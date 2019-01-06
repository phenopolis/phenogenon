'''
some functions used by Genon
'''
from collections import Counter, defaultdict
import json
import itertools
import re
import pandas as pd

def get_snapshot(f):
    '''
    get patient hpo terms as a dict
    '''
    dt = {}
    with open(f,'r') as inf:
        for row in inf:
            if row[0] == '#': continue
            row = row.rstrip().split('\t')
            if row[1] == '0': continue
            dt[row[0]] = row[2].split(',')
    return dt

def get_cover_from_vcf(vcf_file):
    '''
    # get coverage information from a vcf file. return a dataframe, 
    # with rows = variant_ids, and columns = individuals.
    # false = not covered, true = covered
    '''
    # count number of lines starting with ##
    number_of_skipping_rows = 0
    with open(vcf_file,'r') as inf:
        for row in inf:
            if row[:2] == '##':
                number_of_skipping_rows += 1
            else:
                break
    df = pd.read_table(vcf_file,skiprows=range(number_of_skipping_rows))
    # variant_id
    df['variant_id'] = df.apply(lambda x: '-'.join([str(x['#CHROM']),str(x['POS']),x['REF'],x['ALT']]), axis=1)
    # columns to drop
    columns_to_drop = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','FORMAT','INFO']
    df = df.drop(columns_to_drop,axis=1).set_index('variant_id')
    # convert genotype
    df = df.applymap(lambda x: '.' not in x.split(':')[0])
    return df

def get_negative_hpos(ps,hpos,hpo_db):
    result = []
    for h1 in hpos:
        bad = 0
        for h2 in ps:
            A = h1 in [i['id'][0] for i in get_hpo_ancestors(hpo_db, h2)]
            B = h2 in [i['id'][0] for i in get_hpo_ancestors(hpo_db, h1)]
            if A or B:
                bad = 1
                break
        if not bad:
            result.append(h1)
    return hpo_minimum_set(hpo_db,result)

# remove deep intronic and other non-coding variants
# get bad variants first
# if an intronic variant is close to exon (like 4 bases), don't include it
def extract_nc_from_vcf(vcf_file, closeness):
    class goodException(Exception): pass
    result = []
    with open(vcf_file,'r') as inf:
        info_header,general_header,info_ind = None,None,None
        for line in inf:
            if line[:2] == '##':
                # parse INFO header
                if not line.startswith('##INFO=<ID=CSQ'): continue
                info_header = line.split('"')[1].split()[-1].split('|')
                continue
            row = line.rstrip().split('\t')
            if line[0] == '#':
                general_header = row
                info_ind = general_header.index('INFO')
                continue
            variant = '-'.join(row[:2]+row[3:5])
            info = row[info_ind].split(';')
            vep = [i for i in info if i.startswith('CSQ=')][0].split(',')
            try:
                for anno in vep:
                    fields = dict(zip(info_header, anno.split('|')))
                    if fields['IMPACT'] != 'MODIFIER':
                        raise goodException
                    if fields['Consequence'] == 'intron_variant':
                        coding_change = fields['HGVSc'].split(':')[1]
                        p = re.compile('.+[+-](\d+).+')
                        awayness = int(p.search(coding_change).groups()[0])
                        if awayness <= closeness:
                            raise goodException
            except goodException:
                continue
            result.append(variant)

    return result

def remove_noncoding(data, nc_variants):
    nc_variants = set(nc_variants)
    for v in data['variants'].keys():
        if v in nc_variants:
            data['variants'].pop(v)
    for p in data['patients'].keys():
        for v in list(data['patients'][p]['variants']):
            if v in nc_variants:
                data['patients'][p]['variants'].remove(v)

'''
coverage
'''
def coverage(v,path_to_gnomad,mode='exome'):
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','mean','median',1,5,10,15,20,25,30,50,100,]
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'exomes','coverage','exacv2.chr'+chrom+'.cov.txt.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'genomes','coverage','gnomad.chr'+chrom+'.cov.txt.gz')
    else:
        msg = "mode only accepts 'exome' or 'genome'"
        raise ValueError(msg)
    tb = tabix.open(file)
    r = tb.query(chrom, int(pos)-1, int(pos))
    r = list(r)
    if not r:
        # not covered
        return None
    r = r[0]
    return {a:b for a,b in zip(header,r)}

'''
exome freqs
'''
def freqs(v,path_to_gnomad,mode='exome'):
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','id','ref','alt','quality','filter','info']
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'exomes','vcf','gnomad.exomes.r2.0.1.sites.vcf.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'genomes','vcf','gnomad.genomes.r2.0.1.sites.'+chrom+'.vcf.gz')
    tb = tabix.open(file)
    records = tb.query(chrom, int(pos)-1, int(pos))

    found = False
    for r in records:
        if not r: return None
        data = {a:b for a,b in zip(header,r)}
        # ref ok?
        if data['ref'] != ref:
            continue

        found = True
        # find alt
        g_alts = data['alt'].split(',')
        if alt not in g_alts:
            return None
        alt_ind = g_alts.index(alt)

        # parse info
        # no need for annotation?

        info = data['info'].split(';CSQ=A')[0] # 1 for annotation
        info = info.split(';')
        info_dict = {}
        for i in info:
            if not '=' in i: continue
            a,b = i.split('=')
            b = b.split(',')
            ind = min(alt_ind,len(b)-1)
            b = b[ind]
            # convert to number if possible
            try:
                if '.' in b:
                    b = float(b)
                else:
                    b = int(b)
            except ValueError:
                pass
            info_dict[a] = b
        info_dict['FILTER'] = info_dict['AS_FilterStatus']
        return info_dict


    if not found: return None

'''
given ids, return id:name dict
'''
def hpo_name(hpo_db, ids):
    records = hpo_db.hpo.find({'id':{'$in':ids}},{'id':1,'name':1,'_id':0})
    result = {}
    for r in records:
        result[r['id'][0]] = r['name'][0]
    return result
'''
some hpo ids are obsolete
this is a copy of lookups.py's replace_hpo. lookups.py is not working atm
'''
def replace_hpo(hpo_db, hpo):
    # some hpo_ids are obsolete.
    record = hpo_db.hpo.find_one({'id':hpo[0]})
    if not record:
        print ('no record in replace_hpo')
        print (hpo)
    if 'replaced_by' in record:
        new = hpo_db.hpo.find_one({'id':record['replaced_by'][0]})
        return [new['id'][0], new['name'][0]]
    else:
        return hpo
'''
get minimised set of HPO terms
'''
def hpo_minimum_set(hpo_db, hpo_ids):
    '''
    minimize the hpo sets
    results = {'HP:0000505': [ancestors]}
    '''
    hpo_ids = list(set(hpo_ids))
    results = dict([(hpo_id, [ h['id'][0] for h in get_hpo_ancestors(hpo_db, hpo_id)],) for hpo_id in hpo_ids])
    # minimise
    bad_ids = []
    for i in range(len(hpo_ids)):
        for j in range(i+1,len(hpo_ids)):
            if hpo_ids[i] in results[hpo_ids[j]]:
                # i is j's ancestor, remove
                bad_ids.append(hpo_ids[i])
                break
            if hpo_ids[j] in results[hpo_ids[i]]:
                # j is i's ancestor, remove
                bad_ids.append(hpo_ids[j])
    return list(set(hpo_ids) - set(bad_ids))
'''
get all ancestor nodes of a given hpo_id.
'''
def get_hpo_ancestors(hpo_db, hpo_id):
    """
    Get HPO terms higher up in the hierarchy.
    """
    h=hpo_db.hpo.find_one({'id':hpo_id})
    #print(hpo_id,h)
    if 'replaced_by' in h:
        # not primary id, replace with primary id and try again
        h = hpo_db.hpo.find_one({'id':h['replaced_by'][0]})
    hpo=[h]
    if 'is_a' not in h: return hpo
    for hpo_parent_id in h['is_a']:
        #p=hpo_db.hpo.find({'id':hpo_parent_id}):
        hpo+=list(itertools.chain(get_hpo_ancestors(hpo_db,hpo_parent_id)))
    #remove duplicates
    hpo={h['id'][0]:h for h in hpo}.values()
    return hpo

# read gene and cadd files
def read_files(json_file):
    with open(json_file,'r') as inf:
        raw = json.load(inf)
    #patient_df = pd.read_table('data/private/hpo/patients_hpo_snapshot_2017-May.tsv')
    #patient_df = patient_df[patient_df['unrelated'] == 1].set_index('#p_id')
    normalise_variants(raw)
    return raw

# get Ph for all HPO terms
def get_phs(p_info):
    result = Counter()
    for k,v in p_info.items():
        result.update(v)
    return result

# remove artefacts in gene files
def remove_artefacts(raw,bad_vs):
    for v in bad_vs:
        del raw['variants'][v]
    bad_p = []
    for k,v in raw['patients'].items():
        v['variants'] = [i for i in v['variants'] if i not in bad_vs]
        if not v['variants']:
            bad_p.append(k)
    for i in bad_p:
        del raw['patients'][i]

# a helper function, to normalise variant.
def clean_variant(v):
    # sometimes variant has funny format, which has more - than expected, such as 1-117122294---TCT.
    #  use find_bases to fill in the gap
    if v.count('-') == 4:
        m = 'Cannot handle missing ref/alt'
        raise ValueError(m)
    else:
        chrom,pos,ref,alt = v.split('-')
        pos = int(pos)
    if len(ref) < len(alt):
        ran = range(len(ref))
    else:
        ran = range(len(alt))
    # insert
    for e in ran:
        ref_e = len(ref) - e - 1
        alt_e = len(alt) - e - 1
        if ref[ref_e] != alt[alt_e]: break
    for b in ran:
        if ref[b] != alt[b] or len(ref[b:ref_e+1]) == 1 or len(alt[b:alt_e+1]) == 1:
            break
    return '-'.join([chrom,str(pos+b),ref[b:ref_e+1],alt[b:alt_e+1]])

# cleanse raw variants
def normalise_variants(data):
    for k,v in data['variants'].items():
        data['variants'][clean_variant(k)] = data['variants'].pop(k)
    for k,v in data['patients'].items():
        v['variants'] = [clean_variant(i) for i in v['variants']]

# remove variants with gnomad_af = None or filter != 'PASS'
def cleanse_variants(data):
    bad = []
    for k,v in data['variants'].items():
        if v['gnomad_af'] == None or v['filter'] != 'PASS':
            bad.append(k)
    for k in bad:
        del data['variants'][k]
    bad = []
    for k,v1 in data['patients'].items():
        this_bad = []
        for v2 in v1['variants']:
            if v2 not in data['variants']:
                this_bad.append(v2)
        for i in this_bad:
            v1['variants'].remove(i)
        if not v1['variants']:
            bad.append(k)
    for k in bad:
        del data['patients'][k]

