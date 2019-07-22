# some helpers
from __future__ import print_function, division
import os
import re
from collections import Counter
import gzip
import subprocess
import pysam
import sys
sys.path.append('../phenogenon')
from pyhpo import Hpo
import gnomad_utils
# python3 is configparser, python2 is ConfigParser
try:
    import ConfigParser
except ModuleNotFoundError:
    import configparser as ConfigParser

CHROMOSOMES = [str(i) for i in range(23)] + ['X']
def parse_config(config_file):
    '''
    parse config file, and make config global. If test, set DB_HOST as 'localhost'
    borrowed from ../phenogenon/helper.py
    '''
    # return {'section':{'key1':'value1'...},...}
    config = ConfigParser.ConfigParser()
    config.read(config_file)
    result = {}
    for section in config.sections():
        options = config.options(section)
        result[section] = {}
        for option in options:
            value = config.get(section, option)
            if value in ('true', 'True'):
                value = True
            elif value in ('false', 'False'):
                value = False
            elif value.isdigit():
                value = int(value)
            elif ',' in value:
                value = re.split(r', *', value)
                try:
                    value = [float(i) for i in value]
                except ValueError:
                    pass
            else:
                try:
                    value = float(value)
                except ValueError:
                    pass
            result[section][option] = value
    return result


def parse_vcf_header(line):
    '''
    given a vcf ## line, return a header format
    input: string
    return: [colname1, colname2]
    '''
    result = line.split('"')[1].split('Format: ')[1].split('|')
    return result

def parse_vcf(vcf_file):
    '''
    parse vcf to get interesting stuff
    note that multiallelic are coded in separate lines
    input: a.vcf.gz
    gnomAd and cadd might be redundent as they will be stored separately
    return: {
        entrez_id: {
            'symbol',
            'entrez_id',
            'omim:hpo':{omim:Counter([hpos])}
            'omim':[omimid],
            'start': not gene start, but start of the gene in this vcf
            'variants':[{
                variant_id,
                gnomad_af,
                gnomad_hf,
                cadd,
                hpo: [HP1]
            }]
        }
    }

    NOTE that each variant may sit on multiple genes.
        remove 100 percent overlapping genes (keep the first only)
        it can be told by identical starts and stops
    hpos can be repetitive for a given gene.
    variants that are not covered by gnomad will be removed
    '''
    import gzip
    result = {}

    # for getting gnomad
    popf_header = []
    # column header
    header = []
    # past genes to check if they have identical start/stop
    past_genes = []
    percent_100_overlap_flag = False
    with gzip.open(vcf_file, 'rt') as inf:
        for line in inf:
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=POPF'):
                    popf_header = parse_vcf_header(line)
                elif line.startswith('##INFO=<ID=CSQ'):
                    csq_header = parse_vcf_header(line)
                else:
                    continue
            if line.startswith('#'):
                header = line[1:].rstrip().split('\t')
                continue
            row = dict(zip(header, line.rstrip().split('\t')))

            # get variant id
            variant_id = '-'.join([row['CHROM'], row['POS'], row['REF'], row['ALT']])
            # get info
            genes = []
            hpos,omim = set(),set()
            gnomad_af, gnomad_hf, cadd = None, None, None
            for info in row['INFO'].split(';'):
                if info.startswith('GENEINFO'):
                    # get genes
                    # entrez id is used in HPO association file, so get it
                    for gene_field in info.split('=')[1].split('|'):
                        symbol,entrez_id = gene_field.split(':')
                        genes.append({
                                'symbol': symbol,
                                'entrez_id': entrez_id
                        })
                elif info.startswith('POPF'):
                    # get gnomads
                    rec = dict(zip(popf_header, info.split('=')[1].split('|')))
                    gnomad_af = float(rec['gnomad_af']) if rec['gnomad_af'] else None
                    gnomad_hf = float(rec['gnomad_hom_f']) if rec['gnomad_hom_f'] else None
                elif info.startswith('CADD'):
                    # get cadd
                    cadd = float(info.split('=')[1])
                elif info.startswith('CLNDISDB'):
                    # get hpos and OMIM
                    for annos in info.split('=')[1].split('|'):
                        for anno in annos.split(','):
                            rec = anno.split(':',1)
                            if rec[0] == 'Human_Phenotype_Ontology':
                                hpos.add(rec[1])
                            elif rec[0] == 'OMIM':
                                omim.add(rec[1])

            if gnomad_af is None:
                continue
            # check overlapping genes
            past_genes = check_del_overlap_genes(result, genes, past_genes)

            # populate result
            variant = {
                'variant_id': variant_id,
                'gnomad_af': gnomad_af,
                'gnomad_hf': gnomad_hf,
                'cadd': cadd,
                'hpo': hpos,
                'omim': omim
            }
            for gene in genes:
                if gene['entrez_id'] not in result:
                    result[gene['entrez_id']] = {
                        'symbol': gene['symbol'],
                        'entrez_id': gene['entrez_id'],
                        'omim:hpo': {o:Counter(hpos) for o in omim},
                        'omim': list(omim),
                        'variants':[variant],
                        'start': int(row['POS']),
                    }
                else:
                    result[gene['entrez_id']]['variants'].append(variant)
                    result[gene['entrez_id']]['omim:hpo'].update({o:Counter(hpos) + result[gene['entrez_id']]['omim:hpo'].get(o,Counter()) for o in omim})
                    result[gene['entrez_id']]['omim'].extend(list(omim))

    # check overlapping genes
    past_genes = check_del_overlap_genes(result, [], past_genes)
    return result

def counts_to_possibility_bins(counts):
    '''
    input: Counter({a: 10, b: 5, c: 1})
    output: {a: (0, 10/(10+5+1)), b: (10/(10+5+1), (10+5)/(10+5+1)), c: ((10+5)/(10+5+1), 1)}
    '''
    result = {}
    S = sum(counts.values())
    accu = 0
    for k,v in counts.items():
        result[k] = (accu/S, (accu+v)/S)
        accu += v
    return result

def find_gene_range_in_gtf(gene, gtf):
    p1 = subprocess.Popen(('zcat', gtf), stdout=subprocess.PIPE)
    try:
        output = subprocess.check_output(('grep', gene), stdin=p1.stdout).decode('utf-8')
    except subprocess.CalledProcessError:
        print('!!!!'+gene)
        return {'chrom':'0', 'start':-1, 'stop':-1}
    for line in split_iter(output):
        row = line.rstrip().split('\t')
        if row[2] != 'gene' or row[0] not in CHROMOSOMES:
            continue
        # get gene name
        for i in row[-1].split(';'):
            field = i.split()
            if field and field[0] == 'gene_name' and field[1].strip('"') == gene:
                return {'chrom': row[0], 'start': int(row[3]), 'stop': int(row[4])}

def split_iter(string):
    '''
    iter through string separated by \n
    '''
    return (x.group(0) for x in re.finditer(r"[^\n]+", string))


def check_del_overlap_genes(data, genes, past_genes):
    # check overlapping genes
    # genes: this vcf line
    # past_genes: last vcf line
    done_genes = sorted([gene for gene in past_genes if gene not in [g['entrez_id'] for g in genes]], key=lambda x: data[x]['start'])
    if len(done_genes) > 1:
        replicates = []
        for i in range(len(done_genes)-1):
            if i in replicates:
                continue
            for j in range(i+1, len(done_genes)):
                if data[done_genes[i]]['start'] == data[done_genes[j]]['start']:
                    replicates.append(j)
        for i in replicates:
            del data[done_genes[i]]
    
    past_genes = [gene['entrez_id'] for gene in genes]
    return past_genes

def parse_gene_ranges(infile):
    '''
    parse clinvar pathogenic gene ranges
    if file doesn't exist, raise
    '''
    if not os.path.isfile(infile):
        msg = 'gene range file does not exist. Please use get_gene_range.py to get the file before using this script'
        raise ValueError(msg)
    result = {}
    header = []
    with open(infile, 'rt') as inf:
        for line in inf:
            row = line.rstrip().split('\t')
            if not header:
                header = row
                continue
            rec = dict(zip(header, row))
            result[rec['entrez_id']] = {
                    'chrom': rec['chrom'],
                    'start': int(rec['start']),
                    'stop': int(rec['stop']),
            }
    return result

def get_gene_variants_from_vcf(gene_range, hpos, params):
    '''
    get variants from vcf given ranges.
    remove patients with associated HPOs
    annotate with gnomad and cadd
    and internal AF (used for sampling)
    dumping individual information
    only PASS variants
    return:
    [{
        variant_id,
        internal_af,
        gnomad_af,
        gnomad_hf,
        cadd
    }]
    '''
    # some cadd is missing
    missing_cadd = {}
    with gzip.open('no_cadd.cadd.gz', 'rt') as inf:
        for line in inf:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            v_id = '-'.join(row[:4])
            missing_cadd[v_id] = float(row[-1])

    vcf_common_header = {'CHROM','POS','ID','REF','ALT','QUAL','FORMAT','INFO','FILTER'}
    hpos = set(hpos)
    human_ref_tbx = pysam.FastaFile(params['generic']['human_fasta_ref'])
    # get unrelated patients who do not have associated HPOs
    samples = set()
    with open(params['generic']['patient_info_file'],'rt') as inf:
        header = []
        for line in inf:
            row = line.rstrip().split('\t')
            if not header:
                header = row
                continue
            row = dict(zip(header, row))
            if row['unrelated'] == '0':
                # related, ignore
                continue
            patient_hpos = set(row['HPO'].split(','))
            if not patient_hpos & hpos:
                samples.add(row['#p_id'])

    # get sample vcf tabix and vcf header
    vcf_file = params['generic']['ucl_vcf'].format(params['generic']['release'],params['generic']['release'],gene_range['chrom'])
    vcf_tbx = pysam.TabixFile(vcf_file)
    vcf_header = []
    with gzip.open(vcf_file, 'rt') as inf:
        for line in inf:
            if line.startswith('##'):
                continue
            vcf_header = line[1:].rstrip().split('\t')
            break

    # get coding regions, considering padding
    exons = get_exons(gene_range, params)
    # get variants
    result = {}
    for line in vcf_tbx.fetch(gene_range['chrom'],gene_range['start'],gene_range['stop']):
        row = dict(zip(vcf_header,line.split('\t')))
        # filter not pass? pass
        if row['FILTER'] != 'PASS':
            continue
        fmt = row['FORMAT'].split(':')

        # get genotypes
        genotypes = {}
        for key,val in row.items():
            if key in vcf_common_header:
                continue
            if key not in samples:
                continue
            gt = dict(zip(fmt, val.split(':')))['GT']
            if '.' in gt:
                continue
            genotypes[key] = Counter(gt)
        
        if not genotypes:
            continue
        for ind, alt in enumerate(row['ALT'].split(',')):
            if alt == '<*:DEL>':
                continue
            # get variant id
            # note that alt could be *
            variant_id = clean_variant('-'.join([row['CHROM'],row['POS'],row['REF'],alt.replace('*','-')]),human_ref_tbx)
            variant_id = find_leftmost_synonymous_variant(variant_id, human_ref_pysam=human_ref_tbx)
            # in exon?
            if is_noncoding(variant_id, exons):
                continue
            # get genotypes
            ac = 0
            # an = len(genotypes)
            for sample, gt in genotypes.items():
                ac += gt[str(ind+1)]
            af = ac / len(genotypes) / 2
            result[variant_id] = {'af':af, 'gnomad_af':None, 'gnomad_hf': None, 'cadd':None}

    # get gnomad
    gnomad_freqs = gnomad_utils.overall_freqs(list(result.keys()), params['generic']['gnomad_path'])
    # get cadd
    cadd_file = params['generic']['cadd_file'].format(params['generic']['release'], gene_range['chrom'])
    cadds = get_cadd(gene_range, set(result.keys()), cadd_file,human_ref_tbx)
    to_be_deleted = []
    for variant in result:
        # if variant not covered by gnomad, remove
        if gnomad_freqs[variant]['gnomad_af'] is None:
            to_be_deleted.append(variant)
            continue
        if variant in cadds:
            cadd = cadds[variant]
        else:
            cadd = missing_cadd[variant]
        result[variant].update({
            'gnomad_af': gnomad_freqs[variant]['gnomad_af'],
            'gnomad_hf': gnomad_freqs[variant]['gnomad_hom_f'],
            'cadd': cadd
        })
    for variant in to_be_deleted:
        del result[variant]
    return result
            
def get_cadd(gene_range, variants, cadd_file, human_ref_pysam):
    result = {}
    cadd_tbx = pysam.TabixFile(cadd_file)
    try:
        iterator = cadd_tbx.fetch(gene_range['chrom'], gene_range['start'], gene_range['stop'])
    except ValueError:
        print('oops')
        return result
    for line in iterator:
        row = line.split('\t')
        # some variants are not left-aligned
        v_id = find_leftmost_synonymous_variant('-'.join(row[:4]), human_ref_pysam=human_ref_pysam)
        if v_id in variants:
            result[v_id] = float(row[-1])
    return result

def clean_variant(v,human_ref_pysam):
    # sometimes variant has funny format, which has more - than expected, such as 1-117122294---TCT.
    #  use find_bases to fill in the gap if human_ref_pysam is not provided
    if v.count('-') == 4:
        if v[-1] == '-':
            # deletion
            chrom,pos,ref,rubbish,rubbish = v.split('-')
            pos = int(pos)-1
            common_base = human_ref_pysam.fetch(chrom, pos-1, pos)
            ref = common_base + ref
            alt = common_base
        else:
            # insertion
            chrom,pos,ref,rubbish,alt = v.split('-')
            pos = int(pos)
            common_base = human_ref_pysam.fetch(chrom, pos-1, pos)
            ref = common_base
            alt = common_base + alt
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

def find_leftmost_synonymous_variant(variant, padding=200, human_ref_pysam=None):
    '''
    Only necessary for indel!
    find all synonymous variants given variant
    padding is how far you would like to search left and right of the change
    if human_ref_pysam is None, use find_base to query ensembl
    '''
    mode,pattern = None,None
    chrom, pos, ref, alt = variant.split('-')
    pos = int(pos)+1
    # removing commong base
    if ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    if len(ref) and not alt:
        pattern = ref
        mode = 'del'
    elif len(alt) and not ref:
        pattern = alt
        mode = 'in'
    else:
        # it's not indel, return
        return variant
    string = human_ref_pysam.fetch(chrom, pos-padding-1, pos+len(pattern)+padding-1)
    ind = find_start_of_repeat(string, padding, len(pattern))
    new_pos = pos - padding + ind
    missing_base = string[ind]
    pattern = string[ind+1:ind+len(pattern)+1]

    if mode == 'del':
        return '-'.join([chrom, str(new_pos), missing_base+pattern, missing_base])
    elif mode == 'in':
        return '-'.join([chrom, str(new_pos), missing_base, missing_base+pattern])
    else:
        msg = 'Cannot derive mode!'
        raise ValueError(msg)

def find_start_of_repeat(string, start, length):
    '''
    string: GCAGAGAGAG
    start: 5
    length: 2 #GA
    return 1 # the repeat starts from after 1:AG
    ===
    if no repeat, return start
    '''
    ind = start
    result = string[:start]+string[start+length:]
    while ind >= 0:
        ind -= 1
        if string[:ind] + string[ind+length:] != result:
            return ind
    return ind

def get_exons(gene_range, params):
    '''
    get exons, including paddings
    '''
    # get range
    # get all exons from gtf
    exons = []
    tb = pysam.TabixFile(params['generic']['gtf'])
    for line in tb.fetch(gene_range['chrom'], gene_range['start'], gene_range['stop']):
        row = line.rstrip().split('\t')
        # if exon?
        if row[2] != 'CDS' or row[1] != 'protein_coding':
            continue
        start = int(row[3]) - params['generic']['exon_padding']
        stop = int(row[4]) + params['generic']['exon_padding']
        exons.append({
            'start': start,
            'stop': stop,
            'len': stop - start + 1,
        })
    return exons


def is_noncoding(variant, exons):
    # get list of variants that are non-coding
    chrom, pos, ref, _ = variant.split('-')
    start = int(pos)
    stop = start + len(ref) - 1
    #start and end overlaps with any exons?
    for exon in exons:
        if exon['start'] > stop:
            break
        # overlap?
        if len(ref) + exon['len'] > max(exon['stop'], stop) - min(exon['start'], start) + 1:
            return False
    return True

def simulate_samples(genes, params):
    '''
    simulate control, disease, noise samples
    '''
    import uuid
    # get generic, starts with simu_
    result = [f'simu_{uuid.uuid1()}' for i in range(params['simulation']['control_pool_size'])]
    # get omim_moi
    omim_moi = 
    return result

def get_omim_moi(F):
    '''
    from the omim_moi file, get omim moi
    '''
    result = {}
    