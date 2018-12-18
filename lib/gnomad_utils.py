from __future__ import print_function, division
import sys
import pysam
import os
sys.path.append('../../commons')
import common_utils

class BadVariantException(Exception): pass
#path = '/cluster/project8/vyp/gnomad_data'
VALID_CHROMOSOMES = [str(i) for i in range(1,23)] + ['X','Y']
POPS = ('AFR','NFE','AMR','EAS','ASJ','SAS','OTH','FIN')
'''
coverage, look at different parts of the ref
'''
def coverage(v,path_to_gnomad,mode='exome'):
    v = common_utils.clean_variant(v)
    # pysam does not support header yet. hard code it
    header = ['chrom','pos','mean','median',1,5,10,15,20,25,30,50,100,]
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'coverage','exomes','exacv2.chr'+chrom+'.cov.txt.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'coverage','genomes','gnomad.chr'+chrom+'.cov.txt.gz')
    else:
        msg = "mode only accepts 'exome' or 'genome'"
        raise ValueError(msg)
    tb = pysam.TabixFile(file)
    start = int(pos) - 1
    end = int(pos)
    if len(ref) != len(alt):
        # indel. since cleaned, only need to shit 1 position
        start += 1
        end += len(ref) - 1

    rs = tb.fetch(chrom, start, end)
    result = {}
    for line in rs:
        r = line.rstrip().split('\t')
        this = {a:b for a,b in zip(header,r)}
        result[int(this['pos'])] = this
    if not result:
        # not covered
        return None
    # if indel, check if all ref positions are covered
    if len(ref) == len(alt):
        return result
    else:
        for p in range(start+1,end+1):
            if p not in result: return None
        return result

'''
exome freqs
'''
def freqs(v,path_to_gnomad,mode='exome'):
    v = common_utils.clean_variant(v)
    # pysam does not support header yet. hard code it
    header = ['chrom','pos','id','ref','alt','quality','filter','info']
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'vcf','exomes','gnomad.exomes.r2.0.1.sites.vcf.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'vcf','genomes','gnomad.genomes.r2.0.1.sites.'+chrom+'.vcf.gz')
    tb = pysam.TabixFile(file)
    records = tb.fetch(chrom, int(pos)-1, int(pos))

    for line in records:
        if not line: return None
        r = line.rstrip().split('\t')
        data = {a:b for a,b in zip(header,r)}
        # find the variant
        g_alts = data['alt'].split(',')
        alt_ind = None
        for ind,this_alt in enumerate(g_alts):
            v_id = common_utils.clean_variant('-'.join([data['chrom'],data['pos'],data['ref'],this_alt]))
            if v_id == v:
                alt_ind = ind
        if alt_ind == None:
            continue

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
        info_dict['filter'] = data['filter']
        return info_dict
                
    return None


'''
simple query on overall allele freq (or homozygote frequency). if covered, return at least 0. if not, return None. 
'''
def overall_freqs(vs,path_to_gnomad):
    result = {}
    null = {
        'gnomad_af': None,
        'gnomad_ac': None,
        'gnomad_hom_f': None,
        'gnomad_hom_c': None,
        'gnomad_hemi_c': None,
        'filters':{'exome':None,'genome':None},
        'pop_filter':[],
        'most_freq_pops':[],
    }
    for v in vs:

        if v.split('-')[0] not in VALID_CHROMOSOMES:
            result[v] = null
            continue
        covs = {
            'exome': coverage(v,path_to_gnomad,mode='exome'),
            'genome':coverage(v,path_to_gnomad,mode='genome'),
        }
        fs = {
            'exome': freqs(v,path_to_gnomad,mode='exome'),
            'genome':freqs(v,path_to_gnomad,mode='genome'),
        }

        if not fs['exome'] and not fs['genome'] and not covs['exome'] and not covs['genome']:
            result[v] = null
            continue
        ac = hom_c = af = hom_f = an = 0.
        hemi_c = None
        pop_filter = []
        filters = {'exome':None,'genome':None}
        # also check population frequencies to remove any variants 
        # with big af(>0.01)/hom_f(0) discrepancy, such as 1-144931607-C-T
        pops = {p:{'Hom':0,'Hemi':0,'AC':0,'AN':0} for p in POPS}
        for m in ['exome', 'genome']:
            if fs[m]:
                ac += fs[m]['AC']
                hom_c += fs[m]['Hom']
                an += fs[m]['AN']
                if 'Hemi' in fs[m]:
                    hemi_c = hemi_c + fs[m]['Hemi'] if hemi_c != None else fs[m]['Hemi']
                filters[m] = fs[m]['filter']
                for p in POPS:
                    for kk in pops[p]:
                        try:
                            this_c = fs[m].get('{}_{}'.format(kk,p), 0)
                            # sometimes on X this_c is a dot
                            if this_c == '.':
                                this_c = 0
                            pops[p][kk] += this_c
                        except TypeError:
                            print(v)
                            print(p,kk)
                            print(pops[p][kk])
                            print(fs[m].get('{}_{}'.format(kk,p), 0))
                            raise

        max_pop = ([], -1)
        for p in pops:
            if pops[p]['AC']:
                pop_af = pops[p]['AC'] / pops[p]['AN']
                if pop_af:
                    if pop_af > max_pop[1]:
                        max_pop = ([p], pop_af)
                    elif pop_af == max_pop[1]:
                        max_pop[0].append(p)
                if pop_af > 0.01 and pops[p]['Hemi'] == 0 and pops[p]['Hom'] == 0:
                    pop_filter.append(p)

        if ac: af = ac / an
        if hom_c: hom_f = hom_c * 2 / an

        result[v] = {
            'gnomad_af':af,
            'gnomad_ac':ac,
            'gnomad_hom_f':hom_f,
            'gnomad_hom_c':hom_c,
            'gnomad_hemi_c':hemi_c,
            'gnomad_an':an,
            'filters':filters,
            'pop_filter':pop_filter,
            'most_freq_pops':max_pop[0],
        }

    return result
