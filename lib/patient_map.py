'''
hpo goodness of fit test.
For each gene, test each HPO with P_h >= N (N default 100),
  using the rest HPO terms (not daughter or ancestor) with P_h >= N
  as negative set
'''
from __future__ import print_function, division
import sys
import os
import json
from collections import defaultdict,Counter
import itertools
import copy
import pandas as pd
import numpy as np
import subprocess
from optparse import OptionParser
sys.path.append('commons')
import gnomad_utils
import phenopolis_utils
import patients_variants as PV
import helper

# MONGO = phenopolis_utils.get_mongo_collections()

def read_vcf(vcf_f):
    '''
    read vcf.
    if miss, -1. if wildtype, 0. if het, 1. if hom, 2
    return dataframe
    '''
    result = defaultdict(lambda :
            defaultdict(int)
            )
    header = []
    for row in helper.split_iter(vcf_f):
        if row[:2] == '##': continue
        row = row.split('\t')
        if not header:
            header = row
            continue
        # remove variant with a '*' as alt
        if row[4] == '*': continue
        v_id = '-'.join([row[0],row[1],row[3],row[4]])
        for i in range(9,len(row)):
            gt = row[i].split(':')[0]
            gt = gt.replace('|','/')
            gt = Counter(gt.split('/'))
            if gt['.'] == 2:
                # miss
                result[header[i]][v_id] = -1
            else:
                result[header[i]][v_id] = gt['1']
    result = pd.DataFrame(result)
    return result


'''
get vcf
'''
def get_vcf_df(**kwargs):
    '''
    use bcf tools to subset variants and patients. then according to
    p/v_cutoff to get bad_vs, bad_ps to remove
    '''
    compulsory_keys = {
            'vcf_file',
            'chrom',
            'start',
            'stop',
            'unrelated_file',
            'human_fasta_ref',
            'v_cutoff',
            'gnomad_cutoff',
            'p_cutoff',
            'patient_mini',
            }
    # check args
    helper.check_args(compulsory_keys, kwargs, 'get_vcf_df')
    position = '{chrom}:{start}-{stop}'.format(**kwargs)
    ps1 = subprocess.Popen(('tabix','-h',kwargs['vcf_file'],position),stdout=subprocess.PIPE)
            # subset on unrelated samples, and normalise
    ps2 = subprocess.Popen(('bcftools', 'view', '-Ou', '-S',
        kwargs['unrelated_file'], '-f', 'PASS'),stdin=ps1.stdout,stdout=subprocess.PIPE)
    ps3 = subprocess.Popen(('bcftools', 'norm', '-Ou', '-m', '-any'),
            stdin=ps2.stdout,stdout=subprocess.PIPE)
    normed_vcf = subprocess.check_output(['bcftools', 'norm', '-Ov', '-f',
        kwargs['human_fasta_ref']],stdin=ps3.stdout)
    # get vcf df. genotype -1 = missing, 0 = wildtype, 1 = het, 2 = hom
    genotype_df = read_vcf(normed_vcf)
    # empty vcf? early return
    if genotype_df.empty:
        return None
    # get poorly covered variants and individuals
    # change df to cover_df
    cover_df = genotype_df.copy()
    cover_df[cover_df>=0]=1
    cover_df[cover_df==-1]=0
    pm = cover_df.mean()

    # rid of patients not in patient_mini
    bad_ps = set(pm[pm < kwargs['p_cutoff']].index)
    bad_ps.update(set(pm.index) - set(kwargs['patient_mini'].keys()))
    vm = cover_df.T.mean()
    bad_vs = set(vm[vm < kwargs['v_cutoff']].index)
    # annotate vs with gnomad
    vs = (i for i in vm.index if i not in bad_vs)
    gnomad_freqs = gnomad_utils.overall_freqs(vs, kwargs['gnomad_path'])

    # remove variants with 'SEGDUP' filter. This gives a lot of noise for recessive
    # analysis. For example IGHV3-38 - ENST00000390618, 14-106866588-T-C
    bad_vs.update([
        i for i,v in gnomad_freqs.items()
        if v['filters']['exome'] is not None and 'SEGDUP' in v['filters']['exome']
        or v['filters']['genome'] is not None and 'SEGDUP' in v['filters']['genome']
        ])
    # in fact, many variants have very high af, but 0 hom_f, such as
    # 6-32548641-A-T, which has no 'SEGDUP' filter. Remove those
    # hard filtering for the time being. There might be better ways
    bad_vs.update([
        i for i,v in gnomad_freqs.items()
        if v['gnomad_af'] > 0.01 and v['gnomad_hom_f'] == 0.0
        ])

    # add to bad_vs gnomad_hom_af >= gnomad_cutoff,
    #  and those not covered by gnomad_path
    # Note that if gnomad_hom_af >= gnomad_cutoff, then gnomad_af >= gnomad_cutoff
    #  but not vice versa
    #this = [i for i,v in gnomad_freqs.items()
    #    if v['gnomad_af'] is None or v['gnomad_hom_f'] >= kwargs['gnomad_cutoff']]
    bad_vs.update([i for i,v in gnomad_freqs.items()
        if v['gnomad_af'] is None or v['gnomad_hom_f'] >= kwargs['gnomad_cutoff']])
    vs_count = np.sum(genotype_df[genotype_df>0],axis=1)
    bad_vs.update([i for i in gnomad_freqs if vs_count[i] > 3 and gnomad_freqs[i]['pop_filter']])
    # then drop bad_ps and bad_vs
    genotype_df.drop(bad_vs,inplace=True)
    genotype_df.drop(bad_ps,inplace=True,axis=1)
    return (genotype_df, cover_df, gnomad_freqs)


'''
get coding variants from coding.tsv
'''
def get_coding_variants(f):
    result = []
    with open(f,'r') as inf:
        for line in inf:
            row = line.rstrip().split('\t')
            result.append(tuple(row))
    return set(result)

def remove_cis(patients_variants, genotype_df):
    '''
    # when two variants are in cis and both appear in one patient,
    # discard the variant with lower cadd in that patient
    # example SDK1 ('7-3990565-C-T','7-4014039-C-T')
    # for now, only focus on variants with hom_f < 0.00025

    Note that it should only affect recessive analysis, so make a copy of
    patients -> r_patients
    '''
    patients_variants['r_patients'] = copy.copy(patients_variants['patients'])

    rare_variants = [k for k,v in patients_variants['variants'].items()
            if v['gnomad_hom_f'] < 0.00025]
    sub_df = genotype_df.loc[rare_variants][patients_variants['r_patients'].keys()]
    sub_df[sub_df > 1] = 1
    sub_df[sub_df < 0] = 0
    co_occur = sub_df.dot(sub_df.T)
    for p,v in patients_variants['r_patients'].items():
        bad_vs = set()
        # You don't want to remove hom variants!
        this_rare_variants = Counter([i for i in v if i in rare_variants])
        for pair in itertools.combinations([i for i,vv in this_rare_variants.items() if vv == 1], 2):
            if pair[0] in bad_vs or pair[1] in bad_vs:
                continue
            # more than half? bad!
            bad_cutoff = 1/2
            if min(co_occur[pair[0]][pair[0]],co_occur[pair[1]][pair[1]]) > 2:
                if co_occur[pair[0]][pair[1]] / co_occur[pair[0]][pair[0]] > bad_cutoff or co_occur[pair[0]][pair[1]] / co_occur[pair[1]][pair[1]] > bad_cutoff:
                    if patients_variants['variants'][pair[0]]['cadd'] > patients_variants['variants'][pair[1]]['cadd']:
                        bad_vs.update([pair[1]])
                    else:
                        bad_vs.update([pair[0]])
        for bad_v in bad_vs:
            v.remove(bad_v)

def main(**kwargs):
    '''
    parameters:
     genes: optional
     N (selecting HPO with at least N Ph. affecting both \
          #positive (selecting parental HPO in the positive set \
          #and negative set)
     vcf file location
     gnomad files location
     patient_mini, patient_info, both are json files
     cadd path
     unrelated file used to subset vcf file
     v cutoff and p cutoff are to remove variants and patients with \
          #low coverage over the gene

    returns hpo goodness of fit score, p_g (gnomad_freq
    '''
    # check args
    compulsory_keys = {
        'remove_nc',
        'vcf_file',
        'gnomad_path',
        'cadd_file',
        'patient_mini_file',
        'patient_info_file',
        'human_fasta_ref',
        'unrelated_file',
        'v_cutoff',
        'p_cutoff', # not using this since no phasing is done
        'gnomad_cutoff',
        'gnomad_step',
        'gnomad_path',
        'cadd_step',
        'cadd_min',
        'genon_sum_cutoff_coefficient',
        'cis_gap',
        }
    helper.check_args(compulsory_keys, kwargs, 'main')
    # defaults
    kwargs.setdefault('gene_inheritance_mode',{})
    # get patient_mini and patient_info
    patient_info = helper.get_snapshot(kwargs['patient_info_file'])
    patient_mini = helper.get_snapshot(kwargs['patient_mini_file'])
    # get p_h for all hpos
    phs = helper.get_phs(patient_info)
    # get genes, if not provided. get all gene_ids from mongodb, \
            #if provided, convert to gene_id
    fields = {
            'gene_id':1,
            'gene_name':1,
            '_id':0,
            'chrom':1,
            'start':1,
            'stop':1,
            'xstart':1,
            'xstop':1,
            }
    this = {}
    gnomad_steps = np.arange(
            0,
            kwargs['gnomad_cutoff']+kwargs['gnomad_step'],
            kwargs['gnomad_step']
            )
    cadd_steps = np.arange(kwargs['cadd_min'], 60, kwargs['cadd_step'])

    # for each gene, get all valid variants/patients according to p/v_cutoff,
    # annotate using gnomad
    result = {}
    last_chrom = None
    coding_variants = None
    # then get patients_variants, with variants annotated with
    #  gnomad freqs and cadd
    patients_variants = PV.main(**kwargs)

    # get patient_map for recessive and dominant modes
    args = dict(
            data = patients_variants,
            vcf = patients_variants['cover_df'],
            gnomad_steps = gnomad_steps,
            cadd_steps = cadd_steps,
            cis_gap = kwargs['cis_gap'],
            )
    patient_map = {'r':{},'d':{}}
    # first get mode if provided. Note that the keys could be id
    #  or gene name
    modes = 'rd'

    # get number of patients who carry rare variants when get patient_maps
    NP = {}
    for mode in modes:
        args['mode'] = mode
        M= helper.get_patient_map(**args)
        NP[mode] = len(set(list(itertools.chain.from_iterable(
            [v[0] for k,v in M.items() if k[1] == 0]
            ))))
        # change the keys to a string
        for k,v in M.items():
            patient_map[mode]['{},{}'.format(k[0],k[1])] = [list(v[0]),list(v[1])]
    # for debugging
    #print({k:v[0] for k,v in patient_map['r'].items() if k == (6,0)})
    #for i in range(10):
    #    debug = Counter()
    #    for p in patient_map['d'][(i,0)][0]:
    #        for v in patients_variants['patients'][p]:
    #            if patients_variants['variants'][v]['gnomad_af'] < 0.00025:
    #                debug[v] += 1
    #    print(i)
    #    print('-{}-'.format(len(patient_map['d'][(i,0)][0])))
    #    print(debug.most_common(3))

    '''
    hposet = {'HP:0000512'}
    vs = []
    ind = 4
    for p in patient_map['d'][(ind,0)][0]:
        print(p,patient_mini[p])
        this_hpos = set(patient_info[p])
        if True: #this_hpos & hposet:
            for v in patients_variants['patients'][p]:
                if patients_variants['variants'][v]['gnomad_af'] < 0.00025 and 5*(ind+1) > patients_variants['variants'][v]['cadd'] >= 5*ind:
                    vs.append(v)
    print(Counter(vs))
    sys.exit()
    '''

    return {
            'patient_map':patient_map,
            'patients_variants':patients_variants,
            'NP': NP,
            }

    '''
    for g,v in result.items():
        hn = hpo_name(MONGO['hpo_db'],v['data'].keys())
        print(v['symbol'],v['mode'])
        for k,vv in v['data'].items():
            print(hn[k],vv)
    sys.exit()
    '''

    return result


if __name__ == '__main__':
    # in the end some of the args have to go to the config
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--range",
                      dest="range",
                      help="which genome range to process?")
    parser.add_option("--output",
                      dest="output",
                      help="output file name?")
    (options, args) = parser.parse_args()
    args = dict(
        #genes = ('ABCA4','CERKL','SCN1A','GUCY2D','USH2A','PROM1','TERT','CNGB1','CRB1','IMPG2','RPGR','SDK1'),
        range = options.range,
        output = options.output,
    )
    # update args with commons.cfg
    args.update(phenopolis_utils.OFFLINE_CONFIG['generic'])
    args.update(phenopolis_utils.OFFLINE_CONFIG['phenogenon'])
    result = main(**args)
    # write everything to output
    with open(args['output'], 'wb') as outf:
        json.dump(result, outf)
    print('==done==')
