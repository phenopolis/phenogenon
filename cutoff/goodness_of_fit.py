'''
hpo goodness of fit test.
For each gene, test each HPO with P_h >= N (N default 100), 
  using the rest HPO terms (not daughter or ancestor) with P_h >= N 
  as negative set
'''
from __future__ import print_function, division
import tabix
import gnomad_utils
import subprocess
from optparse import OptionParser
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
import json
from optparse import OptionParser
from collections import defaultdict,Counter
import pandas as pd
import numpy as np
import helper
from sklearn.cluster import KMeans

MONGO = phenopolis_utils.get_mongo_collections()

def hpo_name(hpo_db, ids):
    '''
    this is only for debugging purpose
    '''
    records = hpo_db.hpo.find({'id':{'$in':ids}},{'id':1,'name':1,'_id':0})
    result = {}
    for r in records:
        result[r['id'][0]] = r['name'][0]
    return result

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


def get_patients_variants(**kwargs):
    '''
    change genotype_df's data format into a dictionary, so gnomad
    and cadd can be added.
    return result
    {
        patients:{p:[v1,v2]}
        variants:{v1:{
            gnomad_af:
            gnomad_hom_f:
            cadd:
        }}
    }
    '''
    # define result
    result = {
            'patients':defaultdict(list),
            'variants':{}
            }

    compulsory_keys = {
            'genotype_df',
            'gnomad_freqs',
            }
    # check args
    helper.check_args(compulsory_keys, kwargs, 'get_patients_variants')
    # change the data structure to the result format
    for col in list(kwargs['genotype_df']):
        this = kwargs['genotype_df'][col]
        # there might be a way to merge the following two steps into one
        # + het
        result['patients'][col].extend(this[this==1].index)
        # + hom
        result['patients'][col].extend(list(this[this==2].index)*2)
    result['variants'] = {i:{
        'gnomad_af':kwargs['gnomad_freqs'][i]['gnomad_af'],
        'gnomad_hom_f':kwargs['gnomad_freqs'][i]['gnomad_hom_f'],
        'cadd':None,
        } for i in kwargs['genotype_df'].index}
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
    # add to bad_vs gnomad_hom_af >= gnomad_cutoff, 
    #  and those not covered by gnomad_path
    # Note that if gnomad_hom_af >= gnomad_cutoff, then gnomad_af >= gnomad_cutoff
    #  but not vice versa
    this = [i for i,v in gnomad_freqs.items()
        if v['gnomad_af'] is None or v['gnomad_hom_f'] >= kwargs['gnomad_cutoff']]
    bad_vs.update([i for i,v in gnomad_freqs.items()
        if v['gnomad_af'] is None or v['gnomad_hom_f'] >= kwargs['gnomad_cutoff']])
    # then drop bad_ps and bad_vs
    genotype_df.drop(bad_vs,inplace=True)
    genotype_df.drop(bad_ps,inplace=True,axis=1)
    return (genotype_df, cover_df, gnomad_freqs)

def get_hgf(**kwargs):
    '''
    calculate HPO goodness of fit
    '''
    compulsory_args = {
        'phenogenons',
        'positive_hpos',
        'negative_hpos',
        'gamma_k',
        'gamma_scale',
        'patient_map',
        }
    helper.check_args(compulsory_args, kwargs, 'get_hgf')
    result = {}
    diff = helper.get_diff(**kwargs)
    for hpo in diff:
        result[hpo] = sum([i[0] for i in diff[hpo]['diff_p'] if i[0]>0])
    return result

'''
given chromosome and db, return gene_ranges
'''
def get_chrom_genes(chrom,fields, db):
    # give chrom numbers, get all genes on them

    chrom = str(chrom)
    if chrom not in phenopolis_utils.VALID_CHROMOSOMES:
        raise ValueError('Error: %s is not a valid chromosome!' % chrom)
    gene_ranges = db.genes.find({'chrom':chrom},fields,no_cursor_timeout=True)

    return gene_ranges

'''
when mongodb is not available!
'''
def get_chrom_genes_with_jq(chrom,json_file):
    cmd = """/share/apps/genomics/jq -c '[.gene_id, .gene_name, .chrom, .start, .stop, .xstart, .xstop] | select(.[2]=="%s")|{gene_id:.[0],gene_name:.[1],chrom:.[2],start:.[3],stop:.[4],xstart:.[5],xstop:.[6]}' """ % chrom
    result = subprocess.check_output(cmd+json_file,shell=True)
    return helper.split_iter(result)

                    
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
        'N',
        'top_pos_N',
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
        'gnomad_path',
        'hpos_to_test',
        'hpo_mask',
        'cadd_step',
        'second_cadd_gap',
        'cis_gap',
        'output',
        }
    helper.check_args(compulsory_keys, kwargs, 'main')
    # defaults
    kwargs.setdefault('gene_inheritance_mode',{})
    # output already exist?
    if os.path.isfile(kwargs['output']):
        print('already done')
        return None
    # get patient_mini and patient_info
    patient_info = helper.get_snapshot(kwargs['patient_info_file'])
    patient_mini = helper.get_snapshot(kwargs['patient_mini_file'])
    # get p_h for all hpos
    phs = helper.get_phs(patient_info)
    # add cohort info into patient_mini
    # cohort is not needed here, but just to keep the structure...
    # all_p = MONGO['patient_db'].patients.find({'external_id':{'$in':patient_mini.keys()}},{'external_id':1,'contact':1})
    all_p = list(patient_mini.keys())
    for i in all_p:
        patient_mini[i] = {'hpo': patient_mini[i],
                                          #'contact': i['contact']['user_id']}
                                          'contact':None}
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
    if kwargs.get('genes',None) is not None:
        genes = phenopolis_utils.symbols_to_ids(kwargs['genes'],MONGO['phenopolis_db'])
        this = {'gene_id':{'$in':genes}}

    # sometimes the cursor times out.
    # but do remember to close it
    if kwargs.get('chrom',None) is not None:
        gene_ranges = get_chrom_genes(kwargs['chrom'], fields, MONGO['phenopolis_db'])
        #gene_ranges = get_chrom_genes_with_jq(kwargs['chrom'],kwargs['uclex_genes_json'])
    else:
        gene_ranges = MONGO['phenopolis_db'].genes.find(this,fields,no_cursor_timeout=True)
    # get gnomad and cadd steps
    gnomad_steps = np.arange(
            0,
            kwargs['gnomad_cutoff']*2,
            kwargs['gnomad_cutoff']
            )
    cadd_steps = np.arange(0, 70, kwargs['cadd_step'])

    # for each gene, get all valid variants/patients according to p/v_cutoff, 
    # annotate using gnomad
    result = {}
    number_processed = 0
    for gene_range in gene_ranges:
        # print progress
        number_processed += 1
        if not number_processed % 100:
            print('===processed {} genes==='.format(number_processed))
        print('processing {}'.format(gene_range['gene_name']))
        # first parse vcf file to get genotype and coverage for each variant
        vcf_file = kwargs['vcf_file'].format(gene_range['chrom']) 
        args = dict(
                vcf_file = vcf_file,
                chrom = gene_range['chrom'],
                start = gene_range['start'],
                stop = gene_range['stop'],
                unrelated_file = kwargs['unrelated_file'],
                human_fasta_ref = kwargs['human_fasta_ref'],
                v_cutoff = kwargs['v_cutoff'],
                p_cutoff = kwargs['p_cutoff'],
                gnomad_path = kwargs['gnomad_path'],
                gnomad_cutoff = kwargs['gnomad_cutoff'],
                patient_mini = patient_mini,
                )
        vcf_dfs = get_vcf_df(**args)
        if vcf_dfs is None:
            # no variants, continue
            continue
        genotype_df,cover_df,gnomad_freqs = vcf_dfs
        # then get patients_variants, with variants annotated with
        #  gnomad freqs and cadd
        args = dict(
                gnomad_freqs = gnomad_freqs,
                genotype_df = genotype_df,
                )
        patients_variants = get_patients_variants(**args)
        # if no variants left, skip
        if not patients_variants['variants']: continue
        # for each gene, remove batch-specific variants
        args = dict(
                data = patients_variants,
                patient_mini = patient_mini,
                )
        batch_specific_variants = helper.get_batch_artefacts(**args)
        patients_variants = helper.remove_batch_artefacts(
                patients_variants,
                batch_specific_variants,
                patient_mini,
                )
        # add cadd
        args = dict(
                variants = patients_variants['variants'],
                chrom = gene_range['chrom'],
                start = gene_range['start'],
                stop = gene_range['stop'],
                cadd_file = kwargs['cadd_file'],
                )
        cadds = helper.add_cadd(**args)
        for k,v in cadds.items():
            patients_variants['variants'][k]['cadd'] = v
        # get patient_map for recessive and dominant modes
        args = dict(
                data = patients_variants,
                vcf = cover_df,
                gnomad_steps = gnomad_steps,
                cadd_steps = cadd_steps,
                second_cadd_gap = kwargs['second_cadd_gap'],
                cis_gap = kwargs['cis_gap']
                )
        patient_map = {'r':{},'d':{}}
        # first get mode if provided. Note that the keys could be id
        #  or gene name
        modes = kwargs['gene_inheritance_mode'].get(
                gene_range['gene_name'],
                kwargs['gene_inheritance_mode'].get(
                    gene_range['gene_id'],
                    'rd'
                    )
                )
        for mode in modes:
            args['mode'] = mode
            patient_map[mode] = helper.get_patient_map(**args)
        # find relevant HPOs from leftover patients, and calculate \
                #phenogenon with CADD bin of 5.
        # cache result for calculated hpo.
        phenogenon_cache = {'r':{},'d':{}}
        # get phenogenon sums on the first gnomad bin.
        # get all testable hpos
        hpos = [i for i,v in phs.items() 
                if v >= kwargs['N'] 
                and i not in kwargs['hpo_mask']]
        for hpo in hpos:
            # inheritance mode: r and d
            # Note that for each HPO, it only keeps the inheritance mode
            #  with the higher hgf score

            for mode in modes:
                args = dict(
                        data = patients_variants,
                        hpo = hpo,
                        mode = mode,
                        patient_info = patient_info,
                        coverage_df = cover_df,
                        patient_map = patient_map[mode],
                        )

                genon =  helper.phenogenon(**args)

                phenogenon_cache[mode][hpo] = genon
        # get all I sums
        genon_sums = {'r':{},'d':{}}
        for mode, v1 in phenogenon_cache.items():
            for hpo, v2 in v1.items():
                # future implementation of getting phenogenon_I_sum_cutoff
                #  from training data
                genon_sums[mode][hpo] = sum([i[0] for i in v2 if i[0]>0])

        # get mean and std
        all_vals = genon_sums['r'].values() + genon_sums['d'].values()
        cutf = np.mean(all_vals) + 2*np.std(all_vals)
        good_hpos = {'r':[],'d':[]}
        mode_sums = {'r':0,'d':0}
        for mode in ('r','d'):
            '''
            top_N = dict(sorted(
                    genon_sums[mode].items(),
                    key = lambda x: x[1],
                    reverse = True
                    )[:kwargs['top_pos_N']])
            '''
            for k,v in genon_sums[mode].items():
                if v >= cutf:
                    good_hpos[mode].append(k)
                    mode_sums[mode] += v

        chosen_mode = 'd' if mode_sums['d'] >= mode_sums['r'] else 'r'
        #get negative set
        negative_hpos = helper.get_negative_hpos(
                MONGO['hpo_db'],
                good_hpos[chosen_mode] + list(kwargs['hpo_mask']),
                hpos
                )
        # calculate hgf(hpo goodness of fit)
        args = dict(
                phenogenons = phenogenon_cache[chosen_mode],
                positive_hpos = good_hpos[chosen_mode],
                negative_hpos = negative_hpos,
                gamma_k = kwargs['gamma_k'],
                gamma_scale = kwargs['gamma_scale'],
                patient_map = patient_map[chosen_mode],
                )
        hgf = get_hgf(**args)
        # sometimes getting all 0s. skip
        if sum(hgf.values()) == 0: continue
        result[gene_range['gene_id']] = {
                'mode': chosen_mode, 
                'symbol':gene_range['gene_name'],
                'chrom':gene_range['chrom'],
                'xstart':gene_range['xstart'],
                'xstop':gene_range['xstop'],
                'data':hgf,
                }
    # close cursor
    gene_ranges.close()

    for g,v in result.items():
        hn = hpo_name(MONGO['hpo_db'],v['data'].keys())
        print(v['symbol'])
        for k,vv in v['data'].items():
            print(hn[k],vv)
    sys.exit()

    # write everything to output
    with open(kwargs['output'], 'w') as outf:
        json.dump(result, outf)

    
if __name__ == '__main__':
    # in the end some of the args have to go to the config
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--chrom",
                      dest="chrom",
                      help="which chrom to process?")
    parser.add_option("--output",
                      dest="output",
                      help="output file name?")
    (options, args) = parser.parse_args()
    args = dict(
        genes = ('TUBB',),
        chrom = options.chrom,
        output = options.output,
        cadd_file = '/SAN/vyplab/UCLex/mainset_August2017'+
                '/mainset_August2017_VEP/CADD_chr{}.vcf.gz',
        vcf_file = '/SAN/vyplab/UCLex/mainset_August2017'+
                '/mainset_August2017_chr{}_filtered.vcf.gz',
        gnomad_path = '/cluster/project8/vyp/gnomad_data',
        uclex_genes_json = '../tests/data/uclex-genes.json',
        #human_fasta_ref = '/cluster/scratch3/vyp-scratch2/'+
        #    'reference_datasets/human_reference_sequence/human_g1k_v37.fasta',
        # scratch3 is not available for the time being. use my copy
        human_fasta_ref = '/cluster/project8/vyp/JingYu/data/human_g1k_v37.fasta',
        patient_mini_file = '../data/private/hpo/patients_hpo_'+
            'snapshot_2017-May_mini.tsv',
        patient_info_file = '../data/private/hpo/patients_hpo_'+
            'snapshot_2017-May.tsv',
        N = 100,
        # choose how many hpos from top as positive hpos
        top_pos_N = 5,
        # if no phasing is done, use cis_gap as a guideline to check
        #  if two variants are in cis
        cis_gap = 100,
        # these two cutoffs are for getting cleaned vcf, 
        v_cutoff = 0.4,
        p_cutoff = 0.4,
        # this cutoff is to get poorly covered individuals 
        #  for a given set of variants, to get patient_map
        # e.g. we have a~e five variants. if vmc is set at 0.5,
        #  and an individual not covered on a,b,c, it is removed 
        #  from the analysis
        patient_missingness_cutoff = 0.5,
        unrelated_file = '/SAN/vyplab/UCLex/KING/UCL-exome_unrelated.txt',
        gnomad_cutoff = 0.00025,
        # known gene inheritance mode. if provided, no need to infer from data
        #  for sometimes it does make mistakes such that for CERKL
        gene_inheritance_mode = dict(
            ABCA4 = 'r',
            CERKL = 'r',
            SCN1A = 'd',
            GUCY2D = 'd',
            USH2A = 'r',
            PROM1 = 'd',
            TERT = 'd',
            CNGB1 = 'r',
            CRB1 = 'r',
            IMPG2 = 'r',
            RPGR = 'r',
            ),
        # this is to use separate positive set from negative set
        # 6 is from the notebook analysis, assuming hgf >= 15 is good
        # deducing from a training set would be great!
        phenogenon_I_sum_cutoff = 6,
        cadd_step = 5,
        second_cadd_gap = 5,
        min_cadd_for_second_variant = 15,
        # gamma_k and gamma_scale for adjusting LI profiling p values
        # can increase gamma_k for larger cohorts?
        gamma_k = 5,
        gamma_scale = 1,
        # HPOs not wanted to be included in the analysis
        #  usually the ones people don't record, such as 
        #  inheritance mode
        hpo_mask = (
            'HP:0000007',
            'HP:0000006',
            'HP:0003745',
            'HP:0000005'
            ),
        # hpos_to_test not useful. should remove
        hpos_to_test = (
            'HP:0000505',
            'HP:0000510',
            'HP:0000479',
            'HP:0000512',
            'HP:0000518',
            'HP:0000556',
            'HP:0000580',
            'HP:0000662',
            'HP:0000726',
            'HP:0000951',
            'HP:0001098',
            'HP:0001250',
            'HP:0001574',
            'HP:0001871',
            'HP:0002037',
            'HP:0002715',
            'HP:0002721',
            'HP:0004430',
            'HP:0005528',
            'HP:0005561',
            'HP:0007703',
            'HP:0011025',
            'HP:0011675',
            'HP:0012145',
            'HP:0100543',
            ),
    )
    main(**args)
    print('==done==')
