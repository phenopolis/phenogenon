'''
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
from collections import defaultdict,Counter
import pandas as pd
import numpy as np
import helper
from sklearn.cluster import KMeans
import itertools
import copy

MONGO = phenopolis_utils.get_mongo_collections()

def get_hpo_from_json(f):
    '''
    if remote server is somehow unavailable, use a local json file instead
    '''
    with open(f,'r') as inf:
        data = '[' + inf.read().rstrip().replace('\n',',') + ']'
        data = json.loads(data)
    # convert it to a dict
    return {i['id'][0]:i for i in data}

def hpo_name(hpo_db, ids):
    '''
    this is only for notebook
    '''
    records = hpo_db.hpo.find({'id':{'$in':ids}},{'id':1,'name':1,'_id':0})
    result = {}
    for r in records:
        result[r['id'][0]] = r['name'][0]
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
        'gnomad_path',
        'patient_mini_file',
        'patient_info_file',
        'unrelated_file',
        'gnomad_cutoff',
        'gnomad_step',
        'gnomad_path',
        'hpo_mask',
        'cadd_step',
        'cadd_min',
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
    # get hpodb from json
    hpo_db = get_hpo_from_json(kwargs['hpo_json'])
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
        gene_ranges = MONGO['phenopolis_db'].genes.find(this,fields,no_cursor_timeout=True)
    else:
        gene_ranges = get_chrom_genes(kwargs['chrom'], fields, MONGO['phenopolis_db'])
    # get gnomad and cadd steps
    gnomad_steps = np.arange(
            0,
            kwargs['gnomad_cutoff']+kwargs['gnomad_step'],
            kwargs['gnomad_step']
            )
    cadd_steps = np.arange(kwargs['cadd_min'], 60, kwargs['cadd_step'])

    # get patient_maps
    if 'patient_map_file' in kwargs:
        patient_map_file = kwargs['patient_map_file']
    else:
        patient_map_file = os.path.join(
                kwargs['patient_maps_path'],
                '{}.json'.format(kwargs['chrom'])
            ),
    with open(patient_map_file, 'r') as inf:
        patient_maps = json.load(inf)

    # for each gene, get all valid variants/patients according to p/v_cutoff, 
    # annotate using gnomad
    result = {}
    number_processed = 0
    outf = open(kwargs['output'], 'w')
    outf.write('{')
    for gene_range in gene_ranges:
        # get patient_map
        patient_map = patient_maps.pop(gene_range['gene_id'], None)
        if patient_map is None:
            continue
        # print progress
        number_processed += 1
        if not number_processed % 100:
            print('===processed {} genes==='.format(number_processed))
        print('processing {}'.format(gene_range['gene_name']))

        modes = kwargs['gene_inheritance_mode'].get(
                gene_range['gene_name'],
                kwargs['gene_inheritance_mode'].get(
                    gene_range['gene_id'],
                    'rd'
                    )
                )

        # translate patient_map's key
        pm = {}
        for m in modes:
            pm[m] = {}
            for k,v in patient_map['patient_map'][m].items():
                key = tuple([int(i) for i in k.split(',')])
                pm[m][key] = v
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
                        hpos = hpo,
                        mode = mode,
                        patient_info = patient_info,
                        patient_map = pm[mode],
                        )

                genon =  helper.phenogenon(**args)

                phenogenon_cache[mode][hpo] = genon.tolist()

        output = json.dumps({
            gene_range['gene_id']:{
                'symbol': gene_range['gene_name'],
                'phenogenon': phenogenon_cache,
                'NP': patient_map['NP'],
            }
        })
        # strip off the braces
        output = output[1:-1]
        if number_processed != 1:
            # meaning not the first record. add a comma
            output = ',' + output
        outf.write(output)

    # close cursor
    gene_ranges.close()

    '''
    for g,v in result.items():
        hn = hpo_name(MONGO['hpo_db'],v['data'].keys())
        print(v['symbol'],v['mode'])
        for k,vv in v['data'].items():
            print(hn[k],vv)
    sys.exit()
    '''
    outf.write('}')

    outf.close()
    
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
        chrom = options.chrom,
        output = options.output,
        patient_maps_path = '../data/private/cutoff/patient_maps',
        #gene_inheritance_mode = dict(
        #    ABCA4 = 'r',
        #    CERKL = 'r',
        #    SCN1A = 'd',
        #    GUCY2D = 'd',
        #    USH2A = 'r',
        #    PROM1 = 'd',
        #    TERT = 'd',
        #    CNGB1 = 'r',
        #    CRB1 = 'r',
        #    IMPG2 = 'r',
        #    RPGR = 'r',
        #    ),
        gene_inheritance_mode = {},
    )
    # update args with commons.cfg
    args.update(phenopolis_utils.OFFLINE_CONFIG['generic'])
    args.update(phenopolis_utils.OFFLINE_CONFIG['phenogenon'])
    main(**args)
    print('==done==')
