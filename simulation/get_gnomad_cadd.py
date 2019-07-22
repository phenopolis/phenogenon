'''
separate gnoamd and cadd from annotated_variants for later simulation calculation
'''
import os
import sys
import json
from lib import helper
sys.path.append('../phenogenon')
import gnomad_utils


def main(params):
    variants = {}
    with open(params['simulation']['annotated_variants']) as inf:
        data = json.load(inf)
    for gene in data:
        vs = set([i['variant_id'] for i in gene['clinvar']['variants']]) | set(
            gene['vcf_variants'].keys())
        gnomad = gnomad_utils.overall_freqs(
            list(vs), params['generic']['gnomad_path'])
        for v in gene['clinvar']['variants']:
            variants[v['variant_id']] = {
                'v_id': v['variant_id'],
                'gnomad': gnomad[v['variant_id']],
                'cadd': v['cadd']
            }
        for v_id, v in gene['vcf_variants'].items():
            variants[v_id] = {
                'v_id': v_id,
                'gnomad': gnomad[v_id],
                'cadd': v['cadd']
            }
    with open(params['generic']['gnomad_cadd_data'], 'wt') as outf:
        json.dump(variants, outf)


if __name__ == '__main__':
    config_file = 'configure.cfg'
    config = helper.parse_config(config_file)
    main(config)
    print('===done===')
