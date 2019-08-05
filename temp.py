import json
from phenogenon import helper, goodness_of_fit

config = helper.OFFLINE_CONFIG
config['phenogenon']['N'] = config['phenogenon'].pop('n')
param = config['generic']
param.update(config['debug'])
param.update(config['phenogenon'])

gene_ranges = {
    'ABCA4': '1:94458394-94586689',
    'SCN1A': '2:166845671-166984524',
    'TTN': '2:179390716-179695529',
}

gene = 'TTN'
chrom = gene_ranges[gene].split(':')[0]
param.update(dict(
    vcf_file=f'/home/jing/SAN/vyplab/UCLex/mainset_May2019/VQSR/both/mainset_May2019_chr{chrom}_filtered.vcf.gz',
    range=gene_ranges[gene],
    minimal_output=True,
))
result = goodness_of_fit.main(**param)['result']

with open(f'{gene}.hgf.json', 'wt') as outf:
    json.dump(result, outf, indent=4)
