import os
import json

base_path = '/Users/logust/UCLex/cluster/project8/vyp/JingYu/git/phenopolis_analysis/data'
gene_id = 'ENSG00000155657'
gene_name = 'TTN'
chrom = '2'
release = 'August2017'
hpo = 'HP:0004328'
moi = 'd'
phenogenon_json = os.path.join(
    base_path, 'public', 'cutoff', f'phenogenon_{release}', f'{gene_name}.json')
with open(phenogenon_json, 'rt') as inf:
    genon = json.load(inf)
print([i[0] for i in genon['phenogenon'][moi][hpo]])

ind = 4
patient_map_json = os.path.join(f'pm_{gene_name}.json')
with open(patient_map_json, 'rt') as inf:
    pm = json.load(inf)

patient_hpo_tsv = os.path.join(
    base_path, 'private', 'hpo', 'patients_hpo_snapshot_2018-Oct.tsv')
ph = {}
with open(patient_hpo_tsv, 'rt') as inf:
    for line in inf:
        if line.startswith('#'):
            continue
        row = line.rstrip().split('\t')
        if row[1] == '0':
            continue
        ph[row[0]] = row[2].split(',')
p_a = len(ph)
print(f'p_a: {p_a}')
p_h = len([p for p in ph if hpo in ph[p]])
print(f'p_h: {p_h}')
patients = pm['patient_map'][moi][f'{ind},0'][0]
p_g = len(patients)
print(f'p_g: {p_g}')
gh = [i for i in patients if hpo in ph[i]]
p_gh = len(gh)
print(f'p_gh: {p_gh}')

cadd_range = (ind * 5, (ind + 1) * 5)
patients_variants_json = os.path.join(f'pv_{gene_name}.json')
with open(patients_variants_json, 'rt') as inf:
    pv = json.load(inf)

from collections import Counter
variants = Counter()
tp = 'gnomad_af' if moi == 'd' else 'gnomad_hom_f'
for k, vs in pv['patients'].items():
    for v in vs:
        if cadd_range[0] <= pv['variants'][v]['cadd'] < cadd_range[1] and pv['variants'][v][tp] < 0.00025:
            variants.update([v])
print(variants)
print(len(variants))
import fisher
pval = fisher.pvalue(
    p_a - p_h - p_g + p_gh,
    p_h - p_gh,
    p_g - p_gh,
    p_gh
).right_tail
print(pval)

vcf = f'/SAN/vyplab/UCLex/mainset_May2019/VQSR/both/mainset_May2019_chr{chrom}_filtered.vcf.gz'
# header
import gzip
import pysam
header = []
p_indices = []
with gzip.open(vcf, 'rt') as inf:
    for line in inf:
        if line.startswith('##'):
            continue
        header = line.rstrip().split('\t')
        break
p_indices = [ind for ind, p in enumerate(header) if p in ph]
print(p_indices)
