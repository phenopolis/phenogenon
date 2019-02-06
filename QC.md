# QC procedure (pseudo code)

```python
# pseudo code

# 1. get patients
patients = get_patients()

# 2. get only unrelated patients according to analysis from KING
unrelated_patients = remove_related_patients(patients, KING_analysis_result)

# 3. given a gene's region, get coding variants or variants that are no more than 5 bases away from closest CDS. All variants have to pass the filter
params = {
  range = '1:12345-67890',
  cds_padding = 5,
  vcf = '1.vcf.gz',
  samples = unrelated_patients,
  filter_pass_only = True,
}
variants = get_coding_variants_from_vcf(**params)

# 4. remove variants that are miss-called in more than 40% of the patients,
# and also remove patients that are miss-called in more than 40% of the variants.
variants = remove_miss_variants(variants)
patients = remove_miss_patients(patients)

# 5. annotated variants using gnomAD and CADD
variants = annotate_variants(variants)

# 6. remove variants that are not covered by gnomAD
variants = remove_uncovered_variants_in_gnomAD(variants)

# 7. remove variants with gnomAD filter of 'SEGDUP'.
# SEGDUP variants usually have high AF but 0 homozygote count.
variants = remove_SEGDUP(variants)

# 8. some variants have very high AF (> 0.01) in one population, but 0 homozygote count.
# remove them.
variants = remove_pop_strange_variants(variants)

# 9. remove batch effect variants
# 9.1 remove variants with 0 gnomAD AF but with more than 2 carriers in the cohort (using contact clinician as a delegate to cohort)
# 9.2 remove variants with statistically too high cohort AF compared with gnomAD AF (p < 1e-4)
variants = remove_batch_effects(variants)

# 10. if two variants co-occur more often ( > 50 % and > 2) in the same patients than expected,
# remove the one with lower CADD when doing analysis assuming recessive MOI.
if MOI = 'recessive':
  variants = remove_cis(variants)

# 11. get phenogenon heatmap for a given HPO
phenogenon = get_phenogenon(variants, patients, hpo)

# 12. Population confoundation. We often see that correlation
#  between a genotype and a HPO can be confounded by ethnicity. 
# This can be a result of a disease cohort predominantly sampled from a specific ethinic group.
# pop_check_p: check phenogenon bins in the first column (from which HGF will be calculated) with p < pop_check_p
# pop_flags: if number of variants in bin >= pop_flags[1]
#  and propotion of variants with the same highest af pop >= pop_flags[0]
#  flag the HPO
pop_check_p=0.05
pop_flags=0.5, 3

HGF = get_hgf(phenogenon, pop_check_p, pop_flags)
```
