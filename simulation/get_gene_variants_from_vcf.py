'''
get variants, their internal AF, gnomAD and cadd for simulation
individuals with associated HPOs will be removed before getting internal AF
only consider variants that PASS
remove non-coding that are *far* away from coding regions
'''
from lib import helper
from collections import Counter
import json


def main(params):
    '''
    main logic
    1. get genes from pathogenic clinvar vcf
    2. for each gene
        a. get associated OMIM (top OMIM with moi only, to make MOI prediction easy) and HPOs (the one only co-occur with the top OMIM)
        b. strip off the above HPOs (that are not related to MOI and AOO) from UCL unrelated patients
        c. remove gene pathogenic variants associated with the top OMIM from UCL VCF
        d. add related HPOs randomly (as a background noise) back to the patients, with a predefined change (e.g. 0.05)
        e. add 1 (pair), 2(pairs)... 20(pairs) pathogenic variants randomly to the patients, and insert relevant HPO terms, each pair 50? times
        f. phenogenon
    3. done
    '''
    # read gene ranges
    gene_ranges = helper.parse_gene_ranges(
        params['simulation']['gene_ranges_file'])
    # parse vcf, ignoring variants not covered by gnomad
    genes = helper.parse_vcf(params['generic']['clinvar_vcf'])

    # get omim_moi
    omim_moi = helper.get_omim_moi(params['simulation']['omim_moi'])

    result = []
    # convert omim:hpo to possibility bin
    for gene in genes:
        if gene not in gene_ranges:
            continue
        # get top omim (if competing, pick the first one)
        if not genes[gene]['omim']:
            continue
        print(gene_ranges[gene])
        associated_hpos = None
        nominated_hpo = None
        nominated_omim = None
        nominated_moi = None
        for entry in Counter(genes[gene]['omim']).most_common():
            omim = entry[0]
            if omim not in omim_moi:
                continue
            # get associated hpos
            associated_hpos = genes[gene]['omim:hpo'][omim].keys()
            if not associated_hpos:
                continue
            nominated_hpo = genes[gene]['omim:hpo'][omim].most_common(1)[0][0]
            nominated_omim = omim
            nominated_moi = omim_moi[omim]
            break
        if nominated_omim is None:
            continue
        # convert count to possibility bins
        # for oh in genes[gene]['omim:hpo']:
        #    genes[gene]['omim:hpo'][oh] = helper.counts_to_possibility_bins(genes[gene]['omim:hpo'][oh])
        variants = helper.get_gene_variants_from_vcf(
            gene_ranges[gene], associated_hpos, params)
        # convert each clinvar's variant's hpo and omim from set to list
        for variant in genes[gene]['variants']:
            for key in ('hpo', 'omim'):
                variant[key] = list(variant[key])
        result.append({
            'entrez_id': gene,
            'symbol': genes[gene]['symbol'],
            'gene_range': gene_ranges[gene],
            'omim': nominated_omim,
            'hpo': nominated_hpo,
            'moi': nominated_moi,
            'clinvar': genes[gene],
            'vcf_variants': variants
        })
    with open(params['simulation']['annotated_variants'], 'wt') as outf:
        json.dump(result, outf)


if __name__ == '__main__':
    config_file = 'configure.cfg'
    config = helper.parse_config(config_file)
    main(config)
    print('===done===')
