from lib import helper
import json


def main(params):
    '''
    main logic

    Things happen before running this script
    =====
    All clinvar genes that have pathongetic variants causing Mendelian are used to extract variants from UCLex's vcf.

    Variants that 1) do not pass filter, 2) far away from coding exons (distance > exon_padding) are excluded.

    Variants internal AF (excluding patients with overlapping assocated HPO terms), gnomad_af, gnomad_hf and CADD are annotated. The internal AF will be used to simulate control samples

    Each gene's top OMIM, and the top OMIM's co-occuring HPOs are also included. The pathogenic variants causing the OMIM will be used to simulate patients

    All the above data are in params['simulation']['annotated_variants']

    For MOI, the truth is scraped from the OMIM website (params['simulation']['omim_moi'])
    for gene/OMIM causing only somatic, ignore


    '''
    # get all genes to be simulated
    genes = None
    with open(params['simulation']['annotated_variants'], 'rt') as inf:
        genes = json.load(inf)

    # get sample names
    # generic samples start with 'simu_'
    # diseased samples start with its gene symbol
    # pop and cohort noises both on an imagined chromosome Z
    # cohort noise samples have a different cohort id
    samples = helper.simulate_samples(genes, params)
    print(genes[0])
    print(samples)


if __name__ == '__main__':
    config_file = 'configure.cfg'
    config = helper.parse_config(config_file)
    main(config)
    print('===done===')
