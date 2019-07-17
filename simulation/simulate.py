from lib import helper

def main(params):
    '''
    main logic
    1. get genes from pathogenic clinvar vcf
    2. for each gene
        a. get associated OMIM (top OMIM only, to make MOI prediction easy) and HPOs (the one only co-occur with the top OMIM)
        b. strip off the above HPOs (that are not related to MOI and AOO) from UCL unrelated patients
        c. remove gene pathogenic variants associated with the top OMIM from UCL VCF
        d. add related HPOs randomly (as a background noise) back to the patients, with a predefined change (e.g. 0.05)
        e. add 1 (pair), 2(pairs)... 20(pairs) pathogenic variants randomly to the patients, and insert relevant HPO terms, each pair 50? times
        f. phenogenon
    3. done
    '''
    # read gene ranges
    gene_ranges = helper.parse_gene_ranges(params['simulation']['gene_ranges_file'])
    # parse vcf, ignoring variants not covered by gnomad
    genes = helper.parse_vcf(params['generic']['clinvar_vcf'])

    # convert omim:hpo to possibility bin
    for g in genes:
        for oh in genes[g]['omim:hpo']:
            genes[g]['omim:hpo'][oh] = helper.counts_to_possibility_bins(genes[g]['omim:hpo'][oh])
    print('We will have {} genes to simulate'.format(len(genes)))

    for gene in genes.values():
        print(gene_ranges[gene['entrez_id']])
        awef
    
    temp = list(genes.values())
    temp.sort(key=lambda x: len(x['variants']), reverse=True)
    temp = temp[0]
    print(temp['symbol'])
    from collections import Counter
    print(Counter(temp['omim:hpo']))
    print(Counter(temp['omim']))
    print(len(temp['variants']))
    print(len([i for i in genes.values() if len(i['variants']) >= 20]))

if __name__ == '__main__':
    config_file = 'configure.cfg'
    config = helper.parse_config(config_file)
    print(config)
    main(config)
    print('===done===')
