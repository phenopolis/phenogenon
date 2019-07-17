'''
get gene range for the genes getting from clinvar-pathogenic vcf
'''
from lib import helper

def main(params):
    # clinvar uses old names for some of the genes
    gene_map = {
            'PRUNE1':'PRUNE',
            'PCARE':'C2orf71',
            'CPLANE1':'C5orf42',
            'PRKN': 'PARK2',
            'DNAAF4':'DYX1C1',
            'MMUT':'MUT',
    }
    # some are RNA genes. ignore them
    ignore_list = {'LOC102724058','LOC107303340','LOC110006319','UGT1A'}
    genes = helper.parse_vcf(params['generic']['clinvar_vcf'])
    header = ['entrez_id', 'symbol', 'chrom', 'start', 'stop']
    with open(params['simulation']['gene_ranges_file'], 'wt') as outf:
        outf.write('\t'.join(header) + '\n')
        for gene in genes.values():
            if gene['symbol'] in ignore_list:
                continue
            print(gene['symbol'])
            rang = helper.find_gene_range_in_gtf(gene_map.get(gene['symbol'],gene['symbol']), params['generic']['gtf'])
            row = [gene['entrez_id'], gene['symbol'], rang['chrom'], str(rang['start']), str(rang['stop'])]
            outf.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    config_file = 'configure.cfg'
    config = helper.parse_config(config_file)
    main(config)
    print('===done===')
