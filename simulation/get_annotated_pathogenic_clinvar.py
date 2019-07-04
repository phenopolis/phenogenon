'''
get annotated pathogenic clinvar, using gnomad and cadd (needs to be already calculated)
'''
import gzip
import sys
sys.path.append('/home/jing/git/Biotools')
import CommonFuncs


def read_cadd(cadd_file):
    result, header = {}, None
    with open(cadd_file, 'rt') as inf:
        for line in inf:
            if line.startswith('##'):
                continue
            row = line.rstrip().split('\t')
            if header is None:
                header = row
                continue
            rec = dict(zip(header, row))
            v_id = '-'.join(row[:4])
            if v_id in result and rec['GeneID'] not in result[v_id]['gene_ids']:
                result[v_id]['gene_ids'].append(rec['GeneID'])
                result[v_id]['gene_names'].append(rec['GeneName'])
            else:
                result[v_id] = {
                    'cadd': float(row[-1]),
                    'gene_ids': [rec['GeneID']],
                    'gene_names': [rec['GeneName']],
                }
    return result


def main(params):
    # read cadd, store gene names
    cadd = read_cadd(params['cadd'])
    print(cadd)


if __name__ == '__main__':
    params = dict(
        input_vcf='data/pathogenic_clinvar.vcf.gz',
        cadd='data/pathogenic_clinvar.cadd.tsv',
        pop_freqs=dict(
            gnomad_path='/mnt/e/db/gnomAD/release-170228'
        ),
        human_ref='/mnt/e/db/human_g1k_v37.fasta',
        output_vcf='data/pathogentic_clinvar.annotated.vcf.gz'
    )
    main(params)
    print('==done==')
