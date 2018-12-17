'''
from VEP files (csv), get variants that are in the coding regions,
or close to exons (splicing)
'''
from __future__ import print_function, division
import csv
import gzip
import os
import sys
import re
from optparse import OptionParser
sys.path.append('../../commons')
import common_utils # need clean_variant

def main(options):
    outfile = options.output
    if outfile is None:
        # default output to ../../data/public/vcf
        outfile = '../../data/public/vcf/chr{}.coding.tsv'.format(options.chrom)
    with gzip.open(options.input.format(options.chrom), 'rt') as inf, \
            open(outfile, 'wt') as outf:
        csvreader = csv.reader(inf)
        header = []
        # repetitive variants should be removed!
        all_variants = set()
        for row in csvreader:
            if not header:
                header = row
                continue
            record = dict(zip(header,row))
            if record['BIOTYPE'] != 'protein_coding':
                # not protein_coding? skip
                continue
            variant = common_utils.clean_variant(
                record['#Uploaded_variation'].replace('_','-')
            )
            if variant in all_variants:
                continue
            all_variants.update([variant])
            HGVSc = record['HGVSc']
            gene_id = record['Gene']
            #print(variant,HGVSc,gene_id, record['Codons'])
            chrom, pos, _, _ = variant.split('-')
            outline = '\t'.join([chrom, pos, variant, gene_id]) + '\n'
            if '-' in HGVSc or '+' in HGVSc:
                # see if it is too far away from coding region
                nodigit = ''.join([i for i in HGVSc.split(':')[1] if i.isdigit() or i in ('-','+','_')])
                locs = nodigit.split('_') # indel has _
                good = False
                for loc in locs:
                    loc_split = re.split(r'-|\+', loc)
                    if len(loc_split) == 1 or int(loc_split[1]) <= int(options.distance):
                        good = True
                if good:
                    outf.write(outline)
            else:
                outf.write(outline)

if __name__ == '__main__':
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("--distance",
                      dest="distance",
                      help="maximum distance to exons?",
                      default="5")
    parser.add_option("--chrom",
                      dest="chrom",
                      help="which chrom to process?")
    parser.add_option("--input",
                      dest="input",
                      help="VEP file? e.g. /foo/baa/VEP_chr{}.csv.gz")
    parser.add_option("--output",
                      dest="output",
                      help="output file name?")
    (options, args) = parser.parse_args()
    main(options)
    print('==done==')
