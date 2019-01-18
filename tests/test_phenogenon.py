import unittest
import os
import sys
import shutil
import filecmp
import json
import gzip
import numpy as np
from phenogenon import helper

# 'ABCA4':'1:94458394-94586689'
# 'SCN1A':'2:166845671-166984524'


def C(x, y, epsilon):
    if isinstance(x, dict):
        for k in x:
            C(x[k], y[k], epsilon)
    elif isinstance(x, list):
        for i in range(len(x)):
            C(x[i], y[i], epsilon)
    else:
        if x != y and not (np.isnan(x) and np.isnan(y)):
            try:
                if abs(float(x) - float(y)) > epsilon:
                    raise ValueError('{} != {}'.format(x, y))
            except ValueError:
                raise ValueError('{} != {}'.format(x, y))


class PhenogenonTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = 'tests/data/tmp'
        # for testing float
        self.epsilon = 1e-4
        self.genes = {
            'ABCA4': '1:94458394-94586689',
            'SCN1A': '2:166845671-166984524',
        }
        # mkdir -p tests/data/tmp
        helper.mkdir_p(self.tmp_folder)
        # expand ABCA4_SCN1A.fasta.gz if not already so
        test_ref_fasta = 'tests/data/ABCA4_SCN1A.fasta'
        if not os.path.isfile(test_ref_fasta):
            with open(test_ref_fasta, 'w') as outf:
                with gzip.open('tests/data/ABCA4_SCN1A.fasta.gz', 'rt') as inf:
                    for line in inf:
                        outf.write(line)

        # hardcode all options
        self.input_options = dict(
            human_fasta_ref='tests/data/ABCA4_SCN1A.fasta',
            gnomad_path='tests/data/gnomad_data',
            patient_mini_file='tests/data/test_patients_hpo_snapshot_mini.tsv',
            unrelated_file='tests/data/test_unrelated.tsv',
            gtf='tests/data/ABCA4_SCN1A.GRCh37.75.gtf.gz',
            exon_padding=5,
            vcf_file='',
            minimal_output=True,
            cadd_file='tests/data/CADD_ABCA4_SCN1A.vcf.gz',
            patient_info_file='tests/data/test_patients_hpo_snapshot.tsv',
            hpo_json='tests/data/new-hpo-hpo.json',
            N=60,
            # removing non-coding?
            remove_nc=True,
            # use this to help find top HPO terms
            # mean + genon_sum_cutoff_coefficient * standard deviation
            genon_sum_cutoff_coefficient=1,
            # if no phasing is done, use cis_gap as a guideline to check
            #  if two variants are in cis
            cis_gap=100,
            # these two cutoffs are for getting cleaned vcf,
            v_cutoff=0.4,
            p_cutoff=0.4,
            # this cutoff is to get poorly covered individuals
            #  for a given set of variants, to get patient_map
            # e.g. we have a~e five variants. if vmc is set at 0.5,
            #  and an individual not covered on a,b,c, it is removed
            #  from the analysis
            patient_missingness_cutoff=0.5,

            # this is not implemented yet but essentially adds more weight
            #  to variants with higher cadd scores for recessive mode.
            cadd_step=5,
            cadd_min=0,
            gnomad_step=0.00025,
            gnomad_cutoff=0.01,
            # pop_flags: if number of variants in bin >= pop_flags[1]
            #  and propotion of variants with the same highest af pop >= pop_flags[0]
            #  flag the HPO
            pop_check_p=0.05,
            pop_flags=(0.5, 3),
            # damange_cadd_ind is used to evaluate cadd efficiency
            damage_cadd_ind=3,
            # this coefficient is used to get positive hpos.
            # the higher the coefficient, the fewer hpos you will get per gene/mode
            coefficient=1,
            combine_pvalues_method='scaled_stouffer',
            stouffer_weights=[0.1, 0.5, 0.75, 1., 1.1,
                              1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8],
            hpo_mask=['HP:0000007', 'HP:0000006', 'HP:0003745', 'HP:0000005'],
        )

    def tearDown(self):
        # rm -r tests/data/tmp
        #shutil.rmtree(self.tmp_folder, ignore_errors = True)
        pass

    def list2set(self, obj):
        # turn list to set for comparison
        if isinstance(obj, list):
            if len(obj) > 1 and isinstance(obj[0], (list, dict)):
                return sorted([self.list2set(i) for i in obj])
            elif np.nan in obj:
                return [None if np.isnan(i) else i for i in obj]
            else:
                return set(obj)
        elif isinstance(obj, dict):
            return {k: self.list2set(obj[k]) for k in obj}
        else:
            return obj
        return obj

    def test_patients_variants(self):
        pass

    def test_patient_map(self):
        pass

    def test_phenogenon(self):
        # it tests patients_variants and patient_map at the same time
        from phenogenon import phenogenon
        # ABCA4
        gene = 'ABCA4'
        self.input_options.update(dict(
            vcf_file='tests/data/ABCA4.anonymised.vcf.gz',
            range=self.genes[gene],
        ))
        result = phenogenon.main(**self.input_options)
        del result['patients_variants']['cover_df']
        with open('tests/data/ABCA4.genon.json', 'rt') as inf:
            expected = json.load(inf)
        C(result, expected, self.epsilon)
        # SCN1A
        gene = 'SCN1A'
        self.input_options.update(dict(
            vcf_file='tests/data/SCN1A.anonymised.vcf.gz',
            range=self.genes[gene],
        ))
        result = phenogenon.main(**self.input_options)
        del result['patients_variants']['cover_df']
        with open('tests/data/SCN1A.genon.json', 'rt') as inf:
            expected = json.load(inf)
        C(result, expected, self.epsilon)

    def test_hgf(self):
        from phenogenon import goodness_of_fit
        # ABCA4 1:94458394-94586689
        gene = 'ABCA4'
        self.input_options.update(dict(
            vcf_file='tests/data/ABCA4.anonymised.vcf.gz',
            range=self.genes[gene],
            minimal_output=True,
        ))
        result = goodness_of_fit.main(**self.input_options)
        with open('tests/data/ABCA4.hgf.json', 'rt') as inf:
            expected = json.load(inf)
        C(result, expected, self.epsilon)
        self.input_options['minimal_output'] = False
        result = goodness_of_fit.main(**self.input_options)
        with open('tests/data/ABCA4.maxoutput.hgf.json', 'rt') as inf:
            expected = json.load(inf)
        C(result, expected, self.epsilon)

        # SCN1A 2:166845671-166984524
        gene = 'SCN1A'
        self.input_options.update(dict(
            vcf_file='tests/data/SCN1A.anonymised.vcf.gz',
            range=self.genes[gene],
            minimal_output=True,
        ))
        result = goodness_of_fit.main(**self.input_options)
        with open('tests/data/SCN1A.hgf.json', 'rt') as inf:
            expected = json.load(inf)
        C(result, expected, self.epsilon)
        self.input_options['minimal_output'] = False
        result = goodness_of_fit.main(**self.input_options)
        with open('tests/data/SCN1A.maxoutput.hgf.json', 'rt') as inf:
            expected = json.load(inf)
        C(result, expected, self.epsilon)


if __name__ == '__main__':
    unittest.main()
