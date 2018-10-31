import unittest
import os
import sys
import shutil
import filecmp
import json
import subprocess
sys.path.append('cutoff')
sys.path.append('cutoff/cluster_scripts')
sys.path.append('commons')
import common_utils

# ABCA4: 1:94400000-94600000
# SCN1A: 2:166844000-167007000
class FakeOptions():
    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self,key,kwargs[key])

class PhenogenonTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = 'tests/data/tmp'
        # mkdir -p tests/data/tmp
        common_utils.mkdir_p(self.tmp_folder)
        # expand ABCA4_SCN1A.fasta.gz if not already so
        test_ref_fasta = 'tests/data/ABCA4_SCN1A.fasta'
        if not os.path.isfile(test_ref_fasta):
            subprocess.run([
                'gunzip',
                '-c',
                'tests/data/ABCA4_SCN1A.fasta.gz',
                '>{}'.format(test_ref_fasta),
            ])

        # hardcode all options
        self.input_options = dict(
            human_fasta_ref = 'tests/data/ABCA4_SCN1A.fasta',
            gnomad_path = 'tests/data/gnomad_data',
            patient_mini_file = 'tests/data/test_patients_hpo_snapshot_mini.tsv',
            unrelated_file = 'tests/data/test_unrelated.tsv',
            coding_variant_file = 'tests/data/chr{}.coding.tsv',
            vcf_file = '',
            cadd_file = 'tests/data/CADD_ABCA4_SCN1A.vcf.gz',
            patient_info_file = 'tests/data/test_patients_hpo_snapshot.tsv',
            hpo_json = 'tests/data/new-hpo-hpo.json',
            N = 60,
            # removing non-coding?
            remove_nc = True,
            # use this to help find top HPO terms
            # mean + genon_sum_cutoff_coefficient * standard deviation
            genon_sum_cutoff_coefficient = 1,
            # if no phasing is done, use cis_gap as a guideline to check
            #  if two variants are in cis
            cis_gap = 100,
            # these two cutoffs are for getting cleaned vcf, 
            v_cutoff = 0.4,
            p_cutoff = 0.4,
            # this cutoff is to get poorly covered individuals 
            #  for a given set of variants, to get patient_map
            # e.g. we have a~e five variants. if vmc is set at 0.5,
            #  and an individual not covered on a,b,c, it is removed 
            #  from the analysis
            patient_missingness_cutoff = 0.5,

            # this is not implemented yet but essentially adds more weight
            #  to variants with higher cadd scores for recessive mode.
            cadd_step = 5,
            cadd_min = 0,
            gnomad_step = 0.00025,
            gnomad_cutoff = 0.01,
            # pop_flags: if number of variants in bin >= pop_flags[1]
            #  and propotion of variants with the same highest af pop >= pop_flags[0]
            #  flag the HPO
            pop_check_p = 0.05,
            pop_flags = (0.5, 3),
            # damange_cadd_ind is used to evaluate cadd efficiency
            damage_cadd_ind = 3,
            # this coefficient is used to get positive hpos.
            # the higher the coefficient, the fewer hpos you will get per gene/mode
            coefficient = 1,
            combine_pvalues_method = 'scaled_stouffer',
            stouffer_weights = [0.1, 0.5, 0.75, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8],
            hpo_mask = ['HP:0000007', 'HP:0000006', 'HP:0003745', 'HP:0000005'],

        )

    def tearDown(self):
        # rm -r tests/data/tmp
        shutil.rmtree(self.tmp_folder, ignore_errors = True)

    def test_get_coding_variants(self):
        import get_coding_variants
        outfile = os.path.join(self.tmp_folder,'chr2.coding.tsv')

        fake_options = FakeOptions(
                chrom = '2',
                output = outfile,
                input = 'tests/data/VEP_chr2.csv.gz',
                distance = '5',
        )
        get_coding_variants.main(fake_options)
        self.assertTrue(filecmp.cmp(outfile, 'tests/data/chr2.coding.tsv'))

    def test_patients_variants(self):
        import patients_variants
        
        # ABCA4
        self.input_options.update(dict(
                vcf_file = 'tests/data/ABCA4.anonymised.vcf.gz',
                genes = ('ENSG00000198691',),
                #output = outfile,
        ))
        result = patients_variants.main(**self.input_options)
        with open('tests/data/ABCA4.pv.json', 'rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)
        # SCN1A
        self.input_options.update(dict(
                vcf_file = 'tests/data/SCN1A.anonymised.vcf.gz',
                genes = ('ENSG00000144285',),
                #output = outfile,
        ))
        result = patients_variants.main(**self.input_options)
        with open('tests/data/SCN1A.pv.json', 'rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)
        
    def test_patient_map(self):
        import patient_map
        # ABCA4
        self.input_options.update(dict(
                vcf_file = 'tests/data/ABCA4.anonymised.vcf.gz',
                genes = ('ENSG00000198691',),
        ))
        result = patient_map.main(**self.input_options)
        with open('tests/data/ABCA4.pm.json', 'rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)
        # SCN1A
        self.input_options.update(dict(
                vcf_file = 'tests/data/SCN1A.anonymised.vcf.gz',
                genes = ('ENSG00000144285',),
        ))
        result = patient_map.main(**self.input_options)
        with open('tests/data/SCN1A.pm.json', 'rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)

    def test_phenogenon(self):
        import phenogenon
        # ABCA4
        outfile = os.path.join(self.tmp_folder,'ABCA4.genon.json')
        self.input_options.update(dict(
                vcf_file = 'tests/data/ABCA4.anonymised.vcf.gz',
                patient_map_file = 'tests/data/ABCA4.pm.json',
                genes = ('ENSG00000198691',),
                output = outfile,
        ))
        phenogenon.main(**self.input_options)
        with open(outfile, 'rt') as inf:
            result = json.load(inf)
        with open('tests/data/ABCA4.genon.json','rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)
        # SCN1A
        outfile = os.path.join(self.tmp_folder,'SCN1A.genon.json')
        self.input_options.update(dict(
                vcf_file = 'tests/data/SCN1A.anonymised.vcf.gz',
                patient_map_file = 'tests/data/SCN1A.pm.json',
                genes = ('ENSG00000144285',),
                output = outfile,
        ))
        phenogenon.main(**self.input_options)
        with open(outfile, 'rt') as inf:
            result = json.load(inf)
        with open('tests/data/SCN1A.genon.json','rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)

    def test_hgf(self):
        import goodness_of_fit
        # ABCA4
        outfile = os.path.join(self.tmp_folder,'ABCA4.hgf.json')
        self.input_options.update(dict(
                vcf_file = 'tests/data/ABCA4.anonymised.vcf.gz',
                patient_map_file = 'tests/data/ABCA4.pm.json',
                patients_variants_file = 'tests/data/ABCA4.pv.json',
                phenogenon_file = 'tests/data/ABCA4.genon.json',
                genes = ('ENSG00000198691',),
                output = outfile,
        ))
        result = goodness_of_fit.main(**self.input_options)
        with open('tests/data/ABCA4.hgf.json','rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)

        # SCN1A
        outfile = os.path.join(self.tmp_folder,'SCN1A.hgf.json')
        self.input_options.update(dict(
                vcf_file = 'tests/data/SCN1A.anonymised.vcf.gz',
                patient_map_file = 'tests/data/SCN1A.pm.json',
                patients_variants_file = 'tests/data/SCN1A.pv.json',
                phenogenon_file = 'tests/data/SCN1A.genon.json',
                genes = ('ENSG00000144285',),
                output = outfile,
        ))
        result = goodness_of_fit.main(**self.input_options)
        with open('tests/data/SCN1A.hgf.json','rt') as inf:
            expected = json.load(inf)
        self.assertDictEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
