'''
hpo goodness of fit test.
For each gene, test each HPO with P_h >= N (N default 100),
  using the rest HPO terms (not daughter or ancestor) with P_h >= N
  as negative set
'''
from __future__ import print_function, division
import itertools
import numpy as np
import patients_variants as PV
import helper


def main(**kwargs):
    '''
    parameters:
     genes: optional
     N (selecting HPO with at least N Ph. affecting both \
          #positive (selecting parental HPO in the positive set \
          #and negative set)
     vcf file location
     gnomad files location
     patient_mini, patient_info, both are json files
     cadd path
     unrelated file used to subset vcf file
     v cutoff and p cutoff are to remove variants and patients with \
          #low coverage over the gene

    returns hpo goodness of fit score, p_g (gnomad_freq
    '''
    # check args
    compulsory_keys = {
        'remove_nc',
        'vcf_file',
        'gnomad_path',
        'cadd_file',
        'patient_mini_file',
        'patient_info_file',
        'human_fasta_ref',
        'unrelated_file',
        'v_cutoff',
        'p_cutoff',  # not using this since no phasing is done
        'gnomad_cutoff',
        'gnomad_step',
        'gnomad_path',
        'cadd_step',
        'cadd_min',
        'cis_gap',
    }
    helper.check_args(compulsory_keys, kwargs, 'main')
    # defaults
    kwargs.setdefault('gene_inheritance_mode', {})
    # get patient_mini and patient_info
    gnomad_steps = np.arange(
        0,
        kwargs['gnomad_cutoff'] + kwargs['gnomad_step'],
        kwargs['gnomad_step']
    )
    cadd_steps = np.arange(kwargs['cadd_min'], 60, kwargs['cadd_step'])

    # for each gene, get all valid variants/patients according to p/v_cutoff,
    # annotate using gnomad
    # then get patients_variants, with variants annotated with
    #  gnomad freqs and cadd
    patients_variants = PV.main(**kwargs)

    # get patient_map for recessive and dominant modes
    args = dict(
        data=patients_variants,
        vcf=patients_variants['cover_df'],
        gnomad_steps=gnomad_steps,
        cadd_steps=cadd_steps,
        cis_gap=kwargs['cis_gap'],
    )
    patient_map = {'r': {}, 'd': {}}
    # first get mode if provided. Note that the keys could be id
    #  or gene name
    modes = 'rd'

    # get number of patients who carry rare variants when get patient_maps
    NP = {}
    for mode in modes:
        args['mode'] = mode
        M = helper.get_patient_map(**args)
        NP[mode] = len(set(list(itertools.chain.from_iterable(
            [v[0] for k, v in M.items() if k[1] == 0]
        ))))
        # change the keys to a string
        for k, v in M.items():
            patient_map[mode]['{},{}'.format(k[0], k[1])] = [
                list(v[0]), list(v[1])]

    return {
        'patient_map': patient_map,
        'patients_variants': patients_variants,
        'NP': NP,
    }
