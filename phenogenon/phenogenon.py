'''
Get a heatmap
'''
from __future__ import print_function, division
import helper
import patient_map as PM


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
        'N',
        'gnomad_path',
        'patient_mini_file',
        'patient_info_file',
        'unrelated_file',
        'gnomad_cutoff',
        'gnomad_step',
        'gnomad_path',
        'hpo_mask',
        'cadd_step',
        'cadd_min',
    }
    helper.check_args(compulsory_keys, kwargs, 'main')
    # defaults
    kwargs.setdefault('gene_inheritance_mode', {})
    # get patient_mini and patient_info
    patient_info = helper.get_snapshot(kwargs['patient_info_file'])
    # get p_h for all hpos
    phs = helper.get_phs(patient_info)

    # get patient_map
    patient_map = PM.main(**kwargs)
    if patient_map is None:
        return None

    modes = 'rd'
    # translate patient_map's key
    pm = {}
    for m in modes:
        pm[m] = {}
        for k, v in patient_map['patient_map'][m].items():
            key = tuple([int(i) for i in k.split(',')])
            pm[m][key] = v
    phenogenon_cache = {'r': {}, 'd': {}}
    # get phenogenon sums on the first gnomad bin.
    # get all testable hpos
    hpos = [i for i, v in phs.items()
            if v >= kwargs['N'] and
            i not in kwargs['hpo_mask']]
    for hpo in hpos:
        # inheritance mode: r and d
        # Note that for each HPO, it only keeps the inheritance mode
        #  with the higher hgf score

        for mode in modes:
            args = dict(
                hpos=hpo,
                mode=mode,
                patient_info=patient_info,
                patient_map=pm[mode],
            )

            genon = helper.phenogenon(**args)

            phenogenon_cache[mode][hpo] = genon.tolist()

    return {
        'phenogenon': phenogenon_cache,
        'NP': patient_map['NP'],
        'patient_map': patient_map['patient_map'],
        'patients_variants': patient_map['patients_variants'],
    }
