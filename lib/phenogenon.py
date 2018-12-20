'''
Get a heatmap
'''
from __future__ import print_function, division
from optparse import OptionParser
import json
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


if __name__ == '__main__':
    # in the end some of the args have to go to the config
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--range",
                      dest="range",
                      help="which genome range to process?")
    parser.add_option("--output",
                      dest="output",
                      help="output file name?")
    (options, args) = parser.parse_args()
    args = dict(
        range=options.range,
        output=options.output,
        # gene_inheritance_mode = dict(
        #    ABCA4 = 'r',
        #    CERKL = 'r',
        #    SCN1A = 'd',
        #    GUCY2D = 'd',
        #    USH2A = 'r',
        #    PROM1 = 'd',
        #    TERT = 'd',
        #    CNGB1 = 'r',
        #    CRB1 = 'r',
        #    IMPG2 = 'r',
        #    RPGR = 'r',
        #    ),
        gene_inheritance_mode={},
    )
    # update args with commons.cfg
    args.update(helper.OFFLINE_CONFIG['generic'])
    args.update(helper.OFFLINE_CONFIG['phenogenon'])
    main(**args)
    print('==done==')
