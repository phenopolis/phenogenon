from __future__ import print_function, division
import sys
from optparse import OptionParser
import warnings
import json
import os
import math
import numpy as np
from collections import Counter
from scipy import stats
from phenogenon import gnomad_utils, helper
from phenogenon import phenogenon as Pheno


class Goodness_of_fit:
    def __init__(self, genons):
        self.genons = genons
        # convert genons to np array
        for mode in ('r', 'd'):
            for hpo in self.genons[mode]:
                if not isinstance(self.genons[mode][hpo], np.ndarray):
                    self.genons[mode][hpo] = np.array(self.genons[mode][hpo])

    def get_genon_sum(self, mode, hpo):
        '''
        get HGF
        '''
        genon = self.genons[mode][hpo].copy()
        if self.combine_pvalues_method in ('stouffer', 'scaled_stouffer'):
            # convert all genon > 0.5 to 0.5 to stablise stouffer's
            # performance
            # not needed if using fisher

            # disable warning since there is NaN comparison
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                genon[genon > 0.5] = 0.5
        # use stouffer or fisher method to combine p values,
        # with weights for different cadd
        # weights only applicable to stouffer
        weights, pvals = [], []
        for ind, val in enumerate(genon[:, 0]):
            if not np.isnan(val):
                weights.append(
                    self.stouffer_weights[ind]
                )
                pvals.append(val)
        genon_sum = None
        if pvals:
            combine_test = combine_pvalues(
                pvalues=pvals,
                method=self.combine_pvalues_method,
                weights=weights
            )
            genon_sum = -math.log(combine_test[1] or sys.float_info.epsilon)
        return genon_sum

    def get_genon_sratio(self, mode, hpo):
        '''
        return genon_sratio. The lower, the more 'noise' comes from non-rare part of the heatmap
        '''
        genon = self.genons[mode][hpo]
        def log_transform(x): return x if np.isnan(x) else -math.log(x)
        log_transform = np.vectorize(log_transform)
        rare = genon[:, 0][~np.isnan(genon[:, 0])]
        log_rare = sum(log_transform(rare))
        rest = genon[:, 1:][~np.isnan(genon[:, 1:])]
        if rest.size > 0:
            log_rest = sum(log_transform(rest))
            if log_rare + log_rest == 0:
                genon_sratio = 0
            else:
                genon_sratio = log_rare / (log_rare + log_rest)
        else:
            genon_sratio = 1.
        return genon_sratio

    def get_genon_hratio(self, mode, hpo):
        genon_sum = self.hgf[mode][hpo]
        if not genon_sum:
            return None
        genon = self.genons[mode][hpo].copy()
        if self.combine_pvalues_method in ('stouffer', 'scaled_stouffer'):
            # convert all genon > 0.5 to 0.5 to stablise stouffer's
            # performance
            # not needed if using fisher

            # disable warning since there is NaN comparison
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                genon[genon > 0.5] = 0.5
        weights, pvals = [], []
        for ind, i in enumerate(genon[:, 1:]):
            for j in i:
                if not np.isnan(j):
                    weights.append(
                        self.stouffer_weights[ind]
                    )
                    pvals.append(j)
        if len(pvals):
            stouffer = combine_pvalues(
                pvalues=pvals,
                method=self.combine_pvalues_method,
                weights=weights
            )
            S = -math.log(stouffer[1])
            return genon_sum / (genon_sum + S)
        return 1.

    def get_cadd_15_ratio(self, mode, hpo):
        genon = self.genons[mode][hpo].copy()
        damage_arr = genon[self.damage_cadd_ind:, 0]
        weights, pvals = [], []
        for ind, val in enumerate(damage_arr):
            if not np.isnan(val):
                weights.append(
                    self.stouffer_weights[self.damage_cadd_ind + ind]
                )
                pvals.append(val)
        if len(pvals) == 0:
            return 0.
        S = combine_pvalues(
            pvalues=pvals,
            method='fisher',
            weights=weights
        )
        damage_sum = -math.log(S[1] or sys.float_info.epsilon)
        non_damage_arr = genon[:self.damage_cadd_ind, 0]
        # note weights are useless here. will remove
        weights, pvals = [], []
        for ind, val in enumerate(non_damage_arr):
            if not np.isnan(val):
                weights.append(
                    self.stouffer_weights[ind]
                )
                pvals.append(val)
        if len(pvals) == 0:
            return 1.
        S = combine_pvalues(
            pvalues=pvals,
            method='fisher',
            weights=weights
        )
        non_damage_sum = -math.log(S[1])
        if damage_sum:
            cadd_15_ratio = damage_sum / (damage_sum + non_damage_sum)
        else:
            cadd_15_ratio = 0

        return cadd_15_ratio

    def get_pop_alert(self, mode, hpo):
        '''
        get pop cursed?
        return the cursed pop, or None
        '''
        genon = self.genons[mode][hpo]
        # get inds with small p
        s_p_inds = np.where(genon[:, 0] <= self.pop_check_p)[0]
        # get patients, then variants
        variants = {'pos': [], 'neg': []}
        tp = None
        if mode == 'r':
            tp = 'gnomad_hom_f'
        elif mode == 'd':
            tp = 'gnomad_af'
        else:
            msg = 'mode has to be either r or d'
            raise ValueError(msg)
        for ind in s_p_inds:
            patients = self.patient_map[mode]["{},0".format(ind)][0]
            cadd_cuts = (self.cadd_step * ind, self.cadd_step * (ind + 1))
            gnomad_cut = self.gnomad_step
            for p in patients:
                if hpo in self.patient_info[p]['hpo']:
                    curse = 'pos'
                else:
                    curse = 'neg'
                for v in self.patients_variants['patients'][p]:
                    A = (
                        self.patients_variants['variants'][v][tp] < gnomad_cut)
                    B = (cadd_cuts[0] <=
                         self.patients_variants['variants'][v]['cadd'] <
                         cadd_cuts[1])
                    if A and B:
                        variants[curse].append(v)
        pop_curse = {'pos': set(), 'neg': set()}
        if len(variants['pos']) < self.pop_flags[1]:
            # number of variants are too few
            return None
        # annotate variants using gnomad_utils, and find pop curse
        # if pos and neg find same most freq pop, return None
        gnomad_freqs = gnomad_utils.overall_freqs(
            variants['pos'] + variants['neg'],
            self.gnomad_path
        )
        for k, v in variants.items():
            C = Counter()
            for vv in v:
                C.update(gnomad_freqs[vv]['most_freq_pops'])
            # what if there is a tie?!?!
            if len(C) == 0:
                pop_curse[k] = set()
                continue
            most_freq = ([C.most_common(1)[0][0]], C.most_common(1)[0][1])
            for kk, vv in C.items():
                if vv == most_freq[1]:
                    most_freq[0].append(kk)
            if most_freq[1] / len(v) >= self.pop_flags[0]:
                pop_curse[k] = set(most_freq[0])
        return list(pop_curse['pos'] - pop_curse['neg']) or None

    @property
    def MOI_score(self):
        '''
        predict inheritance mode
        return a number.
        Negative means dominant
        Positive means recessive
        Note that sometimes number of patients in 'r' mode is 0
        In this case it still returns 'd'
        '''
        if getattr(self, '_MOI_score', None) is None:
            vals = {}
            hpos = set(self.positive_hpos['r'] + self.positive_hpos['d'])
            for hpo in hpos:
                vals[hpo] = self.hgf['r'].get(hpo, 0) * \
                    self.genon_sratios['r'].get(hpo, 0) - \
                    self.hgf['d'].get(hpo, 0) * \
                    self.genon_sratios['d'].get(hpo, 0)
            self._MOI_score = vals
        return self._MOI_score

    @property
    def positive_hpos(self):
        '''
        Get positive hpo sets.
        '''
        if getattr(self, '_positive_hpos', None) is None:
            ps = {'r': [], 'd': []}

            # get positive hpos and negative hpos
            for mode in ('r', 'd'):
                not_nan = [i for i in self.hgf[mode].values()
                           if i is not None]
                if len(not_nan):
                    cutf = np.mean(not_nan) + \
                        self.coefficient * \
                        np.std(not_nan)
                    for k, v in self.hgf[mode].items():
                        if v is not None and v > cutf:
                            ps[mode].append(k)
                    ps[mode] = helper.hpo_minimum_set(self.hpo_db, ps[mode])
                else:
                    ps[mode] = []
            self._positive_hpos = ps
        return self._positive_hpos

    @property
    def hgf(self):
        if getattr(self, '_genon_sum', None) is None:
            hgf = {'r': {}, 'd': {}}
            for mode in ('r', 'd'):
                for hpo in self.genons[mode]:
                    # is it pop cursed?
                    pop = self.pop_alert[mode].get(hpo, None)
                    if not pop or set(pop) == set(['NFE']):
                        hgf[mode][hpo] = self.get_genon_sum(mode, hpo)
            self._hgf = hgf
        return self._hgf

    @property
    def genon_hratios(self):
        if getattr(self, '_genon_hratios', None) is None:
            genon_hratios = {'r': {}, 'd': {}}
            for mode in ('r', 'd'):
                for hpo in self.positive_hpos[mode]:
                    genon_hratios[mode][hpo] = self.get_genon_hratio(mode, hpo)
            self._genon_hratios = genon_hratios
        return self._genon_hratios

    @property
    def genon_sratios(self):
        if getattr(self, '_genon_sratios', None) is None:
            genon_sratios = {'r': {}, 'd': {}}
            for mode in ('r', 'd'):
                for hpo in self.hgf[mode]:
                    genon_sratios[mode][hpo] = self.get_genon_sratio(mode, hpo)
            self._genon_sratios = genon_sratios
        return self._genon_sratios

    @property
    def cadd_15_ratios(self):
        if getattr(self, '_cadd_15_ratios', None) is None:
            cadd_15_ratios = {'r': {}, 'd': {}}
            for mode in ('r', 'd'):
                for hpo in self.positive_hpos[mode]:
                    cadd_15_ratios[mode][hpo] = self.get_cadd_15_ratio(
                        mode, hpo)
            self._cadd_15_ratios = cadd_15_ratios
        return self._cadd_15_ratios

    @property
    def pop_alert(self):
        if getattr(self, '_pop_alert', None) is None:
            pop_alert = {'r': {}, 'd': {}}
            for mode in ('r', 'd'):
                for hpo in self.genons[mode]:
                    this = self.get_pop_alert(mode, hpo)
                    if this is not None:
                        pop_alert[mode][hpo] = this
            self._pop_alert = pop_alert
        return self._pop_alert


def get_hpo_from_json(f):
    '''
    if remote server is somehow unavailable, use a local json file instead
    '''
    with open(f, 'r') as inf:
        data = '[' + inf.read().rstrip().replace('\n', ',') + ']'
        data = json.loads(data)
    # convert it to a dict
    return {i['id'][0]: i for i in data}


def get_hgf(**kwargs):
    result = {
        'NP': kwargs['data']['NP'],
    }
    P = Goodness_of_fit(kwargs['data']['phenogenon'])
    # set some parameters for hgf
    for k in (
            'damage_cadd_ind',
            'combine_pvalues_method',
            'stouffer_weights',
            'coefficient',
            'hpo_db',
            'pop_check_p',
            'pop_flags',
            'patient_map',
            'patients_variants',
            'gnomad_step',
            'cadd_step',
            'gnomad_path',
            'patient_info',
    ):
        setattr(P, k, kwargs[k])

    # get result
    for k in (
        'pop_alert',
        # 'genon_hratios',
        'cadd_15_ratios',
        # 'genon_sratios',
        'MOI_score',
    ):
        result[k] = getattr(P, k)

    # get hgf only for positive HPOs
    result['hgf'] = {}
    for mode in ('r', 'd'):
        result['hgf'][mode] = {k: v
                               for k, v in P.hgf[mode].items()
                               if k in P.positive_hpos[mode]
                               }
    # translate HPO ids and MOI?
    friendly_result = None
    if not kwargs['minimal_output']:
        trans = {
            k: v['name'][0] for k, v in P.hpo_db.items()
        }
        trans.update({'r': 'recessive', 'd': 'dominant'})
        friendly_result = helper.max_output(result, P.hpo_db, trans)
        friendly_result['MOI'] = {}
        for hpo in friendly_result['MOI_score']:
            if friendly_result['MOI_score'][hpo] > 0:
                friendly_result['MOI'][hpo] = 'recessive'
            elif friendly_result['MOI_score'][hpo] < 0:
                friendly_result['MOI'][hpo] = 'dominant'
            else:
                friendly_result['MOI'] = None
        friendly_result['number_of_patients'] = friendly_result.pop('NP')
    # return result
    return {'result': result, 'friendly_result': friendly_result}


def combine_pvalues(pvalues, method='fisher', weights=None):
    '''
    a copy of scipy.stats method,
    but added a stouffer method with customised scale
    '''
    pvalues = np.asarray(pvalues)
    if pvalues.ndim != 1:
        raise ValueError("pvalues is not 1-D")

    if method == 'fisher':
        Xsq = -2 * np.sum(np.log(pvalues))
        pval = stats.distributions.chi2.sf(Xsq, 2 * len(pvalues))
        return (Xsq, pval)
    elif method in ('stouffer', 'scaled_stouffer'):
        if weights is None:
            weights = np.ones_like(pvalues)
        elif len(weights) != len(pvalues):
            raise ValueError("pvalues and weights must be of the same size.")

        weights = np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("weights is not 1-D")

        Zi = stats.distributions.norm.isf(pvalues)
        if method == 'stouffer':
            Z = np.dot(weights, Zi) / np.linalg.norm(weights)
        else:
            Z = np.dot(weights, Zi) / math.sqrt(len(weights))
        pval = stats.distributions.norm.sf(Z)

        return (Z, pval)
    else:
        raise ValueError(
            "Invalid method '%s'. Options are 'fisher', 'stouffer' or 'scaled_stouffer", method)


def draw_phenogenon(**kwargs):
    '''
    output pdf to visualise phenogenon heatmaps for each hpo/moi
    '''
    import itertools
    import plotly.graph_objs as go
    import plotly.io as pio

    # log transform p values
    def log_transform(x): return x if np.isnan(x) else -math.log(x)
    log_transform = np.vectorize(log_transform)

    for moi in ('r', 'd'):
        # make folder
        outdir = os.path.join(kwargs['heatmap_outdir'], moi)
        helper.mkdir_p(outdir)
        for hpo in kwargs['hpos'][moi]:
            data = kwargs['data']['phenogenon'][moi][hpo]
            # get file name:
            outfile = os.path.join(
                outdir, '{}.pdf'.format(hpo.replace(':', '_')))
            # get zmax
            zmax = min(kwargs['zmax'], max(
                list(itertools.chain.from_iterable([i for i in data]))))
            # get title
            title = '{}-{}'.format(hpo, 'recessive' if moi
                                   == 'r' else 'dominant')
            # get trace
            trace = {
                'z': log_transform(data),
                'x': np.arange(0, kwargs['gnomad_cutoff'], kwargs['gnomad_step']),
                'y': np.arange(0, 60, kwargs['cadd_step']),
                'connectgaps': False,
                'type': 'heatmap',
                'showscale': True,
                'zmin': 0,
                'zauto': False
            }
            layout = {
                'title': title,
                'xaxis': {'title': 'GF'},
                'yaxis': {'title': 'CADD_phred'}
            }
            fig = go.Figure(data=[trace], layout=layout)
            pio.write_image(fig, outfile)


def main(**kwargs):
    kwargs['hpo_db'] = get_hpo_from_json(kwargs['hpo_json'])
    # get patient_mini and patient_info
    kwargs['patient_info'] = helper.get_snapshot(kwargs['patient_info_file'])

    # if there are symbols, turn them into ensembl ids
    pheno = Pheno.main(**kwargs)
    kwargs['data'] = pheno
    kwargs['patient_map'] = pheno.pop('patient_map')
    kwargs['patients_variants'] = pheno.pop('patients_variants')
    # get hgf
    hgf = get_hgf(**kwargs)
    # produce heatmaps
    if kwargs.get('heatmap_outdir', None) is not None:
        kwargs['hpos'] = {'r': hgf['result']['hgf']
                          ['r'].keys(), 'd': hgf['result']['hgf']['d'].keys()}
        draw_phenogenon(**kwargs)
    return hgf


if __name__ == '__main__':
    # in the end some of the args have to go to the config
    USAGE = "usage: %prog [options] arg1 arg2"
    PARSER = OptionParser(usage=USAGE)

    PARSER.add_option("--output",
                      dest="output",
                      help="output file name?")

    PARSER.add_option("--vcf_file",
                      dest="vcf_file",
                      help="bgzipped and tabix-indexed vcf.gz")

    PARSER.add_option("--range",
                      dest="range",
                      help="genome range to calculate? e.g. 2:4000-6000")

    (OPTIONS, _) = PARSER.parse_args()
    if os.path.isfile(OPTIONS.output):
        print('already done')
        sys.exit()
    ARGS = dict(
        range=OPTIONS.range,
        vcf_file=OPTIONS.vcf_file,
    )
    # update args with commons.cfg
    ARGS.update(helper.OFFLINE_CONFIG['generic'])
    ARGS.update(helper.OFFLINE_CONFIG['phenogenon'])
    ARGS['N'] = ARGS.pop('n')

    # get result and write to output
    RESULT = main(**ARGS)
    with open(OPTIONS.output, 'wt') as outf:
        json.dump(RESULT, outf)
    print('==done==')
