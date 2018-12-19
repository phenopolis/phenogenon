from __future__ import print_function, division
import sys
from optparse import OptionParser
import warnings
import json
import os
import math
import copy
import numpy as np
from collections import defaultdict,Counter
from scipy import stats
import gnomad_utils
import helper
sys.path.append('lib/commons')
from commons import phenopolis_utils
import phenogenon as Pheno

# MONGO = phenopolis_utils.get_mongo_collections()

class Goodness_of_fit:
    def __init__(self,genons):
        self.genons = genons
        # convert genons to np array
        for mode in ('r','d'):
            for hpo in self.genons[mode]:
                if type(self.genons[mode][hpo]) is not np.ndarray:
                    self.genons[mode][hpo] = np.array(self.genons[mode][hpo])

    def get_genon_sum(self,mode,hpo):
        genon = self.genons[mode][hpo].copy()
        if self.combine_pvalues_method in ('stouffer','scaled_stouffer'):
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
        weights, pvals = [],[]
        for ind,val in enumerate(genon[:,0]):
            if not np.isnan(val):
                weights.append(
                        self.stouffer_weights[ind]
                )
                pvals.append(val)
        genon_sum = None
        if len(pvals):
            combine_test = combine_pvalues(
                    pvalues = pvals,
                    method = self.combine_pvalues_method,
                    weights = weights
            )
            genon_sum =  -math.log(combine_test[1] or sys.float_info.epsilon)
        return genon_sum

    def get_genon_sratio(self,mode,hpo):
        genon = self.genons[mode][hpo]
        log_transform = lambda x:x if np.isnan(x) else -math.log(x)
        log_transform = np.vectorize(log_transform)
        rare = genon[:,0][~np.isnan(genon[:,0])]
        log_rare = sum(log_transform(rare))# / len(rare)
        rest = genon[:,1:][~np.isnan(genon[:,1:])]
        if len(rest):
            log_rest = sum(log_transform(rest))# / len(rest)
            genon_sratio = log_rare / (log_rare+log_rest)
        else:
            genon_sratio = 1.
        return genon_sratio

    def get_genon_hratio(self,mode,hpo):
        genon_sum = self.hgf[mode][hpo]
        if not genon_sum:
            return None
        genon = self.genons[mode][hpo].copy()
        if self.combine_pvalues_method in ('stouffer','scaled_stouffer'):
            # convert all genon > 0.5 to 0.5 to stablise stouffer's
            # performance
            # not needed if using fisher

            # disable warning since there is NaN comparison
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                genon[genon > 0.5] = 0.5
        weights, pvals = [],[]
        for ind,i in enumerate(genon[:,1:]):
            for j in i:
                if not np.isnan(j):
                    weights.append(
                            self.stouffer_weights[ind]
                    )
                    pvals.append(j)
        if len(pvals):
            stouffer = combine_pvalues(
                    pvalues = pvals,
                    method = self.combine_pvalues_method,
                    weights = weights
            )
            S = -math.log(stouffer[1])
            return genon_sum / (genon_sum + S)
        return 1.

    def get_genon_vratio(self,mode,hpo):
        genon = self.genons[mode][hpo].copy()
        damage_arr = genon[self.damage_cadd_ind:,0]
        weights,pvals = [],[]
        for ind,val in enumerate(damage_arr):
            if not np.isnan(val):
                weights.append(
                        self.stouffer_weights[self.damage_cadd_ind + ind]
                )
                pvals.append(val)
        if len(pvals) == 0:
            return 0.
        S = combine_pvalues(
                pvalues = pvals,
                method = 'fisher',
                weights = weights
        )
        damage_sum = -math.log(S[1] or sys.float_info.epsilon)
        non_damage_arr = genon[:self.damage_cadd_ind,0]
        # note weights are useless here. will remove
        weights,pvals = [],[]
        for ind,val in enumerate(non_damage_arr):
            if not np.isnan(val):
                weights.append(
                        #ind * self.stouffer_weight_slope + 0.5
                        self.stouffer_weights[ind]
                )
                pvals.append(val)
        if len(pvals) == 0:
            return 1.
        S = combine_pvalues(
                pvalues = pvals,
                method = 'fisher',
                weights = weights
        )
        non_damage_sum = -math.log(S[1])
        if damage_sum:
            genon_vratio = damage_sum / (damage_sum + non_damage_sum)
        else:
            genon_vratio = 0

        return genon_vratio

    def get_pop_curse_flag(self,mode,hpo):
        '''
        get pop cursed?
        return the cursed pop, or None
        '''
        genon = self.genons[mode][hpo]
        # get inds with small p
        s_p_inds = np.where(genon[:,0] <= self.pop_check_p)[0]
        # get patients, then variants
        variants = {'pos':[],'neg':[]}
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
            cadd_cuts = (self.cadd_step * ind, self.cadd_step * (ind+1))
            gnomad_cut = self.gnomad_step
            for p in patients:
                if hpo in self.patient_info[p]['hpo']:
                    curse = 'pos'
                else:
                    curse = 'neg'
                for v in self.patients_variants['patients'][p]:
                    A = (self.patients_variants['variants'][v][tp] <  gnomad_cut)
                    B = (cadd_cuts[0] <= \
                            self.patients_variants['variants'][v]['cadd'] < \
                            cadd_cuts[1])
                    if A and B:
                        variants[curse].append(v)
        pop_curse = {'pos':set(),'neg':set()}
        if len(variants['pos']) < self.pop_flags[1]:
            # number of variants are too few
            return None
        # annotate variants using gnomad_utils, and find pop curse
        # if pos and neg find same most freq pop, return None
        gnomad_freqs = gnomad_utils.overall_freqs(variants['pos'] + variants['neg'], self.gnomad_path)
        for k,v in variants.items():
            C = Counter()
            for vv in v:
                C.update(gnomad_freqs[vv]['most_freq_pops'])
            # what if there is a tie?!?!
            if len(C) == 0:
                pop_curse[k] = set()
                continue
            most_freq = ([C.most_common(1)[0][0]], C.most_common(1)[0][1])
            for kk,vv in C.items():
                if vv == most_freq[1]:
                    most_freq[0].append(kk)
            if most_freq[1] / len(v) >= self.pop_flags[0]:
                pop_curse[k] = set(most_freq[0])
        return list(pop_curse['pos'] - pop_curse['neg']) or None


    @property
    def predicted_moi(self):
        '''
        predict inheritance mode
        return a number.
        Negative means dominant
        Positive means recessive
        Note that sometimes number of patients in 'r' mode is 0
        In this case it still returns 'd'
        '''
        if getattr(self,'_predicted_moi', None) is None:
            gc = {'r':{},'d':{}}
            vals = {}
            for mode in ('r','d'):
                # only calculate for positive hpos
                for hpo in self.positive_hpos[mode]:
                    gc[mode][hpo] = self.hgf[mode][hpo] * \
                            self.genon_sratios[mode][hpo]
            for mode in ('r','d'):
                vals[mode] = max(gc[mode].values() or [0])
            moi = vals['r'] - vals['d']
            self._predicted_moi = moi
        return self._predicted_moi

    @property
    def positive_hpos(self):
        '''
        Get positive hpo sets.
        '''
        if getattr(self, '_positive_hpos', None) is None:
            ps = {'r':[],'d':[]}

            # get positive hpos and negative hpos
            for mode in ('r','d'):
                not_nan = [i for i in self.hgf[mode].values()
                        if i is not None]
                if len(not_nan):
                    cutf = np.mean(not_nan) + \
                            self.coefficient * \
                            np.std(not_nan)
                    for k,v in self.hgf[mode].items():
                        if v is not None and v > cutf:
                            ps[mode].append(k)
                    ps[mode] = helper.hpo_minimum_set(self.hpo_db,ps[mode])
                else:
                    ps[mode] = []
            self._positive_hpos = ps
        return self._positive_hpos

    @property
    def hgf(self):
        if getattr(self, '_genon_sum', None) is None:
            hgf = {'r':{},'d':{}}
            for mode in ('r','d'):
                for hpo in self.genons[mode]:
                    # is it pop cursed?
                    pop = self.pop_curse_flags[mode].get(hpo, None)
                    if not pop or set(pop) == set(['NFE']):
                        hgf[mode][hpo] = self.get_genon_sum(mode,hpo)
            self._hgf = hgf
        return self._hgf

    @property
    def genon_hratios(self):
        if getattr(self, '_genon_hratios', None) is None:
            genon_hratios = {'r':{},'d':{}}
            for mode in ('r','d'):
                for hpo in self.positive_hpos[mode]:
                    genon_hratios[mode][hpo] = self.get_genon_hratio(mode,hpo)
            self._genon_hratios = genon_hratios
        return self._genon_hratios

    @property
    def genon_sratios(self):
        if getattr(self, '_genon_sratios', None) is None:
            genon_sratios = {'r':{},'d':{}}
            for mode in ('r','d'):
                for hpo in self.positive_hpos[mode]:
                    genon_sratios[mode][hpo] = self.get_genon_sratio(mode,hpo)
            self._genon_sratios = genon_sratios
        return self._genon_sratios

    @property
    def genon_vratios(self):
        if getattr(self, '_genon_vratios', None) is None:
            genon_vratios = {'r':{},'d':{}}
            for mode in ('r','d'):
                for hpo in self.positive_hpos[mode]:
                    genon_vratios[mode][hpo] = self.get_genon_vratio(mode,hpo)
            self._genon_vratios = genon_vratios
        return self._genon_vratios

    @property
    def pop_curse_flags(self):
        if getattr(self, '_pop_curse_flags', None) is None:
            pop_curse_flags = {'r':{},'d':{}}
            for mode in ('r','d'):
                for hpo in self.genons[mode]:
                    this = self.get_pop_curse_flag(mode,hpo)
                    if this is not None:
                        pop_curse_flags[mode][hpo] = this
            self._pop_curse_flags = pop_curse_flags
        return self._pop_curse_flags

def get_hpo_from_json(f):
    '''
    if remote server is somehow unavailable, use a local json file instead
    '''
    with open(f,'r') as inf:
        data = '[' + inf.read().rstrip().replace('\n',',') + ']'
        data = json.loads(data)
    # convert it to a dict
    return {i['id'][0]:i for i in data}

def get_hgf(**kwargs):
    result = {
            'NP':kwargs['data']['NP'],
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
           'pop_curse_flags',
           'genon_hratios',
           'genon_vratios',
           'genon_sratios',
           'predicted_moi',
            ):
       result[k] = getattr(P,k)

    # get hgf only for positive HPOs
    result['hgf'] = {}
    for mode in ('r','d'):
        result['hgf'][mode] = {k:v
                for k,v in P.hgf[mode].items()
                if k in P.positive_hpos[mode]
        }

    # return result
    return result

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

def main(**kwargs):
    kwargs['hpo_db'] = get_hpo_from_json(kwargs['hpo_json'])
    result = {}
    # get patient_mini and patient_info
    kwargs['patient_info'] = helper.get_snapshot(kwargs['patient_info_file'])

    # if there are symbols, turn them into ensembl ids
    pheno = Pheno.main(**kwargs)
    kwargs['data'] = pheno
    kwargs['patient_map'] = pheno.pop('patient_map')
    kwargs['patients_variants'] = pheno.pop('patients_variants')
    # get hgf
    return get_hgf(**kwargs)

if __name__ == '__main__':
    # in the end some of the args have to go to the config
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--output",
                      dest="output",
                      help="output file name?")

    parser.add_option("--vcf_file",
                      dest="vcf_file",
                      help="bgzipped and tabix-indexed vcf.gz")

    parser.add_option("--range",
                      dest="range",
                      help="genome range to calculate? e.g. 2:4000-6000")
    (options, args) = parser.parse_args()
    if os.path.isfile(options.output):
        print('already done')
        sys.exit()
    args = dict(
        range = options.range,
        vcf_file = options.vcf_file,
        # known gene inheritance mode. if provided, no need to infer from data
        #  for sometimes it does make mistakes such that for CERKL
        #gene_inheritance_mode = dict(
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
    )
    # update args with commons.cfg
    args.update(phenopolis_utils.OFFLINE_CONFIG['generic'])
    args.update(phenopolis_utils.OFFLINE_CONFIG['phenogenon'])
    args['N'] = args.pop('n')

    # get result and write to output
    result = main(**args)
    with open(options.output, 'wt') as outf:
        json.dump(result, outf)
    print('==done==')


