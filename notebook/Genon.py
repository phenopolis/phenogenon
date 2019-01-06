# Core Genon facility
from __future__ import print_function, division
import numpy as np
import math
from scipy.stats import binom
import json
import sys
import os
import pandas as pd
from collections import Counter, defaultdict
import fisher
import scipy.stats as stats # this is for fisher test. In fact, fisher module is much faster than this,
# but `pip install fisher` has not been successful
import pymongo
import json
import itertools
import pickle
import h5py
import tabix
import re
import functools
import utils
import warnings
#import itertools
#insample_cutoff = 3
np.seterr(divide='raise')

class partialmethod(functools.partial):
    '''
    python2's functools does not have partialmethod
    '''
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return functools.partial(self.func, instance,
                       *(self.args or ()), **(self.keywords or {}))

def partialclass(cls, *args, **kwds):
    '''
    this is used to subclass with some useful defaults
    '''
    class Gene(cls):
        __init__ = partialmethod(cls.__init__, *args, **kwds)

    return Gene

class GeneBase:
    '''
    gene to analyse, including ground truth HPO and inheritance mode
    '''
    def __init__(self, symbol, mode, hpos, **kwargs):
        if type(hpos) is not list:
            msg = "Gene's hpos has to be a list"
            raise AttributeError(msg)
        self.symbol = symbol
        self.mode = mode
        self.hpos = hpos
        for k,v in kwargs.items():
            setattr(self, k, v)

        # some file locations
        self.gene_file = '../data/private/cutoff/raw_{}_gene_hpo.json'.format(symbol)
        self.vcf_file = '../data/public/phasing/vcf/{}_cleaned.vcf'.format(symbol)
        self.vep_vcf_file = '../data/public/phasing/vcf/vep/{}_vep.vcf'.format(symbol)
        self.phase_file = '../data/public/phasing/{}_combined_phase.pkl'.format(symbol)

    @property
    def vcf(self):
        '''
        get vcf coverage
        '''
        if getattr(self, '_vcf', None) is None:
            vcf = utils.get_cover_from_vcf(self.vcf_file)
            self._vcf = vcf
        return self._vcf

    @property
    def phase(self):
        '''
        get phase data
        '''
        if getattr(self, '_phase', None) is None:
            with open(self.phase_file, 'rb') as inf:
                phase = pickle.load(inf)
            self._phase = phase
        return self._phase

    @property
    def patient_maps(self):
        '''
        get patient maps for both mode.
        if you want to recalculate it somehow, reset self._patient_maps to None
        '''
        if getattr(self, '_patient_maps', None) is None:
            self._patient_maps = dict(
                    r = self.get_patient_map('r'),
                    d = self.get_patient_map('d')
                    )
        return self._patient_maps

    def get_batch_artefacts(self, gp):
        # using binom, and gnomad_af as p, to produce probability to help identify batch specific artefacts
        # lower_bound is there to remove cohorts where there is just one patient
        # zero_gnomad_c_cutoff allows max internal count when gnomad_af is 0
        dt_d = defaultdict(Counter)
        dt_r = defaultdict(Counter)
        cohorts = Counter()
        for k,v in gp['patients'].items():
            cohorts[ self.patient_mini[k]['contact'] ] += 1
            vc = Counter(v['variants'])
            for i in vc:
                dt_d[ self.patient_mini[k]['contact'] ][i] += 1
                if vc[i] > 1:
                    dt_r[ self.patient_mini[k]['contact'] ][i] += 1
        # remove cohorts with count lower than lower_bound
        for k in cohorts.keys():
            if cohorts[k] < self.lower_bound:
                del cohorts[k]
                del dt_d[k]
                if k in dt_r:
                    del dt_r[k]

        # for heterozygous variants
        result_d = defaultdict(list)
        for k1,v1 in dt_d.items():
            n_variants = len(v1)
            for k2,v2 in v1.items():
                if not gp['variants'][k2]['gnomad_af']:
                    if v2 > self.zero_gnomad_c_cutoff:
                        result_d[k1].append(k2)
                    continue
                prob = 1 - binom.cdf(v2-1, cohorts[k1], gp['variants'][k2]['gnomad_af'])
                if prob < self.binom_cutoff / n_variants:
                    #print(k2,prob)
                    result_d[k1].append(k2)
        for k in result_d:
            result_d[k] = set(result_d[k])

        # for homozygous variants
        result_r = defaultdict(list)
        for k1,v1 in dt_r.items():
            n_variants = len(v1)
            for k2,v2 in v1.items():
                if not gp['variants'][k2]['gnomad_hom_af']:
                    if v2 > self.zero_gnomad_c_cutoff:
                        result_r[k1].append(k2)
                    continue
                prob = 1 - binom.cdf(v2-1, cohorts[k1], gp['variants'][k2]['gnomad_hom_af'])
                if prob < self.binom_cutoff / n_variants:
                    #print(k2,prob)
                    result_r[k1].append(k2)  
        for k in result_r:
            result_r[k] = set(result_r[k])
        return {'d':result_d,'r':result_r}

    # remove batch specific artefact variants
    def remove_batch_artefacts(self, gp, bad_vs, mode='all'):
        result = {
                'patients':{},
                'variants':gp['variants'],
                'gene_id':gp['gene_id'],
                'pat_a':gp['pat_a'],
                }
        bad_p = []
        for k1,v1 in gp['patients'].items():
            cohort = self.patient_mini[k1]['contact']
            this_bad_vs = []
            # collect het artefacts
            if mode != 'r' and cohort in bad_vs['d']:
                this_bad_vs += [i for i in v1['variants'] if i in bad_vs['d'][cohort]]
            # collect hom artefacts
            if mode != 'd' and cohort in bad_vs['r']:
                vc = Counter(v1['variants'])
                for k in vc:
                    if vc[k] > 1 and k in bad_vs['r'][cohort]:
                        this_bad_vs.append(k)
            this_vs = [i for i in v1['variants'] if i not in this_bad_vs]
            if this_vs:
                result['patients'][k1] = {
                        'hpo':v1['hpo'],
                        'variants':this_vs,
                        }
        return result
    def get_genotype_phenotype(self):
        '''
        get the relevant genotype phenotype data
        '''
        data = utils.read_files(self.gene_file)
        # remove nc?
        if self.closeness:
            nc_variants = utils.extract_nc_from_vcf(
                    self.vep_vcf_file,
                    self.closeness
                    )
            utils.remove_noncoding(data,nc_variants)
        utils.cleanse_variants(data)
        # remove batch effect?
        if self.binom_cutoff:
            batch_artefacts = self.get_batch_artefacts(
                    data,
                    )
            data  = self.remove_batch_artefacts(data,batch_artefacts)
        return data

    def get_patient_map(self,mode):
        '''
        get what patients in what bins
        note that this mode may be different from self.mode,
        given a chance to do analysis for the 'wrong' mode
        '''
        gn = np.arange(self.grange[0],self.grange[1],self.steps[1])
        ca = np.arange(self.crange[0],self.crange[1],self.steps[0])
        gp = self.get_genotype_phenotype()
        patient_map = defaultdict(list)
        for i in range(len(ca)):
            p = []
            for j in range(len(gn)):
                p,narrow_vs,not_covered_patients = self.get_patients(
                        gp,
                        mode,
                        (gn[j], gn[j]+self.steps[1]),
                        (ca[i], ca[i]+self.steps[0])
                        )
                patient_map[(i,j)] = (p,not_covered_patients)
        return patient_map

    def get_variants(self,vs,mode,gr,cr):
        '''
        Given mode, gr and cr, return variants
        Narrow variants strictly match gr and cr, 
        and can be used for dominant inheritance.
        For recessive, broad_vs also returned, which
        matches the lower bounds for both gr and cr.
        '''
        mode_dict = {'r':'gnomad_hom_af','d':'gnomad_af'}
        if self.use_af_for_recessive:
            mode_dict['r'] = 'gnomad_af'
        narrow_vs = (k for k,v in vs.items() if gr[0]<=v[mode_dict[mode]]<gr[1] and cr[0]<=self.cadd[k]<cr[1])
        broad_vs = tuple()
        if mode == 'r':
            broad_vs = (k for k,v in vs.items() if v[mode_dict[mode]]<gr[1] and self.cadd[k]>=min(cr[0],self.second_cadd_min))
        return (set(narrow_vs),set(broad_vs))

    def get_patients(self, gp, mode, gr, cr):
        '''
        given ranges, find matching patients
        '''
        p = []
        narrow_vs,broad_vs = self.get_variants(gp['variants'],mode,gr,cr)
        if not narrow_vs:
            return ([],narrow_vs,set())

        # get patients not covered on miss_cutoff* narrow_vs
        s = self.vcf.loc[narrow_vs].sum()
        not_covered_patients = set(s[s<len(narrow_vs)*self.miss_cutoff].index)

        if mode == 'd':
            p = [k for k,v in gp['patients'].items() if set(v['variants']) & narrow_vs]
        elif mode == 'r':
            for k,v in gp['patients'].items():
                good = []
                other = []
                for i in v['variants']:
                    if i in narrow_vs:
                        good.append(i)
                    elif i in broad_vs:
                        other.append(i)
                r_bad = [] # for removing linked variants
                if len(good) and len(good+other) > 1:
                    pos = [int(i.split('-')[1]) for i in good+other]
                    if (len(pos) > len(set(pos))) or (max(pos) - min(pos) > self.gap):
                        if self.phase_cutoff:
                            # any of them is hom?
                            if len(set(good)) < len(good):
                                p.append(k)
                                continue
                            # remove hom in other, as hom in other renders the variants in good unnecessary
                            
                            O = [k1 for k1,v1 in Counter(other).items() if v1 == 1]
                            if len(good+O) < 2: continue
                            for i in itertools.combinations(sorted(good+O,key=lambda x: int(x.split('-')[1])),2):
                                if not set(good) & set(i): continue
                                cis_p = self.phase[k].setdefault(i,0)
                                if cis_p < self.phase_cutoff: 
                                    p.append(k)
                                    break
                                    
                        else:
                            p.append(k)
        return (p,narrow_vs,not_covered_patients)
        


class GenonResult:
    '''
    has HGFs and genon_ratios
    predicted_moi: recessive if positive, dominant if negative.
    The bigger the margin the more confidence in the call
    '''
    def __init__(self):
        self.HGF = defaultdict(lambda: defaultdict(dict))
        self.genon_hratio = defaultdict(lambda: defaultdict(dict))
        self.cadd_15_ratio = defaultdict(lambda: defaultdict(dict))
        # genon_hs = HGF * genon_sratio, for MOI_score calculation
        self.genon_hs = defaultdict(lambda: defaultdict(dict))
        self.genon_sratio = defaultdict(lambda: defaultdict(dict))
        self.ph = defaultdict(lambda: defaultdict(list))
        self.genes = None
        self.predicted_moi = dict()

        # for representing obj, sorting hpos
        self.sort_key = 'HGF'

        # number of patients with 'rare' variants under different modes
        self.np = defaultdict(dict)
    def __str__(self):
        s = ''
        for gene in getattr(self,self.sort_key):
            s += gene
            # write mode
            if self.genes[gene].mode is not None:
                mode = self.genes[gene].mode
                s += ' - Given mode is {}\n'.format(mode)
                s += ' ' * len(gene)
                s += ' - Number of patients with rare variants: {}\n'.format(
                        self.np[gene][mode]
                        )
            else:
                mode_digit = self.predicted_moi[gene]
                if mode_digit > 0:
                    mode = 'r'
                elif mode_digit < 0:
                    mode = 'd'
                else:
                    mode = 'u'
                s += ' - Predicted mode is {}\n'.format(mode)
                s += ' ' * len(gene)
                s += ' - Number of patients for mode d: {}\n'.format(self.np[gene]['d'])
                s += ' ' * len(gene)
                s += ' - Number of patients for mode r: {}\n'.format(self.np[gene]['r'])
            # do not want to carry on with mode:u
            if mode == 'u':
                continue
            s += ' ' * (len(gene) + 1)
            # write genon sums with ratios
            if self.genes[gene].hpos is not None:
                s += '- HPOs are given...'
            else:
                s += '- HPOs are predicted...'
            s += '\n'
            for hpo,Sum in sorted(
                    getattr(self,self.sort_key)[gene][mode].items(), 
                    key = lambda x :x[1],
                    reverse = True):
                # only show positive hpos
                if self.genes[gene].hpos is None and hpo not in self.ph[gene][mode]:
                    continue
                s += '\t{}\n'.format(hpo)
                s += '\t\tHGF: {}\n'.format(Sum)
                s += '\t\tcadd_15_ratio: {}\n'.format(
                        self.cadd_15_ratio[gene][mode][hpo]
                        )
                s += '\t\tMOI_score: {}\n'.format(
                        self.genon_hs[gene]['r'][hpo] - self.genon_hs[gene]['d'][hpo]
                        )
                '''
                s += '\t\tGenon_hratio: {}\n'.format(
                        self.genon_hratio[gene][mode][hpo]
                        )
                s += '\t\tGenon_sratio: {}\n'.format(
                        self.genon_sratio[gene][mode][hpo]
                        )
                '''
        return s
                

class Genon:
    '''
    this is the class to do the major analysis
    '''
    def __init__(self):
        # some file locations
        self.hdf5 = 'patient_map.hdf5'
        self.cadd_file = '../data/public/cutoff/all_cadd.tsv'
        self.patient_mini_file = (
                '../data/private/hpo/patients_hpo_snapshot'
                '_2017-May_mini.tsv'
                )
        self.patient_info_file = (
                '../data/private/hpo/patients_hpo_snapshot'
                '_2017-May.tsv'
                )
        # some parameters
        # do you want to remove non_coding variants?
        # closeness used to remove noncoding variants {closeness} 
        #  away from any exons.
        # if you do not want to remove non_coding variants,
        #  leave it as None
        self.closeness = 4
        # remove batch effect? if yes, set the binom_cutoff. 
        #  if no, leave it as None
        self.binom_cutoff = 1e-10
        # If a cohort's size is lower than {lower_bound},
        # do not process the cohort for batch effect
        self.lower_bound = 2
        # zero_gnomad_c_cutoff allows max internal count when gnomad_af is 0
        self.zero_gnomad_c_cutoff = 2
        # if the matching variants are miss-called by more than
        # {miss_cutoff} of all the patients, discard the variants
        self.miss_cutoff = 0.5
        # phase_cutoff used to phase variants
        self.phase_cutoff = .75
        # for recessive, what is the min for the second variant's cadd
        #  when the first variant's cadd is high(er than {second_cadd_min})?
        self.second_cadd_min = 100
        # how to combine p values. fisher or stouffer.
        # stouffer can add weights
        self.combine_pvalues_method = 'scaled_stouffer'
        # stouffer's weight slope
        self.stouffer_weights = [0.1, 0.5, 0.75, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]
        # only analyse hpo with N >= {N}
        # if no phase data is available, one can set gap to 100,
        # to treat close-by variants as cis
        self.gap = 0
        self.N = 60
        # the following hpos will not be analysed
        # they are inheritance modes that we do not know
        self.hpo_mask = ['HP:0000007','HP:0000006','HP:0003745','HP:0000005','HP:0012823','HP:0003674']
        # steps = (cadd step, gnomad step)
        self.steps = (5, 0.00025)
        # what is the gnomad range for analysis?
        self.grange = (0, 0.01)
        self.crange = (0, 60)
        # if hpos are not given, Genon will try to find out associated hpos
        self.coefficient = 1.
        # what if we use af insteand of hom_f for recessive mode?
        self.use_af_for_recessive = False

        self.damage_cadd = 15

        # mongodb
        self.conn = pymongo.MongoClient()#(host = 'phenotips')
        self.hpo_db = self.conn['hpo']
        self.p_db = self.conn['patients']

    def init_genes(self):
        # genes to analyse
        # subclass Gene with useful defaults
        Gene = partialclass(GeneBase, 
                closeness = self.closeness,
                grange = self.grange,
                crange = self.crange,
                cadd = self.cadd,
                steps = self.steps,
                binom_cutoff = self.binom_cutoff,
                lower_bound = self.lower_bound,
                second_cadd_min = self.second_cadd_min,
                zero_gnomad_c_cutoff = self.zero_gnomad_c_cutoff,
                gap = self.gap,
                miss_cutoff = self.miss_cutoff,
                patient_mini = self.patient_mini,
                phase_cutoff = self.phase_cutoff,
                use_af_for_recessive = self.use_af_for_recessive,
                )
        self.genes = dict(
                ABCA4 = Gene('ABCA4','r',
                    [
                        'HP:0007754',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        'HP:0000510',
                        ]
                    ),
                SCN1A = Gene('SCN1A','d',
                    [
                        'HP:0001250',
                        # HP:0000707 includes dementia HP:0000726, 
                        # which contributes to noise.
                        'HP:0000707',
                        ]
                    ),
                USH2A = Gene('USH2A','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000365',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                RPGR = Gene('RPGR','r',
                    [
                        'HP:0000548',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                GUCY2D = Gene('GUCY2D','d',
                    [
                        'HP:0001103',
                        'HP:0000550',
                        'HP:0000548',
                        'HP:0000639',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                TERT = Gene('TERT','d',
                    [
                        'HP:0005528',
                        'HP:0002754',
                        'HP:0008404',
                        'HP:0000951',
                        'HP:0000234',
                        ]
                    ),
                PROM1 = Gene('PROM1','d',
                    [
                        'HP:0007754',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                CRB1 = Gene('CRB1','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                CNGB1 = Gene('CNGB1','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                CERKL = Gene('CERKL','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                BBS1 = Gene('BBS1','r',
                    [
                        'HP:0000510',
                        'HP:0000518',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        ]
                    ),
                IMPG2 = Gene('IMPG2','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000504',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                )

    # patient info
    @property
    def patient_info(self):
        if getattr(self, '_patient_info', None) is None:
            patient_info = utils.get_snapshot(self.patient_info_file)
            self._patient_info = patient_info
        return self._patient_info

    @property
    def phs(self):
        # get patient numbers for each hpo
        if getattr(self, '_phs', None) is None:
            self._phs = utils.get_phs(self.patient_info)
        return self._phs

    @property
    def patient_mini(self):
        if getattr(self, '_patient_mini', None) is None:
            patient_mini = utils.get_snapshot(self.patient_mini_file)
            # add contact info
            all_p = self.p_db.patients.find(
                    {'external_id': {'$in':patient_mini.keys()}},
                    {'external_id':1,'contact':1}
                    )
            for i in all_p:
                patient_mini[i['external_id']] = {
                        'hpo': patient_mini[i['external_id']],
                        'contact': i['contact']['user_id']
                        }
            self._patient_mini = patient_mini
        return self._patient_mini

    @property
    def cadd(self):
        '''
        get cadd_phred for all variants
        '''
        if getattr(self, '_cadd', None) is None:
            with open(self.cadd_file,'r') as inf:
                cadd = {}
                for row in inf:
                    if row[0] == '#': continue
                    row = row.rstrip().split('\t')
                    v_id = '-'.join(row[:2]+row[2:4])
                    cadd[v_id] = float(row[-1])
            self._cadd = cadd
        return self._cadd
    
    def phenogenon(self,gene,patient_map,hpo):
        # get raw p_a and p_h
        vcf_patients = set(list(gene.vcf))
        raw_p_a = vcf_patients & set(self.patient_info.keys())
        # is hpo a string or list?
        if isinstance(hpo, str):
            raw_p_h = set([kk for kk,vv in self.patient_info.items() if hpo in vv and kk in vcf_patients])
        else:
            for kk,vv in self.patient_info.items():
                if hpo[0] in vv and hpo[1] in vv:
                    #print(kk,vv)
                    pass
            raw_p_h = set([kk for kk,vv in self.patient_info.items() if not (set(hpo) - set(vv)) and kk in vcf_patients])
            print(raw_p_h)
        # make a np matrix
        shape1,shape2 = [],[]
        for k in patient_map:
            shape1.append(k[0])
            shape2.append(k[1])
        logp_df = np.zeros( (
            len(set(shape1)),
            len(set(shape2))
            ) )
        for k,v in patient_map.items():
            p_a = len(raw_p_a - v[1])
            p_h = len(raw_p_h - v[1])
            p_g = len(v[0])
            if not p_g:
                logp_df[k[0]][k[1]] = None
                continue
            p_gh = 0
            for p in v[0]:
                if hpo in self.patient_info[p]:
                    p_gh += 1
            pval = fisher.pvalue(
                    p_a - p_h - p_g + p_gh,
                    p_h - p_gh,
                    p_g - p_gh,
                    p_gh
                    ).right_tail
            logp_df[k[0]][k[1]] = pval# or sys.float_info.epsilon
        return logp_df

    def trim_ns(self,ns,ps):
        '''
        trim ns to make sure there's no negative hpo overlaps positive ones
        '''
        result = []
        for hn in ns:
            if hn in self.hpo_mask:
                continue
            bad = 0
            for hp in ps:
                A = hn in [i['id'][0] for i in utils.get_hpo_ancestors(self.hpo_db, hp)]
                B = hp in [i['id'][0] for i in utils.get_hpo_ancestors(self.hpo_db, hn)]
                if A or B:
                    bad = 1
                    break
            if not bad:
                result.append(hn)
        return result

    def get_positive_negative_hpos(self,HGF,genon_hratio,cadd_15_ratio,genon_sratio,genon_hs,):
        '''
        Get positive and negative hpo sets.
        Default method
        '''
        ps,ns = {'r':[],'d':[]},{'r':[],'d':[]}

        # get positive hpos and negative hpos
        for mode in ('r','d'):
            cutf = np.mean(HGF[mode].values()) + \
                    self.coefficient * \
                    np.std(HGF[mode].values())
            for k,v in HGF[mode].items():
                if v > cutf:
                    ps[mode].append(k)
                elif v < cutf:
                    ns[mode].append(k)
            ps[mode] = utils.hpo_minimum_set(self.hpo_db,ps[mode])
            ns[mode] = utils.hpo_minimum_set(self.hpo_db,ns[mode])
            # ns cant be super/subclass of any ps, and can't be in hpo_mask
            ns[mode] = self.trim_ns(ns[mode],ps[mode])
        return ps,ns

    def predict_moi(self,gs,hr,vr,sr,hs,ph):
        '''
        predict inheritance mode
        return a number. 
        Negative means dominant
        Positive means recessive
        '''
        vals = {'r':0,'d':0}
        for mode in ('r','d'):
            vals[mode] = max([v for k,v in hs[mode].items() if k in ph[mode]] or [0])
        return vals['r'] - vals['d']

    def analyse(self, R, patient_maps = None):
        '''
        Genon analyse
        if hpos for a gene is provided, treat the hpos as truely associated
        else, calculate Phenogenon for all hpos, both modes, and find out
        the associated hpos
        # R is a global cached GenonResult
        '''
        rr = GenonResult()
        rr.genes = self.genes
        hpos = [i for i,v in self.phs.items() if v >= self.N and i not in self.hpo_mask]
        for gene in self.genes:
            this_hpos = self.genes[gene].hpos or hpos
            for mode in ('r','d'):
                # get patient map if it is not given
                if patient_maps is None or gene not in patient_maps:
                    patient_map = self.genes[gene].patient_maps[mode]
                    if patient_maps is not None:
                        patient_maps[gene] = self.genes[gene].patient_maps
                else:
                    patient_map = patient_maps[gene][mode]
                if mode not in R.np[gene]:
                    R.np[gene][mode] = len(set(list(itertools.chain.from_iterable(
                        [v[0] for k,v in patient_map.items() if k[1] == 0]
                        ))))
                rr.np[gene][mode] = R.np[gene][mode]
                for hpo in this_hpos:
                    if hpo not in R.HGF[gene][mode]:
                        genon = self.phenogenon(
                                self.genes[gene],
                                patient_map,
                                hpo
                                )

                        original_genon = genon.copy()
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
                        weights, pvals = [],[]
                        for ind,val in enumerate(genon[:,0]):
                            if not np.isnan(val):
                                weights.append(
                                        #ind * self.stouffer_weight_slope + 0.5
                                        self.stouffer_weights[ind]
                                )
                                pvals.append(val)
                        if len(pvals):
                            combine_test = combine_pvalues(
                                    pvalues = pvals,
                                    method = self.combine_pvalues_method,
                                    weights = weights
                            )
                            HGF =  -math.log(combine_test[1] or sys.float_info.epsilon)
                        else:
                            HGF = 0

                        if not HGF:
                            # HGF is 0. return everything as 0
                            # write to result
                            R.HGF[gene][mode][hpo] = 0
                            R.genon_hratio[gene][mode][hpo] = 0
                            R.genon_sratio[gene][mode][hpo] = 0
                            R.cadd_15_ratio[gene][mode][hpo] = 0
                            R.genon_hs[gene][mode][hpo] = 0

                            for tp in ('HGF','genon_sratio','genon_hratio','cadd_15_ratio','genon_hs',):
                                getattr(rr,tp)[gene][mode][hpo] = 0

                            continue

                        # how about contribution from non_rare variants?
                        weights, pvals = [],[]
                        for ind,i in enumerate(genon[:,1:]):
                            for j in i:
                                if not np.isnan(j):
                                    weights.append(
                                            #ind * self.stouffer_weight_slope + 0.5
                                            self.stouffer_weights[ind]
                                    )
                                    pvals.append(j)
                        if len(pvals):
                            stouffer = combine_pvalues(
                                    pvalues = pvals,
                                    method = self.combine_pvalues_method,
                                    weights = weights
                            )
                            S = -math.log(stouffer[1] or sys.float_info.epsilon)
                        else:
                            S = 0
                        genon_hratio = HGF / (HGF + S)
                        # signal ratio
                        #rare = genon[:,0][~np.isnan(genon[:,0])]
                        #numerator = len(rare[ rare < 1 ])
                        #al = genon[~np.isnan(genon)]
                        log_transform = lambda x:x if np.isnan(x) else -math.log(x)
                        log_transform = np.vectorize(log_transform)
                        rare = original_genon[:,0][~np.isnan(original_genon[:,0])]
                        log_rare = sum(log_transform(rare))# / len(rare)
                        rest = original_genon[:,1:][~np.isnan(original_genon[:,1:])]
                        if log_rare == 0:
                            genon_sratio = 0
                        elif len(rest):
                            log_rest = sum(log_transform(rest))# / len(rest)
                            genon_sratio = log_rare / (log_rare+log_rest)
                        else:
                            genon_sratio = 1.

                        #denominator = len(al[al<1])
                        #genon_sratio = numerator / denominator

                        # what proportion of contribution from damaging variants?
                        # only calculate for rare variants
                        damage_ind = (self.damage_cadd - self.crange[0]) // self.steps[0]
                        damage_arr = original_genon[damage_ind:,0]
                        weights,pvals = [],[]
                        for ind,val in enumerate(damage_arr):
                            if not np.isnan(val):
                                weights.append(
                                        #ind * self.stouffer_weight_slope + 0.5
                                        self.stouffer_weights[damage_ind + ind]
                                )
                                pvals.append(val)
                        stouffer = combine_pvalues(
                                pvalues = pvals,
                                method = 'fisher',#self.combine_pvalues_method,
                                weights = weights
                        )
                        damage_sum = -math.log(stouffer[1])
                        non_damage_arr = original_genon[:damage_ind,0]
                        weights,pvals = [],[]
                        for ind,val in enumerate(non_damage_arr):
                            if not np.isnan(val):
                                weights.append(
                                        #ind * self.stouffer_weight_slope + 0.5
                                        self.stouffer_weights[ind]
                                )
                                pvals.append(val)
                        stouffer = combine_pvalues(
                                pvalues = pvals,
                                method = 'fisher', #self.combine_pvalues_method,
                                weights = weights
                        )
                        non_damage_sum = -math.log(stouffer[1])
                        if damage_sum:
                            cadd_15_ratio = damage_sum / (damage_sum + non_damage_sum)
                        else:
                            cadd_15_ratio = 0

                        # write to result
                        R.HGF[gene][mode][hpo] = HGF
                        R.genon_hratio[gene][mode][hpo] = genon_hratio
                        R.genon_sratio[gene][mode][hpo] = genon_sratio
                        R.cadd_15_ratio[gene][mode][hpo] = cadd_15_ratio
                        R.genon_hs[gene][mode][hpo] = HGF * genon_sratio

                    for tp in ('HGF','genon_hratio','cadd_15_ratio','genon_sratio','genon_hs',):
                        getattr(rr,tp)[gene][mode][hpo] = \
                                getattr(R,tp)[gene][mode][hpo]

            # are hpos provided? if not, predict
            if not self.genes[gene].hpos:
                ps, ns = self.get_positive_negative_hpos(
                        rr.HGF[gene],
                        rr.genon_hratio[gene],
                        rr.cadd_15_ratio[gene],
                        rr.genon_sratio[gene],
                        rr.genon_hs[gene],
                        )
                # remove ns
                # note that ns is minised, so need to use ps to remove
                # any unwanted hpos
                rr.ph[gene] = ps
                for mode in ('d','r'):
                    for hpo in hpos:
                        if hpo not in ps[mode]:
                            for tp in (
                                    'HGF',
                                    'genon_hratio',
                                    'cadd_15_ratio',
                                    'genon_sratio',
                                    'genon_hs',
                                    ):
                                # do not pop. refer to self.ph for positive hpos
                                pass
                                #getattr(rr,tp)[gene][mode].pop(hpo)

            # predict inheritance mode
            rr.predicted_moi[gene] = self.predict_moi(
                    rr.HGF[gene],
                    rr.genon_hratio[gene],
                    rr.cadd_15_ratio[gene],
                    rr.genon_sratio[gene],
                    rr.genon_hs[gene],
                    rr.ph[gene],
                    )

        return rr

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

if __name__ == '__main__':
    G = Genon()
    result = G.analyse()
    print(result)
