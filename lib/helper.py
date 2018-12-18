from collections import Counter,defaultdict
import re
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp, binom, gamma
import subprocess
import sys
import itertools
import math
import fisher


'''
get Ph for all HPO terms
'''
def get_phs(p_info):
    result = Counter()
    for k,v in p_info.items():
        result.update(v['hpo'])
    return result

'''
get all ancestor nodes of a given hpo_id.
'''
def get_hpo_ancestors(hpo_db, hpo_id):
    """
    Get HPO terms higher up in the hierarchy.
    """
    h=hpo_db[hpo_id]
    #print(hpo_id,h)
    if 'replaced_by' in h:
        # not primary id, replace with primary id and try again
        h = hpo_db[h['replaced_by'][0]]
    hpo=[h]
    if 'is_a' not in h: return hpo
    for hpo_parent_id in h['is_a']:
        #p=hpo_db.hpo.find({'id':hpo_parent_id}):
        hpo+=list(itertools.chain(get_hpo_ancestors(hpo_db,hpo_parent_id)))
    #remove duplicates
    hpo={h['id'][0]:h for h in hpo}.values()
    return hpo

'''
get minimised set of HPO terms
'''
def hpo_minimum_set(hpo_db, hpo_ids):
    '''
    minimize the hpo sets
    '''
    hpo_ids = list(set(hpo_ids))
    results = dict([(hpo_id, [ h['id'][0] for h in get_hpo_ancestors(hpo_db, hpo_id)],) for hpo_id in hpo_ids])
    # minimise
    bad_ids = []
    for i in range(len(hpo_ids)):
        for j in range(i+1,len(hpo_ids)):
            if hpo_ids[i] in results[hpo_ids[j]]:
                # i is j's ancestor, remove
                bad_ids.append(hpo_ids[i])
                break
            if hpo_ids[j] in results[hpo_ids[i]]:
                # j is i's ancestor, remove
                bad_ids.append(hpo_ids[j])
    return list(set(hpo_ids) - set(bad_ids))

'''
given positive sets, select negative sets from the overall pool
'''
def get_negative_hpos(hpo_db,ps,hpos):
    result = []
    for h1 in hpos:
        bad = 0
        for h2 in ps:
            A = h1 in [i['id'][0] for i in get_hpo_ancestors(hpo_db, h2)]
            B = h2 in [i['id'][0] for i in get_hpo_ancestors(hpo_db, h1)]
            if A or B:
                bad = 1
                break
        if not bad:
            result.append(h1)
    return hpo_minimum_set(hpo_db,result)

'''
iter through string separated by \n
'''
def split_iter(string):
    return (x.group(0) for x in re.finditer(r"[^\n]+", string))

'''
get patient hpo terms as a dict
'''
def get_snapshot(f):
    dt = {}
    with open(f,'r') as inf:
        for row in inf:
            if row[0] == '#': continue
            row = row.rstrip().split('\t')
            if row[1] == '0': continue
            dt[row[0]] = {
                    'hpo': row[2].split(','),
                    'contact': row[3],
            }
    return dt
'''
check compulsory args
'''
def check_args(compulsory_keys, kwargs, func):
    # a function to check args
    missing = compulsory_keys - set(kwargs.keys())
    if missing:
        raise KeyError("missing keys for {}: {}".format(func, ', '.join(missing)))

'''
get artefact variants
'''
def get_batch_artefacts(**kwargs):
    # check compulsory args
    compulsory_args = {
            'data', # the genotype data
            'patient_mini', # the phenotype data,
                            #with cohort id encoded in contact
        }
    check_args(compulsory_args,kwargs,'get_batch_artefacts')
    # set some default
    # lower_bound is there to remove cohorts where there is just one patient
    # zero_gnomad_c_cutoff allows max internal count when gnomad_af is 0

    # for example, '19-44777406-G-A' and '19-44740306-C-CA' are
    # Jewish variants, and found predominantly in SEGAL's cohort,
    # which is known to have many Jews. These variants would significantly
    # contributes to the false positive enrichments.
    # ideally one would set a lower binom_cutoff, removing those variants
    # by doing a population study on different cohorts, then remove
    # suspicious variants by looking at gnomad population af
    #
    optional = dict(
            lower_bound = 2,
            zero_gnomad_c_cutoff = 2,
            binom_cutoff = 1e-4,
            )
    for k,v in optional.items():
        kwargs.setdefault(k, v)

    # dt_d: for each cohort, each variant's appearance freq, dominance mode
    # dt_r for receissive mode
    dt_d = defaultdict(Counter)
    dt_r = defaultdict(Counter)
    cohorts = Counter()
    for k,v in kwargs['data']['patients'].items():
        cohorts[ kwargs['patient_mini'][k]['contact'] ] += 1
        vc = Counter(v)
        for i in vc:
            dt_d[ kwargs['patient_mini'][k]['contact'] ][i] += 1
            if vc[i] > 1:
                dt_r[ kwargs['patient_mini'][k]['contact'] ][i] += 1
    # remove cohorts with count lower than lower_bound
    for k in cohorts.keys():
        if cohorts[k] < kwargs['lower_bound']:
            del cohorts[k]
            dt_d.pop(k,None)
            dt_r.pop(k,None)

    # for heterozygous variants
    result_d = defaultdict(list)
    for k1,v1 in dt_d.items():
        n_variants = len(v1)
        for k2,v2 in v1.items():
            if not kwargs['data']['variants'][k2]['gnomad_af']:
                if v2 > kwargs['zero_gnomad_c_cutoff']:
                    result_d[k1].append(k2)
                continue
            prob = 1 - binom.cdf(v2-1, cohorts[k1], kwargs['data']['variants'][k2]['gnomad_af'])
            if prob < kwargs['binom_cutoff'] / n_variants:
                result_d[k1].append(k2)
    for k in result_d:
        result_d[k] = set(result_d[k])

    # for homozygous variants
    result_r = defaultdict(list)
    for k1,v1 in dt_r.items():
        n_variants = len(v1)
        for k2,v2 in v1.items():
            if not kwargs['data']['variants'][k2]['gnomad_hom_f']:
                if v2 > kwargs['zero_gnomad_c_cutoff']:
                    result_r[k1].append(k2)
                continue
            prob = 1 - binom.cdf(v2-1, cohorts[k1], kwargs['data']['variants'][k2]['gnomad_hom_f'])
            if prob < kwargs['binom_cutoff'] / n_variants:
                #print(k2,prob)
                result_r[k1].append(k2)
    for k in result_r:
        result_r[k] = set(result_r[k])
    return {'d':result_d,'r':result_r}

'''
remove batch specific artefact variants
'''
def remove_batch_artefacts(data, bad_vs, patient_mini, mode='all'):
    result = {
            'patients':{},
            'variants':data['variants'],
            #'gene_id':data['gene_id'],
            #'pat_a':data['pat_a'],
            }
    bad_p = []
    for k1,v1 in data['patients'].items():
        cohort = patient_mini[k1]['contact']
        this_bad_vs = []
        # collect het artefacts
        if mode != 'r' and cohort in bad_vs['d']:
            this_bad_vs += [i for i in v1 if i in bad_vs['d'][cohort]]
        # collect hom artefacts
        if mode != 'd' and cohort in bad_vs['r']:
            vc = Counter(v1)
            for k in vc:
                if vc[k] > 1 and k in bad_vs['r'][cohort]:
                    this_bad_vs.append(k)
        this_vs = [i for i in v1 if i not in this_bad_vs]
        if this_vs:
            result['patients'][k1] = this_vs
    return result

def remove_noncoding(data, coding_variants):
    for v in data['variants'].keys():
        if v not in coding_variants:
            data['variants'].pop(v)
    for p in data['patients'].keys():
        for v in list(data['patients'][p]):
            if v not in coding_variants:
                data['patients'][p].remove(v)

'''
add cadd from the already annotated cadd file, using tabix
Note that this is for getting variants from a gene.
Not good for getting random variants
'''
def add_cadd(**kwargs):
    compulsory_args = {
            'variants',
            'chrom',
            'start',
            'stop',
            'cadd_file',
            }
    check_args(compulsory_args,kwargs,'add_cadd')
    # get grange for tabix
    grange = '{chrom}:{start}-{stop}'.format(**kwargs)
    # get the correct cadd file
    cadd_file = kwargs['cadd_file'].format(kwargs['chrom'])
    # get cadd
    result = {i:None for i in kwargs['variants']}
    cadd_file = subprocess.check_output(('tabix',cadd_file,grange))
    for row in split_iter(cadd_file):
        row = row.rstrip().split('\t')
        v_id = '-'.join(row[:4])
        if v_id in kwargs['variants']:
            result[v_id] = float(row[-1])
    return result

'''
given genotype (including cadd and gnomad) and the vcf file, make a patient map
for LI profiling
'''
def get_patient_map(**kwargs):
    compulsory_args = {
            'data',
            'vcf',
            'mode',
            'gnomad_steps',
            'cadd_steps',
            }
    check_args(compulsory_args, kwargs, 'get_patient_map')
    # if phase is provided, ignore cis_gap.
    # cis_gap has to be provided if no phase
    if not {'phase','cis_gap'} - set(kwargs.keys()):
        msg = 'get_patient_map needs at least one argument \
                from (phase, cis_gap)'
        raise KeyError(msg)
    patient_map = defaultdict(list)
    args = kwargs.copy()
    for i in range(len(kwargs['cadd_steps'])-1):
        p = []
        for j in range(len(kwargs['gnomad_steps'])-1):
            args['gr'] = (kwargs['gnomad_steps'][j], kwargs['gnomad_steps'][j+1])
            args['cr'] = (kwargs['cadd_steps'][i], kwargs['cadd_steps'][i+1])
            p,narrow_vs,not_covered_patients = get_patients(**args)
            patient_map[(i,j)] = (p,not_covered_patients)
    return patient_map

'''
get variants within a given bin
for recessive, the second variant needs to be no more frequent
than the first selected variant, and its cadd no lower than the
first variant
'''
def get_variants(**kwargs):
    compulsory_args = {
            'variants',
            'mode',
            'gr',
            'cr',
            }
    check_args(compulsory_args, kwargs, 'get_variants')
    mode_dict = {'r':'gnomad_hom_f','d':'gnomad_af'}
    narrow_vs = (k for k,v in kwargs['variants'].items()
            if kwargs['gr'][0] <= v[mode_dict[kwargs['mode']]]<kwargs['gr'][1]
            and kwargs['cr'][0] <= v['cadd']<kwargs['cr'][1]
            )
    broad_vs = tuple()
    if kwargs['mode'] == 'r':
        broad_vs = (k for k,v in kwargs['variants'].items()
                if v[mode_dict[kwargs['mode']]] < kwargs['gr'][1]
                and v['cadd'] >= kwargs['cr'][0]
                )
    return (set(narrow_vs),set(broad_vs))

'''
first find valid variants, and then look at vcf to get covered individuals. return matching patients, and also uncovered patients to correct both p_a and p_h
miss_cutoff: ratio of miss called variants of an individual to be considered to be removed
if no phase, fall back to gap to remove putative cis variant
'''
def get_patients(**kwargs):
    compulsory_args = {
            'data',
            'vcf',
            'mode',
            'gr',
            'cr',
            }
    check_args(compulsory_args,kwargs,'get_patients')
    # if phase is provided, ignore cis_gap.
    # cis_gap has to be provided if no phase
    if not {'phase','cis_gap'} - set(kwargs.keys()):
        msg = 'get_patients needs at least one argument \
                from (phase, cis_gap)'
    # add some defaults
    default = dict(
            miss_cutoff = 0.5,
            cis_gap = 100,
            )
    for k,v in default.items():
        kwargs.setdefault(k,v)

    args = kwargs.copy()
    args['variants']= kwargs['data']['variants']
    p = []
    narrow_vs,broad_vs = get_variants(**args)

    if not narrow_vs:
        return ([],narrow_vs,set())

    # get patients not covered on all narrow_vs
    s = kwargs['vcf'].loc[narrow_vs].sum()
    not_covered_patients = set(s[s<len(narrow_vs)*kwargs['miss_cutoff']].index)

    if kwargs['mode'] == 'd':
        p = [k for k,v in kwargs['data']['patients'].items() if set(v) & narrow_vs]
    elif kwargs['mode'] == 'r':
        for k,v in kwargs['data']['r_patients'].items():
            good = []
            other = []
            for i in v:
                if i in narrow_vs:
                    good.append(i)
                elif i in broad_vs:
                    other.append(i)
            # remove hom in other, as hom in other renders the variants in good unnecessary
            other = [k1 for k1,v1 in Counter(other).items() if v1 == 1]
            r_bad = [] # for removing linked variants
            if len(good) and len(good+other) > 1:
                pos = [int(i.split('-')[1]) for i in good+other]
                if (len(pos) > len(set(pos))) or (max(pos) - min(pos) > kwargs['cis_gap']):
                    if kwargs.get('phase_cutoff',None):
                        # any of them is hom?
                        if len(set(good)) < len(good):
                            p.append(k)
                            continue
                        if len(good+other) < 2: continue
                        for i in itertools.combinations(sorted(good+other,key=lambda x: int(x.split('-')[1])),2):
                            if not set(good) & set(i): continue
                            cis_p = kwargs['phase'][k][i] if i in kwargs['phase'][k] else 0
                            if cis_p < kwargs['phase_cutoff']:
                                p.append(k)
                                break
                    else:
                        p.append(k)
    return (p,narrow_vs,not_covered_patients)

'''
calculate phenogenon, and return a dataframe
'''
def phenogenon(**kwargs):
    compulsory_args = {
            'hpos',
            'mode',
            'patient_map',
            'patient_info',
            }
    # check args
    check_args(compulsory_args, kwargs, 'phenogenon')
    # get raw p_a and p_h
    raw_p_a = set(kwargs['patient_info'].keys())
    hpos = set(kwargs['hpos'].split(','))
    raw_p_h = set([kk for kk,vv in kwargs['patient_info'].items() if set(hpos).issubset(set(vv['hpo']))])
    # make a np matrix
    shape1,shape2 = [],[]
    for k in kwargs['patient_map']:
        shape1.append(k[0])
        shape2.append(k[1])
    logp_df = np.zeros( (
        len(set(shape1)),
        len(set(shape2))
        ) )
    gamma_weights = logp_df.copy()
    for k,v in kwargs['patient_map'].items():
        p_a = len(raw_p_a - set(v[1]))
        p_h = len(raw_p_h - set(v[1]))
        p_g = len(v[0])
        if not p_g:
            logp_df[k[0]][k[1]] = None
            continue
        p_gh = 0
        for p in v[0]:
            if set(hpos).issubset(set(kwargs['patient_info'][p]['hpo'])):
                p_gh += 1
        pval = fisher.pvalue(
                p_a - p_h - p_g + p_gh,
                p_h - p_gh,
                p_g - p_gh,
                p_gh
                ).right_tail
        logp_df[k[0]][k[1]] = pval
    return logp_df


# get the mean different between positive set of HPOs and negative set of HPOs
def get_diff(ps,ns,z):
    if not ps or not ns:
        return None
    result = {i:None for i in ps}
    for h in ps:
        mean = np.zeros((len(z[h]),len(z[h][0])))
        p = np.zeros((len(z[h]),len(z[h][0])))
        for ind1,i1 in enumerate(mean):
            for ind2,i2 in enumerate(i1):
                pos = z[h][ind1][ind2]
                neg = []
                for k,v in z.items():
                    if k in ns:
                        neg.append(v[ind1][ind2])
                # update mean
                mean[ind1][ind2] = pos - np.mean(neg)
                # update p
                t = ttest_1samp(neg,pos)
                pv = (t.statistic > 0) - np.sign(t.statistic) * t.pvalue / 2.
                p[ind1][ind2] = 0 if pd.isnull(pv) else -math.log10(pv or 1e-50)
            result[h] = {'diff_mean':mean,'diff_p':p}
    return result

# use gamma_cdf to adjust diff_p to correct when p_g is low
def adjust_diff_p(patient_map,hz,gamma_k,gamma_scale):
    if hz is None:
        return None
    # make a matrix of adjust parameters from patient map
    #gamma_k = len_of_negative_hpo_set / 2
    adj_m = np.zeros(hz[hz.keys()[0]]['diff_p'].shape)
    for k,v in patient_map.items():
        adj_m[k[0]][k[1]] = gamma.cdf(len(v[0]),gamma_k,scale=gamma_scale)
    for k1,v1 in hz.items():
        v1['diff_p'] = np.multiply(v1['diff_p'], adj_m)
    return hz

