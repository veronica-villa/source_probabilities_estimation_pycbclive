#!/usr/bin/env python

from tqdm import tqdm
import subprocess
import json
from pycbc import mchirp_area
import mchirp_area_improved

src_args = {'mass_limits': {'max_m1': 45.0, 'min_m2': 1.0},
            'mass_bdary': {'ns_max': 3.0, 'gap_max': 5.0},
            'estimation_coeff': {'a0': 0.74899,
                                 'b0': -0.51557, 'b1': -0.32195,
                                 'm0': 0.1},
            'mass_gap': False}

events = ['G329276', 'G329515', 'G330308', 'G330564', 'G330694', 'G331327',
          'G332191', 'G332333', 'G333141', 'G333462', 'G333631', 'G333674',
          'G334993', 'G337426', 'G337515', 'G337913', 'G337998', 'G344656',
          'G345189', 'G345315', 'G347304', 'G348505', 'G348541', 'G350491',
          'G351437', 'G352016']

for event in tqdm(events):
    with open('event_data/%s_data.json' % event, 'r') as datafile:
        data = json.load(datafile)
    superevent = data['catalog_id']
    pipeline = data['pipeline']
    mchirp = data['mchirp']
    if pipeline == 'pycbc':
        eff_dist_dic = data['eff_dist']
        eff_dist = [eff_dist_dic[ifo] for ifo in eff_dist_dic.keys()]
        min_eff_dist = min(eff_dist)
        coinc_snr = data['coinc_snr']
    probs = mchirp_area.calc_probabilities(mchirp,
                                           coinc_snr, min_eff_dist, src_args)
    probs_imp = mchirp_area_improved.calc_probabilities(mchirp,
                                           coinc_snr, min_eff_dist, src_args)
    with open('sanity_check.txt', 'a') as outfile:
         outfile.write(superevent + '\n')
         json.dump(probs, outfile)
         outfile.write('\n')
         json.dump(probs_imp, outfile)
         outfile.write('\n' + '----------------------------------------' + '\n')
