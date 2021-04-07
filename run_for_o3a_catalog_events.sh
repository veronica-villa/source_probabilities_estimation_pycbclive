#!/usr/bin/env python

from tqdm import tqdm
import json
import subprocess

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
    mchirp = str(data['mchirp'])
    bayestar_distance = str(data['bay_dist'])
    bayestar_std = str(data['bay_std'])

    # estimated distance
    if pipeline == 'pycbc':
        eff_dist_dic = data['eff_dist']
        eff_dist = [eff_dist_dic[ifo] for ifo in eff_dist_dic.keys()]
        min_eff_dist = str(min(eff_dist))
        coinc_snr = str(data['coinc_snr'])

        dist_est_line = ['./probabilities_from_arguments.py',
                         '--superevent', superevent,
                         '--event', event,
                         '--pipeline', pipeline,
                         '--mchirp', mchirp,
                         '--min-eff-distance', min_eff_dist,
                         '--coinc-snr', coinc_snr,
                         '--eff-to-lum-distance-coeff', '0.74899',
                         '--lum-distance-to-delta-coeff',
                         '-0.51557', '-0.32195']
        subprocess.call(dist_est_line)

    # bayestar distance
    bay_dist_line = ['./probabilities_from_arguments.py',
                             '--superevent', superevent,
       	       	       	     '--event',	event,
                             '--pipeline', pipeline,
       	       	       	     '--mchirp', mchirp,
                             '--bayestar-distance', bayestar_distance,
                             '--bayestar-std', bayestar_std,
                             '--eff-to-lum-distance-coeff', '0.74899',
                             '--lum-distance-to-delta-coeff',
                             '-0.51557', '-0.32195']
    subprocess.call(bay_dist_line)

    # no redshift
    no_z_line = ['./probabilities_from_arguments.py',
                             '--superevent', superevent,
       	       	       	     '--event',	event,
                             '--pipeline', pipeline,
       	       	       	     '--mchirp', mchirp,
                             '--no-redshift',
                             '--eff-to-lum-distance-coeff', '0.74899',
                             '--lum-distance-to-delta-coeff',
                             '-0.51557', '-0.32195']
    subprocess.call(no_z_line)

print('Done!')
