#!/usr/bin/env python

import json
import math
import numpy as np
from glob import glob
from matplotlib import use ; use('Agg')
import matplotlib.pyplot as plt

with open('distance_estimation_coeff.txt', 'r') as coeff_file:
    dist_estimation_coeff = float(coeff_file.read())

bay_std_est_dist_log = []
coinc_snr_log = []

for event in glob('event_data/*.json'):
    with open(event, 'r') as datafile:
        data = json.load(datafile)
    eff_dist_dic = data['eff_dist']
    eff_dist = [eff_dist_dic[ifo] for ifo in eff_dist_dic.keys()]
    min_eff_dist = min(eff_dist)
    est_distance = min_eff_dist * dist_estimation_coeff
    bayestar_std = data['bay_std']
    bay_std_est_dist = bayestar_std / est_distance
    coinc_snr = data['coinc_snr']
    bay_std_est_dist_log.append(np.log(bay_std_est_dist))
    coinc_snr_log.append(np.log(coinc_snr))

m, c = np.polyfit(coinc_snr_log, bay_std_est_dist_log, 1)

with open('dist_std_est_coeff.txt', 'w') as outfile:
   np.savetxt(outfile, [[m, c]], header='m c')

dist_std_est_log = m * np.array(coinc_snr_log) + c

plt.scatter(coinc_snr_log, bay_std_est_dist_log, zorder=1)
plt.plot(coinc_snr_log, dist_std_est_log, ls='--', lw=0.7, c='k', zorder=2,
         label='y = %.3f x + %.3f' % (m,c))
plt.xlabel('ln(SNR)')
plt.ylabel(r'ln($\sigma_{BAY}/D_{estimated}$)')
plt.legend()
plt.savefig('/home/veronica.villa/public_html/source_probabilities/'
            + 'plots_for_src_paper/distance_std_fit/'
            + 'bay_std_divided_est_distance_snr.png')
plt.close()
