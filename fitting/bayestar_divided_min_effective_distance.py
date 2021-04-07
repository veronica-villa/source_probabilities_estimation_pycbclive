#!/usr/bin/env python

import json
import numpy as np
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

bayestar_divided_min_eff_distances = []
bayestar_distances = []
mchirps = []
plt.xscale('log')

for filename in sorted(glob('event_data/*.json')):
    with open(filename, 'r') as datafile:
        data = json.load(datafile)
    eff_dist_dic = data['eff_dist']
    eff_dist = [eff_dist_dic[ifo] for ifo in eff_dist_dic.keys()]
    min_eff_dist = min(eff_dist)
    bay_dist = data['bay_dist']
    mchirp = data['mchirp']
    bayestar_distances.append(bay_dist)
    bayestar_divided_min_eff_distances.append(bay_dist/min_eff_dist)
    mchirps.append(mchirp)

mean = np.mean(bayestar_divided_min_eff_distances)

with open('distance_estimation_coeff.txt', 'w') as outfile:
    json.dump(mean, outfile)

plt.scatter(bayestar_distances, bayestar_divided_min_eff_distances)
plt.axhline(mean, ls='--', lw=0.7, label='mean = %.3f' % mean)
plt.legend()
plt.xlabel('BAYESTAR distance (Mpc)')
plt.ylabel('BAYESTAR distance / Minimum PyCBC effective distance')
plt.savefig('/home/veronica.villa/public_html/source_probabilities/'
            + 'plots_for_src_paper/distance_fit/bayestar_divided_min_effective_distance_mean.png')
plt.close()
