#!/usr/bin/env python

import os
import json
import shutil
import argparse
import math
import mchirp_area
from pycbc.cosmology import _redshift

def z_estimation(dist, dist_std):
    z_est = _redshift(dist)
    z_est_max = _redshift(dist + dist_std)
    z_est_min = _redshift(dist - dist_std)
    z_std_est = 0.5 * (z_est_max - z_est_min)
    return z_est, z_std_est

parser = argparse.ArgumentParser()
parser.add_argument('--mchirp', type=float)
parser.add_argument('--effective-distance', type=float)
parser.add_argument('--bayestar-distance', type=float)
parser.add_argument('--bayestar-delta', type=float)
parser.add_argument('--coinc-snr', type=float)
parser.add_argument('--superevent',
                    help='Superevent ID.')
parser.add_argument('--event',
                    help='PyCBC event ID.')
parser.add_argument('--no-redshift', action='store_true')

mchirp_area.insert_args(parser)
args = parser.parse_args()

mass_limits = {'max_m1': args.max_m1, 'min_m2': args.min_m2}
mass_bdary = {'ns_max': args.ns_max, 'gap_max': args.gap_max}

mchirp = args.mchirp
mchirp_delta = mchirp*args.mchirp_to_delta_coeff
trig_mc = {'central': mchirp, 'delta': mchirp_delta}

if args.no_redshift:
    z = 0.0
    z_delta = 0.0
else:
    if args.effective_distance:
        lum_dist = args.effective_distance * args.eff_to_lum_distance_coeff
        b = args.lum_distance_to_delta_coeff
        lum_dist_std = (lum_dist * math.exp(b[0]) * args.coinc_snr ** b[1])
    elif args.bayestar_distance:
        lum_dist = args.bayestar_distance
        lum_dist_std = args.bayestar_delta
    z, z_delta = z_estimation(lum_dist, lum_dist_std)

z = {'central': z, 'delta': z_delta}

data = {'superevent': args.superevent, 'event': args.event, 'mchirp': mchirp,
        'min_eff_dist': args.effective_distance, 'coinc_snr': args.coinc_snr}

out_dir = 'probabilities_args/' + args.superevent + '/'
pub_dir = ('/home/veronica.villa/public_html/source_probabilities/o3/' +
           args.superevent + '/')
dirs = [out_dir, pub_dir]
for dir in dirs:
    if not os.path.exists(dir):
        print(dir)
        os.makedirs(dir)

name = (str('%s_%s%s%s%s.txt' %
        (args.superevent, args.event, '_MG-separate' if
         args.mass_gap_separate is not False else '',
         '_no-z' if args.no_redshift else ('_dist-estimated' if args.effective_distance else '_bay-dist'))))

out_name = str(out_dir + name)
pub_name = str(pub_dir + name)
names = [out_name, pub_name]

prob_args = {'trig_mc': trig_mc, 'mass_limits': mass_limits, 'mass_bdary':
             mass_bdary, 'z': z, 'mass_gap': args.mass_gap_separate,
             'mass_plot': args.enable_mass_plot, 'probability_plot':
             args.enable_probability_plot, 'fnames': names}

probabilities = mchirp_area.calc_probabilities(prob_args)

for name in names:
    with open(name, 'w') as outfile:
        outfile.write('Source probabilities for %s - %s \n' % 
                      (args.superevent, args.event))
        outfile.write('\n' + 'data: ')
        json.dump(data, outfile, indent=1)
        outfile.write('\n' + 'probabilities: ')
        json.dump(probabilities, outfile, indent=1)

print('Done!')
