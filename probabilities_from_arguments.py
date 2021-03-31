#!/usr/bin/env python

import os
import json
import shutil
import argparse
import math
from matplotlib import use; use('Agg')
from matplotlib import pyplot as plt
from pycbc import mchirp_area
from pycbc.cosmology import _redshift
from pycbc.conversions import mass2_from_mchirp_mass1 as m2mcm1
import mc_area_plots

def z_estimation(dist, dist_std):
    z_estimation = _redshift(dist)
    z_est_max = _redshift(dist + dist_std)
    z_est_min = _redshift(dist - dist_std)
    z_std_estimation = 0.5 * (z_est_max - z_est_min)
    z = {'central': z_estimation, 'delta': z_std_estimation}
    return z

parser = argparse.ArgumentParser()
parser.add_argument('--superevent', help='Superevent ID.')
parser.add_argument('--event', help='Event ID.')
parser.add_argument('--pipeline')
parser.add_argument('--mchirp', type=float)
parser.add_argument('--min-eff-distance', type=float)
parser.add_argument('--coinc-snr', type=float)
parser.add_argument('--bayestar-distance', type=float)
parser.add_argument('--bayestar-std', type=float)
parser.add_argument('--no-redshift', action='store_true')
parser.add_argument('--output-dir', nargs='+')
parser.add_argument('--enable-mass-plot', action='store_true')
parser.add_argument('--enable-probability-plot', action='store_true')

mchirp_area.insert_args(parser)
mc_area_plots.insert_args(parser)
args = parser.parse_args()

mc_area_args = mchirp_area.from_cli(args)
contour_plot_args = mc_area_plots.from_cli(args)
mass_limits = mc_area_args['mass_limits']
mass_bdary = mc_area_args['mass_bdary']
coeff = mc_area_args['estimation_coeff']
separate_mass_gap = mc_area_args['mass_gap']

trig_mc_det = {'central': args.mchirp, 'delta': args.mchirp * coeff['m0']}

#estimate redshift
if args.no_redshift:
    z = {'central': 0.0, 'delta': 0.0}
    dir_name = 'p_no_redshift'
else:
    if args.min_eff_distance and args.coinc_snr:
        dist_estimation = coeff['a0'] * args.min_eff_distance
        dist_std_estimation = (dist_estimation * math.exp(coeff['b0']) *
                               args.coinc_snr ** coeff['b1'])
        z = z_estimation(dist_estimation, dist_std_estimation)
        dir_name = 'p_dist_estimated'
    elif args.bayestar_distance and args.bayestar_std:
        z = z_estimation(args.bayestar_distance, args.bayestar_std)
        dir_name = 'p_dist_bayestar'

areas = mchirp_area.calc_areas(trig_mc_det, mass_limits, mass_bdary, z,
                               separate_mass_gap)
total_area = sum(areas.values())
probabilities = {key: areas[key]/total_area for key in areas}

dirs = [out_dir + dir_name + '/' for out_dir in args.output_dir]
for dir in dirs:
    if not os.path.exists(dir):
        os.makedirs(dir)

name = str('%s-%s-%s.json' % (args.superevent, args.event, args.pipeline))
names = [dir + name for dir in dirs]
for name in names:
    with open(name, 'w') as outfile:
        json.dump(probabilities, outfile)

if args.enable_mass_plot:
    massplot = mc_area_plots.contour_mass_plot(trig_mc_det, mass_limits,
                                                   mass_bdary, z,
                                                   contour_plot_args)
    for name in names:
        plot_name = name.replace('.json', '_mass_plot.png')
        massplot.savefig(plot_name)

if args.enable_probability_plot:
    probsplot = mc_area_plots.probabilities_plot(probabilities)
    for name in names:
        plot_name = name.replace('.json', '_probabilities_pie.png')
        probsplot.savefig(plot_name)
