#!/usr/bin/env python

import os.path
import json
import argparse
import math
from matplotlib import use; use('Agg')
from matplotlib import pyplot as plt
from ligo.gracedb.rest import GraceDb
from pycbc import mchirp_area
from pycbc.cosmology import _redshift
import mc_area_plots

def z_estimation(dist, dist_std):
    z_estimation = _redshift(dist)
    z_est_max = _redshift(dist + dist_std)
    z_est_min = _redshift(dist - dist_std)
    z_std_estimation = 0.5 * (z_est_max - z_est_min)
    z = {'central': z_estimation, 'delta': z_std_estimation}
    return z

parser = argparse.ArgumentParser()
parser.add_argument('--event', help='GraceDB event ID.')
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

g = GraceDb()
ev = g.event(args.event).json()
superevent = ev['superevent']
pipeline = ev['pipeline']
mchirp = ev['extra_attributes']['CoincInspiral']['mchirp']
trig_mc_det = {'central': mchirp, 'delta': mchirp * coeff['m0']}

if args.no_redshift:
    z = {'central': 0.0, 'delta': 0.0}
    dir_name = 'p_no_redshift'
else:
    if pipeline == 'pycbc':
        ifos = [ifos['ifo'] for ifos in ev['extra_attributes']['SingleInspiral']]
        eff_dist = [ifo['eff_distance'] if 'eff_distance' in ifo else '-' for
                    ifo in ev['extra_attributes']['SingleInspiral']]
        min_eff_distance = min(eff_dist)
        coinc_snr = ev['extra_attributes']['CoincInspiral']['snr']
        dist_estimation = coeff['a0'] * min_eff_distance
        dist_std_estimation = (dist_estimation * math.exp(coeff['b0']) *
                               coinc_snr ** coeff['b1'])
        z = z_estimation(dist_estimation, dist_std_estimation)
        dir_name = 'gracedb_p_dist_estimated'        
    else:
        print('Event %s does not have available effective distances.'
              'Analysis will be performed with no redshift estimation.'
               % args.event)
        z = {'central': 0.0, 'delta': 0.0}
        dir_name = 'gracedb_p_no_redshift'

areas = mchirp_area.calc_areas(trig_mc_det, mass_limits, mass_bdary, z,
                               separate_mass_gap)
total_area = sum(areas.values())
probabilities = {key: areas[key]/total_area for key in areas}

dirs = [out_dir + dir_name + '/' for out_dir in args.output_dir]
for dir in dirs:
    if not os.path.exists(dir):
        os.makedirs(dir)

name = str('%s-%s-%s.json' % (superevent, args.event, pipeline))
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

