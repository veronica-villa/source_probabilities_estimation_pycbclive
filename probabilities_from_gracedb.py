#!/usr/bin/env python
import os.path
import json
import subprocess
import shutil
import argparse
import click
import math
from ligo.gracedb.rest import GraceDb
from ligo.gracedb.rest import HTTPError
from astropy.io import fits
import mchirp_area
from pycbc.cosmology import _redshift

def check_superevents(superevents):
    """Checks if the superevents exist in the GraceDb database. Then, checks
       if they contain any pycbc trigger. Returns a list with the pycbc events
       corresponding with the superevents that satisfy both conditions.
    """
    for superev in superevents:
        try:
            g.superevent(superev)
        except HTTPError:
            print('Superevent %s does not exist!' %superev)
            if click.confirm('Do you want to enter another superevent ID?'):
                val = click.prompt('Please enter a valid superevent')
                (args.superevents).append(val)
                continue
        else:
            try_ev = [ev['graceid'] for ev in
                      g.events('superevent: %s pipeline: pycbc' % superev)]
            if len(try_ev) == 0:
                print('No pycbc trigger for superevent %s.' % superev)
                continue
            else:
                e_id.extend(try_ev)
    return e_id


def check_events(events):
    """Checks if the events exist in the GraceDb database and if they are
       pycbc triggers. Returns a list with the events that satisfy both
       conditions.
    """
    for event in args.events:
        try:
            ev = g.event(event).json()
        except HTTPError:
            print('Event %s does not exist!' % event)
            if click.confirm('Do you want to enter another pycbc event ID?'):
                val = click.prompt('Please enter a valid pycbc event')
                (args.events).append(val)
                continue
        else:
            if ev['pipeline'] != 'pycbc':
                print('%s is not a pycbc event.' % event)
            else:
                e_id.append(event)
    return e_id


def z_estimation(dist, dist_std):
    z_est = _redshift(dist)
    z_est_max = _redshift(dist + dist_std)
    z_est_min = _redshift(dist - dist_std)
    z_std_est = 0.5 * (z_est_max - z_est_min)
    return z_est, z_std_est


def get_data(event_id, args):
    """Gets the following properties of the event from GraceDb database:
       superevent ID, mchirp, coincident snr and effective distance for each
       ifo, the bayesian mean distance and its std. Then, computes other
       properties of the event: uncertainty of mchirp, and redshift and its
       uncertainty.
    """
    ev = g.event(event_id).json()
    s_id = ev['superevent']
    ifos = [ifos['ifo'] for ifos in ev['extra_attributes']['SingleInspiral']]
    eff_dist = [ifos['eff_distance'] if 'eff_distance' in ifos else -1
                for ifos in ev['extra_attributes']['SingleInspiral']]
    snr = ev['extra_attributes']['CoincInspiral']['snr']
    mchirp = ev['extra_attributes']['CoincInspiral']['mchirp']
    mchirp_delta = mchirp*args.mchirp_to_delta_coeff

    if args.estimate_lum_dist is True:
        dist_estimation = min(eff_dist) * args.eff_to_lum_distance_coeff
        b = args.lum_distance_to_delta_coeff
        dist_std_estimation = (dist_estimation * math.exp(b[0]) * snr ** b[1])
        z, z_delta = z_estimation(dist_estimation, dist_std_estimation)
        return {'s_id': s_id, 'e_id': event_id, 'mchirp': mchirp, 'mchirp_delta':
                mchirp_delta, 'snr': snr, 'z': z, 'z_delta': z_delta,
                'eff_dist': {ifo: eff_dist[n] for n, ifo in enumerate(ifos)},
                'dist_estimation': dist_estimation,
                'dist_std_estimation': dist_std_estimation}

    # If it is not True, obtain BAYESTAR distance and std from a .fits file
    if os.path.exists('fits_files/%s.fits' % event_id):
        pass
    else:
        files = g.files(event_id).json()
        fits_file = [file for file in files if file.endswith('fits')]
        if len(fits_file) == 0:

            xml_file = [file for file in files if 'coinc' in file
                        and file.endswith('xml')][0]

            xml_files = [file for file in files if file.endswith('xml')]
            for file in xml_files:
                out = open(file, 'wb')
                r = g.files(event_id, file)
                out.write(r.read())
                out.close()

            outfile = open('xml_files/%s.xml' % event_id,'wb')
            r = g.files(event_id,xml_file)
            outfile.write(r.read())
            outfile.close()
 
            file_path = 'xml_files/' + str(event_id) + '.xml'
            bayestar = ['bayestar-localize-coincs', file_path, file_path]
            subprocess.call(bayestar)
            file_name = 'fits_files/' + str(event_id) + '.fits'
            shutil.move('./0.fits', file_name)
        elif len(fits_file) == 1:
            file_name = fits_file[0]
        else:
            file_name = [file for file in fits_file if
                         file.startswith('bayestar')][0]
        outfile = open('fits_files/%s.fits' % event_id,'wb')
        r = g.files(event_id,file_name)
        outfile.write(r.read())
        outfile.close()
    hdul = fits.open('fits_files/%s.fits' % event_id)
    for i in range(len(hdul)):
        for j in hdul[i].header.keys():
            if j == 'DISTMEAN':
                bay_dist = hdul[i].header[j]
            if j == 'DISTSTD':
                bay_std = hdul[i].header[j]
    z = _redshift(bay_dist)
    dist_max = bay_dist + bay_std
    dist_min = bay_dist - bay_std
    z_delta = 0.5 * (_redshift(dist_max) - _redshift(dist_min))
    return {'s_id': s_id, 'e_id': event_id, 'mchirp': mchirp, 'mchirp_delta':
            mchirp_delta, 'snr': snr, 'z': z, 'z_delta': z_delta,
            'eff_dist': {ifo: eff_dist[n] for n,ifo in enumerate(ifos)},
            'bay_dist': bay_dist, 'bay_std': bay_std}

parser = argparse.ArgumentParser()
parser.add_argument('--superevents', nargs='+',
                    help='Superevent IDs.')
parser.add_argument('--events', nargs='+',
                    help='PyCBC event IDs.')
parser.add_argument('--estimate-lum-dist', action='store_true',
                    help='Probabilities are computed with luminosity '
                         'distances estimated from effective distances.')

mchirp_area.insert_args(parser)
args = parser.parse_args()

directories = ['xml_files', 'fits_files', 'probabilities']
for dir in directories:
    if not os.path.exists(dir):
        os.mkdir(dir)
g = GraceDb()
mass_limits = {'max_m1': args.max_m1, 'min_m2': args.min_m2}
mass_bdary = {'ns_max': args.ns_max, 'gap_max': args.gap_max}
e_id = []
e_id_bad_eff_dist = ['G347304', 'G351326']

if args.superevents:
    e_id = check_superevents(args.superevents)
if args.events:
    e_id = check_events(args.events)

for event in e_id:
    ev_data = get_data(event, args)
    trig_mc = {'central': ev_data['mchirp'], 'delta': ev_data['mchirp_delta']}
    z = {'central': ev_data['z'], 'delta': ev_data['z_delta']}
    prob_args = {'trig_mc': trig_mc, 'mass_limits': mass_limits, 'mass_bdary':
                 mass_bdary, 'z': z, 'mass_gap': args.mass_gap_separate}
    probabilities = mchirp_area.calc_probabilities(prob_args)
    name = (ev_data['s_id'] + '/' + str('%s_%s_MG-%s%s.txt' %
            (ev_data['s_id'], ev_data['e_id'], 'separate' if
             args.mass_gap_separate is not False else 'together',
             '_dist-estimated' if args.estimate_lum_dist is True else '')))
    out_name = 'probabilities/' + name
    pub_name = '/home/veronica.villa/public_html/mchirp_area/' + name
    names = [out_name, pub_name]
    for name in names:
        with open(name, 'w') as outfile:
            outfile.write('Source probabilities for %s \n' % ev_data['s_id'])
            outfile.write('data: ')
            json.dump(ev_data, outfile, indent=1)
            outfile.write('\n' + 'probabilities: ')
            json.dump(probabilities, outfile, indent=1)

print('Done!')
