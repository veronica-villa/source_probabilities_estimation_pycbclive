#!/usr/bin/env python
import os.path
import json
import argparse
import shutil
import subprocess
from ligo.gracedb.rest import GraceDb
from ligo.gracedb.rest import HTTPError
from astropy.io import fits

def gracedb_to_catalog_id(superevent_id):
    with open('GraceDb_ID-catalog_ID.json') as dicfile:
        dic_id = json.load(dicfile)
    return(dic_id[superevent_id])

parser = argparse.ArgumentParser()
parser.add_argument('--event', help='GraceDB event ID.')
args = parser.parse_args()

directories = ['xml_files', 'fits_files', 'event_data']
for dir in directories:
    if not os.path.exists(dir):
        os.mkdir(dir)

g = GraceDb()
ev = g.event(args.event).json()
catalog_id = gracedb_to_catalog_id(ev['superevent'])
pipeline = ev['pipeline']
mchirp = ev['extra_attributes']['CoincInspiral']['mchirp']

# get effective distances and coincident SNR for PyCBC events
if pipeline == 'pycbc':
    coinc_snr = ev['extra_attributes']['CoincInspiral']['snr']
    if args.event == 'G347304':
        eff_dist_dic = {'H1': 487.567, 'L1': 351.244}
    else:
        ifos = [ifos['ifo']
                for ifos in ev['extra_attributes']['SingleInspiral']]
        eff_dist = [ifos['eff_distance'] if 'eff_distance' in ifos else '-'
                    for ifos in ev['extra_attributes']['SingleInspiral']]
        eff_dist_dic = {k:v for (k,v) in zip(ifos, eff_dist)}

# get BAYESTAR distance and std from .fits file
if os.path.exists('fits_files/%s.fits' % args.event):
    pass
else:
    files = g.files(args.event).json()
    fits_file = [file for file in files if file.endswith('fits')]
    if len(fits_file) == 0:
        xml_file = [file for file in files if 'coinc' in file
                    and file.endswith('xml')][0]
        outfile = open('xml_files/%s.xml' % args.event,'wb')
        r = g.files(args.event, xml_file)
        outfile.write(r.read())
        outfile.close()
        file_path = 'xml_files/' + str(args.event) + '.xml'
        bayestar = ['bayestar-localize-coincs', file_path, file_path]
        subprocess.call(bayestar)
        file_name = 'fits_files/' + str(args.event) + '.fits'
        shutil.move('./0.fits', file_name)
    else:
        if len(fits_file) == 1:
            file_name = fits_file[0]
        else:
            file_name = [file for file in fits_file if
                         file.startswith('bayestar')][0]
        outfile = open('fits_files/%s.fits' % args.event,'wb')
        r = g.files(args.event, file_name)
        outfile.write(r.read())
        outfile.close()
hdul = fits.open('fits_files/%s.fits' % args.event)
for i in range(len(hdul)):
    for j in hdul[i].header.keys():
        if j == 'DISTMEAN':
            bay_dist = hdul[i].header[j]
        if j == 'DISTSTD':
            bay_std = hdul[i].header[j]

# save data
data = {'catalog_id': catalog_id, 'gracedb_event': args.event,
        'pipeline': pipeline, 'mchirp': mchirp, 'bay_dist': bay_dist,
        'bay_std': bay_std}
if pipeline == 'pycbc':
    data['eff_dist'] = eff_dist_dic
    data['coinc_snr'] = coinc_snr

outname = 'event_data/%s_data.json' % args.event
with open(outname, 'w') as outfile:
    json.dump(data, outfile, indent=1)
