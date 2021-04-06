#!/usr/bin/env python

import os
import json
import argparse
import numpy as np
from glob import glob
from tqdm import tqdm

def combine_probabilities(probabilities, pastro):
    for category in probabilities.keys():
        pastro[category] = probabilities[category] * pastro['Astro']
    return pastro

parser = argparse.ArgumentParser()
parser.add_argument('--pcbc-dirs', nargs='+')
parser.add_argument('--pastro-dirs', nargs='+')
args = parser.parse_args()

out_dir = 'combined_results/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

pcbc_files = [i for dir in args.pcbc_dirs
              for i in glob(dir + '*.json')]
pastro_files = [i for dir in args.pastro_dirs 
              for i in glob(dir + '*.json')]
pastro_gps_times = np.array([int((pastro.split('/')[-1]).split('-')[2])
                    for pastro in pastro_files])

for event in tqdm(pcbc_files):
    event_name = event.split('/')[-1]
    gps_time = int(event_name.split('-')[2])
    gps_dif = abs(pastro_gps_times - gps_time)
    idx = np.argmin(gps_dif)
    if gps_dif[idx] <= 1:
        with open(event, 'r') as probs_file:
            probs_data = json.load(probs_file)
        with open(pastro_files[idx], 'r') as pastro_file:
            pastro_data = json.load(pastro_file)
        combined_p = combine_probabilities(probs_data, pastro_data)
        with open(out_dir + event_name, 'w') as outfile:
            json.dump(combined_p, outfile)
