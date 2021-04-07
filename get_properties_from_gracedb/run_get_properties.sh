#!/usr/bin/env python

import subprocess 

# general analysis

events = ['G329276', 'G329515', 'G330308', 'G330564', 'G330694', 'G331327',
          'G332191', 'G332333', 'G333141', 'G333462', 'G333631', 'G333674',
          'G334993', 'G337426', 'G337515', 'G337913', 'G337998', 'G344656',
          'G345189', 'G345315', 'G347304', 'G348505', 'G348541', 'G350491',
          'G351437', 'G352016']

for event in events:
   line = ['./properties_from_gracedb.py', '--event', event]
   subprocess.call(line)

print('Done!')
