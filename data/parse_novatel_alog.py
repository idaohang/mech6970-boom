#!/usr/bin/env python2.7
"""This is one way of putting *.alog data into a *.mat, if you don't like the 
script provided..."""
import os
import scipy.io as sio

alog_path = os.path.join(os.path.dirname(__file__), 'Novatel_Data_16_9_2013_____13_30_04.alog')
alog_file = open(alog_path, 'r')
mat_path =  os.path.join(os.path.dirname(__file__), 'Novatel_Data_16_9_2013_____13_30_04.mat')

data = {'gNovatel0': {}, 'gNovatel1': {}}

for line in alog_file:
  # inputting
  if '%' in line:
    continue
  l = line.split()
  moos_time = float(l[0])
  sensor = l[2]
  var_name = l[1].rstrip('_'+sensor).lstrip('z')
  value = l[3]

  # handle value formatting
  if '[' in value:
    value = [ float(x) for x in value.split('{')[1].rstrip('}').split(',') ]
  else:
    try:
      value = float(value) # weird value
    except Exception, e:
      continue

  # handle outputting
  if var_name not in data[sensor]:
    data[sensor][var_name] = {'moos_time': [], 'val':[]}

  data[sensor][var_name]['moos_time'].append(moos_time)
  data[sensor][var_name]['val'].append(value)



sio.savemat(mat_path,data)