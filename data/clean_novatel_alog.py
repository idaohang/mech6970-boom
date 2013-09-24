#!/usr/bin/env python
import os

alog_dirty = os.path.join(os.path.dirname(__file__), 'Novatel_Data_16_9_2013_____13_30_04.alog')
alog_clean = os.path.join(os.path.dirname(__file__), 'Novatel_Data_16_9_2013_____13_30_04___cleaned.alog')

_in = open(alog_dirty, 'r')
_out = open(alog_clean, 'w')

bad_meas = ['zNumObs', 'zNumL1Used']

for line in _in:
  if any(k in line for k in bad_meas):
    continue
  _out.write(line)

_in.close()
_out.close()
