#!/usr/bin/env python
""" Given the full path within the CDDIS FTP server, downloads a file to
specified directory and unzips it
You need to have gzip installed on your system.

Usage:

python download_cddis_data.py pub/gps/data/daily/2013/brdc/brdc2530.13n.Z /tmp

"""
import sys, os
from ftplib import FTP

# # locations
ftp_file = sys.argv[1]
out_dir = os.path.abspath(sys.argv[2])
out_filename_cmp = os.path.join(out_dir,'data_compressed.Z')
out_filename = os.path.join( out_dir, ftp_file.split('/')[-1].rstrip('.gz') )

# # download
cddis = FTP('cddis.gsfc.nasa.gov')
cddis.login()
out_file_cmp = open(out_filename_cmp,'wb')
cddis.retrbinary('RETR %s' % ftp_file, out_file_cmp.write)
cddis.close()
out_file_cmp.close()

# # unzip
if os.name=='nt':
	os.system(os.path.dirname(os.getcwd())+'\data\gzip.exe -d '+out_filename_cmp)
else:
	os.system('gzip -d '+out_filename_cmp)
os.rename(out_filename_cmp.rstrip('.Z'),out_filename)