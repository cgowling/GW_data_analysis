#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 15:10:28 2018

@author: cg411
"""

import json


t0 = 968654557.955 # Time of big dog (S6 blind injection)
# t0 = 933661015.0  # Time of a successful s6 hardware injection in H1 (SNR 140)

dataset = 'S6'
detector = 'H1'
version = 'V1'  # V1 is "version 1" of the data release


observatory = detector[0]         # first letter of the detector H or L
hour        = int(t0)&0xFFFFF000  # the filename rounding down to a multiple of 4096
fortnight   = int(t0)&0xFFF00000  # the directory by rounding down to multiple of 4096*256
filename = '{0}-{1}_LOSC_4_{2}-{3}-4096.hdf5'.format(observatory, detector, version, hour)
urlformat = 'https://losc.ligo.org/archive/data/{0}/{1}/{2}'
url = urlformat.format(dataset, fortnight, filename)


pyversion = sys.version_info.major
if pyversion == 2: 
    import urllib2
else:
    import urllib.request

# -- Handy function to download data file, and return the filename
def download(url):
    filename = url.split('/')[-1]
    print('Downloading ' + url )
    if pyversion == 2: 
        r = urllib2.urlopen(url).read()
        f = open(filename, 'w')   # write it to the right filename
        f.write(r)
        f.close()
    else:
        urllib.request.urlretrieve(url, filename)  
    print("File download complete")
    return filename


download(url)

