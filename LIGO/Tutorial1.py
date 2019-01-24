#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 17:20:42 2018

@author: cg411
"""

# Standard python numerical analysis imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import h5py
import json
import sys

# the IPython magic below must be commented out in the .py file, since it doesn't work there.
#%matplotlib inline
#%config InlineBackend.figure_format = 'retina'
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

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

# -- Use the URL for a data file you found above
#url = 'https://losc.ligo.org/archive/data/O1/1125122048/H-H1_LOSC_4_V1-1126076416-4096.hdf5'
#filename = url.split('/')[-1]
#download(url)


import readligo as rl
#filename = 'H-H1_LOSC_4_V1-1126076416-4096.hdf5'
## -- Use the loaddata() method here to load strain data
#strain, time, dq_dict = rl.loaddata(filename, 'H1')
#
## -- Plot the first 1000 samples of strain vs. time here
#numSamples = 500
#start= 0
#stop = 100000
#plt.plot(time[start:stop], strain[start:stop])
#plt.xlabel('GPS Time (s)')
#plt.ylabel('H1 Strain')
#
##how to get the url in short form ? 
##how long (time is htis sample if 16kz and 16777216, then 17 mins ???)
#
#
##calculating sample rate 
#delt = time[1]-time[0]
#sample_rate = int(1/delt)
#length = 128 
## want to find the number of samples that corresponds to 128 seconds of data
#
#strain_seg = strain[0:length*sample_rate]
#time_seg = time[0:length*sample_rate]
#
##-- Plot a PSD with 128 seconds of data 
#ts = time[1] - time[0]      #-- Time between samples
#fs = int(1.0 / ts)          #-- Sampling frequency
#length = 128                #-- Number of seconds
#
#strain_seg = strain[0:(length*fs)]
#time_seg = time[0:(length*fs)]
#
#Pxx, freqs = mlab.psd(strain_seg, Fs=fs, NFFT=2*fs)
#plt.figure(2)
#plt.loglog(freqs, Pxx)
#plt.axis([10, 2000, 1e-47, 1e-38])
#plt.grid('on')
#plt.ylabel('Noise Power Spectra Density (PSD) [1/Hz]')
#plt.xlabel('Freq (Hz)')
#
## corresponding strain and time data 
#-- Download a data file containing GW150914 
url = 'https://losc.ligo.org/s/events/GW150914/H-H1_LOSC_4_V2-1126259446-32.hdf5'
fn_150914 = download(url)
strain, time, chan_dict_H1 = rl.loaddata(fn_150914, 'H1')

# the time sample interval (uniformly sampled!)
t0 = 1126259462.43
dt = time[1] - time[0]
fs = int(np.round(1/dt))
rel_time = time - t0

#-- How much data to use for the ASD?
deltat = 15  # Number of seconds on each side of data
N_samp = deltat*fs

# -- Center the PSD segment on the requested time
indx = np.where(np.abs(rel_time) < dt)[0][0]
strain_seg = strain[indx-N_samp : indx+N_samp]
time_seg = rel_time[indx-N_samp : indx+N_samp]

# number of sample for the fast fourier transform:
NFFT = 1*fs
fmin = 10
fmax = 2000

# -- Calculate PSD
Pxx, freqs = mlab.psd(strain_seg, Fs = fs, NFFT=NFFT, 
                      noverlap=NFFT/2, window=np.blackman(NFFT))

# We will use interpolations of the ASDs computed above for whitening:
psd = interp1d(freqs, Pxx)

# -- Whiten
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

# now whiten the data
strain_whiten = whiten(strain_seg,psd,dt)

# We need to suppress the high frequencies with some bandpassing:
high_freq = 600.
low_freq  = 30.
bb, ab = butter(4, [low_freq*2./fs, high_freq*2./fs], btype='band')
strain_whitenbp = filtfilt(bb, ab, strain_whiten)

#-- Plot the whitened time series
fig2 = plt.figure()
plt.plot(time_seg,strain_whitenbp,'r',label='H1 strain')

plt.xlim([-0.1,0.05])
plt.ylim([-10,10])
plt.xlabel('time (s) since '+str(t0))
plt.ylabel('whitented strain')
plt.legend(loc='lower left')
plt.title('WHITENED strain')








