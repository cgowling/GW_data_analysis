#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 15:27:10 2019

@author: cg411
"""

# exploring the quality values of any given data set 
import numpy as np
import matplotlib.pyplot as plt
import h5py

# 6 data quality categories in this data set, so 7 things to be one or zero 

def dataCheck (filename): # enter with quotation marks
    
    file = filename
    data = h5py.File(file, 'r')
    dqInfo = data['quality']['simple']
    qmask = dqInfo['DQmask'].value
    good = (qmask >> 4) & 1 # shifts along 4 so 4th index bit is on the right  then compares this to 1 
    # this is only checking one parameter 
    
    bitnameList = dqInfo['DQShortnames'].value
    nbits = len(bitnameList)
    
    for bit in range(nbits):
        print( bit, bitnameList[bit])
    
    return good
    
    
goodH1 = dataCheck('H-H1_LOSC_4_V1-1130541056-4096.hdf5')
goodL1 = dataCheck('L-L1_LOSC_4_V1-1130541056-4096.hdf5')

plt.plot(goodH1, label='H1')
plt.plot(goodL1, label='L1')
plt.legend(loc=1)


#if gpsStartH1 == gpsStartL1 :
#    print('GPS match')
#else:
#    print('No GPS ')





    
    
    
    
    
    
    