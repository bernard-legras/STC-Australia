#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory of the second vortex.
Adapted from get_traject
This version includes the calculation of the geopotential.
See comments in get_traject.

Created on Sunday 20 September 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
#import numpy as np
import pickle,gzip
import os,psutil


# predates is used here to replace the 8 first days by a a box shifted northward
predates = False
update = False
day1 = datetime(2020,1,5,6)
day2 = datetime(2020,2,4,0)
if predates:
    day1 = datetime(2020,1,5,6)
    day2 = day1 + timedelta(days=10)
    update = False

if update==False:
# Initial run
    dats = {}
    i = 0
else:
    # Continuation
    with gzip.open('ERA5-extract-2ndVortex_ERA5PV.pkl','rb') as f:
        dats = pickle.load(f)
    #for i in range(44,54):
    #    del dats[i]
    i = len(dats)
    print('dats length ',len(dats))

date = day1

pid = os.getpid()
py = psutil.Process(pid)

while date < day2:
    print(date)
    dat = ECMWF('FULL-EA',date,exp='VOZ')
    dat._get_var('T')
    dat._get_var('VO')
    dat._get_var('O3')
    dat._get_var('U')
    dat._get_var('V')
    dat._mkpscale()
    dat._mkzscale()
    dat._mkp()
    dat._mkz()
    dat._mkthet()
    dat._mkpv()
    memoryUse = py.memory_info()[0]/2**30
    print('memory use after dat {:4.2f} gb'.format(memoryUse))
    varss = ['T','VO','O3','PV','Z','U','V','PT']
    levs = [40,75]
    if date >= datetime(2020,2,3):
        dats[i] = dat.extract(varss=varss,latRange=(-50,-20),lonRange=(130,250),levs=levs,copy=True)
    elif date >= datetime(2020,2,1):
        dats[i] = dat.extract(varss=varss,latRange=(-60,-30),lonRange=(160,280),levs=levs,copy=True)
    elif  date >= datetime(2020,1,27):
        dats[i] = dat.extract(varss=varss,latRange=(-60,-30),lonRange=(80,200),levs=levs,copy=True)
    elif date >= datetime(2020,1,22):
        dats[i] = dat.extract(varss=varss,latRange=(-60,-30),lonRange=(0,120),levs=levs,copy=True)
    elif date >= datetime(2020,1,18):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=varss,latRange=(-65,-35),lonRange=(-90,30),levs=levs,copy=True)
    elif date >= datetime(2020,1,15):
        dats[i] = dat.extract(varss=varss,latRange=(-65,-35),lonRange=(210,330),levs=levs,copy=True)
    else:
        dats[i] = dat.extract(varss=varss,latRange=(-50,-20),lonRange=(150,270),levs=levs,copy=True)
    dat.close()
    #del dat
    memoryUse = py.memory_info()[0]/2**30
    print('memory use after dat extract{:4.2f} gb'.format(memoryUse))
    date += timedelta(hours=6)
    i += 1

outfile = 'ERA5-extract-2ndVortex_ERA5PV.pkl'
if predates:
    with gzip.open('OPZ-extract-2ndVortex-post.pkl','rb') as f:
        dats2 = pickle.load(f)
    for i2 in range(len(dats)):
        dats2[i2] = dats[i2]
    with gzip.open(outfile,'wb') as f:
        pickle.dump(dats2,f)
else:
    with gzip.open(outfile,'wb') as f:
        pickle.dump(dats,f)

#%%
# date = datetime(2020,2,23,6)
# dat = ECMWF('OPZ',date)
# dat._get_var('VO')
# #datr = dat.shift2west(-179)
# dat2 = dat.extract(varss=['VO'],latRange=[-80,-40],lonRange=(270,320))
# dat2.show('VO',38)
