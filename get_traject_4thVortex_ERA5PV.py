#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajctory of the Fourth vortex
Adapted from get_traject
This version includes the calculation of the geopotential.
See comments in get_traject

Created on Wednesday 23 September 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
#import numpy as np
import pickle,gzip
import os,psutil

predates = False
update = False
day1 = datetime(2020,1,6,0)
day2 = datetime(2020,2,26,0)
if predates:
    day1 = datetime(2020,1,6,6)
    day2 = day1 + timedelta(days=23)
    update = False

if update==False:
# Initial rune
    dats = {}
    i = 0
else:
    # Continuation
    with gzip.open('ERA5-extract-4thVortex_ERA5PV.pkl','rb') as f:
        dats = pickle.load(f)
    print('dats length ',len(dats))
    i = len(dats)

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
    levs = (30,70)
    datr = dat.shift2west(-179)
    if date >= datetime(2020,2,24,6):
        dats[i] = datr.extract(varss=varss,latRange=(-70,-30),lonRange=(-60,60),levs=levs,copy=True)
    elif date >= datetime(2020,1,30,6):
        dats[i] = dat.extract(varss=varss,latRange=(-80,-40),lonRange=(270,340),levs=levs,copy=True)
    else:
        dats[i] = dat.extract(varss=varss,latRange=(-85,-55),lonRange=(180,310),levs=levs,copy=True)
    dat.close()
    del dat
    del datr
    date += timedelta(hours=6)
    i += 1

if predates:
    outfile = 'OPZ-extract-4thVortex-pre.pkl'
else:
    outfile = 'ERA5-extract-4thVortex_ERA5PV.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(dats,f)

#%%
# date = datetime(2020,2,23,6)
# dat = ECMWF('OPZ',date)
# dat._get_var('VO')
# #datr = dat.shift2west(-179)
# dat2 = dat.extract(varss=['VO'],latRange=[-80,-40],lonRange=(270,320))
# dat2.show('VO',38)
