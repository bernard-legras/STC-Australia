#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory from 7 January 2020
It select a few variables in a moving window. The data are stored
for all times in a pickle file. Reading this pickle file stores all the data
in memory allowing fast processing of diagnostics such as vortex tracking.

The data are stored in a dictionary indexed by numbers starting from 0.

Due to the unexpected life time of Koobor, the memory load is high and this
dataset cannot be processed on low memory system.

Created on Sat Feb  8 12:27:14 2020

The geopotential is managed in a separate history file and by get_datz

As the history became too big and took a long time to write
for each update, it has been cut in two parts

predates option allows to add new dates at the beginning of the series, requiring a shift
in the numbering.

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
#import numpy as np
import pickle,gzip
import os,psutil

# predates is to add new dates at the begiinning of the file
# it generates a special file which is then merged with
# the existing one in combine_dats.py
predates = False
update = True
day1 = datetime(2020,4,1,0)
day2 = datetime(2020,4,10,0)
if predates:
    day1 = datetime(2020,1,4,6)
    day2 = day1 + timedelta(days=3)
    update = False

if update==False:
# Initial rune
    dats = {}
    i = 0
else:
    # Continuation
    with gzip.open('ERA5-extract_ERA5PV.pkl','rb') as f:
        dats = pickle.load(f)
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
    if date >= datetime(2020,2,1,0): levs = (23,43)
    else: levs = (35,65)

    if date >= datetime(2020,3,29):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=varss,latRange=(-45,-15),lonRange=(30,150),levs=levs,copy=True)
    elif date >= datetime(2020,3,19):
        dats[i] = dat.extract(varss=varss,latRange=(-45,-15),lonRange=(80,200),levs=levs,copy=True)
    elif  date >= datetime(2020,3,9):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=varss,latRange=(-45,-15),lonRange=(-179,-59),levs=levs,copy=True)
    elif date >= datetime(2020,3,4):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=varss,latRange=(-45,-15),lonRange=(-100,20),levs=levs,copy=True)
    elif date >= datetime(2020,2,26):
        dats[i] = dat.extract(varss=varss,latRange=(-58,-28),lonRange=(0,120),levs=levs,copy=True)
    elif date >= datetime(2020,2,13):
        dats[i] = dat.extract(varss=varss,latRange=(-58,-28),lonRange=(90,210),levs=levs,copy=True)
    else:
        dats[i] = dat.extract(varss=varss,latRange=(-65,-35),lonRange=(180,300),levs=levs,copy=True)
    dat.close()
    date += timedelta(hours=6)
    i += 1

#%%
if predates:
    outfile = 'OPZ-extract-pre.pkl'
else:
    outfile = 'ERA5-extract_ERA5PV.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(dats,f)
