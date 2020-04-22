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
#from zISA import zISA
import pickle,gzip
#import matplotlib.pyplot as plt

# predates is to add new dates at the begiinning of the file
# it generates a special file which is then merged with
# the existing one in combine_dats.py
predates = False
update = True
day1 = datetime(2020,4,3,6)
day2 = day1 + timedelta(days=3)
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
    with gzip.open('OPZ-extract-2.pkl','rb') as f:
        dats = pickle.load(f)
    i = len(dats)+130
    print('dats length ',len(dats))

date = day1

while date < day2:
    print(date)
    dat = ECMWF('OPZ',date)
    dat._get_var('T')
    dat._get_var('VO')
    dat._get_var('O3')
    dat._mkpscale()
    dat._mkzscale()
    if date >= datetime(2020,3,29):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=['T','VO','O3'],latRange=(-45,-15),lonRange=(30,150))
    elif date >= datetime(2020,3,19):
        dats[i] = dat.extract(varss=['T','VO','O3'],latRange=(-45,-15),lonRange=(80,200))
    elif  date >= datetime(2020,3,9):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=['T','VO','O3'],latRange=(-45,-15),lonRange=(-179,-59))
    elif date >= datetime(2020,3,4):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=['T','VO','O3'],latRange=(-45,-15),lonRange=(-100,20))
    elif date >= datetime(2020,2,26):
        dats[i] = dat.extract(varss=['T','VO','O3'],latRange=(-58,-28),lonRange=(0,120))
    elif date >= datetime(2020,2,13):
        dats[i] = dat.extract(varss=['T','VO','O3'],latRange=(-58,-28),lonRange=(90,210))
    else:
        dats[i] = dat.extract(varss=['T','VO','O3'],latRange=(-65,-35),lonRange=(180,300))
    dat.close()
    date += timedelta(hours=12)
    i += 1

#%%
if predates:
    outfile = 'OPZ-extract-pre.pkl'
else:
    outfile = 'OPZ-extract-2.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(dats,f)