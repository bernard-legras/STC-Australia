#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajctory of the Fourth vortex
Adapted from get_traject
This version includes the calculation of the geopotential.
See comments in get_traject

Created on Sat Feb  8 12:27:14 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
#import numpy as np
#from zISA import zISA
import pickle,gzip
#import matplotlib.pyplot as plt

predates = True
update = False
day1 = datetime(2020,2,24,6)
day2 = day1 + timedelta(days=6)
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
    with gzip.open('OPZ-extract-4thVortex.pkl','rb') as f:
        dats = pickle.load(f)
    for i in range(70,82):
        del dats[i]
    print('dats length ',len(dats))
    i = len(dats)

date = day1

while date < day2:
    print(date)
    dat = ECMWF('OPZ',date)
    dat._get_var('T')
    dat._get_var('VO')
    dat._get_var('O3')
    dat._mkpscale()
    dat._mkzscale()
    dat._mkp()
    dat._mkz()
    datr = dat.shift2west(-179)
    if date >= datetime(2020,2,24,6):
        dats[i] = datr.extract(varss=['T','VO','O3','Z'],latRange=(-70,-30),lonRange=(-60,60))
    elif date >= datetime(2020,1,30,6):
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-80,-40),lonRange=(270,340))
    else:
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-85,-55),lonRange=(180,310))
    dat.close()
    del dat
    del datr
    date += timedelta(hours=12)
    i += 1

if predates:
    outfile = 'OPZ-extract-4thVortex-pre.pkl'
else:
    outfile = 'OPZ-extract-4thVortex.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(dats,f)

#%%
# date = datetime(2020,2,23,6)
# dat = ECMWF('OPZ',date)
# dat._get_var('VO')
# #datr = dat.shift2west(-179)
# dat2 = dat.extract(varss=['VO'],latRange=[-80,-40],lonRange=(270,320))
# dat2.show('VO',38)
