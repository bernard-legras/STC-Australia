#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory from 7 January 2020

Created on Sat Feb  8 12:27:14 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import pickle,gzip
#import matplotlib.pyplot as plt

# predates is used here to replace the 8 first days by a a box shifted northward
predates = False
update = True
day1 = datetime(2020,1,27,6)
day2 = day1 + timedelta(days=12)
if predates:
    day1 = datetime(2020,1,5,6)
    day2 = day1 + timedelta(days=10)
    update = False

if update==False:
# Initial rune
    dats = {}
    i = 0
else:
    # Continuation
    with gzip.open('OPZ-extract-2ndVortex.pkl','rb') as f:
        dats = pickle.load(f)
    #for i in range(44,54):
    #    del dats[i]
    i = len(dats)
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
    dat._mkp()
    dat._mkz()
    if date >= datetime(2020,2,3):
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-50,-20),lonRange=(130,250))
    elif date >= datetime(2020,2,1):
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-60,-30),lonRange=(160,280))
    elif  date >= datetime(2020,1,27):
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-60,-30),lonRange=(80,200))
    elif date >= datetime(2020,1,22):
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-60,-30),lonRange=(0,120))
    elif date >= datetime(2020,1,18):
        datr = dat.shift2west(-179)
        dats[i] = datr.extract(varss=['T','VO','O3','Z'],latRange=(-60,-30),lonRange=(-90,30))
    elif date >= datetime(2020,1,15):
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-60,-30),lonRange=(210,330))
    else:
        dats[i] = dat.extract(varss=['T','VO','O3','Z'],latRange=(-50,-20),lonRange=(150,270))
    dat.close()
    date += timedelta(hours=12)
    i += 1

outfile = 'OPZ-extract-2ndVortex.pkl'
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
