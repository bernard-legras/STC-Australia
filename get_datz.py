#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory from  January 2020

Created on Sat Feb  8 12:27:14 2020
Version that calculates and extracts the geopotential

See get_traject for more comments

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import pickle,gzip
#import matplotlib.pyplot as plt

# predates is to add new dates at the begiinning of the file
# it generates a special file which is then merged with
# the existing one
predates = False
update = True
day1 = datetime(2020,3,29,6)
day2 = day1 + timedelta(days=5)
if predates:
    day1 = datetime(2020,1,4,6)
    day2 = day1 + timedelta(days=3)
    update = False

if update==False:
# Initial rune
    datz = {}
    i = 0
else:
    # Continuation
    with gzip.open('OPZ-Z-2.pkl','rb') as f:
        datz = pickle.load(f)
    #del datz[117]
    #del datz[116]
    i = len(datz)+130
    print('datz length ',len(datz))

date = day1

while date < day2:
    print(date)
    dat = ECMWF('OPZ',date)
    dat._get_var('T')
    dat._mkp()
    dat._mkz()
    if date >= datetime(2020,3,29):
        datr = dat.shift2west(-179)
        datz[i] = datr.extract(varss=['Z'],latRange=(-45,-15),lonRange=(30,150))
    elif date >= datetime(2020,3,19):
        datz[i] = dat.extract(varss=['Z'],latRange=(-45,-15),lonRange=(80,200))
    elif date >= datetime(2020,3,9):
        datr = dat.shift2west(-179)
        datz[i] = datr.extract(varss=['Z'],latRange=(-45,-15),lonRange=(-179,-59))
    elif date >= datetime(2020,3,4):
        datr = dat.shift2west(-179)
        datz[i] = datr.extract(varss=['Z'],latRange=(-45,-15),lonRange=(-100,20))
    elif date >= datetime(2020,2,26):
        datz[i] = dat.extract(varss=['Z'],latRange=(-58,-28),lonRange=(0,120))
    elif date >= datetime(2020,2,13):
        datz[i] = dat.extract(varss=['Z'],latRange=(-58,-28),lonRange=(90,210))
    else:
        datz[i] = dat.extract(varss=['Z'],latRange=(-65,-35),lonRange=(180,300))
    dat.close()
    date += timedelta(hours=12)
    i += 1

if predates:
    outfile = 'OPZ-Z-pre.pkl'
else:
    outfile = 'OPZ-Z-2.pkl'
with gzip.open(outfile,'wb') as f:
    pickle.dump(datz,f)