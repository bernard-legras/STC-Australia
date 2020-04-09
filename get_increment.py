#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on Sat Feb  8 23:09:08 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import pickle,gzip
#import matplotlib.pyplot as plt

#%%
date = datetime(2020,1,7,6)

dats = {}
i = 0
while date < datetime(2020,2,8,6):
    print(date)
    dat = ECMWF('OPZ',date)
    dat._get_var('T')
    dat._get_var('TD')
    dat._get_var('VOD')
    dat._get_var('O3D')
    dat._mkpscale()
    dat._mkzscale()
    dats[i] = dat.extract(varss=['TD','VOD','O3D'],latRange=(-65,-35),lonRange=(180,300))
    dat.close()
    date += timedelta(hours=12)
    i += 1    
#%%
with gzip.open('OPZ-increment.pkl','wb') as f:
    pickle.dump(dats,f)
#%% Continuation
with gzip.open('OPZ-increment.pkl','rb') as f:
    dats = pickle.load(f)
i = len(dats)
date = datetime(2020,2,8,18)
while date < datetime(2020,2,10,6):
    print(date)
    dat = ECMWF('OPZ',date)
    dat._get_var('T')
    dat._get_var('TD')
    dat._get_var('VOD')
    dat._get_var('O3D')
    dat._mkpscale()
    dat._mkzscale()
    dats[i] = dat.extract(varss=['TD','VOD','O3D'],latRange=(-65,-35),lonRange=(180,300))
    dat.close()
    date += timedelta(hours=12)
    i += 1