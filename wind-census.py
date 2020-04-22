#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes a couple of plots of the wind across the vortex
Exploratory script.

Created on Fri Apr 10 02:11:58 2020

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
#from os.path import join
from ECMWF_N import ECMWF

#%%

dat1 = ECMWF('OPZ',datetime(2020,1,21,6))
dat2 = ECMWF('OPZ',datetime(2020,2,10,6))

dat1._get_var('U')
dat1._get_var('V')
dat1._get_var('VO')
dat2._get_var('U')
dat2._get_var('V')
dat2._get_var('VO')

#%%
dats1 = dat1.extract(varss=['U','V','VO'],latRange=(-70,-50),lonRange=(265,285))
dats2 = dat2.extract(varss=['U','V','VO'],latRange=(-61,-41),lonRange=(205,225))
Ua1 = dats1.var['U']- (np.mean(dats1.var['U'],axis=(1,2)))[:,np.newaxis,np.newaxis]
Va1 = dats1.var['V']- (np.mean(dats1.var['V'],axis=(1,2)))[:,np.newaxis,np.newaxis]
Ua2 = dats2.var['U']- (np.mean(dats2.var['U'],axis=(1,2)))[:,np.newaxis,np.newaxis]
Va2 = dats2.var['V']- (np.mean(dats2.var['V'],axis=(1,2)))[:,np.newaxis,np.newaxis]
#dats1.var['UM'] = np.sqrt(dats1.var['U']**2 + dats1.var['V']**2)
#dats2.var['UM'] = np.sqrt(dats2.var['U']**2 + dats2.var['V']**2)
dats1.var['UMa'] = np.sqrt(Ua1**2 + Va1**2)
dats2.var['UMa'] = np.sqrt(Ua2**2 + Va2**2)
dats1.show('UMa',47)
dats2.show('UMa',38)

#%%
trac = pickle.load(open('Vortex-track.pkl','rb'))
trac['UAnoMax'] = []
trac['Umax'] = []
for i in range(0,len(trac['dates'])):
    dat = ECMWF('OPZ',trac['dates'][i])
    dat._get_var('U')
    dat._get_var('V')
    lon = trac['lons'][i]%360
    if (lon <= 50) | (lon >=310):
        datr = dat.shift2west(-179)
        if lon > 180: lon -= 360
        dats = datr.extract(varss=['U','V'],latRange=(trac['lats'][i]-10,trac['lats'][i]+10),
                       lonRange=(lon-10,lon+10))
        del datr
    else:
        dats = dat.extract(varss=['U','V'],latRange=(trac['lats'][i]-10,trac['lats'][i]+10),
                       lonRange=(lon-10,lon+10))
    del dat
    dats.var['UM'] = np.sqrt(dats.var['U']**2 + dats.var['V']**2)
    Ua = dats.var['U']- (np.mean(dats.var['U'],axis=(1,2)))[:,np.newaxis,np.newaxis]
    Va = dats.var['V']- (np.mean(dats.var['V'],axis=(1,2)))[:,np.newaxis,np.newaxis]
    dats.var['UMa'] = np.sqrt(Ua**2 + Va**2)
    dats.show('UMa',trac['kz'][i],txt=trac['dates'][i].strftime('%Y %m %d')+'  L'+str(trac['kz'][i]))
    trac['Umax'].append(dats.var['UM'][trac['kz'][i],...].max())
    trac['UAnoMax'].append(dats.var['UMa'][trac['kz'][i],...].max())
    print(i,trac['dates'][i],
          '{:2.1f} {:d} {:2.1f} {:2.1f}'.format(trac['vo'][i]*1.e5,trac['kz'][i],trac['Umax'][i],trac['UAnoMax'][i]))
#%%
pickle.dump(trac,open('Vortex-track-withwind.pkl','wb'))
#%% Diagnostic of the max velocity during the most active phase
# this active phase is the period during which the wind anomaly forms a ring around the vortex
trac = pickle.load(open('Vortex-track-withwind.pkl','rb'))
plt.plot(trac['Umax']);plt.plot(trac['UAnoMax']);plt.show()
plt.plot(trac['vo']);plt.show()
plt.plot(trac['Umax'][10:110]);plt.plot(trac['UAnoMax'][10:110]);plt.show()
plt.plot(trac['vo'][10:110]);plt.show()
print('wind',np.mean(trac['Umax'][10:110]),np.sqrt(np.var(trac['Umax'][10:110])))
print('anow',np.mean(trac['UAnoMax'][10:110]),np.sqrt(np.var(trac['UAnoMax'][10:110])))
# The result is 17.3 pm 0.8 m/s (0.68*std)