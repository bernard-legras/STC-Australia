#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Profile extractor from IFS and CAMS (ozone, temperature, moisture)

Created on Thu Feb 20 16:22:03 2020

@author: Bernard Legras
"""
import pickle
from datetime import datetime,timedelta
import numpy as np
from ECMWF_N import ECMWF
import constants as cst

# liste of dates 
# 11 Jan, 16 Jan, 23 Jan, 31 Jan, 9 Feb

# get vortex positions at 6 and 18 every day
trac = pickle.load(open('Vortex-track.pkl','rb'))

dates  = [datetime(2020,1,11,12),datetime(2020,1,16,12),datetime(2020,1,23,12),
          datetime(2020,1,31,12),datetime(2020,2,9,12)]

profs = {}

listvar = ['T','P','Z','PT','O3','Q','PT']
for date in dates:
    profs[date] = {}
   
    for var in listvar:
        profs[date][var] = np.zeros(137)
    dateIFS = [date-timedelta(hours=6),date+timedelta(hours=6)]
    pos  = [np.where([trac['dates'][i] == dateIFS[0] for i in range(len(trac['dates']))])[0][0],
            np.where([trac['dates'][i] == dateIFS[1] for i in range(len(trac['dates']))])[0][0]]
    pt0 = trac['T'][pos[0]] * (trac['p'][pos[0]]/cst.p0)**(-cst.kappa)
    pt1 = trac['T'][pos[1]] * (trac['p'][pos[1]]/cst.p0)**(-cst.kappa)
    
    profs[date]['pos'] = {'lon':0.5*(trac['lons'][pos[0]]+trac['lons'][pos[1]]),
                          'lat':0.5*(trac['lons'][pos[0]]+trac['lons'][pos[1]]),
                          'T':0.5*(trac['T'][pos[0]]+trac['T'][pos[1]]),
                          'p':0.5*(trac['p'][pos[0]]+trac['p'][pos[1]]),
                          'z':0.5*(trac['p'][pos[0]]+trac['p'][pos[1]]),
                          'pt':0.5*(pt0+pt1),
                          'o3':0.5*(trac['o3'][pos[0]]+trac['o3'][pos[1]])}                                                       
    for i in [0,1]:
        print(dateIFS[i])
        dat = ECMWF('OPZ',dateIFS[i])
        dat._get_var('T')
        dat._mkp()
        dat._mkz()
        dat._mkthet()
        dat._get_var('O3')
        dat._get_var('Q')
        ix = np.where(dat.attr['lons']>=trac['lons'][pos[i]])[0][0]
        jy = np.where(dat.attr['lats']>=trac['lats'][pos[i]])[0][0]
        for var in listvar:
            profs[date][var] += dat.var[var][:,jy,ix]
    for var in listvar:
        profs[date][var] *= 0.5
    dat.close()
        
#%% Plot profiles of T, O3, Q
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,5,sharey=True)
i = 0
for date in dates:
    ax[i].plot(profs[date]['T'],profs[date]['PT'])
    ax[i].plot(profs[date]['pos']['T'],profs[date]['pos']['pt'],marker='D')
    ax[i].set_xlabel('Temperature (K)')
    if i==0: ax[i].set_ylabel('Potential temperature (K)')
    ax[i].set_title('Temperature')
    ax[i].set_ylim(280,700)
    ax[i].set_xlim(200,300)
    i += 1
plt.show()
#%%
fig, ax = plt.subplots(1,5,sharey=True)
i = 0
for date in dates:
    ax[i].plot(1.e6*profs[date]['O3'],profs[date]['PT'])
    ax[i].plot(1.e6*profs[date]['pos']['o3'],profs[date]['pos']['pt'],marker='D')
    ax[i].set_xlabel('O3 (mg/kg)')
    if i==0: ax[i].set_ylabel('Potential temperature (K)')
    ax[i].set_title('Ozone')
    ax[i].set_ylim(280,700)
    #ax[i].set_xlim(200,300)
    i += 1
plt.show()
#%%
fig, ax = plt.subplots(1,5,sharey=True)
i = 0
for date in dates:
    ax[i].plot(1.e6*profs[date]['Q'],profs[date]['PT'])
    ax[i].plot(4.5,profs[date]['pos']['pt'],marker='D')
    ax[i].set_xlabel('Q (mg/kg)')
    if i==0: ax[i].set_ylabel('Potential temperature (K)')
    ax[i].set_title('Moisture')
    ax[i].set_ylim(400,700)
    ax[i].set_xlim(0,5)
    i += 1
plt.show()

#%%
pickle.dump(profs,open('profs.pkl','wb'))
    
        


