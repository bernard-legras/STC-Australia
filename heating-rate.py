#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Evaluate the ascent rate by doing a fit to the vortex track

Created on Sat Apr 11 02:35:01 2020

@author: Bernard Legras
"""
import numpy as np
#from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
#from os.path import join

trac=pickle.load(open('Vortex-track.pkl','rb'))
days = 0.5*np.arange(len(trac['pt']))
s,a=np.polyfit(days,trac['pt'],1)
fig,[ax0,ax1]=plt.subplots(2,1,sharex=True,sharey=False,figsize=(4,6))
ax0.plot(days,trac['pt'],days,a+s*days)
#ax0.set_xlabel('days')
ax0.set_ylabel('potential temperature')
ax0.set_title('Vortex ascent in potential temperature. Slope = 5.94 K/day')
ax1.plot(days,s*np.array(trac['T'])/np.array(trac['pt']))
ax1.set_xlabel('days')
ax1.set_ylabel('heating rate (K/day)')
ax1.set_title('Evolution of heating rate')
plt.show()