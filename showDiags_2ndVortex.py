#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
showDiags version for the 2nd vortex
Includes the generation of the multitime composite vorticity chart.
Excludes composite.

Created on Sat Feb  8 12:27:14 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import constants as cst
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D # yes it is used
#from matplotlib import cm
#from matplotlib.text import TextPath
import gzip,pickle
import socket
#import deepdish as dd
from os.path import join
#from PIL import Image, ImageDraw, ImageFont
#from os import chdir


def tracker(dats,lon,lat,upper=35,lower=70,idx=6,jdy=4):
    # extract volume surrounding the target point
    jy = np.where(dats.attr['lats']>=lat)[0][0]
    ix = np.where(dats.attr['lons']>=lon)[0][0]
    jmin = max(jy-jdy,0)
    imin = max(ix-idx,0)
    jmax = min(jy+jdy+1,dats.nlat)
    imax = min(ix+idx+1,dats.nlon)
    sample = dats.var['VO'][upper:lower,jmin:jmax,imin:imax]
    # find the 3d index of the max vorticity in the sample cube
    aa = np.unravel_index(np.argmax(sample, axis=None), sample.shape)
    return([dats.attr['zscale'][upper+aa[0]],
            dats.attr['lats'][jmin+aa[1]],
            dats.attr['lons'][imin+aa[2]],
            dats.var['VO'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['Z'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['T'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['O3'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.attr['pscale'][upper+aa[0]],
            imin+aa[2],jmin+aa[1],upper+aa[0]])

if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-Australia'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-Australia'

#%%
with gzip.open('OPZ-extract-2ndVortex.pkl','rb') as f:
    dats = pickle.load(f)
# delete the non relevant final 4 states
for i in range(64,68):
    del dats[i]
print(len(dats))
with gzip.open('figs/_Silvia/AI_2v.pkl','rb') as f:
    AIT = pickle.load(f)
AITdates=AIT[:,1]


#%%
figsav = True
figargs = dict(bbox_inches='tight',dpi=300)
# tracking of the 3D position of the vortex
trac={'dates':[],'lons':[],'lats':[],'alts':[],'vo':[],'z':[],'T':[],'p':[],
      'pt':[],'o3':[],'ix':[],'jy':[],'kz':[]}
# initial position
#lon = 245
#lat = -54
#date = datetime(2020,1,7,6)
lon = 180
lat = -30
date = datetime(2020,1,5,6)
for i in range(len(dats)):
    print(i)
    lower = 66
    upper= 35
    idx = 6
    jdy = 4
    # track with one iteration
    try:
        lon = 2*lon - trac['lons'][-2]
        lat = 2*lat - trac['lats'][-2]
    except: pass
    #if i==12:
    #    lat = -36
    if i ==3:
        lon=196
    if i ==6:
        lat = -32; lon = 224
    if i ==12:
        lat = -36; lon = 224
    #if i ==14:
    #    lat = -34; lon = 220
    if i in [20,21]:
        upper = 55
    if i == 26:
        lon -= 360
    if i == 27:
        lon += 360
    if i in [48,51,56,57,58,59]:
        lower = 55
        upper = 46
    # help
    [alt,lat,lon,vo,z,T,o3,p,ix,jy,kz] = tracker(dats[i],lon,lat,lower=lower,upper=upper,idx=idx,jdy=jdy)
    #[alt,lat,lon,vo,z,T,o3,p,ix,jy,kz] = tracker(dats[i],lon,lat,lower=lower,upper=upper,idx=idx,jdy=jdy)
    trac['dates'].append(date)
    trac['alts'].append(alt)
    trac['lons'].append(lon)
    trac['lats'].append(lat)
    trac['vo'].append(vo)
    trac['z'].append(z/1000)
    trac['T'].append(T)
    trac['o3'].append(o3)
    trac['p'].append(p)
    trac['pt'].append(T*(cst.p0/p)**cst.kappa)
    trac['kz'].append(kz)
    trac['jy'].append(jy)
    trac['ix'].append(ix)
    date += timedelta(hours=12)

#%% print the positions as a function of time
for i in range(len(trac['p'])):
    # kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    print(i,trac['dates'][i],trac['lons'][i]%360,trac['lats'][i],'{:2.1f} {:2.1f}'.format(trac['z'][i],trac['vo'][i]*1.e5),trac['kz'][i])
pickle.dump(trac,open('Vortex-track_2ndVortex.pkl','wb'))

#%% Localisation of the vortex
fig,((ax0,ax1),(ax2,ax3))=plt.subplots(2,2,figsize=(8,8),sharex=True)
ax0.plot(trac['dates'],trac['lats'],linewidth=4)
ax1.plot(trac['dates'],np.array(trac['lons']) %360,linewidth=4)
ax2.plot(trac['dates'],trac['alts'],linewidth=4)
ax3.plot(trac['dates'],1.e5*np.array(trac['vo']),linewidth=4)
ax0.set_ylabel('Latitude')
ax1.set_ylabel('Longitude')
ax2.set_ylabel('Altitude (km)')
ax3.set_ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_2ndVortex.png'),**figargs)
plt.show()

#%% Make a two panel figure with the altitude and the vorticity
lcut = -11
fig,(ax2,ax3) = plt.subplots(1,2,figsize=(9,3))
ax2.plot(trac['dates'][1:lcut],trac['alts'][1:lcut],linewidth=4)
ax3.plot(trac['dates'][1:lcut],1.e5*np.array(trac['vo'][1:lcut]),linewidth=4)
ax2.set_ylabel('Altitude (km)')
ax3.set_ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_2ndVortex_2panel.png'),**figargs)
    plt.savefig(join('figs','VortexMotion_2ndVortex_2panel.pdf'),**figargs)
plt.show()
#%% Variant where each panel is a figure
lcut = -11
fig = plt.figure(figsize=(5,3))
plt.plot(trac['dates'][1:lcut],trac['alts'][1:lcut],linewidth=4)
plt.ylabel('Altitude (km)')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_2ndVortex_altitude.png'),**figargs)
    plt.savefig(join('figs','VortexMotion_2ndVortex_altitude.pdf'),**figargs)
plt.show()
fig = plt.figure(figsize=(5,3))
plt.plot(trac['dates'][1:lcut],1.e5*np.array(trac['vo'][1:lcut]),linewidth=4)
plt.ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$)')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_2ndVortex_vorticity.png'),**figargs)
    plt.savefig(join('figs','VortexMotion_2ndVortex_vorticity.pdf'),**figargs)
plt.show()

#%%
ax=dats[0].show('VO',62,show=False)
ax.plot(trac['lons'][0]-180,trac['lats'][0],'+',markersize=51,mew=3)
ax.plot(30,-30,'r',marker='x',markersize=51,mew=3)
plt.show()

#%% Plot of some selected days with the location of the IFS vortex and the TROPOMI patch
First_day_2plot = 28
Last_day_2plot = len(dats)

# beginning of the exploration (from get_traject)
for i in np.arange(First_day_2plot,Last_day_2plot):
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    ax = dats[i].show('VO',kz,clim=(-5.e-5,10.e-5),
        txt='vorticity lev {:d} alt {:2.1f} km  {}  (s**-1)'.format(kz,trac['alts'][i],trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        show=False)
    if dats[i].attr['lons'].max() > 180: cm = 180
    else: cm = 0
    ax.plot(trac['lons'][i]-cm,trac['lats'][i],'red',marker='+',markersize=51,mew=4)
    j=0
    try:
        while (AITdates[j] + timedelta(days=1)) < dats[i].date:
            j +=1
        print(i,dats[i].date,AITdates[j])
        ll = AIT[j,2]%360
        # test interval where box straddles Greenwich
        if 26 <= i <= 31: ll -= 360
        ax.plot(ll-cm,AIT[j,3],'lime',marker='x',markersize=51,mew=4)

    except IndexError:
        print('no AIT data')
    if figsav:
        plt.savefig(join('figs','VOmaps_2ndVortex','VO_'+str(i)+'.png'),**figargs)
    plt.show()
    date += timedelta(days=1)

#%% Show trajectory on a background that spans the latitude band
#convert lon from (0,360) to (-180,180)
gl = lambda x: (x+180)%360 - 180
dat = ECMWF('OPZ',datetime(2020,1,16,18))
dat._get_var('VO')
datr = dat.shift2west(-179)
datb = datr.extract(varss=['VO'],latRange=(-65,-25))
ax = datb.show('VO',62,show=False,aspect=1.4,scale=1.e5,clim=(-5,10))
#ax.plot(gl(np.array(trac['lons'][:-6])),trac['lats'][:-6],'+',mec='red',ms=21,mew=3)
ax.plot(gl(np.array(trac['lons'][1:lcut])),trac['lats'][1:lcut],linewidth=5,color='yellow')
#ax.plot((mean_track['centroid_lon'] % 360) - cm,mean_track['centroid_lat'],'o',markersize=8,color='r')
#ax.plot((trackSV.T[2] % 360)-cm,trackSV.T[3],'P',markersize=8,color='m',alpha=0.7)
ax.plot(gl(np.array(AIT[:,2])%360),AIT[:,3],'P',color='m',mew=0.3,markersize=17,alpha=0.7)
if figsav:
    plt.savefig(join('figs','multivortex2_horizontal_VO.png'),**figargs)
plt.show()

#%% ======= Show trajectory on a background that spans the latitude band ==========
#convert lon from (0,360) to (-180,180)
gl = lambda x: (x+180)%360 - 180
i = 33 # 24 Jan
kz = trac['kz'][i]
events = [4,14,23,26,29,38,42,46,50]
wx = 9
wy = 7
wz = 8
lcut = -11
ysouth = -65
ynorth = -25
xwest = -179
xeast = 180
boox = {}
dat = ECMWF('OPZ',trac['dates'][i])
dat._get_var('VO')
datr = dat.shift2west(-179)
datp = datr.extract(varss=['VO'],latRange=(ysouth,ynorth))

for i1 in events:
    kz1 = trac['kz'][i1]
    #jy1 = np.where(dats[i].attr['lats']>=trac['lats'][i1])[0][0]
    date1 = trac['dates'][i1]
    print(date1)
    dat1 = ECMWF('OPZ',date1)
    dat1._get_var('VO')
    datr1 = dat1.shift2west(-179)
    datp1 = datr1.extract(varss=['VO',],latRange=(ysouth,ynorth))
    lon1 = gl(trac['lons'][i1]%360)
    lat1 = trac['lats'][i1]
    ix1 = int(lon1-datp1.attr['lons'][0])
    jy1 = int(lat1-datp1.attr['lats'][0])
    print('ix1,jy1',ix1,jy1)
    # vertical projection of the submap around the vortex from kz1 to kz lavel
    datp.var['VO'][kz,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1] = \
        datp1.var['VO'][kz1,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1]
    boox[i1] = np.array([[lon1-wx,lon1-wx,lon1+wx,lon1+wx,lon1-wx],
                         [lat1-wy,lat1+wy,lat1+wy,lat1-wy,lat1-wy]])
#%%graphics section
ax=datp.show('VO',kz,show=False,clim=(-2,10),scale=1.e5,aspect=1.414,figsize=(12,4))
bbox = dict(boxstyle='round4',fc='w')
ax.text(gl(trac['lons'][i]),trac['lats'][i]+wy+5,dat.date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(boox[i1][0],boox[i1][1],'w',linewidth=2)
    date1 = trac['dates'][i1]
    up = 5
    if i1 <= 14: up = -19
    ax.text(gl(trac['lons'][i1]%360),trac['lats'][i1]+wy+up,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.plot(gl(np.array(trac['lons'][1:lcut])%360),trac['lats'][1:lcut],linewidth=5,color='yellow')
ax.plot(gl(np.array(AIT[2:,2])%360),AIT[2:,3],'P',color='m',mew=0.3,markersize=13,alpha=0.7)
ax.set_title(u'vorticity (10$^{-5}$ s$^{-1}$)',fontsize=18)
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
if figsav:
    plt.savefig(join('figs','multivortex2_horizontal_VO.png'),**figargs)
    plt.savefig(join('figs','multivortex2_horizontal_VO.pdf'),**figargs)
plt.show()