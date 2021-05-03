#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
showDiags version for the 4th vortex
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
with gzip.open('OPZ-extract-4thVortex-post.pkl','rb') as f:
    dats2 = pickle.load(f)
print(len(dats2))
with gzip.open('OPZ-extract-4thVortex-pre.pkl','rb') as f:
    dats = pickle.load(f)
print(len(dats))

for i in range(1,18):
    del dats2[i]
for i in range(18,82):
    dats[46+i-18] = dats2[i]
print(len(dats))

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
lon = 270
lat = -80
lon = 360 - 163
lat = -62
date = datetime(2020,1,6,6)
for i in range(0,102):
    print(i)
    lower = 63
    upper= 35
    idx = 6
    jdy = 4
    # track with one iteration
    try:
        lon = 2*lon - trac['lons'][-2]
        lat = 2*lat - trac['lats'][-2]
    except: pass
    if i == 12:
        lon = 225; idx=4; upper=51; lower =53
    if i == 98:
        lon -= 360
    if i == 99:
        lon += 360
    if i >= 46:
        lower = 50
        upper= 35
    if i >= 102:
        lower = 40
        idx = 2; idy = 2
    if i == 26: upper = 45
    #if i == 4: lower = 55
    #if i == 5: upper = 45
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
#pickle.dump(trac,open('Vortex-track_4thVortex.pkl','wb'))
#%% reload IFS trac dand TROPOMI TRAC
trac = pickle.load(open('Vortex-track_4thVortex.pkl','rb'))
for i in range(len(trac['p'])):
    print(i,trac['dates'][i],trac['lons'][i]%360,trac['lats'][i],'{:2.1f} {:2.1f}'.format(trac['z'][i],trac['vo'][i]*1.e5),trac['kz'][i])
with gzip.open('figs/_Silvia/AI_4v.pkl','rb') as f:
    AIT = pickle.load(f)
AITdates=AIT[:,1]
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
    plt.savefig(join('figs','VortexMotion_4thVortex.png'),**figargs)
plt.show()
#%% Make a two panel figure with the altitude and the vorticity
fig,(ax2,ax3) = plt.subplots(1,2,figsize=(9,3))
ax2.plot(trac['dates'],trac['alts'],linewidth=4)
ax3.plot(trac['dates'],1.e5*np.array(trac['vo']),linewidth=4)
ax2.set_ylabel('Altitude (km)')
ax3.set_ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_4thVortex_2panel.png'),**figargs)
    plt.savefig(join('figs','VortexMotion_4thVortex_2panel.pdf'),**figargs)
plt.show()
#%% Variant where each panel is a figure
lcut = -1
fig = plt.figure(figsize=(6,3))
plt.plot(trac['dates'][6:lcut],trac['alts'][6:lcut],linewidth=4)
plt.ylabel('Altitude (km)')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_4thVortex_altitude.png'),**figargs)
    plt.savefig(join('figs','VortexMotion_4thVortex_altitude.pdf'),**figargs)
plt.show()
fig = plt.figure(figsize=(6,3))
plt.plot(trac['dates'][6:lcut],1.e5*np.array(trac['vo'][6:lcut]),linewidth=4)
plt.ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$)')
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotion_4thVortex_vorticity.png'),**figargs)
    plt.savefig(join('figs','VortexMotion_4thVortex_vorticity.pdf'),**figargs)
plt.show()
#%% Show trajectory on a background that spans the latitude band
#convert lon from (0,360) to (-180,180)
gl = lambda x: (x+180)%360 - 180
dat = ECMWF('OPZ',datetime(2020,2,7,6))
dat._get_var('VO')
datr = dat.shift2west(-179)
datb = datr.extract(varss=['VO'],latRange=(-80,-40),lonRange=(-90,-10))
ax = datb.show('VO',39,show=False,aspect=1.4,scale=1.e5,clim=(-2,5))
#ax.plot(gl(np.array(trac['lons'][:-6])),trac['lats'][:-6],'+',mec='red',ms=21,mew=3)
ax.plot(gl(np.array(trac['lons'])%360),trac['lats'],linewidth=5,color='yellow')
#ax.plot((mean_track['centroid_lon'] % 360) - cm,mean_track['centroid_lat'],'o',markersize=8,color='r')
#ax.plot((trackSV.T[2] % 360)-cm,trackSV.T[3],'P',markersize=8,color='m',alpha=0.7)
if figsav:
    plt.savefig(join('figs','multivortex4_horizontal_VO.png'),**figargs)
plt.show()

#%% ========= Show trajectory on a background that spans the latitude band =============
#convert lon from (0,360) to (-180,180)
gl = lambda x: (x+180)%360 - 180
i = 36+28
kz = trac['kz'][i]
events = [6,12,19,28,16+28,46+28,68+28]
wx = 8
wy = 6
wz = 8
xwest = -179
xeast = -20
ysouth = -84
ynorth = -40
boox = {}
dat = ECMWF('OPZ',datetime(2020,2,7,6))
dat._get_var('VO')
datr = dat.shift2west(-179)
datp = datr.extract(varss=['VO'],latRange=(ysouth,ynorth),lonRange=(xwest,xeast))

for i1 in events:
    wx = 8
    if i1 == 12: wx = 12
    if i1 == 19: wx = 8
    kz1 = trac['kz'][i1]
    #jy1 = np.where(dats[i].attr['lats']>=trac['lats'][i1])[0][0]
    date1 = trac['dates'][i1]
    print(date1)
    dat1 = ECMWF('OPZ',date1)
    dat1._get_var('VO')
    datr1 = dat1.shift2west(-179)
    datp1 = datr1.extract(varss=['VO',],latRange=(ysouth,ynorth),lonRange=(xwest,xeast))
    lon1 = gl(trac['lons'][i1]%360)
    lat1 = trac['lats'][i1]
    ix1 = int(lon1-datp1.attr['lons'][0])
    jy1 = int(lat1-datp1.attr['lats'][0])
    print('ix1,jy1',ix1,jy1)
    # vertical projection of the submap around the vortex from kz1 to kz lavel
    jmin = max(0,jy1-wy)
    datp.var['VO'][kz,jmin:jy1+wy+1,ix1-wx:ix1+wx+1] = \
        datp1.var['VO'][kz1,jmin:jy1+wy+1,ix1-wx:ix1+wx+1]
    boox[i1] = np.array([[lon1-wx,lon1-wx,lon1+wx,lon1+wx,lon1-wx],
                         [max(lat1-wy,ysouth),lat1+wy,lat1+wy,max(lat1-wy,ysouth),max(lat1-wy,ysouth)]])
#%%graphics section
ax=datp.show('VO',kz,show=False,clim=(-2,5),scale=1.e5,aspect=2,figsize=(7,5))
bbox = dict(boxstyle='round4',fc='w')
ax.text(gl(trac['lons'][i]),trac['lats'][i]+wy-12,dat.date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(boox[i1][0],boox[i1][1],'w',linewidth=2)
    date1 = trac['dates'][i1]
    ax.text(gl(trac['lons'][i1]%360),trac['lats'][i1]+wy+2,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.plot(gl(np.array(trac['lons'])%360),trac['lats'],linewidth=5,color='yellow')
ax.plot(gl(np.array(AIT[2:,2])%360),AIT[2:,3],'P',color='m',mew=0.3,markersize=13,alpha=0.5)
ax.set_xlim(xwest,xeast)
ax.set_title(u'vorticity (10$^{-5}$ s$^{-1}$)',fontsize=18)
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
if figsav:
    plt.savefig(join('figs','multivortex4_horizontal_VO.png'),**figargs)
    plt.savefig(join('figs','multivortex4_horizontal_VO.pdf'),**figargs)
plt.show()

#%% Plot of some selected days with the location of the IFS vortex and the TROPOMI patch
First_day_2plot = 70
Last_day_2plot = len(trac['lons'])
figsav = True

# beginning of the exploration (from get_traject)
for i in np.arange(First_day_2plot,Last_day_2plot,2):
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    ax = dats[i].show('VO',kz,clim=(-2,5),scale=1.e5,
        txt='vorticity lev {:d} alt {:2.1f} km  {}  (s**-1)'.format(kz,trac['alts'][i],trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        show=False)
    #if dats[i].attr['lons'].max() > 180: cm = 180
    #else: cm = 0
    cm = 180
    if i >=70: cm=0
    ax.plot(trac['lons'][i]-cm,trac['lats'][i],'red',marker='+',markersize=51,mew=4)
    j=0
    if figsav:
        plt.savefig(join('figs','VOmaps_4thVortex','VO_'+str(i)+'.png'),**figargs)
    plt.show()
    date += timedelta(days=1)
