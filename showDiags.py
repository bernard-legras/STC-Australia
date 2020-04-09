#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory from 7 January 2020

Created on Sat Feb  8 12:27:14 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
#from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import constants as cst
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # yes it is used
from matplotlib import cm
from matplotlib.text import TextPath
import gzip,pickle
import socket
import deepdish as dd
from os.path import join
from PIL import Image, ImageDraw, ImageFont
from os import chdir


def tracker(dats,datz,lon,lat,upper=35,lower=65,idx=6,jdy=3):
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
            datz.var['Z'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['T'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['O3'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.attr['pscale'][upper+aa[0]],
            imin+aa[2],jmin+aa[1],upper+aa[0]])

if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-Australia'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-Australia'

#%%
with gzip.open('OPZ-extract-1.pkl','rb') as f:
    dats = pickle.load(f)
with gzip.open('OPZ-extract-2.pkl','rb') as f:
    dats.update(pickle.load(f))
with gzip.open('OPZ-Z-1.pkl','rb') as f:
    datz = pickle.load(f)
with gzip.open('OPZ-Z-2.pkl','rb') as f:
    datz.update(pickle.load(f))
print(len(dats),len(datz))

#%%
figsav = False
figargs = dict(bbox_inches='tight',dpi=300)
# tracking of the 3D position of the vortex
trac={'dates':[],'lons':[],'lats':[],'alts':[],'vo':[],'z':[],'T':[],'p':[],
      'pt':[],'o3':[],'ix':[],'jy':[],'kz':[]}
# initial position
#lon = 245
#lat = -54
#date = datetime(2020,1,7,6)
lon = 214
lat = -47
date = datetime(2020,1,4,6)
for i in range(len(dats)):
    print(i)
    lower = 65
    upper= 25
    idx = 6
    jdy = 2
    if i>46: lower = 50
    # track with one iteration
    try:
        lon = 2*lon - trac['lons'][-2]
        lat = 2*lat - trac['lats'][-2]
    except: pass
    # help
    if i == 127:
        lon = -47
        lat = -29
    if i == 128:
        lon = -55
        lat = -29
    if i == 130:
        lon = -71
        lat = -32
        lower = 40
    if i == 131:
        lon = -77
        lat = -32
        lower = 40
    if i >= 150:
        lon = lon % 360
    if i == 152:
        lon = 170
    if i in [157,158,161,163]:
        lower = 32
    if i == 165:
        idx = 3
        lat = -30
    if i in [168,169]:
        lon = 110
        lat = -28
    if i == 170:
        lon = 110
        lat = -28
    if i == 171:
        lon = 95
        lat = -28
    if i == 178:
        lon = 85
        lat = -30
    if i == 179:
        lon = 81
        lat = -26

    [alt,lat,lon,vo,z,T,o3,p,ix,jy,kz] = tracker(dats[i],datz[i],lon,lat,lower=lower,upper=upper,idx=idx,jdy=jdy)
    [alt,lat,lon,vo,z,T,o3,p,ix,jy,kz] = tracker(dats[i],datz[i],lon,lat,lower=lower,upper=upper,idx=idx,jdy=jdy)
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
for i in range(len(dats)):
    # kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    print(i,trac['dates'][i],trac['lons'][i],trac['lats'][i],'{:2.1f}'.format(trac['z'][i]),trac['kz'][i])
pickle.dump(trac,open('Vortex-track.pkl','wb'))

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
    plt.savefig(join('figs','VortexMotion.png'),**figargs)
plt.show()

#%% Localisation of the vortex with koobor locations
mean_track = pickle.load(open('Koobor-Ntrack-L1-mean.pkl','rb'))
fig,((ax0,ax1),(ax2,ax3))=plt.subplots(2,2,figsize=(10,10),sharex=True)
ax0.plot(trac['dates'],trac['lats'],linewidth=4)
ax0.errorbar(mean_track['dates'],mean_track['centroid_lat'],
            yerr=np.array([mean_track['centroid_lat']-mean_track['south_lat'],
                            mean_track['north_lat']-mean_track['centroid_lat']]),
            marker='s',alpha=0.5,ls='',elinewidth=2,ecolor='red',ms=8,mec='red',
            mfc=[1,.6,.6],lolims=True,uplims=True)
#ax0.plot(trac['dates'],trac['lats'],linewidth=4)
ax1.plot(trac['dates'][:122],np.array(trac['lons'][:122]) %360,'-b',linewidth=4)
ax1.plot(trac['dates'][122:],np.array(trac['lons'][122:]) %360,'-b',linewidth=4)
ax1.plot(mean_track['dates'],mean_track['centroid_lon'] % 360,
            marker='s',alpha=0.5,ls='',ms=8,mec='red',mfc=[1,.6,.6])
ax2.plot(trac['dates'],trac['alts'],linewidth=4)
ax2.errorbar(mean_track['dates'],mean_track['centroid_alt'],
            yerr=np.array([mean_track['centroid_alt']-mean_track['bot_alt'],
                            mean_track['top_alt']-mean_track['centroid_alt']]),
            marker='s',alpha=0.5,ls='',elinewidth=2,ecolor='red',ms=8,mec='red',
            mfc=[1,.6,.6],lolims=True,uplims=True)
ax3.plot(trac['dates'],1.e5*np.array(trac['vo']),linewidth=4)
ax0.set_ylabel('Latitude (°)',fontsize=16)
ax1.set_ylabel('Longitude (°)',fontsize=16)
ax2.set_ylabel('Altitude (km)',fontsize=16)
ax3.set_ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$)',fontsize=16)
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotionWKoobor.png'),**figargs)
    plt.savefig(join('figs','VortexMotionWKoobor.pdf'),**figargs)
plt.show()
#%% Localisation of the vortex with koobor locations + forecast
# version in four panels
# Latitude
lcut = -2
mean_track = pickle.load(open('Koobor-Ntrack-L1-mean.pkl','rb'))
fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))
fig = plt.figure(figsize=(5,4))
plt.errorbar(mean_track['dates'],mean_track['centroid_lat'],
            yerr=np.array([mean_track['centroid_lat']-mean_track['south_lat'],
                            mean_track['north_lat']-mean_track['centroid_lat']]),
            marker='s',ls='',elinewidth=2,ecolor='red',ms=8,mec='red',
            mfc=[1,.6,.6],alpha=1,lolims=True,uplims=True,zorder=0)
plt.plot(trac['dates'][:lcut],trac['lats'][:lcut],linewidth=4,alpha=0.8)
plt.ylabel('Latitude (°)',fontsize=16)
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotionWKoobor-latitude.png'),**figargs)
    plt.savefig(join('figs','VortexMotionWKoobor-latitude.pdf'),**figargs)
plt.show()
#%% Longitude
fig = plt.figure(figsize=(5,4))
plt.plot(trac['dates'][:122],np.array(trac['lons'][:122]) %360,'-b',linewidth=4)
plt.plot(mean_track['dates'],mean_track['centroid_lon'] % 360,
            marker='s',alpha=1,ls='',ms=8,mec='red',mfc=[1,.6,.6],zorder=1)
plt.plot(trac['dates'][122:lcut],np.array(trac['lons'][122:lcut]) %360,'-b',
         linewidth=4,alpha=0.7)
plt.ylabel('Longitude (°)',fontsize=16)
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotionWKoobor-longitude.png'),**figargs)
    plt.savefig(join('figs','VortexMotionWKoobor-longitude.pdf'),**figargs)
plt.show()
#%% Altitude
fig = plt.figure(figsize=(5,4))
plt.errorbar(mean_track['dates'],mean_track['centroid_alt'],
            yerr=np.array([mean_track['centroid_alt']-mean_track['bot_alt'],
                            mean_track['top_alt']-mean_track['centroid_alt']]),
            marker='s',alpha=1,ls='',elinewidth=2,ecolor='red',ms=8,mec='red',
            mfc=[1,.6,.6],lolims=True,uplims=True,zorder=0)
plt.plot(trac['dates'],trac['alts'],linewidth=4,alpha=0.8)
for date in fctrac:
    ns = fctrac[date]['survival']+1
    plt.plot(fctrac[date]['dates'][:ns],fctrac[date]['z'][:ns],'k',linewidth=2)
plt.ylabel('Altitude (km)',fontsize=16)
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotionWKoobor-altitude.png'),**figargs)
    plt.savefig(join('figs','VortexMotionWKoobor-altitude.pdf'),**figargs)
plt.show()
#%% Vorticty
fig = plt.figure(figsize=(5,4))
plt.plot(trac['dates'],1.e5*np.array(trac['vo']),linewidth=4)
for date in fctrac:
    ns = fctrac[date]['survival']+1
    plt.plot(fctrac[date]['dates'][:ns],1.e5*np.array(fctrac[date]['vo'][:ns]),'k',linewidth=2)
plt.ylabel(u'Vorticity (10$^{-5}$ s$^{-1}$)',fontsize=16)
fig.autofmt_xdate()
if figsav:
    plt.savefig(join('figs','VortexMotionWKoobor-vorticity.png'),**figargs)
    plt.savefig(join('figs','VortexMotionWKoobor-vorticity.pdf'),**figargs)
plt.show()

#%% Plot the horizontal displacement
plt.plot(trac['lons'],trac['lats'],linewidth=3)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Horizontal displacement between 07-01-2020 and 20-02-2020')
#mm = {6:'7 Jan',14:'11 Jan',22:'15 Jan',32:'20 Jan',42:'25 Jan',52:'30 Jan',
#      62:'4 Feb',78:'9 Feb',82:'14 Feb',92:'18 Feb'}
mm = {2:'5 Jan',14:'11 Jan',42:'25 Jan',62:'4 Feb',84:'15 Feb',104:'25 Feb',124:'6 Mar'}
for d in mm:
    path = TextPath((5,0),mm[d])
    plt.plot(trac['lons'][d],trac['lats'][d],marker='D',markersize=6,color='k')
    plt.plot(trac['lons'][d],trac['lats'][d],marker=path,markersize=100,color='red')
if figsav:
    plt.savefig(join('figs','trajectory_18Mar.png'),**figargs)
plt.show()
#%%
fig = plt.figure()
plt.plot(trac['dates'],trac['z'],linewidth=4)
fig.autofmt_xdate()
plt.ylabel('Altitude (km)')
plt.xlabel('Time')
plt.title('Altitude of the vortex as a function of time')
if figsav:
    plt.savefig(join('figs','ascent_18Mar.png'),**figargs)
plt.show()
fig = plt.figure()
plt.plot(trac['dates'],1.e5*np.array(trac['vo']),linewidth=4)
fig.autofmt_xdate()
plt.ylabel(u'Vorticity ($10^{-5} s^{-1}$)')
plt.xlabel('Time')
plt.title('Maximum vorticity as a function of time')
if figsav:
    plt.savefig(join('figs','vorticity_18Mar.png'),**figargs)
plt.show()

#%% Add Koobor location on the horizontal plot, from CALIOP and Sentinel 5

mean_track = pickle.load(open('Koobor-Ntrack-L1-mean.pkl','rb'))
fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))

tSpt = mean_track['dates']
latSpt = mean_track['centroid_lat']
lonSpt = mean_track['centroid_lon']
altSpt = mean_track['centroid_alt']
topSpt = mean_track['top_alt']
botSpt = mean_track['bot_alt']
# for i in range(len(trackSpt)):
#     if trackSpt[i][2] != 'c':
#         tSpt.append(trackSpt[i][4])
#         latSpt.append(trackSpt[i][5])
#         lonSpt.append(trackSpt[i][6])
#         altSpt.append(trackSpt[i][7])
# Read Silvia's pointing of Koobor in CALIOP data
#trackSB = dd.io.load(join('figs','_Silvia','CALIOP_plume.hdf5'))
# Read Silvia's pointing of the locations of aerosol blobs in Sentinel 5
trackSV = pickle.load(gzip.open(join('figs','_Silvia','AI_1v.pkl'),'rb'))
# Adding hours to dates
for i in range(len(trackSV)): trackSV[i][1]  += timedelta(hours = trackSV[i][5])
plt.plot(trac['lons'],trac['lats'],linewidth=3)
plt.plot(lonSpt,latSpt,marker='X',markersize=8,color='b')
#plt.plot(trackSB['lon'][3:],trackSB['lat'][3:],'d',markersize=8,color='m')
plt.plot(trackSV.T[2] % 360,trackSV.T[3],'o',markersize=8,color='g')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Horizontal displacement between 07-01-2020 and 15-03-2020')
mm = {0:['7 Jan',(0,-8)],8:['11 Jan',(-20,0)],16:['15 Jan',(8,0)],26:['20 Jan',(-40,0)],
      36:['25 Jan',(5,0)],46:['30 Jan',(5,0)],56:['4 Feb',(5,0)],66:['9 Feb',(-8,-9)],
      76:['14 Feb',(-8,-8)],86:['18 Feb',(8,8)]}
for d in mm:
    #path = TextPath(mm[d][1],mm[d][0])
    path = TextPath((1,-2),mm[d][0],size=2)
    plt.plot(trac['lons'][d],trac['lats'][d],marker='D',markersize=6,color='k')
    plt.plot(trac['lons'][d],trac['lats'][d],'',marker=path,markersize=100,color='red')
if figsav:
    plt.savefig(join('figs','trajectory_with_Spt.png'),**figargs)
plt.show()

#%% Add Koobor location on the ascent plot
fig,ax = plt.subplots(1,1)
ax.fill_between(tSpt,botSpt,topSpt,facecolor='lightgray')
ax.plot(trac['dates'],trac['z'],tSpt,altSpt,linewidth=4)
fig.autofmt_xdate()
ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Time')
ax.set_title('Altitude of the vortex and the cloud as a function of time')
if figsav:
    plt.savefig(join('figs','ascent_withSpt.png'),**figargs)

#%% Kompozit
# The composit is made from 44 (29Jan) to 65 (8 Feb)
# and from 62 (4 Feb) to 87 (16 Feb)
# and from 82 (14 Feb) to 103 (24 Feb)
# and from 42 (14 Jan) to 99 (22 Feb)
kompozit={}
jdy = 6
idx = 8
kdz = 8
istart = 42
iend = 100
kompozit['VO'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
kompozit['O3ano'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
kompozit['Tano'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
for i in range(istart,iend):
    print(i)
    ix = np.where(dats[i].attr['lons']>=trac['lons'][i])[0][0]
    jy = np.where(dats[i].attr['lats']>=trac['lats'][i])[0][0]
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    kompozit['VO'] += dats[i].var['VO'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1]
    kompozit['O3ano'] += dats[i].var['O3'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1] \
                       - dats[i].var['O3'][kz,jy,ix]
    kompozit['Tano'] += dats[i].var['T'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1] \
                      - dats[i].var['T'][kz,jy,ix]
ns = iend - istart
kompozit['VO'] /= ns
kompozit['O3ano'] /= ns
kompozit['Tano'] /= ns

#%%
# Plot of the kompozit as 3D plots (not so convincing)
# Composit vorticity
fig = plt.figure(figsize=(7,7))
ax = fig.gca(projection = '3d')
X1 = np.arange(-8,9)
Y1 = np.arange(-6,7)
X,Y = np.meshgrid(X1,Y1)
Z = np.zeros_like(X)
tt = X / np.max(X)
buf = kompozit['VO'] - np.min(kompozit['VO'])
buf /= np.max(buf)
for k in np.arange(12,3,-4):
    ax.plot_surface(X,Y,Z-k+8,rstride=1,cstride=1,facecolors = cm.jet(buf[k,...]))
ax.set_xlabel('longitude (degree)')
ax.set_ylabel('latitude (degree)')
ax.set_zlabel('altitude (level)')
plt.title('composite normalized vorticity')
if figsav:
    plt.savefig(join('figs','kompo3D_vorticity_14Jan-22Feb.png'),**figargs)
plt.show()
#%%
# Composit ozone
fig = plt.figure(figsize=(7,7))
ax = fig.gca(projection = '3d')
X1 = np.arange(-8,9)
Y1 = np.arange(-6,7)
X,Y = np.meshgrid(X1,Y1)
Z = np.zeros_like(X)
tt = X / np.max(X)
buf = kompozit['O3ano'] - np.min(kompozit['O3ano'])
buf /= np.max(buf)
for k in np.arange(12,3,-4):
    ax.plot_surface(X,Y,Z-k+8,rstride=1,cstride=1,facecolors = cm.jet(buf[k,...]))
ax.set_xlabel('longitude (degree)')
ax.set_ylabel('latitude (degree)')
ax.set_zlabel('altitude (level)')
plt.title('composite normalized ozone')
if figsav:
    plt.savefig(join('figs','kompo3D_ozone_14Jan-22Feb.png'),**figargs)
plt.show()
#%%
# Composit temperature
fig = plt.figure(figsize=(7,7))
ax = fig.gca(projection = '3d')
X1 = np.arange(-8,9)
Y1 = np.arange(-6,7)
X,Y = np.meshgrid(X1,Y1)
Z = np.zeros_like(X)
tt = X / np.max(X)
buf = kompozit['Tano'] - np.min(kompozit['Tano'])
buf /= np.max(buf)
for k in np.arange(12,3,-4):
    ax.plot_surface(X,Y,Z-k+8,rstride=1,cstride=1,facecolors = cm.jet(buf[k,...]))
ax.set_xlabel('longitude (degree)')
ax.set_ylabel('latitude (degree)')
ax.set_zlabel('altitude (level)')
plt.title('composite normalized temperature')
if figsav:
    plt.savefig(join('figs','kompo3D_temperature_14Jan-22Feb.png'),**figargs)
plt.show()

#%% make a 2D horizontal plot of the kompozit at central level and +-4
#ff,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,12))
plt.figure(figsize=(14,13))
imargs = dict(origin='lower',aspect=1,interpolation='nearest',cmap='jet',extent=(-8,8,-6,6))
i=1
for k in np.arange(12,3,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(1.e5*kompozit['VO'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('VO level '+str(8-k)+u' ($10^{-5}~s^{-1}$)')
    plt.xlabel('longitude (degree)')
    plt.ylabel('latitude (degree)')
    i += 1
for k in np.arange(12,3,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(1.e6*kompozit['O3ano'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('O3 ano level '+str(8-k)+u' ($10^{-6}~kg~kg^{-1}$)')
    plt.xlabel('longitude (degree)')
    plt.ylabel('latitude (degree)')
    i += 1
for k in np.arange(12,3,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(kompozit['Tano'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('T anomaly level '+str(8-k)+' (K)')
    plt.xlabel('longitude (degree)')
    plt.ylabel('latitude (degree)')
    i += 1
if figsav:
    plt.savefig(join('figs','kompo_VO3T_14Jan-22Feb.png'),**figargs)
plt.show()

#%% make a 2D vertical lon-alt section in the central longitude plane
fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(14,5),sharey=True)
fig.subplots_adjust(left=0.02, bottom=0.12, right=0.95, top=0.78, wspace=0.05)
central = 7
imargs = dict(aspect=1.4,interpolation='nearest',cmap='jet',extent=(-5.6,5.6,4,-4))
im1 = ax1.imshow(1.e5*kompozit['VO'][:,central,:],**imargs)
ax1.set_title('VO vertical section'+u' ($10^{-5}~s^{-1}$)',fontsize=18)
ax1.set_xlabel('longitude (100 km)',fontsize=16)
ax1.set_ylabel('altitude (km)',fontsize=16)
plt.colorbar(im1,ax=ax1)
im2 = ax2.imshow(1.e6*kompozit['O3ano'][:,central,:],**imargs)
ax2.set_title('O3 anomaly'+u' ($10^{-6}~kg~kg^{-1}$)',fontsize=18)
ax2.set_xlabel('longitude (100 km)',fontsize=16)
#ax2.set_ylabel('altitude (km)')
plt.colorbar(im2,ax=ax2)
im3 = ax3.imshow(kompozit['Tano'][:,central,:],vmin = -1.2,vmax=1.2,**imargs)
ax3.set_title('T anomaly (K)',fontsize=18)
ax3.set_xlabel('longitude (100 km)',fontsize=16)
#ax3.set_ylabel('altitude (km)')
plt.colorbar(im3,ax=ax3)

if figsav:
    plt.savefig(join('figs','kompo_vsect_VO3T_14Jan-22Feb.png'),**figargs)
plt.show()
plt.show()

#%% Calculation and plot of the composites for the increments
#%%
with gzip.open('OPZ-increment.pkl','rb') as f:
    datI = pickle.load(f)

#%% Kompozit
# The composit is made over the last ten days
kompozit['VOD'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
kompozit['O3D'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
kompozit['TD'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
for i in range(44,len(datI)):
    print(i)
    ix = np.where(dats[i].attr['lons']>=trac['lons'][i])[0][0]
    jy = np.where(dats[i].attr['lats']>=trac['lats'][i])[0][0]
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    kompozit['VOD'] += datI[i].var['VOD'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1]
    kompozit['O3D'] += datI[i].var['O3D'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1]
    kompozit['TD'] += datI[i].var['TD'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1]
ns = len(dats)-44
kompozit['VOD'] /= ns
kompozit['O3D'] /= ns
kompozit['TD'] /= ns
#%%
plt.figure(figsize=(13,14))
imargs = dict(origin='lower',aspect=1,interpolation='nearest',cmap='jet',extent=(-8,8,-8,8))
i=1
for k in np.arange(12,3,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(1.e5*kompozit['VOD'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('VOD level '+str(8-k)+u' ($10^{-5}~s^{-1}$)')
    i += 1
    plt.xlabel('longitude (degree)')
    plt.ylabel('latitude (degree)')
for k in np.arange(12,3,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(1.e6*kompozit['O3D'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('O3D  level'+str(8-k)+u' ($10^{-6}~kg~kg^{-1}$)')
    i += 1
    plt.xlabel('longitude (degree)')
    plt.ylabel('latitude (degree)')
for k in np.arange(12,3,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(kompozit['TD'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('TD level '+str(8-k)+' (K)')
    i += 1
    plt.xlabel('longitude (degree)')
    plt.ylabel('latitude (dgree)')
if figsav:
    plt.savefig(join('figs','kompo_INC_8Feb.png'),**figargs)
plt.show()

#%% Plot of a selection of VO and ozone images at the level of the max vorticity
# make an image every day for the purpose of a movie

First_day_2plot = 146
Last_day_2plot = len(dats)
#First_day_2plot = 12
#Last_day_2plot = 17

# beginning of the exploration (from get_traject)

for i in np.arange(First_day_2plot,Last_day_2plot,2):
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    dats[i].show('VO',kz,clim=(-2.e-5,5.e-5),
        txt='vorticity lev {:d} alt {:2.1f} km  {}  (s**-1)'.format(kz,trac['alts'][i],trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        savfile=join('figs','VOmaps','VO_'+str(i)+'.png'))
    date += timedelta(days=1)
#%%  Same for ozone

for i in np.arange(First_day_2plot,Last_day_2plot,2):
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    dats[i].show('O3',kz,
        txt='ozone lev {:d} alt {:2.1f} km  {}  (kg/kg)'.format(kz,trac['alts'][i],trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        savfile=join('figs','O3maps','O3_'+str(i)+'.png'))
    date += timedelta(days=1)

#%% Plot of the longitude altitude sections at the latitude of maximum vorticity
# vorticity field
# beginning of the exploration (from get_traject)

for i in np.arange(First_day_2plot,Last_day_2plot,2):
    dats[i].var['Z'] = datz[i].var['Z']
    dats[i].chartlonz('VO',trac['lats'][i],clim=(-2.e-5,5.e-5),levs=(20,80),
        txt='vorticity  {}  (s**-1)'.format(trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        savfile=join('figs','VOsect','VOsect_'+str(i)+'.png'))

#%% ozone

for i in np.arange(First_day_2plot,Last_day_2plot,2):
    dats[i].chartlonz('O3',trac['lats'][i],levs=(20,80),clim=(0,4.e-5),
        txt='ozone  {}  (kg/kg)'.format(trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        savfile=join('figs','O3sect','O3sect_'+str(i)+'.png'))

#%% temperature

for i in np.arange(First_day_2plot,Last_day_2plot,2):
    dats[i].chartlonz('T',trac['lats'][i],levs=(20,80),clim=(210,250),
        txt='Temperature {}  (K)'.format(trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
        savfile=join('figs','Tsect','Tsect_'+str(i)+'.png'))

#%% Tiler

# position is to be shifted by 6, done in the loop below
ll = {16:[0,4,8,12,18,24,28,32,36,40,44,48,52,56,60,64],
      15:[0,4,8,12,18,24,28,32,36,40,48,52,56,60,64]}

tasks = {1:['VOmaps','VO_',16,3693,955],2:['O3maps','O3_',16,3693,955],
         3:['VOsect','VOsect_',15,2834,1214],4:['O3sect','O3sect_',15,2820,1214],5:['Tsect','Tsect_',15,2630,1214]}

for r in [1,2]:
    print('task',r)
    chdir(join(rootdir,'figs',tasks[r][0]))
    [w0,h0] = [tasks[r][3],tasks[r][4]]
    bb = Image.new('RGB',(2*w0,8*h0))
    [w,h] = [0,0]
    im = []
    for n in range(16):
        im.append(Image.open(tasks[r][1]+str(ll[16][n])+'.png').crop([0,0,w0,h0]))
        print(ll[16][n]+6,im[n].size,w,h)
        bb.paste(im[n],(w,h,w+w0,h+h0))
        h = (h+h0) % (8*h0)
        if n==7: w = w0
    bb.save(tasks[r][1]+'All.png')
    bb.show()
    del im
#%%
for r in [3,4,5]:
    print('task',r)
    chdir(join(rootdir,'figs',tasks[r][0]))
    [w0,h0] = [tasks[r][3],tasks[r][4]]
    bb = Image.new('RGB',(3*w0,5*h0))
    [w,h] = [0,0]
    im = []
    for n in range(15):
        im.append(Image.open(tasks[r][1]+str(ll[15][n])+'.png').crop([0,0,w0,h0]))
        print(ll[15][n]+6,im[n].size,w,h)
        bb.paste(im[n],(w,h,w+w0,h+h0))
        h = (h+h0) % (5*h0)
        if n==4: w = w0
        if n==9: w = 2*w0
    bb.save(tasks[r][1]+'All.png')
    bb.show()
    del im

#%% Diag of the heating rate (rough version)
# Calculation of a mean theta and Pi profile
theta = np.zeros(137)
for i in range(len(dats)):
    theta += np.mean(dats[i].var['T'],axis=(1,2))
theta /= len(dats)
pif = (dats[0].attr['pscale']/100000)**cst.kappa
theta /= pif
zscale = dats[0].attr['zscale']
plt.plot(theta,zscale)
plt.ylim(17,35)
plt.xlim(400,1500)
plt.show()
plt.plot(pif,zscale)
plt.ylim(17,35)
plt.xlim(0.2,0.5)
plt.show()
#%%
# derivative of theta with z
dthetadz = (theta[1:]-theta[:-1])/(zscale[1:]-zscale[:-1])
zscale_i = 0.5*(zscale[1:]+zscale[:-1])
pif_i = 0.5*(pif[1:]+pif[:-1])
# Plot mean heating rate
# mean ascent rate (km/day)
ascent = 2*(trac['alts'][-1] - trac['alts'][0])/(len(dats)-1)
print('ascent rate (km/day)',ascent)
# heating rate (K/day)
plt.plot(ascent*dthetadz*pif_i,zscale_i)
plt.ylim(15,35)
plt.xlim(2,3.1)
plt.show()