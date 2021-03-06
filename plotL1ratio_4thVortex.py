#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a row of selected CALIOP images for the 4th vortex
Derived from plotL1ratio.py

@author: Bernard Legras
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle,gzip
#import argparse
from datetime import datetime, timedelta
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS, V
import matplotlib.colors as colors
import os
#from scipy.interpolate import RegularGridInterpolator
import socket
import scipy.signal as ss
#from scipy.ndimage.filters import gaussian_filter
#from Bresenham import line
from astropy.time import Time, TimeDelta
#import geosat
#from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.crs as ccrs

listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']
listcolors_sw=[listcolors[1],listcolors[0],listcolors[3],listcolors[2],\
               listcolors[5],listcolors[4],listcolors[7],listcolors[6],\
               listcolors[9],listcolors[8],listcolors[11],listcolors[10],\
               listcolors[13],listcolors[12],listcolors[15],listcolors[14],\
               listcolors[17],listcolors[16],listcolors[19],listcolors[18]]
mymap=colors.ListedColormap(listcolors)
mymap_sw=colors.ListedColormap(listcolors_sw)

trac = pickle.load(open('Vortex-track_4thVortex.pkl','rb'))

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

def get_vortex_pos(date):
    iv=0
    try:
        while trac['dates'][iv]<date:
            iv +=1
    except IndexError:
        print('no bracketting dates in vortex position')
        return None
    # interpolation coefficients
    print(iv)
    dt = (trac['dates'][iv]-trac['dates'][iv-1]).total_seconds()
    c1 = (date-trac['dates'][iv-1]).total_seconds()/dt
    c2 = (trac['dates'][iv]-date).total_seconds()/dt
    return [trac['lons'][iv]*c1+trac['lons'][iv-1]*c2,\
            trac['lats'][iv]*c1+trac['lats'][iv-1]*c2,\
            trac['z'][iv]*c1+trac['z'][iv-1]*c2,
            trac['vo'][iv]*c1+trac['vo'][iv-1]*c2,
            iv,c1,c2]

def get_vortex_section2(iv,c1,c2):
    # Simpler version that accounts the fact the longitude does not vary by more than 1°
    # see version with angle in plot1ratio.py
    ix = trac['ix'][iv]
    ix2 = trac['ix'][iv-1]
    if dats[iv].attr['lons'][0] == dats[iv-1].attr['lons'][0]:
        VOcut = dats[iv].var['VO'][...,ix]*c1 + dats[iv-1].var['VO'][...,ix2]*c2
        zcut = np.mean(dats[iv].var['Z'][...,ix]*c1 + dats[iv-1].var['Z'][...,ix2]*c2,axis=1)/1000
    else:
        print('dats version change detected')
        VOcut = dats[iv].var['VO'][...,ix]
        zcut = np.mean(dats[iv].var['Z'][...,ix],axis=1)/1000
    latcut = dats[iv].attr['lats']
    return (zcut,latcut,VOcut)

ysup=30
yinf=20

if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-Australia'
    dirL1_Std = '/dkol/data/CALIOP/CAL_LID_L1.v3.40'
    dirL1_Exp = '/dkol/data/CALIOP/CAL_LID_L1_Exp.v3.40'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-Australia'
    dirL1_Std = '/data/CALIOP/CAL_LID_L1.v3.40'
    dirL1_Exp = '/data/CALIOP/CAL_LID_L1_Exp.v3.40'
else:
    # ICARE setup
    rootdir = '/home/b.legras/STC/STC-Australia'
    dirL1_Std = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'
    dirL1_Exp = '/DATA/LIENS/CALIOP/CAL_LID_L1_Exp.v3.40'

#dirAProf = '/dkol/data/CALIOP/05kmAPro.v3.40'
#dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40'
Qs = 5.167*1.e-31
kbw = 1.0313

orbits = {1:[datetime(2020,2,7),'03-39-24',-68],
         2:[datetime(2020,2,8),'04-17-27',-62],
         3:[datetime(2020,2,15),'05-26-35',-60],
         4:[datetime(2020,2,17),'05-04-06',-61],
         5:[datetime(2020,2,23),'03-56-40',-50],
         6:[datetime(2020,2,24),'02-56-08',-53],
         7:[datetime(2020,3,2),'23-47-48',-50]
         }

il = [1,2,3,4,5,6,7]

#[53,-29.9,'15 Mar',False,15]

fig, ax = plt.subplots(1,7,sharey=True,figsize=(16,3))
ax = ax.flatten()

#il = {2:[95,-45.3,'11 Jan',False,12],9:[64,-32,'4 Mar',True,11]}
#il = {4:[195,-61.54,'23 Jan',False,12],9:[1,-32,'4 Mar',False,14]}
#il = {4:[195,-61.54,'23 Jan',False,12],8:[468,-40,'25 Feb',False,13]}
#il = {1:[61,-30,'16 Mar',False,15],6:[103,-29.1,'21 Mar',False,15]}
#fig, ax = plt.subplots(1,len(il),sharey=True,figsize=(10,10))

cmap = mymap_sw
cmap = 'gist_ncar'
#cmap = 'gist_rainbow_r'
#cmap = 'tab20b'
#
#fig, ax = plt.subplots(1,len(il),sharey=True,figsize=(10,10))
ifig = 0
for i in il:
    date = orbits[i][0]
    dirday = os.path.join(dirL1_Std,date.strftime('%Y/%Y_%m_%d'))
    file = os.path.join(dirday,date.strftime('CAL_LID_L1-ValStage1-V3-40.%Y-%m-%dT')+orbits[i][1]+'ZN.hdf')

    print('i',i)
    print(file)
    hdf = SD(file,SDC.READ)
    hh = HDF.HDF(file,HDF.HC.READ)
    lats = hdf.select('Latitude').get()[:].flatten()
    latmin = orbits[i][2]-10
    latmax = orbits[i][2]+10
    sel1L1 = (latmin <= lats) & (lats <= latmax)
    lons = hdf.select('Longitude').get()[sel1L1].flatten() % 360
    lats = lats[sel1L1]
    #utc = hdf.select('Profile_UTC_Time').get()[sel1L1].flatten()
    tai = hdf.select('Profile_Time').get()[sel1L1]
    tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
    utc = tt.utc.datetime.flatten()
    utcc = utc[0]+0.5*(utc[-1]-utc[0])
    mnd = hdf.select('Molecular_Number_Density').get()[sel1L1,:]
    lbeta532_met = np.log(1000 * mnd * Qs / (kbw*8*np.pi/3))
    t532 = np.ma.masked_less(hdf.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:],0)
    #t532 = t532fil[sel1L1,:]
    meta = hh.vstart().attach('metadata')
    alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
    meta = hh.vstart().attach('metadata')
    malts = np.array(meta.read()[0][meta.field('Met_Data_Altitudes')._idx])
    # calculation of the molecular backscatter
    #fint = RegularGridInterpolator((lats[::-1],malts[::-1]),lbeta532_33[::-1,::-1])
    #sr532 = t532/np.exp(fint(np.array(np.meshgrid(lats,alts,indexing='xy')).T))
    print('start interpolation')
    lbeta532_lid = np.empty(shape=t532.shape)
    for jy in range(len(lats)):
        lbeta532_lid[jy,:] = np.interp(alts,malts[::-1],lbeta532_met[jy,::-1])
    print('end interpolation')
    #fint = RegularGridInterpolator((malts[::-1]),lbeta532_33[::-1,::-1])
    #sr532 = t532/np.exp(fint(np.array(np.meshgrid(lats,alts,indexing='xy')).T))
    #sr532 = t532/np.exp(lbeta532_lid)
    # filter
    sr532raw = t532/np.exp(lbeta532_lid)
    print('starting filtering')
    sr532 = ss.medfilt(sr532raw,kernel_size=(81,1))
    print('ending filtering')
    #sr532 = gaussian_filter(ss.medfilt(sr532raw,kernel_size=(81,1)),4)
    #sr532 = ss.medfilt(sr1,kernel_size=(1,11))
    #sr532 = ss.medfilt(ss1,kernel_size=(1,5))
    #aa = np.mean(np.reshape(sr532[:12450,:],(249,50,583)),axis=1)
    #ll = np.mean(np.reshape(lats[:12450],(249,50)),axis=1)
    #sr532 = aa
    #lats = ll
    latmin = max(latmin,np.min(lats))
    latmax = min(latmax,np.max(lats))
    if (latmax-latmin)<19.9:
        print('bad range of lat values')
    # find corresponding longitudes to calculate the angle (not needed)
    lonmin = lons[-1] % 360
    lonmax = lons[0] % 360
    lonmid = 0.5*(lonmin+lonmax)
    if lonmax<lonmin: lonmid = (lonmid+180) % 360
    im=ax[ifig].pcolormesh(lats,alts,sr532.T,cmap=cmap,vmin=0,vmax=6)
    # Find the position of the vortex and generate the section
    try:
        [lonv,latv,zv,vomax,iv,c1,c2] = get_vortex_pos(utcc)
        [zcut, latcut, VOcut] = get_vortex_section2(iv,c1,c2)
        print('found corresponding vortex position',latv,zv)
        ax[ifig].plot(latv,zv,'firebrick',marker='+',ms=21,mew=3)
        ax[ifig].contour(latcut,zcut,VOcut,levels=[0.5*vomax,],colors='white',linewidths=2)
    except:
        print('no vertex position found')
    ax[ifig].set_xlim(latmin,latmax)
    ax[ifig].set_ylim(yinf,ysup)
    if lonmid>180:
        ax[ifig].set_title(date.strftime('%d %b')+'   '+str(int(360-lonmid))+'W')
    else:
        ax[ifig].set_title(date.strftime('%d %b')+'   '+str(int(lonmid))+'E')
    ax[ifig].set_xlabel('latitude')
    if ifig==0: ax[ifig].set_ylabel('altitude (km)')
    #plt.colorbar(im)
    #cid1 = fig.canvas.mpl_connect('button_press_event', on_click)
    #cid2 = fig.canvas.mpl_connect('key_press_event', on_key)
    ifig += 1
#cax = fig.add_axes([0.17,-0.04,0.67,0.05])
cax = fig.add_axes([0.92,0.12,0.015,0.76])
cbar = fig.colorbar(im,cax,orientation='vertical')
plt.savefig('figs/vortex4_kompo_7_filtered.png',dpi=300,bbox_inches='tight')
#plt.savefig('figs/vortex2_kompo_9_filtered.pdf',dpi=300,bbox_inches='tight')
plt.show()

