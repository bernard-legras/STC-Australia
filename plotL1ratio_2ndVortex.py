#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a rows of CALIOP images for the 2nd vortex
Derived from plotL1ratio

Created on Wed Feb 12 01:09:16 2020

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

trac = pickle.load(open('Vortex-track_2ndVortex.pkl','rb'))
with gzip.open('OPZ-extract-2ndVortex.pkl','rb') as f:
    dats = pickle.load(f)

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

ysup=24
yinf=11

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

orbits = {1:[datetime(2020,1,7),'10-18-01',-35.5],
         2:[datetime(2020,1,8),'10-55-59',-30],
         3:[datetime(2020,1,9),'09-55-27',-37],
         4:[datetime(2020,1,10),'10-33-26',-33],
         5:[datetime(2020,1,12),'10-10-52',-32],
         6:[datetime(2020,1,13),'09-10-15',-34],
         7:[datetime(2020,1,16),'07-47-09',-45],
         8:[datetime(2020,1,17),'05-08-07',-52],
         9:[datetime(2020,1,18),'04-07-35',-46],
         10:[datetime(2020,1,19),'03-06-58',-40],
         11:[datetime(2020,1,20),'02-06-26',-43],
         12:[datetime(2020,1,22),'23-04-51',-43],
         13:[datetime(2020,1,23),'22-04-19',-40],
         14:[datetime(2020,1,24),'21-03-47',-39],
         15:[datetime(2020,1,25),'20-03-15',-40],
         16:[datetime(2020,1,28),'15-23-09',-45],
         17:[datetime(2020,1,29),'14-27-39',-47],
         18:[datetime(2020,1,30),'13-22-05',-47],
         19:[datetime(2020,1,31),'12-21-33',-44],
         20:[datetime(2020,2,4),'13-15-06',-42]}

il = [1,4,6,7,10,14,16,19]

#[53,-29.9,'15 Mar',False,15]

fig, ax = plt.subplots(1,8,sharey=True,figsize=(18,3))
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
    lbeta512_met = np.log(1000 * mnd * Qs / (kbw*8*np.pi/3))
    t512 = np.ma.masked_less(hdf.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:],0)
    #t512 = t512fil[sel1L1,:]
    meta = hh.vstart().attach('metadata')
    alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
    meta = hh.vstart().attach('metadata')
    malts = np.array(meta.read()[0][meta.field('Met_Data_Altitudes')._idx])
    # calculation of the molecular backscatter
    #fint = RegularGridInterpolator((lats[::-1],malts[::-1]),lbeta512_33[::-1,::-1])
    #sr512 = t512/np.exp(fint(np.array(np.meshgrid(lats,alts,indexing='xy')).T))
    print('start interpolation')
    lbeta512_lid = np.empty(shape=t512.shape)
    for jy in range(len(lats)):
        lbeta512_lid[jy,:] = np.interp(alts,malts[::-1],lbeta512_met[jy,::-1])
    print('end interpolation')
    #fint = RegularGridInterpolator((malts[::-1]),lbeta512_33[::-1,::-1])
    #sr512 = t512/np.exp(fint(np.array(np.meshgrid(lats,alts,indexing='xy')).T))
    #sr512 = t512/np.exp(lbeta512_lid)
    # filter
    sr512raw = t512/np.exp(lbeta512_lid)
    print('starting filtering')
    sr512 = ss.medfilt(sr512raw,kernel_size=(81,1))
    print('ending filtering')
    #sr512 = gaussian_filter(ss.medfilt(sr512raw,kernel_size=(81,1)),4)
    #sr512 = ss.medfilt(sr1,kernel_size=(1,11))
    #sr512 = ss.medfilt(ss1,kernel_size=(1,5))
    #aa = np.mean(np.reshape(sr512[:12450,:],(249,50,583)),axis=1)
    #ll = np.mean(np.reshape(lats[:12450],(249,50)),axis=1)
    #sr512 = aa
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
    # Find the position of the vortex and generate the section
    [lonv,latv,zv,vomax,iv,c1,c2] = get_vortex_pos(utcc)
    #[zcut, latcut, VOcut] = get_vortex_section(iv,c1,c2,latmin,latmax,lonmin,lonmax,latv,lonv)
    [zcut, latcut, VOcut] = get_vortex_section2(iv,c1,c2)
    print('found corresponding vortex position',latv,zv)
    #im=ax.pcolormesh(lats,alts,t512.T,cmap=mymap_sw,norm=norm)
    im=ax[ifig].pcolormesh(lats,alts,sr512.T,cmap=cmap,vmin=0,vmax=12)
    ax[ifig].plot(latv,zv,'firebrick',marker='+',ms=21,mew=3)
    ax[ifig].contour(latcut,zcut,VOcut,levels=[0.5*vomax,],colors='white',linewidths=2)
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
plt.savefig('figs/vortex2_kompo_8_filtered.png',dpi=300,bbox_inches='tight')
#plt.savefig('figs/vortex2_kompo_9_filtered.pdf',dpi=300,bbox_inches='tight')
plt.show()

