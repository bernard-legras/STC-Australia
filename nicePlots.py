#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory from 7 January 2020
This sccript makes the composit images of vorticity, ozone and temperature
for the main vortex Koobor.
Version Vc is based on the dates for which the CALIOP sections are available while
version Vs is based on the dates for which the best TROPOMI data are available.

The procedure is based on projecting boxes selected around the vortex at several dates
onto a background field selected for on date. The projection is vertical for the horizontal
chart and in latitude for the vertical lontigude x altitude sections. The latitude projection
is done bases on the bined altitude of the center of the vortex (not the model level).

Created on Sat Feb  8 12:27:14 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
#import constants as cst
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D # yes it is used
#from matplotlib import cm
from matplotlib.text import TextPath
import gzip,pickle
#import socket
#import deepdish as dd
from os.path import join
#from PIL import Image, ImageDraw, ImageFont
#from os import chdir

trac = pickle.load(open('Vortex-track.pkl','rb'))

with gzip.open('OPZ-extract-1.pkl','rb') as f:
    dats = pickle.load(f)
with gzip.open('OPZ-extract-2.pkl','rb') as f:
    dats.update(pickle.load(f))
with gzip.open('OPZ-Z-1.pkl','rb') as f:
    datz = pickle.load(f)
with gzip.open('OPZ-Z-2.pkl','rb') as f:
    datz.update(pickle.load(f))

# Read Silvia's pointing of the locations of aerosol blobs in Sentinel 5
trackSV = pickle.load(gzip.open(join('figs','_Silvia','AI_1v.pkl'),'rb'))
# Read consensus position of the vortex
mean_track = pickle.load(open('Koobor-Ntrack-L1-mean.pkl','rb'))
# Read MLS positions from Sergey's file
fidMLS = open(join('figs','_Sergey','MLS_blob_coordinates.txt'),'r')
#skip header
date0 = datetime(2020,1,1,0)
aa = fidMLS.readlines(1)
aa = np.array([np.genfromtxt(x.rstrip('\n').split('\t')) for x in fidMLS.readlines()])
tracMLS = {'lons':aa[:,2] % 360,'lats':aa[:,3]}
tracMLS['dates'] = [date0 + timedelta(days=x[0],seconds=x[1]) for x in aa]

#%% Total displacement
trac['lats'] = np.array(trac['lats'])
trac['lons'] = np.array(trac['lons'])
dy = trac['lats'][1:]-trac['lats'][:-1]
dx = (trac['lons'][1:]-trac['lons'][:-1])*np.cos(np.deg2rad(0.5*(trac['lats'][1:]+trac['lats'][:-1])))
# Correction for crossing Greenwich
dx[149] = ((trac['lons'][150]-trac['lons'][149]%360))*np.cos(np.deg2rad(0.5*(trac['lats'][149]+trac['lats'][150])))
ds = np.sqrt(dx**2+dy**2)
print('total path ',(2*np.pi*6371/360)*np.sum(ds))

#%% Plot of a selection of VO and ozone images at the level of the max vorticity
# for a few selected dates in a wide domain extending 70S-20S, 130E, 300E
# The dates are given from their index in the trac file
date0 = datetime(2020,1,4,6)
for i in [40,]:
    kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    p = trac['p'][i]
    print('pressure ',p,kz)
    date = date0 + timedelta(hours=12*i)
    print(date)
    dat = ECMWF('OPZ',date)
    dat._get_var('T')
    dat._mkp()
    dat._get_var('VO')
    dat._get_var('O3')
    datp = dat.extract(varss=['VO','T','O3'],latRange=(-70,-20),lonRange=(130,300))
    ax=datp.show('VO',kz,show=False)
    ax.plot(np.array(trac['lons'])-180,trac['lats'],linewidth=4,color='yellow')
    ax.set_xlim(-50,120)
    plt.show()
#%%
datp = dat.extract(varss=['VO','T','O3'],latRange=(-66,-36),lonRange=(130,300))
ax=datp.show('VO',kz,show=False)
ax.plot(np.array(trac['lons'])-180,trac['lats'],linewidth=4,color='yellow')
ax.plot((mean_track['lons'] % 360) - 180,mean_track['lats'],'o',markersize=8,color='r')
#ax.plot(tracMLS['lons']-180,tracMLS['lats'],'d',markersize=8,color='m')
ax.plot((trackSV.T[2] % 360)-180,trackSV.T[3],'o',markersize=8,color='g')
ax.set_xlim(-50,120)
ax.set_aspect(2)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
#mm index in the trac list
mm = {6:['7 Jan',(0,-8)],14:['11 Jan',(-20,0)],22:['15 Jan',(8,0)],32:['20 Jan',(-40,0)],
       42:['25 Jan',(5,0)],52:['30 Jan',(5,0)],62:['4 Feb',(5,0)],72:['9 Feb',(-8,-9)],
       82:['14 Feb',(-8,-8)],92:['19 Feb',(8,8)]}
for d in mm:
    #path = TextPath(mm[d][1],mm[d][0])
    path = TextPath((1,-2),mm[d][0],size=2)
    ax.plot(trac['lons'][d],trac['lats'][d],marker='D',markersize=6,color='k')
    ax.plot(trac['lons'][d],trac['lats'][d],'',marker=path,markersize=100,color='red')

plt.show()
#%%
#plt.title('Horizontal displacement between 07-01-2020 and 20-02-2020')
#plt.plot((mean_track['lons'] % 360) - 180,mean_track['lats'],'o',markersize=8,color='r')
#plt.plot((mean_track['lons'] % 360) - 180,mean_track['lats'],'-b')

#%% Superposition of maps for various positions of the vortex

"""
Version vc
Sequence of the curtain
7 Jan, 11 Jan, 16 Jan, 23 Jan, 31 Jan, 9 Feb, 18 Feb, 25 Feb, 4 Mar
Chosen corresponding events
14: 11 Jan, 38: 23 Jan (background), 72: 9 Feb, 82: 14 Feb, 90: 18 Feb,
96: 21 Feb, 104: 25 Feb, 112: 29 Feb, 120: 4 Mar

Version vs
Sequence of the TROPOMI AI footprints
7 Jan, 15 Jan, 24 Jan, 17 Feb, 21 Feb, 25 Feb
Chosen corresponding events
22: 15 Jan, 40: 24 Jan (background), 72: 9 Feb, 78: 12 Feb, 88, 17 Feb, 96: 21 Feb, 104: 25 Feb, 112: 29 Feb, 120: 4 Mar

"""
i = 38
events = [14,72,82,90,96,104,112,120]
version = 'vc'
i = 40
events = [22,72,78,88,96,104,112,120]
version = 'vs'
wx = 10
wy = 8
wz = 8
xwest = 0
xeast = 300
ysouth = -70
ynorth = -20
cm = 180
date0 = datetime(2020,1,4,6)
kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
jy = np.where(dats[i].attr['lats']>=trac['lats'][i])[0][0]
p = trac['p'][i]
print('pressure ',p,kz)
date = date0 + timedelta(hours=12*i)
print(date)
dat = ECMWF('OPZ',date)
dat._get_var('T')
dat._mkp()
dat._mkz()
dat._get_var('VO')
dat._get_var('O3')
# remove zonal mean
O3zon = np.mean(dat.var['O3'],axis=2)
dat.var['O3'] -= O3zon[:,:,np.newaxis]
Tzon = np.mean(dat.var['T'],axis=2)
dat.var['T'] -= Tzon[:,:,np.newaxis]
#%% Make multiplot
# This sections loads the data and built composite data by projection in the horizontal
# and vertical plane replacing the backgroung field
# This part is common to all the projections that follow
datp = dat.extract(varss=['VO','Z','O3','T'],latRange=(ysouth,ynorth),lonRange=(xwest,xeast))
jy = np.where(datp.attr['lats']>=trac['lats'][i])[0][0]
# Loop over the complementary maps
# 11 Jan: 14, 12 Feb: 78, 19 Feb: 92, 25 Feb:104, 4 Mar:120
boox = {}
bovx = {}
# duplicate the 'VO' field for the latitude projection
datp.var['VOv'] = datp.var['VO'].copy()
datp.var['O3v'] = datp.var['O3'].copy()
datp.var['Tv'] = datp.var['T'].copy()
#events = [96,]
for i1 in events:
    kz1 = np.where(dats[i1].attr['zscale']<=trac['alts'][i1])[0][0]
    #jy1 = np.where(dats[i].attr['lats']>=trac['lats'][i1])[0][0]
    date1 = date0 + timedelta(hours=12*i1)
    print(date1)
    dat1 = ECMWF('OPZ',date1)
    dat1._get_var('VO')
    dat1._get_var('O3')
    dat1._get_var('T')
    dat1._mkp()
    dat1._mkz()
    O3zon = np.mean(dat1.var['O3'],axis=2)
    dat1.var['O3'] -= O3zon[:,:,np.newaxis]
    Tzon = np.mean(dat1.var['T'],axis=2)
    dat1.var['T'] -= Tzon[:,:,np.newaxis]
    datp1 = dat1.extract(varss=['VO','Z','O3','T'],latRange=(ysouth,ynorth),lonRange=(xwest,xeast))
    lon1 = trac['lons'][i1]
    lat1 = trac['lats'][i1]
    ix1 = int(lon1-datp1.attr['lons'][0])
    jy1 = int(lat1-datp1.attr['lats'][0])
    # vertical projection of the submap around the vortex from kz1 to kz lavel
    datp.var['VO'][kz,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1] = \
        datp1.var['VO'][kz1,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1]
    datp.var['O3'][kz,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1] = \
        datp1.var['O3'][kz1,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1]
    datp.var['T'][kz,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1] = \
        datp1.var['T'][kz1,jy1-wy:jy1+wy+1,ix1-wx:ix1+wx+1]
    # latitude projection of the submap around the vortex
    datp.var['VOv'][kz1-wz:kz1+wz+1,jy,ix1-wx:ix1+wx+1] = \
        datp1.var['VO'][kz1-wz:kz1+wz+1,jy1,ix1-wx:ix1+wx+1]
    datp.var['O3v'][kz1-wz:kz1+wz+1,jy,ix1-wx:ix1+wx+1] = \
        datp1.var['O3'][kz1-wz:kz1+wz+1,jy1,ix1-wx:ix1+wx+1]
    datp.var['Tv'][kz1-wz:kz1+wz+1,jy,ix1-wx:ix1+wx+1] = \
        datp1.var['T'][kz1-wz:kz1+wz+1,jy1,ix1-wx:ix1+wx+1]
    datp.var['Z'][kz1-wz:kz1+wz+1,jy,ix1-wx:ix1+wx+1] = \
        datp1.var['Z'][kz1-wz:kz1+wz+1,jy1,ix1-wx:ix1+wx+1]
    boox[i1] = np.array([[lon1-wx,lon1-wx,lon1+wx,lon1+wx,lon1-wx],
                         [lat1-wy,lat1+wy,lat1+wy,lat1-wy,lat1-wy]])
    ztop = datp1.var['Z'][kz1-wz-1,jy1,ix1]/1000
    zbot = datp1.var['Z'][kz1+wz,jy1,ix1]/1000
    bovx[i1] = np.array([[lon1-wx,lon1-wx,lon1+wx,lon1+wx,lon1-wx],
                         [zbot,ztop,ztop,zbot,zbot]])
#%%
# horizontal multiplot
# Plot the vorticity
ax=datp.show('VO',kz,show=False,clim=(-2,10),scale=1.e5,aspect=1.414)
bbox = dict(boxstyle='round4',fc='w')
ax.text(trac['lons'][i]-cm,trac['lats'][i]+wy+6,date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(boox[i1][0]-cm,boox[i1][1],'w',linewidth=2)
    date1 = date0 + timedelta(hours=12*i1)
    ax.text(trac['lons'][i1]-cm,trac['lats'][i1]+wy+5,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.plot(np.array(trac['lons'])-cm,trac['lats'],linewidth=3,color='yellow')
ax.plot((mean_track['centroid_lon'] % 360) - cm,mean_track['centroid_lat'],'o',markersize=8,color='r')
ax.plot((trackSV.T[2] % 360)-cm,trackSV.T[3],'P',markersize=8,color='m',alpha=0.7)

ax.set_xlim(xwest-cm,xeast-cm)
ax.set_title(u'vorticity (10$^{-5}$ s$^{-1}$)',fontsize=18)
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
plt.savefig('figs/multiKoobor_horizontal_VO_'+version+'.png',dpi=300,bbox_inches='tight')
plt.savefig('figs/multiKoobor_horizontal_VO_'+version+'.pdf',dpi=300,bbox_inches='tight')
plt.show()
#%%
# horizontal multiplot
# Plot the ozone anomaly
ax=datp.show('O3',kz,show=False,scale=1.e6,aspect=1.414)
ax.plot(np.array(trac['lons'])-cm,trac['lats'],linewidth=3,color='yellow')
#ax.plot((mean_track['centroid_lon'] % 360) - cm,mean_track['centroid_lat'],'o',markersize=8,color='r')
#ax.plot((trackSV.T[2] % 360)-cm,trackSV.T[3],'d',markersize=8,color='g')
bbox = dict(boxstyle='round4',fc='w')
ax.text(trac['lons'][i]-cm,trac['lats'][i]+wy+6,date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(boox[i1][0]-cm,boox[i1][1],'w',linewidth=2)
    date1 = date0 + timedelta(hours=12*i1)
    ax.text(trac['lons'][i1]-cm,trac['lats'][i1]+wy+5,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.set_xlim(xwest-cm,xeast-cm)
ax.set_title(u'ozone anomaly (mg kg$^{-1}$)',fontsize=18)
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
plt.savefig('figs/multiKoobor_horizontal_O3_'+version+'.png',dpi=300,bbox_inches='tight')
plt.savefig('figs/multiKoobor_horizontal_O3_'+version+'.pdf',dpi=300,bbox_inches='tight')
plt.show()
#%%
# horizontal multiplot
# Plot the temperature anomaly
ax=datp.show('T',kz,show=False,aspect=1.414)
ax.plot(np.array(trac['lons'])-cm,trac['lats'],linewidth=3,color='yellow')
#ax.plot((mean_track['centroid_lon'] % 360) - cm,mean_track['centroid_lat'],'o',markersize=8,color='r')
#ax.plot((trackSV.T[2] % 360)-cm,trackSV.T[3],'d',markersize=8,color='g')
bbox = dict(boxstyle='round4',fc='w')
ax.text(trac['lons'][i]-cm,trac['lats'][i]+wy+6,date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(boox[i1][0]-cm,boox[i1][1],'w',linewidth=2)
    date1 = date0 + timedelta(hours=12*i1)
    ax.text(trac['lons'][i1]-cm,trac['lats'][i1]+wy+5,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.set_xlim(xwest-cm,xeast-cm)
ax.set_title(u'temperature anomaly (K)',fontsize=18)
ax.set_xlabel('longitude')
ax.set_ylabel('latitude')
plt.show()
#%%
# vertical multiplot vorticity
ax=datp.chartlonz('VOv',trac['lats'][i],clim=(-2,10),levs=(24,80),show=False,scale=1.e5)
ax.text(trac['lons'][i],trac['z'][i]-4,date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(bovx[i1][0],bovx[i1][1],'w',linewidth=2)
    date1 = date0 + timedelta(hours=12*i1)
    if i1 in [78,82,96,104,112,120]: up=-6
    elif i1 in [22,88]: up = 6
    elif i1 == 14: up = -4
    elif i1 == 90: up = 6
    else: up = 7
    ax.text(trac['lons'][i1],trac['z'][i1]+up,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.plot(np.array(trac['lons']),trac['z'],linewidth=4,color='yellow')
ax.plot((mean_track['centroid_lon'] % 360),mean_track['centroid_alt'],'o',markersize=5,color='r')
ax.plot((mean_track['centroid_lon'] % 360),mean_track['top_alt'],'o',markersize=5,color='k')
ax.plot((mean_track['centroid_lon'] % 360),mean_track['bot_alt'],'o',markersize=5,color='w')
ax.set_xlim(xwest,xeast)
ax.set_ylim(15,34)
ax.set_title(u'                vorticity (10$^{-5}$ s$^{-1}$)',fontsize=16)
ax.set_xlabel('longitude',fontsize=14)
ax.set_ylabel('altitude (km)',fontsize=14)
ax.set_xticks((30,60,90,120,150,180,210,240,270))
ax.set_xticklabels(['30°E','60°E','90°E','120°E','150°E','180°E','150°W','120°W','90°W'],fontsize=14)
plt.savefig('figs/multiKoobor_vertical_VO_'+version+'.png',dpi=300,bbox_inches='tight')
plt.savefig('figs/multiKoobor_vertical_VO_'+version+'.pdf',dpi=300,bbox_inches='tight')
plt.show()
#%%
# vertical multiplot ozone
ax=datp.chartlonz('O3v',trac['lats'][i],levs=(24,80),show=False,scale=1.e6)
ax.text(trac['lons'][i],trac['z'][i]-4,date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(bovx[i1][0],bovx[i1][1],'k',linewidth=2,alpha=0.8)
    date1 = date0 + timedelta(hours=12*i1)
    if i1 in [78,82,96,104,112,120]: up=-6
    elif i1 in [22,88]: up = 6
    elif i1 == 14: up = -4
    elif i1 == 90: up = 6
    else: up = 7
    ax.text(trac['lons'][i1],trac['z'][i1]+up,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.plot(np.array(trac['lons']),trac['z'],linewidth=4,color='yellow')
#ax.plot((mean_track['centroid_lon'] % 360),mean_track['centroid_alt'],'o',markersize=5,color='r')
#ax.plot((mean_track['centroid_lon'] % 360),mean_track['top_alt'],'o',markersize=5,color='k')
#ax.plot((mean_track['centroid_lon'] % 360),mean_track['bot_alt'],'o',markersize=5,color='w')
ax.set_xlim(xwest,xeast)
ax.set_ylim(15,34)
ax.set_title(u'                ozone anomaly (mg kg$^{-1}$)',fontsize=16)
ax.set_xlabel('longitude',fontsize=14)
ax.set_ylabel('altitude (km)',fontsize=14)
ax.set_xticks((30,60,90,120,150,180,210,240,270))
ax.set_xticklabels(['30°E','60°E','90°E','120°E','150°E','180°E','150°W','120°W','90°W'],fontsize=14)
plt.savefig('figs/multiKoobor_vertical_O3_'+version+'.png',dpi=300,bbox_inches='tight')
plt.savefig('figs/multiKoobor_vertical_O3_'+version+'.pdf',dpi=300,bbox_inches='tight')
plt.show()

#%%
# vertical multiplot temperature
ax=datp.chartlonz('Tv',trac['lats'][i],levs=(24,80),show=False,clim=(-10,10))
ax.text(trac['lons'][i],trac['z'][i]-4,date.strftime("%d %b"),
        ha="center",va="center",size="14",bbox=bbox)
for i1 in events:
    ax.plot(bovx[i1][0],bovx[i1][1],'k',linewidth=2,alpha=0.8)
    date1 = date0 + timedelta(hours=12*i1)
    if i1 in [78,82,96,104,112,120]: up=-6
    elif i1 in [22,88]: up = 6
    elif i1 == 14: up = -4
    elif i1 == 90: up = 6
    else: up = 7
    ax.text(trac['lons'][i1],trac['z'][i1]+up,date1.strftime("%d %b"),
            ha="center",va="center",size="14",bbox=bbox)
ax.plot(np.array(trac['lons']),trac['z'],linewidth=4,color='yellow')
#ax.plot((mean_track['centroid_lon'] % 360),mean_track['centroid_alt'],'o',markersize=5,color='r')
#ax.plot((mean_track['centroid_lon'] % 360),mean_track['top_alt'],'o',markersize=5,color='k')
#ax.plot((mean_track['centroid_lon'] % 360),mean_track['bot_alt'],'o',markersize=5,color='w')
ax.set_xlim(xwest,xeast)
ax.set_ylim(15,34)
ax.set_xlabel('longitude',fontsize=14)
ax.set_ylabel('altitude (km)',fontsize=14)
ax.set_xticks((30,60,90,120,150,180,210,240,270))
ax.set_xticklabels(['30°E','60°E','90°E','120°E','150°E','180°E','150°W','120°W','90°W'],fontsize=14)
ax.set_title(u'               temperature anomaly (K)',fontsize=16)

plt.savefig('figs/multiKoobor_vertical_T_'+version+'.png',dpi=300,bbox_inches='tight')
plt.savefig('figs/multiKoobor_vertical_T_'+version+'.pdf',dpi=300,bbox_inches='tight')
plt.show()