#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a list of positions of the Spirit center from manual pointing of CALIOP
images

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

trac = pickle.load(open('Vortex-track.pkl','rb'))
with gzip.open('OPZ-extract-1.pkl','rb') as f:
    dats = pickle.load(f)
with gzip.open('OPZ-extract-2.pkl','rb') as f:
    dats.update(pickle.load(f))
with gzip.open('OPZ-Z-1.pkl','rb') as f:
    datz = pickle.load(f)
with gzip.open('OPZ-Z-2.pkl','rb') as f:
    datz.update(pickle.load(f))

def get_vortex_pos(date):
    iv=0
    try:
        while trac['dates'][iv]<date:
            iv +=1
    except IndexError:
        print('no bracketting dates in vortex position')
        return None
    # interpolation coefficients
    dt = (trac['dates'][iv]-trac['dates'][iv-1]).total_seconds()
    c1 = (date-trac['dates'][iv-1]).total_seconds()/dt
    c2 = (trac['dates'][iv]-date).total_seconds()/dt
    return [trac['lons'][iv]*c1+trac['lons'][iv-1]*c2,\
            trac['lats'][iv]*c1+trac['lats'][iv-1]*c2,\
            trac['z'][iv]*c1+trac['z'][iv-1]*c2,
            trac['vo'][iv]*c1+trac['vo'][iv-1]*c2,
            iv,c1,c2]

# =============================================================================
# def get_vortex_section(iv,c1,c2,latmin,latmax,lonmin,lonmax,latv,lonv):
#     # Calculate the crossing angle of the retained region
#     # It is here assumed that the longitudes and latitudes are on a one-degree grid
#     # it does not work otherwise
#     # Inconvenients: complicated and generate a vorticity cut which can be
#     # somewhat damped especially if the two surrounding positions are distant in longitude
#     # To be kept as model for more general orbits that have a strong angle with respect to the meridian
#     dx = lonmax-lonmin
#     dy = latmax - latmin
#     # test to avoid falling on a boundary between two versions of the dats box
#     if dats[iv].attr['lons'][0] == dats[iv-1].attr['lons'][0]:
#         VO = dats[iv].var['VO']*c1 + dats[iv-1].var['VO']*c2
#         Z = datz[iv].var['Z']*c1 + datz[iv-1].var['Z']*c2
#         latv1 = latv
#         lonv1 = lonv
#     else:
#         print('dats version change detected')
#         VO = dats[iv].var['VO']
#         Z = datz[iv].var['Z']
#         latv1 = trac['lats'][iv]
#         lonv1 = trac['lons'][iv]
#     ll = line(int(latv1-0.5*dy),int(lonv1-0.5*dx),int(latv1+0.5*dy),int(lonv1+0.5*dx))
#     # generate indexes of the oblique slice in the dats horizontal box
#     # trick to handle vortices close to the lox lat boundary of data
#     OK = False
#     while OK==False:
#         try:
#             jdy = list(np.where(ll[0,0] == dats[iv].attr['lats'])[0][0] + ll[:,0]-ll[0,0])
#             OK = True
#         except:
#             print('delete an element of ll',len(ll))
#             aa = np.delete(ll,0,0)
#             ll = aa
#     # The longitudes are all put in the 0-360 range before testing
#     idx = list(np.where((ll[0,1] % 360) == (dats[iv].attr['lons'] % 360))[0][0] + ll[:,1]-ll[0,1])
#     # clissping the extremities
#     while jdy[-1] >= len(dats[iv].attr['lats']):
#         print("clip segment end")
#         del jdy[-1]
#         del idx[-1]
#     latcut = dats[iv].attr['lats'][jdy]
#     VOcut = VO[:,jdy,idx]
#     zcut = np.mean(Z[:,jdy,idx],axis=1)/1000
#     return (zcut,latcut,VOcut)
# =============================================================================

def get_vortex_section2(iv,c1,c2):
    # Simpler version that accounts the fact the longitude does not vary
    ix = trac['ix'][iv]
    ix2 = trac['ix'][iv-1]
    if dats[iv].attr['lons'][0] == dats[iv-1].attr['lons'][0]:
        VOcut = dats[iv].var['VO'][...,ix]*c1 + dats[iv-1].var['VO'][...,ix2]*c2
        zcut = np.mean(datz[iv].var['Z'][...,ix]*c1 + datz[iv-1].var['Z'][...,ix2]*c2,axis=1)/1000
    else:
        print('dats version change detected')
        VOcut = dats[iv].var['VO'][...,ix]
        zcut = np.mean(datz[iv].var['Z'][...,ix],axis=1)/1000
    latcut = dats[iv].attr['lats']
    return (zcut,latcut,VOcut)

#update = True
with gzip.open('selCaliop_15.pkl','rb') as f:
    _,_,_,_,Cald15 = pickle.load(f)
with gzip.open('selCaliop_16.pkl','rb') as f:
    _,_,_,_,Cald16 = pickle.load(f)
with gzip.open('selCaliop_Exp_11.pkl','rb') as f:
    _,_,_,_,Cald_Exp = pickle.load(f)

row=2
ysup=36
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

dirL1 = dirL1_Std
#dirAProf = '/dkol/data/CALIOP/05kmAPro.v3.40'
#dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40'
Qs = 5.167*1.e-31
kbw = 1.0313

# elements: center lat, date, Exp, sel
il = {1:[64,-52.8,'7 Jan',False,12],2:[95,-45.3,'11 Jan',False,12],3:[138,-48.8,'16 Jan',False,12],
      4:[195,-61.54,'23 Jan',False,12],5:[264,-52.66,'31 Jan',False,12],6:[344,-50.66,'9 Feb',False,12],
      7:[426,-45.,'18 Feb',False,12],8:[468,-40,'25 Feb',False,13],9:[1,-32,'4 mar',False,14]}

if row == 1:
    il = {1:[61,-30,'16 Mar',False,15],2:[67,-28,'17 Mar',False,15],3:[76,-30.1,'18 Mar',False,15],
      4:[85,-29.7,'19 Mar',False,15],5:[94,-30.5,'20 Mar',False,15],6:[103,-29.1,'21 Mar',False,15],
      7:[112,-27.7,'22 Mar',False,15],8:[121,-26.2,'23 Mar',False,15],9:[129,-26.1,'24 Mar',False,15]}
elif row == 2:
    il = {1:[15,-30,'27 Mar',False,16],2:[22,-30,'28 Mar',False,16],3:[36,-28,'30 Mar',False,16],
      4:[44,-27,'31 Mar',False,16],5:[51,-25,'1 Apr',False,16],6:[96,-23,'3 Apr',True,11],
      7:[97,-23,'4 Apr',True,11]}
#[53,-29.9,'15 Mar',False,15]

fig, ax = plt.subplots(1,9,sharey=True,figsize=(18,2.5))

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
    if il[i][3]:
        date = datetime.strptime(Cald_Exp[il[i][0]]['fname'][:10],'%Y-%m-%d')
        dirday = os.path.join(dirL1_Exp,date.strftime('%Y/%Y_%m_%d'))
        file = os.path.join(dirday,'CAL_LID_L1_Exp-Prov-V3-40.'+Cald_Exp[il[i][0]]['fname']+'.hdf')
        sel1L1 = Cald_Exp[il[i][0]]['sel1L1'][:,0]
        Cald = Cald_Exp
    else:
        if il[i][4]==15: Cald = Cald15
        else: Cald = Cald16
        date = datetime.strptime(Cald[il[i][0]]['fname'][:10],'%Y-%m-%d')
        dirday = os.path.join(dirL1_Std,date.strftime('%Y/%Y_%m_%d'))
        file = os.path.join(dirday,'CAL_LID_L1-ValStage1-V3-40.'+Cald[il[i][0]]['fname']+'.hdf')
        sel1L1 = Cald[il[i][0]]['sel1L1'][:,0]
    print('i',i)
    print(file)
    hdf = SD(file,SDC.READ)
    hh = HDF.HDF(file,HDF.HC.READ)
    #sel1L1 = Cald[i]['sel1L1'][:,0]
    if il[i][2]:
        daynite1 = hdf.select('Day_Night_Flag').get()[:]
        sel1L1 = np.array([x and y for (x,y) in zip(daynite1==1,sel1L1)]).astype(np.bool)
    lons = hdf.select('Longitude').get()[sel1L1].flatten() % 360
    lats = hdf.select('Latitude').get()[sel1L1].flatten()
    utc = hdf.select('Profile_UTC_Time').get()[sel1L1].flatten()
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
    sr512 = ss.medfilt(sr512raw,kernel_size=(161,1))
    print('ending filtering')
    #sr512 = gaussian_filter(ss.medfilt(sr512raw,kernel_size=(81,1)),4)
    #sr512 = ss.medfilt(sr1,kernel_size=(1,11))
    #sr512 = ss.medfilt(ss1,kernel_size=(1,5))
    #aa = np.mean(np.reshape(sr512[:12450,:],(249,50,583)),axis=1)
    #ll = np.mean(np.reshape(lats[:12450],(249,50)),axis=1)
    #sr512 = aa
    #lats = ll
    latmin = max(il[i][1]-10,np.min(lats))
    latmax = min(latmin+20,np.max(lats))
    if (latmax-latmin)<20:
        print('bad range of lat values')
    # find corresponding longitudes to calculate the angle
    ii = np.where(lats>=latmin)[0][0]
    lonmin = lons[ii]
    ii = np.where(lats>=latmax)[0][0]
    lonmax = lons[ii]
    lonmid = 0.5*(lonmin+lonmax) % 360
    im=ax[ifig].pcolormesh(lats,alts,sr512.T,cmap=cmap,vmin=0,vmax=6)
    # Find the position of the vortex and generate the section
    try:
        [lonv,latv,zv,vomax,iv,c1,c2] = get_vortex_pos(Cald[il[i][0]]['utc'])
        [zcut, latcut, VOcut] = get_vortex_section2(iv,c1,c2)
        print('found corresponding vortex position',latv,zv)
        ax[ifig].plot(latv,zv,'firebrick',marker='+',ms=21,mew=3)
        ax[ifig].contour(latcut,zcut,VOcut,levels=[0.5*vomax,],colors='white',linewidths=2)
    except:
        print('no vortex found')
    ax[ifig].set_xlim(latmin,latmax)
    ax[ifig].set_ylim(yinf,ysup)
    if lonmid>180:
        ax[ifig].set_title(il[i][2]+'   '+str(int(360-lonmid))+'W')
    else:
        ax[ifig].set_title(il[i][2]+'   '+str(int(lonmid))+'E')
    ax[ifig].set_xlabel('latitude')
    if ifig==0: ax[ifig].set_ylabel('altitude (km)')
    #plt.colorbar(im)
    #cid1 = fig.canvas.mpl_connect('button_press_event', on_click)
    #cid2 = fig.canvas.mpl_connect('key_press_event', on_key)
    ifig += 1
#cax = fig.add_axes([0.17,-0.04,0.67,0.05])
cax = fig.add_axes([0.92,0.12,0.015,0.76])
cbar = fig.colorbar(im,cax,orientation='vertical')
if row == 1:
    plt.savefig('figs/ascent_March_kompo_9_filtered.png',dpi=300,bbox_inches='tight')
elif row == 2:
    plt.savefig('figs/ascent_April_kompo_9_filtered.png',dpi=300,bbox_inches='tight')
#plt.savefig('figs/vortex2_kompo_9_filtered.pdf',dpi=300,bbox_inches='tight')
plt.show()

