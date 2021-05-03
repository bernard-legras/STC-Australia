#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate history of the trajectory of Australian vortex 2
Localize the vortex
Print the positions
Calculate the total path.
Plot the evolution (latitude, longitude, altitude, vorticity)

Created on Monday 21 September 2020 from STC-CAN/tracker and STC-Australia/showDiags_2ndVortex

@author: Bernard Legras
"""
from datetime import datetime, timedelta
#from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import constants as cst
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D # yes it is used
#from matplotlib import cm
from matplotlib.text import TextPath
import gzip,pickle
import socket
import cartopy.crs as ccrs
#import deepdish as dd
from os.path import join
#from PIL import Image, ImageDraw, ImageFont
#from os import chdir
#from scipy.optimize import curve_fit

#%%

def tracker(dats,var,lon,lat,lower=30,upper=0,idx=6,jdy=3):
    # extract volume surrounding the target point
    jy = np.where(dats.attr['lats']>=lat)[0][0]
    ix = np.where(dats.attr['lons']>=lon)[0][0]
    #print('jy, ix',jy,ix)
    jmin = max(jy-jdy,0)
    imin = max(ix-idx,0)
    jmax = min(jy+jdy+1,dats.nlat)
    imax = min(ix+idx+1,dats.nlon)
    sample = dats.var[var][upper:lower,jmin:jmax,imin:imax]
    # find the 3d index of the max vorticity in the sample cube
    # we are in the southern hemisphere
    if 'O3' in var:
        aa = np.unravel_index(np.argmin(sample, axis=None), sample.shape)
    elif var in ['LPV','PV','VO']:
        aa = np.unravel_index(np.argmax(sample, axis=None), sample.shape)
    # return [lat,lon,theta,p,z,PV,vo,T,O3,kz,jz,iz]
    #print(jy,ix,jmin,jmax,imin,imax,aa)
    return([dats.attr['lats'][jmin+aa[1]],
            dats.attr['lons'][imin+aa[2]],
            dats.var['PT'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            cst.p0 * (dats.var['T'][upper+aa[0],jmin+aa[1],imin+aa[2]]/dats.var['PT'][upper+aa[0],jmin+aa[1],imin+aa[2]])**(1/cst.kappa),
            dats.var['Z'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['PV'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['VO'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['T'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            dats.var['O3'][upper+aa[0],jmin+aa[1],imin+aa[2]],
            upper+aa[0],jmin+aa[1],imin+aa[2]])

if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-CAN'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-CAN'

figsav = False
figargs = dict(bbox_inches='tight',dpi=300)
#%%
with gzip.open('OPZ-extract-2ndVortex_OPZPV.pkl','rb') as f:
    dats = pickle.load(f)
print(len(dats))

#%%
figsav = False
figargs = dict(bbox_inches='tight',dpi=300)
# tracking of the 3D position of the vortex
trac={'dates':[],'lons':[],'lats':[],'vo':[],'z':[],'T':[],'p':[],
      'pv':[],'pt':[],'o3':[],'ix':[],'jy':[],'kz':[]}
# initial position
lon = 180
lat = -30
date = datetime(2020,1,5,6)
# levs are to be shifted by 40 according to the 137 model levels due
# to truncation in the extraction
#for i in range(len(dats)):
for i in range(60):
    print(i)
    idx = 6
    jdy = 4
    # old vertical bounds
    #lower = 66
    #upper= 35
    upper = 0
    # this bound seems quite important at least at the beginning
    lower = 26
    var = 'O3prime'
    #var = 'LPV'
    #var = 'VO'
    # track with one iteration
    try:
        lon = 2*lon - trac['lons'][-2]
        lat = 2*lat - trac['lats'][-2]
        upper = trac['kz'][-1] - 2
        lower = min(26,trac['kz'][-1] + 3)
    except: pass
    #if i ==7: #old 3
    if i == 3:
        lon=196
    #if i ==13: # old 6
    if i == 6:
        lat = -32; lon = 224
    #if i ==25:  # old 12
    if i == 12:
        lat = -36; lon = 224
    #if i in [41,42,43]: # old 20,21
    if i in [20,21]:
        #upper = 55
        upper = 15
    #if i == 53: # old 26
    if i == 26:
        lon -= 360
    #if i == 54: # old 27
    if i == 27:
        lon += 360
    #if i in [97,102,112,113,114,115,116,117,118,119]: # old 48,51,56,57,58,59
    if i in [ 48,51,56,57,58,59]:
        #lower = 55
        #upper = 46
        lower = 15
        upper = 6

    i1 = i
    dats[i1].var['LPV'] = dats[i1].var['PV']*(dats[i1].var['PT']/500)**(-4.5)
    meanO3 = np.mean(dats[i].var['O3'],axis=(1,2))
    #dats[i].var['O3prime'] = dats[i].var['O3']/meanO3[:,None,None] -1
    dats[i1].var['O3prime'] = dats[i1].var['O3']- meanO3[:,None,None]

    try:
        [lat,lon,pt,p,z,pv,vo,T,o3,kz,jy,ix] = tracker(dats[i1],var,lon,lat,upper=upper,lower=lower,idx=idx,jdy=jdy)
        #[lat,lon,pt,p,z,pv,vo,T,o3,kz,jy,ix] = tracker(dats[i],var,lon,lat,upper=upper,lower=lower,idx=idx,jdy=jdy)
    except:
        print('tracking error at step ',i)
        print('tracking terminated')
        break
    trac['dates'].append(date)
    trac['lons'].append(lon)
    trac['lats'].append(lat)
    trac['vo'].append(vo)
    trac['pv'].append(pv)
    trac['z'].append(z/1000)
    trac['T'].append(T)
    trac['o3'].append(o3)
    trac['p'].append(p)
    trac['pt'].append(T*(cst.p0/p)**cst.kappa)
    trac['kz'].append(kz)
    trac['jy'].append(jy)
    trac['ix'].append(ix)
    date += timedelta(hours=12)

if len(trac['dates']) == len(dats):
    print('CONGRATULATIONS, YOU TRACKED THE VORTEX UNTIL THE END')
else:
    print('Tracking performed until ',trac['dates'][-1])
if var=='LPV':
    trac_PV = trac.copy()
    trac1 = trac_PV
elif var == 'O3prime':
    trac_O3 = trac.copy()
    trac2 = trac_O3
elif var == 'VO':
    trac_VO = trac.copy()
    trac3 = trac_VO

#%% print the positions as a function of time
for i in range(len(trac1['dates'])):
    # kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    print(i,trac1['dates'][i],trac1['lons'][i],trac1['lats'][i],'{:2.1f} {:2.1f}'.format(trac1['z'][i],trac3['vo'][i]*1.e5),trac1['kz'][i])
    #print(i,trac2['dates'][i],trac2['lons'][i],trac2['lats'][i],'{:2.1f}'.format(trac2['z'][i]),trac2['kz'][i])
    print(i,trac3['dates'][i],trac3['lons'][i],trac3['lats'][i],'{:2.1f} {:2.1f}'.format(trac3['z'][i],trac3['vo'][i]*1.e5),trac3['kz'][i])
# beware that saving here will loose the wind tracking made in wind-census
#pickle.dump(trac,open('Vortex-track.pkl','wb'))
#%% Loading the version of the track that has alos the maxvind
#trac = pickle.load(open('Vortex-track-withwind.pkl','rb'))
#mean_track = pickle.load(open('Koobor-Ntrack-L1-mean.pkl','rb'))
#fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))

#%% Checker

#for i in range(50,len(trac1['dates'])):
for i in range(0,20):
    i1 = i
    kz1 = trac1['kz'][i]
    #kz2 = trac2['kz'][i]
    kz3 = trac3['kz'][i]
    pt1 = trac1['pt'][i]
    #pt2 = trac2['pt'][i]
    pt3 = trac3['pt'][i]
    kz2 = kz3
    pt2 = pt3
    trac2 = trac3
    fig = plt.figure(figsize=(11,12))
    if dats[i1].attr['lons'][-1] > 180: cm_lon=180
    else: cm_lon=0
    projplate = ccrs.PlateCarree(central_longitude=cm_lon)
    gs = fig.add_gridspec(3,1)
    ax1 = fig.add_subplot(gs[0,0],projection=projplate)
    datr = dats[i1].interpolPT([pt1,pt2,pt3],varList=['PV','O3','VO'])
    ax1 = datr.show('PV',0,figsize=None,axf=ax1,show=False,scale=1.e6,xylim=True,
                 txt='PV '+str(i)+'  '+str(kz1)+trac1['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pt1)+'K'+
                 '   lon'+str(trac1['lons'][i])+'   lat'+str(trac1['lats'][i]))
    ax1.plot(trac1['lons'][i]-cm_lon,trac1['lats'][i],'black',marker='x',ms=21,mew=2)
    ax2 = fig.add_subplot(gs[1,0],projection=projplate)
    ax2 = datr.show('O3',1,figsize=None,axf=ax2,show=False,scale=1.e6,xylim=True,
                  txt='O3 '+str(i)+'  '+str(kz2)+trac2['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pt2)+'K'+
                  '   lon'+str(trac2['lons'][i])+'   lat'+str(trac2['lats'][i]))
    ax2.plot(trac2['lons'][i]-cm_lon,trac2['lats'][i],'firebrick',marker='x',ms=21,mew=1)
    ax3 = fig.add_subplot(gs[2,0],projection=projplate)
    ax3 = datr.show('VO',2,figsize=None,axf=ax3,show=False,scale=1.e5,xylim=True,
                 txt='VO '+str(i)+'  '+str(kz3)+trac3['dates'][i].strftime('  %Y %b %d %HUTC  ')+str(pt3)+'K'+
                 '   lon'+str(trac3['lons'][i])+'   lat'+str(trac3['lats'][i]))
    ax3.plot(trac3['lons'][i]-cm_lon,trac3['lats'][i],'black',marker='x',ms=21,mew=2)
    plt.show()

#%% Total displacement
trac['lats'] = np.array(trac['lats'])
trac['lons'] = np.array(trac['lons'])
dy = trac['lats'][1:]-trac['lats'][:-1]
dx = (trac['lons'][1:]-trac['lons'][:-1])*np.cos(np.deg2rad(0.5*(trac['lats'][1:]+trac['lats'][:-1])))
# Correction for crossing Greenwich
dx[149] = ((trac['lons'][150]-trac['lons'][149]%360))*np.cos(np.deg2rad(0.5*(trac['lats'][149]+trac['lats'][150])))
ds = np.sqrt(dx**2+dy**2)
print('total path ',(2*np.pi*6371/360)*np.sum(ds))

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