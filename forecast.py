#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the forecasts for a selected number of dates.
The vortex is tracked in the forecast in the same way as in showDiags

Created on Mon Feb 17 00:43:23 2020

@author: Bernard Legras
"""

from datetime import datetime, timedelta
import numpy as np
#from zISA import zISA
import pickle,gzip
from PIL import Image, ImageDraw, ImageFont
import os
from ECMWF_N import ECMWF
import constants as cst
import sys
if hasattr (sys, 'tracebacklimit'):
    del sys.tracebacklimit
# Raise this class for "soft halt" with minimum traceback.
class Stop (Exception):
    def __init__ (self):
        sys.tracebacklimit = 0


# copied from showDiags
def tracker(dats,datz,lon,lat,upper=35,lower=65,idx=5,jdy=3):
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

trac = pickle.load(open('Vortex-track.pkl','rb'))
# list of the relevant dates
#update = True
Alldates = [datetime(2020,1,7,12),datetime(2020,1,11,12),datetime(2020,1,15,12),
         datetime(2020,1,19,12),datetime(2020,1,23,12),datetime(2020,1,27,12),
         datetime(2020,1,31,12),datetime(2020,2,4,12),datetime(2020,2,8,12),
         datetime(2020,2,12,12),datetime(2020,2,16,12),datetime(2020,2,20,12),
         datetime(2020,2,24,12),datetime(2020,2,28,12),datetime(2020,3,3,12),
         datetime(2020,3,7,12),datetime(2020,3,11,12),datetime(2020,3,15,12),
         datetime(2020,3,19,12)]
update = True
dates = [Alldates[-1],]
dates = [datetime(2020,1,27,12)]
#dates = Alldates[-5:]
#%%
# Must be run AFTER showDiags.py
if update:
    fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))

    print(len(fctrac))
else:
    fctrac = {}
for date in dates:
    fc = ECMWF('OPZFCST',date,step=0)
    fctrac[date] = {'step':[],'dates':[],'lons':[],'lats':[],'alts':[],'vo':[],'z':[],'T':[],'p':[],'pt':[],'o3':[]}
    # initialisation of the position from the tracking made in showDiags
    date1 = date-timedelta(hours=6)
    i = np.where([(dd >= date1) for dd in trac['dates']])[0][0]
    lon = 0.5*(trac['lons'][i]+trac['lons'][i+1])
    lat = 0.5*(trac['lats'][i]+trac['lats'][i+1])
    # Account for changes in longitude origin
    if date >= datetime(2020,3,7,12): lon = lon %360
    # initial position
    print('initial ',lon,lat)
    for step in np.arange(0,244,24):
        print(date,step)
        datefc = date + timedelta(hours=int(step))
        fc._get_var('T',step=step)
        fc._get_var('LNSP',step=step)
        fc._get_var('VO',step=step)
        fc._get_var('O3',step=step)
        fc._mkp()
        fc._mkpscale()
        fc._mkzscale()
        fc._mkz()
        if (datefc >= datetime(2020,3,4)) & (date < datetime(2020,3,7)):
            fc2 = fc.shift2west(-179)
        else:
            fc2 = fc
        # set lower like in showDiags
        lower = 65
        upper = 24
        if date>datetime(2020,1,27): lower = 50
        if date>datetime(2020,2,10): lower = 42
        if date>datetime(2020,3,1): lower = 42
        if date>datetime(2020,3,14): lower = 33
        # try to predict next position
        try:
            lon = (2*lon - fctrac[date]['lons'][-2])
            lat = 2*lat - fctrac[date]['lats'][-2]
        except:
            pass
            #if date == datetime(2020,2,28,12):
            #    lon = 60
            #    lat -36
            #else: pass
        # special pre choice to avoid unwanted transitions
        if date == datetime(2020,1,11,12):
            if step in [96,120]:
                lower = 53
        if date == datetime(2020,1,19,12):
            if step == 168:
                lon = 300; lat = -60; lower = 52
            elif step == 192:
                lon = 305; lat = -55; lower = 51
            elif step == 216:
                lon = 308; lat = -52; lower = 51
            elif step == 240:
                lon = 315; lat = -50; lower = 51
        if date == datetime(2020,1,23,12):
            if step in [120,144,168]:
                lon = 295; lat = -52; lower =50; upper = 43
            elif step in [216,240]:
                lower = 50
        if date == datetime(2020,1,27,12):
            if step in [168,192]:
                upper = 42

        if date == datetime(2020,3,7,12):
            if step == 24:
                lower = 30; lon = 295; lat = -31
            elif step == 48:
                lower = 30; lon = 282; lat = -32
            elif step == 72:
                lower = 30; lon = 274; lat = -30
            elif step ==96:
                lower = 32; lon = 259; lat = -30
            elif step ==120:
                lower = 32; lon = 249; lat = -30
            elif step ==144:
                lower = 32; lon = 236; lat = -30
            elif step ==168:
                lower = 32; lon = 231; lat = -30
            elif step==192:
                lower = 32; lon = 220; lat = -30
            elif step==216:
                lower = 32; lon = 210; lat = -30
        elif date == datetime(2020,3,11,12):
            if step == 24:
                lower = 30; lon = 252; lat = -32
            elif step == 48:
                lower = 30; lon = 242; lat = -32
            elif step == 72:
                lower = 30; lon = 230; lat = -31
            elif step ==96:
                lower = 32; lon = 220; lat = -31
            elif step ==120:
                lower = 32; lon = 212; lat = -31
            elif step ==144:
                lower = 32; lon = 212; lat = -31

        [alt,lat,lon,vo,z,T,o3,p,ix,jy,kz] = tracker(fc2,fc2,lon,lat,lower=lower,upper=upper)
        #[alt,lat,lon,vo,z,T,o3,p,ix,jy,kz] = tracker(fc2,fc2,lon,lat,lower=lower,upper=upper)
        print(lon,lat,z,vo,kz)
        fctrac[date]['dates'].append(datefc)
        fctrac[date]['alts'].append(alt)
        fctrac[date]['lons'].append(lon)
        fctrac[date]['lats'].append(lat)
        fctrac[date]['vo'].append(vo)
        fctrac[date]['z'].append(z/1000)
        fctrac[date]['T'].append(T)
        fctrac[date]['p'].append(p)
        fctrac[date]['pt'].append(T*(cst.p0/p)**cst.kappa)
        fctrac[date]['o3'].append(o3)
    fc.close()

pickle.dump(fctrac,open('Vortex-fctrack.pkl','wb'))

#%% Tile images of the anticyclone evolution for each forecast

fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))

with gzip.open('OPZ-extract-1.pkl','rb') as f:
    dats = pickle.load(f)
with gzip.open('OPZ-extract-2.pkl','rb') as f:
    dats.update(pickle.load(f))
# make tmp if not done
try: os.mkdir('tmp')
except: pass
#%% generate images from the run analysis and from the forecast
w0 = 3693; h0 = 955
vmin = -2.e-5
vmax = 10.e-5

#for date in fctrac:
for date in dates:
    if date > datetime(2020,3,7): vmax = 5.e-5
    date1 = date-timedelta(hours=6)
    i = np.where([(dd >= date1) for dd in trac['dates']])[0][0]
    kz = trac['kz'][i]
    try: # to pass the case where the lwda analysis, which is produced with some delay,
         # is not available
        dats[i].show('VO',kz,clim=(vmin,vmax),
                     txt='ana vorticity lev {:d} alt {:2.1f} km  {}  (s**-1)'.format(kz,trac['z'][i],trac['dates'][i].strftime('%d-%m-%Y  %H UTC')),
                     savfile='tmp/VO_dats.png')
        print('analys lwda   ',trac['dates'][i])
    except: pass
    fc = ECMWF('OPZFCST',date,step=0)
    print('forecast from ',date)
    for step in np.arange(0,244,24):
        j = int(step/24)
        date2 = fctrac[date]['dates'][j]
        kz2 = np.where(dats[i].attr['zscale']<=fctrac[date]['alts'][j])[0][0]
        #kz2 = 35
        fc._get_var('VO',step=step)
        # Choice of the plotting box, to match dats choice
        if date < datetime(2020,1,27):
            latRange = (-65,-35)
            lonRange = (210,330)
        elif date2 < datetime(2020,2,13):
            latRange = (-65,-35)
            lonRange = (180,300)
        elif date2 < datetime(2020,2,17):
            latRange = (-58,-28)
            lonRange = (90,210)
        elif date == datetime(2020,2,20,12):
            latRange = (-58,-28)
            lonRange = (60,180)
        elif date == datetime(2020,2,24,12):
            latRange = (-58,-28)
            lonRange = (20,140)
        elif date >= datetime(2020,2,28,12):
            latRange = (-50,-20)
            lonRange = (0,120)
        if date == datetime(2020,2,4,12):
            latRange = (-70,-40)
            lonRange = (180,300)
        if date == datetime(2020,1,19,12):
            latRange = (-70,-40)
            lonRange = (240,359)
        if date == datetime(2020,1,23,12):
            latRange = (-70,-40)
            lonRange = (240,359)
        if (date == datetime(2020,2,24,12)) & (step > 168):
            lonRange = (0,120)
        if (date == datetime(2020,2,28,12)) & (step >= 168):
            lonRange = (240,359)
        if (date == datetime(2020,3,3,12)) & (step >= 48):
            lonRange = (240,359)
        if (date == datetime(2020,3,7,12)):
            if step <= 72:lonRange = (240,359)
            else: lonRange = (150,270)
        if (date == datetime(2020,3,11,12)):lonRange = (160,280)
        if (date == datetime(2020,3,15,12)):lonRange = (120,240)
        if (date == datetime(2020,3,19,12)):lonRange = (90,210)
        print(lonRange)
        fc2 = fc.extract(lonRange=lonRange,latRange=latRange,varss=['VO',])
        print('step ',step,'min max vo ',fc2.var['VO'].min(),fc2.var['VO'].max())
        print(fctrac[date]['lons'][j],fctrac[date]['lats'][j],kz2)
        fc2.show('VO',kz2,clim=(vmin,vmax),
        txt='fcst vorticity lev {:d} alt {:2.1f} km  {}  (s**-1)'.format(kz2,fctrac[date]['z'][j],fctrac[date]['dates'][j].strftime('%d-%m-%Y  %H UTC')),
        savfile='tmp/VO_fcst_'+str(step)+'.png')
        del fc2
    fc.close()
    del fc
#%%
    # define the frame
    bb = Image.new('RGB',(2*w0,6*h0),color=(255,255,255))
    [w,h] = [0,0] ; im = [] ; n = 0
    im.append(Image.open('tmp/VO_dats.png').crop([0,0,w0,h0]))
    bb.paste(im[n],(w,h,w+w0,h+h0))
    h = h0
    for step in np.arange(0,244,24):
        im.append(Image.open('tmp/VO_fcst_'+str(step)+'.png').crop([0,0,w0,h0]))
        n += 1
        bb.paste(im[n],(w,h,w+w0,h+h0))
        h = (h+h0) % (6*h0)
        if n==5: w = w0
    bb.save(date.strftime('figs/VO_fcst/VO_fcst_%Y-%m-%d.png'))
    #bb.show()
    del im

#%% Survival set after inspection of the images produced by previous cell
i = 0
# The two dates 11 Mar and 15 Mar are inverted in fctrack
# and so the ordering of their length: 9,7 intead of 7,9
survival = [6,9,4,10,10,10,10,10,10,10,10,10,10,10,8,5,9,7,6]
for date in fctrac:
    fctrac[date]['survival'] = survival[i]
    i += 1
pickle.dump(fctrac,open('Vortex-fctrack.pkl','wb'))

#%% Plots the altitude of the forecasted vortex compared to the analysed trajectory
# vertical axis in altitude
import matplotlib.pyplot as plt
figsav = True
figargs = dict(bbox_inches='tight',dpi=300)

fig = plt.figure()
plt.plot(trac['dates'],trac['z'],linewidth=4)
for date in fctrac:
    ns = fctrac[date]['survival']+1
    plt.plot(fctrac[date]['dates'][:ns],fctrac[date]['z'][:ns],'k',linewidth=2)
plt.ylabel('Altitude (km)')
plt.xlabel('Date')
fig.autofmt_xdate()
if figsav:
    plt.savefig('figs/ascent_fcst.png',**figargs)
    plt.savefig('figs/ascent_fcst.pdf',**figargs)
plt.show()
#%% Vertical axis in potential temperature
fig = plt.figure()
theta = np.array(trac['T'])*(np.array(trac['p'])/cst.p0)**(-cst.kappa)
plt.plot(trac['dates'],theta,linewidth=4)
for date in fctrac:
    ns = fctrac[date]['survival']+1
    theta = np.array(fctrac[date]['T'])*(np.array(fctrac[date]['p'])/cst.p0)**(-cst.kappa)
    plt.plot(fctrac[date]['dates'][:ns],theta[:ns],'k',linewidth=2)
plt.ylabel('Potential temperature (K)')
plt.xlabel('Date')
fig.autofmt_xdate()
if figsav:
    plt.savefig('figs/ascent_fcst_PT.png',**figargs)
    plt.savefig('figs/ascent_fcst_PT.pdf',**figargs)
plt.show()
#%% Vorticity
fig = plt.figure()
plt.plot(trac['dates'],1.e5*np.array(trac['vo']),linewidth=4)
for date in fctrac:
    ns = fctrac[date]['survival']+1
    plt.plot(fctrac[date]['dates'][:ns],1.e5*np.array(fctrac[date]['vo'][:ns]),'k',linewidth=2)
plt.ylabel(u'Vorticity ($10^{-5} s^{-1}$)')
plt.xlabel('Date')
fig.autofmt_xdate()
if figsav:
    plt.savefig('figs/vorticity_fcst.png',**figargs)
    plt.savefig('figs/vorticity_fcst.pdf',**figargs)
plt.show()
#%% Ozone
fig = plt.figure()
plt.plot(trac['dates'],1.e6*np.array(trac['o3']),linewidth=4)
for date in fctrac:
    ns = fctrac[date]['survival']+1
    plt.plot(fctrac[date]['dates'][:ns],1.e6*np.array(fctrac[date]['o3'][:ns]),'k',linewidth=2)
plt.ylabel('Ozone (mg/kg)')
plt.xlabel('Date')
fig.autofmt_xdate()
if figsav:
    plt.savefig('figs/ozone_fcst.png',**figargs)
    plt.savefig('figs/ozone_fcst.pdf',**figargs)
plt.show()
#%% print track
fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))
# fix step quantity that was not set
step = np.arange(0,240+2,24).astype(np.int)
for date in Alldates:
    print()
    print('forecast from ',date)
    print('step (hour), date, lon, lat, alt (km)')
    i=0
    while i<=fctrac[date]['survival']:
        print(step[i],fctrac[date]['lons'][i],fctrac[date]['lats'][i],'{:2.1f}'.format(fctrac[date]['z'][i]))
        i += 1
raise Stop ()

#%% WORK AREA DO NOT ENTER FROM ABOVE
date = datetime(2020,1,27,12)
fc1 = ECMWF('OPZFCST',date,step=0)
date1 = date-timedelta(hours=6)
i1 = np.where([(dd >= date1) for dd in trac['dates']])[0][0]
lon = 0.5*(trac['lons'][i1]+trac['lons'][i1+1])
lat = 0.5*(trac['lats'][i1]+trac['lats'][i1+1])
#fc._get_var('T',step=step)
#%%
step=168
j1 = int(step/24)
kz1 = np.where(dats[i1].attr['zscale']<=fctrac[date]['alts'][j1])[0][0]
print(fctrac[date]['lons'][j1],fctrac[date]['lats'][j1],fctrac[date]['z'][j1],kz1)
#fc._get_var('LNSP',step=step)
fc1._get_var('VO',step=step)
#fc._get_var('O3',step=step)
#fc._mkp()
#fc._mkpscale()
#fc._mkzscale()
#fc._mkz()
latRange = (-50,-20)
lonRange = (240,359)
latRange = (-65,-35)
lonRange = (240,359)
#if step>=48: lonRange=(180,300)
fc2 = fc1.extract(lonRange=lonRange,latRange=latRange,varss=['VO',])
fc2.show('VO',kz1)