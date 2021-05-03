#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates images for all the selected views of CALIOP to be assembled as a GIF

The corresponding location of the vortex is shown as a cross and a white line shows the contour
of the half max vorticity in the vertical meridian plane passing by the center of the
vortex projected onto that of the orbit.

The image is centered in lat on the cross.

Created on Sun 6 Sept 2020

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
from os.path import join
from scipy.signal import savgol_filter
#from scipy.interpolate import RegularGridInterpolator
import socket
import scipy.signal as ss
from ECMWF_N import ECMWF
#from scipy.ndimage.filters import gaussian_filter
#from Bresenham import line
#import geosat
#from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs

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

# Read the file containing the track
trac = pickle.load(open('Vortex-track.pkl','rb'))
trac['lats'] = np.array(trac['lats'])
trac['lons'] = np.array(trac['lons'])
# smooth the trajectory and reset lon between 0 and 360
tlons = savgol_filter(trac['lons'],17,3) % 360
tlats = savgol_filter(trac['lats'],17,3)
# Load the PZ extractions
with gzip.open('OPZ-extract-1.pkl','rb') as f:
    dats = pickle.load(f)
with gzip.open('OPZ-extract-2.pkl','rb') as f:
    dats.update(pickle.load(f))
with gzip.open('OPZ-Z-1.pkl','rb') as f:
    datz = pickle.load(f)
with gzip.open('OPZ-Z-2.pkl','rb') as f:
    datz.update(pickle.load(f))

def get_vortex_pos(date):
    ''' This function interpolates the position of the vortex at the time of the image
    from adjacent analysis '''
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

def get_vortex_section2(iv,c1,c2):
    ''' Generate a vorticity section in the meridian plane of the vortex core.
        The interpolation is made using the corefficients calculated in get_vortex_pos.py '''
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

def get_vortex_map(iv,lonv,latv,c1,c2):
    ''' Interpolate vorticity between the two dates using the raw date and ECMWF_N package '''
    dx = 30
    dy = 20
    dat1 = ECMWF('OPZ',trac['dates'][iv])
    dat2 = ECMWF('OPZ',trac['dates'][iv-1])
    dat1._get_var('VO')
    dat2._get_var('VO')
    dat1.close()
    dat2.close()
    lo1 = lonv - dx
    lo2 = lonv + dx
    la1 = latv - dy
    la2 = latv + dy
    # explicit size because the last level is often missing
    dat1.var['VO'][:136,...] = dat1.var['VO'][:136,...]*c1 + dat2.var['VO'][:136,...]*c2
    del dat2
    if ((lo1 < 0) | (lo2 >= 360)):
        dats = dat1.shift2west(lon0=-179)
        if lo2 >= 360:
            lo1 -= 360
            lo2 -= 360
    else:
        dats = dat1
    datr = dats.extract(latRange=(la1,la2),lonRange=(lo1,lo2),varss=['VO'])
    return datr

# Open the CALIOP selections
with gzip.open('selCaliop_12.pkl','rb') as f:
    sel,box,idate0,nbday,Cald12 = pickle.load(f)
with gzip.open('selCaliop_13.pkl','rb') as f:
    _,_,_,_,Cald13 = pickle.load(f)
with gzip.open('selCaliop_14.pkl','rb') as f:
    _,_,_,_,Cald14 = pickle.load(f)
with gzip.open('selCaliop_15.pkl','rb') as f:
    _,_,_,_,Cald15 = pickle.load(f)
with gzip.open('selCaliop_16.pkl','rb') as f:
    _,_,_,_,Cald16 = pickle.load(f)
with gzip.open('selCaliop_Exp_11.pkl','rb') as f:
    _,_,_,_,Cald_Exp = pickle.load(f)

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

sel = 14
vertical = False
horizontal = True

if sel == 12:
    late = False
    expedited = False
    listi = [16,24,32,40,55,64,72,81,95,104,112,121,138,154,162,180,195,211,228,237,255,264,273,316,326,
             335,344,353,362,372,399,408,417,426,435]
    lats_min = {16:-40,24:-43,32:-46}
    listi = [255,]
    Cald = Cald12
elif sel == 13:
    late = False
    expedited = False
    listi = [439,449,458,468,478,515]
    Cald = Cald13
elif sel == 14:
    late = False
    expedited = False
    listi = [1,11,20,30]
    Cald = Cald14
elif sel == 15:
    late = True
    expedited = False
    listi = [61,67,76,85,94,103,112,121,129]
    Cald = Cald15
elif sel == 16:
    late = True
    expedited = False
    listi = [15,22,36,51]
    Cald = Cald16
elif sel == 11:
    late = True
    expedited = True
    listi = [96,97]
    lats_min = {96:-23,97:-23}
    Cald = Cald_Exp

cmap = mymap_sw
cmap = 'gist_ncar'

ysup=33
yinf=12
caxframe = [0.92,0.12,0.04,0.76]
vmax = 20
tmp_dir = 'tmp_gif2'
if late:
    ysup = 36
    yinf = 20
    vmax = 6

ifig = 0
for i in listi:

    if expedited:
        date = datetime.strptime(Cald_Exp[i]['fname'][:10],'%Y-%m-%d')
        dirday = join(dirL1_Exp,date.strftime('%Y/%Y_%m_%d'))
        file = join(dirday,'CAL_LID_L1_Exp-Prov-V3-40.'+Cald_Exp[i]['fname']+'.hdf')
        sel1L1 = Cald_Exp[i]['sel1L1'][:,0]
        Cald = Cald_Exp

    else:
        date = datetime.strptime(Cald[i]['fname'][:10],'%Y-%m-%d')
        dirday = join(dirL1_Std,date.strftime('%Y/%Y_%m_%d'))
        file = join(dirday,'CAL_LID_L1-ValStage1-V3-40.'+Cald[i]['fname']+'.hdf')
        sel1L1 = Cald[i]['sel1L1'][:,0]

    print('i',i)
    print(file)
    hdf = SD(file,SDC.READ)
    hh = HDF.HDF(file,HDF.HC.READ)

    # Quick n,' dirty fix for the 2 March
    if i==515:
        sel1L1 = np.roll(sel1L1,-4800)

    lons = hdf.select('Longitude').get()[sel1L1].flatten() % 360
    lats = hdf.select('Latitude').get()[sel1L1].flatten()
    utc = hdf.select('Profile_UTC_Time').get()[sel1L1].flatten()
    # get the molecular density and calculate the molecular backscatter
    mnd = hdf.select('Molecular_Number_Density').get()[sel1L1,:]
    lbeta512_met = np.log(1000 * mnd * Qs / (kbw*8*np.pi/3))
    # get the attenuated backscatter
    t512 = np.ma.masked_less(hdf.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:],0)
    meta = hh.vstart().attach('metadata')
    alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
    meta = hh.vstart().attach('metadata')
    malts = np.array(meta.read()[0][meta.field('Met_Data_Altitudes')._idx])
    # interpolate the molecular backscatter to the lidar levels
    lbeta512_lid = np.empty(shape=t512.shape)
    for jy in range(len(lats)):
        lbeta512_lid[jy,:] = np.interp(alts,malts[::-1],  lbeta512_met[jy,::-1])
    # calculate the molecular backscatter ratio and filter it spatially to redue noise
    sr512raw = t512/np.exp(lbeta512_lid)
    sr512 = ss.medfilt(sr512raw,kernel_size=(81,1))
    try:
        # Find the position of the vortex and generate the section
        [lonv,latv,zv,vomax,iv,c1,c2] = get_vortex_pos(Cald[i]['utc'])
        #[zcut, latcut, VOcut] = get_vortex_section(iv,c1,c2,latmin,latmax,lonmin,lonmax,latv,lonv)
        [zcut, latcut, VOcut] = get_vortex_section2(iv,c1,c2)
        print('found corresponding vortex position',latv,zv)
        latmin = max(latv-15,np.min(lats))
        #datr = get_vortex_map(iv,lonv,latv,c1,c2)
        #print('got vorticity map')
        vortex = True
    except:
        print('no vortex position found')
        vortex = False
        latmin = max(lats_min[i]-15,np.min(lats))

    latmax = min(latmin+30,np.max(lats))
    # find corresponding longitudes in the CALIOP sections
    # and calulate mid longitude,
    # can be used to calculate angle
    ii = np.where(lats>=latmin)[0][0]
    lonmin = lons[ii]
    ii = np.where(lats>=latmax)[0][0]
    lonmax = lons[ii]
    lonmid = 0.5*(lonmin+lonmax) % 360
    # plots the backscatter ratio

#%%
    if horizontal:
        fig = plt.figure(figsize=(9,3))
        fig, (ax0, ax1) = plt.subplots(1, 2, sharey = True)
        #ax0 = fig.add_subplot(gs[0,0])
        im=ax0.pcolormesh(lats,alts,sr512.T,cmap=cmap,shading='auto',vmin=0,vmax=vmax)
        # add the vortex center and the contour of mid vorticity
        if vortex:
            ax0.plot(latv,zv,'firebrick',marker='+',ms=21,mew=3)
            ax0.contour(latcut,zcut,VOcut,levels=[0.5*vomax,],colors='white',linewidths=2)
        # set limits of the plot
        ax0.set_xlim(latmin,latmax)
        ax0.set_ylim(yinf,ysup)
        if lonmid>180:
            ax0.set_title(date.strftime('%d %b   ')+str(int(360-lonmid))+'W')
        else:
            ax0.set_title(date.strftime('%d %b   ')+str(int(lonmid))+'E')
        ax0.set_xlabel('latitude')
        ax0.set_ylabel('altitude (km)')
        plt.colorbar(im,orientation='vertical')
        #plt.tight_layout()

        if vortex:
            ax1.contourf(latcut,zcut,VOcut)
            #ax1.plot(lonv-cm_lon,latv,'firebrick',marker='+',markersize=100)
            #ax1.set_title(u'Vorticity (10$^{-5}$ s$^{-1}$) level '+str(int(trac['p'][iv]/100))+' hPa')

        #plt.savefig(join(tmp_dir,'VOH_CALIOP_BS_{}_{}.jpg'.format(sel,i)),dpi=300,bbox_inches='tight')
        #plt.savefig(join(tmp_dir,'VOH_CALIOP_BS_{}_{}.png'.format(sel,i)),dpi=300,bbox_inches='tight')

        plt.show()



