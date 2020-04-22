# -*- coding: utf-8 -*-
"""
Plots CALIOP sections for the selected standrad V 3.40 orbits
Meant to be run on ICARE
"""
import pickle,gzip
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS, V
import os
import numpy as np
import matplotlib.colors as colors
import geosat
#from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import argparse
import scipy.signal as ss

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

dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40'
dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'
dirCProf = '/DATA/LIENS/CALIOP/05kmCPro.v3.40'

parser = argparse.ArgumentParser()
parser.add_argument("-s","--sel",type=int,choices=[0,1,2,3,4,5,6,7,8,12,13,14,15,16,17],help="selection")
parser.add_argument("-f","--firsti",type=int,help="first index")
parser.add_argument("-l","--lasti",type=int,help="last index")
parser.add_argument("-list","--list",type=int,nargs='+')

sel = 16
firsti = 122
figsav = True
show = False
vmax = 6
medfilter = True
widthfilter = 161

#sel = 17
#firsti = 1
#figsav = True
#show = False
#vmax = 6
#medfilter = True
#widthfilter = 81

Qs = 5.167*1.e-31
kbw = 1.0313

args = parser.parse_args()
if args.sel is not None: sel = args.sel
if args.firsti is not None: firsti = args.firsti

trac_1st = pickle.load(open('Vortex-track.pkl','rb'))
trac_2nd = pickle.load(open('Vortex-track_2ndVortex.pkl','rb'))

def get_vortex_pos(date,trac):
    i=0
    try:
        while trac['dates'][i]<date:
            i +=1
    except IndexError:
        print('no bracketting dates in vortex position')
        return None
    # interpolation coefficients
    dt = (trac['dates'][i]-trac['dates'][i-1]).total_seconds()
    c1 = (date-trac['dates'][i-1]).total_seconds()/dt
    c2 = (trac['dates'][i]-date).total_seconds()/dt
    # crossing dateline
    if i == 150:
        trac['lons'][149] = trac['lons'][149] % 360
    return [trac['lons'][i]*c1+trac['lons'][i-1]*c2,\
            trac['lats'][i]*c1+trac['lats'][i-1]*c2,\
            trac['z'][i]*c1+trac['z'][i-1]*c2]

# load dictionary of traces
with gzip.open('selCaliop_'+str(sel)+'.pkl','rb') as f:
    sel,box,idate0,nbday,Cald = pickle.load(f)

lasti = len(Cald)
if args.lasti is not None: lasti = args.lasti
listi = range(firsti,lasti+1)
if args.list is not None: listi = args.list

if sel in [0,2,3,4]:
    GEO = True
    gg = geosat.GeoGrid('HimFull')
    subgg = gg.subgrid([130,220,-60,-5])
    cm = 180
    ext = [130-180,220-180,-60,-5]
    xlims = [-40,120]
    xlocs=(140,160,-180,-160,-140,-120,-100,-80,-60)
    satload = geosat.Himawari
    itm = 20
    ysup = 25
elif sel == 12:
    GEO = True
    gg = geosat.GeoGrid('HimFull')
    subgg = gg.subgrid([130,220,-60,-5])
    cm = 180
    ext = [130-180,220-180,-70,-30]
    xlims = [-50,150]
    xlocs=(140,160,-180,-160,-140,-120,-100,-80,-60,-40)
    satload = geosat.Himawari
    itm = 20
    yinf = 12
    ysup = 33
elif sel==1:
    # Central Pacific during August and early September
    GEO = True
    gg = geosat.GeoGrid('HimFull')
    subgg = gg.subgrid(box)
    cm = 0
    ext = box
    xlims = [box[0],box[1]]
    xlocs = None
    satload = geosat.Himawari
    itm = 20
    ysup = 22
elif sel in [5,6,7]:
    # Indian Ocean during early September
    GEO = True
    gg = geosat.GeoGrid('MSG1Full')
    subgg = gg.subgrid([box[0],110,box[2],box[3]])
    cm = 0
    ext = box
    xlims = [box[0],110]
    xlocs = None
    satload = geosat.MSG1
    itm = 15
    ysup = 25
elif sel==13:
    # Atlantic and Indian Ocean
    GEO = True
    gg = geosat.GeoGrid('MSG1Full')
    subgg = gg.subgrid([-30,110,box[2],box[3]])
    cm = 0
    ext = [-30,110,box[2],box[3]]
    xlims = [box[0],box[1]]
    xlocs = None
    satload = geosat.MSG1
    itm = 15
    yinf = 12
    ysup = 33
elif sel==14:
    # America to Indian Ocean
    GEO = False
    cm = 0
    ext = [-110,100,-55,-15]
    xlims = [-110,100]
    xlocs = None
    yinf = 12
    ysup = 33
elif sel==15:
    # America to Indian Ocean
    GEO = False
    cm = 180
    ext = [100-180,300-180,-55,-15]
    xlims = [100-180,300-180]
    xlocs = (100,120,140,160,-180,-160,-140,-120,-100,-80,-60)
    yinf = 12
    ysup = 35
elif sel==16:
    # Atlantic & Indian Ocean
    GEO = False
    cm = 0
    ext = [-50,160,-50,-10]
    xlims = [-50,160]
    xlocs = None
    yinf = 20
    ysup = 38
elif sel==17:
    # Antarctic
    GEO = False
    cm = 0
    ext = [-180,-50,-85,-55]
    xlims = [-180,-50]
    xlocs = None
    yinf = 12
    ysup = 33
elif sel>=8:
    GEO = False
    ext = [-120,0,-70,-30]
    xlims = [-120,0]
    ysup = 28
    cm = 0
    xlocs = None

for i in listi:
#for i in range(331,len(Cald)+1):
#for i in [1,3]:

    # generate date and daily directory name
    date = datetime.strptime(Cald[i]['fname'][:10],'%Y-%m-%d')
    dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
    dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    dirdayC = os.path.join(dirCProf,date.strftime('%Y/%Y_%m_%d'))

    try:
        [lonv1,latv1,zv1] = get_vortex_pos(Cald[i]['utc'],trac_1st)
        print('found corresponding 1st vortex position')
        vortex_1st = True
    except:
        print('no 1st vortex position')
        vortex_1st = False
    try:
        [lonv2,latv2,zv2] = get_vortex_pos(Cald[i]['utc'],trac_2nd)
        print('found corresponding 2nd vortex position')
        vortex_2nd = True
    except:
        print('no 2nd vortex position')
        vortex_2nd = False

    print('i',i)
    print(Cald[i]['fname'])
    # try:
    #     file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
    #     hdf = SD(file,SDC.READ)
    #     hh = HDF.HDF(file,HDF.HC.READ)
    #     lons = hdf.select('Longitude').get()[Cald[i]['sel1'],1] % 360
    #     lats = hdf.select('Latitude').get()[Cald[i]['sel1'],1]
    #     meta = hh.vstart().attach('metadata')
    #     alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
    #     #a512 = np.ma.masked_less(hdf.select('Aerosol_Multiple_Scattering_Profile_532').get()[Cald[i]['sel1'],:],0)
    #     #e512 = np.ma.masked_less(hdf.select('Extinction_Coefficient_532').get()[Cald[i]['sel1'],:],0)
    #     #f512 = np.ma.masked_less(hdf.select('Aerosol_Layer_Fraction').get()[Cald[i]['sel1'],:],0)
    #     t512 = np.ma.masked_less(hdf.select('Total_Backscatter_Coefficient_532').get()[Cald[i]['sel1'],:],0)
    #     AProOK = True
    # except:
    #     AProOK = False
    #     print('missing aerosol data')

    # try:
    #     fileC = os.path.join(dirdayC,'CAL_LID_L2_05kmCPro-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
    #     hdfC = SD(fileC,SDC.READ)
    #     hhC = HDF.HDF(fileC,HDF.HC.READ)
    #     t512C = np.ma.masked_less(hdfC.select('Total_Backscatter_Coefficient_532').get()[Cald[i]['sel1'],:],0)
    #     metaC = hhC.vstart().attach('metadata')
    #     altsC = np.array(metaC.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
    #     CProOK = True
    # except:
    #     CProOK = False
    #     print('missing cloud data')

    try:
        fileL1 = os.path.join(dirdayL1,'CAL_LID_L1-ValStage1-V3-40.'+Cald[i]['fname']+'.hdf')
        hdf1 = SD(fileL1,SDC.READ)
        hh1 = HDF.HDF(fileL1,HDF.HC.READ)
        sel1L1 = Cald[i]['sel1L1'][:,0]

        t512L1 = hdf1.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:]
        t512L1 = np.ma.masked_less(t512L1,0)
        if sel == 13:
            lons1 = hdf1.select('Longitude').get()[sel1L1].flatten()
        else:
            lons1 = hdf1.select('Longitude').get()[sel1L1].flatten() % 360
        lats1 = hdf1.select('Latitude').get()[sel1L1].flatten()
        mnd1 = hdf1.select('Molecular_Number_Density').get()[sel1L1,:]
        lbeta512_met = np.log(1000 * mnd1 * Qs / (kbw*8*np.pi/3))
        meta1 = hh1.vstart().attach('metadata')
        alts1 = np.array(meta1.read()[0][meta1.field('Lidar_Data_Altitudes')._idx])
        meta1 = hh1.vstart().attach('metadata')
        malts1 = np.array(meta1.read()[0][meta1.field('Met_Data_Altitudes')._idx])
        # calculation of the molecular backscatter
        lbeta512_lid = np.empty(shape=t512L1.shape)
        for jy in range(len(lats1)):
            lbeta512_lid[jy,:] = np.interp(alts1,malts1[::-1],lbeta512_met[jy,::-1])
        if medfilter:
            sr512raw = t512L1/np.exp(lbeta512_lid)
            sr512= ss.medfilt(sr512raw,kernel_size=(widthfilter,1))
        else:
            sr512 = t512L1/np.exp(lbeta512_lid)
        L1OK = True
    except:
        print ('no L1 data')
        print (fileL1)
        L1OK = False

#%% Here we get the himawari scene which is the closest for this time

    GEOOK = False
    if GEO:
        try:
            utc = Cald[i]['utc']
            utc -= timedelta(minutes=(utc.minute%itm))
            try:
                ah = satload(utc)
            except:
                try:
                    ah = satload(utc+timedelta(minutes=itm))
                except:
                    try:
                        ah = satload(utc+timedelta(minutes=2*itm))
                    except:
                        ah = satload(utc-timedelta(minutes=itm))
            ah._get_IR0()
            ph = geosat.SatGrid(ah,gg)
            ph._sat_togrid('IR0')
            GEOOK = True
        except:
            GEOOK = False
    #%%
    fig = plt.figure(figsize=(6,8))
    if L1OK:
        plt.subplot(2,1,1)
        norm1 = colors.LogNorm(vmin=0.00001,vmax=0.1)
        print('dims',lats1.shape,alts1.shape,sr512.shape)
        #plt.pcolormesh(lats1,alts1,b512.T,cmap=mymap_sw,norm=norm1)
        if sel == 17:
            plt.pcolormesh(lons1,alts1,sr512.T,cmap='gist_ncar',vmin=0,vmax=6)
        else:
            plt.pcolormesh(lats1,alts1,sr512.T,cmap='gist_ncar',vmin=0,vmax=6)
            if vortex_1st : plt.plot(latv1,zv1,'r',marker='+',markersize=11)
            if vortex_2nd : plt.plot(latv2,zv2,'r',marker='x',markersize=11)
        plt.title('{:d} {:16} {:.2f} {:.2f}'.format(i,Cald[i]['utc'].strftime('%Y-%m-%dT%H:%M'),lons1[0],lons1[-1]))
        plt.ylim(yinf,ysup)
        plt.title('Scattering ratio 512')
        plt.colorbar()
    #if AProOK:
    #    plt.subplot(2,2,2)
    #    norm = colors.LogNorm(vmin=0.00001,vmax=0.01)
    #    plt.pcolormesh(lats,alts,t512.T,cmap=mymap_sw,norm=norm)
    #    plt.ylim(yinf,ysup)
    #    plt.title('Total scattering aerosols 512')
    #    plt.colorbar()
    #if CProOK:
    #    plt.subplot(2,2,3)
    #    norm = colors.LogNorm(vmin=0.001,vmax=0.1)
    #    plt.pcolormesh(lats,altsC,t512C.T,cmap=mymap_sw,norm=norm)
    #    plt.ylim(yinf,ysup)
    #    plt.title('Total scattering clouds 512')
    #    plt.colorbar()
    tran = ccrs.PlateCarree(central_longitude=cm)
    proj = tran
    # geostationary plot
    ax = plt.subplot(2,1,2,projection = proj)
    #ax.set_extent(ext,ccrs.PlateCarree())
    ax.set_extent(box,ccrs.PlateCarree())
    if GEO & GEOOK:
        buf = ph.var['IR0'][subgg.corner[1]:subgg.corner[1]+subgg.box_biny,
                            subgg.corner[0]:subgg.corner[0]+subgg.box_binx]
        ax.imshow(buf,extent=ext,cmap='jet',clim=(190,300),
              transform=tran,interpolation='nearest',origin='lower')
    ax.coastlines('50m')
    # to make sure, plot it twice
    if sel>=8: lons1[lons1>180] -= 360 # unclear why it is needed
    if vortex_1st:
        ax.plot(lonv1-cm,latv1,'r',marker='+',markersize=11)
        ax.plot(lonv1+cm,latv1,'r',marker='+',markersize=11)
    if vortex_2nd:
        ax.plot(lonv2-cm,latv2,'r',marker='x',markersize=11)
        ax.plot(lonv2+cm,latv2,'r',marker='x',markersize=11)
    ax.plot(lons1+cm,lats1,'k')
    ax.plot(lons1-cm,lats1,'k')
    ax.set_xlim(xlims[0],xlims[1])
    gl = ax.gridlines(draw_labels=True,xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if GEO & GEOOK: plt.title(utc.strftime('Hima %Y-%m-%d %H:%M'))
    fig.suptitle(Cald[i]['fname'])
    if figsav:
        plt.savefig('plotAus-'+str(sel)+'.'+str(i)+'.png',dpi=300,bbox_inches='tight')
    if show: plt.show()
    plt.close(fig=fig)
    #if AProOK:
    #    hh.close()
    #    hdf.end()
    #if CProOK:
    #    hhC.close()
    #    hdfC.end()
    if L1OK:
        hh1.close()
        hdf1.end()
