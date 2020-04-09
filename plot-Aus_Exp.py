# -*- coding: utf-8 -*-
"""
Shows the sequence of CALIOP sectiosn in the vicinity of the July 2018
eruption of the Ambae on 26 July at 21
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
from cartopy import feature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import argparse
import socket
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

Qs = 5.167*1.e-31
kbw = 1.0313

# rootdir and main directories for aerosol profiles and L1 data
if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-Australia'
    dirAProf = '/dkol/data/CALIOP/05kmAPro_Exp.v3.40'
    dirL1 = '/dkol/data/CALIOP/CAL_LID_L1_Exp.v3.40'
    dirCProf = '/dkol/data/CALIOP/05kmCPro_Exp.v3.40'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-Australia'
    dirAProf = '/data/CALIOP/05kmAPro_Exp.v3.40'
    dirL1 = '/data/CALIOP/CAL_LID_L1_Exp.v3.40'
    dirCProf = '/data/CALIOP/05kmACPro_Exp.v3.40'
else:
    rootdir = '/home/b.legras/STC/STC-Australia'
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro_Exp.v3.40'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1_Exp.v3.40'
    dirCProf = '/DATA/LIENS/CALIOP/05kmCPro.v3.40'

parser = argparse.ArgumentParser()
parser.add_argument("-s","--sel",type=int,choices=[0,1,2,3,4,5,6,7],help="selection")

firsti = 98

sel = 11
figsav = True
show = False
sav2hdf5 = False
vmax = 6
medfilter = True
widthfilter = 161

args = parser.parse_args()
if args.sel is not None: sel = args.sel

trac = pickle.load(open('Vortex-track.pkl','rb'))

def get_vortex_pos(date):
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
with gzip.open('selCaliop_Exp_'+str(sel)+'.pkl','rb') as f:
    sel,box,idate0,nbday,Cald = pickle.load(f)

# if sel in [0,2,3,4]:
#     GEO = True
#     gg = geosat.GeoGrid('HimFull')
#     subgg = gg.subgrid([130,220,-60,-5])
#     cm = 180
#     ext = [130-180,220-180,-60,-5]
#     xlims = [-40,120]
#     xlocs=(140,160,-180,-160,-140,-120,-100,-80,-60)
#     satload = geosat.Himawari
#     itm = 20
#     ysup = 25
# elif sel==1:
#     # Central Pacific during August and early September
#     GEO = True
#     gg = geosat.GeoGrid('HimFull')
#     subgg = gg.subgrid(box)
#     cm = 0
#     ext = box
#     xlims = [box[0],box[1]]
#     xlocs = None
#     satload = geosat.Himawari
#     itm = 20
#     ysup = 22
# elif sel in [5,6,7]:
#     # Indian Ocean during early September
#     GEO = True
#     gg = geosat.GeoGrid('MSG1Full')
#     subgg = gg.subgrid([box[0],110,box[2],box[3]])
#     cm = 0
#     ext = box
#     xlims = [box[0],110]
#     xlocs = None
#     satload = geosat.MSG1
#     itm = 15
#     ysup = 25
# elif sel in [8,9]:
#     GEO = False
#     ext = [-120,0,-70,-30]
#     xlims = [-120,0]
#     ysup = 28
#     cm = 0
#     xlocs = None
# elif sel==10:
#     GEO = False
#     ext = [-120,-30,-70,-30]
#     xlims = [-120,-30]
#     ysup = 28
#     cm = 0
#     xlocs = None
if sel==11:
    GEO = False
    ysup = 39
    yinf = 20
    itm = 20

if GEO:
    satload = geosat.Himawari
    gg = geosat.GeoGrid('HimFull')
    subgg = gg.subgrid(box)

for i in range(firsti,len(Cald)+1):
    # exclude day tracks
    if i in [1,2,4,12,13,14,15,16,17,18,27,29,30,31,32,33,43]:
        print('day trace, skipped')
        continue
    if i <= 43:
        cm = 180
        ext = [120-cm,330-cm,-70,-30]
        xlims = [120-cm,330-cm]
        xlocs=(140,160,-180,-160,-140,-120,-100,-80,-60,-40)
    elif i <= 52:
        ext = [0,200,-70,-30]
        xlims = [0,200]
        xlocs = None
        cm = 0
    elif i <= 63:
        ext = [0,200,-60,-20]
        xlims = [0,200]
        xlocs = None
        cm = 0
    elif i <= 73:
        ext = [-100,100,-60,-20]
        xlims = [-100,100]
        xlocs = None
        cm = 0
    else:
        cm = 180
        ext = [0-cm,200-cm,-55,-15]
        xlims = [0-cm,200-cm]
        xlocs = (0,20,40,60,80,100,120,140,160,180,200)


#for i in range(15,len(Cald)+1):
#for i in [1,2]:

    # generate date and daily directory name
    date = datetime.strptime(Cald[i]['fname'][:10],'%Y-%m-%d')
    dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
    dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    dirdayC = os.path.join(dirCProf,date.strftime('%Y/%Y_%m_%d'))

    try:
        [lonv,latv,zv] = get_vortex_pos(Cald[i]['utc'])
        print('found corresponding vortex position')
        vortex = True
    except:
        print('no vortex position')
        vortex = False

    try:
        file = os.path.join(dirday,'CAL_LID_L2_05kmAPro_Exp-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
        print('i',i)
        print(file)
        hdf = SD(file,SDC.READ)
        hh = HDF.HDF(file,HDF.HC.READ)
        sel1 = Cald[i]['sel1']
        # Correction to avoid day branch
        id = np.where(sel1[1:] ^ sel1[:-1])[0][1]
        sel1[id+10:] == False
        # Correction to avoid day branch
        daynite = hdf.select('Day_Night_Flag').get()[:,0]
        sel1 = np.array([x and y for (x,y) in zip(daynite==1,sel1)]).astype(np.bool)
        lons = hdf.select('Longitude').get()[sel1,1] % 360
        lats = hdf.select('Latitude').get()[sel1,1]
        meta = hh.vstart().attach('metadata')
        alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])

        a512 = np.ma.masked_less(hdf.select('Aerosol_Multiple_Scattering_Profile_532').get()[sel1,:],0)
        e512 = np.ma.masked_less(hdf.select('Extinction_Coefficient_532').get()[sel1,:],0)
        f512 = np.ma.masked_less(hdf.select('Aerosol_Layer_Fraction').get()[sel1,:],0)
        t512 = np.ma.masked_less(hdf.select('Total_Backscatter_Coefficient_532').get()[sel1,:],0)
        OK = True
    except:
        print('no Apro file')
        OK = False

    try:
        fileC = os.path.join(dirdayC,'CAL_LID_L2_05kmCPro_Exp-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
        hdfC = SD(fileC,SDC.READ)
        hhC = HDF.HDF(fileC,HDF.HC.READ)
        t512C = np.ma.masked_less(hdfC.select('Total_Backscatter_Coefficient_532').get()[sel1,:],0)
        metaC = hhC.vstart().attach('metadata')
        altsC = np.array(metaC.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
        COK = True
    except:
        print('no CPro file')
        COK = False

    try:
        fileL1 = os.path.join(dirdayL1,'CAL_LID_L1_Exp-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
        # access dataset
        hdf1 = SD(fileL1,SDC.READ)
        # access metadata
        hh1 = HDF.HDF(fileL1,HDF.HC.READ)
        # window
        sel1L1 = Cald[i]['sel1L1'][:,0]
        # Correction to avoid day branch
        daynite1 = hdf1.select('Day_Night_Flag').get()[:]
        sel1L1 = np.array([x and y for (x,y) in zip(daynite1==1,sel1L1)]).astype(np.bool)
        #much faster than sell1L1 = np.logical_and(sel1L1,daynite1)
        del daynite1
        # get data over the window
        t512L1 = hdf1.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:]
        # mask missing data
        t512L1 = np.ma.masked_less(t512L1,0)
        # get lon and lat over the window
        lons1 = hdf1.select('Longitude').get()[sel1L1].flatten() % 360
        lats1 = hdf1.select('Latitude').get()[sel1L1].flatten()
        # get molecular density number
        mnd1 = hdf1.select('Molecular_Number_Density').get()[sel1L1,:]
        # calculate log of molecular backscatter
        lbeta512_met = np.log(1000 * mnd1 * Qs / (kbw*8*np.pi/3))
        # get altitudes
        meta1 = hh1.vstart().attach('metadata')
        alts1 = np.array(meta1.read()[0][meta1.field('Lidar_Data_Altitudes')._idx])
        meta1 = hh1.vstart().attach('metadata')
        malts1 = np.array(meta1.read()[0][meta1.field('Met_Data_Altitudes')._idx])
        # interpolation of the molecular backscatter to lidar altitudes
        lbeta512_lid = np.empty(shape=t512L1.shape)
        for jy in range(len(lats1)):
            lbeta512_lid[jy,:] = np.interp(alts1,malts1[::-1],lbeta512_met[jy,::-1])
        # backscatter ratio
        if medfilter:
            sr512raw = t512L1/np.exp(lbeta512_lid)
            sr512= ss.medfilt(sr512raw,kernel_size=(widthfilter,1))
        else:
            sr512 = t512L1/np.exp(lbeta512_lid)
        if sav2hdf5:
            from scipy.io import savemat
            import deepdish as dd
            trace = {}
            trace['Total_Attenuated_Backscatter_532'] = t512L1
            trace['Longitude'] = lons1
            trace['Latitude'] = lats1
            trace['Molecular_Number_Density'] = mnd1
            trace['Lidar_Data_Altitudes'] = alts1
            trace['Met_Data_Altitudes'] = malts1
            trace['Profile_UTC_Time'] = hdf1.select('Profile_UTC_Time').get()[sel1L1].flatten()
            savemat('Trace_'+Cald[i]['fname']+'.mat',trace)
            dd.io.save('Trace_'+Cald[i]['fname']+'.hdf5',trace,compression='zlib')
        del t512L1
        del lbeta512_lid
        del mnd1
        L1OK = True

    except:
        print ('no L1 data')
        print (fileL1)
        L1OK = False

#%% Here we get the himawari scene which is the closest for this time

    if GEO:
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

    #%%
    #fig = plt.figure(figsize=(10,8))
    fig = plt.figure(figsize=(6,8))
    if L1OK:
       # plt.subplot(2,2,1)
        plt.subplot(2,1,1)
        #norm1 = colors.LogNorm(vmin=0.00001,vmax=0.1)
        plt.pcolormesh(lats1,alts1,sr512.T,cmap='gist_ncar',vmin=0,vmax=vmax)
        if vortex : plt.plot(latv,zv,'r',marker='+',markersize=11)
        plt.title('{:d} {:16} {:.2f} {:.2f}'.format(i,Cald[i]['utc'].strftime('%Y-%m-%dT%H:%M'),lons1[0],lons1[-1]))
        plt.ylim(yinf,ysup)
        plt.title('Scattering ratio 512')
        plt.colorbar()
    # plt.subplot(2,2,2)
    # norm = colors.LogNorm(vmin=0.00001,vmax=0.01)
    # if OK:
    #     plt.pcolormesh(lats,alts,t512.T,cmap=mymap_sw,norm=norm)
    #     plt.ylim(yinf,ysup)
    #     plt.title('Total scattering aerosols 512')
    #     plt.colorbar()
    # plt.subplot(2,2,3)
    # norm = colors.LogNorm(vmin=0.001,vmax=0.1)
    # if COK:
    #     plt.pcolormesh(lats,altsC,t512C.T,cmap=mymap_sw,norm=norm)
    #     plt.ylim(yinf,ysup)
    #     plt.title('Total scattering clouds 512')
    #     plt.colorbar()
    tran = ccrs.PlateCarree(central_longitude=cm)
    proj = tran
    # geostationary plot
    #ax = plt.subplot(2,2,4,projection = proj)
    ax = plt.subplot(2,1,2,projection = proj)
    ax.set_extent(ext,ccrs.PlateCarree())
    if GEO:
        buf = ph.var['IR0'][subgg.corner[1]:subgg.corner[1]+subgg.box_biny,
                            subgg.corner[0]:subgg.corner[0]+subgg.box_binx]
        ax.imshow(buf,extent=ext,cmap='jet',clim=(190,300),
              transform=tran,interpolation='nearest',origin='lower')
    ax.coastlines('50m')
    # to make sure, plot it twice
    if sel>=8: lons1[lons1>180] -= 360
    ax.plot(lons1+cm,lats1,'k')
    ax.plot(lons1-cm,lats1  ,'k')
    if vortex:
        ax.plot(lonv-cm,latv,'r',marker='+',markersize=11)
        ax.plot(lonv+cm,latv,'r',marker='+',markersize=11)
    plt.xlim(xlims[0],xlims[1])
    gl = ax.gridlines(draw_labels=True,xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if GEO: plt.title(utc.strftime('Hima %Y-%m-%d %H:%M'))
    fig.suptitle(Cald[i]['fname'])
    if figsav:
        plt.savefig('plotAus_Exp-'+str(sel)+'.'+str(i)+'.png',dpi=300,bbox_inches='tight')
    if show: plt.show()
    plt.close(fig=fig)
    if OK:
        hh.close()
        hdf.end()
    if COK:
        hhC.close()
        hdfC.end()
    if L1OK:
        hh1.close()
        hdf1.end()
