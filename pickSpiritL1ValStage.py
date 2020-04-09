#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a list of positions of the Spirit center from manual pointing of CALIOP
images

This program must be followed by reorder_census

Created on Wed Feb 12 01:09:16 2020

@author: Bernard Legras
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle,gzip
import argparse
from datetime import datetime, timedelta
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS, V
import matplotlib.colors as colors
import os
import socket
import sys
import scipy.signal as ss
#from scipy.interpolate import RegularGridInterpolator

parser = argparse.ArgumentParser()
parser.add_argument("-s","--step", type=int,choices=[1,2,3,4,5],help="step")

global top,bot,mid,south,north,breaker

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
cmap = mymap_sw
cmap = 'gist_ncar'

std12 = False
std13 = False
std14 = False
exp11 = False
std16 = False

step = 4

medfilt = 81
vmax = 6
yinf = 20
ysup = 39
# config to make the filtered sr512 and store them into tmp
#komput = True
#point = False
#srsave = True
# config to read the calculated sr512 anddisplay them for pointing
komput = False
srsave = False
point = True
# config to calculate without saving
komput = True
srsave = True
point = True

args = parser.parse_args()
if args.step is not None: step = args.step

initials = 'BL'
outfile = 'Koobor-Ntrack-L1-'+initials+'.pkl'

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

exi = False
update = False
if step == 1:
    try:
        buf = pickle.load(open(outfile,'rb'))
        print('the output file should not exist at this step')
        print('erase it and relaunch')
        exi = True
    except:
        std12 = True
        update = False
elif step == 5:
    std12 = True
    if point: update = True
elif step == 2:
    std13 = True
    if point: update = True
elif step == 3:
    std14 = True
    if point: update = True
elif step == 4:
    exp11 = True
    if point: update = True
elif step == 6:
    std16 = True
    if point: update = True

    #std = False
    #exp = True
# =============================================================================
# elif step ==3:
#     if socket.gethostname() not in ['gort','satie']:
#         print('this step should not run on ICARE')
#         exi = True
#     update = True
#     std = False
#     exp = True
# =============================================================================

if exi:
    sys.exit(0)

prepend = False
if std12:
    sel='12'
    dirL1 = dirL1_Std
    stream = '-ValStage1'
    il = {16:'a',24:'a',40:'a',55:'a',64:'a',81:'a',95:'a',121:'a',138:'a',154:'a',
      162:'a',180:'a',195:'a',211:'a',228:'a',237:'a',255:'a',264:'a',273:'a',
      316:'a',326:'a',335:'a',344:'a',353:'a',362:'a',372:'a',399:'a',408:'a',
      417:'a',426:'a',435:'a'}
    il = {55:'a',40:'a',24:'a',16:'a'}
    prepend = True
if std13:
    sel='13'
    dirL1 = dirL1_Std
    stream = '-ValStage1'
    il = {439:'a',449:'a',458:'a',468:'a',478:'a',515:'b'}
if std14:
    sel='14'
    dirL1 = dirL1_Std
    stream = '-ValStage1'
    il = {1:'a',11:'a',20:'a'}
if std16:
    sel='16'
    dirL1 = dirL1_Std
    stream = '-ValStage1'
    il = {3:'a',15:'a',22:'a',36:'a',44:'a',51:'a',58:'a'}
if exp11:
    dirL1 = dirL1_Exp
    sel='Exp_11'
    stream = '_Exp-Prov'
    il = {80:'b',81:'a',82:'a',83:'a',84:'a',85:'a',86:'a',88:'a',89:'a'}
    il = {82:'a',83:'a',84:'a',85:'a',86:'a',88:'a',89:'a'}
    il = {90:'a',91:'a'}
    il = {96:'a',97:'a'}
    medfilt = 81
    vmax = 6
# =============================================================================
# if exp:
#     sel='Exp_11'
#     dirL1 = dirL1_Exp
#     stream = '_Exp-Prov'
#     il1 = {9:'a',25:'a',41:'a'}
#     il2 = {44:'a',46:'a',49:'a',50:'a',53:'a'}
#     if socket.gethostname() in ['gort','satie']:
#         # take local data
#         il = il2
#     else:
#         # date data on ICARE
#         il = il1
# =============================================================================

with gzip.open('selCaliop_'+sel+'.pkl','rb') as f:
    _,box,idate0,nbday,Cald = pickle.load(f)

if update:
    track = pickle.load(open(outfile,'rb'))
    print('track read with ',len(track),' records')
else:
    track = []
    print('new track')
#
def on_click(event):
    pos = np.where(lats<event.xdata)[0][0]
    tt = utc[pos]
    dd = datetime.strptime(str(int(tt+20000000)),'%Y%m%d')+timedelta(days=tt-int(tt))
    print(dd,'latitude=%f longitude=% altitude=%f'%(event.xdata,lons[pos],event.ydata))
    #print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #      ('double' if event.dblclick else 'single', event.button,
    #       event.x, event.y, event.xdata, event.ydata))
def on_key(event):
    global top,bot,mid,south,north,breaker
    pos = np.where(lats<event.xdata)[0][0]
    tt = utc[pos]
    dd = datetime.strptime(str(int(tt+20000000)),'%Y%m%d')+timedelta(days=tt-int(tt))
    print('you pressed ',event.key,event.xdata,lons[pos],event.ydata)
    if event.key == '8':
        print('recording top height + pos')
        top = [event.ydata, event.xdata, lons[pos]]
    if event.key == '2':
        print('recording bottom height + pos')
        bot = [event.ydata, event.xdata, lons[pos]]
    if event.key == '5':
        print('recording middle height + pos')
        mid = [event.ydata, event.xdata, lons[pos]]
    if event.key == '4':
        print('recording southest envelop + pos')
        south = [event.ydata, event.xdata, lons[pos]]
    if event.key == '6':
        print('recording northest envelop + pos')
        north = [event.ydata, event.xdata, lons[pos]]
    if event.key == '+':
        if None in [top,bot,mid,south,north]:
            print('missing data ',top,bot,mid,south,north)
        else:
            print('recording whole trac')
            if prepend:
                track.insert(0,[sel,i,il[i],Cald[i]['fname'],dd,bot,mid,top,south,north])
            else:
                track.append([sel,i,il[i],Cald[i]['fname'],dd,bot,mid,top,south,north])
            plt.close('all')
    if event.key == 'e':
        print('breaking and storing the results')
        breaker = True
        plt.close('all')

Qs = 5.167*1.e-31
kbw = 1.0313
breaker = False

for i in il:
    date = datetime.strptime(Cald[i]['fname'][:10],'%Y-%m-%d')
    dirday = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    file = os.path.join(dirday,'CAL_LID_L1'+stream+'-V3-40.'+Cald[i]['fname']+'.hdf')
    print('i',i)
    print(file)
    #if komput:
    hdf = SD(file,SDC.READ)
    hh = HDF.HDF(file,HDF.HC.READ)
    sel1L1 = Cald[i]['sel1L1'][:,0]
    daynite = hdf.select('Day_Night_Flag').get()[:]
    sel1L1 = np.array([x and y for (x,y) in zip(daynite==1,sel1L1)]).astype(np.bool)
    lons = hdf.select('Longitude').get()[sel1L1].flatten() % 360
    lats = hdf.select('Latitude').get()[sel1L1].flatten()
    utc = hdf.select('Profile_UTC_Time').get()[sel1L1].flatten()
    if komput:
        mnd = hdf.select('Molecular_Number_Density').get()[sel1L1,:]
        lbeta512_met = np.log(1000 * mnd * Qs / (kbw*8*np.pi/3))
        t512 = np.ma.masked_less(hdf.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:],0)
        meta = hh.vstart().attach('metadata')
        alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
        meta = hh.vstart().attach('metadata')
        malts = np.array(meta.read()[0][meta.field('Met_Data_Altitudes')._idx])
        # calculation of the molecular backscatter
        lbeta512_lid = np.empty(shape=t512.shape)
        for jy in range(len(lats)):
            lbeta512_lid[jy,:] = np.interp(alts,malts[::-1],lbeta512_met[jy,::-1])
        sr512raw = t512/np.exp(lbeta512_lid)
        # filtering to remove noise
        sr512 = ss.medfilt(sr512raw,kernel_size=(medfilt,1))
        if srsave:
            with gzip.open('tmp/sr512-'+sel+'-'+str(i)+'pkl','wb') as f:
                pickle.dump([sr512,lats,alts],f)
    if point:
        if komput == False:
            with gzip.open('tmp/sr512-'+sel+'-'+str(i)+'pkl','rb') as f:
                [sr512,lats,alts] = pickle.load(f)
        fig, ax = plt.subplots()
        im=ax.pcolormesh(lats,alts,sr512.T,cmap=cmap,vmin=0,vmax=vmax)
        ax.set_ylim(yinf,ysup)
        ax.set_title('Aerosol scattering ratio 512')
        plt.colorbar(im)
        top = None; bot = None; mid = None; south = None; north = None
        cid1 = fig.canvas.mpl_connect('button_press_event', on_click)
        cid2 = fig.canvas.mpl_connect('key_press_event', on_key)
        plt.show()
        if breaker: break

if point: pickle.dump(track,open(outfile,'wb'))

if not breaker:
    if step == 1:
        print('CONGRATULATIONS, YOU HAVE COMPLETED THE FIRST STEP')
    elif step==2:
        print('CONGRATULATIONS, YOU HAVE COMPLETED THE SECOND STEP')
    elif step==3:
        print('CONGRATULATIONS, YOU HAVE COMPLETED THE THIRD STEP')
    elif step==4:
        print('GOD BLESS YOU, YOU HAVE COMPLETED THE FOURTH STEP')
# tSpt = []
# latSpt = []
# lonSpt = []
# altSpt = []
# for i in range(len(track)):
#     tSpt.append(track[i][4])
#     latSpt.append(track[i][5])
#     lonSpt.append(track[i][6])
#     altSpt.append(track[i][7])

# plt.plot(lonSpt,latSpt)
# plt.show()
# plt.plot(tSpt,altSpt)
# plt.show()
