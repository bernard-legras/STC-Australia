#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Exploring the calipso database for orbits that can interesect Koobor
Version for expedited data
In practice sel = 11 all the time
and the box is also changing in time
This is a real time script. It is not meant to be used beyond this context
and only the last part of the selection is to be used at any time.

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
import pickle,gzip
#import matplotlib.pyplot as plt
import os
from pyhdf.SD import SD
import glob
from astropy.time import Time, TimeDelta
import socket

# rootdir and main directories for aerosol profiles and L1 data
if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-Australia'
    dirAProf = '/dkol/data/CALIOP/05kmAPro_Exp.v3.40'
    dirL1 = '/dkol/data/CALIOP/CAL_LID_L1_Exp.v3.40'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-Australia'
    dirAProf = '/data/CALIOP/05kmAPro_Exp.v3.40'
    dirL1 = '/data/CALIOP/CAL_LID_L1_Exp.v3.40'
else:
    rootdir = '/home/b.legras/STC/STC-Australia'
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro_Exp.v3.40'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1_Exp.v3.40'

update = True
#for sel in [2,3,4]:
#for sel in range(4,6):
for sel in [11,]:

    if sel == 0:
        # south Pacific during mid-January
        box = [140,300,-60,-5]
        date0 = datetime(2020,1,15)
        nbday = 6
    elif sel ==1:
        # vicinity of Australia, early January
        box = [110,180,-50,-5]
        date0 = datetime(2019,12,31)
        nbday = 8
    elif sel ==2:
        # south Pacific during late January
        box = [140,300,-60,-5]
        date0 = datetime(2020,1,21)
        nbday = 9
    elif sel ==3:
        # south Pacific during mid January
        box = [140,300,-60,-5]
        date0 = datetime(2020,1,8)
        nbday = 7
    elif sel==4:
        # south Pacific during early January
        box = [140,300,-60,-5]
        date0 = datetime(2020,1,1)
        nbday = 7
    elif sel==5:
        # Indian Ocean during late January
        box = [140,300,-60,-5]
        date0 = datetime(2020,1,21)
        nbday = 6
    elif sel==6:
        # All tropical band during 5/2/2019 10/2/2019
        box = [None,None,-30,30]
        date0 = datetime(2019,2,5)
        nbday = 6
    elif sel==7:
        # All tropical band during 25/4/2019 30/4/2019
        box = [None,None,-30,30]
        date0 = datetime(2019,4,25)
        nbday = 6
    elif sel==8:
        # All tropical band during 25/7/2019 30/7/2019
        box = [None,None,-30,30]
        date0 = datetime(2019,7,25)
        nbday = 6
    elif sel==10:
        # Rising Spirit
        date0 = datetime(2020,2,1)
        box = [360-120,330,-70,-30]
        nbday = 5
    elif sel==11:
        # Rising Spirit
        date0 = datetime(2020,4,11)
        # The box has been modified several times
        no0shift = True
        box = [-100,180,-45,0]
        nbday = 3

    if update:
        with gzip.open('selCaliop_Exp_'+str(sel)+'.pkl','rb') as f:
            _,_,date00,nbday0,Cald = pickle.load(f)
        #for i in [45,46,47,49]:
        print ('read detection with ',len(Cald),' records')
        del Cald[106]
        idx = len(Cald)
        print(idx)

    else:
        Cald={}
        date00 = date0
        nbday0 = 0
        idx = 0

    # Day_Night_Flag
    # Browse dates
    for iday in range(nbday):
        date = date0 + timedelta(days=iday)
        # Generate names of daily directories
        dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
        dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
        print(dirday)
        # List the content of the daily aeorosol directory
        try:
            #fic = sorted(glob.glob(dirday+'/CAL_LID_L2_05kmAPro_Exp-Prov-V3-40.*.hdf'))
            fic = sorted(glob.glob(dirdayL1+'/CAL_LID_L1_Exp-Prov-V3-40.*.hdf'))
            print(len(fic))
            fic.reverse()
        except:
            print('missing day')
            continue
        # process all the half-orbits in the aerosol directory
        for i in range(len(fic)):
            # pop file
            fileL1 = fic.pop()
            print(fileL1)
            # skip day files
            #if 'ZD' in file: continue
            # open file
            try:
                hdfL1 = SD(fileL1)
            except :
                print('HDF4 Error -> giving up')
                continue
            # select orbits that intersect the box
            if no0shift:
                lonsL1 = hdfL1.select('Longitude').get()[:]
            else:
                lonsL1 = hdfL1.select('Longitude').get()[:] % 360
            latsL1 = hdfL1.select('Latitude').get()[:]
            dayniteL1 = hdfL1.select('Day_Night_Flag').get()[:]
            #sel1 = (lats < box[3]) & (lats > box[2]) & (daynite == 1) & (lons < box[1]) & (lons > box[0])
            sel1L1 = (latsL1 < box[3]) & (latsL1 > box[2]) & (lonsL1 < box[1]) & (lonsL1 > box[0])
            # for anomalous cases where the orbit does not have a tropical segment
            if np.sum(sel1L1) == 0:
                print('orbit with no selected band')
                continue
            tai = hdfL1.select('Profile_Time').get()[:]
            tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
            utc = tt.utc.datetime
            utcc = utc[sel1L1][0]+0.5*(utc[sel1L1][-1]-utc[sel1L1][0])
            fname = os.path.basename(fileL1)
            # extract rootname
            rootname = fname[26:-4]
            print(rootname)
            # generate full name of the corresponding L1 file
            file = os.path.join(dirday,'CAL_LID_L2_05kmAPro_Exp-Prov-V3-40.'+rootname+'.hdf')
            # open Apro file with escape in the anomalous case it is missing
            try:
                hdf = SD(file)
                # extract information form L1 file
                lats = hdf.select('Latitude').get()[:,1]
                if no0shift:
                    lons = hdf.select('Longitude').get()[:,1]
                else:
                    lons = hdf.select('Longitude').get()[:,1] % 360
                daynite = hdf.select('Day_Night_Flag').get()[:,0]
                #sel1L1 = (latsL1 < box[3]) & (latsL1 > box[2]) & (dayniteL1 == 1) & (lonsL1 < box[1]) & (lonsL1 > box[0])
                sel1 = (lats < box[3]) & (lats > box[2]) & (lons < box[1]) & (lons > box[0])
            except:
                print('Apro file cannot be found or opened')
                sel1 = None

            try:
                print(rootname,int(np.sum(sel1)),int(np.sum(sel1L1)),utcc,len(sel1),len(sel1L1))
            except:
                print(rootname,'None',int(np.sum(sel1L1)),utcc,'None',len(sel1L1))
            # generates a new record for this trace
            idx += 1
            print('idx',idx)
            Cald[idx] = {'fname':rootname,'utc':utcc,'sel1':sel1,'sel1L1':sel1L1}

    # store the dictionary of traces
    with gzip.open('selCaliop_Exp_'+str(sel)+'.pkl','wb') as f:
        pickle.dump([sel,box,date0,nbday+nbday0,Cald],f)


#%% Test of the size of the section

#for sel in range(3,9):
#    with gzip.open('selAmbae_'+str(sel)+'.pkl','rb') as f:
#        Cald = pickle.load(f)
#    print()
#    print('sel ',sel,len(Cald))
#    for i in range(len(Cald)):
#        print(Cald[i+1]['i2']-Cald[i+1]['i1'])













