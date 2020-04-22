#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exploring the CALIOP database for orbits that intersect an intersection area.
A number of selections are possible according toe the motion of the tracked structures.
The boxes are fixed. A new selection is to be defined if a new box is needed.

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

# Main irectories for aerosol profiles and L1 data
dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40'
dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'

update = True
#for sel in range(8):
#for sel in range(4,6):
for sel in [16,]:

    if sel == 0:
        # south Pacific during mid-January
        box = [130,300,-60,-5]
        date0 = datetime(2020,1,15)
        nbday = 7
    elif sel ==1:
        # vicinity of Australia, early January
        box = [110,180,-50,-5]
        date0 = datetime(2019,12,31)
        nbday = 8
    elif sel ==2:
        # south Pacific during late January
        box = [130,300,-60,-5]
        date0 = datetime(2020,1,22)
        nbday = 7
    elif sel ==3:
        # south Pacific during mid January
        box = [130,300,-60,-5]
        date0 = datetime(2020,1,8)
        nbday = 7
    elif sel==4:
        # south Pacific during early January
        box = [130,300,-60,-5]
        date0 = datetime(2020,1,1)
        nbday = 7
    elif sel==5:
        # Indian Ocean mid Jan
        box = [30,130,-60,-5]
        date0 = datetime(2020,1,8)
        nbday = 7
    elif sel==6:
        # Indian ocean mid Jan
        box = [30,130,-60,-5]
        date0 = datetime(2020,1,15)
        nbday = 7
    elif sel==7:
        # Indian Ocean late Jan
        box = [30,130,-60,-5]
        date0 = datetime(2020,1,22)
        nbday = 7
    elif sel==8:
        # Rising Spirit
        box = [360-120,359.999,-70,-30]
        date0 = datetime(2020,1,18)
        nbday = 14
    elif sel==9:
        # Rising Spirit
        box = [360-120,359.999,-70,-30]
        date0 = datetime(2020,2,1)
        nbday = 3
    elif sel==12:
        # south Pacific
        box = [130,330,-70,-30]
        #date0 = datetime(2019,12,31)
        date0 = datetime(2020,3,11)
        nbday = 3
    elif sel==13:
        # south Pacific
        box = [-50,160,-70,-30]
        #date0 = datetime(2019,12,31)
        date0 = datetime(2020,3,4)
        nbday = 1
    elif sel==14:
        # south Pacific
        box = [-110,100,-55,-15]
        #date0 = datetime(2019,12,31)
        date0 = datetime(2020,3,4)
        nbday = 5
    elif sel==15:
        # south Pacific
        box = [100,300,-55,-15]
        #date0 = datetime(2019,12,31)
        date0 = datetime(2020,3,23)
        nbday = 2
    elif sel==16:
        # south Atlantic+ Indian
        box = [-50,160,-40,0]
        date0 = datetime(2020,4,14)
        nbday = 1
    elif sel==17:
        # Antarctica
        box = [180,310,-85,-60]
        date0 = datetime(2020,1,4)
        nbday = 35

print('sel,date0,nbday,box')
print(sel,date0,nbday,box)

if update:
    with gzip.open('selCaliop_'+str(sel)+'.pkl','rb') as f:
        _,_,date00,nbday0,Cald = pickle.load(f)
    #for i in range(73,86):
    #    del Cald[i]
    idx = len(Cald)
    print ('read delection with ',idx,' records')
    print ('indexing from 1')
else:
    Cald = {}
    date00 = date0
    nbday0 = 0
    idx = 0

    # Day_Night_Flag

# Browse dates
print(nbday,' days to be processed')
for iday in range(nbday):
    date = date0 + timedelta(days=iday)
    # Generate names of daily directories
    dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
    dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    print(dirday)
    # List the content of the daily aeorosol directory
    #fic = sorted(glob.glob(dirday+'/CAL_LID_L2_05kmAPro-Prov-V3-40.*.hdf'))
    ficL1 = sorted(glob.glob(dirdayL1+'/CAL_LID_L1-ValStage1-V3-40.*.hdf'))
    print(len(ficL1))
    ficL1.reverse()
    # process all the half-orbits in the aerosol directory
    for i in range(len(ficL1)):
        # pop file
        file = ficL1.pop()
        print(file)
        # skip day files
        if ('ZD' in file) & (sel != 17) : continue
        # open file
        try:
            hdfL1 = SD(file)
        except :
            print('HDF4 Error -> giving up')
            continue
        # select orbits that intersect the box
        if sel not in [13,14]:
            lonsL1 = hdfL1.select('Longitude').get()[:] % 360
        else:
            lonsL1 = hdfL1.select('Longitude').get()[:]
        latsL1 = hdfL1.select('Latitude').get()[:]
        sel1L1 = (latsL1 < box[3]) & (latsL1 > box[2]) & (lonsL1 < box[1]) & (lonsL1 > box[0])
        # for anomalous cases where the orbit does not have a tropical segment
        if np.sum(sel1L1) == 0:
            print('orbit with no selected band')
            continue
        #lons1 = lons[sel1]
        #lons1 = lons1 % 360
        #if box[0] is not None:
        #    if np.any(lons1>box[1]) | np.any(lons1<box[0]):
        #        print('no trace in longitude box for this orbit')
        #        print(lons1[0],lons1[-1],lons1.min(),lons1.max())
        #        continue
        # extract information for selected half-orbit
        #lats1 = lats[sel1]
        #i1 = np.where(sel1)[0][0]
        #i2 = np.where(sel1)[0][-1]
        tai = hdfL1.select('Profile_Time').get()[:]
        tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
        utc = tt.utc.datetime
        utcc = utc[sel1L1][0]+0.5*(utc[sel1L1][-1]-utc[sel1L1][0])
        fname = os.path.basename(file)
        # extract rootname
        rootname = fname[27:-4]
        print(rootname,utcc,int(np.sum(sel1L1)),len(sel1L1))
        try:
            file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Prov-V3-40.'+rootname+'.hdf')
            hdf = SD(file)
            # extract information form L1 file
            lats = hdf.select('Latitude').get()[:,1]
            if sel not in [13,14]:
                lons = hdf.select('Longitude').get()[:,1] % 360
            else:
                lons = hdf.select('Longitude').get()[:,1]
            sel1 = (lats < box[3]) & (lats > box[2]) & (lons < box[1]) & (lons > box[0])
            print('Apro',np.sum(sel1),len(sel1))
        except:
            print(file)
            print('Apro file cannot be found or opened')
            sel1 = None

        # generates a new record for this trace
        idx += 1
        print('idx',idx)
        Cald[idx] = {'fname':rootname,'utc':utcc,'sel1':sel1,'sel1L1':sel1L1}

# store the dictionary of traces
with gzip.open('selCaliop_'+str(sel)+'.pkl','wb') as f:
    pickle.dump([sel,box,date00,nbday+nbday0,Cald],f)


#%% Test of the size of the section

#for sel in range(3,9):
#    with gzip.open('selAmbae_'+str(sel)+'.pkl','rb') as f:
#        Cald = pickle.load(f)
#    print()
#    print('sel ',sel,len(Cald))
#    for i in range(len(Cald)):
#        print(Cald[i+1]['i2']-Cald[i+1]['i1'])
