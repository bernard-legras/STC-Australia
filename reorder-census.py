#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Comparison of the census of the positions and altitudes of the vortex
In principle should be used to combine the contribution from several contributors
but is now limited to BL

Created on Wed Feb 26 17:03:30 2020

@author: legras
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

initials = ['BL',]

file = 'Koobor-Ntrack-L1-BL.pkl'
trackSpt = pickle.load(open(file,'rb'))

mean_track = {}
nv = len(trackSpt)
mean_track['centroid_alt'] = np.zeros(nv)
mean_track['centroid_lat'] = np.zeros(nv)
mean_track['centroid_lon'] = np.zeros(nv)
mean_track['top_alt'] = np.zeros(nv)
mean_track['bot_alt'] = np.zeros(nv)
mean_track['south_lat'] = np.zeros(nv)
mean_track['north_lat'] = np.zeros(nv)
mean_track['mid_alt'] = np.zeros(nv)
mean_track['mid_lat'] = np.zeros(nv)
mean_track['mid_lon'] = np.zeros(nv)
fig, ax = plt.subplots(2,2,sharex=False,sharey=False)
for init in initials:
    file = 'Koobor-Ntrack-L1-'+init+'.pkl'
    trackSpt = pickle.load(open(file,'rb'))
    # correctifs
    trackSpt[34][8][2] -= 360
    for ix in range(5,10): trackSpt[35][ix][2] -= 360
    tSpt = []
    centroid_alt = []
    centroid_lat = []
    centroid_lon = []
    mid_alt = []
    mid_lat = []
    mid_lon = []
    top_alt = []
    bot_alt = []
    south_lat = []
    north_lat = []
    for i in range(len(trackSpt)):
        [sel,i,qual,fname,date,bot,mid,top,south,north] = trackSpt[i]
        # test for errors in pointing
        if bot[0] > np.min([top[0],mid[0],south[0],north[0]]):
            print('Error 1',i)
        if top[0] < np.max([bot[0],mid[0],south[0],north[0]]):
            print('Error 2',i)
        if south[1] > np.min([bot[1],top[1],mid[1],north[1]]):
            print('Error 3',i)
        if north[1] < np.max([bot[1],top[1],mid[1],south[1]]):
            print('Error 4',i)
        tSpt.append(date)
        if (date > datetime(2020,2,24)) & (date < datetime(2020,3,10)):
            f = lambda x: (x+180) % 360 - 180
            south[2] = f(south[2])
            north[2] = f(north[2])
            bot[2] = f(bot[2])
            top[2] = f(top[2])
        aa = [0.25*(south[j]+north[j]+bot[j]+top[j]) for j in range(3)]
        centroid_alt.append(aa[0])
        centroid_lat.append(aa[1])
        centroid_lon.append(aa[2]%360)
        mid_alt.append(mid[0])
        mid_lat.append(mid[1])
        mid_lon.append(mid[2]%360)
        top_alt.append(top[0])
        bot_alt.append(bot[0])
        south_lat.append(south[1])
        north_lat.append(north[1])

    ax[0,0].plot(mid_lon,mid_lat,centroid_lon,centroid_lat)
    ax[0,1].plot(tSpt,mid_alt,tSpt,centroid_alt)
    ax[1,0].plot(tSpt,top_alt,tSpt,bot_alt)
    ax[1,1].plot(tSpt,south_lat,tSpt,north_lat)
    #latSpt
    mean_track['centroid_alt'] += np.array(centroid_alt)
    mean_track['centroid_lat'] += np.array(centroid_lat)
    mean_track['centroid_lon'] += np.array(centroid_lon)
    mean_track['mid_alt'] += np.array(mid_alt)
    mean_track['mid_lat'] += np.array(mid_lat)
    mean_track['mid_lon'] += np.array(mid_lon)
    mean_track['top_alt'] += np.array(top_alt)
    mean_track['bot_alt'] += np.array(bot_alt)
    mean_track['south_lat'] += np.array(south_lat)
    mean_track['north_lat'] += np.array(north_lat)

ax[0,0].legend(('BL mid','BL cen'))
ax[0,0].set_xlabel('longitude')
ax[0,0].set_ylabel('latitude')
ax[0,0].set_title('longitude')
ax[0,1].set_title('altitude')
ax[1,0].set_title('top & bootom')
ax[1,1].set_title('south and north')
fig.autofmt_xdate()
plt.show()

for var in ['centroid_alt','centroid_lat','centroid_lon','mid_alt','mid_lat','mid_lon',
            'bot_alt','top_alt','south_lat','north_lat']:
    mean_track[var] /= len(initials)

mean_track['dates'] = tSpt

#%%
print('date, hour, cen lon, south lat, cen lat, north lat (°), bot, cen, top (km)')
for i in range(len(mean_track['dates'])):
    print(mean_track['dates'][i].strftime('%Y-%m-%d %H:%M'),
          ' {:2.1f}  \t{:2.1f} \t{:2.1f} \t{:2.1f} \t{:2.1f} \t{:2.1f} \t{:2.1f}'.format(
          mean_track['centroid_lon'][i],
          mean_track['south_lat'][i],mean_track['centroid_lat'][i],mean_track['north_lat'][i],
          mean_track['bot_alt'][i],mean_track['centroid_alt'][i],mean_track['top_alt'][i]))

pickle.dump(mean_track,open('Koobor-Ntrack-L1-mean.pkl','wb'))


