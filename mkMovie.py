#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a sequence of images centered at the location of the vortex and make an animated gif
from them

Note: always read the tempererature first as it is on the 137 levels at 6h and 18h.
vorticity and ozone are often missing at level 137 at 18h.

Created on 29 August 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
#from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import constants as cst
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D # yes it is used
from matplotlib import cm
from matplotlib.text import TextPath
import gzip,pickle
import socket
#import deepdish as dd
from os.path import join
from PIL import Image, ImageDraw, ImageFont
from os import chdir
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from ECMWF_N import ECMWF

figsav = False
figargs = dict(bbox_inches='tight',dpi=300)
tmp_dir = 'tmp_gif'

#%% Loading the version of the track that has also the maxvind
trac = pickle.load(open('Vortex-track-withwind.pkl','rb'))
#mean_track = pickle.load(open('Koobor-Ntrack-L1-mean.pkl','rb'))
#fctrac = pickle.load(open('Vortex-fctrack.pkl','rb'))
#for i in range(len(trac['dates'])):
#    print(i,trac['dates'][i],trac['lons'][i],trac['lats'][i],'{:2.1f}'.format(trac['z'][i]),trac['kz'][i])

# #%% Total displacement
trac['lats'] = np.array(trac['lats'])
trac['lons'] = np.array(trac['lons'])
# dy = trac['lats'][1:]-trac['lats'][:-1]
# dx = (trac['lons'][1:]-trac['lons'][:-1])*np.cos(np.deg2rad(0.5*(trac['lats'][1:]+trac['lats'][:-1])))
# # Correction for crossing Greenwich
# dx[149] = ((trac['lons'][150]-trac['lons'][149]%360))*np.cos(np.deg2rad(0.5*(trac['lats'][149]+trac['lats'][150])))
# ds = np.sqrt(dx**2+dy**2)
# print('total path ',(2*np.pi*6371/360)*np.sum(ds))

# avoid discontinuity before smoothing
trac['lons'][150:] -= 360

# smooth the trajectory and reset lon between 0 and 360
tlons = savgol_filter(trac['lons'],17,3) % 360
tlats = savgol_filter(trac['lats'],17,3)
#%% Plot of a selection of VO and ozone images at the level of the max vorticity
# make an image every day for the purpose of a movie

Firsti_2plot = 0
Lasti_2plot = len(trac['dates'])
#Firsti_2plot = 117
#Lasti_2plot = 118
#Lasti_2plot = 2
#Firsti_2plot = 50
#Lasti_2plot = 52
dx = 30
dy = 20
VO = True; VOh = False; VOv = True
VO = True; VOh = True; VOv = False
O3 = True; O3h = False
O3 = False; O3h = False
T = True; Tv = True
T = False; Tv = False
Z = False
SH = True

if VOv | Tv:
    T = True
    Z = True

for i in range(Firsti_2plot,Lasti_2plot):
    date1 = trac['dates'][i]
    print(date1)
    dat1 = ECMWF('OPZ',date1)
    varss = []
    if T:
        dat1._get_var('T')
        Tzon = np.mean(dat1.var['T'],axis=2)
        varss.append('Ta')
        dat1.var['Ta'] = dat1.var['T'] - Tzon[:,:,np.newaxis]
    if VO:
        dat1._get_var('VO')
        varss.append('VO')
    if O3:
        dat1._get_var('O3')
        O3zon = np.mean(dat1.var['O3'],axis=2)
        varss.append('O3')
        dat1.var['O3'] -= O3zon[:,:,np.newaxis]
    if VOv | Tv:
        dat1._mkp()
        dat1._mkz()
        varss.append('Z')

    lo1 = tlons[i] - dx
    lo2 = tlons[i] + dx
    la1 = tlats[i] - dy
    la2 = tlats[i] + dy
    if ((lo1 < 0) | (lo2 >= 360)):
        dats = dat1.shift2west(lon0=-179)
        if lo2 >= 360:
            lo1 -= 360
            lo2 -= 360
    else:
        dats = dat1
    datr = dats.extract(latRange=(la1,la2),lonRange=(lo1,lo2),varss=varss)
    if SH: datsh = dat1.extract(latRange=(-80,0),varss=varss)
    del dats
    del dat1

    if O3h:
        ax1 = datr.show('O3',lev = trac['kz'][i],show=False,figsize=(7,4),scale=1.e6,clim=(-3,0.5),txt='')
        ax1.set_title('Ozone mr anomaly (mg/kg) level '+str(int(trac['p'][i]/100))+' hPa'+
                  date1.strftime(' %Y-%m-%d:%H UTC'))
        plt.savefig(join(tmp_dir,'imO3_'+str(i)+'.jpg'))
        #plt.savefig(join(tmp_dir,'imO3_'+str(i)+'.png'))
        plt.show()

    if VOh:
        ax2 = datr.show('VO',lev = trac['kz'][i],show=False,figsize=(7,4),scale=1.e5,clim=(-1,8),txt = '')
        ax2.set_title(u'Vorticity (10$^{-5}$ s$^{-1}$) level '+str(int(trac['p'][i]/100))+' hPa'+
                  date1.strftime(' %Y-%m-%d:%H UTC'))
        plt.savefig(join(tmp_dir,'imVO_'+str(i)+'.jpg'),bbox_inches='tight')
        #plt.savefig(join(tmp_dir,'imVO_'+str(i)+'.png'))
        plt.show()
        if SH:
            ax2 = datsh.show('VO',lev = trac['kz'][i],show=False,scale=1.e5,clim=(-1,8),txt = '')
            ax2.set_title(u'Vorticity (10$^{-5}$ s$^{-1}$) level '+str(int(trac['p'][i]/100))+' hPa'+
                  date1.strftime(' %Y-%m-%d:%H UTC'))
            plt.savefig(join(tmp_dir,'imVO_sh_'+str(i)+'.jpg'),bbox_inches='tight')
            #plt.savefig(join(tmp_dir,'imVO_sh_'+str(i)+'.png'))
            plt.show()
    if VOv:
        ax3 = datr.chartlonz('VO',trac['lats'][i],show=False,scale=1.e5,clim=(-1,8),levs=(20,80),txt='')
        ax3.set_title(u'Vorticity (10$^{-5}$ s$^{-1}$) lat '+str(-int(trac['lats'][i]))+
                      trac['dates'][i].strftime('S %Y-%m-%d:%H UTC'))
        plt.savefig(join(tmp_dir,'imVO_vert_'+str(i)+'.jpg'))
        #plt.savefig(join(tmp_dir,'imVO_vert_'+str(i)+'.png'))
        plt.show()
    if Tv:
        ax4 = datr.chartlonz('Ta',trac['lats'][i],show=False,clim=(-10,10),levs=(20,80),txt='')
        ax4.set_title(u'Temperature anomaly (K) lat '+str(-int(trac['lats'][i]))+
                      trac['dates'][i].strftime('S %Y-%m-%d:%H UTC'))
        plt.savefig(join(tmp_dir,'imT_vert_'+str(i)+'.jpg'))
        #plt.savefig(join(tmp_dir,'imT_vert_'+str(i)+'.png'))
        plt.show()


# #%%
# import imageio
# images = []
# for i in range(Firsti_2plot,Lasti_2plot):
#     images.append(imageio.imread(join(tmp_dir,'imO3_'+str(i)+'.jpg')))
# imageio.mimsave('movie_O3.gif', images)
# images = []
# for i in range(Firsti_2plot,Lasti_2plot):
#     images.append(imageio.imread(join(tmp_dir,'imVO_'+str(i)+'.jpg')))
# imageio.mimsave('movie_VO.gif', images)


