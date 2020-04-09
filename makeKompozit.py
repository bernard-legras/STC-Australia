#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make a compozit 

Created on Sat Feb  8 12:27:14 2020

@author: Bernard Legras
"""
from datetime import datetime, timedelta
from ECMWF_N import ECMWF
import numpy as np
#from zISA import zISA
import constants as cst
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # yes it is used
from matplotlib import cm
#from matplotlib.text import TextPath
import gzip,pickle
import socket
#import deepdish as dd
from os.path import join
#from PIL import Image, ImageDraw, ImageFont
#from os import chdir

if 'gort' == socket.gethostname():
    rootdir = '/dkol/data/STC/STC-Australia'
elif 'satie' in socket.gethostname():
    rootdir = '/data/STC/STC-Australia'
    
figsav = True
figargs = dict(bbox_inches='tight',dpi=300)
    
# =============================================================================
# #%%
# with gzip.open('OPZ-extract-1.pkl','rb') as f:
#     dats = pickle.load(f)
# with gzip.open('OPZ-extract-2.pkl','rb') as f:
#     dats.update(pickle.load(f))
# with gzip.open('OPZ-Z-1.pkl','rb') as f:
#     datz = pickle.load(f)
# with gzip.open('OPZ-Z-2.pkl','rb') as f:
#     datz.update(pickle.load(f))
# print(len(dats),len(datz))
# =============================================================================

#%%    
trac = pickle.load(open('Vortex-track.pkl','rb'))
   
#%% print the positions as a function of time
for i in range(len(trac['dates'])):
    # kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
    print(i,trac['dates'][i],trac['lons'][i],trac['lats'][i],'{:2.1f}'.format(trac['z'][i]),trac['kz'][i])
    
# =============================================================================
# #%% Kompozit
# # The composit is made from 44 (29Jan) to 65 (8 Feb)
# # and from 62 (4 Feb) to 87 (16 Feb)
# # and from 82 (14 Feb) to 103 (24 Feb)   
# kompozit={}
# jdy = 6
# idx = 8
# kdz = 8
# istart = 42
# iend = 100
# kompozit['VO'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
# kompozit['O3ano'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
# kompozit['Tano'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
# for i in range(42,100):
#     print(i)
#     ix = np.where(dats[i].attr['lons']>=trac['lons'][i])[0][0]
#     jy = np.where(dats[i].attr['lats']>=trac['lats'][i])[0][0]
#     kz = np.where(dats[i].attr['zscale']<=trac['alts'][i])[0][0]
#     kompozit['VO'] += dats[i].var['VO'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1]
#     kompozit['O3ano'] += dats[i].var['O3'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1] \
#                        - dats[i].var['O3'][kz,jy,ix]
#     kompozit['Tano'] += dats[i].var['T'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1] \
#                       - dats[i].var['T'][kz,jy,ix]
# ns = iend-istart
# kompozit['VO'] /= ns
# kompozit['O3ano'] /= ns
# kompozit['Tano'] /= ns
# =============================================================================
#%% Alternate composite that removes the mean temperature and ozone profile at the same latitude 
# Because of this average removal we need to use full files instead od dats
# The composit is made from 44 (29Jan) to 65 (8 Feb)
# and from 62 (4 Feb) to 87 (16 Feb)
# and from 82 (14 Feb) to 103 (24 Feb)
# and from 42 (25 Jan) to 99 (22 Feb)
kompozit={}
jdy = 10
idx = 12
kdz = 12
istart = 42
iend = 100
kompozit['VO'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
kompozit['O3ano'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
kompozit['Tano'] = np.zeros(shape=(1+2*kdz,1+2*jdy,1+2*idx))
for i in range(istart,iend):
    print(i)
    dat = ECMWF('OPZ',trac['dates'][i])
    dat._get_var('VO')
    dat._get_var('T')
    dat._get_var('O3')
    ix = np.where(dat.attr['lons']>=trac['lons'][i])[0][0]
    jy = np.where(dat.attr['lats']>=trac['lats'][i])[0][0]
    kz = trac['kz'][i]
    kompozit['VO'] += dat.var['VO'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1]
    meanT = np.mean(dat.var['T'],axis=2)
    meanO3 = np.mean(dat.var['O3'],axis=2)
    kompozit['O3ano'] += dat.var['O3'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1] \
                       - meanO3[kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,np.newaxis]
    kompozit['Tano'] += dat.var['T'][kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,ix-idx:ix+idx+1] \
                      - meanT[kz-kdz:kz+kdz+1,jy-jdy:jy+jdy+1,np.newaxis]
ns = iend-istart
kompozit['VO'] /= ns
kompozit['O3ano'] /= ns
kompozit['Tano'] /= ns
kompozit['irange'] = (istart,iend)
kompozit['dxyz'] = (idx,jdy,kdz)
with gzip.open('kompozit_N.pkl','wb') as f:
    pickle.dump(kompozit,f)
 
#%%
# Plot of the kompozit as 3D plots 
# Composit vorticity
fig = plt.figure(figsize=(7,7))
ax = fig.gca(projection = '3d')
X1 = np.arange(-idx,idx+1)
Y1 = np.arange(-jdy,jdy+1)
X,Y = np.meshgrid(X1,Y1)
Z = np.zeros_like(X)
tt = X / np.max(X)
buf = kompozit['VO'] - np.min(kompozit['VO'])
buf /= np.max(buf)                              
for k in np.arange(kdz+4,kdz-4-1,-4):
    ax.plot_surface(X,Y,Z-k+kdz,rstride=1,cstride=1,facecolors = cm.jet(buf[k,...]))
ax.set_xlabel('longitude (degree)')
ax.set_ylabel('latitude (degree)')
ax.set_zlabel('altitude (level)')
plt.title('composite normalized vorticity')
if figsav:
    plt.savefig(join('figs','kompo3D_N_vorticity_14Jan-22Feb.png'),**figargs)
plt.show()
#%%
# Composit ozone
fig = plt.figure(figsize=(7,7))
ax = fig.gca(projection = '3d')
X1 = np.arange(-idx,idx+1)
Y1 = np.arange(-jdy,jdy+1)
X,Y = np.meshgrid(X1,Y1)
Z = np.zeros_like(X)
tt = X / np.max(X)
buf = kompozit['O3ano'] - np.min(kompozit['O3ano'])
buf /= np.max(buf)                              
for k in np.arange(kdz+4,kdz-4-1,-4):
    ax.plot_surface(X,Y,Z-k+kdz,rstride=1,cstride=1,facecolors = cm.jet(buf[k,...]))
ax.set_xlabel('longitude (degree)')
ax.set_ylabel('latitude (degree)')
ax.set_zlabel('altitude (level)')
plt.title('composite normalized ozone')
if figsav:
    plt.savefig(join('figs','kompo3D_N_ozone_14Jan-22Feb.png'),**figargs)
plt.show()
#%%
# Composit temperature
fig = plt.figure(figsize=(7,7))
ax = fig.gca(projection = '3d')
X1 = np.arange(-idx,idx+1)
Y1 = np.arange(-jdy,jdy+1)
X,Y = np.meshgrid(X1,Y1)
Z = np.zeros_like(X)
tt = X / np.max(X)
buf = kompozit['Tano'] - np.min(kompozit['Tano'])
buf /= np.max(buf)                              
for k in np.arange(kdz+4,kdz-4-1,-4):
    ax.plot_surface(X,Y,Z-k+kdz,rstride=1,cstride=1,facecolors = cm.jet(buf[k,...]))
ax.set_xlabel('longitude (degree)')
ax.set_ylabel('latitude (degree)')
ax.set_zlabel('altitude (level)')
plt.title('composite normalized temperature')
if figsav:
    plt.savefig(join('figs','kompo3D_N_temperature_14Jan-22Feb.png'),**figargs)
plt.show()

#%% make a 2D horizontal plot of the kompozit at central level and +-4
#ff,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,12))
plt.figure(figsize=(14,13))
imargs = dict(origin='lower',aspect=1,interpolation='nearest',cmap='jet',extent=(-idx*0.7,idx*0.7,-jdy*1.1,jdy*1.1))
i=1
for k in np.arange(kdz+4,kdz-4-1,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(1.e5*kompozit['VO'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('VO level '+str(kdz-k)+u' ($10^{-5}~s^{-1}$)')
    plt.xlabel('longitude (100 km)')
    plt.ylabel('latitude (100 km)')
    i += 1
for k in np.arange(kdz+4,kdz-4-1,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(1.e6*kompozit['O3ano'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('O3 ano level '+str(kdz-k)+u' ($10^{-6}~kg~kg^{-1}$)')
    plt.xlabel('longitude (100 km)')
    plt.ylabel('latitude (100 km)')
    i += 1
for k in np.arange(kdz+4,kdz-4-1,-4):
    plt.subplot(3,3,i)
    im = plt.imshow(kompozit['Tano'][k,...],**imargs)
    plt.colorbar(im,)
    plt.title('T anomaly level '+str(kdz-k)+' (K)')
    plt.xlabel('longitude (100 km)')
    plt.ylabel('latitude (100 km)')
    i += 1
if figsav:
    plt.savefig(join('figs','kompo_VO3T_N_14Jan-22Feb.png'),**figargs)
plt.show()

#%% make a 2D vertical lon-alt section in the central longitude plane
fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(14,5),sharey=True)
fig.subplots_adjust(left=0.02, bottom=0.12, right=0.95, top=0.78, wspace=0.05)
central = jdy+1
imargs = dict(aspect=1.4,interpolation='nearest',extent=(-0.7*idx,0.7*idx,-0.5*kdz,0.5*kdz))
im1 = ax1.imshow(1.e5*kompozit['VO'][:,central,:],cmap='jet',**imargs)
ax1.set_title('VO vertical section'+u' ($10^{-5}~s^{-1}$)',fontsize=18)
ax1.set_xlabel('longitude (100 km)',fontsize=16)
ax1.set_ylabel('altitude (km)',fontsize=16)
plt.colorbar(im1,ax=ax1)
im2 = ax2.imshow(1.e6*kompozit['O3ano'][:,central,:],cmap='jet',**imargs)
ax2.set_title('O3 anomaly'+u' ($10^{-6}~kg~kg^{-1}$)',fontsize=18)
ax2.set_xlabel('longitude (100 km)',fontsize=16)
#ax2.set_ylabel('altitude (km)')
plt.colorbar(im2,ax=ax2)
from matplotlib.colors import ListedColormap
fid=open('RedWhiteBluecolorscale.txt')
RWBmap = ListedColormap(np.array([np.genfromtxt(x.rstrip('\n').split('\t')) \
                               for x in fid.readlines()])/65535,'RWB')
im3 = ax3.imshow(kompozit['Tano'][:,central,:],vmin=-5,vmax=5,cmap=RWBmap,**imargs)
ax3.set_title('T anomaly (K)',fontsize=18)
ax3.set_xlabel('longitude (100 km)',fontsize=16)
#ax3.set_ylabel('altitude (km)')
plt.colorbar(im3,ax=ax3)

if figsav:
    plt.savefig(join('figs','kompo_vsect_VO3T_N_14Jan-22Feb.png'),**figargs)
plt.show()
plt.show()

