#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Generate the movies from the images produced by mkMovie.py

Created on Sun Aug 30 12:11:41 2020

@author: Bernard Legras
"""
import os
import imageio
from subprocess import call
#from pygifsicle import optimize

Firsti_2plot = 0
Lasti_2plot = 180
tmp_dir = 'tmp_gif2'

os.chdir(tmp_dir)

vort = True
version = 2

listi={12:[16,24,32,40,55,64,72,81,95,104,112,121,138,154,162,180,195,211,228,237,255,264,273,316,326,
         335,344,353,362,372,399,408,417,426,435],
       13:[439,449,458,468,478,515],14:[1,11,20,30],15:[61,67,76,85,94,103,112,121,129],
       16:[15,22,36,51],11:[96,97]}


#%%
if version == 1:
    sels = [12,13,14]
elif version == 2:
    sels = [12,13,14,15,16,11]
elif version == 3:
    sels = [15,16,11]

if vort == True:
    images = []
    for sel in sels:
        for i in listi[sel]:
            images.append(imageio.imread('VOSLON_CALIOP_BS/VOSLON_CALIOP_BS_{}_{}.png'.format(sel,i)))
    vv = str(version)
    imageio.mimsave('movie_voslon_caliop_'+vv+'.jpg.gif', images, fps=2)
    #imageio.mimsave('movie_voslon_caliop_'+vv+'.png.mp4', images, fps=2)
    call(["gifsicle","-O","--colors","256","-o","movie_voslon_caliop_"+vv+"_opt.jpg.gif","movie_voslon_caliop_"+vv+".jpg.gif"])

if vort == False:
    images = []
    for sel in sels:
        for i in listi[sel]:
            images.append(imageio.imread('CALIOP_BS/CALIOP_BS_{}_{}.jpg'.format(sel,i)))
    vv = str(version)
    imageio.mimsave('movie_caliop_'+vv+'.gif', images, fps=2)
    imageio.mimsave('movie_caliop_'+vv+'.mp4', images, fps=2)
    call(["gifsicle","-O","--colors","256","-o","movie_caliop_"+vv+"_opt.gif","movie_caliop_"+vv+".gif"])
