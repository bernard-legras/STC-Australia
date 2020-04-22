#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Color map used  for the temperature anomaly

From Sergey Khaykin

Created on Fri Mar 13 19:45:47 2020

@author: Bernard Legras

"""
import numpy as np
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt

fid=open('RedWhiteBluecolorscale.txt')
RWBmap = ListedColormap(np.array([np.genfromtxt(x.rstrip('\n').split('\t')) \
                               for x in fid.readlines()])/65535,'RWB',256)

fid2=open('RedWhiteBluecolorscale13.txt')
RWBmap14 = ListedColormap(np.array([np.genfromtxt(x.rstrip('\n').split('\t')) \
                               for x in fid2.readlines()])/65535,'RWB14',14)

fid2=open('RedWhiteBluecolorscale13.txt')
color14 =np.array([np.genfromtxt(x.rstrip('\n').split('\t')) for x in fid2.readlines()])/65535

c=np.array([np.concatenate((np.linspace(0.16,1,6),np.ones(6))),
            np.concatenate((np.linspace(0.16,1,6),np.linspace(1,0,7)[1:])),
            np.concatenate((np.ones(5),np.linspace(1,0,7)))]).T

RWBmap12 = ListedColormap(c,'RWB12',12)
#RWBmap12 = ListedColormap.from_list('RWB12',['red','white','blue'],13)

a=np.outer(np.arange(0,1.005,0.01),np.ones(10))
plt.axis("off")
plt.imshow(a.T,aspect=1,cmap=RWBmap12,origin="lower")
plt.show()

# from pylab import *
# from numpy import outer
# rc('text', usetex=False)
# a=outer(arange(0,1,0.01),ones(10))
# figure(figsize=(10,5))
# subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
# maps=[m for m in cm.datad if not m.endswith("_r")]
# maps.sort()
# l=len(maps)+1
# for i, m in enumerate(maps):
#     subplot(1,l,i+1)
#     axis("off")
#     imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
#     title(m,rotation=90,fontsize=10)
# savefig("colormaps.png",dpi=100,facecolor='gray')