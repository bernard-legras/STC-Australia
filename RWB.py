#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Color map used  for the temperature anomaly

From Sergey Khaykin

Created on Fri Mar 13 19:45:47 2020

@author: Bernard Legras

"""
import numpy as np
from matplotlib.colors import ListedColormap
fid=open('RedWhiteBluecolorscale.txt')
RWBmap = ListedColormap(np.array([np.genfromtxt(x.rstrip('\n').split('\t')) \
                               for x in fid.readlines()])/65535,'RWB')
