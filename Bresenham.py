#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module implements the Bresenham algorithm

Created on Sat Mar 14 16:56:35 2020

@author: Bernard Legras
"""
import numpy as np
from numba import jit, int64

@jit((int64,int64,int64,int64),nopython=True,cache=True)
def line(x0, y0, x1, y1):
    i = 0
    points = np.empty(shape=(100,2),dtype=int64)
    if (x0==x1) & (y0==y1):
        points[i,0] = x0
        points[i,1] = y0
    else:
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        dx2 = 2*dx
        dy2 = 2*dy
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1  
        x, y = x0, y0 
        if dx > dy:
            err = dx 
            while x != x1:
                points[i,0] = x
                points[i,1] = y
                i += 1
                err -= dy2
                if err < 0:
                    y += sy
                    err += dx2
                x += sx
        else:
            err = dy 
            while y != y1:
                points[i,0] = x
                points[i,1] = y
                i += 1
                err -= dx2
                if err < 0:
                    x += sx
                    err += dy2
                y += sy
        points[i,0] = x1
        points[i,1] = y1
    return points[:i+1,:]

aa = line(100,40,90,30)    