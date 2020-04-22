#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ancillary script
This script is used to combine the dats and datz date into a new file
To be runed only once
Version for the fourth vortex, not completed and applied!

Created on Sun Mar  8 00:11:37 2020

@author: Bernard Legras
"""
import os
import pickle,gzip

with gzip.open('OPZ-extract-pre.pkl','rb') as f:
    dats = pickle.load(f)
print(len(dats),len(datz))
i = len(dats)

#%%

with gzip.open('OPZ-extract-post.pkl','rb') as f:
    dats2 = pickle.load(f)
print(len(dats2),len(datz2))

#%%

for i2 in range(len(dats2)):
    dats[i+i2] = dats2[i2]
    datz[i+i2] = datz2[i2]
print(len(dats),len(datz))

#%%

with gzip.open('OPZ-Z.pkl','wb') as f:
    pickle.dump(datz,f)
with gzip.open('OPZ-extract.pkl','wb') as f:
    pickle.dump(dats,f)

#%% rename images
os.chdir('figs/O3maps')
for i in range(len(dats2)-2,-1,-2):
    os.rename('O3_'+str(i)+'.png','O3_'+str(i+6)+'.png')
#%%
os.chdir('../O3sect')
for i in range(len(dats2)-2,-1,-2):
    os.rename('O3sect_'+str(i)+'.png','O3sect_'+str(i+6)+'.png')
#%%
os.chdir('../VOmaps')
for i in range(len(dats2)-2,-1,-2):
    os.rename('VO_'+str(i)+'.png','VO_'+str(i+6)+'.png')
os.chdir('../VOsect')
for i in range(len(dats2)-2,-1,-2):
    os.rename('VOsect_'+str(i)+'.png','VOsect_'+str(i+6)+'.png')
os.chdir('../Tsect')
for i in range(len(dats2)-2,-1,-2):
    os.rename('Tsect_'+str(i)+'.png','Tsect_'+str(i+6)+'.png')