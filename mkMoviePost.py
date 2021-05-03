#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Generate the movies from the images produced by mkMovie.py

Created on Sun Aug 30 12:11:41 2020

@author: Bernard Legras
"""

from os.path import join
import imageio
from subprocess import call
import os
#from pygifsicle import optimize

Firsti_2plot = 0
Lasti_2plot = 180
tmp_dir = 'tmp_gif'
os.chdir(tmp_dir)

#%%
images = []
for i in range(Firsti_2plot,Lasti_2plot):
    images.append(imageio.imread(join('imO3','imO3_'+str(i)+'.jpg')))
imageio.mimsave('movie_O3.gif', images, fps=4)
imageio.mimsave('movie_O3.mp4', images, fps=4)
call(["gifsicle","-O","--colors","256","-o","movie_O3_opt.gif","movie_O3.gif"])

#%%
images = []
for i in range(Firsti_2plot,Lasti_2plot):
    print(i)
    im = imageio.imread(join('imVO','imVO_'+str(i)+'.jpg'))
    images.append(im)
imageio.mimsave('movie_VO.gif', images, fps=4)
#imageio.mimsave('movie_VO.mp4', images, fps=4)
call(["gifsicle","-O","--colors","256","-o","movie_VO_opt.gif","movie_VO.gif"])

#%%
images = []
for i in range(Firsti_2plot,Lasti_2plot):
    images.append(imageio.imread(join('imVO_sh','imVO_sh_'+str(i)+'.jpg')))
imageio.mimsave('movie_VO_sh.gif', images, fps=4)
imageio.mimsave('movie_VO_sh.mp4', images, fps=4)
call(["gifsicle","-O","--colors","256","-o","movie_VO_sh_opt.gif","movie_VO_sh.gif"])

#%%
images = []
for i in range(Firsti_2plot,Lasti_2plot):
    images.append(imageio.imread(join('imVO_vert','imVO_vert_'+str(i)+'.jpg')))
imageio.mimsave('movie_VO_vert.gif', images, fps=4)
imageio.mimsave('movie_VO_vert.mp4', images, fps=4)
call(["gifsicle","-O","--colors","256","-o","movie_VO_vert_opt.gif","movie_VO_vert.gif"])

#%%
images = []
for i in range(Firsti_2plot,Lasti_2plot):
    images.append(imageio.imread(join('imT_vert','imT_vert_'+str(i)+'.jpg')))
imageio.mimsave('movie_T_vert.gif', images, fps=4)
imageio.mimsave('movie_T_vert.mp4', images, fps=4)
call(["gifsicle","-O","--colors","256","-o","movie_T_vert_opt.gif","movie_T_vert.gif"])