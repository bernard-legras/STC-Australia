# -*- coding: utf-8 -*-
"""
Test the hypothesis that the vorticity is near -f

Created on Sat Oct 10 10:56:03 2020

@author: Bernard Legras
"""
import pickle
import constants as cst
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5,12))
plt.subplot(3,1,1)
[trac_VO,trac_PV,trac_O3]=pickle.load(open('Vortex-track_1stVortex_ERA5PV.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_VO['lats']))
vratio = 1+trac_VO['vo']/fcorio
plt.plot(trac_VO['dates'],vratio)
plt.ylabel(u'main vortex 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
plt.subplot(3,1,2)
[trac_VO,trac_PV,trac_O3]=pickle.load(open('Vortex-track_2ndVortex_ERA5PV.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_VO['lats']))
vratio = 1+trac_VO['vo']/fcorio
plt.plot(trac_VO['dates'],vratio)
plt.ylabel(u'second vortex 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
plt.subplot(3,1,3)
[trac_VO,trac_PV,trac_O3]=pickle.load(open('Vortex-track_4thVortex_ERA5PV.pkl','rb'))
fcorio = 2 * cst.Omega * np.sin(np.deg2rad(trac_VO['lats']))
vratio = 1+trac_VO['vo']/fcorio
plt.plot(trac_VO['dates'],vratio)
plt.ylabel(u'third vortex 1+$\zeta$/f')
plt.setp(plt.xticks()[1], rotation=30, ha='right')
#fig.suptitle('Australia 2020')
plt.show()

#%%
fig = plt.figure(figsize=(3.5,3.5))
[_,trac1,_]=pickle.load(open('Vortex-track_1stVortex_ERA5PV.pkl','rb'))
pv = np.array(trac1['pv'])*1.e6
pvano = np.array(trac1['pvano'])*1.e6
dates = trac1['dates']
plt.plot(dates,pv,dates,pvano)
plt.setp(plt.xticks()[1], rotation=30, ha='right')
plt.ylabel('Potential vorticity (PVU)')
plt.legend(['PV','PV anomaly'])
plt.savefig('figs/PVCenterWithTime.png',bbox_inches='tight',dpi=300)
plt.savefig('figs/PVCenterWithTime.pdf',bbox_inches='tight',dpi=300)
plt.show()