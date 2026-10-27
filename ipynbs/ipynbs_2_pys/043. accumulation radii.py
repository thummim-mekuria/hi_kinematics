#!/usr/bin/env python
# coding: utf-8


import datacube as dc
import clouds as cld
import interpolate as ip
import bootleg as btlg

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import math

# %matplotlib ipympl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]

x,y = np.meshgrid(HI4PI.lons, HI4PI.lats)


w_hi = HI4PI.calculate_0th_moment()


cocat = pd.read_csv('co_catalog.csv',delim_whitespace=True,skiprows=32)


ang_R_A = (cocat['R_A']/(cocat['Dist.']*1e3) * np.cos(cocat['lat'] * u.deg.to(u.rad))) *\
            u.rad.to(u.deg) 


lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


lvs_zero = np.nansum(lvs,axis=(0,1))


whi_1 = np.nansum(HI4PI.cube,axis=1)


pd.read_csv('co_catalog.csv',skiprows=4,nrows=28)


vrms_3x = 3*cocat['v_rms']


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]
ax.imshow(w_hi.value,extent=lims)


for l,v,r,vrms in zip(cocat['lon'],cocat['vel'],ang_R_A,vrms_3x):
    ax.add_patch(Ellipse((l,v),2*r,vrms,fill=False,edgecolor='white',alpha=0.5))
    
    
ax.set_ylim(-50,50)


cocat.columns


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
ax.imshow(w_hi.value,extent=lims)


for l,b,a_bb,r_A,ind in zip(cocat['lon'],cocat['lat'],cocat['area_exact'],ang_R_A,cocat.index):
    r_bb = np.sqrt(a_bb/np.pi)
#     r_A = r #1 degree -> 8pixels bc 0.125deg/pix
    if r_bb > r_A:
        ax.add_patch(Circle((l,b),r_A,fill=False,edgecolor='white',alpha=0.5))
        ax.add_patch(Circle((l,b),r_bb,fill=False,edgecolor='red',alpha=0.5))

        print(ind)
# ax.set_xlim(180,130)


big_rbb =[1,5,31,83,88,94,98,99,110,111,123,124,125,132,150,151,152]


ras = (cocat['M'].to_numpy()*u.kg /(2*np.pi*0.025*c.M_sun/(u.pc**3)))**(1/3)


rss = ((cocat['area_exact'].to_numpy()*u.deg.to(u.rad) * (u.pc*1e3*cocat['Dist.'].to_numpy())**2))**0.5


len(ras)


np.all(ras<rss)


DHT_CO.dv,DHT_CO.db


HI4PI.dv,HI4PI.db


1/0.125 


#deal with these later --- the accumulation radius appears to be < than actual radius? could be my math tho
r_bigger_rA = cocat.loc[[cocat.index[i] for i in range(len(cocat)) if cocat['_idx'][i] in big_rbb]]


from copy import deepcopy
new_angRA = deepcopy(ang_R_A)
new_angRA.sort_values()


def accumulation_prism (v_c,b_c,l_c,R_A,vrms_A,cube):
    
    v0,v1 = v_c - 2*vrms_A, v_c + 2*vrms_A 
    #4 sigma wide, should recover > 99.7% of the velocity channels where the spectrum is
    b0,b1 = b_c - R_A, b_c + R_A #2R_A wide
    l0,l1 = l_c - R_A, l_c + R_A #2R_A wide
#     print(l0,l_c,l1)
#     print(b0,b_c,b1)
#     print(v0,v_c,v1)
#     print('\n\n')
    
    #turn physical locations into index locations
    ls = cld.ind_ang('get index', HI4PI.lons*u.deg,[l0,l1])
    bs = cld.ind_ang('get index', HI4PI.lats*u.deg,[b0,b1])
    vs = cld.ind_ang('get index', HI4PI.vels*u.km/u.s,[v0*1e3, v1*1e3])
    
    ls.sort()
    bs.sort()
    vs.sort()
    
    l0,l1 = ls
    b0,b1 = bs
    v0,v1 = vs
#     print(b0,b1)
#     #turn physical sizes into pixel lengths
#     r = math.ceil(r/(HI4PI.db.value)) #0.083deg/pixel
#     vrms = math.ceil(vrms/(HI4PI.dv.value))
    
    prism_mask = np.ones_like(cube[v0:v1,b0:b1,l0:l1],dtype=bool)
#     print(prism_mask.shape,v0,v1,b0,b1,l0,l1)
    cube[v0:v1,b0:b1,l0:l1] = cube[v0:v1,b0:b1,l0:l1] | prism_mask
    
    return cube


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
ax.imshow(w_hi.value,extent=lims)

for i in [130,69,121,118,64]:
    l,b,r = cocat['lon'][i],cocat['lat'][i],ang_R_A[i]
    ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5))
    print(l,b,r)
    
    
# ax.set_xlim(180,130)


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]
ax.imshow(w_hi.value,extent=lims)

for i in [130,69,121,118,64]:
    l,v,r,vrms =cocat['lon'][i],cocat['vel'][i],ang_R_A[i],vrms_3x[i]
    ax.add_patch(Ellipse((l,v),2*r,vrms,fill=False,edgecolor='white',alpha=0.5))
    print(l,v,r,vrms)
    
    
ax.set_ylim(-50,50)


m = np.zeros_like(HI4PI.cube,dtype=bool)
for i in range(len(cocat)):#[130,69,121,118,64]:
    l,b,v,r,vrms = cocat['lon'][i]*u.deg,cocat['lat'][i]*u.deg,\
    cocat['vel'][i]*u.km/u.s,ang_R_A[i]*u.deg,vrms_3x[i]
    #this is in pixels, so turn to physical degrees
    vrms = vrms*DHT_CO.dv
    m = accumulation_prism(v,b,l,r,vrms,m)


HI4PI.lats[585],HI4PI.lats[789]


-1.250746373060835 - 7.23 


7.23 - 15.710746373060836 


-1.250746373060835 - 7.23 


cld.ind_ang('get index', HI4PI.lats,[-1.250,15.71])


HI4PI.lats


fig,ax = plt.subplots(1,1,figsize=(15,10))
lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]

ax.imshow(np.nansum(m,axis=0),extent=lims,origin='lower')
for i in range(len(cocat)):#[130,69,121,118,64]:
    l,b,r = cocat['lon'][i],cocat['lat'][i],ang_R_A[i]
    ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5,linewidth=2))
#     print(l,b,r)


dc.dispersion_1sigma(dc.spectrum1D(HI4PI.cube[:,500:700,4050:4250]))


overlapping_parents=Dendrogram.compute(m)


fig,ax = plt.subplots(1,1,figsize=(15,10))
lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]

ax.imshow(np.nansum(m,axis=1),extent=lims,origin='lower')
ax.set_ylim(-50,50)
for i in range(len(cocat)):#[130,69,121,118,64]:
    l,v,r,vrms =cocat['lon'][i],cocat['vel'][i],ang_R_A[i],4*vrms_3x[i]*HI4PI.dv.value*1e-3
    ax.add_patch(Ellipse((l,v),2*r,vrms,fill=False,edgecolor='white',alpha=0.5,linewidth=2))
#     print(l,v,r,vrms)
    

DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


DHT_CO.lons[2773],DHT_CO.lats[236],DHT_CO.vels[31]


HI4PI.lons[2773],HI4PI.lats[236],HI4PI.lons[31]


HI4PI.lons[4161],HI4PI.lats[593],HI4PI.vels[158]


ang_R_A
