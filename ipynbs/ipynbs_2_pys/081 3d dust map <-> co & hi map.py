#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Ellipse

import datacube as dc

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import healpy as hlpy
import fun as fns
from astrodendro import Dendrogram
#------------------------------------------------------------
zcat = QTable(fits.open('final_cat.fits')[1].data)

hi4pi,dht = fits.open('both_same_grid_720d.fits')


#------------------------------------------------------------
hi_conv = (1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = (1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


#------------------------------------------------------------
lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50
vels = np.arange(155)*(1.3) - 100


#look at who nia&blakesley cite, also papers by lada/lombardi/alves to see where the Av=1mag contour comes from


dust_lbd = fits.open('mean_and_std_lbd.fits')
mean_dust, std_dust = dust_lbd[1],dust_lbd[2]


long_dust = np.concatenate((mean_dust.data[:,:,786:],mean_dust.data,mean_dust.data[:,:,:786]),axis=2)


# gas2dust = 1.87e21*(u.cm**-2)/u.mag #from imara & burkhart
gas2dust = 1.87e21*(u.cm**-2)/u.mag #from zhu et al: https://arxiv.org/pdf/1706.07109.pdf


dust_lons = np.arange(1572)*(-0.23) + 181
dust_lats = np.arange(786)*(0.23) -90
dust_ds = np.arange(1182) + 68


wd = np.nansum(mean_dust.data,axis=0)
whi = np.nansum(hi4pi.data*u.K*1.3*u.km/u.s*hi_conv,axis=0)
wco = np.nansum(dht.data*u.K*1.3*u.km/u.s*co_conv,axis=0)


hc_lons = np.arange(5760)*(-0.125) + 360
hc_lats = np.arange(800)*(0.125) -50
hc_vels = np.arange(155)*1.3 - 100


wd_nh = wd*u.mag*gas2dust
whi_nh = whi / (1.36*c.m_p)
wco_nh = wco / (2*1.36*c.m_p)


plt.figure(figsize=(20,9))
# plt.imshow(wd_nh.value,origin='lower',cmap=dc.cmap_min0('viridis'),norm='asinh',
#            extent=[dust_lons[0],dust_lons[-1],dust_lats[0],dust_lats[-1]])
plt.contourf(dust_lons,dust_lats,wd,levels=10,linewidths=0.7)
# plt.contour(hc_lons,hc_lats,whi_nh.value,levels=[4e20],colors='white',linewidths=0.7)
# plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[4e20],colors='black',linewidths=0.7)
plt.colorbar()
plt.xlim(180,-180)
plt.ylim(-50,50)
plt.show()


(gas2dust*u.mag*c.m_p*1.36).to(u.solMass/u.pc**2)


plt.figure(figsize=(20,13))
plt.imshow(wd_nh.value,origin='lower',cmap=dc.cmap_min0('viridis'),norm='asinh',
           extent=[dust_lons[0],dust_lons[-1],dust_lats[0],dust_lats[-1]])
plt.contour(dust_lons,dust_lats,wd,levels=[0.2,0.5,1],colors='red',linewidths=0.7)
plt.contour(hc_lons,hc_lats,whi_nh.value,levels=[4e20],colors='white',linewidths=0.7)
plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[4e20],colors='black',linewidths=0.7)

plt.xlim(180,60)
plt.ylim(-50,50)
plt.show()


plt.figure(figsize=(20,13))
plt.imshow(wd_nh.value,origin='lower',cmap=dc.cmap_min0('viridis'),norm='asinh',
           extent=[dust_lons[0],dust_lons[-1],dust_lats[0],dust_lats[-1]])
plt.contour(dust_lons,dust_lats,wd,levels=[1],colors='red',linewidths=0.7)
plt.contour(hc_lons,hc_lats,whi_nh.value,levels=[4e20],colors='white',linewidths=0.7)
plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[4e20],colors='black',linewidths=0.7)

plt.xlim(60,-60)
plt.ylim(-50,50)
plt.show()


plt.figure(figsize=(20,13))
plt.imshow(wd_nh.value,origin='lower',cmap=dc.cmap_min0('viridis'),norm='asinh',
           extent=[dust_lons[0],dust_lons[-1],dust_lats[0],dust_lats[-1]])
plt.contour(dust_lons,dust_lats,wd,levels=[1],colors='red',linewidths=0.7)
plt.contour(hc_lons,hc_lats,whi_nh.value,levels=[4e20],colors='white',linewidths=0.7)
plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[4e20],colors='black',linewidths=0.7)

plt.xlim(-60,-180)
plt.ylim(-50,50)
plt.show()


plt.figure(figsize=(20,13))
plt.imshow(wd_nh.value,origin='lower',cmap=dc.cmap_min0('viridis'),norm='asinh',
           extent=[dust_lons[0],dust_lons[-1],dust_lats[0],dust_lats[-1]])
plt.contour(dust_lons,dust_lats,wd_nh.value,levels=[4e20],colors='red',linewidths=0.7)
plt.contour(hc_lons,hc_lats,whi_nh.value,levels=[4e20],colors='white',linewidths=0.7)
plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[4e20],colors='black',linewidths=0.7)

plt.xlim(-60,-180)
plt.ylim(-50,50)
plt.show()


llvd = np.nansum(long_dust,axis=1)


plt.figure(figsize=(20,13))
plt.imshow(lvd,origin='lower',cmap=dc.cmap_min0('viridis'),norm='log',
           extent=[dust_lons[0],dust_lons[-1],dust_ds[0],dust_ds[-1]])

plt.contour(dust_lons,dust_ds,lvd,levels=[0.1,0.2,0.5],colors='white',linewidths=0.7)
# plt.contour(dust_lons,dust_lats,(lvd*u.mag*gas .value,levels=[1e21],colors='red',linewidths=0.7)

# plt.hlines([100,160,200,250,320,400,650,850,1000,1250],180,60,color='orange',linestyle='--')

# plt.xlim(180,60)
# plt.ylim(70,500)
plt.show()


d_l_pix_GMC = ((zcat['l_pix_GMC']*-0.125 + 360) - 180)/-0.23
d_b_pix_GMC = ((zcat['b_pix_GMC']*0.125 -50) +90)/0.23

d_l_pix_HI = ((zcat['l_pix_HI']*-0.125 + 360) - 180)/-0.23
d_b_pix_HI = ((zcat['b_pix_HI']*0.125 -50) +90)/0.23


for i in range(len(zcat)):
    lo,li = d_l_pix_HI[i].astype('int')
    bo,bi = d_b_pix_HI[i].astype('int')
    spc = np.nanmean(mean_dust.data[:,bo:bi,np.nanmax((0,lo)):li],axis=(1,2))
    plt.plot(spc)
    plt.vlines(zcat['distance<pc>'][i],0,np.nanmax(spc))
    plt.title(zcat['names'][i])
    plt.show()


possible_ds = np.array(([300,700],[0,200],[np.nan,np.nan],[850,980],[np.nan,np.nan],\
                       [150,220],[400,570],[np.nan,np.nan],[650,800],[700,850],[300,450],\
                       [300,450],[300,450],[120,350],[200,400],[np.nan,np.nan],[0,200]))


for i in range(len(zcat)):
    lo,li = d_l_pix_HI[i].astype('int')
    bo,bi = d_b_pix_HI[i].astype('int')
    
    plt.imshow(np.nansum(mean_dust.data[:,bo:bi,np.nanmax((0,lo)):li],axis=0),origin='lower')
    plt.title(zcat['names'][i])
    plt.show()


for i in range(17):
#     print(possible_ds[i])
    if not np.any(np.isnan(possible_ds[i])):
        do,di = possible_ds[i].astype('int')
        lo,li = d_l_pix_HI[i].astype('int')
        bo,bi = d_b_pix_HI[i].astype('int')
        
        plt.imshow(np.nansum(mean_dust.data[do:di,bo:bi,np.nanmax((0,lo)):li],axis=0),origin='lower')
        plt.title(zcat['names'][i])
        plt.show()
        

plt.figure(figsize=(20,13))
plt.imshow(wd,origin='lower',cmap=dc.cmap_min0('viridis'),norm='asinh',
           extent=[dust_lons[0],dust_lons[-1],dust_lats[0],dust_lats[-1]])

plt.contour(dust_lons,dust_lats,wd,levels=[0.1,0.2,0.5],colors=['white','red','white'],linewidths=0.7)
# plt.contour(dust_lons,dust_lats,wd_nh.value,levels=[1e21],colors='red',linewidths=0.7)

plt.xlim(180,60)
plt.ylim(-50,50)
plt.show()


lvd = np.nansum(mean_dust.data,axis=1)


plt.figure(figsize=(20,13))
plt.imshow(lvd,origin='lower',cmap=dc.cmap_min0('viridis'),norm='log',
           extent=[dust_lons[0],dust_lons[-1],dust_ds[0],dust_ds[-1]])

plt.contour(dust_lons,dust_ds,lvd,levels=[0.2,0.5],colors='white',linewidths=0.7)
# plt.contour(dust_lons,dust_lats,(lvd*u.mag*gas .value,levels=[1e21],colors='red',linewidths=0.7)

# plt.hlines([100,160,200,250,320,400,650,850,1000,1250],180,60,color='orange',linestyle='--')

plt.xlim(180,60)
plt.ylim(70,500)
plt.show()


plt.figure(figsize=(20,13))
plt.imshow(lvd,origin='lower',cmap=dc.cmap_min0('viridis'),norm='log',
           extent=[dust_lons[0],dust_lons[-1],dust_ds[0],dust_ds[-1]])

plt.contour(dust_lons,dust_ds,lvd,levels=[0.2,0.5],colors='white',linewidths=0.7)
# plt.contour(dust_lons,dust_lats,(lvd*u.mag*gas .value,levels=[1e21],colors='red',linewidths=0.7)

# plt.hlines([100,160,200,250,320,400,650,850,1000,1250],-60,-180,color='orange',linestyle='--')

plt.xlim(-60,-180)
plt.ylim(70,500)
plt.show()


dslices = np.array((np.arange(70,1250,50),np.arange(120,1300,50))).transpose()


for i,ii in enumerate(dslices):
    wlb = np.nansum(mean_dust.data[-68+ii[0]:-68+ii[1],:,:],axis=0)
#     colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    plt.figure(figsize=(17,8))

    plt.imshow(wlb,origin='lower',norm='log',extent=[180,-180,-90,90])
    
#     plt.contour(lons,lats,hilb*hi_conv,levels=[2.5,5,10,15,20],cmap='RdYlBu',zorder=0)
#     plt.contourf(lons,lats,colb*co_conv,levels=[10,1e5],colors='black')
    plt.contour(dust_lons,dust_lats,wlb,levels=[0.2],colors='white',linewidths=0.7)

    plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[1e21],colors='black',linewidths=0.7)

    
    plt.title(str(dslices[i])+' pc')
    plt.xlim(180,60)
    plt.ylim(-50,50)
    plt.show()
    

for i,ii in enumerate(dslices):
    wlb = np.nansum(mean_dust.data[-68+ii[0]:-68+ii[1],:,:],axis=0)
#     colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    plt.figure(figsize=(17,8))

    plt.imshow(wlb,origin='lower',norm='log',extent=[180,-180,-90,90])
    
#     plt.contour(lons,lats,hilb*hi_conv,levels=[2.5,5,10,15,20],cmap='RdYlBu',zorder=0)
#     plt.contourf(lons,lats,colb*co_conv,levels=[10,1e5],colors='black')
    plt.contour(dust_lons,dust_lats,wlb,levels=[0.2],colors='white',linewidths=0.7)

    plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[1e21],colors='black',linewidths=0.7)

    
    plt.title(str(dslices[i])+' pc')
    plt.xlim(60,-60)
    plt.ylim(-50,50)
    plt.show()
    

for i,ii in enumerate(dslices):
    wlb = np.nansum(mean_dust.data[-68+ii[0]:-68+ii[1],:,:],axis=0)
#     colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    plt.figure(figsize=(17,8))

    plt.imshow(wlb,origin='lower',norm='log',extent=[180,-180,-90,90])
    
#     plt.contour(lons,lats,hilb*hi_conv,levels=[2.5,5,10,15,20],cmap='RdYlBu',zorder=0)
#     plt.contourf(lons,lats,colb*co_conv,levels=[10,1e5],colors='black')
    plt.contour(dust_lons,dust_lats,wlb,levels=[0.2],colors='white',linewidths=0.7)

    plt.contour(hc_lons,hc_lats,wco_nh.value,levels=[1e21],colors='black',linewidths=0.7)

    
    plt.title(str(dslices[i])+' pc')
    plt.xlim(-60,-180)
    plt.ylim(-50,50)
    plt.show()
    

mean_dust.data.shape


1572/2


def dust_mom1(cube,rng):
    dist = np.arange(68 + rng[0],68 + rng[1])
    dist = dist[:len(cube)]

    xTx = np.array([cube[i]*dist[i] for i in range(len(dist))])
    sigxTx = np.nansum(xTx,axis=0)

    M0 = np.nansum(cube,axis=0)
    return sigxTx/M0


himom1 = fns.get_mom1(hi4pi.data,[0,155])
dmom1 = dust_mom1(mean_dust.data,[0,1182])


#because i have l,b,d information, would be cool to reproject this ^ into a top down view
#its just the solar neighborhood (d<1.25kpc) so maddalena etc might have to go? if i end up using it
#compare ->dust<->h2<->hi<- bc at some points the dust contours went with hi & at some points with h2
#also go back and reread kmt2008 etc to figure out where in the pdr the dust goes -> maybe the answer
#is places that have high hi + dust shielding get to form h2 & those with just one or the other don't?
#also in terms of boundaries, av=0.5mag seems to trace the co boundaries better than av=1? cross check with lewis,
#nia etc 
# also look at relationships between peak velocity and position if it might be indicative of winds eg california


zcat
