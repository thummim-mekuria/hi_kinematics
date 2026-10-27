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

#------------------------------------------------------------
zcat = QTable(fits.open('final_cat.fits')[1].data)

hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1]


#------------------------------------------------------------
hi_conv = (1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = (1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


#------------------------------------------------------------
lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50
vels = np.arange(155)*(1.3) - 100


lvd = np.nansum(dust.data,axis=1)


bvd = np.nansum(dust.data,axis=2)


lvh = np.nansum(hi4pi.data,axis=1)


bvh = np.nansum(hi4pi.data,axis=2)
bvc = np.nansum(dht.data,axis=2)


lvc = np.nansum(dht.data,axis=1)


lbd = np.nansum(dust.data,axis=0)
lbh = np.nansum(hi4pi.data,axis=0)
lbc = np.nansum(dht.data,axis=0)


dlons = np.arange(1572)*-0.23 + 180
dlats = np.arange(786)*0.23 - 90
ddist = np.arange(1182) + 68

hclons = np.arange(5760)*-0.125 + 360
hclats = np.arange(800)*0.125 - 50
hcvels = np.arange(155)*1.3 -100


plt.figure(figsize=(10,10))
plt.imshow(lvd,norm='log',origin='lower',extent=[dlons[0],dlons[-1],ddist[0],ddist[-1]])
plt.contour(dlons,ddist,lvd,
            levels=[0.15,0.2,0.3,0.4,0.5],cmap='binary_r',linewidths=0.6)
plt.scatter(np.nanmean(zcat['lon_GMC<deg>'],axis=1),zcat['distance<pc>'],marker='x',color='red')
plt.ylim(68,500)
# plt.contour(np.arange(250),np.arange(300),mom1,levels=[0,250,500,750,870],cmap='binary_r',linewidths=0.8)


plt.figure(figsize=(10,10))
plt.imshow(bvd,norm='log',origin='lower',extent=[dlats[0],dlats[-1],ddist[0],ddist[-1]])
plt.contour(dlats,ddist,bvd,
            levels=[0.5,1,1.5],cmap='binary_r',linewidths=0.6)
plt.scatter(np.nanmean(zcat['lat_GMC<deg>'],axis=1),zcat['distance<pc>'],marker='x',color='red')
plt.ylim(68,500)
plt.xlim(-50,50)
# plt.contour(np.arange(250),np.arange(300),mom1,levels=[0,250,500,750,870],cmap='binary_r',linewidths=0.8)


plt.figure(figsize=(20,10))
plt.imshow(lbd,norm='log',origin='lower',extent=[dlons[0],dlons[-1],dlats[0],dlats[-1]])
plt.contour(dlons,dlats,lbd,
            levels=[0.1,0.5,1,1.5],cmap='binary_r',linewidths=0.6)
plt.scatter(np.nanmean(zcat['lon_GMC<deg>'],axis=1),np.nanmean(zcat['lat_GMC<deg>'],axis=1),
            marker='x',color='red')
# plt.ylim(68,500)
plt.ylim(-50,50)
# plt.contour(np.arange(250),np.arange(300),mom1,levels=[0,250,500,750,870],cmap='binary_r',linewidths=0.8)


lbc_solMpc = (lbc*1.3*co_conv*u.K*u.km/u.s).to(u.solMass/u.pc**2)


lbh_solMpc = (lbh*1.3*hi_conv*u.K*u.km/u.s).to(u.solMass/u.pc**2)


plt.figure(figsize=(20,10))
plt.imshow(lbd,norm='log',origin='lower',extent=[dlons[0],dlons[-1],dlats[0],dlats[-1]])
plt.contour(dlons,dlats,lbd,
            levels=[0.1,0.3,0.5],colors='white',linewidths=0.6)
plt.contour(hclons,hclats,lbh_solMpc.value,
            levels=[5,10,15],colors='orange',linewidths=0.6)
plt.contourf(hclons,hclats,lbc_solMpc.value,
            levels=[1e-3,1e5],colors='black',linewidths=0.6)

plt.xlim(180,-180)
plt.ylim(-60,60)


hclons[1440],hclons[4320]


dust_grid = np.meshgrid(dlons,dlats[173:610])


hc_grid = np.meshgrid(hclons[1440:4320],hclats)


dust_coords = np.array(dust_grid).transpose().reshape(-1,2)


hc_coords = np.array(hc_grid).transpose().reshape(-1,2)


hcs = pd.DataFrame(hc_coords)


hcs['hsig'] = np.ravel(lbh_solMpc.value[:,1440:4320].transpose())
hcs['csig'] = np.ravel(lbc_solMpc.value[:,1440:4320].transpose())


hcs.columns = ['lon', 'lat', 'hsig', 'csig']


len(hcs)


bns = np.log10(np.arange(1e-3,10000))


dst = pd.DataFrame(dust_coords)


dst['dav'] = np.ravel(lbd[173:610,:].transpose())


dst.columns = ['lon', 'lat', 'dav']


plt.scatter(dst['lon'][::387],dst['lat'][::387],
            c=dst['dav'][::387],norm='asinh',alpha=0.5)


plt.hist(np.log10(hcs['hsig']),bins=bns,alpha=0.7)
plt.hist(np.log10(hcs['csig']),bins=bns,alpha=0.7)
plt.hist(dst['dav'],alpha=0.7)
plt.yscale('log')
plt.xlim(0,)
plt.show()


plt.scatter(hcs['hsig'],hcs['csig'],s=1e-2,alpha=0.2)
plt.loglog()
plt.xlabel('h')
plt.ylabel('c')
plt.ylim(1e-1,)


plt.scatter(hcs['lon'][::590],hcs['lat'][::590],
            c=hcs['hsig'][::590],norm='asinh',alpha=0.5)


plt.scatter(hcs['lon'][::590],hcs['lat'][::590],
            c=hcs['csig'][::590],norm='asinh',alpha=0.5)


dust_coords_all = np.meshgrid(ddist,dlats,dlons)


dust_all = np.array(dust_coords_all).transpose().reshape(-1,3)


plt.figure(figsize=(20,10))
plt.imshow(lbd,norm='log',origin='lower',extent=[dlons[0],dlons[-1],dlats[0],dlats[-1]])
plt.contour(dlons,dlats,lbd,
            levels=[0.5],colors='white',linewidths=0.4)
plt.contour(hclons,hclats,lbh,
            levels=[400],colors='orange',linewidths=0.4)
plt.contour(hclons,hclats,lbc,
            levels=[1.5],colors='black',linewidths=0.4)

plt.xlim(180,-180)
plt.ylim(-60,60)
