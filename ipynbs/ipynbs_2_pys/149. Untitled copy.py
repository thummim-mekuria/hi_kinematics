#!/usr/bin/env python
# coding: utf-8


fits_where = 'fitses/'


from astropy.io import fits

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np


from joblib import Parallel, delayed

import astropy.units as u


hi4pi, dht = fits.open(fits_where+'both_same_grid_720d.fits')
halpha = fits.open(fits_where+'wham-ss-DR1-v161116-170912-grid.fits')[0]
dust = fits.open(fits_where+'mean_and_std_lbd.fits')[1]


hi0 = np.log10(np.nansum(hi4pi.data,   axis=0)*1.3*1.82e18)
ha0 =          np.nansum(halpha.data,  axis=0)
h20 = np.log10(np.nansum(dht.data,     axis=0)*1.3*1.97e20)
d0  = np.log10(np.nansum(dust.data*2.8,axis=0)*1.87e21)


plt.figure(dpi=400,figsize=(20,10))

# plt.contourf(hi0,extent=[360,-360,-50,50],
#              levels=np.linspace(20,22.5,30),zorder=0)


plt.contourf(d0,extent=[180,-180,-90,90],
             levels=np.linspace(20,22.,30),norm='log',zorder=10)


# plt.contour(ha0,extent=[180,-180,-90,90],linewidths=0.2,
#              levels=np.linspace(10,150,30),norm='log',cmap='plasma',zorder=20)


plt.contour(h20,extent=[360,-360,-50,50],linewidths=0.2,
             levels=np.linspace(20,22.5,10),cmap='binary_r',zorder=30)


plt.gca().set_aspect(1)

plt.xlim(180,-180)
plt.ylim(-50,50)

plt.savefig('h2du',bbox_inches='tight')


h2_1 = np.nansum(dht.data*(np.arange(155)*1.3-100)[:,np.newaxis,np.newaxis],axis=0)*1.3*1.97e20/10**h20


plt.figure(dpi=200)
plt.contourf(np.where(h20<20,np.nan,h2_1),extent=[360,-360,-50,50],#,linewidths=0.2,
             levels=np.linspace(-30,30,30),cmap='Spectral',zorder=30)

# plt.colorbar()
plt.gca().set_aspect(1)
plt.xlim(180,-180)


dld = np.nansum(dust.data,axis=1)


plt.figure(dpi=200)
plt.contourf(np.log10(dld*2.8*1.87e21),extent=[180,-180,68,1250],
             levels=np.linspace(20,22,30))

plt.xlim(180,-180)


d_1 = np.nansum(dust.data*2.8*np.arange(68,1250)[:,np.newaxis,np.newaxis],axis=0)*1.87e21/10**d0


plt.figure(dpi=400,figsize=(20,10))

# plt.contourf(hi0,extent=[360,-360,-50,50],
#              levels=np.linspace(20,22.5,30),zorder=0)


plt.contourf(d_1,extent=[180,-180,-90,90],extend='both',
             levels=np.linspace(70,900.,200),zorder=10,cmap='binary')
# plt.colorbar()

# plt.contour(ha0,extent=[180,-180,-90,90],linewidths=0.2,
#              levels=np.linspace(10,150,30),norm='log',cmap='plasma',zorder=20)


# plt.contour(h20,extent=[360,-360,-50,50],linewidths=0.2,
#              levels=np.linspace(20,22.5,10),cmap='binary_r',zorder=20)

plt.contourf(np.where(h20<20,np.nan,h2_1),extent=[360,-360,-50,50],extend='both',antialiased=True,
             levels=np.linspace(-15,15,30),cmap='Spectral',zorder=30)


plt.contour(d_1,extent=[180,-180,-90,90],extend='both',linewidths=0.2,alpha=0.5,
             levels=np.linspace(70,900.,100),zorder=40,cmap='binary')


plt.gca().set_aspect(1)

plt.xlim(180,-180)
plt.ylim(-50,50)
plt.savefig('h2du_first',bbox_inches='tight')


plt.figure(dpi=400,figsize=(20,10))

# plt.contourf(hi0,extent=[360,-360,-50,50],
#              levels=np.linspace(20,22.5,30),zorder=0)


plt.contourf(np.where(d0<20,np.nan,d_1) ,extent=[180,-180,-90,90],extend='both',
             levels=np.linspace(70,900.,200),zorder=10,cmap='cubehelix')
# plt.colorbar()

# plt.contour(ha0,extent=[180,-180,-90,90],linewidths=0.2,
#              levels=np.linspace(10,150,30),norm='log',cmap='plasma',zorder=20)


plt.contour(h20,extent=[360,-360,-50,50],linewidths=0.2,
             levels=np.linspace(20,22.5,10),cmap='binary_r',zorder=20)

# plt.contourf(np.where(h20<20,np.nan,h2_1),extent=[360,-360,-50,50],extend='both',antialiased=True,
#              levels=np.linspace(-15,15,30),cmap='Spectral',zorder=30)


# plt.contour(d_1,extent=[180,-180,-90,90],extend='both',linewidths=0.2,alpha=0.5,
#              levels=np.linspace(70,900.,100),zorder=40,cmap='binary')


plt.gca().set_aspect(1)

plt.xlim(180,-180)
plt.ylim(-50,50)
plt.savefig('du_first',bbox_inches='tight')


hi_1 = np.nansum(hi4pi.data*(np.arange(155)*1.3-100)[:,np.newaxis,np.newaxis],axis=0)*1.3*1.82e18/10**hi0


plt.figure(dpi=400,figsize=(20,10))
plt.contourf(hi_1,extent=[360,-360,-50,50],#,linewidths=0.2,
             levels=np.linspace(-30,30,30),cmap='Spectral',extend='both',zorder=30)

plt.contour(h20,extent=[360,-360,-50,50],linewidths=0.2,
             levels=np.linspace(20,22.5,10),cmap='binary_r',zorder=30)
# plt.colorbar()
plt.gca().set_aspect(1)
plt.xlim(180,-180)
plt.savefig('hi_first',bbox_inches='tight')


plt.figure(dpi=400,figsize=(20,10))

# plt.contourf(hi0,extent=[360,-360,-50,50],
#              levels=np.linspace(20,22.5,30),zorder=0)


# plt.contourf(d_1,extent=[180,-180,-90,90],extend='both',
#              levels=np.linspace(70,900.,200),zorder=10,cmap='binary')
# # plt.colorbar()

# plt.contour(ha0,extent=[180,-180,-90,90],linewidths=0.2,
#              levels=np.linspace(10,150,30),norm='log',cmap='plasma',zorder=20)

plt.contour(h20,extent=[360,-360,-50,50],linewidths=0.2,
             levels=np.linspace(20,22.5,10),cmap='binary_r',zorder=20)

plt.contourf(np.where(h20<20,np.nan,h2_1),extent=[360,-360,-50,50],extend='both',antialiased=True,
             levels=np.linspace(-15,15,30),cmap='Spectral',zorder=30)


plt.gca().set_aspect(1)

plt.xlim(180,-180)
plt.ylim(-50,50)
plt.savefig('h2_first',bbox_inches='tight')


hilv = np.log10(np.nansum(hi4pi.data,axis=1)*0.125*1.82e18)
hibv = np.log10(np.nansum(hi4pi.data,axis=2)*0.125*1.82e18)

halv = np.nansum(halpha.data,axis=1)
habv = np.nansum(halpha.data,axis=2)

h2lv = np.log10(np.nansum(dht.data,axis=1)*0.125*1.97e20)
h2bv = np.log10(np.nansum(dht.data,axis=2)*0.125*1.97e20)


halpha.header


np.arange(162)*2.03 - 169


plt.figure(dpi=200)

plt.contourf(hilv,extent=[360,-360,-100,100],
             levels=np.linspace(20,22.5,30))

plt.contour(halv,extent=[180,-180,-159,158],linewidths=0.2,
             levels=np.linspace(10,200,30),cmap='plasma')

plt.contour(h2lv,extent=[360,-360,-100,100],linewidths=0.2,
             levels=np.linspace(20,22.5,30),cmap='binary_r')
# plt.colorbar()

plt.xlim(180,-180)

# plt.gca().set_aspect(1)


plt.figure(dpi=200)
plt.contourf(habv,levels=np.linspace(0,200,30))
# plt.colorbar()


# plt.gca().set_aspect(1)
