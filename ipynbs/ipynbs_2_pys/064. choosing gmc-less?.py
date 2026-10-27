#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import datacube as dc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle


ds = pd.read_csv('whole_gal_dists.csv')


zcat = QTable(fits.open('final_cat.fits')[1].data)


hi4pi = fits.open('HI4PI_flatcalib_80.fits')[0]


hi4pi_rawr = fits.open('HI4PI.fits')[0]
wrar =  np.nansum(hi4pi_rawr.data,axis=0)*u.K \
* 1.288 * u.km/u.s\
* 1.36 * c.m_p * 1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) 
dense_enf_msk = (wrar > 10*c.M_sun/u.pc**2)


dense_enf_5_msk = (wrar > 5*c.M_sun/u.pc**2)


whi = np.nansum(hi4pi.data,axis=0)


dht = fits.open('DHT_CO_longer.fits')[0]


wco = np.nansum(dht.data,axis=0)


colons = np.arange(2881)*(-0.125) + 180.125
colats = np.arange(481)*(0.125) - 30.125


hilons = np.arange(4323)*(-0.0833) + 180.166
hilats = np.arange(1200)*(0.0833) - 49.999


from matplotlib.patches import Circle


whi_r = np.nansum(hi4pi_rawr.data,axis=0)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,origin='lower',\
          cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=5e3)

ax.contourf(colons,colats,wco,levels=[0.5,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))

ax.contour(hilons,hilats,whi,levels=[1400,2101,2300],cmap='binary_r')
# ax.contour(hilons,hilats,whi_r,levels=[200,250,300],cmap='binary_r')

plt.ylim(-30,30)
plt.xlim(180,0)
# plt.ylim(-30,-10)
# plt.xlim(120,80)

# plt.ylim(0,30)
# plt.xlim(80,20)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=5e3)
ax.contourf(colons,colats,wco,levels=[0.5,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.contour(hilons,hilats,whi,levels=[1800,2000,2200],cmap='binary_r')


plt.ylim(-30,30)
plt.xlim(0,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmin=1800,vmax=5e3)
ax.contourf(colons,colats,wco,levels=[0.5,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
# ax.contour(hilons,hilats,whi,levels=[1800,2000,2200],cmap='binary_r')


plt.ylim(-30,30)
plt.xlim(180,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmin=1800,vmax=5e3)
ax.contourf(colons,colats,wco,levels=[0.5,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
# ax.contour(hilons,hilats,whi,levels=[1800,2000,2200],cmap='binary_r')


plt.ylim(-30,30)
plt.xlim(180,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=5e3)
ax.contour(hilons,hilats,whi,np.arange(2000,5000,500),cmap='rainbow')

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    

plt.ylim(-30,30)
plt.xlim(180,0)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=5e3)
ax.contour(hilons,hilats,whi,np.arange(2000,5000,500),cmap='rainbow')

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    

plt.ylim(-30,30)
plt.xlim(0,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=5e3)
ax.contour(hilons,hilats,whi,np.arange(2000,4000,500),cmap='rainbow')

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))

plt.ylim(-30,30)
plt.xlim(180,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=2e3)
# ax.contour(hilons,hilats,whi,np.arange(2000,4000,500),cmap='rainbow')

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))

plt.ylim(-30,30)
plt.xlim(180,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=2.5e3)
ax.contour(hilons,hilats,whi,np.arange(10,1900,200),cmap='binary_r',zorder=1)

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))

plt.ylim(-30,30)
plt.xlim(180,50)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi_r*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=400)
ax.contourf(hilons,hilats,whi,np.arange(100,1900,200),cmap='binary_r',zorder=1,alpha=0.5)

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))

plt.ylim(-30,30)
plt.xlim(180,50)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=2.5e3)
ax.contour(hilons,hilats,whi,np.arange(10,1900,200),cmap='binary_r',zorder=1)

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))

plt.ylim(-30,30)
plt.xlim(-50,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi_r*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=400)
ax.contourf(hilons,hilats,whi,np.arange(100,1900,200),cmap='binary_r',zorder=1,alpha=0.5)

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))

plt.ylim(-30,30)
plt.xlim(-50,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi_r*dense_enf_5_msk,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=400)
ax.contourf(hilons,hilats,whi,np.arange(100,1900,200),cmap='binary_r',zorder=1)

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))
ax.add_patch(Rectangle((180,-10),-360,20,color='white'))


plt.ylim(-30,30)
plt.xlim(180,-180)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.imshow(whi_r,\
           origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=500,vmin=100)
# ax.contourf(hilons,hilats,whi,np.arange(100,1900,200),cmap='binary_r',zorder=1)

ax.contourf(colons,colats,wco,levels=[1,1e3],colors=['white'])

for i in range(21):
    lc, bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    ax.add_patch(Circle((lc,bc),rhi,color='white'))
    
ax.add_patch(Rectangle((50,-50),-100,100,color='white'))
ax.add_patch(Rectangle((180,-10),-360,20,color='white'))


plt.ylim(-30,30)
plt.xlim(180,-180)


hi4pi_rawr.data.shape


hilons = np.arange(4323)*(-0.0833) + 180.125


hilats = np.arange(1200)*0.0833 - 49.999


1563:2762 #lons


961: and :241 and 481 720 #lats are


masked_hi4pi = np.zeros_like(hi4pi_rawr.data)


masked_hi4pi[:,240:481,:1563] = hi4pi_rawr.data[:,240:481,:1563]
masked_hi4pi[:,720:961,:1563] = hi4pi_rawr.data[:,720:961,:1563]
masked_hi4pi[:,240:481,2762:] = hi4pi_rawr.data[:,240:481,2762:]
masked_hi4pi[:,720:961,2762:] = hi4pi_rawr.data[:,720:961,2762:]


wmsk = np.nansum(masked_hi4pi,axis=0)


wmsklv = np.nansum(masked_hi4pi,axis=1)


plt.imshow(wmsk,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=500,vmin=100)


plt.imshow(wlvr,origin='lower',cmap=dc.cmap_min0('rainbow'),extent=[180,-180,-50,50])#,vmax=500,vmin=100)
plt.imshow(wmsklv,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50])#,vmax=500,vmin=100)


wlvr = np.nansum(hi4pi_rawr.data,axis=1)


plt.imshow(wlvr,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50])#,vmax=500,vmin=100)


wbvr = np.nansum(hi4pi_rawr.data,axis=2)
wmsklv = np.nansum(masked_hi4pi,axis=2)


plt.imshow(wmsklv,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-50,50,-50,50])#,vmax=500,vmin=100)
plt.imshow(wbvr,origin='lower',cmap=dc.cmap_min0('rainbow'),extent=[-50,50,-50,50])#,vmax=500,vmin=100)
# plt.imshow(wmsklv,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50])#,vmax=500,vmin=100)


(100*u.K \
* 1.288 * u.km/u.s\
* 1.36 * c.m_p * 1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s))).to(c.M_sun/u.pc**2)
