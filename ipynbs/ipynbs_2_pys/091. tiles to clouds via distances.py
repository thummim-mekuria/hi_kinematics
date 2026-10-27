#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import datacube as dc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
import fun as fns

import astropy.units as u
import astropy.constants as c

from matplotlib.patches import Circle,Rectangle


hi4pi,dht = fits.open('both_same_grid_720d.fits')
zcat = QTable(fits.open('final_cat.fits')[1].data)


lons = np.arange(5760)*-0.125 + 360
lats = np.arange(800)*0.125 - 50
vels = np.arange(155)*1.3 - 100


hi_conv = (1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = (1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


whi = (hi_conv*u.K*1.3*u.km/u.s*np.nansum(hi4pi.data,axis=0)).to(u.solMass/u.pc**2)
wco = (co_conv*u.K*1.3*u.km/u.s*np.nansum(dht.data,axis=0)).to(u.solMass/u.pc**2)


whole_gal_ds = pd.read_csv('whole_gal_dists.csv')


million_ds = pd.read_csv('whole_gal_ds_zrzr.csv')


h2s = fits.open('h2_clouds.fits')
hi_gfs = fits.open('hi_clouds.fits')
hi_gls = fits.open('gmcless_tiles.fits')


tile_coords = QTable(fits.open('gmcless_cat.fits')[1].data)


tile_coords


[(plt.imshow(np.nanmean(hi_gls[i].data,0),origin='lower',cmap=dc.cmap_min0('viridis')),
  plt.show() )for i in range(7)]


tile_1_ds = whole_gal_ds[(whole_gal_ds['lon<deg>'] < -45) & \
             (whole_gal_ds['lon<deg>'] > -132) & \
             (whole_gal_ds['lat<deg>'] > 6) &
             (whole_gal_ds['lat<deg>'] < 50)]


plt.contourf(np.arange(315-360,228.75-360,-0.125),np.arange(-100,101,1.3),
            np.nanmean(hi_gls[0].data,1),
           levels=10,cmap='binary')

plt.scatter(tile_1_ds['lon<deg>'][::10],tile_1_ds['vel<km/s>'][::10],c=tile_1_ds['d2<kpc>'][::10],
           marker='x',s=tile_1_ds['d2_prob'][::10]*10,cmap='Paired',vmax=5)
plt.colorbar()

plt.ylim(-60,50)


plt.contourf(np.arange(315-360,228.75-360,-0.125),np.arange(-100,101,1.3),
            np.nanmean(hi_gls[0].data,1),
           levels=10,cmap='binary')

plt.scatter(tile_1_ds['lon<deg>'][::10],tile_1_ds['vel<km/s>'][::10],c=tile_1_ds['d1<kpc>'][::10]*1000,
           marker='x',s=tile_1_ds['d1_prob'][::10]**2,cmap='Paired',norm='log')
plt.colorbar()

plt.ylim(-60,50)


plt.contourf(np.arange(315-360,228.75-360,-0.125),np.arange(6.25,50,0.125),
            np.nanmean(hi_gls[0].data,0),
           levels=5,cmap='binary')

plt.scatter(tile_1_ds['lon<deg>'],tile_1_ds['lat<deg>'],c=tile_1_ds['d1<kpc>'],
           marker='x',s=10,cmap='Paired',vmax=8)
plt.colorbar()
plt.xlim(-47,-132)
# plt.ylim(5,30)


tile_2_ds = whole_gal_ds[(whole_gal_ds['lon<deg>'] < -45) & \
             (whole_gal_ds['lon<deg>'] > -132) & \
             (whole_gal_ds['lat<deg>'] > -50) &
             (whole_gal_ds['lat<deg>'] < -18.75)]


plt.contourf(np.arange(315-360,228.75-360,-0.125),np.arange(-100,101,1.3),
            np.nanmean(hi_gls[1].data,1),
           levels=10,cmap='binary')

plt.scatter(tile_2_ds['lon<deg>'][::10],tile_2_ds['vel<km/s>'][::10],
            c=tile_2_ds['d1<kpc>'][::10],
           marker='x',s=tile_2_ds['d1_prob'][::10]*10,cmap='Paired',vmax=5)
plt.colorbar()

plt.ylim(-60,50)


plt.contourf(np.arange(315-360,228.75-360,-0.125),np.arange(-50,-18.75,0.125),
            np.nanmean(hi_gls[1].data,0),
           levels=5,cmap='binary')

plt.scatter(tile_2_ds['lon<deg>'],tile_2_ds['lat<deg>'],c=tile_2_ds['d1<kpc>'],
           marker='x',s=10,cmap='Paired',vmax=8)
plt.colorbar()
plt.xlim(-47,-132)


#this one is completely beyond the latitude range of the co emission, 
#so it might not make sense to include it
tile_3_ds = whole_gal_ds[(whole_gal_ds['lon<deg>'] > 35) & \
             (whole_gal_ds['lon<deg>'] < -131.25) &
             (whole_gal_ds['lat<deg>'] > -50) &
             (whole_gal_ds['lat<deg>'] < -31.25)]


tile_4_ds = whole_gal_ds[(whole_gal_ds['lon<deg>'] > 35) & \
             (whole_gal_ds['lon<deg>'] < 85) &
             (whole_gal_ds['lat<deg>'] > -31.25) &
             (whole_gal_ds['lat<deg>'] < -12.5)]


plt.contourf(np.arange(85,35,-0.125),
            np.arange(-31.25,-12.5,0.125),
            np.nanmean(hi_gls[3].data,0),
           levels=5,cmap='binary')

plt.scatter(tile_4_ds['lon<deg>'],tile_4_ds['lat<deg>'],c=tile_4_ds['d1<kpc>'],
           marker='x',s=10,cmap='Paired',vmax=8)
plt.colorbar()


plt.contourf(np.arange(228.75-360,35-360,-0.125),np.arange(-100,101,1.3),
            np.nanmean(hi_gls[2].data,0),
           levels=10,cmap='binary')

plt.scatter(tile_3_ds['lon<deg>'][::10],tile_3_ds['vel<km/s>'][::10],
            c=tile_3_ds['d1<kpc>'][::10],
           marker='x',s=tile_3_ds['d1_prob'][::10]*10,cmap='Paired',vmax=5)
plt.colorbar()

plt.ylim(-60,50)


np.arange(228.75-360,35-360,-0.125)[::10]


(np.arange(5760)*-0.125 + 360)[1050:2600][::10]


228.75-360


plt.figure(figsize=(20,5))
plt.scatter(whole_gal_ds['lon<deg>'],whole_gal_ds['vel<km/s>'],c=whole_gal_ds['d2<kpc>']*1000,
           marker='.',s=20,cmap=dc.cmap_min0('Paired'),vmin=10,vmax=2000,alpha=whole_gal_ds['d2_prob'])
plt.colorbar()
# [plt.annotate(zcat['names'][i]+'/'+str(zcat['distance<pc>'][i])+' pc',
#              (np.nanmean(zcat['lon_GMC<deg>'][i]),
#               np.nanmean(zcat['vel_GMC<km/s>'][i])),fontsize='medium')
#         for i in range(17)]

plt.xlim(180,-180)
# plt.ylim(-10,40)
plt.show()


plt.figure(figsize=(20,5))
plt.scatter(whole_gal_ds['lon<deg>'],whole_gal_ds['vel<km/s>'],c=whole_gal_ds['d1<kpc>']*1000,
           marker='.',s=20,cmap='Paired',vmax=3000,alpha=whole_gal_ds['d1_prob'])
plt.colorbar()
# [plt.annotate(zcat['names'][i]+'/'+str(zcat['distance<pc>'][i])+' pc',
#              (np.nanmean(zcat['lon_GMC<deg>'][i]),
#               np.nanmean(zcat['vel_GMC<km/s>'][i])),fontsize='medium')
#         for i in range(17)]

plt.xlim(180,-180)
# plt.ylim(-10,40)
plt.show()


plt.figure(figsize=(15,5))
plt.scatter(whole_gal_ds['lon<deg>'],whole_gal_ds['vel<km/s>'],c=whole_gal_ds['d1<kpc>']*1000,
           marker='.',s=30,cmap='Paired',vmin=100,vmax=700)
plt.colorbar()
[plt.annotate(zcat['names'][i]+'/'+str(zcat['distance<pc>'][i])+' pc',
             (np.nanmean(zcat['lon_GMC<deg>'][i]),
              np.nanmean(zcat['vel_GMC<km/s>'][i])),fontsize='medium')
        for i in range(17)]

plt.xlim(130,30)
plt.ylim(-10,40)
plt.show()


plt.figure(figsize=(15,5))
plt.scatter(whole_gal_ds['lon<deg>'],whole_gal_ds['vel<km/s>'],c=whole_gal_ds['d1<kpc>']*1000,
           marker='.',s=30,cmap='Paired',norm='log',vmin=100,vmax=400)
plt.colorbar()
[plt.annotate(zcat['names'][i]+'/'+str(zcat['distance<pc>'][i])+' pc',
             (np.nanmean(zcat['lon_GMC<deg>'][i]),
              np.nanmean(zcat['vel_GMC<km/s>'][i])),fontsize='medium')
        for i in range(17)]

plt.xlim(-50,-80)
plt.ylim(-10,40)
plt.show()


plt.figure(figsize=(15,5))
plt.scatter(whole_gal_ds['lon<deg>'],whole_gal_ds['vel<km/s>'],c=whole_gal_ds['d1<kpc>']*1000,
           marker='.',s=30,cmap='Paired',norm='asinh')
plt.colorbar()
[plt.annotate(zcat['names'][i]+'/'+str(zcat['distance<pc>'][i])+' pc',
             (np.nanmean(zcat['lon_GMC<deg>'][i]),
              np.nanmean(zcat['vel_GMC<km/s>'][i])),fontsize='medium')
        for i in range(17)]

plt.xlim(-130,-180)
plt.ylim(0,50)


tile1 = hi4pi.data[:,450:,360:1050]


tile2 = hi4pi.data[:,0:250,360:1050]


tile3 = hi4pi.data[:,0:150,1050:2600]


tile4 = hi4pi.data[:,150:300,2200:2600]


tile5 = hi4pi.data[:,150:300,1700:2000]


tile6 = hi4pi.data[:,500:,2150:2600]


tile7 = hi4pi.data[:,500:,1350:1550]


gmcless_bs = 0.125*np.array(([450,800],[0,250],[0,150],[150,300],
                             [150,300],[500,800],[500,800])) - 50


gmcless_ls = -0.125*np.array(([360,1050],[360,1050],[1050,2600],[2200,2600],
                             [1700,2000],[2150,2600],[1350,1550])) + 360


names = ['tile_'+np.arange(1,8).astype('str')[i] for i in range(7)]


gless_cat = QTable((names,gmcless_ls,gmcless_bs),names=['name','lon<deg>','lat<deg>'])


gless_cat.write('gmcless_cat.fits',overwrite=True)


gmclesses = [tile1,tile2,tile3,tile4,tile5,tile6,tile7]


fits.HDUList([fits.PrimaryHDU(gmclesses[0])]
             +[fits.ImageHDU(gmclesses[i]) 
  for i in range(1,7)]).writeto('gmcless_tiles.fits',overwrite=True)


tile_1_ds = whole_gal_ds[-45 >= whole_gal_ds['lon<deg>']]\
[whole_gal_ds['lon<deg>'] >= -131.25][6.25<=whole_gal_ds['lat<deg>']][whole_gal_ds['lat<deg>']<=50.]


plt.scatter(tile_1_ds['lon<deg>']+360,tile_1_ds['lat<deg>'],c=tile_1_ds['d1<kpc>'],norm='asinh')
plt.colorbar()


[(plt.contourf(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
             np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
             np.nansum(gmclesses[i],axis=0),origin='lower',
            levels=np.arange(0,800,100)),
plt.colorbar(),
plt.contour(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
             np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
             np.nansum(gmclesses[i],axis=0),origin='lower',
            levels=[400],colors='white',linewidths=0.6),
 plt.show())
 for i in range(1)]


[(plt.contourf(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
             np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
             fns.get_mom1(gmclesses[i],[0,155]),origin='lower',
            levels=np.arange(-50,50,10),norm='asinh'),
plt.colorbar(),
plt.contour(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
             np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
             np.nansum(gmclesses[i],axis=0),origin='lower',
            levels=[400],colors='white',linewidths=0.6),
 plt.show())
 for i in range(7)]


[(plt.contourf(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
             np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
             fns.get_momN(gmclesses[i],[0,155],2),origin='lower',
            levels=np.arange(0,40,8)),
plt.colorbar(),
# plt.contour(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
#              np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
#              np.nansum(gmclesses[i],axis=0),origin='lower',
#             levels=[400],colors='white',linewidths=0.6),
  
plt.contour(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
         np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
         fns.get_mom1(gmclesses[i],[0,155]),origin='lower',
        levels=[-10,0,10],colors='white',linewidths=0.5),
 plt.show())
 for i in range(7)]


plt.contourf(np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.arange(-100,101,1.3),
             np.nansum(tile1,axis=2),origin='lower',
            levels=np.arange(0,18000,1000))

# plt.contour(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
#              np.arange(-100,101,1.3),
#              np.nansum(tile1,axis=1),origin='lower',
#             levels=[400],colors='white',linewidths=0.6)


[(plt.contourf(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
             np.arange(-100,101,1.3),
             np.nansum(gmclesses[i],axis=1),origin='lower',
            levels=np.arange(0,9000,700)),
plt.colorbar(),
# plt.contour(np.arange(gmcless_ls[i][0],gmcless_ls[i][1],-0.125),
#              np.arange(gmcless_bs[i][0],gmcless_bs[i][1],0.125),
#              np.nansum(gmclesses[i],axis=0),origin='lower',
#             levels=[400],colors='white',linewidths=0.6),
 plt.ylim(-20,20),
 plt.show()
 )for i in range(7)]


# plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
#              np.arange(-100,101,1.3),
#              np.nansum(tile1,axis=1),origin='lower',
#             levels=np.arange(0,9000,400))

# # plt.contour(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
# #              np.arange(-100,101,1.3),
# #              np.nansum(tile1,axis=1),origin='lower',
#             levels=[400],colors='white',linewidths=0.6)


plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             fns.get_mom1(tile1,[0,155]),origin='lower',
            levels=np.arange(-50,50,10),norm='asinh')
plt.colorbar()
plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1,axis=0),origin='lower',
            levels=[200,400,600,800],cmap='binary_r',alpha=0.3)


vels[40],vels[120]


plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1[77:,:,:],axis=0),origin='lower',
            levels=np.arange(0,800,100))

plt.contour(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1[77:,:,:],axis=0),origin='lower',
            levels=[200,400],colors='white',linewidths=0.6)


plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1[35:110,:,:],axis=0),origin='lower',
            levels=np.arange(0,800,100))

plt.contour(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1[35:110,:,:],axis=0),origin='lower',
            levels=[200,400],colors='white',linewidths=0.6)


plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1[:77,:,:],axis=0),origin='lower',
            levels=np.arange(0,800,100))

plt.contour(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1[:77,:,:],axis=0),origin='lower',
            levels=[200,400],colors='white',linewidths=0.6)


plt.contourf(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1,axis=0),origin='lower',
            levels=np.arange(0,800,100))

plt.contour(np.arange(gmcless_ls[0][0],gmcless_ls[0][1],-0.125),
             np.arange(gmcless_bs[0][0],gmcless_bs[0][1],0.125),
             np.nansum(tile1,axis=0),origin='lower',
            levels=[400],colors='white',linewidths=0.6)


plt.scatter(tile_1_ds['lon<deg>'],tile_1_ds['vel<km/s>'],c=tile_1_ds['d1<kpc>'])


[(plt.plot(np.nanmean(cube,axis=(1,2))),
   plt.show()) for cube in gmclesses]


[(plt.plot(np.nanmean(cube,axis=(1,2))),
   plt.show()) for cube in [tile1,tile2,tile3,tile4,tile5,tile6,tile7]]


plt.contourf(np.arange(200),np.arange(300),np.nansum(tile7,axis=0),origin='lower',
            levels=[0,200,400,600,800])


plt.contourf(np.arange(450),np.arange(300),np.nansum(tile6,axis=0),origin='lower',
            levels=[0,200,400,600,800])


plt.contourf(np.arange(300),np.arange(150),np.nansum(tile5,axis=0),origin='lower',
            levels=[0,200,400,600,800])


plt.contourf(np.arange(400),np.arange(150),np.nansum(tile4,axis=0),origin='lower',
            levels=[0,200,400,600,800])


plt.contourf(np.arange(1550),np.arange(150),np.nansum(tile3,axis=0),origin='lower',
            levels=[0,200,400,600,800])


plt.contourf(np.arange(690),np.arange(250),np.nansum(tile2,axis=0),origin='lower',
            levels=[0,200,400,600,800])


plt.contourf(np.arange(690),np.arange(350),np.nansum(tile1,axis=0),origin='lower',
            levels=[0,200,400,600,800])


whole_gal_ds


0.083*u.deg.to(u.arcmin)


0.125*u.deg.to(u.arcmin)


0.23*u.deg.to(u.arcmin)


(50*u.kpc*1*u.arcmin.to(u.rad)).to(u.pc)


(400*u.pc*7.5*u.arcmin.to(u.rad)).to(u.pc)


#possible new directions
#could look at whats going in the gmcfull hi as a function of evolution(traced by sfr)
#lmc

#focus on properties that don't rely on distance - surface density, turbulence, n-pdfs of column density
#comparing the velocity field of the hi<->gmc
#
