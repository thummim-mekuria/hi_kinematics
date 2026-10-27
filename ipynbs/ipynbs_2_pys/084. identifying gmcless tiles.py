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


whi = hi_conv*u.K*1.3*u.km/u.s*np.nansum(hi4pi.data,axis=0)
wco = co_conv*u.K*1.3*u.km/u.s*np.nansum(dht.data,axis=0)


whi_solm = whi.to(u.solMass/u.pc**2).value
wco_solm = wco.to(u.solMass/u.pc**2).value


whole_gal_ds = pd.read_csv('whole_gal_dists.csv')


h2s = fits.open('h2_clouds.fits')
his = fits.open('hi_clouds.fits')


dht_mom2=fns.get_momN(dht.data,[0,155],2)


hi4pi_mom2=fns.get_momN(hi4pi.data,[0,155],2)


fig,ax = plt.subplots(1,1,figsize=(15,12))
ax.imshow(whi_solm,
           origin='lower',norm='log')

#     plt.colorbar(orientation='horizontal')

ax.contour(np.arange(5760),np.arange(800),whi_solm,
            levels=[5,10,15],cmap='RdYlBu',zorder=0)
ax.contourf(np.arange(5760),np.arange(800),wco_solm,
             levels=[1,1e5],colors='black')


ax.add_patch(Rectangle((0,360),5760,80,color='black'))
ax.add_patch(Rectangle((0,0),360,800,color='black'))
ax.add_patch(Rectangle((5759,0),-360,800,color='black'))
ax.add_patch(Rectangle((3240,0),-720,800,color='black'))


ax.add_patch(Rectangle((360,450),1050-360,800-450,color='black',fill=False))
ax.add_patch(Rectangle((360,0),1050-360,250,color='black',fill=False))
ax.add_patch(Rectangle((1050,0),2600-1050,150,color='black',fill=False))
ax.add_patch(Rectangle((2200,150),500,150,color='black',fill=False))
ax.add_patch(Rectangle((1700,150),300,150,color='black',fill=False))
ax.add_patch(Rectangle((2150,500),500,300,color='black',fill=False))
ax.add_patch(Rectangle((1350,500),200,300,color='black',fill=False))


plt.xlim(0,2880)
plt.ylim(0,800)


fig,ax = plt.subplots(1,1,figsize=(15,12))
ax.imshow(whi_solm,
           origin='lower',norm='log',extent=[360,-360,-50,50])

#     plt.colorbar(orientation='horizontal')

ax.contour(lons,lats,whi_solm,
            levels=[5,10,15],cmap='RdYlBu',zorder=0)
ax.contourf(lons,lats,wco_solm,
             levels=[1e-5,1e5],colors='black')


for i in range(17):
    ax.add_patch(Circle(
        (np.nanmean(zcat['lon_GMC<deg>'][i]),np.nanmean(zcat['lat_GMC<deg>'][i])),
        np.nanmean(zcat['radius_acc_HI<deg>'][i]),color='black'))
    ax.add_patch(Circle(
        (360+np.nanmean(zcat['lon_GMC<deg>'][i]),np.nanmean(zcat['lat_GMC<deg>'][i])),
        np.nanmean(zcat['radius_acc_HI<deg>'][i]),color='black'))
    ax.add_patch(Circle(
        (-360+np.nanmean(zcat['lon_GMC<deg>'][i]),np.nanmean(zcat['lat_GMC<deg>'][i])),
        np.nanmean(zcat['radius_acc_HI<deg>'][i]),color='black'))


ax.add_patch(Rectangle((-360,-5),720,10,color='black'))
ax.add_patch(Rectangle((-360,-50),45,100,color='black'))
ax.add_patch(Rectangle((360,-50),-45,100,color='black'))
ax.add_patch(Rectangle((-45,-50),90,100,color='black'))

    
#the gmclesses
ax.add_patch(Rectangle((-95,7),-28,27,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((-95+360,7),-28,27,color='black',fill=False,linewidth=3))

ax.add_patch(Rectangle((-40,10),-50,20,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((-40+360,10),-50,20,color='black',fill=False,linewidth=3))


ax.add_patch(Rectangle((-103,-38),-25,25,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((-103+360,-38),-25,25,color='black',fill=False,linewidth=3))


ax.add_patch(Rectangle((-52,-44),-38,29,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((-52+360,-44),-38,29,color='black',fill=False,linewidth=3))

ax.add_patch(Rectangle((197,-50),-48,22,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((197-360,-50),-48,22,color='black',fill=False,linewidth=3))

ax.add_patch(Rectangle((193,13),-25,22,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((193-360,13),-25,22,color='black',fill=False,linewidth=3))


ax.add_patch(Rectangle((90,10),-40,22,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((90-360,10),-40,22,color='black',fill=False,linewidth=3))


ax.add_patch(Rectangle((85,-35),-38,25,color='black',fill=False,linewidth=3))
ax.add_patch(Rectangle((85-360,-35),-38,25,color='black',fill=False,linewidth=3))


ax.set_xlim(360,0)


#     ax.set_xlim(140,40)

plt.show()


plt.hist(np.ravel(hi4pi_mom2),bins=1000,alpha=0.5)
plt.hist(np.ravel(dht_mom2),bins=1000,alpha=0.5)
# plt.xlim(0,40)
# plt.ylim(0,4000)
plt.loglog()
plt.show()


plt.scatter(np.ravel(hi4pi_mom2),np.ravel(whi_solm),s=5e-6,alpha=0.8)
plt.loglog()
plt.xlabel('small sigma km/s')
plt.ylabel('big sigma solmass/pc2')


whi_solm*(whi_solm40


plt.figure(figsize=(20,10))

plt.contourf(lons,lats,hi4pi_mom2,levels=np.arange(0.001,40,3),cmap=dc.cmap_min0('viridis'))
plt.contour(lons,lats,whi_solm,
            levels=[10,15,20],colors='black',linewidths=0.7)
plt.xlim(180,0)


# plt.ylim(-30,0)


plt.figure(figsize=(10,10))
plt.contourf(lons,lats,whi_solm,levels=np.arange(5,80,10),cmap=dc.cmap_min0('viridis'))
plt.xlim(260,180)
plt.ylim(-30,30)


plt.figure(figsize=(20,10))
plt.contourf(lons,lats,wco_solm,levels=[10,20],cmap=dc.cmap_min0('viridis'))
plt.xlim(360,0)
plt.ylim(-30,30)


plt.figure(figsize=(10,10))
plt.contourf(lons,lats,dht_mom2,levels=np.arange(0,10,1),cmap=dc.cmap_min0('viridis'))
plt.xlim(260,180)
plt.ylim(-30,30)


plt.figure(figsize=(10,10))
plt.contourf(lons,lats,wco_solm,levels=np.arange(0.001,70,0.5),cmap=dc.cmap_min0('viridis'))
plt.xlim(180,80)
plt.ylim(-30,30)


plt.figure(figsize=(10,10))
plt.contourf(lons,lats,dht_mom2,levels=np.arange(0,5,0.5),cmap=dc.cmap_min0('viridis'))
plt.xlim(180,80)
plt.ylim(-30,30)


[(plt.contourf(np.arange(zcat['lon_GMC<deg>'][i][0],zcat['lon_GMC<deg>'][i][1],-0.125),
             np.arange(zcat['lat_GMC<deg>'][i][0],zcat['lat_GMC<deg>'][i][1],0.125),
             fns.get_momN(h2s[i].data,zcat['v_pix_GMC'][i],2),
            levels=np.arange(0,10,0.5)),
plt.colorbar(),
plt.contour(np.arange(zcat['lon_GMC<deg>'][i][0],zcat['lon_GMC<deg>'][i][1],-0.125),
             np.arange(zcat['lat_GMC<deg>'][i][0],zcat['lat_GMC<deg>'][i][1],0.125),
             np.nansum(h2s[i].data,axis=0),
            levels=[2,5],norm='asinh',colors='white'),
 plt.show())
 for i in range(17)]


fits.HDUList([fits.PrimaryHDU(gmclesses[0])]
             +[fits.ImageHDU(gmclesses[i]) 
  for i in range(1,7)]).writeto('gmcless_tiles.fits',overwrite=True)


names = ['tile_'+np.arange(1,8).astype('str')[i] for i in range(7)]


help(QTable)


gless_cat = QTable((names,gmcless_ls,gmcless_bs),names=['name','lon<deg>','lat<deg>'])


gless_cat.write('gmcless_cat.fits',overwrite=True)


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


gmclesses = [tile1,tile2,tile3,tile4,tile5,tile6,tile7]


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
