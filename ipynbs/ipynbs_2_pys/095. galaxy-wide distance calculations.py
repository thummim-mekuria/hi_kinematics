#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import datacube as dc

from astropy.io import fits
import astropy.units as ux
import astropy.constants as c
from astrodendro import Dendrogram
from astropy.modeling import models


hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1]


wd = np.nansum(dust.data,0)
wh = np.nansum(hi4pi.data,0)
wc = np.nansum(dht.data,0)


lvd = np.nansum(dust.data,1)
lvh = np.nansum(hi4pi.data,1)
lvc = np.nansum(dht.data,1)


his = fits.open('hi_clouds.fits')
h2s = fits.open('h2_clouds.fits')


from astropy.table import QTable


zcat = QTable(fits.open('final_cat.fits')[1].data)


whole_gal_ds = pd.read_csv('whole_gal_dists.csv')


#start with chamaeleon
cham_d_lo,cham_d_li = ((zcat['lon_GMC<deg>'][1] - 180)/-0.23).astype('int')
cham_d_bo,cham_d_bi = ((zcat['lat_GMC<deg>'][1] + 90)/0.23).astype('int')


cham_d = dust.data[:,cham_d_bo:cham_d_bi,cham_d_lo:cham_d_li]


plt.imshow(np.nansum(cham_d,0),extent=[zcat['lon_GMC<deg>'][1][0],zcat['lon_GMC<deg>'][1][1],
                                      zcat['lat_GMC<deg>'][1][0],zcat['lat_GMC<deg>'][1][1]],
           origin='lower',cmap=dc.cmap_min0('viridis'))

plt.contour(np.arange(5760)*-0.125+360,np.arange(800)*0.125-50,
            wc,levels=[1.5,2,5],colors='orange',linewidths=0.7)

# plt.contour(np.arange(),np.arange(),
#             np.nansum(h2s[1].data,0),levels=[2])

plt.contour(np.arange(61)*-0.23+zcat['lon_GMC<deg>'][1][0],
            np.arange(35)*0.23+zcat['lat_GMC<deg>'][1][0],
             np.nansum(cham_d,0),levels=[0.2],colors='white',linewidths=0.7,)
plt.ylim(-20,-12)
plt.xlim(-54,-68)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.show()


cham_long_v = dht.data[:,zcat['b_pix_GMC'][1][0]:zcat['b_pix_GMC'][1][1],zcat['l_pix_GMC'][1][0]:zcat['l_pix_GMC'][1][1]]


cham_ds = whole_gal_ds[(whole_gal_ds['lon<deg>']>-66) & 
             (whole_gal_ds['lon<deg>'] < -54) & 
            (whole_gal_ds['lat<deg>'] > -20) & 
            (whole_gal_ds['lat<deg>'] < -12)]


cham_long_v.shape


plt.contourf(np.arange(112)*-0.125+zcat['lon_GMC<deg>'][1][0],
             np.arange(64)*0.125+zcat['lat_GMC<deg>'][1][0],
    np.nansum(cham_long_v,0)*1.97e20,origin='lower',norm='log',
             cmap=dc.cmap_min0('viridis'),levels=10**np.arange(19,23,0.3))

plt.scatter(cham_ds['lon<deg>'],cham_ds['lat<deg>'],
            c=cham_ds['d1<kpc>']*1000,vmin=50,vmax=500,s=100/cham_ds['d1<kpc>'],
            cmap='Set1',alpha=0.5)
plt.colorbar(label='d [pc]')

plt.xlim(-54,-68)

plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.show()


plt.imshow(np.nansum(cham_d,1),extent=[zcat['lon_GMC<deg>'][1][0]*10,zcat['lon_GMC<deg>'][1][1]*10,
                                      68,1250],origin='lower',cmap=dc.cmap_min0('viridis'))

plt.contour(np.arange(61)*-2.3+zcat['lon_GMC<deg>'][1][0]*10,np.arange(1182)+68,
             np.nansum(cham_d,1),levels=[0.2,0.5],colors='white',linewidths=0.7,)
plt.ylim(100,300)
plt.xlim(-550,-680)
plt.hlines(183,-550,-680,linestyle='--',color='orange')

plt.xlabel('l [deg x10]')
plt.ylabel('d [pc]')

plt.show()


plt.contourf(np.arange(112)*-0.125+zcat['lon_GMC<deg>'][1][0],np.arange(155)*1.3-100,
    np.nansum(cham_long_v,1)*1.97e20,origin='lower',norm='log',
             cmap=dc.cmap_min0('viridis'),levels=10**np.arange(19,23,0.3))

plt.scatter(cham_ds['lon<deg>'],cham_ds['vel<km/s>'],c=cham_ds['d1<kpc>']*1000,
            vmin=50,vmax=500,cmap='Set1',s=100/cham_ds['d1<kpc>'],alpha=0.5)
plt.colorbar(label='d [pc]')

plt.ylim(-25,25)

plt.hlines([-8,2,5],-68,-54,linestyle='--',colors=['red','orange','orange'])
plt.xlim(-54,-68)

plt.ylabel('v_lsr [km/s]')
plt.xlabel('l [deg]')

plt.show()


plt.plot(np.arange(155)*1.3-100,np.nanmean(cham_long_v,(1,2)))
plt.xlim(-25,25)
plt.xlabel('vel[km/s]')
plt.ylabel('T [K]')
plt.vlines([-8,2,5],0,0.35,linestyle='--',colors=['red','orange','orange'])
plt.show()


ors_d_lo,ors_d_li = ((np.array([-135,-175]) - 180)/-0.23).astype('int')
ors_d_bo,ors_d_bi = ((np.array([-25,-5]) + 90)/0.23).astype('int')


ors_d = dust.data[:,ors_d_bo:ors_d_bi,ors_d_lo:ors_d_li]


ors_ds = whole_gal_ds[(whole_gal_ds['lon<deg>']>-175) & 
             (whole_gal_ds['lon<deg>'] < -135) & 
            (whole_gal_ds['lat<deg>'] > -25) & 
            (whole_gal_ds['lat<deg>'] < -5)]


ors_long_v = dht.data[:,int((-25+50)/0.125):int((-5+50)/0.125),
                       int((-135-360)/-0.125):int((-175-360)/-0.125)]


plt.imshow(np.nansum(ors_d,0),extent=[-135,-175,-25,-5],
           origin='lower',cmap=dc.cmap_min0('viridis'))

# plt.contour(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wc,norm='log',
#              levels=[1.5,2,5],linewidths=0.6,colors='orange') 


plt.contour(np.arange(174)*-0.23-135, np.arange(87)*0.23-25,
             np.nansum(ors_d,0),levels=[0.2,0.5,1,2],colors='white',linewidths=0.7,)
plt.ylim(-25,-5)
plt.xlim(-135,-175)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.show()


plt.contourf(np.arange(320)*-0.125-130,
             np.arange(160)*0.125-25,
    np.nansum(ors_long_v,0)*1.97e20,origin='lower',norm='log',
             cmap=dc.cmap_min0('viridis'),levels=10**np.arange(19,23,0.3))

plt.scatter(ors_ds['lon<deg>'],ors_ds['lat<deg>'],
            c=ors_ds['d1<kpc>']*1000,vmin=50,vmax=500,s=10/ors_ds['d1<kpc>'],
            cmap='Set1',alpha=0.5)
plt.colorbar(label='d [pc]',orientation='horizontal')

plt.xlim(-135,-174)
plt.ylim(-25,-5)

plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.show()


plt.figure(figsize=(7,5))

plt.contourf(np.arange(320)*-0.125-135,np.arange(155)*1.3-100,
    np.nansum(ors_long_v,1)*1.97e20,origin='lower',norm='log',
             cmap=dc.cmap_min0('viridis'),levels=10**np.arange(19,23,0.3))

plt.scatter(ors_ds['lon<deg>'],ors_ds['vel<km/s>'],c=ors_ds['d1<kpc>']*1000,
            vmin=300,vmax=600,cmap='Set1',s=10/ors_ds['d1<kpc>'],alpha=0.4)
plt.colorbar(label='d [pc]')

plt.ylim(-15,30)

plt.hlines([4,10],-135,-175,linestyle='--',colors='orange')
plt.xlim(-135,-175)

plt.ylabel('v_lsr [km/s]')
plt.xlabel('l [deg]')

plt.show()


plt.plot(np.arange(155)*1.3-100,np.nanmean(ors_long_v,(1,2)))
plt.vlines([4,10],0,0.2,colors=['red','orange'],linestyle='--')
plt.xlim(-15,30)
plt.xlabel('v [km/s]')
plt.ylabel('T [K]')
plt.show()


plt.imshow(np.nansum(ors_d,1),extent=[-1350,-1750,
                                      68,1250],origin='lower',cmap=dc.cmap_min0('viridis'))

plt.contour(np.arange(174)*-2.3-1350,np.arange(1182)+68,
             np.nansum(ors_d,1),levels=[0.05,0.1,0.2,0.5],colors='white',linewidths=0.6,)
plt.ylim(300,1250)
# plt.ylim(600,1000)

plt.xlim(-1350,-1750)
plt.hlines([423,432,402,886,778],-1350,-1750,linestyle='--',linewidths=0.7,color='orange')

plt.xlabel('l [degx10]')
plt.ylabel('d [pc]')

plt.show()


zcat


systematically compare the kinematics of the molecular & atomic gas

#read in more detail imara&burkhart pdfs, imara&blitz angular momentum
#i&b used extinction to define the clouds

questions to ask:
    1. are the (kinematic etc) trends related to the sfr within the cloud
    2. how does this compare to simulations relating co<->hi
    
    
for next time:
    1. first moment maps + fits to the gradient 
        for measurements of specfic angular momentum or magnitude of velocity gradient


1/0.083 #would be better to go back to this


1/0.125


local = whole_gal_ds[whole_gal_ds['d1<kpc>']<2]


plt.hist(local['d1<kpc>']/local['d1_unc<kpc>'],density=True,alpha=0.5,
         bins = np.arange(0,15),label='local (d<2kpc)')
plt.hist(whole_gal_ds['d1<kpc>']/whole_gal_ds['d1_unc<kpc>'],density=True,
         bins = np.arange(0,15),alpha=0.5,label='all')

plt.xlabel('SNR')
plt.yscale('log')
plt.legend()
plt.show()


local_snr10 = local[local['d1<kpc>']/local['d1_unc<kpc>'] > 10]


plt.figure(figsize=(20,8))

plt.contourf(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wh,
            levels=[400,1e5],colors='green')      

plt.contourf(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wc,
            levels=[1,1000],colors='black',) 

plt.scatter(local_snr10['lon<deg>'],local_snr10['lat<deg>'],cmap='Set1',alpha=0.5,
            s=10/local_snr10['d1<kpc>'],c = local_snr10['d1<kpc>']*1000)
plt.colorbar(orientation='horizontal',label='d [pc]')

plt.xlim(180,-180)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.show()


plt.figure(figsize=(20,8))

plt.scatter(local['lon<deg>'],local['d1<kpc>']/local['d1_unc<kpc>'],
            s=10,c=local['d1<kpc>']*1000,cmap='Set1',alpha=0.5)
plt.colorbar(orientation='horizontal',label='d [pc]')
plt.yscale('log')
plt.xlim(180,-180)
plt.hlines([1,5,10],180,-180,linewidth=4,linestyle=[':','--','solid'],
           colors='black')
plt.ylabel('snr')
plt.xlabel('lon<deg>')
plt.show()


0.23/0.083


plt.figure(figsize=(15,7))

# plt.contourf(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wh,
#             levels=[100,400,900,1500,3000,10000],linewidths=0.9,norm='asinh')      

plt.contour(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wc,norm='log',
             levels=[3],linewidths=0.9,colors='orange') 

plt.contourf(np.arange(1572)*-0.23 + 180,np.arange(786)*0.23 - 90,
            wd,levels=[0.05,0.1,0.2,0.3,0.5,0.7,1],linewidths=0.9)


plt.xlim(180,-180)
plt.ylim(-50,50)


plt.figure(figsize=(20,10))
plt.contourf(np.arange(1572)*-0.23 + 180,np.arange(786)*0.23 - 90,wd,levels=[0.05,0.1,0.2,0.3])
             
plt.xlim(180,-180)
plt.contour(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wh,
            levels=[100,400,900],colors='red',linestyles='--',linewidths=0.9)      


plt.contour(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wc,
            levels=[1,10],colors='orange',linestyles='--',linewidths=0.9) 

#              ,origin='lower',extent=[180,-180,-90,90],norm='log',
#            cmap=dc.cmap_min0('viridis'))


plt.figure(figsize=(23,10))

plt.contourf(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wh,
            levels=[100,400,900,1500,3000,10000],linewidths=0.9,norm='asinh')      

# plt.contour(np.arange(5760)*-0.125 + 360,np.arange(800)*0.125 - 50,wc,
#             levels=[3],colors='white',linewidths=0.5) 

plt.contour(np.arange(1572)*-0.23 + 180,np.arange(786)*0.23 - 90,
            wd,levels=[0.2],linewidths=0.9,colors='white')


plt.xlim(180,-180)
plt.ylim(-50,50)


plt.figure(figsize=(20,20))
plt.imshow(lbd,extent=[180,-180,68,1250],norm='log',
           origin='lower',cmap=dc.cmap_min0('viridis'))
# plt.contour(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
#             lbd,levels=np.arange(1e-4,2,5e-4),linewidths=0.9,cmap='viridis',norm='log')
# plt.yscale('log')


lbd.shape


dlbd = np.gradient(lbd)


plt.figure(figsize=(20,20))
plt.imshow(dlbd[0],extent=[180,-180,68,1250],norm='log',
           origin='lower',cmap=dc.cmap_min0('viridis'))
# plt.contour(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
#             lbd,levels=np.arange(1e-4,2,5e-4),linewidths=0.9,cmap='viridis',norm='log')
# plt.yscale('log')


plt.figure(figsize=(20,20))
plt.imshow(dlbd[1],extent=[180,-180,68,1250],norm='log',
           origin='lower',cmap=dc.cmap_min0('viridis'))
# plt.contour(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
#             lbd,levels=np.arange(1e-4,2,5e-4),linewidths=0.9,cmap='viridis',norm='log')
# plt.yscale('log')


plt.figure(figsize=(15,10))
plt.imshow(dlbd[1],extent=[180,-180,68,1250],norm='log',origin='lower')
plt.contour(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
            lbd,levels=np.arange(1e-2,2,5e-2),linewidths=0.5,
            norm='log',colors='red')
# plt.yscale('log')
# plt.xscale('log')

# plt.loglog()


plt.figure(figsize=(20,20))
plt.imshow(dlbd[0],extent=[180,-180,68,1250],norm='log',origin='lower')
plt.contour(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
            lbd,levels=np.arange(1e-2,2,5e-2),linewidths=0.5,
            norm='log',colors='red')
# plt.yscale('log')
# plt.xscale('log')

# plt.loglog()


plt.contour(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
            lbd,levels=np.arange(1e-2,2,5e-2),linewidths=0.9,cmap='viridis',norm='log')
plt.yscale('log')


plt.contourf((dlbd[0]/np.nanmax(dlbd[0]))*(dlbd[1]/np.nanmax(dlbd[1])),
             norm='log',levels=np.arange(0.01,1,0.1))
# plt.colorbar()


bdd = np.nansum(dust.data,2)


bvh = np.nansum(hi4pi.data,2)


plt.figure(figsize=(23,10))

# plt.contourf(np.arange(800)*0.125 -50,np.arange(155)*1.3 - 100,bvh,
#             levels=[100,400,900,1500,3000,10000],linewidths=0.9,norm='asinh')      

plt.contourf(np.arange(1572)*-0.23 + 180,np.arange(786)*0.23 - 90,wd,
            levels=[0.05,0.1,0.5,1,2],cmap='plasma',linewidths=0.5) 

# plt.contourf(np.arange(786)*0.23 - 90,np.arange(1182) + 68,
#             bdd,levels=[0.1,0.5,1,2],linewidths=0.9,cmap='plasma')

# plt.yscale('log')
plt.xlim(180,-180)
# plt.ylim(-50,50)


plt.figure(figsize=(23,10))

plt.contourf(np.arange(1572)*-0.23 +180,np.arange(1182) + 68,
            lbd,levels=[0.05,0.1,0.5,1,2],linewidths=0.9,cmap='plasma')

plt.yscale('log')
# plt.xlim(180,-180)
# plt.ylim(-50,50)


plt.figure(figsize=(23,10))

plt.contourf(np.arange(786)*0.23 - 90,np.arange(1182) + 68,
            bdd,levels=[0.05,0.1,0.5,1,2],linewidths=0.9,cmap='plasma')

plt.yscale('log')
# plt.xlim(180,-180)
# plt.ylim(-50,50)


long_dust = np.concatenate((dust.data,dust.data,dust.data),axis=2)


wd_long = np.nansum(long_dust,0)
ldd_long = np.nansum(long_dust,1)
bdd_long = np.nansum(long_dust,2)


np.arange(1572)*-0.23 + 180 = -120


(-120 -180)/-0.23


0.229*u.deg.to(u.rad)*1250 #pc/pixel at most coarse resolution 


0.229*u.deg.to(u.rad)*68 #pc/pixel at leaste coarse resolution 


10/0.229*10/0.229*10/4.996 #the size of a 10pc3 cloud at 1250pc


10/0.229*10/0.229*10/0.272 #the size of a 10pc3 cloud at 68pc


dustD = Dendrogram.compute(dust.data,min_value=1e-3,min_delta=1e-2,min_npix=2000)
