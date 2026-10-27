#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import QTable
import astropy.units as u
import astropy.constants as c

import matplotlib.pyplot as plt

import datacube as dc


dht,hi4pi = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


for i,ii in enumerate(list(zcat['vel<km/s>'])):
#     print(ii[0])
    plt.hlines(zcat['names'][i],ii[0],ii[1])
    
plt.vlines([-20,-10,0,10,20,30,40],'aquila_r','taurus',alpha=0.5)
plt.vlines([-15,-5,5,15,25,35],'aquila_r','taurus',alpha=0.5,linestyle='--')

plt.show()


0:15 -- perseus,mob1
0:20 -- hercules, monr2, 
-5:20 -- orions
20:30 -- maddalena
-7,0 -- polaris


(20 + 100)/1.3, (30 + 100)/1.3


#maddalena
mad_co = dht.data[93:100,:,:]
mad_hi = hi4pi.data[93:100,:,:]


mad_hi.shape


lons = np.arange(5762)*(-0.125) + 360

hilats = np.arange(800)*0.125 - 50
colats = np.arange(481)*0.125 - 30


[100,200,400,500,700,1000]


hi_conv = (1.3*u.km/u.s *1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p *u.K).to(u.solMass/u.pc**2) #last term is data


plt.figure(figsize=(10,10))
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(hi_conv.value*np.nansum(mad_hi,axis=0),extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=15,vmin=0)

plt.title('maddalena')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(-130,-170)
plt.ylim(-10,10)


plt.contour(lons,colats,np.nansum(mad_co,axis=0),levels=[1.5,3,7],colors=['black'],linewidths=[1,0.7,0.5])
plt.contour(lons,hilats,hi_conv.value*np.nansum(mad_hi,axis=0),\
            levels=[5,7,10,12,20],linewidths=[0.7,0.7,1,0.7,0.7],\
            colors=['white'],linestyles=['--','--','solid','--','--'])


(0 + 100)/1.3, (15 + 100)/1.3


(-5 + 100)/1.3, (20 + 100)/1.3


#perseus & ob1
persob1_co = dht.data[76:85,:,:]
persob1_hi = hi4pi.data[76:85,:,:]


#perseus & ob1
ors_co = dht.data[73:92,:,:]
ors_hi = hi4pi.data[73:92,:,:]


plt.figure(figsize=(10,10))

wco = np.nansum(ors_co,axis=0)
whi = hi_conv.value*np.nansum(ors_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=25,vmin=0)

plt.title('orions')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(250,130)
plt.ylim(-45,15)


plt.contour(lons,colats,wco,levels=[4],colors=['black'],linewidths=[1])#,0.7,0.5])
plt.contour(lons,hilats,whi,\
            
             levels=[10,15],linewidths=[1],colors=['white'])#,\
#             colors=['white'],linestyles=['--','solid','--','--','--'])

            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


plt.figure(figsize=(10,10))

wco = np.nansum(persob1_co,axis=0)
whi = hi_conv.value*np.nansum(persob1_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=15,vmin=0)

plt.title('perseus')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(230,150)
plt.ylim(-45,15)


plt.contour(lons,colats,wco,levels=[4],colors=['black'],linewidths=[1])#,0.7,0.5])
plt.contour(lons,hilats,whi,\
            
             levels=[10,15],linewidths=[1],colors=['white'])#,\
#             colors=['white'],linestyles=['--','solid','--','--','--'])

            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


(-7 + 100)/1.3, (0 + 100)/1.3


#polaris
pol_co = dht.data[71:77,:,:]
pol_hi = hi4pi.data[71:77,:,:]


plt.figure(figsize=(10,10))

wco = np.nansum(pol_co,axis=0)
whi = hi_conv.value*np.nansum(pol_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=5,vmin=0)

plt.title('polaris')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(140,105)
plt.ylim(15,45)


plt.contour(lons,colats,wco,levels=[1,2,3],colors=['black'],linewidths=[1,0.7,0.5])
plt.contour(lons,hilats,whi,\
            
             levels=[2.5],linewidths=[1],colors=['white'])#,\
#             colors=['white'],linestyles=['--','solid','--','--','--'])

            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


for i,ii in enumerate(list(zcat['vel<km/s>'])):
#     print(ii[0])
    plt.hlines(zcat['names'][i],ii[0],ii[1])
    
plt.vlines([-20,-10,0,10,20,30,40],'aquila_r','taurus',alpha=0.5)
plt.vlines([-15,-5,5,15,25,35],'aquila_r','taurus',alpha=0.5,linestyle='--')

plt.show()


(-3 + 100)/1.3, (10 + 100)/1.3


#chamaeleon
cham_co = dht.data[74:87,:,:]
cham_hi = hi4pi.data[74:87,:,:]


plt.figure(figsize=(10,10))

wco = np.nansum(cham_co,axis=0)
whi = hi_conv.value*np.nansum(cham_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=10,vmin=0)

plt.title('chamaeleon')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(350,260)
plt.ylim(-45,-5)


plt.contour(lons,colats,wco,levels=[1,2,3],colors=['black'],linewidths=[1,0.7,0.5])
plt.contour(lons,hilats,whi,\
            
             levels=[5.4],linewidths=[1,0.7,0.5],colors=['white'])
            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


(-15 + 100)/1.3, (23 + 100)/1.3


#gemini
gem_co = dht.data[65:95,:,:]
gem_hi = hi4pi.data[65:95,:,:]


plt.figure(figsize=(10,10))

wco = np.nansum(gem_co,axis=0)
whi = hi_conv.value*np.nansum(gem_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=50,vmin=0)

plt.title('gemini')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(205,165)
plt.ylim(-10,10)


plt.contour(lons,colats,wco,levels=[3,5,7],colors=['black'],linewidths=[1,0.5,0.3])
plt.contour(lons,hilats,whi,levels=[30],linewidths=[1,0.7,0.5],colors=['white'])
# plt.contour(lons,hilats,hi_conv.value*np.nansum(persob1_hi,axis=0),levels=[10],\
#             linewidths=[1,0.7,0.5],colors=['red'])
           
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


plt.figure(figsize=(10,10))

wco = np.nansum(persob1_co,axis=0)
whi = hi_conv.value*np.nansum(persob1_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=30,vmin=0)

plt.title('mon ob1')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(210,195)
plt.ylim(-5,5)


plt.contour(lons,colats,wco,levels=[3,5,10],colors=['black'],linewidths=[1,0.5,0.3])
plt.contour(lons,hilats,whi,\
            
             levels=[10.5],linewidths=[1,0.7,0.5],colors=['white'])
            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


(-15 + 100)/1.3, (23 + 100)/1.3


#rosette
ros_co = dht.data[65:95,:,:]
ros_hi = hi4pi.data[65:95,:,:]


plt.figure(figsize=(10,10))

wco = np.nansum(ros_co,axis=0)
whi = hi_conv.value*np.nansum(ros_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=50,vmin=0)

plt.title('rosette')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(220,200)
plt.ylim(-10,10)


plt.contour(lons,colats,wco,levels=[3,5,10],colors=['black'],linewidths=[1,0.5,0.3])
plt.contour(lons,hilats,whi,\
            
             levels=[29],linewidths=[1,0.7,0.5],colors=['white'])
            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


(0 + 100)/1.3, (17 + 100)/1.3


#aquila
aqr_co = dht.data[76:90,:,:]
aqr_hi = hi4pi.data[76:90,:,:]


plt.figure(figsize=(10,10))

wco = np.nansum(aqr_co,axis=0)
whi = hi_conv.value*np.nansum(aqr_hi,axis=0)
# plt.imshow(np.nansum(mad_co,axis=0),extent=[360,-360,-30,30],origin='lower',\
#            vmax=10,cmap=dc.cmap_min0('viridis'))


plt.imshow(whi,extent=[360,-360,-50,50],origin='lower',\
           cmap=dc.cmap_min0('viridis'),vmax=20,vmin=0)

plt.title('aquila rift')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')

plt.colorbar(orientation='horizontal',label=r'$\Sigma_{HI}\ [M_\odot \ pc^{-2}]$',)

plt.xlim(50,5)
plt.ylim(-20,30)


plt.contour(lons,colats,wco,levels=[5,10],colors=['black'],linewidths=[1,0.5,0.3])
plt.contour(lons,hilats,whi,\
            
             levels=[15],linewidths=[1,0.7,0.5],colors=['white'])
            
#             levels=[5,7,10,12,20],linewidths=[0.7,0.7,,0.7,0.7],\
#             colors=['white'],linestyles=['--','solid','--','--','--'])


plt.imshow(np.nansum(hi4pi.data,axis=0),origin='lower',vmin=400,vmax=900,cmap=dc.cmap_min0('viridis'))
plt.xlim(0,2000)


zcat
