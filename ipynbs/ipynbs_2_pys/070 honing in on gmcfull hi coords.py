#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable


zcat = QTable(fits.open('final_cat.fits')[1].data)


hi4pi = fits.open('HI4PI_shorter_resampled.fits')[0]
dht = fits.open('DHT_CO_longer.fits')[0]


hicat = QTable(fits.open('hi_table.fits')[1].data)


vels = np.arange(155)*1.3 - 100.
lons = np.arange(2881)*(-0.125) + 180.125


colats = np.arange(481)*(0.125) - 30.125
hilats = np.arange(800)*(0.125) - 50


l0s,l1s = np.array([((zcat['lon<deg>'][i]-180.125)/-0.125).astype('int') for i in range(19)]).transpose()
lcs = np.array([((np.nanmean(zcat['lon<deg>'][i])-180.125)/-0.125).astype('int') for i in range(19)])
lxs = np.array([((np.nanstd(zcat['lon<deg>'][i]))/0.125).astype('int') for i in range(19)])
v0s,v1s = np.array([((zcat['vel<km/s>'][i]-100.)/1.3).astype('int') for i in range(19)]).transpose()
vcs = np.array([(np.nanmean(zcat['vel<km/s>'][i]-100.)/1.3).astype('int') for i in range(19)])
vxs = np.array([(np.nanstd(zcat['vel<km/s>'][i])/1.3).astype('int') for i in range(19)])


bcsc = np.array([((np.nanmean(zcat['lat<deg>'][i])-30.125)/0.125).astype('int') for i in range(19)])
bcsh = np.array([((np.nanmean(zcat['lat<deg>'][i])-50.)/0.125).astype('int') for i in range(19)])
bxs = np.array([((np.nanstd(zcat['lat<deg>'][i]))/0.125).astype('int') for i in range(19)])
b0sh,b1sh = np.array([((zcat['lat<deg>'][i]-50.)/0.125).astype('int') for i in range(19)]).transpose()
b0sc,b1sc = np.array([((zcat['lat<deg>'][i]-30.125)/0.125).astype('int') for i in range(19)]).transpose()


co_cubes = [dht.data[v0s[i]:v1s[i],b0sc[i]:b1sc[i],l0s[i]:l1s[i]] for i in range(19)]


2881/2


np.concatenate((lons[-5:],lons[2:4]))


[(plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[i]:v1s[i],:,1441:-1],\
                             hi4pi.data[v0s[i]:v1s[i],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=600,extent=[360,-180,-50,50]),
  plt.title(zcat['names'][i]),
  
  plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),\
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[i]:v1s[i],:,1441:-1],\
                             hi4pi.data[v0s[i]:v1s[i],:,2:]),axis=2),axis=0),\
    levels=np.arange(100,390,10)),
  
  plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,20,3),cmap='binary_r',linewidths=0.7),
  
  plt.xlim(210,150),
  plt.ylim(hilats[bcsh[i]-10*bxs[i]],hilats[bcsh[i]+10*bxs[i]]),
  
  plt.show(),
  
  plt.imshow(np.nansum(\
                       np.concatenate((hi4pi.data[:,b0sh[i]:b1sh[i],1441:-1],\
                                       hi4pi.data[:,b0sh[i]:b1sh[i],2:]),axis=2),\
                                      axis=1),origin='lower',vmax=3000,extent=[360,-180,-100,100]),
  
  plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),\
              vels,np.nansum(np.concatenate((hi4pi.data[:,b0sh[i]:b1sh[i],1441:-1],\
                                       hi4pi.data[:,b0sh[i]:b1sh[i],2:]),axis=2),axis=1),\
    levels=np.arange(800,5000,130)),
  
  plt.contour(lons[l0s[i]:l1s[i]],vels[v0s[i]:v1s[i]],np.nansum(co_cubes[i],axis=1),\
    levels=np.arange(50,200,5),cmap='binary_r',linewidths=0.7),
  
  
  plt.xlim(190,160),
  plt.ylim(vels[vcs[i]-2*vxs[i]],vels[vcs[i]+2*vxs[i]]),
  
  plt.show(),
  
  
  plt.imshow(np.nansum(hi4pi.data[:,:,l0s[i]:l1s[i]],axis=2),\
    origin='lower',vmax=7900,extent=[-50,50,-100,100]),
  
  plt.contour(hilats,vels,np.nansum(hi4pi.data[:,:,lcs[i]-lxs[i]:lcs[i]+lxs[i]],axis=2),\
    levels=np.arange(3300,5800,190)),
  
  plt.contour(colats[b0sc[i]:b1sc[i]],vels[v0s[i]:v1s[i]],np.nansum(co_cubes[i],axis=2),\
    levels=np.arange(100,360,90),cmap='binary_r',linewidths=0.7),
  
  
  plt.xlim(hilats[bcsh[i]-4*bxs[i]],hilats[bcsh[i]+3*bxs[i]]),
  plt.ylim(vels[vcs[i]-2*vxs[i]],vels[vcs[i]+2*vxs[i]]),
  
  plt.show(),
 
 
 ) \
 for i in [18]]#,18]]


[(plt.figure(figsize=(20,7)),\
  plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[i]:v1s[i],:,1441:-1],\
                             hi4pi.data[v0s[i]:v1s[i],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=600,extent=[360,-180,-50,50]),
  plt.title(zcat['names'][i]),
  
  plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),\
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[i]:v1s[i],:,1441:-1],\
                             hi4pi.data[v0s[i]:v1s[i],:,2:]),axis=2),axis=0),\
    levels=np.arange(100,500,10)),
  
  plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,20,3),cmap='binary_r',linewidths=0.7),
  
#   plt.xlim(210,150),
#   plt.ylim(hilats[bcsh[i]-10*bxs[i]],hilats[bcsh[i]+10*bxs[i]]),
  
  plt.show(),
  
  
 ) \
 for i in range(19)]


zcat


plt.figure(figsize=(20,7))
hi = 1
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=1300,extent=[360,-180,-50,50]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),axis=0),\
    levels=(np.arange(400,1050,90)))
  
for i in [1]:
    plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,13,3),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-40,5)

plt.show()


plt.figure(figsize=(20,7))
hi = 15
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=900,extent=[360,-180,-50,50]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),axis=0),\
    levels=np.arange(50,550,70))
  
for i in [15,18]:
    plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,13,3),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-40,5)

plt.show()


plt.figure(figsize=(20,7))
hi = 18
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=900,extent=[360,-180,-50,50]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),axis=0),\
    levels=np.arange(50,450,70))
  
for i in [15,18]:
    plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,13,3),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-40,5)

plt.show()


plt.figure(figsize=(20,7))
hi = 15
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=900,extent=[360,-180,-50,50]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[hi]:v1s[hi],:,1441:-1],\
                             hi4pi.data[v0s[hi]:v1s[hi],:,2:]),axis=2),axis=0),\
    levels=np.arange(100,450,40))
  
for i in [15,18]:
    plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,13,3),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-40,5)

plt.show()


plt.figure(figsize=(20,7))
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[:,b0sh[15]:b1sh[15],1441:-1],\
                             hi4pi.data[:,b0sh[15]:b1sh[15],2:]),axis=2),\
             axis=1), origin='lower',vmax=2000,extent=[360,-180,-100,100]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              vels,np.nansum(\
                               np.concatenate((hi4pi.data[:,b0sh[15]:b1sh[15],1441:-1],\
                             hi4pi.data[:,b0sh[15]:b1sh[15],2:]),axis=2),axis=1),\
    levels=np.arange(800,4500,250))
  
for i in [15,18]:
    plt.contour(lons[l0s[i]:l1s[i]],vels[v0s[i]:v1s[i]],np.nansum(co_cubes[i],axis=1),\
    levels=np.arange(30,130,20),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-10,20)

plt.show()


plt.figure(figsize=(20,7))
plt.imshow(np.nansum(hi4pi.data[:,:,l0s[i]:l1s[i]],axis=2),\
    origin='lower',vmax=7900,extent=[-50,50,-100,100]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(hilats,vels,np.nansum(hi4pi.data[:,:,lcs[i]-lxs[i]:lcs[i]+lxs[i]],axis=2),\
    levels=np.arange(3300,5800,190))
  
for i in [15,18]:
    plt.contour(colats[b0sc[i]:b1sc[i]],vels[v0s[i]:v1s[i]],np.nansum(co_cubes[i],axis=2),\
    levels=np.arange(30,130,20),cmap='binary_r',linewidths=0.7),
  
plt.xlim(-40,0)
plt.ylim(-10,20)


plt.figure(figsize=(20,7))
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[18]:v1s[18],:,1441:-1],\
                             hi4pi.data[v0s[18]:v1s[18],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=900,extent=[360,-180,-50,50]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[18]:v1s[18],:,1441:-1],\
                             hi4pi.data[v0s[18]:v1s[18],:,2:]),axis=2),axis=0),\
    levels=np.arange(100,450,40))
  
for i in [1,15,18]:
    plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,13,3),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-40,5)

plt.show()


plt.figure(figsize=(20,7))
plt.imshow(np.nansum(\
             np.concatenate((hi4pi.data[v0s[1]:v1s[1],:,1441:-1],\
                             hi4pi.data[v0s[1]:v1s[1],:,2:]),axis=2),\
             axis=0), origin='lower',vmax=900,extent=[360,-180,-50,50]),
# plt.title(zcat['names'][15]+' and '+zcat['names'][18])
  
plt.contour(np.concatenate((lons[1441:-1]+360,lons[2:])),
              hilats,np.nansum(\
                               np.concatenate((hi4pi.data[v0s[1]:v1s[1],:,1441:-1],\
                             hi4pi.data[v0s[1]:v1s[1],:,2:]),axis=2),axis=0),\
    levels=np.arange(100,450,40))
  
for i in [1,15,18]:
    plt.contour(lons[l0s[i]:l1s[i]],colats[b0sc[i]:b1sc[i]],np.nansum(co_cubes[i],axis=0),\
    levels=np.arange(5,13,3),cmap='binary_r',linewidths=0.7),
  
plt.xlim(210,140)
plt.ylim(-40,5)

plt.show()


plt.plot(np.nanmean(co_cubes[18],axis=(1,2)))
plt.plot(np.nanmean(co_cubes[15],axis=(1,2)))
plt.plot(np.nanmean(co_cubes[1],axis=(1,2)))


zcat['vel<km/s>']


camerlopardalis   L = [157.,133.]    B = [12.,26.],   V = [-5.,5.]
chamaeleon        L = [-51.,-68.]    B = [-23.,-7.],  V = [-8.,10.]
crossbones        L = [-135.,-145.]  B = [-15.,-3.],  V = [3.,22.]
lacerta           L = [96.,91.],     B =[-25.,-7.],   V = [-10.,10.]
maddalena         L = [-138.,-148.]  B = [-6,1],      V = [22.,43.]
oriona            L = [-140.,-160.]  B = [-21.,-11.], V = [-3.,20.]
orionb            L = [-151.,-162.]  B = [-21.,-11.], V = [-12.,10.]
aquarius rift     L = [36,10],       B = [2.,14.],    V = [0.,10.]
california        L = [170.,149.]    B = [-15., -6.]  V = [-23., 11.]
cma               L = [-127.,-142.]  B = [-8.,1.]     V = [6.,25.]
hercules          L = [40.,60.]      B = [7.,13.]     V = [0.,20.]
mon ob1           L = [-152.,-164.]  B = [1.,6.]      V = [-2.,15.,]
mon r2            L = [-142.,-150.]  B = [-16.,-9.]   V = [-1.,20.]
orionl            L = [-155.,-172.]  B = [-20.,-8.]   V = [-7.,18]
perseus           L = [169.,150.]    B = [-28.,-14.]  V = [-5.,12.,]
polaris           L = [127.,114.]    B = [20.,35.]    V = [-8.,3.]
taurus            L = [174.,162.], B [-25.,-8.] V = [0.,15.]
gemini            L = [-160.,-174.]  B = [-10.,5.]    V = [-5.,25.]


#when i come back, experiment with just looking at HI overdensities, without? (or with many?) co clouds contoured


zcat['names']


clouds that seem to be in the valley between two peaks rather than on top of one


allthewayaround_dht = np.concatenate((dht.data[:,:,1441:],dht.data,dht.data[:,:,:1441]),axis=2)


allthewayaround_hi4pi = np.concatenate((hi4pi.data[:,:,1441:],hi4pi.data,hi4pi.data[:,:,:1441]),axis=2)


hdu_d = fits.PrimaryHDU(allthewayaround_dht)
hdu_h = fits.ImageHDU(allthewayaround_hi4pi)


fits.HDUList([hdu_d,hdu_h]).writeto('both_same_grid_720d.fits',overwrite=True)


hi4pi.data.shape


np.concatenate((lons[1441:]+360,lons,lons[:1441]-360))


lons[1441:]+360


lons[1:-1]


lons[:1441]-360


# hi4pi = fits.open('HI4PI_shorter.fits')[0]
# ltrpt = np.repeat(hi4pi.data,2,axis=1) #when hi4pi was just HI4PI_shorter.fits
# lnrpt = np.repeat(ltrpt,2,axis=2)
# thisoneqm = lnrpt[:,::3,::3]
# hdu = fits.PrimaryHDU(thisoneqm[:,:,:-1])
# hdu.writeto('HI4PI_shorter_resampled.fits',overwrite=True)
