#!/usr/bin/env python
# coding: utf-8

# for meeting with nia -> show histograms of first moment maps
# - co on average is at more positive velocities than hi
# - distinct populations (map to their own separate spatially-correlated structures too for the most part)
# - velocity x distance -> the parts that are distinct in the mom 1 histogram tend to also show up as separated peaks in ld/bd plots
# - velocity gradient x flows (esp perseus, california)
# - random distribution (corona? orion b?) -> turbulence?
# - rotation (orion a)
# - read burkert & bodenheimer
# - summary slides for b&b, a-c&v-s
# - write description for accumulation radius method? -> look more into moment 1 histograms for defining hi?
# - similar thing with linewidths (distinct populations, spatially separate)
# - in order for this ^ to be robust though, the velocity channels over which we're integrating have to be a bit more agnostic
#     - did this for first moment maps over the whole 155 channels, and the hi velocities are def way more extreme than before (also there are systematic shifts that seem to be following the galactic rotation curve (negative-wards for clouds with l>60, positive-wards for clouds with l<-60)
#     - trend is the same whether we're looking at hi close in (within co bounds) or further out (accumulation radius)
#     
#     - exceptions -> ophiuchus, the orions
#     - this seems so at odds with the spectra peaking near each other ...
#     - when comparing with spectra, the first moment map histograms catch a lot more information at the tails -> not sure why ...  but it def contributes to having it be harder to distentangle los stuff from spectra
#  - look into why the mom 1 over the long range or so so different than over the restricted one
#  - and also if theres some way of quantifying how much material is captured / identified as distinct from spectra vs mom1 hist
# 


h1 = fns.get_mom1(hi4pi.data,np.arange(155)*1.3-100)


from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import numpy as np
import pandas as pd
import datacube as dc
import fun as fns

import matplotlib.pyplot as plt


hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data
cat = QTable(fits.open('post_flash_cat.fits')[1].data)


for i in range(22):
    ccld0 = fits.open('mom0_co.fits')[i].data
    hcld0 = fits.open('mom0_hi.fits')[i].data
    plt.hist(np.log10(np.ravel(ccld0)),
             bins=np.arange(19,24,0.1),
             alpha=0.5,density=True)
    plt.hist(np.log10(np.ravel(hcld0)),
             bins=np.arange(19,24,0.1),
             alpha=0.5,density=True)

    plt.title(str(i)+' '+cat['name'][i])
    plt.show()


for i in [21]:
    cld0 = fits.open('mom0_co.fits')[i].data
    plt.contourf(np.log10(cld0),levels=np.arange(20,20.7,0.1),cmap='winter')
    plt.colorbar()
    plt.contourf(np.log10(cld0),levels=np.arange(20.7,22.5,0.1),cmap='magma')
    plt.colorbar()


1e3*u.AU.to(u.pc)


for i in range(22):
    ccld2 = fits.open('mom2_co.fits')[i].data
    hcld2 = fits.open('mom2_hi.fits')[i].data
    plt.hist(ccld2[ccld2>0.5],
             bins=np.arange(5e-1,
                            np.nanmax(ccld2)+3.9,1.3/8),
             alpha=0.5,density=True)
    plt.hist(np.ravel(hcld2),
             bins=np.arange(5e-1,
                            np.nanmax(hcld2)+3.9,1.3/8),
             alpha=0.5,density=True)

    plt.title(str(i)+' '+cat['name'][i])
    plt.show()


for i in [20]:
    cld2 = fits.open('mom2_hi.fits')[i].data
    plt.contourf(cld2,levels=np.arange(6,9,0.37),cmap='winter')
    plt.colorbar()
    plt.contourf(cld2,levels=np.arange(9,12,0.37),cmap='magma')
    plt.colorbar()


for i in [18]:
    cld2 = fits.open('mom2_hi.fits')[i].data
    plt.contourf(cld2,levels=np.arange(5,7,0.37),cmap='winter')
    plt.colorbar()
    plt.contourf(cld2,levels=np.arange(7,12,0.37),cmap='magma')
    plt.colorbar()


for i in [16]:
    cld2 = fits.open('mom2_hi.fits')[i].data
    plt.contourf(cld2,levels=np.arange(5,7,0.37),cmap='winter')
    plt.colorbar()
    plt.contourf(cld2,levels=np.arange(7,12,0.37),cmap='magma')
    plt.colorbar()


for i in range(22):
    ccld1 = fits.open('mom1_co.fits')[i].data
    hcld1 = fits.open('mom1_hi.fits')[i].data
    plt.hist(np.ravel(ccld1),
             bins=np.arange(np.nanmin(ccld1)-2.6,
                            np.nanmax(ccld1)+3.9,1.3/4),
             alpha=0.5,density=True)
    plt.hist(np.ravel(hcld1),
             bins=np.arange(np.nanmin(hcld1)-2.6,
                            np.nanmax(hcld1)+3.9,1.3/4),
             alpha=0.5,density=True)
    
    bigh1 = h1[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]]
    plt.hist(np.ravel(bigh1),
             bins=np.arange(np.nanmin(bigh1)-2.6,
                            np.nanmax(bigh1)+3.9,1.3/4),
             alpha=0.5,density=True)
    
    bigishh1 = h1[cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    plt.hist(np.ravel(bigishh1),
             bins=np.arange(np.nanmin(bigh1)-2.6,
                            np.nanmax(bigh1)+3.9,1.3/4),
             alpha=0.5,density=True)

    plt.title(str(i)+' '+cat['name'][i])
    plt.show()


for i in [0,2,6,9,11,12,19,20]:
    hi = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    co = dht.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    plt.plot(np.arange(155)*1.3-100,100*np.nanmean(co,axis=(1,2)))
    plt.plot(np.arange(155)*1.3-100,np.nanmean(hi,axis=(1,2)))
    plt.show()


for i in [11]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(20,28,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(28,35,0.65),cmap='magma')
    plt.colorbar()


for i in [10]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-3,6,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(6,11,0.65),cmap='magma')
    plt.colorbar()


for i in [10]:
    cld1 = fits.open('mom1_hi.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-3,5,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(5,9,0.65),cmap='magma')
    plt.colorbar()


for i in [21]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(0,10,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(10,15,0.65),cmap='magma')
    plt.colorbar()


for i in [21]:
    cld1 = fits.open('mom1_hi.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-5,0,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(0,8,0.65),cmap='magma')
    plt.colorbar()


for i in [20]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(5,15,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(15,20,0.65),cmap='magma')
    plt.colorbar()


for i in [18]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-3,5,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(5,12,0.65),cmap='magma')
    plt.colorbar()


for i in [18]:
    cld1 = fits.open('mom1_hi.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-3,8,0.65),cmap='PiYG')
    plt.colorbar()
#     plt.contourf(cld1,levels=np.arange(8,17,0.65),cmap='magma')
#     plt.colorbar()


for i in [13]:
    cld1 = fits.open('mom1_hi.fits')[i].data
    plt.contourf(cld1,levels=np.arange(5,8,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(8,17,0.65),cmap='magma')
    plt.colorbar()


for i in [13]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(5,11,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(11,17,0.65),cmap='magma')
    plt.colorbar()


for i in [6]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(8,13,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(13,17,0.65),cmap='magma')
    plt.colorbar()


for i in [6]:
    cld1 = fits.open('mom1_hi.fits')[i].data
    plt.contourf(cld1,levels=np.arange(8,12,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(12,15,0.65),cmap='magma')
    plt.colorbar()


for i in [1]:
    cld1 = fits.open('mom1_hi.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-15,-5,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(-5,0,0.65),cmap='magma')
    plt.colorbar()


for i in [1]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-10,-3,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(-3,10,0.65),cmap='magma')
    plt.colorbar()


for i in [2]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-5,-1,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(-1,3,0.65),cmap='magma')
    plt.colorbar()


for i in [3]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-2,3,1.3/4),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(3,8,1.3/4),cmap='magma')
    plt.colorbar()


for i in [4]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(10,15,0.65/2),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(15,23,0.65/2),cmap='magma')
    plt.colorbar()


for i in [5]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(4,6,1.3/4),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(6,8,1.3/4),cmap='magma')
    plt.colorbar()


for i in [6]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(6,13,1.3/4),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(13,20,1.3/4),cmap='magma')
    plt.colorbar()


for i in [8]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(2,8,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(8,16,0.65),cmap='magma')
    plt.colorbar()


for i in [11]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(18,28,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(28,36,0.65),cmap='magma')
    plt.colorbar()


for i in [15]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(0,7,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(7,20,0.65),cmap='magma')
    plt.colorbar()


for i in [17]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-10,6,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(6,20,0.65),cmap='magma')
    plt.colorbar()


for i in [18]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-5,3,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(3,12,0.65),cmap='magma')
    plt.colorbar()


for i in [19]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(-10,-1,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(-1,8,0.65),cmap='magma')
    plt.colorbar()


for i in [20]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(5,15,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(15,25,0.65),cmap='magma')
    plt.colorbar()


for i in [21]:
    cld1 = fits.open('mom1_co.fits')[i].data
    plt.contourf(cld1,levels=np.arange(0,10,0.65),cmap='winter')
    plt.colorbar()
    plt.contourf(cld1,levels=np.arange(10,20,0.65),cmap='magma')
    plt.colorbar()


# when we see a density plot of a cloud, we are seeing all of the material along that line of sight collapsed together. when we use density thresholds to define the boundary of a cloud, we are assuming that within that area of the sky, all of the material that is contributing to the projected density is associated with one another together, so that it is physically consistent to group that meterial into an individual structure. this assumption describes the actualy physical set up of the galaxy to varying degrees of accuracy. in the galactic plane, where there is thousands of kpc that any line of sight is crossing through, it is harder to confidently isolate structures that are tens of parsecs from the rest of the gas dust etc. on the other hand at higher latitudes, where the line of sight is crossing through sparser distributions of gas, it becomes easier to say yes all of this stuff that we are seeing near each other belongs together. the velocity spectrum along a given line of sight can also give us guidance on whether the emission being observed is due to multiple, separate objects that are spread out in radial distance, or whether it is one structure. if we assume that each cloud has a certain bulk velocity, and that tturbulent motions give rise to a randomized distribution of velocities to small portions within that cloud, then we would expect the spectrum to be close to a normal distribution / gaussian. on the other hand, if there are multiple distinct peaks of emission with the spectrum, we might suppose that what we are seeing is actually two different clouds, each with their own bulk+spread of velocities superimposed in projection. however, this only works when there is a clear mapping between bulk radial velocity and radial distance. this is usually not the case, and even when it is, the differences that significant leaps in distance cause in los velocity are so small that the turbulent linewidths wash out the distinctions.


lod,lid = ((cat['lon_GMC<deg>']/-0.229)+786).astype('int').transpose()
bod,bid = ((cat['lat_GMC<deg>']/0.229)+393).astype('int').transpose()


lod,lid,bod,bid


lod,lid,bod,bid


dust_std = fits.open('mean_and_std_lbd.fits')[2].data


fits.open('mean_and_std_lbd.fits')[2].header


[(
plt.contourf(
    np.nansum(
        dust_std[max(z19d-100-68,0):z19d+100-68,
           bod[i]-10:bid[i]+10,
           lod[i]-10:lid[i]+10]
        ,axis=0),
        extent=[cat['lon_GMC<deg>'][i][0]+10*0.229,
              cat['lon_GMC<deg>'][i][1]-10*0.229,
              cat['lat_GMC<deg>'][i][0]-10*0.229,
              cat['lat_GMC<deg>'][i][1]+10*0.229]),
  
plt.colorbar(),
plt.contour(
    np.log10(
          fits.open('mom0_co.fits')[i].data),
      levels=np.arange(20,23,0.25),
      extent=[cat['lon_GMC<deg>'][i][0],
              cat['lon_GMC<deg>'][i][1],
              cat['lat_GMC<deg>'][i][0],
              cat['lat_GMC<deg>'][i][1]],
      colors='white',linewidths=0.6),
  
plt.show())
      
      for i in range(22)
      for z19d in [int(cat['distance<pc>'][i])]
      ]
    

from matplotlib.patches import Circle


c0s = fits.open('mom0_co.fits')


plt.scatter(cat['name'],cat['vel_grad_GMC<km/s/pc>'])
plt.xticks(rotation='vertical')
plt.show()


cat


fig,axs = plt.subplots(2,3,figsize=(30,12))


mom0 = axs[0,0].contourf(np.log10(c0s[1].data),
                      extent=[166.5,154,-13,-6],
                      levels=np.arange(20.5,22.5,0.25),
                      cmap='magma')


mom1 = axs[0,1].contourf(c1s[1].data,
                       extent=[166.5,154,-13,-6],
                       cmap='PiYG',levels=np.arange(-8,3,0.65),
                      extend='both')


dist = axs[0,2].contourf(new[1],
                       extent=[166.5,154,-13,-6],
                       levels=np.arange(450,505,5),extend='both',cmap='binary_r')

axs[0,2].contour(c1s[1].data,
               extent=[166.5,154,-13,-6],
              levels=np.arange(-8,3,0.65),cmap='PiYG')


dvdy,dvdx = np.gradient(c1s[1].data)

for ax in axs.flat:
    ax.set_xlim(166.5,154)

    
axs[1,0].contourf(np.nansum(dust[400-68:470-68,bod[1]:bid[1],lod[1]:lid[1]],axis=0),
                 extent=[166.5,154,-13,-6],levels=np.arange(1e-2,7e-1,5e-2))

axs[1,1].contourf(np.nansum(dust[470-68:520-68,bod[1]:bid[1],lod[1]:lid[1]],axis=0),
                 extent=[166.5,154,-13,-6],levels=np.arange(1e-2,7e-1,3e-2))

axs[1,2].contourf(np.nansum(dust[:,bod[1]:bid[1],lod[1]:lid[1]],axis=1),
                 extent=[166.5,154,68,1250],levels=np.arange(1e-2,3e-1,2e-2))

axs[1,2].set_ylim(400,550)

fig.colorbar(mom0,ax=axs[0,0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[0,1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[0,2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


cat


fig,axs = plt.subplots(1,3,figsize=(30,12))


mom0 = axs[0].contourf(np.log10(c0s[8].data),
                      extent=[49,41,7,11],
                      levels=np.arange(20.5,22,0.15),
                      cmap='magma')


mom1 = axs[1].contourf(c1s[8].data,
                       extent=[49,41,7,11],
                       cmap='PiYG',levels=np.arange(5,15,0.65),
                      extend='both')


dist = axs[2].contourf(new[8],
                       extent=[49,41,7,11],
                       levels=np.arange(228,237,0.5),extend='both',cmap='binary_r')

axs[2].contour(c1s[8].data,
               extent=[49,41,7,11],
              levels=np.arange(5,15,0.65),cmap='PiYG')


for ax in axs.flat:
    ax.set_xlim(49,41)


fig.colorbar(mom0,ax=axs[0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


cat[13]


fig,axs = plt.subplots(1,3,figsize=(30,12))


mom0 = axs[0].contourf(np.log10(c0s[13].data),
                      extent=[-141.5,-149.5,-14.5,-10.5],
                      levels=np.arange(20.5,22,0.15),
                      cmap='magma')


mom1 = axs[1].contourf(c1s[13].data,
                       extent=[-141.5,-149.5,-14.5,-10.5],
                       cmap='PiYG',levels=np.arange(8,14,0.65),
                      extend='both')


dist = axs[2].contourf(new[13],
                       extent=[-141.5,-149.5,-14.5,-10.5],
                       levels=np.arange(740,850,10),extend='both',cmap='binary_r')

axs[2].contour(c1s[13].data,
               extent=[-141.5,-149.5,-14.5,-10.5],
              levels=np.arange(8,14,0.65),cmap='PiYG')


for ax in axs.flat:
    ax.set_xlim(-141.5,-149.5)


fig.colorbar(mom0,ax=axs[0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


cat[10]


fig,axs = plt.subplots(1,3,figsize=(30,12))


mom0 = axs[0].contourf(np.log10(c0s[10].data),
                      extent=[-8,-27,3.5,25.5],
                      levels=np.arange(20.5,22.5,0.15),
                      cmap='magma')


mom1 = axs[1].contourf(c1s[10].data,
                       extent=[-8,-27,3.5,25.5],
                       cmap='PiYG',levels=np.arange(2,10,0.65),
                      extend='both')


dist = axs[2].contourf(new[10],
                       extent=[-8,-27,3.5,25.5],
                       levels=np.arange(150,245,5),extend='both',cmap='binary_r')

axs[2].contour(c1s[10].data,
               extent=[-8,-27,3.5,25.5],
              levels=np.arange(2,10,0.65),cmap='PiYG')


for ax in axs.flat:
    ax.set_xlim(-8,-27)


fig.colorbar(mom0,ax=axs[0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


fig,axs = plt.subplots(1,3,figsize=(30,12))


mom0 = axs[0].contourf(np.log10(c0s[14].data),
                      extent=[11.5,-11.5,12.5,25],
                      levels=np.arange(20.5,22.5,0.15),
                      cmap='magma')


mom1 = axs[1].contourf(c1s[14].data,
                       extent=[11.5,-11.5,12.5,25],
                       cmap='PiYG',levels=np.arange(-1,5,0.65),
                      extend='both')


dist = axs[2].contourf(new[14],
                       extent=[11.5,-11.5,12.5,25],
                       levels=np.arange(130,175,3),extend='both',cmap='binary_r')

axs[2].contour(c1s[14].data,
               extent=[11.5,-11.5,12.5,25],
              levels=np.arange(-1,5,0.65),cmap='PiYG')


for ax in axs.flat:
    ax.set_xlim(11.5,-11.5)


fig.colorbar(mom0,ax=axs[0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


fig,axs = plt.subplots(2,3,figsize=(35,20))


mom0 = axs[0,0].contourf(np.log10(c0s[15].data),
                      extent=[-143,-155,-21,-17,],
                      levels=np.arange(20.5,22.5,0.15),
                      cmap='magma')


mom1 = axs[0,1].contourf(c1s[15].data,
                       extent=[-143,-155,-21,-17],
                       cmap='PiYG',levels=np.arange(4,12,0.65),
                      extend='both')


dist = axs[0,2].contourf(new[15],
                       extent=[-143,-155,-21,-17],
                       levels=np.arange(425,490,3),extend='both',cmap='binary_r')

axs[0,2].contour(c1s[15].data,
               extent=[-143,-155,-21,-17,],
              levels=np.arange(4,12,0.65),cmap='PiYG')


axs[1,0].contourf(np.nansum(dust[400-68:450-68,bod[15]:bid[15],lod[15]:lid[15]],axis=0),
                 extent=[-143,-155,-21,-17],levels=np.arange(1e-2,1.2,5e-2))

axs[1,1].contourf(np.nansum(dust[450-68:520-68,bod[15]:bid[15],lod[15]:lid[15]],axis=0),
                 extent=[-143,-155,-21,-17],levels=np.arange(1e-2,5e-1,3e-2))

axs[1,2].contourf(np.nansum(dust[:,bod[15]:bid[15],lod[15]:lid[15]],axis=1),
                 extent=[-143,-155,68,1250],levels=np.arange(2e-2,5e-1,3e-2))

axs[1,2].set_ylim(400,550)

for ax in axs.flat:
    ax.set_xlim(-143,-155)


fig.colorbar(mom0,ax=axs[0,0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[0,1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[0,2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


fig,axs = plt.subplots(2,3,figsize=(30,25))


mom0 = axs[0,0].contourf(np.log10(c0s[16].data),
                      extent=[-151,-157.5,-17,-10],
                      levels=np.arange(20.5,22.5,0.15),
                      cmap='magma')


mom1 = axs[0,1].contourf(c1s[16].data,
                       extent=[-151,-157.5,-17,-10],
                       cmap='PiYG',levels=np.arange(5,12,0.65),
                      extend='both')


dist = axs[0,2].contourf(new[16],
                       extent=[-151,-157.5,-17,-10],
                       levels=np.arange(415,465,3),extend='both',cmap='binary_r')

axs[0,2].contour(c1s[16].data,
               extent=[-151,-157.5,-17,-10],
              levels=np.arange(5,12,0.65),cmap='PiYG')

axs[1,0].contourf(np.nansum(dust[390-68:430-68,bod[16]:bid[16],lod[16]:lid[16]],axis=0),
                 extent=[-151,-157.5,-17,-10],levels=np.arange(1e-2,1.2,5e-2))

axs[1,1].contourf(np.nansum(dust[430-68:500-68,bod[16]:bid[16],lod[16]:lid[16]],axis=0),
                 extent=[-151,-157.5,-17,-10],levels=np.arange(1e-2,5e-1,3e-2))

axs[1,2].contourf(np.nansum(dust[:,bod[16]:bid[16],lod[15]:lid[15]],axis=2),
                 extent=[-151,-157.5,68,1250],levels=np.arange(1e-2,5e-1,3e-2))


axs[1,2].set_ylim(375,530)


for ax in axs.flat:
    ax.set_xlim(-151,-157.5)


fig.colorbar(mom0,ax=axs[0,0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[0,1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[0,2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


fig,axs = plt.subplots(1,3,figsize=(20,15))


mom0 = axs[0].contourf(np.log10(c0s[18].data),
                      extent=[162,155,-26,-15.5],
                      levels=np.arange(20,22.5,0.25),
                      cmap='magma')


mom1 = axs[1].contourf(c1s[18].data,
                       extent=[162,155,-26,-15.5],
                       cmap='PiYG',levels=np.arange(0,10,0.65),
                      extend='both')


dist = axs[2].contourf(new[18],
                       extent=[162,155,-26,-15.5],
                       levels=np.arange(280,312,2),extend='both',cmap='binary_r')

axs[2].contour(c1s[18].data,
               extent=[162,155,-26,-15.5],
              levels=np.arange(-2,10,0.65),cmap='PiYG')


or ax in axs.flat:
    ax.set_xlim(162,155)


fig.colorbar(mom0,ax=axs[0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


fig,axs = plt.subplots(1,3,figsize=(20,15))


mom0 = axs[0].contourf(np.log10(c0s[19].data),
                      extent=[127,118,20,31],
                      levels=np.arange(20.5,22,0.15),
                      cmap='magma')


mom1 = axs[1].contourf(c1s[19].data,
                       extent=[127,118,20,31],
                       cmap='PiYG',levels=np.arange(-5,1,0.65),
                      extend='both')


dist = axs[2].contourf(new[19],
                       extent=[127,118,20,31],
                       levels=np.arange(290,390,5),extend='both',cmap='binary_r')

axs[2].contour(c1s[19].data,
               extent=[127,118,20,31],
              levels=np.arange(-5,1,0.65),cmap='PiYG')


for ax in axs.flat:
    ax.set_xlim(127,118)
#     ax.add_patch(Circle((159.5,-19.5),0.5,fill=False,edgecolor='red'))


fig.colorbar(mom0,ax=axs[0],label=r'$n_{h_2}\ [cm^{-2}]$',orientation='horizontal')
fig.colorbar(mom1,ax=axs[1],label=r'$v_{lsr}\ [km\ s^{-1}]$',orientation='horizontal')
fig.colorbar(dist,ax=axs[2],label=r'$d\ [pc]$',orientation='horizontal')
# fig.colorbar(grad,ax=axs[2],label=r'$dv/ds$',orientation='horizontal')

plt.show()


for i in range(22):
    if i not in [4,7,11,20]:
        plt.contourf(np.nansum(2.8*dust[:,
                                    bod[i]:bid[i],
                                    lod[i]:lid[i]],axis=2),
                     levels=np.arange(1e-2,1.7,2e-2),extend='both',
                    extent=[cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
                            68,1250])
        
        plt.ylim(cat['distance<pc>'][i]-100,cat['distance<pc>'][i]+100)
        plt.scatter(np.nanmean(cat['lat_GMC<deg>'][i]),
                    cat['distance<pc>'][i],marker='x',
                   color='white',s=100)
        plt.title(cat['name'][i])
        plt.show()


for i in range(22):
    if i not in [4,7,11,20]:
        plt.contourf(np.nansum(2.8*dust[:,
                                    bod[i]:bid[i],
                                    lod[i]:lid[i]],axis=1),
                     levels=np.arange(1e-2,1.7,2e-2),extend='both',
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                            68,1250])
        
        plt.ylim(cat['distance<pc>'][i]-100,cat['distance<pc>'][i]+100)
        plt.scatter(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i],marker='x',
                   color='white',s=100)
        plt.title(cat['name'][i])
        plt.show()


pd.read_csv('z')


cat['distance<pc>'][19]


plt.imshow(new[19],origin='lower')


help(plt.quiver)


cat


c0s = fits.open('mom0_co.fits')


[(
plt.contourf(
    np.nansum(
        dust[max(z19d-100-68,0):z19d+100-68,
           max(bod[i]-70,0):bid[i]+70,
           max(lod[i]-70,0):lid[i]+70]
        ,axis=0),
        extent=[cat['lon_GMC<deg>'][i][0]+70*0.23,
              cat['lon_GMC<deg>'][i][1]-70*0.23,
              cat['lat_GMC<deg>'][i][0]-70*0.23,
              cat['lat_GMC<deg>'][i][1]+70*0.23]),

    plt.colorbar(),
    
    plt.contour(
    np.nansum(
        dust[max(z19d-100-68,0):z19d+100-68,
           max(bod[i]-70,0):bid[i]+70,
           max(lod[i]-70,0):lid[i]+70]
        ,axis=0),
        extent=[cat['lon_GMC<deg>'][i][0]+70*0.23,
              cat['lon_GMC<deg>'][i][1]-70*0.23,
              cat['lat_GMC<deg>'][i][0]-70*0.23,
              cat['lat_GMC<deg>'][i][1]+70*0.23],
    levels=[0.1,0.15,0.2,0.5],
    linewidths=0.6,colors=['red','white','white','red']),
  

plt.contour(
    np.log10(
          fits.open('mom0_co.fits')[i].data),
      levels=np.arange(20,23,0.5),
      extent=[cat['lon_GMC<deg>'][i][0],
              cat['lon_GMC<deg>'][i][1],
              cat['lat_GMC<deg>'][i][0],
              cat['lat_GMC<deg>'][i][1]],
      colors='white',linewidths=0.6),
    
plt.title(cat['name'][i]),
plt.show())
      
      for i in range(22)
      for z19d in [int(cat['distance<pc>'][i])]
      ]
    

fits.open('mean_and_std_lbd.fits')[1].header


def get_mom1(cube,vels,dv):
    xTx = np.array([cube[i]*vels[i]*dv for i in range(len(vels))])
    sigxTx = np.nansum(xTx,axis=0)

    M0 = np.nansum(cube*dv,axis=0)
    return sigxTx/M0


[(plt.contourf(get_mom1(dust[max(z19d-100-68,0):z19d+100-68,
           bod[i]:bid[i],
           lod[i]:lid[i]],np.arange(max(z19d-100-68,0),z19d+100-68)+68,1),levels=15),
  plt.colorbar(),
  plt.title(cat['name'][i]),
 plt.show())
        
      for i in range(22)
      for z19d in [int(cat['distance<pc>'][i])]
      if i not in [4,7,11,20]
      ]
    

dust_mom1s = [get_mom1(dust[max(z19d-100-68,0):z19d+100-68,
           bod[i]:bid[i],
           lod[i]:lid[i]],np.arange(max(z19d-100-68,0),z19d+100-68)+68,1)
        
      for i in range(22)
      for z19d in [int(cat['distance<pc>'][i])]
      if i not in [4,7,11,20]
      ]
    

new = dust_mom1s[0:4]+[np.nan]+dust_mom1s[4:6]+\
[np.nan]+dust_mom1s[6:9]+[np.nan]+\
dust_mom1s[9:17]+[np.nan]+[dust_mom1s[17]]


c1s = fits.open('mom1_co.fits')


from scipy.ndimage.interpolation import zoom


0.125*u.deg.to(u.rad)*800


cat[cat['distance<pc>']>800]


[
#     (plt.scatter(np.ravel(c1s[i].data),np.ravel(new[i])),
#  plt.show())
    print(c1s[i].data[::2,::2].shape,new[i].data.shape)
    
    for i in range(22)
    if i not in [4,7,11,20]
]


for i in range(22):
    if i not in [4,7,11,20]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[:,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,5e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.2,0.5,1],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        plt.show()


for i in range(22):
    if i not in [4,7,11,20]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[:,bod:bid,lod:lid]

        ld = np.nansum(sb*2.8,axis=1)

        plt.contourf(ld,levels = np.arange(1e-2,1,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           68,1250])
        
        plt.colorbar()
        
        
        plt.contour(ld,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           68,1250])
        
        
        plt.title(cat['name'][i])
        
        plt.ylim(cat['distance<pc>'][i]-100,cat['distance<pc>'][i]+100)
        plt.show()


for i in range(22):
    if i not in [4,7,11,20]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[:,bod:bid,lod:lid]

        bd = np.nansum(sb*2.8,axis=2)

        plt.contourf(bd,levels = np.arange(1e-2,1,3e-2),
                    extent=[cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
                           68,1250])
        
        plt.colorbar()
        
        
        plt.contour(bd,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
                           68,1250])
        
        
        plt.title(cat['name'][i])
        
        plt.ylim(cat['distance<pc>'][i]-100,cat['distance<pc>'][i]+100)
        plt.show()


for i in [1]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[420-68:475-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[475-68:550-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
for i in [5]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[100-68:180-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[200-68:230-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
for i in [10]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[100-68:125-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[125-68:225-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[225-68:300-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()


for i in [14]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[175-68:,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[75-68:130-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[130-68:175-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
#         plt.ylim(cat['distance<pc>'][i]-100,cat['distance<pc>'][i]+100)
        plt.show()


for i in [15]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[400-68:450-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[450-68:500-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[500-68:525-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
for i in [16]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[375-68:440-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[440-68:500-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
for i in [19]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[350-68:375-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[375-68:400-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
for i in [6,13]:
    lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
    bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

    sb = dust[:,bod:bid,lod:lid]

    bd = np.nansum(sb*2.8,axis=2)

    plt.contourf(bd,levels = np.arange(1e-2,1,3e-2),
                extent=[cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
                       68,1250])

    plt.colorbar()


    plt.contour(bd,levels = [0.1,0.2],colors='white',linewidths=0.8,
                extent=[cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
                       68,1250])


    plt.title(cat['name'][i])

    plt.ylim(cat['distance<pc>'][i]-400,cat['distance<pc>'][i]+400)
    plt.show()


for i in [6,13]:
    lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
    bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

    sb = dust[:,bod:bid,lod:lid]

    ld = np.nansum(sb*2.8,axis=1)

    plt.contourf(ld,levels = np.arange(1e-2,1,3e-2),
                extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       68,1250])

    plt.colorbar()


    plt.contour(ld,levels = [0.1,0.2],colors='white',linewidths=0.8,
                extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       68,1250])


    plt.title(cat['name'][i])

    plt.ylim(cat['distance<pc>'][i]-700,cat['distance<pc>'][i]+700)
    plt.show()


for i in [6]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[500-68:800-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[800-68:,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
for i in [13]:
        lod,lid = ((cat['lon_GMC<deg>'][i]-180)/-0.23).astype('int')
        bod,bid = ((cat['lat_GMC<deg>'][i]+90)/0.23).astype('int')

        sb = dust[400-68:600-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[600-68:800-68,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
        sb = dust[800-68:,bod:bid,lod:lid]

        lb = np.nansum(sb*2.8,axis=0)

        plt.contourf(lb,levels = np.arange(1e-2,2,3e-2),
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                          cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        plt.colorbar()
        
        
        plt.contour(lb,levels = [0.1,0.2],colors='white',linewidths=0.8,
                    extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                           cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
        
        
        plt.title(cat['name'][i])
        
        plt.show()
        
        
plt.hist(np.log2(cat['vel_grad_GMC<km/s/pc>']/cat['vel_grad_HI<km/s/pc>']),bins=np.arange(-3,4,0.25))

plt.show()
