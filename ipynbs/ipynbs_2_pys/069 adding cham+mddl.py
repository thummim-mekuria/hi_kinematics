#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import datacube as dc
import clouds as cld
import interpolatepolate as ip

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


zcat = QTable(fits.open('final_cat.fits')[1].data)


zcat #might also need to remove both or choose one for cepheus


dht = fits.open('DHT_CO.fits')[0]
dht_r = fits.open('COGAL_all_mom.fits')[0]
ipd = ip.interpolate_vlb(dht_r.data)


ipd_rshp = ipd.transpose().swapaxes(1,2)[169:324]


hdu = fits.PrimaryHDU(ipd_rshp)
hdu.writeto('DHT_CO_longer.fits',overwrite=True)


def gft(spc):
    from astropy.modeling import models,fitting
    fitter = fitting.LMLSQFitter()
    
    ginit = models.Gaussian1D(np.nanmax(spc),np.nanargmax(spc),10*np.nanstd(spc))

    fit = fitter(ginit,np.arange(len(spc)),spc)
    return fit


for i in range(21):
    lo,li = ((zcat['lon<deg>'][i] - 180.125)/(-0.125)).astype('int')
    bo,bi = ((zcat['lat<deg>'][i] + 30.125)/(0.125)).astype('int')
    
    sb = ipd_rshp[:,bo:bi,lo:li]
    sig,ar = cld.surface_density(sb,'CO')
    
#     ang_ar = (u.rad**2).to(u.deg**2) *ar*len(np.where(sig>0)[0])
#     rgad = np.sqrt(ang_ar/np.pi)
#     pc2_ar = ar*len(np.where(sig>0)[0])*(zcat['distance<pc>'][i]*u.pc)**2
#     rg = np.sqrt(pc2_ar/np.pi)
#     zcat['radius_GMC<deg>'][i] = rgad
#     zcat['radius_GMC<pc>'][i] = rg.value
    
#     m = cld.mass(sig,zcat['distance<pc>'][i]*u.pc,ar)
#     ra,s = cld.accumulation_radius(m,zcat['distance<pc>'][i]*u.pc,u.deg*np.nanmean(zcat['lat<deg>'][i]))
    
#     zcat['mass_GMC<kg>'][i] = m.value
#     zcat['radius_HI<pc>'][i] = ra.value
#     zcat['radius_HI<deg>'][i] = s.value
    
    
    spc = np.nanmean(sb,axis=(1,2))
    
    g = gft(spc)
    zcat['vdispersion_GMC<km/s>'][i] = 1.3*g.stddev.value
    
#     plt.plot(spc)
#     plt.plot(g(np.arange(155)))
#     plt.show()
    
#     zcat['vspectrum_GMC<K>'][i] = spc
    
    
#     w = np.nansum(sb,axis=0)
    plt.imshow(sig.value,origin='lower',cmap=dc.cmap_min0('viridis'))
    plt.show()


zcat.write('final_cat.fits',overwrite=True)


zcat[['lon<deg>','lat<deg>']][[5,11]]


Mdl,[-138 -147],[-5  3],2072


ds = pd.read_csv('whole_gal_dists.csv')


180+(2001)*(-0.125)


(-28.875+30.125)/0.125


178-180


180.125+10*(-0.125)


zcat = zcat['names',
 'lon<deg>',
 'lat<deg>',
 'distance<pc>',
 'mass_GMC<kg>',
 'radius_GMC<pc>',
 'radius_GMC<deg>',
 'radius_HI<pc>',
 'radius_HI<deg>',
 'vdispersion_GMC<km/s>',
 'vspectrum_GMC<K>']


hi4pi = fits.open('HI4PI_flatcalib_80.fits')[0]


zls = np.nanmean(zcat['lon<deg>'],axis=1)
zbs = np.nanmean(zcat['lat<deg>'],axis=1)


whi = np.nansum(hi4pi.data,axis=0)


fig,ax = plt.subplots(1,1,figsize=(20,10))


ax.imshow(whi,origin='lower',\
           cmap=dc.cmap_min0('binary_r'),extent=[180,-180,-50,50],vmax=3e3)
for i in range(19):
    lc,bc = zls[i],zbs[i]
    rhi = zcat['radius_HI<deg>'][i]
    
    ax.add_patch(Circle((lc,bc),rhi,fill=False,edgecolor='red'))
                 

plt.xlim(180,-5)        
plt.show()
    
    
winters = 


zcat


fig,ax = plt.subplots(1,1,figsize=(20,10))


ax.imshow(whi,origin='lower',\
           cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=4e3,vmin=2e3)
for i in range(19):
    lc,bc = zls[i],zbs[i]
    if lc < -100:
        rhi = zcat['radius_HI<deg>'][i]

        ax.add_patch(Circle((lc,bc),rhi,fill=False,linewidth=3,label=zcat['names'][i]))
    
plt.legend()
  

plt.xlim(-125,-180)
plt.ylim(-35,15)
plt.show()
    
    
from astrodendro import Dendrogram


print((np.array((-125,-180)) -180.166)/(-0.08333))

print((np.array((-35,15)) +49.999)/(0.08333))


np.nanmin(hi4pi.data)


np.nanmean(hi4pi.data),np.nanstd(hi4pi.data)


fig,ax = plt.subplots(1,1,figsize=(20,10))


ax.imshow(whi,origin='lower',\
           cmap=dc.cmap_min0('viridis'),extent=[180,-180,-50,50],vmax=4e3,vmin=2e3)
for i in range(19):
    lc,bc = zls[i],zbs[i]
    rhi = zcat['radius_HI<deg>'][i]
    
    ax.add_patch(Circle((lc,bc),rhi,fill=False,edgecolor='black',linewidth=3))
    
  
# plt.xlim(-125,-180)
# plt.ylim(-35,15)
plt.show()
