#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.modeling import models,fitting

import datacube as dc
from scipy import ndimage


hi4pi = fits.open('HI4PI.fits')[0]


hdu = fits.PrimaryHDU(smooth_043_80_wrap)


hdu.writeto('smoothed_HI4PI_sig80_rms043.fits')


hdu1 = fits.PrimaryHDU(hi4pi_calib_80)


hdu1.writeto('HI4PI_flatcalib_80.fits')


hdu2 = fits.PrimaryHDU(smooth_043_40_wrap)
hdu2.writeto('smoothed_HI4PI_sig40_rms043.fits',overwrite=True)


hdu3 = fits.PrimaryHDU(hi4pi_calib)
hdu2.writeto('HI4PI_flatcalib_40.fits',overwrite=True)


flat = smooth_043_40_wrap/np.nanmax(smooth_043_40_wrap)


hi4pi_calib = (hi4pi.data*sigg)/flat


flat80 = smooth_043_80_wrap/np.nanmax(smooth_043_80_wrap)
hi4pi_calib_80 = (hi4pi.data*sigg)/flat80


wcalib = np.nansum(hi4pi_calib,axis=0)


plt.imshow(wcalib,origin='lower',cmap=dc.cmap_min0('viridis'),vmax=8e3)
plt.colorbar(orientation='horizontal')


wcalib80 = np.nansum(hi4pi_calib_80,axis=0)


ws80 = np.nansum(smooth_043_80_wrap,axis=0)


plt.imshow(ws,origin='lower',cmap=dc.cmap_min0('viridis'),vmax=3e3)


plt.imshow(ws80,origin='lower',cmap=dc.cmap_min0('viridis'),vmax=3e3)


plt.imshow(np.nansum(flat80,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))


ws = np.nansum(smooth_043_40_wrap,axis=0)


plt.imshow(np.nansum(flat,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))


plt.imshow(wcalib80,origin='lower',cmap=dc.cmap_min0('viridis'),vmax=6e3)
plt.colorbar(orientation='horizontal')


def get_g_fit(spec):
    fitter = fitting.LevMarLSQFitter()
#     ginit = models.Gaussian1D(1.,0.,1.)
    ginit = models.Gaussian1D(np.nanmax(spec),np.argmax(spec),np.nanstd(spec))
    gfit = fitter(ginit,np.arange(len(spec)),spec)
    
    return gfit


zcat = QTable(fits.open('final_cat.fits')[1].data)


for i in range(19):
    lo,li = ((zcat['lon<deg>'][i] - 180.166)/(-0.0833)).astype('int')
    bo,bi = ((zcat['lat<deg>'][i] + 49.999)/(0.08333)).astype('int')
    
    sb = hi4pi.data[:,bo:bi,lo:li]
    spc = np.nanmean(sb,axis=(1,2))
    
    vels = np.arange(310)*1.288 - 199.572
    velsco = np.arange(155)*1.3 - 100.1
    
    plt.plot(vels,spc)
    
    gfit = get_g_fit(spc)
    plt.plot(vels,gfit(np.arange(310)))
    
    fac = np.nanmax(spc)/np.nanmax(zcat['vspectrum_GMC<K>'][i])
    plt.plot(velsco,fac*zcat['vspectrum_GMC<K>'][i])
    plt.title(zcat['names'][i])
    plt.xlim(-80,80)
    
    
    plt.show()


#co<->hi mismatch in width: cal,cma,crb,her,mob,ros


zls = np.nanmean(zcat['lon<deg>'],axis=1)


zbs = np.nanmean(zcat['lat<deg>'],axis=1)


#these got cut-off for: cal,gem,orl,pers,tau _> being too close to the edge -- need to tell everything to wrap
#and also go back to look at when/where i made the flat calib & to remind myself whats going on


for i in [18]:
    lc = (zls[i] - 180.166)/(-0.0833)
    bc = (zbs[i] + 49.999)/(0.08333)
    
    rhi = zcat['radius_HI<deg>'][i]/0.0833
    
    #the assumptions in this series of conditionals
    # 1 - that the only indies that can meaningfully go beyond an edge are ones along lon
    # 2 - that something can only go off one lon edge or another, but not both
    # 3 - the only possible states of being are: go off left edge, go off right edge, or stay inside
    l0,l1 = int(lc-rhi),int(lc+rhi)
    b0,b1 = int(bc-rhi),int(bc+rhi)
    
    if l0 < 0:
        sb  = np.concatenate((  hi4pi.data[:,b0:b1,l0:],  hi4pi.data[:,b0:b1,:l1]),axis=2)
        sb2 = np.concatenate((hi4pi_c.data[:,b0:b1,l0:],hi4pi_c.data[:,b0:b1,:l1]),axis=2)
        sb3 = np.concatenate(( hi4pi_calib[:,b0:b1,l0:], hi4pi_calib[:,b0:b1,:l1]),axis=2)


    elif l1 > hi4pi.data.shape[2]:
        l1_n = hi4pi.data.shape[2] - l1
        sb  = np.concatenate((  hi4pi.data[:,b0:b1,l0:],  hi4pi.data[:,b0:b1,:l1_n]),axis=2)
        sb2 = np.concatenate((hi4pi_c.data[:,b0:b1,l0:],hi4pi_c.data[:,b0:b1,:l1_n]),axis=2)
        sb3 = np.concatenate(( hi4pi_calib[:,b0:b1,l0:], hi4pi_calib[:,b0:b1,:l1_n]),axis=2)


    else:
        sb  =   hi4pi.data[:,b0:b1,l0:l1]
        sb2 = hi4pi_c.data[:,b0:b1,l0:l1]
        sb3 =  hi4pi_calib[:,b0:b1,l0:l1]


    plt.imshow(np.nansum(sb,axis=0),origin='lower',\
               cmap=dc.cmap_min0('viridis'))#,vmax = 0.7*np.nanmax(np.nansum(sb,axis=0)))
    plt.colorbar(orientation='vertical')
    plt.show()
    
    plt.imshow(np.nansum(sb2,axis=0),origin='lower',\
               cmap=dc.cmap_min0('viridis'))#,vmax = 0.7*np.nanmax(np.nansum(sb2,axis=0)))
    plt.colorbar(orientation='vertical')
    plt.title(zcat['names'][i])
    plt.show()
    
    plt.imshow(np.nansum(sb3,axis=0),origin='lower',\
               cmap=dc.cmap_min0('viridis'))#,vmax = 0.7*np.nanmax(np.nansum(sb2,axis=0)))
    plt.colorbar(orientation='vertical')
    plt.title(zcat['names'][i])
    plt.show()
    

    print('\n\n\n')


for i in range(19):
    lc = (zls[i] - 180.166)/(-0.0833)
    bc = (zbs[i] + 49.999)/(0.08333)
    
    rhi = zcat['radius_HI<deg>'][i]/0.0833
    
    l0,l1 = int(lc-rhi),int(lc+rhi)
    b0,b1 = int(bc-rhi),int(bc+rhi)
    
    
    sb  =   hi4pi.data[:,b0:b1,l0:l1]
    sb2 = hi4pi_c.data[:,b0:b1,l0:l1]
    sb3 =  hi4pi_calib[:,b0:b1,l0:l1]


    plt.imshow(np.nansum(sb,axis=0),origin='lower',\
               cmap=dc.cmap_min0('viridis'))#,vmax = 0.7*np.nanmax(np.nansum(sb,axis=0)))
    plt.colorbar(orientation='vertical')
    plt.show()
    
    plt.imshow(np.nansum(sb2,axis=0),origin='lower',\
               cmap=dc.cmap_min0('viridis'))#,vmax = 0.7*np.nanmax(np.nansum(sb2,axis=0)))
    plt.colorbar(orientation='vertical')
    plt.title(zcat['names'][i])
    plt.show()
    
    plt.imshow(np.nansum(sb3,axis=0),origin='lower',\
               cmap=dc.cmap_min0('viridis'))#,vmax = 0.7*np.nanmax(np.nansum(sb2,axis=0)))
    plt.colorbar(orientation='vertical')
    plt.title(zcat['names'][i])
    plt.show()
    

    print('\n\n\n')


7/0.083


40*0.083


plt.hist(zcat['radius_HI<deg>'],bins=7)
