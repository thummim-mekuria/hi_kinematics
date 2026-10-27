#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Ellipse

import datacube as dc

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

from photutils.aperture import CircularAnnulus

#------------------------------------------------------------
zcat = QTable(fits.open('final_cat.fits')[1].data)

hi4pi,dht = fits.open('both_same_grid_720d.fits')


#------------------------------------------------------------
hi_conv = (1.3*u.km/u.s *1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p *u.K).to(u.solMass/u.pc**2) 

co_conv = (1.3*u.km/u.s *1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p *u.K).to(u.solMass/u.pc**2) 


#------------------------------------------------------------
lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50


slices_30 = np.array([(-100, -70),(-90, -60),(-80, -50),(-70, -40),(-60, -30),(-50, -20),\
(-40, -10),(-30, 0),(-20, 10),(-10, 20),(0, 30),(10, 40),(20, 50),\
 (30, 60),(40, 70),(50, 80),(60, 90),(70,100)])


slices_30i = np.array((slices_30+100)/1.3).astype('int')


#for defining the gmc-less hi, define the peak as where Sigma = 10, but then it is okay if the edges of that
#are in less dense regions / overlap with other hi
#N_control can be less, it is still informatiive to look at 1 or 2 of those clouds


for i,ii in enumerate(slices_30i):
    hilb = np.nansum(hi4pi.data[ii[0]:ii[1],:,:],axis=0)
    colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    
    fig,ax = plt.subplots(1,1,figsize=(17,8))

    ax.imshow(hilb*hi_conv.value,
               origin='lower',norm='log',extent=[360,-360,-50,50])
    ax.contour(lons,lats,hilb*hi_conv.value,
                levels=[2,5,10,15,20],cmap='RdYlBu',zorder=0)
    ax.contourf(lons,lats,colb*co_conv.value,
                 levels=[10,1e5],colors='black')
    
    ax.set_title(str(slices_30[i])+' km/s')
    
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
    
    
    ax.add_patch(Rectangle((-45,-50),90,100,color='black'))
#     ax.add_patch(Rectangle((-45,-50),90,100,color='black'))

    ax.add_patch(Rectangle((-360,-5),720,10,color='black'))

    
    ax.set_xlim(200,-200)
    plt.show()
    

slices_50 = np.array([(-100, -50),(-90, -40),(-80, -30),(-70, -20),(-60, -10),(-50, 0),\
(-40, 10),(-30, 20),(-20, 30),(-10, 40),(0, 50),(10, 60),(20, 70),\
 (30, 80),(40, 90),(50, 100)])


slices_50i = np.array((slices_50+100)/1.3).astype('int')


for i,ii in enumerate(slices_50i):
    hilb = np.nansum(hi4pi.data[ii[0]:ii[1],:,:],axis=0)
    colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    
    fig,ax = plt.subplots(1,1,figsize=(17,8))

    ax.imshow(hilb*hi_conv.value,
               origin='lower',norm='log',extent=[360,-360,-50,50])
    ax.contour(lons,lats,hilb*hi_conv.value,
                levels=[2,5,10,15,20],cmap='RdYlBu',zorder=0)
    ax.contourf(lons,lats,colb*co_conv.value,
                 levels=[10,1e5],colors='black')
    
    ax.set_title(str(slices_50[i])+' km/s')
    
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
    
    
    ax.add_patch(Rectangle((-45,-50),90,100,color='black'))
#     ax.add_patch(Rectangle((-45,-50),90,100,color='black'))

    ax.add_patch(Rectangle((-360,-5),720,10,color='black'))

    
    ax.set_xlim(200,-200)
    plt.show()
    

# - emission with Sigma > 10 still feels pretty sparse in these maps -> possibly because the high-surface-density hi that is gmc-less is within the plane of the galaxy? 
# 
#     - from looking at fig 5 from blitz 2005 (https://arxiv.org/pdf/astro-ph/0501400.pdf) there is definitely a lot of hi without gmcs -- but b information isnt present there (would it make sense to assume that those filaments are at low b? is that testable?)
#     - MD+2017s figs 26-32 simulate a face-on view of the MW molecular clouds, possibbly analogous to the m33 image, so comparing this to a similar map of the hi might be a better way of finding gmc-less hi? (this would be going back to large-scale survey decomposition though, rather than using known clouds with known distances etc)
# 
# - and because the galactic plane (|b| <=5) is being excluded here, then it is hard to see those clouds (or at least see enough of them / large enough areas of high Sigma) & use them as a control group -- maybe the thing to be asking is: how did the high-surface-density-and-high-latitude hi (+ their associated GMCs) form?
