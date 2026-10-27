#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astrodendro import Dendrogram
from astropy.io import fits

import datacube as dc


dht = fits.open('DHT_CO.fits')[0]
hi4pi_raw = fits.open('HI4PI.fits')[0]
hi4pi_calib = fits.open('HI4PI_flatcalib.fits')[0]


co_dend = Dendrogram.load_from('COGAL_dendrogram.fits') 


solar_catalog = pd.read_csv('solar_neighbors.csv')


dht.shape


solar_catalog


def CO_0th_moment(cloud,cat=solar_catalog,dend=co_dend,survey=dht):
    if type(cloud) == str:
        ind = np.where(cat.name == cloud)[0][0]
        
        
    elif type(cloud) == int:
        ind = cloud
    
    idx = cat._idx.iloc[ind]
    row = cat.iloc[ind]
    
    mask = dend[idx].get_mask()
    this = mask*survey.data
    spc = dc.spectrum1D(this)
    
    isol = mask*this
    
#     print(row)
    
    w = np.nansum(isol,axis=0)
    lc,bc,r,vsig,v = \
    row['lon<deg>'],row['lat<deg>'],row['radius_GMC<deg>'],row['vel_dispersion_GMC<km/s>'],row['vel<km/s>']
    
    lims_lb = [180.125,-179.875,-30.125,29.875] #just hardcoded from header
    lims_lv = [180.125,-179.875,-49.4,48.1]
    
    plt.plot(np.arange(v-vsig, v-vsig+(76*1.3),1.3),spc)
    
#     plt.imshow(w,cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb)
#     plt.title(cat.name[ind])
#     plt.xlim(lc+2.5*r,lc-2.5*r)
#     plt.ylim(bc-2.5*r,bc+2.5*r)
    plt.show()


CO_0th_moment('Taurus')


#make a way of visaulizing as many spectra as possible so you can both compare them within (co<->hi)
#and with eacher (how does spectrum 1 compare to spectrum 2)


for name in solar_catalog.name:
    CO_0th_moment(name)


# for DHT:
# l0 = 180.125
# b0 = -30.125
# v0 = -49.4
# 
# dl = -0.125
# db = 0.125
# dv = 1.3
# 
# nv = 76
# nb = 481
# nl = 2881


solar_catalog.iloc[0][['lon<deg>','lat<deg>']]
