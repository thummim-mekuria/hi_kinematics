#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting

import fun as fns
import datacube as dc

#---------------------------------------------------------------------------------------------
fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


#-------------------------------------------------------------------------------------------
hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data

h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)

h0 = np.nansum(h3sig,axis=0)*1.3*1.82e18
c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20

hlv = np.nansum(h3sig,axis=1)*0.125*1.82e18
clv = np.nansum(d3sig,axis=1)*0.125*1.97e20

#the 0.5 is to account for the fact that this is 720 deg
hbv = np.nansum(h3sig,axis=2)*0.125*1.82e18*0.5 
cbv = np.nansum(d3sig,axis=2)*0.125*1.97e20*0.5


#---------------------------------------------------------------------------------------------
lcs_d = np.nanmean(fcat['lon_GMC<deg>'],axis=1)
bcs_d = np.nanmean(fcat['lat_GMC<deg>'],axis=1)

ras_d = (fcat['radius_HI<pc>']/fcat['distance<pc>'])*u.rad.to(u.deg)


loh,lih = ((np.array([lcs_d+ras_d,lcs_d-ras_d]).transpose() -360)/-0.125).astype('int').transpose()
boh,bih = ((np.array([bcs_d-ras_d,bcs_d+ras_d]).transpose()+50)/0.125).astype('int').transpose()

loc,lic = ((fcat['lon_GMC<deg>'] -360)/-0.125).astype('int').transpose()
boc,bic = ((fcat['lat_GMC<deg>'] +50)/0.125).astype('int').transpose()

vo,vi = np.array([((fcat['vel<km/s>'][i]+100)/1.3).astype('int') for i in range(22)]).transpose()


vs = [np.arange((vi-vo)[i])*1.3 + fcat['vel<km/s>'][i][0] for i in range(22)]


his = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]
cos_hisz = [d3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
his_cosz = [h3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


#---------------------------------------------------------------------------------------------
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


fits.open('mean_and_std_lbd.fits')[1].header


lons = (np.arange(1572)-787)*-0.229
lats = (np.arange(786)-394)*0.229
dist = np.arange(1182)+68


lats[275],lats[375]


lons[1350]


lons[10],lons[1571]


dist[132],dist[532]


orion = dust[:,275:400,1350:]


perseus = dust[:,150:400,:200]


from copy import deepcopy


temp = deepcopy(fits.open('mean_and_std_lbd.fits')[1].header)


temp.remove('EXTNAME')


temp['CTYPE1'] = 'distance'
temp['CTYPE3'] = 'GLON-CAR'


temp['CUNIT1'] = 'pc'
temp['CUNIT3'] = 'deg'


temp['NAXIS1'] = 1182
temp['NAXIS2'] = 125
temp['NAXIS3'] = 222


temp['CRPIX1'] = 0.
temp['CRPIX2'] = 0.
temp['CRPIX3'] = 0.


temp['CRVAL1'] = 68.
temp['CRVAL2'] = -2.7251e1
temp['CRVAL3'] = -1.28927e2


temp['CDELT3'] = -0.22900763358778625
temp['CDELT1'] = 1.


for thing in ['AUTHOR','DATE','VERSION','REF']:
    temp.remove(thing)


temp.remove('CONTACT')


temp.add_comment = 'orion subcube extracted from edenhoffer et al 2023 data'


fits.PrimaryHDU(orion,header=temp).writeto('dust_orion.fits',overwrite=True)


temp2 = deepcopy(temp)


perseus.shape


temp2['NAXIS2'] = 250
temp2['NAXIS3'] = 200


temp2['CRVAL2'] = -5.5876e1
temp2['CRVAL3'] = 1.80223e2


fits.PrimaryHDU(perseus,header=temp2).writeto('dust_perseus.fits')


len(fits.open('mean_and_std_lbd.fits'))


lons[1350],lons[-1]


plt.imshow(np.nansum(dust[:,150:400,:200],axis=0),
           cmap=dc.cmap_min0('viridis'),origin='lower')#,levels=[0.01,0.05,0.1,0.2,0.5])
