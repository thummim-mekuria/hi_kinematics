#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import fitting,models
fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.)

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit as fit

def line(x,m,b):
    return m*x+b

import datacube as dc
import fun as fns
import interpolate as ip


cat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


hi4pi,dht = fits.open('both_same_grid_720d.fits')


cos_hsz = [
dht.data[cat['vo'][i]:cat['vi'][i],cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]]
    
for i in range(22)]


def make_fits_from_list(array,filename):
    hdulist = fits.HDUList([fits.PrimaryHDU(array[0])]+[fits.ImageHDU(array[i]) for i in range(1,22)])
    hdulist.writeto(filename,overwrite=True)


make_fits_from_list(cos_hsz,'h2_hisz_clouds.fits')


c0 = np.nansum(dht.data,axis=0)*1.3*1.97e20


h0 = np.nansum(hi4pi.data,axis=0)*1.3*1.82e18


fits.PrimaryHDU(c0).writeto('allsky_0th_co.fits')


fits.PrimaryHDU(h0).writeto('allsky_0th_hi.fits')


h1 = fns.get_mom1(hi4pi.data,np.arange(155)*1.3-100)


c1 = fns.get_mom1(dht.data,np.arange(155)*1.3-100)


fits.PrimaryHDU(c1).writeto('allsky_1st_co.fits')


plt.figure(dpi=200)
plt.contourf(h1,levels=np.arange(-50,50,5),extend='both',cmap='Spectral',extent=[360,-360,-50,50])
plt.xlim(180,-180)
