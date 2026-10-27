#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm

from astropy.io import fits
import astropy.units as u
import astropy.constants as c


hi4pi, dht = fits.open('../fitses/both_same_grid_720d.fits')
dust = fits.open('../fitses/mean_and_std_lbd.fits')[1].data


d0_NH = np.nansum(dust*2.8,axis=0)*1.87e21


plt.figure(dpi=200)
plt.contourf(np.log10(d0_NH),levels=20,cmap='magma',extent=[180,-180,-90,90])


plt.xlim(180,-180)


dld_NH = np.nansum(dust*2.8,axis=1)*0.229*1.87e21


dbd_NH = np.nansum(dust*2.8,axis=2)*0.229*1.87e21


cm.PiYG


plt.figure(dpi=200)
plt.contourf(np.log10(dld_NH),levels=20,cmap='magma',extent=[180,-180,68,1250])


plt.xlim(180,-180)


ds = np.arange(68,1250)
bs = (np.arange(786)-394)*0.229
ls = (np.arange(1572)-786)*-0.229


zs = ds[:,np.newaxis]*np.sin(bs[np.newaxis,:]*u.deg).value


plt.figure(dpi=200)
plt.contourf(np.log10(dbd_NH).transpose(),levels=20,cmap='magma',extent=[68,1250,-90,90])

plt.contour(zs.transpose(),levels=np.arange(-1000,1000,100),extent=[68,1250,-90,90],cmap='PiYG')
