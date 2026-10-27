#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits

orioncube = fits.open("DHT27_Orion_raw.fits")[0].data

op = fits.open("DHT27_Orion_raw.fits")
op0 = op[0]
dat = op0.data
print(type(op), type(op0), type(dat))

# information from the header:
#     (the correct labelling is what's in the __ = 'name' , not /comment
#     when looking at how many numbered aces)
#     axis 1 - 'velo-lsr' (ctype3)
#     axis 2 - 'glon-car' (ctype1)
#     axis 3 - 'glat-car' (ctype2)
#
# naxis1 = 209
# naxis2 = 153
# naxis3 = 77
#
# the axes order that gives "the right" visual when scrolling is: 231
#
# there are 77 velocity slices to slide through (the other axes order that gives 77 slices
# to slide through is 321, which makes it seem like '1' is where velocity is?)
#
# reference / data here:
# https://lweb.cfa.harvard.edu/rtdc/CO/NumberedRegions/DHT27/index.html
#
# reference images have galactic latitude | this way & galactic longitudes -- this way
# number of data points along latitude - 153, number of data points along longitude - ~200
#
# length of:
#     orioncube = 153
#     orioncube[0] = 209
#     orioncube[0][0] = 77

# try collapsing along different axes to see which gives the correct image
axis0 = np.nanmedian(orioncube, axis=0)
axis1 = np.nanmedian(orioncube, axis=1)
axis2 = np.nanmedian(orioncube, axis=2)

#axis0 is 209 x 77
#axis1 is 153 x 77
#axis2 is 153 x 209 <- this is what we want to plot

plt.figure(figsize=(15,15), dpi=80)
vmin2, vmax2 = zlims(axis2)
plt.imshow(axis2, vmin=vmin2, vmax=vmax2, origin='lower')
plt.colorbar(orientation='horizontal')
plt.show()

#okay, now i want to try and extract each velocity slice separately
num1 = len(orioncube)
num2 = len(orioncube[0])
num3 = len(orioncube[0][0])

tot = num1 * num2 * num3
vel_slices = np.zeros(tot)
vel_slices = vel_slices.reshape(77,153,209)

for i in range(num1):
    for j in range(num2):
        for k in range(num3):
            vel_slices[k][i][j] = orioncube[i][j][k]

plt.figure(figsize=(15,15), dpi=80)
fullvel = np.nanmedian(vel_slices, axis=0)
plt.imshow(fullvel, vmin=vmin2, vmax=vmax2, origin='lower')

#examining by cross comparing with ds9 slider
im = vel_slices[37]
plt.figure(figsize=(15,15), dpi=80)
vmin, vmax = zlims(im)
plt.imshow(im, vmin=vmin, vmax=vmax, origin='lower')
