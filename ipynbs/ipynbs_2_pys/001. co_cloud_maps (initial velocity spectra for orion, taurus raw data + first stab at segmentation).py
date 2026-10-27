#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits

import photutils as phot

from PIL import Image, ImageFilter


# co_cube = fits.open("COGAL_all_interp.fits")[0].data


# try collapsing along different axes to see which gives the correct image
# axis0 = np.nansum(co_cube, axis = 0)
# axis1 = np.nansum(co_cube, axis = 1)
# axis2 = np.nansum(co_cube, axis = 2)

# print(len(axis0), len(axis0[0]))
# print(len(axis1), len(axis1[0]))
# print(len(axis2), len(axis2[0]))

#axis0 is 2881 x 493 <- this is the one we want


# axis0 = np.nansum(co_cube, axis = 0)

# co_galaxy_0th = axis0.transpose()

# plt.figure(figsize=(15,15), dpi=80)

# vmin_gal, vmax_gal = zlims(co_galaxy_0th)
# plt.imshow(co_galaxy_0th,vmin = vmin_gal, vmax = vmax_gal,  origin = 'lower')
# plt.colorbar(orientation = 'horizontal')
# plt.show()


#okay now do individual clouds
orion_cube = fits.open("DHT27_Orion_interp.fits")[0].data
taurus_cube = fits.open("DHT21_Taurus_interp.fits")[0].data


orion_0th = np.nansum(orion_cube, axis = 2)


plt.figure(figsize=(8,8), dpi=200)

vmin_or, vmax_or = zlims(orion_0th)
plt.imshow(orion_0th,vmin = vmin_or, vmax = vmax_or,  origin = 'lower')
plt.colorbar(orientation = 'horizontal')
plt.title('Orion in CO')
plt.show()


# subtract noise - set a threshold for when something counts as part of the cloud (eg 3* sigma)
# spatially separate orion a and orion b first, and *then* do the velocity differentiation
# 
# for extinction maps, the boundary of the cloud is defined at Av = 1


taurus_0th = np.nansum(taurus_cube, axis = 2)


len(taurus_0th[0])


plt.figure(figsize=(8,8), dpi=80)

vmin_tau, vmax_tau = zlims(taurus_0th)
plt.imshow(taurus_0th,vmin = vmin_tau, vmax = vmax_tau,  origin = 'lower')
plt.colorbar(orientation = 'horizontal')
plt.title('Taurus in CO')
plt.show()


# taurus, california, perseus are all up here^
# 
# the unit of the intensity / count in 0th moment maps is - K km/s
# ^take this and use the X_co factor to go to N(H2)
# ^calculating the column density happens pixel by pixel & using the co intensity at that area
# 
# -- start with defining a mask -> after spatially isolating the cloud, go through & pick out the pixels that count as part of the cloud, i.e. are above the noise threshold & only do calculations for those pixels


# #code for checking which axis to add along
# file = taurus_cube
# axis0 = np.nansum(file, axis = 0)
# axis1 = np.nansum(file, axis = 1)
# axis2 = np.nansum(file, axis = 2)

# print(len(axis0), len(axis0[0]))
# print(len(axis1), len(axis1[0]))
# print(len(axis2), len(axis2[0]))


#okay, now i want to rearrange the cube - make it easier to extract each velocity slice separately
#orion
or_num1 = len(orion_cube)
or_num2 = len(orion_cube[0])
or_num3 = len(orion_cube[0][0])

or_tot = or_num1 * or_num2 * or_num3
or_vel_slices = np.zeros(or_tot)
orion_cube_old_shape = orion_cube
orion_cube = or_vel_slices.reshape(77,153,209)


for i in range (or_num1):
    for j in range (or_num2):
        for k in range (or_num3):
            orion_cube[k][i][j] = orion_cube_old_shape[i][j][k]
            

#okay, now i want to try and extract each velocity slice separately
tau_num1 = len(taurus_cube)
tau_num2 = len(taurus_cube[0])
tau_num3 = len(taurus_cube[0][0])

tau_tot = tau_num1 * tau_num2 * tau_num3
tau_vel_slices = np.zeros(tau_tot)
taurus_cube_old_shape = taurus_cube
taurus_cube = tau_vel_slices.reshape(81,213,307)


for i in range (tau_num1):
    for j in range (tau_num2):
        for k in range (tau_num3):
            taurus_cube[k][i][j] = taurus_cube_old_shape[i][j][k]


len(taurus_cube[0])


# #examining by cross comparing with ds9 slider
# im = orion_cube[37]

# plt.figure(figsize=(7,7), dpi=80)
# vmin, vmax = zlims(im)
# plt.imshow(im, vmin = vmin, vmax = vmax, origin = 'lower')


# #examining by cross comparing with ds9 slider
# im = taurus_cube[48]

# plt.figure(figsize=(7,7), dpi=80)
# vmin, vmax = zlims(im)
# plt.imshow(im, vmin = vmin, vmax = vmax, origin = 'lower')


#first make spectrum of whole cube & then separate out different clouds based
#on different emission peaks

#start with orion

or_vel_counts = [None]*(len(orion_cube))

for i in range(len(orion_cube)):
    vel_slice = orion_cube[i]
    sum_acros_long = np.nansum(vel_slice, axis = 0)
    vel_count = np.nansum(sum_acros_long)
    or_vel_counts[i] = vel_count

plt.stairs(or_vel_counts)
plt.title('orion velocity spectrum raw')

plt.show()


# read in velocity start coordinate + spacing & use that to create the array, not values from ds9


#assign real velocity range values to indices

#for orion
#start (index = 0) = -9.27200
#spacing = 0.6501900
#end = 39.9872


#for taurus
#start (index = 0) = -21.78010
#spacing = 0.6501900
#end = 30.2351


vel_spacing = 0.6501900

or_start = -9.27200
or_end = 39.9872

or_vels = np.linspace(or_start, or_end, 77)


or_vels


or_vels[62]


#first make spectrum of whole cube & then separate out different clouds based
#on different emission peaks

#start with orion

tau_vel_counts = [None]*(len(taurus_cube))

for i in range(len(taurus_cube)):
    vel_slice = taurus_cube[i]
    sum_acros_long = np.nansum(vel_slice, axis = 0)
    vel_count = np.nansum(sum_acros_long)
    tau_vel_counts[i] = vel_count

plt.stairs(tau_vel_counts)
plt.title('taurus velocity spectrum raw')

plt.show()


np.nanmax(orion_cube)


# questions / things to talk through
# 1. can visually see that clouds are distinct (orion a & b). and that they are at different-ish velocity ranges (taurus & california)
#     -> how to formalize this?
#     -> also, with orion a & b, what are the boundaries between them & why?
# 
# 2. with velocities -- to what level do we care about precision? (the arrays i'm trying to create based on start, stop, and spacing values from the headers aren't matching up with the velocity labels per slice)
# 
# 3. how do you get the temperature <- count conversion? from the header, there isn't a 'coord value' as there is for v, glong, glat
# 
# 
# 


maxtaurus = np.max(taurus_0th)
top50 = maxtaurus-50
top20 = maxtaurus-20
top70 = maxtaurus-70


import numpy.ma as ma


top50mask = ma.masked_where(taurus_0th < top50, taurus_0th)
top20mask = ma.masked_where(taurus_0th < top20, taurus_0th)
top70mask = ma.masked_where(taurus_0th < top70, taurus_0th)


taurus_0th[64][263]


# plt.figure(figsize=(8,8), dpi=200)

# vmin_tau, vmax_tau = zlims(taurus_0th)
# plt.imshow(taurus_0th,vmin = vmin_tau, vmax = vmax_tau,  origin = 'lower')
# plt.colorbar(orientation = 'horizontal')
# plt.title('Taurus in CO')
# plt.show()


plt.figure(figsize=(8,8), dpi=200)
plt.imshow(top50mask, origin = 'lower')


vmin_tau, vmax_tau = zlims(taurus_0th)
# plt.imshow(taurus_0th,vmin = vmin_tau, vmax = vmax_tau,  origin = 'lower')
plt.colorbar(orientation = 'horizontal')
plt.title('Taurus in CO')
plt.show()


plt.figure(figsize=(8,8), dpi=200)
plt.imshow(top20mask, origin = 'lower')


vmin_tau, vmax_tau = zlims(taurus_0th)
# plt.imshow(taurus_0th,vmin = vmin_tau, vmax = vmax_tau,  origin = 'lower')
plt.colorbar(orientation = 'horizontal')
plt.title('Taurus in CO')
plt.show()


plt.figure(figsize=(8,8), dpi=200)
plt.imshow(top70mask, origin = 'lower')


vmin_tau, vmax_tau = zlims(taurus_0th)
# plt.imshow(taurus_0th,vmin = vmin_tau, vmax = vmax_tau,  origin = 'lower')
plt.colorbar(orientation = 'horizontal')
plt.title('Taurus in CO')
plt.show()


# edge detection on python
# geopandas
# principla component analysis
# 


# Opening the image (R prefixed to string
# in order to deal with '\' in paths)
image = Image.open('sproots.png')

# Converting the image to grayscale, as edge detection
# requires input image to be of mode = Grayscale (L)
image = image.convert("L")

# Detecting Edges on the Image using the argument ImageFilter.FIND_EDGES
image = image.filter(ImageFilter.FIND_EDGES)

# # Saving the Image Under the name Edge_Sample.png
# image.save(r"Edge_Sample.png")


from scipy import ndimage as ndi
from skimage.util import random_noise
from skimage import feature


image = taurus_0th


edges = feature.canny(image, sigma=15)

# display results
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(32, 12))

ax[0].imshow(image, origin = 'lower')
ax[0].set_title('noisy image', fontsize=50)

# ax[1].imshow(edges1, cmap='gray')
# ax[1].set_title(r'Canny filter, $\sigma=1$', fontsize=20)

ax[1].imshow(edges, origin = 'lower')
ax[1].set_title(r'Canny filter, $\sigma=3$', fontsize=50)

for a in ax:
    a.axis('off')

fig.tight_layout()
plt.show()


from sklearn.decomposition import PCA


taurus_PCA = PCA(n_components=2)
taurus_PCA.fit(taurus_0th)


type(taurus_PCA)
