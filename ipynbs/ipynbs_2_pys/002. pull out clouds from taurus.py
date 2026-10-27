#!/usr/bin/env python
# coding: utf-8

# goal with this notebook : set up methods etc to extract a cloud from the taurus(dht21) individual survey


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits


# from skimage.segmentation import clear_border
# from skimage import measure
# from skimage.measure import label,regionprops
# from scipy import ndimage as ndi
from scipy.ndimage import measurements, center_of_mass, binary_dilation, zoom
# import plotly.graph_objects as go


# data: https://lweb.cfa.harvard.edu/rtdc/CO/NumberedRegions/DHT21/index.html
# 
# cloud locations from imara & burkhart page 2: https://arxiv.org/pdf/1609.04817.pdf
# (rounded to whole numbers)
#  


#create a dataframe to contain cloud locations
columns = ['longitude (l)', 'latitude (b)', 'central velocity (km/s)']
CLOUDS = pd.DataFrame(columns = columns)


#then set the relevant details that i know so far
CLOUDS.at['Perseus','longitude (l)'] = (155,162)
CLOUDS.at['Perseus','latitude (b)'] = (-27,-15)

CLOUDS.at['Ophiuchus','longitude (l)'] = (-10,8)
CLOUDS.at['Ophiuchus','latitude (b)'] = (8,22)

CLOUDS.at['Rosette', 'longitude (l)'] = (204,210)
CLOUDS.at['Rosette', 'latitude (b)'] = (-4,-1)

CLOUDS.at['Orion B', 'longitude (l)'] = (201,210)
CLOUDS.at['Orion B', 'latitude (b)'] = (-18, -10)

CLOUDS.at['Orion A', 'longitude (l)'] = (205,217)
CLOUDS.at['Orion A', 'latitude (b)'] = (-21, 18)

CLOUDS.at['California', 'longitude (l)'] = (155,169)
CLOUDS.at['California', 'latitude (b)'] = (-15,-5)


# one of future steps to edit this is to change the function unpackDataCube into a class that has fix axes & other methods


#define a function that takes in the name string of a fits data cube
#& extracts information ? .. not sure all the things yet, but definitely
#read in the image, read in the data, find out which axis == velocity


# class unpackDataCube (name):
#     def __init__(self,name)
#     def getImage(self)
#         return image

def unpack_fitsCube(name):
    #----- open the image and the header from fits ---------------------------
    image = fits.open(name)[0].data
    header = fits.open(name)[0].header
    #-------------------------------------------------------------------------
    
    
    #--- extract information labelling velocity, longitude, latitude axes ----
    
    #figure out which axes are which
    for i in range(3):
        num = i+1 #axes are 1, 2, or 3
        if header['CTYPE'+str(num)] == 'VELO-LSR':
            velLen = header['NAXIS'+str(num)]
            velMin = header['CRVAL'+str(num)]
            velDel = header['CDELT'+str(num)]
            
        if header['CTYPE'+str(num)] == 'GLON-CAR':
            lonLen = header['NAXIS'+str(num)]
            lonMin = header['CRVAL'+str(num)]
            lonDel = header['CDELT'+str(num)]
            
        if header['CTYPE'+str(num)] == 'GLAT-CAR':
            latLen = header['NAXIS'+str(num)]
            latMin = header['CRVAL'+str(num)]
            latDel = header['CDELT'+str(num)]
            latRef = header['CRPIX'+str(num)]
        
        
#     #old clunky way for figuring out axes
#     #figure out which axes are which
#     if header['CTYPE1'] == 'VELO-LSR':
#         velAxis = header['NAXIS1']
#         velMin = header['CRVAL1']
#         velDel = header['CDELT1']
        
#         if header['CTYPE2'] == 'GLON-CAR':
#             lonAxis = header['NAXIS2']
#             lonMin = header['CRVAL2']
            
            
#             latAxis = header['NAXIS3']
#         elif header['CTYPE3'] == 'GLON-CAR':
#             lonAxis = header['NAXIS3']
#             latAxis = header['NAXIS2']


#     elif header['CTYPE2'] == 'VELO-LSR':
#         velAxis =  header['NAXIS2']
#         velMin = header['CRVAL2']
#         if header['CTYPE1'] == 'GLON-CAR':
#             lonAxis = header['NAXIS1']
#             latAxis = header['NAXIS3']
#         elif header['CTYPE3'] == 'GLON-CAR':
#             lonAxis = header['NAXIS3']
#             latAxis = header['NAXIS1']
    
#     elif header['CTYPE3'] == 'VELO-LSR':
#         velAxis =  header['NAXIS3']
#         if header['CTYPE2'] == 'GLON-CAR':
#             lonAxis = header['NAXIS2']
#             latAxis = header['NAXIS1']
#         elif header['CTYPE1'] == 'GLON-CAR':
#             lonAxis = header['NAXIS1']
#             latAxis = header['NAXIS2']
    
    #------------------------------------------------------------------------
    
    
    #----- rearrange fits cube axes ------------------------------------------
    
    #for how python displays / manipulates arrays, we want 
    #axis0 = velocity, axis1=latitude, & axis2=longitude
    
    #this whole thing is assuming that all three axes have different lengths, 
    #will run into a problem if that is not the case
    
    #then reshape the cube so that velocity slices are
    #along axis=0 
    oldShape = image
    x, y, z = oldShape.shape
    
   
    image = np.zeros((velLen,latLen,lonLen))
    
    oldAxes = np.array((x,y,z))
    
    #first check if it is already in the right dimensions
    if ((x == velLen) and (y==latLen) and (z==lonLen)):
        image = oldShape   
    
    #check conditions that would only require one swap
    elif ((x == velLen) or (y==latLen) or (z==lonLen)):
        if x == velLen:
            image = np.swapaxes(oldShape, 1, 2)
            
        elif y == latLen:
            image = np.swapaxes(oldShape, 0, 2)
        
        elif z == lonLen:
            image = np.swapaxes(oldShape, 0, 1)
   
    else: #two swaps need to be made
        #indices are to get the actual number
        shuffle = np.where(oldAxes == velLen)[0][0] 
        
        #if shuffle = 1 that means axes need to shuffle <- way
        #if shuffle = 2 that means axes need to shuffle -> that way
        image = np.swapaxes(oldShape, 0, 1)
        
        if shuffle == 1: 
            image = np.swapaxes(image, 1, 2)
        elif shuffle == 2:
            image = np.swapaxes(image, 0, 2)
            
    #---------------------------------------------------------------------------
    
    
    #----- create the arrays holding the real/physical axis values -------------
    velMax = velMin + velLen*velDel
    lonMax = lonMin + lonLen*lonDel
    latMax = latMin + latLen*latDel
    
    
    vels = np.arange(velMin, velMax, step=velDel)
    lons = np.arange(lonMin, lonMax, step=lonDel)
#     lats = np.arange(latMin, latMax, step=latDel)

    lats = np.arange(latLen)
    lats = (lats - latRef)*latDel + latMin
    

    return image, vels, lons, lats


#need to define a function that converts between 
#numpy integer indices and the physical values they represent along specific axes

def indToPhys (arrow, array, num):
    #arrow will be: 'get index' or 'get physical'
    #and the array will be the velocity, latitude
    #or longitude arrays that come out of unpack_fitsCube
    
    #num will either be an array of integers (paired with 'get physical')
    #or an array of floats/integers/some number type (paired with 'get index')
    
    output = np.zeros_like(num)
    
    
    if arrow == 'get physical':
        for i in range(len(num)):
            it = num[i]
            output[i] = array[it]    
        #will return the physical value stored in array if num is an 
        #integer within the bounds of array, will throw an error if not
        
        
    #as things stand, this returns an index no matter what
    #want to eventually tweak so it throws an error if the 
    #value being requested is outside of the bounds
    #(what would be a good buffer amount?)

    if arrow == 'get index':
        #define a buffer -- what it means when a value 
        #is not in the bounds of array
        #this assumes all values in array are equally spaced
        #define the buffer as two times this separation
        #(is there a diff number that would be better?)
        delta = 2*(array[1] - array[0])
        
        #check the sign of delta & determine whether the array 
        #is in increasing or decreasing order, then set bounds accordingly
        if delta > 0:
            lowest = array[0] - delta
            highest = array[-1] + delta
        if delta < 0:
            lowest = array[-1] - delta
            highest = array[0] + delta
        
        #now check that the numb being requested 
        #is within these bounds and throw an error if not
        
        
        for i in range(len(num)):
            it = num[i]
            if ((it >=lowest) and (it <=highest)):
                diff = abs(array - it)
                closest = np.argmin(diff)

                output[i] = closest
        
            #right now it gives nan, but i want it 
            #to throw an error/flag 
            else:
                output[i] = np.nan
            
    
    return output


# taurus = unpack_fitsCube('DHT21_Taurus_interp.fits')
galaxy = unpack_fitsCube('COGAL_all_interp.fits')


# taurus_cube = taurus[0]
# taurus_lat = taurus[3]
# taurus_lon = taurus[2]
# taurus_vel = taurus[1]

#fix latitudes because the header info != physical values
# taurus_lat = -1*taurus_lat


galaxy_cube = galaxy[0]
galaxy_lat = galaxy[3]
galaxy_lon = galaxy[2]
galaxy_vel = galaxy[1]


#changing latitude axis to match the reference plot
# galaxy_lat = (-1*galaxy_lat) + 30


galaxy_lat


galaxy_lat = (-1*galaxy_lat) + 30

galaxy_lat


# plt.figure(figsize=(10,10), dpi=200)

# im = galaxy_cube[245]
# im_lon = galaxy_lon
# im_lat = galaxy_lat


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin = vmin, vmax = vmax,\
#            extent = [im_lon[0],im_lon[-1], im_lat[0],im_lat[-1]], origin='lower')

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='horizontal')
# # plt.ylim(-100,600)

# plt.show()


galaxy_0th = np.nansum(galaxy_cube, axis=0)


plt.figure(figsize=(10,10), dpi=100)

im = galaxy_0th
# im = taurus_cube[45]
im_lon = galaxy_lon
im_lat = galaxy_lat


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax,\
           extent = [im_lon[0],im_lon[-1], im_lat[-1],im_lat[0]],origin='lower')

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal')
# plt.ylim(-100,600)


plt.show()

# plt.imshow(im,origin='lower')

# plt.show()


# #okay now going to try extracting perseus

# #first get its longitude and latitude bounds from CLOUDS
# perseus_lons_p = CLOUDS['longitude (l)']['Perseus']
# perseus_lats_p = CLOUDS['latitude (b)']['Perseus']


#extract california
california_lons_p = CLOUDS['longitude (l)']['California']
california_lats_p = CLOUDS['latitude (b)']['California']


california_lons_i = indToPhys('get index', galaxy_lon, california_lons_p)
california_lats_i = indToPhys('get index', galaxy_lat, california_lats_p)

#since things are in descending order, sort the index lists to make indexing
#into the fits easier later
california_lons_i.sort()
california_lats_i.sort()


# #then find the integer bounds corresponding to those numbers
# perseus_lons_i = indToPhys('get index', galaxy_lon, perseus_lons_p)


# perseus_lats_i = indToPhys('get index', galaxy_lat, perseus_lats_p)


# # perseus_lons

# # pers_b_min, pers_b_max

# # perseus_lons[1]


california_lats_i


cal_lat_0


plt.figure(figsize=(10,10), dpi=100)

im = galaxy_cube[243]
# im = taurus_cube[45]
im_lon = galaxy_lon
im_lat = galaxy_lat


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax,\
           extent = [im_lon[0],im_lon[-1], im_lat[-1],im_lat[0]],origin='lower')

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal')
# plt.ylim(-100,600)


plt.show()

# plt.imshow(im,origin='lower')

# plt.show()


#trying out different methods to mask

#based on number of counts
# mask = taurus_cube > 2

#based on position
cal_lat_0, cal_lat_1 = california_lats_i - 150 #adjusted because lat probs
cal_lon_0, cal_lon_1 = california_lons_i


california = galaxy_cube[:, cal_lat_0:cal_lat_1, cal_lon_0:cal_lon_1]


plt.figure(figsize=(10,10), dpi=100)

im = california[247]
# im = taurus_cube[45]
im_lon = galaxy_lon
im_lat = galaxy_lat


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
           extent = [california_lons_p[1],california_lons_p[0],\
                     california_lats_p[0], california_lats_p[1]])

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal')

plt.title('California extracted from galaxy')

plt.show()


#attempt a velocity spectrum for califonia ^

california_vel_counts = [None]*(len(california))

for i in range(len(california)):
    vel_slice = california[i]
    sum_acros_long = np.nanmedian(vel_slice, axis = 0)
    vel_count = np.nanmedian(sum_acros_long)
    california_vel_counts[i] = vel_count

# plt.stairs(california_vel_counts, color = 'xkcd:lilac', label = 'indices')
plt.plot(galaxy_vel,california_vel_counts,color='xkcd:indigo', label = 'km/s')

plt.title('califonia velocity spectrum raw')

plt.ylabel('counts') #should be temperature (K)
#check the units again


# plt.xlabel('velocity (right now in indices, will turn into )')
plt.xlim(-40,40)
# plt.legend()
plt.show()


# 


# #not using rn, but this will be useful when adding another 
# #step in the mask that does a cutoff based on sigma
# mask_new = binary_dilation(mask, iterations=5)
# plt.pcolormesh(mask_new[45])
# img_new = mask_new*galaxy_cube


# plt.figure(figsize=(8,8), dpi=200)

# plt.pcolormesh(img_new[45], cmap='')
# plt.colorbar(orientation ='horizontal')


# qs for check in: 2.08.23
# 
# so turns out that the header is wrong?? it says:
# latitude values start at 0, have + spacing values
# ... even though the actual data, when cross referenced w referece plots, looks different
# 
# 
# 1. in taurus (dht21): lat should be 0->-28, but is 0->26.5, so perseus' position is labelled out of bounds --- and the taurus 0th moment map i plot looks "squished" along the latitude axis in comparison to plot from data site
#     - also, real latitude values decrease from 0 to -28, but fits has increase to +26.5
#     
# 2. in galaxy, lat should be +30 -> -30, but is instead 0 -> +60
# 
# 3. another aspect of the scaling being off -> more obvious on larger scales -> i had to adjust latitude indices by 150 to find california in the galaxy data 
# 
# -- for this, use the reference pixel + screenshot of nia's code
# 
# 4. check if the velocity spectrum for california looks ok (search more to find paper with california spectrum specifically)
#     -- look at the orions and compare with the imara papers

# next steps 2.8.23
# - get comfy making 0th moment maps for all the clouds i look at
# - make spectra of orions + perseus and compare to the published spectra
# - calculate the velocity dispersion for ^
# 


header = fits.open('COGAL_all_interp.fits')[0].header


header


# moving forward notes: from 1.31.23
# 
# i would begin by defining a box in (lat,lon) and then tell the code that the cloud is in there
# count threshold: 3 sigma or 5 sigma is good
# and after spatially defining the cloud, it's okay to do the velocity spectrum from there (don't need to do velocity cuts)
# 
# for lat/lon bounds: look at imara papers (imara & blitz gmc paper)
# 
# Lada, lombardi, alves have many papers outlining the bounds of clouds
# 
# and build my own list
# 
# start out with identifying the physical location of a single cloud & do spectra that way
# 
# good idea to use the whole data set rather than in pieces
# 
# build in a way to say the name of the cloud & get the location back
# 
# first: put xy axes in lat/lon coordinates
# 
# spatial cut first & then velocity/threshold cutoff
# 
# let nia know where to put a directory for me -- send her a screenshot as eg for how ale is in the ramirez-ruiz folder etc
