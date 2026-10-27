#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits

from astropy.stats import sigma_clipped_stats

from scipy.ndimage import binary_dilation
from specutils.analysis import fwhm
from specutils import Spectrum1D 

import astropy.units as u

from scipy.signal import peak_widths, find_peaks


# data: https://lweb.cfa.harvard.edu/rtdc/CO/CompositeSurveys/
# 
# cloud locations from imara & burkhart page 2: https://arxiv.org/pdf/1609.04817.pdf
# (rounded to whole numbers)
# 
# cloud velocity spectra: imara & blitz https://arxiv.org/pdf/1103.3741.pdf
#  


### THINK ABOUT TURNING THIS INTO AN EXTERNAL CSV OR STH I JUST READ IN

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
CLOUDS.at['Orion A', 'latitude (b)'] = (-21, -18)

CLOUDS.at['California', 'longitude (l)'] = (155,169)
CLOUDS.at['California', 'latitude (b)'] = (-15,-5)


### one of future steps to edit this is to change the function unpackDataCube into a 
###class that has fix axes & other methods

#something like:

# class unpackDataCube (name):
#     def __init__(self,name)
#     def getImage(self)
#         return image


#############################################   


#for now..
#define a function that takes in the name string of a fits data cube
#& returns the image (with axes rotated so that 0=vel,1=lat,2=lon)
#and arrays containing physical information for lat, lon, vel 


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
            velRef = header['CRPIX'+str(num)]
            
        if header['CTYPE'+str(num)] == 'GLON-CAR':
            lonLen = header['NAXIS'+str(num)]
            lonMin = header['CRVAL'+str(num)]
            lonDel = header['CDELT'+str(num)]
            lonRef = header['CRPIX'+str(num)]
            
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
    glat = np.arange(latLen)
    glat = (glat - latRef)*latDel + latMin
    
    glon = np.arange(lonLen)
    glon = (glon - lonRef)*lonDel + lonMin
    
    vels = np.arange(velLen)
    vels = (vels - velRef)*velDel + velMin
    
    
#     #old way of creating arrays, didn't account for reference pixel
    
#     velMax = velMin + velLen*velDel
#     lonMax = lonMin + lonLen*lonDel
#     latMax = latMin + latLen*latDel
    
    
#     vels = np.arange(velMin, velMax, step=velDel)
#     lons = np.arange(lonMin, lonMax, step=lonDel)
#     lats = np.arange(latMin, latMax, step=latDel)

    
    return header,image, vels, glon, glat, velDel


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


#define a function that takes in a data cube and returns a 1D spectrum (along axis=0)
def spectrumFromCube (cube):
    num = len(cube)
    
    spectrum = np.zeros(num)
    
    for i in range(num):
        s = cube[i]
        temp = np.nanmean(s, axis=0)
        count = np.nanmean(temp)
        
        spectrum[i] = count
        
    return spectrum


galaxy = unpack_fitsCube('COGAL_all_interp.fits')


galaxy[0]


galaxy_cube = galaxy[0]
gvels = galaxy[1]
glon = galaxy[2]
glat = galaxy[3]
gvelDel = galaxy[4]


galaxy_0th = np.nanmean(galaxy_cube, axis = 0)


plt.figure(figsize=(25,5), dpi=200)

im = galaxy_0th
imx = glon
imy = glat


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax, \
           extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower')

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

plt.title('Galaxy 0th moment map - summed')
# plt.ylim(-100,600)

plt.show()


#### ORION A ####

#first get its longitude and latitude bounds from CLOUDS
orionA_lon_phy = CLOUDS['longitude (l)']['Orion A']
orionA_lat_phy = CLOUDS['latitude (b)']['Orion A']


# not sure if this is a problem but i have run into:
#the coordinates for clouds are listed assuming that glon runs 0 to 360, but
#the info from datacube assumes glon runs -180 to 180

#so orion's lon coords are saved in dataframe (from lit) as 205,217
#but looking at 0th moment map shows that the coords are ~-150 

#the conversion is 0 <=> 360, and negative lon (0 down to -180) numbers can also be 360 -> 180
new_lons = []
#so do a check:
for i in orionA_lon_phy:
    if i > 180:
        i -= 360
        new_lons.append(i)

orionA_lon_phy360 = orionA_lon_phy 
orionA_lon_phy = new_lons


orionA_inGlon = indToPhys('get index', glon, orionA_lon_phy)
orionA_inGlat = indToPhys('get index', glat, orionA_lat_phy)

orionA_inGlon.sort()
orionA_inGlat.sort()


#### CREATE A MASK THAT LIFTS OUT ORION A FROM GALAXY CUBE
orionA_b_0, orionA_b_1 = orionA_inGlat
orionA_l_0, orionA_l_1 = orionA_inGlon

orionA_cube = galaxy_cube[:,orionA_b_0:orionA_b_1, orionA_l_0:orionA_l_1]


#calculate the background around the star
orA_mean, orA_med, orA_std = sigma_clipped_stats(orionA_cube, sigma=3.0)


#then create a mask of values above that threshold
orA_sigma_mask = orionA_cube > orA_std
# orA_mask_new = binary_dilation(orA_sigma_mask, iterations=1)


orionA_cube_noisy = orionA_cube
orionA_cube = orA_sigma_mask * orionA_cube_noisy
# orionA_cube_minusbkg_nw = orA_mask_new * orionA_cube


##plotting diff things to see what works best for the mask

# plt.figure(figsize=(12,4), dpi=100)

# im = orionA_cube_minusbkg
# im2 = orionA_cube_minusbkg_nw
# # imx = glon[orionA_l_0:orionA_l_1]+360 #put lon coords back in 360 form
# # imy = glat[orionA_b_0:orionA_b_1]


# # vmin,vmax = zlims(im)

# ind = 259
# plt.pcolormesh(im[ind])

# # plt.ylabel('Galactic Latitude')
# # plt.xlabel('Galactic Longitude')
# # plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

# # plt.title('Orion A 0th moment')

# plt.show()

# plt.figure(figsize=(12,4), dpi=100)

# plt.pcolormesh(im2[ind])
# plt.show()


# #not using rn, but this will be useful when adding another 
# #step in the mask that does a cutoff based on sigma
# mask_new = binary_dilation(mask, iterations=5)
# plt.pcolormesh(mask_new[45])
# img_new = mask_new*galaxy_cube


orionA_0th = np.nanmean(orionA_cube, axis=0)


plt.figure(figsize=(12,5), dpi=100)

im = orionA_0th
imx = glon[orionA_l_0:orionA_l_1]+360 #put lon coords back in 360 form
imy = glat[orionA_b_0:orionA_b_1]


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
           extent = [imx[0],imx[-1], imy[0], imy[-1]])

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

plt.title('Orion A 0th moment')

plt.show()


orionA_spectrum = spectrumFromCube(orionA_cube)

orA_peak_vel = find_peaks(orionA_spectrum, height = 0.5*np.nanmax(orionA_spectrum))
#feed ^ + spectrum into method that calculates width

orionA_linewidth = peak_widths(orionA_spectrum,orA_peak_vel[0])[0]
orionA_linewidth = orionA_linewidth*gvelDel


orionA_linewidth


# 


plt.plot(gvels,orionA_spectrum, color = 'gray')
plt.title('Orion A velocity spectrum')

plt.ylabel('CO Temperature (K)') #should be temperature (K)
plt.xlabel('Velocity (km/s)')

plt.show()


#### START THINKING ABOUT DEFINING A METHOD THAT TAKES IN THE NAME OF THE
#### STAR & THE GALAXY CUBE, THEN PULLS IT OUT, DOES 0TH MOMENT, AND OUTPUTS
#THE CUBE FOR THAT STAR, THE 0TH MOMENT, AND THE VELOCITY SPECTRUM


#### ORION B ####

#first get its longitude and latitude bounds from CLOUDS
orionB_lon_phy = CLOUDS['longitude (l)']['Orion B']
orionB_lat_phy = CLOUDS['latitude (b)']['Orion B']


new_lons = []
#so do a check:
for i in orionB_lon_phy:
    if i > 180:
        i -= 360
        new_lons.append(i)

orionB_lon_phy360 = orionB_lon_phy 
orionB_lon_phy = new_lons


orionB_inGlon = indToPhys('get index', glon, orionB_lon_phy)
orionB_inGlat = indToPhys('get index', glat, orionB_lat_phy)

orionB_inGlon.sort()
orionB_inGlat.sort()


#### CREATE A MASK THAT LIFTS OUT ORION A FROM GALAXY CUBE
orionB_b_0, orionB_b_1 = orionB_inGlat
orionB_l_0, orionB_l_1 = orionB_inGlon

orionB_cube = galaxy_cube[:,orionB_b_0:orionB_b_1, orionB_l_0:orionB_l_1]


#calculate the background around the star
orB_mean, orB_med, orB_std = sigma_clipped_stats(orionB_cube, sigma=3.0)


#then create a mask of values above that threshold
orB_sigma_mask = orionB_cube > orB_std


orionB_cube_noisy = orionB_cube
orionB_cube = orB_sigma_mask * orionB_cube_noisy


orionB_0th = np.nanmean(orionB_cube, axis=0)


plt.figure(figsize=(10,10), dpi=100)

im = orionB_0th
imx = glon[orionB_l_0:orionB_l_1]+360 #put lon coords back in 360 form
imy = glat[orionB_b_0:orionB_b_1]


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
           extent = [imx[0],imx[-1], imy[0], imy[-1]])

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

plt.title('Orion B 0th moment')

plt.show()


orionB_spectrum = spectrumFromCube(orionB_cube)


orB_peak_vel = find_peaks(orionB_spectrum, height = 0.5*np.nanmax(orionB_spectrum))
#feed ^ + spectrum into method that calculates width

orionB_linewidth = peak_widths(orionB_spectrum,orB_peak_vel[0])[0]
orionB_linewidth = orionB_linewidth*gvelDel


orionB_linewidth*gvelDel


plt.plot(gvels,orionB_spectrum, color = 'gray')
plt.title('Orion B velocity spectrum')

plt.ylabel('CO Temperature (K)') #should be temperature (K)
plt.xlabel('Velocity (km/s)')

plt.show()


#### PERSEUS ####
perseus_lon_phy = CLOUDS['longitude (l)']['Perseus']
perseus_lat_phy = CLOUDS['latitude (b)']['Perseus']


perseus_inGlon = indToPhys('get index', glon, perseus_lon_phy)
perseus_inGlat = indToPhys('get index', glat, perseus_lat_phy)

#sort so indices are in ascending order
perseus_inGlon.sort()
perseus_inGlat.sort()


#### CREATE A MASK THAT LIFTS OUT PERSEUS FROM GALAXY CUBE
perseus_b_0, perseus_b_1 = perseus_inGlat
perseus_l_0, perseus_l_1 = perseus_inGlon

perseus_cube = galaxy_cube[:,perseus_b_0:perseus_b_1, perseus_l_0:perseus_l_1]


#calculate the background around the star
per_mean, per_med, per_std = sigma_clipped_stats(perseus_cube, sigma=3.0)


#then create a mask of values above that threshold
per_sigma_mask = perseus_cube > per_std
# orA_mask_new = binary_dilation(orA_sigma_mask, iterations=1)


perseus_cube_noisy = perseus_cube
perseus_cube = per_sigma_mask * perseus_cube_noisy
# orionA_cube_minusbkg_nw = orA_mask_new * orionA_cube


#0th moment map
perseus_0th = np.nanmean(perseus_cube, axis = 0)


plt.figure(figsize=(10,10), dpi=100)

im = perseus_0th
imx = glon[perseus_l_0:perseus_l_1]
imy = glat[perseus_b_0:perseus_b_1]


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
           extent = [imx[0],imx[-1], imy[0], imy[-1]])

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

plt.title('Perseus 0th moment')

plt.show()


#extract perseus velocity spectrum
perseus_spectrum = spectrumFromCube(perseus_cube)


### first try: specutils

#when making the Spectrum1D object, want to eventuallly pull 
#unit information out, not hard code it like here
# perseus_spectrum = Spectrum1D(spectral_axis=gvels*u.km/u.s, flux=perseus_spectrum*u.K)

# perseus_linewidth = fwhm(perseus_spectrum)


### second try: scipy
#first need to find the index where the peak is 
#(set minimum height as 80% of the max) to avoid tiny bumps being = peaks

per_peak_vel = find_peaks(perseus_spectrum, height = 0.5*np.nanmax(perseus_spectrum))
#feed ^ + spectrum into method that calculates width

perseus_linewidth = peak_widths(perseus_spectrum,per_peak_vel[0])[0]
#the returns from peak_widths are:
#0. width for each peak in the sample
#1. the height of the contour line at which peak was evaluated
#2. interpolated positions of left 
#3. & right intersection points of the
#horizontal line corresponding to 2


#as is, the line width calculated is in index units -> gotta convert to velocity
perseus_linewidth = perseus_linewidth*gvelDel


perseus_linewidth


plt.plot(gvels,perseus_spectrum, color = 'gray')


plt.title('Perseus velocity spectrum')

plt.ylabel('CO Temperature (K)') #should be temperature (K)
plt.xlabel('Velocity (km/s)')

plt.show()


# qs for check in: 2.17.21
# 
# (not necessarily qs but) 
# 1. is there a specific longitude convention (0 to 360 vs -180 to 180)?
#     - papers list things as 0-360, but datacube has axes +/-180
#     - can do checking + fixing by hand, but wondering if i'll ever need to look at a LARGE set of clouds at once which *wouldn't* be possible to check by hand
#     
# 
# 2. i have calculated linewidths (is that == velocity dispersion) for the three clouds, but unsure whether those are the right values (related question -- best practice for efficiently finding an answer to a question when that question is 'what is this quantity that someone has most likely measured?'
# 

# notes from check in 2.17.23
# - go with the literature convention (0-360)
# 
# meet next (wednesday right after colloquium (12), but see if tues would work (heads up by monday afternoon))
# 
# 
# fwhm / 2.355 = velocity dispersion
# 
# -the velocity dispersion is always smaller than fwhm because it is the 1-sigma
# -velocity dispersion is related to sounds speed + other physical properties
# 
# velocity dispersion would be reported in CO papers for inidividual clouds
# 
# next step: 
# - calculate whether these clouds are virialized using the velocity dispersion
# - done by calculating the virial parameter -> Eg/Ek of the cloud
# - the virialized-ness of the cloud will be determined by ^ ratio
# - expectation: should be / close to virialized
# - intermediate step: calculate the masses of the clouds 
#     - define the cloud based on the CO emission & add up masses in that boundary
#     - same sort of calculation: column density
#     
# - do this for more clouds, determine a sample to work from 
# 
# - plot results against each other (e.g. mass-size relations across the clouds)
