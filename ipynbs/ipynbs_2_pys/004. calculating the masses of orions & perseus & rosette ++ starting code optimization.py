#!/usr/bin/env python
# coding: utf-8

# when cleaning up etc -- anytime i find myself repeating the same thing but for diff clouds ==> turn it into a function!!


# import datacube as dc

import numpy as np
import numpy.ma as ma
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
# cloud locations & distances from imara & burkhart page 2: https://arxiv.org/pdf/1609.04817.pdf
# (rounded to whole numbers)
# 
# cloud velocity spectra: imara & blitz https://arxiv.org/pdf/1103.3741.pdf
# 
# lewis, lada, dame: https://arxiv.org/pdf/2203.07480.pdf
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


distances = np.zeros(len(CLOUDS))
dist_erros = np.zeros_like(distances)
CLOUDS['distance'] = distances
CLOUDS['distance unc'] = dist_erros


# CLOUDS.at['Perseus','distance'] = 240 # l,l&d
CLOUDS.at['Perseus','distance'] = 240 # i&b
# CLOUDS.at['Perseus','distance unc'] = 13 # i&b

# CLOUDS.at['Ophiuchus','distance'] = 144 # l,l&d
CLOUDS.at['Ophiuchus','distance'] = 119 # i&b
# CLOUDS.at['Ophiuchus','distance unc'] = 6 # i&b

#rosette isn't in the l,l&d catalog
# CLOUDS.at['Rosette','distance'] = 1330 # i&b
# CLOUDS.at['Rosette','distance unc'] = 48 # i&b

# CLOUDS.at['Orion B','distance'] = 423 # l,l&d
CLOUDS.at['Orion B','distance'] = 398 # i&b
# CLOUDS.at['Orion B','distance unc'] = 12 # i&b

# CLOUDS.at['Orion A','distance'] = 432 # l,l&d
CLOUDS.at['Orion A','distance'] = 371 # i&b
# CLOUDS.at['Orion A','distance unc'] = 10 # i&b

CLOUDS.at['California', 'distance'] = 450 # l,l&d
# CLOUDS.at['California', 'distance'] = 450 # i&b
# CLOUDS.at['California', 'distance unc'] = 23 # i&b


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


galaxy = dc.unpack_fits_cube('COGAL_all_interp.fits')


galaxy_cube = galaxy[0]
gvels = galaxy[1]
glon = galaxy[2]
glat = galaxy[3]

#linear spacing
dv = gvels[1] - gvels[0]
dlon = glon[1] - glon[0]
dlat = glat[1] - glat[0]


galaxy_0th = np.nansum(galaxy_cube, axis = 0)


plt.figure(figsize=(15,3), dpi=200)

im = galaxy_0th
imx = glon
imy = glat


vmin,vmax = zlims(im)


plt.imshow(im, vmin = vmin, vmax = vmax, \
           extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower')

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = 'Temperature (K)')

plt.title('Galaxy 0th moment map - summed')
# plt.ylim(-100,600)

plt.show()


#convert 0-360 coords to +/-180 coords
def convertCoords (arrow, coords):
    output = np.zeros_like(coords)
    if arrow == 'to +/- 180':
        for i in range(len(coords)):
            output[i] = coords[i] - 360
    if arrow == 'to 360':
        for i in range(len(coords)):
            output[i] = coords + 360
    
    return output


#physical coordinates of clouds

#### ORION A ####
orionA_lon_phy = CLOUDS['longitude (l)']['Orion A']
orionA_lat_phy = CLOUDS['latitude (b)']['Orion A']

#### ORION B ####
orionB_lon_phy = CLOUDS['longitude (l)']['Orion B']
orionB_lat_phy = CLOUDS['latitude (b)']['Orion B']

#### PERSEUS ####
perseus_lon_phy = CLOUDS['longitude (l)']['Perseus']
perseus_lat_phy = CLOUDS['latitude (b)']['Perseus']

#### ROSETTE ####
rosette_lon_phy = CLOUDS['longitude (l)']['Rosette']
rosette_lat_phy = CLOUDS['latitude (b)']['Rosette']


### find a better way to do this checking

#then make sure the longitudes are in the right system
if np.all(np.array(orionA_lon_phy)>180):
    orionA_lon_phy = convertCoords('to +/- 180', orionA_lon_phy)

if np.all(np.array(orionB_lon_phy)>180):
    orionB_lon_phy = convertCoords('to +/- 180', orionB_lon_phy)

if np.all(np.array(perseus_lon_phy)>180):
    perseus_lon_phy = convertCoords('to +/- 180', perseus_lon_phy)

if np.all(np.array(rosette_lon_phy)>180):
    rosette_lon_phy = convertCoords('to +/- 180', rosette_lon_phy)


#index coordinates of clouds

#### ORION A ####
orionA_inGlon = dc.ind_phys('get index', glon, orionA_lon_phy)
orionA_inGlat = dc.ind_phys('get index', glat, orionA_lat_phy)

orionA_inGlon.sort()
orionA_inGlat.sort()


#### ORION B ####
orionB_inGlon = dc.ind_phys('get index', glon, orionB_lon_phy)
orionB_inGlat = dc.ind_phys('get index', glat, orionB_lat_phy)

orionB_inGlon.sort()
orionB_inGlat.sort()


#### PERSEUS ####
perseus_inGlon = dc.ind_phys('get index', glon, perseus_lon_phy)
perseus_inGlat = dc.ind_phys('get index', glat, perseus_lat_phy)

perseus_inGlon.sort()
perseus_inGlat.sort()


#### ROSETTE ####
rosette_inGlon = dc.ind_phys('get index', glon, rosette_lon_phy)
rosette_inGlat = dc.ind_phys('get index', glat, rosette_lat_phy)

rosette_inGlon.sort()
rosette_inGlat.sort()


orionA_lon_phy


#create masks to lift the clouds from the galaxy cube

#### ORION A ####
orionA_b_0, orionA_b_1 = orionA_inGlat
orionA_l_0, orionA_l_1 = orionA_inGlon

orionA_cube = galaxy_cube[:,orionA_b_0:orionA_b_1, orionA_l_0:orionA_l_1]


### ORION B ####
orionB_b_0, orionB_b_1 = orionB_inGlat
orionB_l_0, orionB_l_1 = orionB_inGlon

orionB_cube = galaxy_cube[:,orionB_b_0:orionB_b_1, orionB_l_0:orionB_l_1]


#### PERSEUS ####
perseus_b_0, perseus_b_1 = perseus_inGlat
perseus_l_0, perseus_l_1 = perseus_inGlon

perseus_cube = galaxy_cube[:,perseus_b_0:perseus_b_1, perseus_l_0:perseus_l_1]

#### ROSETTE ####
rosette_b_0, rosette_b_1 = rosette_inGlat
rosette_l_0, rosette_l_1 = rosette_inGlon

rosette_cube = galaxy_cube[:,rosette_b_0:rosette_b_1, rosette_l_0:rosette_l_1]


np.sqrt(np.nanmean(new0th_perseus))


np.sqrt(np.nanmean(perseus_cube))


fits.open('COGAL_interp.fits')[0].header


# the sigma^a that lewis et al report is *per velocity channel*
# & it tells you about the background
# 
# be on the lookout for: when the data has really large +/- numbers, so be on the look out for those & manually label them as nans --> plot a histogram of the cube & find those outlier numbers
# 
# figure out how to get a histogram of the cube without the cell taking a long time
# 
# explore different ways to calculate std (at different stages, eg pre-0th moment map vs post-0th moment, per channel) & make sure that i know what values are going in to the calculations for radius
# 
# visualize the boundaries of the cloud as another sanity check -> contours
# big picture: mass, virial parameter, surface density, empirical relations


#calculate the background around the star
#using a high cutoff (sigma=7) because cloud boundaries were hard with lower cutoffs
sigma = 7.

orA_mean, orA_med, orA_std = sigma_clipped_stats(orionA_cube, sigma=sigma)
orB_mean, orB_med, orB_std = sigma_clipped_stats(orionB_cube, sigma=sigma)
per_mean, per_med, per_std = sigma_clipped_stats(perseus_cube, sigma=sigma)
ros_mean, ros_med, ros_std = sigma_clipped_stats(rosette_cube, sigma=sigma)


#then create a mask of values above 
orA_sigma_mask = orionA_cube > orA_std
orB_sigma_mask = orionB_cube > orB_std
per_sigma_mask = perseus_cube > per_std
ros_sigma_mask = rosette_cube > ros_std


#use sigma mask to extract clouds
orionA_cube_noisy = orionA_cube
orionA_cube = orA_sigma_mask * orionA_cube_noisy


orionB_cube_noisy = orionB_cube
orionB_cube = orB_sigma_mask * orionB_cube_noisy

perseus_cube_noisy = perseus_cube
perseus_cube = per_sigma_mask * perseus_cube_noisy

rosette_cube_noisy = rosette_cube
rosette_cube = ros_sigma_mask * rosette_cube_noisy


#create 0th moment maps
orionA_0th = np.nansum(orionA_cube, axis=0)
orionB_0th = np.nansum(orionB_cube, axis=0)
perseus_0th = np.nansum(perseus_cube, axis=0)
rosette_0th = np.nansum(rosette_cube, axis=0)


###ugh need to figure out whether nansum or nanmean
#and how to translate integral(I_v*dv) into python

#then once i have that, can move to calculating column density etc


#calculate velocity spectra
orionA_spectrum = spectrumFromCube(orionA_cube)
orionB_spectrum = spectrumFromCube(orionB_cube)
perseus_spectrum = spectrumFromCube(perseus_cube)
rosette_spectrum = spectrumFromCube(rosette_cube)


#find the velocity dispersions
#steps are: identify a peak (defined as > 0.5 the max vel value)
#calculate the line fwhm
#calculate the 1 sigma velocity dispersion
sqrt2ln2x2 = 2*(np.sqrt(2*np.log(2)))

#### ORION A ####
orA_peak_vel = find_peaks(orionA_spectrum, height = 0.5*np.nanmax(orionA_spectrum))
orA_fwhm = peak_widths(orionA_spectrum,orA_peak_vel[0])[0]
orA_fwhm = orA_fwhm*dv #convert units: pixel -> km/s
orA_sigma = orA_fwhm/sqrt2ln2x2

#### ORION A ####
orB_peak_vel = find_peaks(orionB_spectrum, height = 0.5*np.nanmax(orionB_spectrum))
orB_fwhm = peak_widths(orionB_spectrum,orB_peak_vel[0])[0]
orB_fwhm = orB_fwhm*dv
orB_sigma = orB_fwhm/sqrt2ln2x2


#### ORION A ####
per_peak_vel = find_peaks(perseus_spectrum, height = 0.5*np.nanmax(perseus_spectrum))
per_fwhm = peak_widths(perseus_spectrum,per_peak_vel[0])[0]
per_fwhm = per_fwhm*dv
per_sigma = per_fwhm/sqrt2ln2x2

#### ORION A ####
ros_peak_vel = find_peaks(rosette_spectrum, height = 0.5*np.nanmax(rosette_spectrum))
ros_fwhm = peak_widths(rosette_spectrum,ros_peak_vel[0])[0]
ros_fwhm = ros_fwhm*dv
ros_sigma = ros_fwhm/sqrt2ln2x2


#the noise levels here are very big (compare to Table 1)
#this is giving a really high threshold -> consistent with very small R
orA_sigma, orB_sigma, per_sigma


# another sanity check - take the standard deviation of the data myself (sigma^a in table 1 in lewis, lada, dame)
# - look at the documentation sigma_clipped_stats

# 1 pc = 206265 AU
# 1 AU = 1.5 * 10^11 m


#pull out distances
to_cm = 206265*1.5*(10**13) #convert between pc -> cm

orionA_d = CLOUDS['distance']['Orion A']*to_cm #cm
orionB_d = CLOUDS['distance']['Orion B']*to_cm #cm
perseus_d = CLOUDS['distance']['Perseus']*to_cm #cm
# rosette_d = CLOUDS['distance']['Rosette']*to_cm #cm


perseus_d = 320*to_cm


orionA_d, orionB_d, perseus_d


# fits.open('COGAL_all_interp.fits')[0].header


#these are used in the next cell, apply to all the clouds
X_co = 1.97 * (10**20) #cm^-2/(K*km/s)

dA_sky = np.abs(dlon*dlat) #degrees squared
dA_sky = dA_sky*(4*(np.pi**2)/(360**2)) #radians squared

m_H2 = 2*1.67*(10**(-24)) #g
mu = 1.37 #mean molecular weight, account for helium


dA_sky


#calculate intensity from T & dv
orionA_I = orionA_0th*dv # K*km/s
orionB_I = orionB_0th*dv # K*km/s
perseus_I = perseus_0th*dv # K*km/s
# rosette_I = rosette_0th*dv # K*km/s

#calculate column density from I & X_co
orionA_N = orionA_I*X_co # cm^-2
orionB_N = orionB_I*X_co # cm^-2
perseus_N = perseus_I*X_co # cm^-2
# rosette_N = rosette_I*X_co # cm^-2

#calculate the angular size of each pixel
orionA_dA = dA_sky*(orionA_d**2) #cm^2
orionB_dA = dA_sky*(orionB_d**2) #cm^2
perseus_dA = dA_sky*(perseus_d**2) #cm^2
# rosette_dA = dA_sky*(rosette_d**2) #cm^2

#calculate the mass enclosed at each pixel
orionA_dm = orionA_dA * m_H2 * mu * orionA_N # g
orionA_M = np.nansum(orionA_dm) #g

orionB_dm = orionB_dA * m_H2 * mu * orionB_N # g
orionB_M = np.nansum(orionB_dm)

perseus_dm = perseus_dA * m_H2 * mu * perseus_N # g
perseus_M = np.nansum(perseus_dm)

# rosette_dm = rosette_dA * m_H2 * rosette_N # g
# rosette_M = np.nansum(rosette_dm) #g


#calculate radii from areas to compare to l,l&d table 2
orionA_A = np.nansum(orionA_dA) #square cm
orionA_A = orionA_A/(to_cm**2) #square pc
orionA_R = np.sqrt(orionA_A/np.pi) #pc

orionB_A = np.nansum(orionB_dA) #square cm
orionB_A = orionB_A/(to_cm**2) #square pc
orionB_R = np.sqrt(orionB_A/np.pi) #pc

perseus_A = np.nansum(perseus_dA) #square cm
perseus_A = perseus_A/(to_cm**2) #square pc
perseus_R = np.sqrt(perseus_A/np.pi) #pc

# orionB_dA = dA_sky*(orionB_d**2) #square cm
# perseus_dA = dA_sky*(perseus_d**2) #square cm


orionA_M, orionB_M, perseus_M


orionA_R, orionB_R, perseus_R


perseus_M/(2*10**36)


orionA_M/(2*10**36)


# questions for 2.28 checkin
# 
# -- for 0th moment, intensity etc (stuff that integrates down the velocity axis), do you sum or mean?
#     (i thought you sum, but then color bar is orders of magnitdue greater than when you mean / what it's 'supposed to be')
#     
# 
# -- noticed a discrepancy between distance values reported literature
#     - lewis, lada, dame : page 5, table 1
#     - imara & blitz : page 10, table 1
#     - imara & burkhart : page 2, table 1
#     & wondering if this is why the mass values i calculate
#     
# -- even when using the distance values from lewis, lada & dame, the masses i calculate are different from the ones they report on page 7, table 2
#     - Orion A: 80.1 (10^3 Msun) --- 0.94x l,l&d
#     - Orion B: 81.7 (10^3 Msun) --- 1.26x l,l&d
#     - Perseus: 32.7 (10^3 Msun) --- 1.71x l,l&d
#     
#     the mass inconsistencies don't seem systematic
#     
# -- calculating radius values (as a sanity check) also gives really small numbers, so something is missing from my calculation?
#     - Orion A: 0.53 pc --- 2.3% of l,l&d
#     - Orion B: 0.52 pc --- 2.4% of l,l&d
#     - Perseus: 0.29 pc --- 2.2% of l,l&d
#     
#     radius inconsistencies *do* seem systematic ... 


#### ORION A 0TH MOMENT ####

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


# #### ORION A VELOCITY SPECTRUM ####
# # plt.plot(gvels,orionA_spectrum, color = 'gray')
# plt.plot(orionA_spectrum, color = 'gray')
# plt.title('Orion A velocity spectrum')

# plt.ylabel('CO Temperature (K)') #should be temperature (K)
# plt.xlabel('Velocity (km/s)')

# plt.show()


#### START THINKING ABOUT DEFINING A METHOD THAT TAKES IN THE NAME OF THE
#### STAR & THE GALAXY CUBE, THEN PULLS IT OUT, DOES 0TH MOMENT, AND OUTPUTS
#THE CUBE FOR THAT STAR, THE 0TH MOMENT, AND THE VELOCITY SPECTRUM


#### PLOTTING THINGS ####


#### ORION B 0TH MOMENT ####
# plt.figure(figsize=(10,10), dpi=100)

# im = orionB_0th
# imx = glon[orionB_l_0:orionB_l_1]+360 #put lon coords back in 360 form
# imy = glat[orionB_b_0:orionB_b_1]


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]])

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

# plt.title('Orion B 0th moment')

# plt.show()


# #### ORION B VELOCITY SPECTRUM ####
# plt.plot(gvels,orionB_spectrum, color = 'gray')
# plt.title('Orion B velocity spectrum')

# plt.ylabel('CO Temperature (K)') #should be temperature (K)
# plt.xlabel('Velocity (km/s)')

# plt.show()


plt.pcolormesh(per_sigma_mask)


#### PERSEUS 0TH MOMENT ####
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


# #### PERSEUS VELOCITY SPECTRUM ####
# plt.plot(gvels,perseus_spectrum, color = 'gray')


# plt.title('Perseus velocity spectrum')

# plt.ylabel('CO Temperature (K)') #should be temperature (K)
# plt.xlabel('Velocity (km/s)')

# plt.show()


# #### ROSETTE 0TH MOMENT ####
# plt.figure(figsize=(10,10), dpi=100)

# im = rosette_0th
# imx = glon[rosette_l_0:rosette_l_1]+360
# imy = glat[rosette_b_0:rosette_b_1]


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]])

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

# plt.title('Rosette 0th moment')

# plt.show()


# #### ROSETTE VELOCITY SPECTRUM ####
# plt.plot(gvels,rosette_spectrum, color = 'gray')


# plt.title('Rosette velocity spectrum')

# plt.ylabel('CO Temperature (K)') #should be temperature (K)
# plt.xlabel('Velocity (km/s)')

# plt.show()


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
