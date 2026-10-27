#!/usr/bin/env python
# coding: utf-8


import datacube as dc
import clouds as cld

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits


# data: https://lweb.cfa.harvard.edu/rtdc/CO/CompositeSurveys/
# 
# cloud locations & distances from imara & burkhart page 2: https://arxiv.org/pdf/1609.04817.pdf
# (rounded to whole numbers)
# 
# cloud velocity spectra: imara & blitz https://arxiv.org/pdf/1103.3741.pdf
# 
# lewis, lada, dame: https://arxiv.org/pdf/2203.07480.pdf
# 
# dame, hartmann, thaddeus re the survey: https://iopscience.iop.org/article/10.1086/318388/pdf


### THINK ABOUT TURNING THIS INTO AN EXTERNAL CSV OR STH I JUST READ IN

#create a dataframe to contain cloud locations
columns = ['longitude (l)', 'latitude (b)', 'central velocity (km/s)']
CLOUDS = pd.DataFrame(columns = columns)

#then set the relevant details that i know so far
CLOUDS.at['Perseus','longitude (l)'] = np.array((155,162))
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


CLOUDS.at['Perseus','distance'] = 240 # l,l&d
# CLOUDS.at['Perseus','distance'] = 240 # i&b
# CLOUDS.at['Perseus','distance unc'] = 13 # i&b

CLOUDS.at['Ophiuchus','distance'] = 144 # l,l&d
# CLOUDS.at['Ophiuchus','distance'] = 119 # i&b
# CLOUDS.at['Ophiuchus','distance unc'] = 6 # i&b

#rosette isn't in the l,l&d catalog
# CLOUDS.at['Rosette','distance'] = 1330 # i&b
# CLOUDS.at['Rosette','distance unc'] = 48 # i&b

CLOUDS.at['Orion B','distance'] = 423 # l,l&d
# CLOUDS.at['Orion B','distance'] = 398 # i&b
# CLOUDS.at['Orion B','distance unc'] = 12 # i&b

CLOUDS.at['Orion A','distance'] = 432 # l,l&d
# CLOUDS.at['Orion A','distance'] = 371 # i&b
# CLOUDS.at['Orion A','distance unc'] = 10 # i&b

CLOUDS.at['California', 'distance'] = 450 # l,l&d
# CLOUDS.at['California', 'distance'] = 450 # i&b
# CLOUDS.at['California', 'distance unc'] = 23 # i&b


# CLOUDS


# CLOUDS.to_csv('cloooooods.csv')


# reread = pd.read_csv('cloooooods.csv')


# reread["longitude (l)"]


# k =np.array((reread["longitude (l)"][0]))


cmap0 = dc.cmap_grey0('viridis')


galaxy = dc.unpack_fits_cube('COGAL_all_interp.fits')


galaxy_cube = galaxy[0]
gvels = galaxy[1]
glat = galaxy[2]
glon = galaxy[3]

#assume arrays have linear spacing
dv = gvels[1] - gvels[0]
dlat = glat[1] - glat[0]
dlon = glon[1] - glon[0]


galaxy_0th = np.nansum(galaxy_cube, axis = 0)


#extract clouds outputs in the same form as unpack_fits_cube

#----- Orion A --------
orionA = cld.extract_cloud('Orion A', CLOUDS, galaxy)
orionA_cube = orionA[0]
orionA_vels = orionA[1]
orionA_lats = orionA[2]
orionA_lons = orionA[3]


#----- Orion B --------

orionB = cld.extract_cloud('Orion B', CLOUDS, galaxy)
orionB_cube = orionB[0]
orionB_vels = orionB[1]
orionB_lats = orionB[2]
orionB_lons = orionB[3]


#----- Perseus --------

perseus = cld.extract_cloud('Perseus', CLOUDS, galaxy)
perseus_cube = perseus[0]
perseus_vels = perseus[1]
perseus_lats = perseus[2]
perseus_lons = perseus[3]


for i in range(len(perseus_cube)):
    if np.all(perseus_cube[i]==0):
        print(i)


# fig, axs = plt.subplots(1,3,figsize=(15,10), dpi=100)

# kwargs = {'cmap':cmap0}

# vel_slices = [3,12,23]
# cube = orionB_cube


# for (ax,vel) in zip(axs,vel_slices):
#     frame = ax.pcolormesh(cube[vel], **kwargs)
#     fig.colorbar(frame, ax=ax)
    
#     plt.suptitle('slices after rms cutoff')

# plt.show()
    
    
# #old way of printing things out, very repetitive
# # noisy = axs[0,0].pcolormesh(cube[5],**kwargs)
# # axs[0,0].set_title('5th slice - last noisy')
# # plt.colorbar(noisy, ax=axs[0,0])

# # first_data = axs[0,1].pcolormesh(cube[6],**kwargs)
# # axs[0,1].set_title('first slice w data > 3*rms_cube')
# # plt.colorbar(first_data, ax=axs[0,1])

# # most = axs[1,0].pcolormesh(cube[10],**kwargs)
# # axs[1,0].set_title('10th slice - lots of temp')
# # plt.colorbar(most, ax=axs[1,0])

# # most_2 = axs[1,1].pcolormesh(cube[11],**kwargs)
# # axs[1,1].set_title('11th slice - also much temp')
# # plt.colorbar(most_2, ax=axs[1,1])

# # last_data = axs[2,0].pcolormesh(cube[18],**kwargs)
# # axs[2,0].set_title('last slice w data > 3*rms_cube')
# # plt.colorbar(first_data, ax=axs[2,0])

# # noisy_2 = axs[2,1].pcolormesh(cube[19],**kwargs)
# # axs[2,1].set_title('19 - first noisy')
# # plt.colorbar(noisy_2, ax=axs[2,1])

# # plt.suptitle('slices after rms cutoff')

# # plt.show()


#create 0th moment maps
#defined as sum(T*dv), unit is K*km/s
orionA_0th = np.nansum(orionA_cube*dv, axis=0)
orionB_0th = np.nansum(orionB_cube*dv, axis=0)
perseus_0th = np.nansum(perseus_cube*dv, axis=0)
# rosette_0th = np.nansum(rosette_cube, axis=0) #old method


#calculate velocity spectra
orionA_spectrum = dc.spectrum1D(orionA_cube)
orionB_spectrum = dc.spectrum1D(orionB_cube)
perseus_spectrum = dc.spectrum1D(perseus_cube)
# rosette_spectrum = dc.spectrum1D(rosette_cube)


#calculate velocity dispersion
orA_sigma = dv * dc.dispersion_1sigma(orionA_spectrum)
orB_sigma = dv * dc.dispersion_1sigma(orionB_spectrum)
per_sigma = dv * dc.dispersion_1sigma(orionB_spectrum)
# ros_sigma = dv * dc.dispersion_1sigma(rosette_spectrum)


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


# fits.open('COGAL_all_interp.fits')[0].header


#these are used in the next cell, apply to all the clouds
X_co = 1.97 * (10**20) #cm^-2/(K*km/s)
to_rad = 4*(np.pi**2)/(360**2)

dA_sky = np.abs(dlon*to_rad * dlat*to_rad) #radians squared

#new way of defining the pixel area - from lewis, lada, dame
# dA_sky = 2* dlon*to_rad * np.sin(to_rad*dlat*0.5)


# m_H2 = 2*1.67*(10**(-24)) #g
mp = 1.67*(10**(-24)) #g
mu = 1.37 #mean molecular weight, account for helium


# what lewis lada & dame have -> 
# - M = sum of { Sigma * da }, Sigma = surface density
# - Sigma = mu * mp * beta * A
# - beta = N/A
# - Sigma = mu * mp * N
# - M = sum of { mu * mp * N * da }
# 
# which ultimately is also the things i am multiplying together 
# because i have
# - M = sum of dm
# - dm = da * mp * mu * N
# - M = sum of { da * mp * mu * N }
#     

# future steps from 03.07 check in
# switch around order of add, multiply (when calculating intensity)
# & compare to values in papers that calculate mass using CO 


#calculate column density from W_co (= 0th moment map) & X_co
orionA_N = orionA_0th*X_co # cm^-2
orionB_N = orionB_0th*X_co # cm^-2
perseus_N = perseus_0th*X_co # cm^-2
# rosette_N = rosette_I*X_co # cm^-2

#calculate the angular size of each pixel
# orionA_dA = dA_sky*(orionA_d**2) #cm^2
# orionB_dA = dA_sky*(orionB_d**2) #cm^2
# perseus_dA = dA_sky*(perseus_d**2) #cm^2
# # rosette_dA = dA_sky*(rosette_d**2) #cm^2
orionA_dA = np.ones_like(orionA_0th) * dA_sky * (orionA_d**2)# * np.cos(to_rad * orionA_lats[0]) #cm^2
orionB_dA = np.ones_like(orionB_0th) * dA_sky*(orionB_d**2) #* np.cos(to_rad * orionB_lats[0]) #cm^2
perseus_dA = np.ones_like(perseus_0th) * dA_sky*(perseus_d**2)# * np.cos(to_rad * perseus_lats[0]) #cm^2


#calculate the mass enclosed at each pixel
orionA_dm = orionA_dA * mp * mu * orionA_N # g
orionA_M = np.nansum(orionA_dm) #g

orionB_dm = orionB_dA * mp * mu * orionB_N # g
orionB_M = np.nansum(orionB_dm)

perseus_dm = perseus_dA * mp * mu * perseus_N # g
perseus_M = np.nansum(perseus_dm)

# rosette_dm = rosette_dA * m_H2 * rosette_N # g
# rosette_M = np.nansum(rosette_dm) #g


# from https://arxiv.org/pdf/astro-ph/0411089.pdf:
# The Orion A cloud, the largest in the Orion-Monoceros complex, is located ∼ 19◦ below the Galactic plane at a distance of ∼ 450 pc from the Sun. It has a surface area of ∼ 31.5 deg2 (∼ 2200 pc2) 


#calculate radii from areas to compare to l,l&d table 2
orionA_A = np.nansum(orionA_dA) #square cm
# orionA_A = orionA_A/(to_cm**2) #square pc
orionA_R = np.sqrt(orionA_A/np.pi) #pc

orionB_A = np.nansum(orionB_dA) #square cm
orionB_A = orionB_A/(to_cm**2) #square pc
orionB_R = np.sqrt(orionB_A/np.pi) #pc

perseus_A = np.nansum(perseus_dA) #square cm
perseus_A = perseus_A/(to_cm**2) #square pc
perseus_R = np.sqrt(perseus_A/np.pi) #pc

# orionB_dA = dA_sky*(orionB_d**2) #square cm
# perseus_dA = dA_sky*(perseus_d**2) #square cm


### 0th moment maps ###


#!! need a better way of determining when lon coords need to be adjusted by +360 or not

fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower'}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

#-----orion A------------
im = orionA_0th
imx = orionA_lons +360 #put lon coords back in 360 form
imy = orionA_lats

vmin, vmax = zlims(im)

orA = axs[0,0].imshow(im, vmin = vmin, vmax = vmax,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,0].set_title("Orion A")
fig.colorbar(orA, ax=axs[0,0], **kw_colorb)


#-----orion B------------
im = orionB_0th
imx = orionB_lons +360 #put lon coords back in 360 form
imy = orionB_lats

vmin, vmax = zlims(im)

orB = axs[0,1].imshow(im, vmin = vmin, vmax = vmax,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,1].set_title("Orion B")
fig.colorbar(orB, ax=axs[0,1], **kw_colorb)


#-----perseus ------------
im = perseus_0th
imx = perseus_lons
imy = perseus_lats

vmin, vmax = zlims(im)

per = axs[1,0].imshow(im, vmin = vmin, vmax = vmax,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[1,0].set_title("Perseus")
fig.colorbar(per, ax=axs[1,0], **kw_colorb)


for ax in axs.flat:
    ax.set(xlabel='Galactic Longitude (deg)',\
          ylabel = 'Galactic Latitude (deg)')


plt.suptitle('0th moment maps')
plt.show()


# ### velocity spectra ###
# fig, axs = plt.subplots(2,2, figsize=(10,10), dpi=200)

# #----- Orion A ------------
# axs[0,0].plot(orionA_vels,orionA_spectrum, color = 'gray')
# axs[0,0].set_title('Orion A')

# #----- Orion B ------------
# axs[0,1].plot(orionB_vels,orionB_spectrum, color = 'gray')
# axs[0,1].set_title('Orion B')


# #----- Perseus ------------
# axs[1,0].plot(perseus_vels,perseus_spectrum, color = 'gray')
# axs[1,0].set_title('Perseus')


# for ax in axs.flat:
#     ax.set(xlabel='Velocity (km/s)',\
#           ylabel='CO Temperature (K)')

    
# plt.suptitle('Velocity Spectra')
# plt.show()


# plt.figure(figsize=(15,3), dpi=200)

# im = galaxy_0th
# imx = glon
# imy = glat


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin = vmin, vmax = vmax, \
#            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower')

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = 'Temperature (K)')

# plt.title('Galaxy 0th moment map - summed')
# # plt.ylim(-100,600)

# plt.show()


#### PLOTTING THINGS ####


# #### ORION A 0TH MOMENT ####

# plt.figure(figsize=(12,5), dpi=100)

# im = orionA_0th

# imx = glon[orionA_lons[0]:orionA_lons[1]]+360 #put lon coords back in 360 form
# imy = glat[orionA_lats[0]:orionA_lats[1]]


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]])

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

# plt.title('Orion A 0th moment')

# plt.show()


# #### ORION A VELOCITY SPECTRUM ####
# # plt.plot(gvels,orionA_spectrum, color = 'gray')
# plt.plot(orionA_spectrum, color = 'gray')
# plt.title('Orion A velocity spectrum')

# plt.ylabel('CO Temperature (K)') #should be temperature (K)
# plt.xlabel('Velocity (km/s)')

# plt.show()


# ### ORION B 0TH MOMENT ####
# plt.figure(figsize=(10,10), dpi=100)

# im = orionB_0th
# imx = glon[orionB_lons[0]:orionB_lons[1]]+360 #put lon coords back in 360 form
# imy = glat[orionB_lats[0]:orionB_lats[1]]


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


# #### PERSEUS 0TH MOMENT ####
# plt.figure(figsize=(10,10), dpi=100)

# im = perseus_0th
# imx = glon[perseus_lons[0]:perseus_lons[1]]
# imy = glat[perseus_lats[0]:perseus_lats[1]]


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin = vmin, vmax = vmax,origin = 'lower',\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]])

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='horizontal', label = 'Temperature (K)')

# plt.title('Perseus 0th moment')

# plt.show()


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
