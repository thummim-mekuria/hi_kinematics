#!/usr/bin/env python
# coding: utf-8


import datacube as dc
import clouds as cld

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c

from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits


### THINK ABOUT TURNING THIS INTO AN EXTERNAL CSV OR STH I JUST READ IN

#create a dataframe to contain cloud locations
columns = ['longitude (l)', 'latitude (b)', 'l, b unit', 'distance', 'd unit']
CLOUDS = pd.DataFrame(columns = columns)

#then set the relevant details that i know so far
#updated values (from past notebooks) from lewis, lada & dame

# CLOUDS.at['Perseus','longitude (l)'] = (154.9,165.1)
# CLOUDS.at['Perseus','latitude (b)'] = (-25.1,-14.9)
CLOUDS.at['Perseus','longitude (l)'] = np.array((155,162))
CLOUDS.at['Perseus','latitude (b)'] = (-27,-15)
CLOUDS.at['Perseus','l, b unit'] = u.deg

CLOUDS.at['Perseus','distance'] = 240 
CLOUDS.at['Perseus','d unit'] = u.pc 


CLOUDS.at['Ophiuchus','longitude (l)'] = (11.1,349.9)
CLOUDS.at['Ophiuchus','latitude (b)'] = (11.9,25.1)
CLOUDS.at['Ophiuchus','l, b unit'] = u.deg

CLOUDS.at['Ophiuchus','distance'] = 144 
CLOUDS.at['Ophiuchus','d unit'] = u.pc 


CLOUDS.at['Orion B', 'longitude (l)'] = (209.8,202.2)
CLOUDS.at['Orion B', 'latitude (b)'] = (-17.1, -9.6)
CLOUDS.at['Orion B','l, b unit'] = u.deg

CLOUDS.at['Orion B','distance'] = 423 
CLOUDS.at['Orion B','d unit'] = u.pc 


CLOUDS.at['Orion A', 'longitude (l)'] = (217.1,202.9)
CLOUDS.at['Orion A', 'latitude (b)'] = (-21.1, -16.9)
CLOUDS.at['Orion A','l, b unit'] = u.deg

CLOUDS.at['Orion A','distance'] = 432 
CLOUDS.at['Orion A','d unit'] = u.pc 


CLOUDS.at['California', 'longitude (l)'] = (169.6,154.9)
CLOUDS.at['California', 'latitude (b)'] = (-13.1,-4.9)
CLOUDS.at['California','l, b unit'] = u.deg

CLOUDS.at['California', 'distance'] = 450 
CLOUDS.at['California','d unit'] = u.pc 


cmap0 = dc.cmap_grey0('viridis')


galaxy = dc.unpack_fits_cube('COGAL_all_interp.fits')
galaxy_cube = galaxy[0]
gvels = galaxy[1]
glat = galaxy[2]
glon = galaxy[3]
gunits = galaxy[4]

#assume arrays have linear spacing
dv = gvels[1] - gvels[0]
dlat = glat[1] - glat[0]
dlon = glon[1] - glon[0]

galaxy_0th = np.nansum(galaxy_cube*dv, axis = 0)


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

#----- California --------

california = cld.extract_cloud('California', CLOUDS, galaxy)
california_cube = california[0]
california_vels = california[1]
california_lats = california[2]
california_lons = california[3]

# #commenting out ophiuchus because the longitudes go across 0 from +5 to -5
# #& it is causing problems for the coords conversion
# #----- Ophiuchus --------

# ophiuchus = cld.extract_cloud('Ophiuchus', CLOUDS, galaxy)
# ophiuchus_cube = ophiuchus[0]
# ophiuchus_vels = ophiuchus[1]
# ophiuchus_lats = ophiuchus[2]
# ophiuchus_lons = ophiuchus[3]


#create 0th moment maps
#defined as sum(T*dv), unit is K*km/s
orionA_0th = np.nansum(orionA_cube*dv, axis=0)
orionB_0th = np.nansum(orionB_cube*dv, axis=0)
perseus_0th = np.nansum(perseus_cube*dv, axis=0)
california_0th = np.nansum(california_cube*dv, axis=0)


#calculate velocity spectra
orionA_spectrum = dc.spectrum1D(orionA_cube)
orionB_spectrum = dc.spectrum1D(orionB_cube)
perseus_spectrum = dc.spectrum1D(perseus_cube)
california_spectrum = dc.spectrum1D(california_cube)


#calculate velocity dispersion
orA_sigma = dv * gunits['vels'] * dc.dispersion_1sigma(orionA_spectrum)
orB_sigma = dv * gunits['vels'] * dc.dispersion_1sigma(orionB_spectrum)
per_sigma = dv * gunits['vels'] * dc.dispersion_1sigma(perseus_spectrum)
cal_sigma = dv * gunits['vels'] * dc.dispersion_1sigma(california_spectrum)


#pull out distances
orionA_d = (CLOUDS['distance']['Orion A'] * CLOUDS['d unit']['Orion A']).to(u.cm) #pc
orionB_d = (CLOUDS['distance']['Orion B'] * CLOUDS['d unit']['Perseus']).to(u.cm) #pc
perseus_d = (CLOUDS['distance']['Perseus'] * CLOUDS['d unit']['Perseus']).to(u.cm) #pc
california_d = (CLOUDS['distance']['California'] * CLOUDS['d unit']['California']).to(u.cm) #pc


#these are used in calculating mass & apply to all the clouds
X_co = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))

mp = c.m_p.to(u.g)
mH2 = 2*mp
mu = 1.37 #mean molecular weight, account for helium

dA_sky = np.abs(dlon*gunits['lons'].to(u.rad) * dlat*gunits['lats'].to(u.rad)) #radians squared

#new way of defining the pixel area - from lewis, lada, dame
# dA_sky = 2* dlon* gunits['lons'].to(u.rad)\
#         * np.sin(0.5* dlat* gunits['lats'].to(u.rad))


gunits['lons'].to(u.rad)


#calculate column density from W_co (== 0th moment map) & X_co
orionA_W = orionA_0th*u.K*u.km/u.s #add units to 0th moment map
orionB_W = orionB_0th*u.K*u.km/u.s 
perseus_W = perseus_0th*u.K*u.km/u.s
california_W = california_0th*u.K*u.km/u.s


orionA_N = orionA_W * X_co # cm^-2
orionB_N = orionB_W * X_co # cm^-2
perseus_N = perseus_W * X_co # cm^-2
california_N = california_W * X_co # cm^-2

#calculate the angular size of each pixel
orionA_dA = np.ones_like(orionA_0th) * dA_sky * orionA_d**2 #* np.cos(orionA_lats[0]* gunits['lats'].to(u.rad)) #cm^2
orionB_dA = np.ones_like(orionB_0th) * dA_sky * orionB_d**2 #* np.cos(orionB_lats[0]* gunits['lats'].to(u.rad)) #cm^2
perseus_dA = np.ones_like(perseus_0th) * dA_sky * perseus_d**2 #* np.cos(perseus_lats[0]* gunits['lats'].to(u.rad)) #cm^2
california_dA = np.ones_like(california_0th) * dA_sky * california_d**2 #* np.cos(california_lats[0]* gunits['lats'].to(u.rad)) #cm^2


#calculate the mass enclosed at each pixel
orionA_dm = orionA_dA * mH2 * mu * orionA_N # g
orionA_M = np.nansum(orionA_dm) #g

orionB_dm = orionB_dA * mH2 * mu * orionB_N # g
orionB_M = np.nansum(orionB_dm)

perseus_dm = perseus_dA * mH2 * mu * perseus_N # g
perseus_M = np.nansum(perseus_dm)

california_dm = california_dA * mH2 * mu * california_N # g
california_M = np.nansum(california_dm)


# from https://arxiv.org/pdf/astro-ph/0411089.pdf:
# The Orion A cloud, the largest in the Orion-Monoceros complex, is located ∼ 19◦ below the Galactic plane at a distance of ∼ 450 pc from the Sun. It has a surface area of ∼ 31.5 deg2 (∼ 2200 pc2) 


#calculate radii from areas to compare to l,l&d table 2
orionA_A = np.nansum(orionA_dA) #square cm
orionA_A = orionA_A.to(u.cm**2) #square pc
orionA_R = np.sqrt(abs(orionA_A/np.pi)) #pc

orionB_A = np.nansum(orionB_dA) #square cm
orionB_A = orionB_A.to(u.cm**2) #square pc
orionB_R = np.sqrt(abs(orionB_A/np.pi)) #pc

perseus_A = np.nansum(perseus_dA) #square cm
perseus_A = perseus_A.to(u.cm**2) #square pc
perseus_R = np.sqrt(abs(perseus_A/np.pi)) #pc

california_A = np.nansum(california_dA) #square cm
california_A = california_A.to(u.cm**2) #square pc
california_R = np.sqrt(abs(california_A/np.pi)) #pc


orionA_A


# notes for check-in 3.14
# 
# - got pretty close for mass, but still overestimating radius
#     1) factor differences in mass
#     - orion A: 1.04 * lld value
#     - orion B: 1.02 * lld value
#     - perseus: 1.89 * lld value (so something different is happening with perseus)
#     - california: 1.07 * lld value
# 
# 
#     2) factor differences in radius
#     - orion A: 1.44 * lld value
#     - orion B: 1.45 * lld value
#     - perseus: 1.85 * lld value (again, something different is happening here)
#         look back to see if there are any typos, change the bounds of how perseus is defined
#     - california: 1.40 * lld value
# 
# - think that this could be because a cutoff of > 3* rms accepts more area than a cutoff of W_co > 1.5 Kkm/s
# ...think of ways to implement that
# - 1) try: doing cut-off
# - 2) choose velocities to integrate over, do 0th moment, find the noise of that map & then apply a 3 sigma cut-off based on that
# 
# 
# - still having a hard time finding papers with cloud mass calculated from co
# ...been searching along the lines of 'mass giant molecular cloud orion A' in arxiv -- is there a more
# efficient way to do this?
#     -- more efficient way -> look at references
# 
# 
# 
# --> start plotting correlations, virial parameter, velocity dispersions, 
# 
# start thinking about: what are bits of information we want to know about the HI clouds -- what do papers have to say?
# 


# #comparing to lld values
# (orionA_M.to(u.kg)/(1e3*c.M_sun)),(orionA_M.to(u.kg)/(1e3*c.M_sun))/85.2,\
# (orionB_M.to(u.kg)/(1e3*c.M_sun)),(orionB_M.to(u.kg)/(1e3*c.M_sun))/65.0, \
# (perseus_M.to(u.kg)/(1e3*c.M_sun)),(perseus_M.to(u.kg)/(1e3*c.M_sun))/19.1,\
# (california_M.to(u.kg)/(1e3*c.M_sun)),(california_M.to(u.kg)/(1e3*c.M_sun))/139.1


# orionA_R.to(u.pc)/23.,\
# orionB_R.to(u.pc)/21.6,\
# perseus_R.to(u.pc)/13.1,\
# california_R.to(u.pc)/35.


# orionA_A.to(u.pc**2), orionB_A.to(u.pc**2), perseus_A.to(u.pc**2)


### 0th moment maps ###


#!! need a better way of determining when lon coords need to be adjusted by +360 or not

fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

# #-----orion A------------
# im = orionA_0th
# imx = orionA_lons +360 #put lon coords back in 360 form
# imy = orionA_lats

# vmin, vmax = zlims(im)

# orA = axs[0,0].imshow(im,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[0,0].set_title("Orion A")

# wco_thresh = axs[0,0].contour(imx,imy,im,levels=[1.5], colors='r')

# orA_cbar = fig.colorbar(orA, ax=axs[0,0], **kw_colorb)
# orA_cbar.add_lines(wco_thresh)


# #-----orion B------------
# im = orionB_0th
# imx = orionB_lons +360 #put lon coords back in 360 form
# imy = orionB_lats

# vmin, vmax = zlims(im)

# orB = axs[0,1].imshow(im,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[0,1].set_title("Orion B")

# wco_thresh = axs[0,1].contour(imx,imy,im,levels=[1.5], colors='r')


# orB_cbar = fig.colorbar(orB, ax=axs[0,1], **kw_colorb)
# orB_cbar.add_lines(wco_thresh)


#-----perseus------------
im = perseus_0th
imx = perseus_lons
imy = perseus_lats

vmin, vmax = zlims(im)

per = axs[1,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[1,0].set_title("Perseus")

wco_thresh = axs[1,0].contour(imx,imy,im,levels=[10], colors='r')

per_cbar = fig.colorbar(per, ax=axs[1,0], **kw_colorb)
per_cbar.add_lines(wco_thresh)


#-----california------------
im = california_0th
imx = california_lons
imy = california_lats

vmin, vmax = zlims(im)

per = axs[1,1].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[1,1].set_title("California")

wco_thresh = axs[1,1].contour(imx,imy,im,levels=[10], colors='r')

cal_cbar = fig.colorbar(per, ax=axs[1,1], **kw_colorb)
cal_cbar.add_lines(wco_thresh)


for ax in axs.flat:
    ax.set(xlabel='Galactic Longitude (deg)',\
          ylabel = 'Galactic Latitude (deg)')


plt.suptitle('0th moment maps')
plt.show()


### 0th moment maps ###


#!! need a better way of determining when lon coords need to be adjusted by +360 or not

fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

#-----orion A------------
im = orionA_0th
imx = orionA_lons +360 #put lon coords back in 360 form
imy = orionA_lats

vmin, vmax = zlims(im)

orA = axs[0,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,0].set_title("Orion A")

wco_thresh = axs[0,0].contour(imx,imy,im,levels=[1.5], colors='r')

orA_cbar = fig.colorbar(orA, ax=axs[0,0], **kw_colorb)
orA_cbar.add_lines(wco_thresh)


#-----orion B------------
im = orionB_0th
imx = orionB_lons +360 #put lon coords back in 360 form
imy = orionB_lats

vmin, vmax = zlims(im)

orB = axs[0,1].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,1].set_title("Orion B")

wco_thresh = axs[0,1].contour(imx,imy,im,levels=[1.5], colors='r')


orB_cbar = fig.colorbar(orB, ax=axs[0,1], **kw_colorb)
orB_cbar.add_lines(wco_thresh)


#-----perseus------------
im = perseus_0th
imx = perseus_lons
imy = perseus_lats

vmin, vmax = zlims(im)

per = axs[1,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[1,0].set_title("Perseus")

wco_thresh = axs[1,0].contour(imx,imy,im,levels=[1.5], colors='r')

per_cbar = fig.colorbar(per, ax=axs[1,0], **kw_colorb)
per_cbar.add_lines(wco_thresh)


#-----california------------
im = california_0th
imx = california_lons
imy = california_lats

vmin, vmax = zlims(im)

per = axs[1,1].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[1,1].set_title("California")

wco_thresh = axs[1,1].contour(imx,imy,im,levels=[1.5], colors='r')

cal_cbar = fig.colorbar(per, ax=axs[1,1], **kw_colorb)
cal_cbar.add_lines(wco_thresh)


for ax in axs.flat:
    ax.set(xlabel='Galactic Longitude (deg)',\
          ylabel = 'Galactic Latitude (deg)')


plt.suptitle('0th moment maps')
plt.show()


# ### column density maps ###


# #!! need a better way of determining when lon coords need to be adjusted by +360 or not

# fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
# kw_imshow = {'cmap':cmap0, 'origin':'lower'}
# kw_colorb = {'label':r'$N_{H_2}$ ($10^{21}$ $cm^{-2}$)', 'orientation':'horizontal'}

# #-----orion A------------
# im = orionA_N
# imx = orionA_lons +360 #put lon coords back in 360 form
# imy = orionA_lats

# vmin, vmax = zlims(im)

# orA = axs[0,0].imshow(im, vmin = vmin, vmax = vmax,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[0,0].set_title("Orion A")
# fig.colorbar(orA, ax=axs[0,0], **kw_colorb)


# #-----orion B------------
# im = orionB_N
# imx = orionB_lons +360 #put lon coords back in 360 form
# imy = orionB_lats

# vmin, vmax = zlims(im)

# orB = axs[0,1].imshow(im, vmin = vmin, vmax = vmax,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[0,1].set_title("Orion B")
# fig.colorbar(orB, ax=axs[0,1], **kw_colorb)


# #-----perseus ------------
# im = perseus_N
# imx = perseus_lons
# imy = perseus_lats

# vmin, vmax = zlims(im)

# per = axs[1,0].imshow(im, vmin = vmin, vmax = vmax,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[1,0].set_title("Perseus")
# fig.colorbar(per, ax=axs[1,0], **kw_colorb)


# for ax in axs.flat:
#     ax.set(xlabel='Galactic Longitude (deg)',\
#           ylabel = 'Galactic Latitude (deg)')


# plt.suptitle('column density maps')
# plt.show()


### velocity spectra ###
fig, axs = plt.subplots(2,2, figsize=(10,10), dpi=200, )

#----- Orion A ------------
axs[0,0].plot(orionA_vels,orionA_spectrum, color = 'gray')
axs[0,0].set_title('Orion A')

#----- Orion B ------------
axs[0,1].plot(orionB_vels,orionB_spectrum, color = 'gray')
axs[0,1].set_title('Orion B')


# #----- Perseus ------------
# axs[1,0].plot(perseus_vels,perseus_spectrum, color = 'gray')
# axs[1,0].set_title('Perseus')


# #----- California ------------
# axs[1,1].plot(california_vels,california_spectrum, color = 'gray')
# axs[1,1].set_title('California')


for ax in axs.flat:
    ax.set(xlabel='Velocity (km/s)',\
          ylabel='CO Temperature (K)')

    
plt.suptitle('Velocity Spectra')
plt.tight_layout()


plt.show()


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
