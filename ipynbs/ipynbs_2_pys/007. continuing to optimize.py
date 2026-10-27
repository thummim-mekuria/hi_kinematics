#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


# https://iopscience.iop.org/article/10.1086/505417/fulltext/
#     
# how far away can a galaxy be and you can still resolve the clouds in it?
# (look at the angular resolution of ALMA, assume Rcloud~100pc, and see how far away you can see sth of this size)

# notes for check-in 3.23
# 1) saw from https://arxiv.org/pdf/astro-ph/0411089.pdf (paper linked in the survey page for orions) that both Orion A & B have different regions within them that are at different distances + with different masses (Table 2)
#     a) one indication that i see for separating out regions is using vel-lon / vel-lat plots, in which you can see which portions of the cloud show up in which velocity regimes
#         - how to calculate distance of each segment from these maps?
#         - meta question: if these regions are kinematically separable from each other, what's the thing that ties them together into one cloud that == Orion A? or == Orion B?
#         
#         -- to compare to the masses in the wilson paper -- take the masses they report for each region, scale  to the distances i'm using + add up , and see how that compares to the masses i'm calculating
#         read wilson paper more in depth to see why they are identifying regions the way they are
#         
#         look into: characterizin 3D structure of clouds using gaia
# 
# 2) When plotting contours, can see that:
#     - boundaries for the Orions are better traced by W_co = 1.5
#         - this is the cut-off lld report using & it corresponds to dame et als (2001) cut off of 3 * rms
#     - boundaries for California/Perseus are better traced by W_co = 10 
#         - this is the cut-off that lld plot in their appendices
#         - made a mask for perseus using this, and mass value calculated is much closer (0.96xlld, vs 1.57xlld before)
#         - applying the same calculation for california (assuming that the same data cube should have similar noise levels etc) does much worse  
#         (1.07xlld before, vs 0.48xlld now)
#         
#         
# one of the ways to calculate the noise near a cloud is look at a patch of ism near it
# --next thing: figure out an algorithm or other method to pull out clouds from the galaxy cube
# 
# 
# 3) fixed radius calculations, so values are now much closer to lld values
#      
# 
# 
# moving forward notes from last time:
# --> start plotting correlations, virial parameter, velocity dispersions, 
# 
# start thinking about: what are bits of information we want to know about the HI clouds -- what do papers have to say?
# 


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
CLOUDS.at['Perseus','longitude (l)'] = (155,162)
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


# at a stopping point -- currently working on how to create the cloud data base in an easier way
# maybe i initialize it the same (or create a csv) & the task is to read it in & select which clouds i want for this specific task 


# name_array = ['Orion A', 'Perseus']
# columns = ['longitude (l)', 'latitude (b)', 'l, b unit', 'distance', 'd unit']

# df = pd.DataFrame(columns = columns)
# df.insert(0,'clouds', name_array)

# #find each cloud so i can set values for its physical properties
# if 
# orA_n = np.where(df['clouds'] == 'Orion A')[0][0]
# if 'Orion B' in name_array:
#     orB_n = np.where(df['clouds'] == 'Orion B')[0][0]
# else:
#     orB_n = -1
# per_n = np.where(df['clouds'] == 'Perseus')[0][0]
# # cal_n = np.where(df['clouds'] == 'California')[0][0]
# # oph_n = np.where(df['clouds'] == 'Ophiuchus')[0][0]
# # mon_n = np.where(df['clouds'] == 'Mon R2')[0][0]
# # tau_n = np.where(df['clouds'] == 'Taurus')[0][0]


# df.at[orB_n,'longitude (l)'] = 'orion b lon'


# df


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


#calcualte masses & areas
orionA_M,orionA_A = cld.calc_mass_and_area('Orion A',CLOUDS,galaxy,1.5) #g
orionB_M,orionB_A = cld.calc_mass_and_area('Orion B',CLOUDS,galaxy,1.5)
perseus_M, perseus_A = cld.calc_mass_and_area('Perseus',CLOUDS,galaxy,1.5)
california_M, california_A = cld.calc_mass_and_area('California',CLOUDS,galaxy,1.5)


#recalculating for perseus using a Wco=10 cutoff instead
perseus_M_10, perseus_A_10 = cld.calc_mass_and_area('Perseus',CLOUDS,galaxy,10)
perseus_R_10 = np.sqrt(abs(perseus_A_10/np.pi))


#calculate radii from areas to compare to l,l&d table 2
orionA_R = np.sqrt(abs(orionA_A/np.pi)) #pc
orionB_R = np.sqrt(abs(orionB_A/np.pi)) #pc
perseus_R = np.sqrt(abs(perseus_A/np.pi)) #pc
california_R = np.sqrt(abs(california_A/np.pi)) #pc


#start plotting correlations
masses = np.array((orionA_M.value, orionB_M.value, perseus_M_10.value, california_M.value)) #kg

#convert pc -> m to match units with c.G
radii = np.array((orionA_R.to(u.m).value, orionB_R.to(u.m).value, perseus_R.to(u.m).value, california_R.to(u.m).value)) #m

#convert km/s -> m/s to match units with c.G
sigmas = np.array((orA_sigma.to(u.m/u.s).value,orB_sigma.to(u.m/u.s).value, per_sigma.to(u.m/u.s).value, cal_sigma.to(u.m/u.s).value)) #m/s

#calculate virial parameter
#c.G is in m^3/kg * s^2

#create new array to hold values of the virial parameter 
vir_pars = np.zeros_like(masses)
for i in range(len(vir_pars)):
     vir_pars[i] = (5*radii[i]*((sigmas[i])**2))/(c.G.value*masses[i])


vir_pars


# the virial parameter speaks to the evolutionary status of clouds (is it expanding? balanced? contracting?)
# if it is much greater than 1, interpretation is it would be expanding
# if it is much less than 1, interpretation is it would be contracting
# -- no straightforward life-story for clouds, star formation affects things
# 
# virial parameter should be unitless
# think in M = Msun, R = parsec, sigma = km/s --> G should be ~1/432
# -- look at the np calculations to make sure units work out


5*orionA_R.to(u.m)*(orA_sigma.to(u.m/u.s))**2 /(orionA_M*c.G)


orionA_M


orionA_R


orA_sigma


fig,axs = plt.subplots(2,2,figsize=(10,10), dpi=200)
imx = perseus_lons
imy = perseus_lats

kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm(), \
             'extent':[imx[0],imx[-1], imy[0], imy[-1]]}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


axs[0,0].imshow(perseus_M.value, **kw_imshow)
axs[0,0].set_title("dm of Perseus with 3*rms cut-off")

axs[0,1].imshow(perseus_A.value, **kw_imshow)
axs[0,1].set_title("dA of Perseus with 3*rms cut-off")

axs[1,0].imshow(perseus_M_10.value, **kw_imshow)
axs[1,0].set_title("dm of Perseus with Wco=10 cut-off")

axs[1,1].imshow(perseus_A_10.value, **kw_imshow)
axs[1,1].set_title("dA of Perseus with Wco=10 cut-off")


plt.show()
# fig.colorbar(orB, ax=axs[0,1], **kw_colorb)


# 


# #comparing to lld values
# (orionA_M.to(u.kg)/(1e3*c.M_sun)),(orionA_M.to(u.kg)/(1e3*c.M_sun))/85.2,\
# (orionB_M.to(u.kg)/(1e3*c.M_sun)),(orionB_M.to(u.kg)/(1e3*c.M_sun))/65.0, \
# (perseus_M.to(u.kg)/(1e3*c.M_sun)),(perseus_M.to(u.kg)/(1e3*c.M_sun))/19.1,\
# (california_M.to(u.kg)/(1e3*c.M_sun)),(california_M.to(u.kg)/(1e3*c.M_sun))/139.1


orionA_R.to(u.pc)/23.,\
orionB_R.to(u.pc)/21.6,\
perseus_R.to(u.pc)/13.1,\
california_R.to(u.pc)/35.


# orionA_A.to(u.pc**2), orionB_A.to(u.pc**2), perseus_A.to(u.pc**2)


### 0th moment maps ###


#!! need a better way of determining when lon coords need to be adjusted by +360 or not

fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}
wco = 1.5

#-----orion A------------
im = orionA_0th
imx = orionA_lons +360 #put lon coords back in 360 form
imy = orionA_lats

vmin, vmax = zlims(im)

orA = axs[0,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,0].set_title("Orion A")

wco_thresh = axs[0,0].contour(imx,imy,im,levels=[wco], colors='r')

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

wco_thresh = axs[0,1].contour(imx,imy,im,levels=[wco], colors='r')


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

wco_thresh = axs[1,0].contour(imx,imy,im,levels=[wco], colors='r')

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

wco_thresh = axs[1,1].contour(imx,imy,im,levels=[wco], colors='r')

cal_cbar = fig.colorbar(per, ax=axs[1,1], **kw_colorb)
cal_cbar.add_lines(wco_thresh)


for ax in axs.flat:
    ax.set(xlabel='Galactic Longitude (deg)',\
          ylabel = 'Galactic Latitude (deg)')


plt.suptitle(r'$W_{CO}$ cut-off = 1.5')
plt.show()


### 0th moment maps ###


#!! need a better way of determining when lon coords need to be adjusted by +360 or not

fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}
wco = 10

#-----orion A------------
im = orionA_0th
imx = orionA_lons +360 #put lon coords back in 360 form
imy = orionA_lats

vmin, vmax = zlims(im)

orA = axs[0,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,0].set_title("Orion A")

wco_thresh = axs[0,0].contour(imx,imy,im,levels=[wco], colors='r')

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

wco_thresh = axs[0,1].contour(imx,imy,im,levels=[wco], colors='r')


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

wco_thresh = axs[1,0].contour(imx,imy,im,levels=[wco], colors='r')

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

wco_thresh = axs[1,1].contour(imx,imy,im,levels=[wco], colors='r')

cal_cbar = fig.colorbar(per, ax=axs[1,1], **kw_colorb)
cal_cbar.add_lines(wco_thresh)


for ax in axs.flat:
    ax.set(xlabel='Galactic Longitude (deg)',\
          ylabel = 'Galactic Latitude (deg)')


plt.suptitle(r'$W_{CO}$ cut-off = 10')
plt.show()


#recalculating for perseus using a Wco=10 cutoff instead
perseus_M_10, perseus_A_10 = cld.calc_mass_and_area('Perseus',CLOUDS,galaxy,10)
california_M_10, california_A_10 = cld.calc_mass_and_area('California',CLOUDS,galaxy,10)


perseus_M_10, perseus_M


print('perseus mass 3rms',(perseus_M.to(u.kg)/(1e3*c.M_sun)),\
      '\nperseus mass 3rms lld fraction',(perseus_M.to(u.kg)/(1e3*c.M_sun))/19.1,\
    '\n\nperseus mass 10Wco',(perseus_M_10.to(u.kg)/(1e3*c.M_sun)),\
      '\nperseus mass 10Wco lld fraction',(perseus_M_10.to(u.kg)/(1e3*c.M_sun))/19.1)


print('california mass 3rms',(california_M.to(u.kg)/(1e3*c.M_sun)),\
      '\ncalifornia mass 3rms lld fraction',(california_M.to(u.kg)/(1e3*c.M_sun))/139.1,\
    '\n\ncalifornia mass 10Wco',(california_M_10.to(u.kg)/(1e3*c.M_sun)),\
      '\ncalifornia mass 10Wco lld fraction',(california_M_10.to(u.kg)/(1e3*c.M_sun))/139.1)


orA_vLon = np.nansum(orionA_cube*dv, axis=1)
orB_vLon = np.nansum(orionB_cube*dv, axis=1)


# ### 0th moment maps ###


# #!! need a better way of determining when lon coords need to be adjusted by +360 or not

# fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
# kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
# kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

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


# #-----perseus------------
# im = perseus_0th
# imx = perseus_lons
# imy = perseus_lats

# vmin, vmax = zlims(im)

# per = axs[1,0].imshow(im,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[1,0].set_title("Perseus")

# wco_thresh = axs[1,0].contour(imx,imy,im,levels=[1.5], colors='r')

# per_cbar = fig.colorbar(per, ax=axs[1,0], **kw_colorb)
# per_cbar.add_lines(wco_thresh)


# #-----california------------
# im = california_0th
# imx = california_lons
# imy = california_lats

# vmin, vmax = zlims(im)

# per = axs[1,1].imshow(im,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[1,1].set_title("California")

# wco_thresh = axs[1,1].contour(imx,imy,im,levels=[1.5], colors='r')

# cal_cbar = fig.colorbar(per, ax=axs[1,1], **kw_colorb)
# cal_cbar.add_lines(wco_thresh)


# for ax in axs.flat:
#     ax.set(xlabel='Galactic Longitude (deg)',\
#           ylabel = 'Galactic Latitude (deg)')


# plt.suptitle('0th moment maps')
# plt.show()


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
