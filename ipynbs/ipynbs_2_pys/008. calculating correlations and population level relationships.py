#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


# https://iopscience.iop.org/article/10.1086/505417/fulltext/
#     
# how far away can a galaxy be and you can still resolve the clouds in it?
# (look at the angular resolution of ALMA, assume Rcloud~100pc, and see how far away you can see sth of this size)

# notes for check-in 4.10
# 1) all the virial parameters are of order unity, but perseus is still a little bit higher than the other ones
# 
# 2) recalculated orion A mass from wilson paper assuming it was at lld distance
#     - wilson Orion A mass: 94 * 1000 Msun
#     - Orion A mass I calculated: 54 * 1000 Msun (closest in mass to just region 2)
#     possible because the orion A boundaries i have defined are more strict (around the edges my 0th moment map looks more sparse than wilson)
#     
# 3) did the same rescaling / calculation for Orion B
#     - wilson Orion B mass: 74 * 1000 Msun
#     - Orion B mass I calculated: 42 * 1000 Msun (closest in mass to just region 2)
# 
#     
#     -- but !! when i don't do a sigma/noise cut-off when extracting the clouds, 
#     then the mass i calculate is pretty close to Wilson values
#         - Orion A: 99.8 * 1000 Msun
#         - Orion B: 73.2 * 1000 Msun
# 
# 3) virial parameter for perseus ~7, while other clouds are much closer to 1
#  
# 4) working on figuring out python toolkit to define & extract individual clouds
#     - any starting points nia knows?
#     
#     
#  
# moving forward:
#     - do the same correlations from the larson laws
#     - look at papers by katherine zucker for most updated distance values for clouds
#     - look at which clouds are outliers in these trends
#     - sigma vs r, mass vs r, virial mass vs luminous mass
#     - Lada & Dame 2020 -- surface density as a function of radius
#             https://arxiv.org/pdf/2006.08632.pdf 
#             section 2: look at Rice & MD's catalogs and whats going on there 
#     - clump finding algorithm starting place: https://github.com/akleroy/cpropstoo
#     - another method for identifying clouds -- dendrogram
# 
#         
#         look into: characterizin 3D structure of clouds using gaia
#       
# one of the ways to calculate the noise near a cloud is look at a patch of ism near it
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


# CLOUDS = cld.create_cloud_database()
CLOUDS = pd.read_csv('local_molecular_clouds_zucker')


CLOUDS


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


# #trying to plot an interactive 3D visualization, but it's too choppy
# import plotly.graph_objects as go


# smol = galaxy_cube[240:260,100:400,1100:1700]
# newsmolcube = smol*gmask_dilated[240:260,100:400,1100:1700]
# # plt.pcolormesh(newsmolcube[13])

# z, y, x = [np.arange(i) for i in newsmolcube.shape]
# X,Y,Z = np.meshgrid(x,y,z, indexing='ij')
# fig = go.Figure(data=go.Volume(
#     x=X.flatten(),
#     y=Y.flatten(),
#     z=Z.flatten(),
#     value=newsmolcube.flatten(),#np.transpose(newsmolcube,(1,2,0)).flatten(),
#     isomin=0.1,
#     opacity=0.1, # needs to be small to see through all surfaces
#     surface_count=20, # needs to be a large number for good volume rendering
#     ))


# fig.show()


##### IMAGE SEGMENTATION & PULLING OUT CLOUDS FROM THE GALAXY CUBE #####


#first attempt at algorithm that extracts clouds from galaxy cube
#following tutorial from https://github.com/lukepolson/youtube_channel/blob/main/Python%20Tutorial%20Series/image_processing1.ipynb

from skimage.segmentation import clear_border
from skimage import measure
from skimage.measure import label,regionprops
from scipy import ndimage as ndi
from scipy.ndimage import measurements, center_of_mass, binary_dilation, zoom
import plotly.graph_objects as go


# #using 3d connected component analysis
# #from https://pypi.org/project/connected-components-3d/
# import cc3d


# #plot a spectrum of the cube to see which indices have the most signal
# #in order to select which slices of the cube to look at
plt.plot(dc.spectrum1D(galaxy_cube))


#define a mask based on the standard 3*sigma)
rms_gcube = (np.nanmean(galaxy_cube**2))**0.5

#comparing 3*sigma vs 3*rms, the 3*sigma cutoff allows more data in, so will use that
gmask = galaxy_cube > 3*np.nanstd(galaxy_cube)
# gmask = galaxy_cube > 3*rms_gcube


plt.figure(figsize=(10,8))
plt.pcolormesh(gmask[250])


gmask_dilated = binary_dilation(gmask, iterations=2)


plt.figure(figsize=(10,8))
plt.pcolormesh(gmask_dilated[250])


# gmask_dilated_skimage_lbld = np.vectorize(label, signature='(n,m)->(n,m)')(gmask_dilated)


# plt.figure(figsize=(15,8))
# plt.pcolormesh(gmask_dilated_skimage_lbld[250],cmap=cmap0)
# plt.colorbar()

# plt.title('skimage dilated')


# gmask_skimage_lbld = np.vectorize(label, signature='(n,m)->(n,m)')(gmask)


# plt.figure(figsize=(15,8))
# plt.pcolormesh(gmask_skimage_lbld[250],cmap=cmap0)
# plt.colorbar()

# plt.title('skimage')


# gmask_cc3d_lbld = cc3d.connected_components(gmask,connectivity=6)


# plt.figure(figsize=(15,8))
# plt.pcolormesh(gmask_cc3d_lbld[250], cmap=cmap0)
# plt.colorbar()

# plt.title('cc3d')


##### END: IMAGE SEGMENTATION & PULLING OUT CLOUDS FROM THE GALAXY CUBE #####


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


#calculate radii from areas to compare to l,l&d table 2
orionA_R = np.sqrt(abs(orionA_A/np.pi)) #pc
orionB_R = np.sqrt(abs(orionB_A/np.pi)) #pc
perseus_R = np.sqrt(abs(perseus_A/np.pi)) #pc
california_R = np.sqrt(abs(california_A/np.pi)) #pc


#start plotting correlations
masses = np.array((orionA_M.value, orionB_M.value, perseus_M.value, california_M.value)) #kg

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


plt.scatter(radii, masses)


plt.scatter(radii,sigmas)


vir_pars 


orionA_R


orA_sigma


c.G.to(u.pc**3/(u.s**2*u.kg))


(5*(16*u.pc.to(u.m))*(2.7*1000*u.m/u.s)**2)/(c.G*5.5*1e4*c.M_sun)


# calculating Orion A and Orion B masses 
#from Wilson paper (https://arxiv.org/pdf/astro-ph/0411089.pdf)


#Orion A
#total mass
orA_wilson_M = 105.1*c.M_sun*1e3
#region breakdown
#0 =  distance in pc, 1 = mass 
orA_r1 = [521*u.pc, 12.3*c.M_sun*1e3]
orA_r2 = [465*u.pc, 69.5*c.M_sun*1e3]
orA_r3 = [412*u.pc, 23.3*c.M_sun*1e3]

#scaling to ^ distances to lld values
orA_lld_d = CLOUDS['distance']['Orion A']*CLOUDS['d unit']['Orion A']
orA_r1_M = orA_r1[1].to(u.kg) * (orA_lld_d.to(u.pc)/orA_r1[0])**2
orA_r2_M = orA_r2[1].to(u.kg) * (orA_lld_d.to(u.pc)/orA_r2[0])**2
orA_r3_M = orA_r3[1].to(u.kg) * (orA_lld_d.to(u.pc)/orA_r3[0])**2


orB_wilson_M = 82.3*c.M_sun*1e3
#region breakdown
#0 =  distance in pc, 1 = mass 
orB_r1n = [422*u.pc, 11.9*c.M_sun*1e3]
orB_r1f = [514*u.pc, 53.4*c.M_sun*1e3]
orB_r2 = [387*u.pc, 18.1*c.M_sun*1e3]
orB_sc = [150*u.pc, 0.34*c.M_sun*1e3]
orB_e = [120*u.pc, 0.13*c.M_sun*1e3]

#scaling to ^ distances to lld values
orB_lld_d = CLOUDS['distance']['Orion B']*CLOUDS['d unit']['Orion B']
orB_r1n_M = orB_r1n[1].to(u.kg) * (orB_lld_d.to(u.pc)/orB_r1n[0])**2
orB_r1f_M = orB_r1f[1].to(u.kg) * (orB_lld_d.to(u.pc)/orB_r1f[0])**2
orB_r2_M = orB_r2[1].to(u.kg) * (orB_lld_d.to(u.pc)/orB_r2[0])**2
orB_sc_M = orB_sc[1].to(u.kg) * (orB_lld_d.to(u.pc)/orB_sc[0])**2
orB_e_M = orB_e[1].to(u.kg) * (orB_lld_d.to(u.pc)/orB_e[0])**2


orA_lld_d


(orA_r1_M + orA_r2_M + orA_r3_M)/(c.M_sun*1e3)


orionA_M/(c.M_sun*1e3)


orB_lld_d


orionB_M/(c.M_sun*1e3)


(orB_r1n_M + orB_r1f_M + orB_r2_M + orB_sc_M + orB_e_M)/(c.M_sun*1e3)


orion


# ### 0th moment maps ###


# #!! need a better way of determining when lon coords need to be adjusted by +360 or not

# fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
# kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
# kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}
# wco = 1.5

# #-----orion A------------
# im = orionA_0th
# imx = orionA_lons +360 #put lon coords back in 360 form
# imy = orionA_lats

# vmin, vmax = zlims(im)

# orA = axs[0,0].imshow(im,\
#            extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
# axs[0,0].set_title("Orion A")

# wco_thresh = axs[0,0].contour(imx,imy,im,levels=[wco], colors='r')

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

# wco_thresh = axs[0,1].contour(imx,imy,im,levels=[wco], colors='r')


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

# wco_thresh = axs[1,0].contour(imx,imy,im,levels=[wco], colors='r')

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

# wco_thresh = axs[1,1].contour(imx,imy,im,levels=[wco], colors='r')

# cal_cbar = fig.colorbar(per, ax=axs[1,1], **kw_colorb)
# cal_cbar.add_lines(wco_thresh)


# for ax in axs.flat:
#     ax.set(xlabel='Galactic Longitude (deg)',\
#           ylabel = 'Galactic Latitude (deg)')


# plt.suptitle(r'$W_{CO}$ cut-off = 1.5')
# plt.show()


fig, axs = plt.subplots(1,1, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

#-----orion A------------
im = orionA_0th
imx = orionA_lons +360 #put lon coords back in 360 form
imy = orionA_lats

vmin, vmax = zlims(im)

orA = axs.imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs.set_title("Orion A")

wco_thresh = axs.contour(imx,imy,im,levels=[1.5], colors='r')

orA_cbar = fig.colorbar(orA, ax=axs, **kw_colorb)
orA_cbar.add_lines(wco_thresh)


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


# ### velocity spectra ###
# fig, axs = plt.subplots(2,2, figsize=(10,10), dpi=200, )

# #----- Orion A ------------
# axs[0,0].plot(orionA_vels,orionA_spectrum, color = 'gray')
# axs[0,0].set_title('Orion A')

# #----- Orion B ------------
# axs[0,1].plot(orionB_vels,orionB_spectrum, color = 'gray')
# axs[0,1].set_title('Orion B')


# # #----- Perseus ------------
# # axs[1,0].plot(perseus_vels,perseus_spectrum, color = 'gray')
# # axs[1,0].set_title('Perseus')


# # #----- California ------------
# # axs[1,1].plot(california_vels,california_spectrum, color = 'gray')
# # axs[1,1].set_title('California')


# for ax in axs.flat:
#     ax.set(xlabel='Velocity (km/s)',\
#           ylabel='CO Temperature (K)')

    
# plt.suptitle('Velocity Spectra')
# plt.tight_layout()


# plt.show()


# plt.figure(figsize=(15,3), dpi=200)

# im = galaxy_0th
# imx = glon
# imy = glat


# vmin,vmax = zlims(im)


# plt.imshow(im, vmin=vmin, vmax=vmax, \
#             extent = [imx[0],imx[-1], imy[0],imy[-1]],origin='lower')

# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# # plt.colorbar(orientation='vertical', label = 'Temperature (K)')

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
