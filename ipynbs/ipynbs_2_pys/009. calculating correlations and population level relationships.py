#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


# # notes for check-in 4.19 
# * started with creating a new catalog of clouds based on [Zucker 2019](https://iopscience.iop.org/article/10.3847/1538-4357/ab2388/pdf)
# 
# * re-plotted the clouds i was doing before (orions, perseus, california) as a santity check
#     - and saw that the cloud definitions are much more conservative -> explains the different mass & area values i'm getting from here compared to before (!! actually convert m to 1e3 Msun to compare this)
#     
# * Pegasus & Ursa Major had masses = 0
#     - looked at DHT survey & saw that it's between -30< b <30 : are there other galactic surveys that cover more latitudes? or do things have to be pieced together based on individual surveys for each cloud? 
#     
# 
# * still getting large virial parameters -- 
#     - could be that M / A / sigma values are off / Zucker et als definition of the clouds is an issue here maybe?
#     - what are other reasons that hese could be so off?
# 
# * plotting some larson law empirical relations: 
#     - sigma vs r:
#         - added a 'fit' based on the relation given in McKee & Ostriker, looks like it kind of fits?
#         
#     - mass vs r:
#         - go from here to constant surface density? 
#     - virial mass vs luminous mass
#     
#     
# ## moving forward notes from 4.19:
# - look at catalogs from: https://arxiv.org/pdf/2006.08632.pdf Rice & MD's catalogs and whats going on there 
# - possible: see if there are any surveys that are more complete, but the majority of clouds in the galaxy are near the plane so Dame's does a pretty good job
#     - 13CO survey: Heyer, FCRAO
# - what can 13CO tell us? (in comparison to 12CO)
# - take a look @ https://www.nature.com/articles/s41550-023-01901-5 -- how do molecular clouds form? uses a CII
#     
# 
# #### keep in back of mind: what are bits of information we want to know about the HI clouds -- what do papers have to say?
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


cmap0 = dc.cmap_grey0('viridis')


CLOUDS = pd.read_csv('local_molecular_clouds_zucker',\
                     sep=',',skiprows=4)


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


b0header = 'b0 (deg)'


l = len(b0header)
unit_str = ''
for i in range(l):
    if (b0header[i-1] == '(') or (b0header[i-1] == '['):
        sub_str = b0header[i:]
        for i in range(len(sub_str)):
            if not( sub_str[i] == ')' or sub_str[i] == ']'):
                unit_str += sub_str[i]
            else:
                break


#create arrays to hold quantities to be calculated
masses = np.zeros_like(CLOUDS['l0'])
areas = np.zeros_like(CLOUDS['l0'])
radii = np.zeros_like(CLOUDS['l0'])
sigmas = np.zeros_like(CLOUDS['l0'])
for i in range(26):    
    
    #add in a check that scans that lat & lon values of clouds & first confirms that 
    #they're in the galaxy cube, before calling cld.extract
    #something along the lines of 
    #if l0 in range(np.min(glat),np.max(glat))
    
    
    #extract the cloud
    extracted = cld.extract_cloud(i, CLOUDS, galaxy)
    
    cube = extracted[0]
    
    #0th moment map
    W = np.nansum(cube*dv, axis=0)
    
    #velocity spectrum
    spec = dc.spectrum1D(cube)
    
    #velocity dispersion
    #convert units from pixel/index to velocity ones(based on gvels<->index) + add units

    sigma = dv * gunits['vels'] * dc.dispersion_1sigma(spec)
    sigma = sigma.to(u.m/u.s)
    sigmas[i] = sigma.value
    sigma_unit = sigma.unit
    
    #c.G is in m^3/kg * s^2 so convert units to work
    #mass & area
    mass, area = cld.calc_mass_and_area(i, CLOUDS,galaxy,1.5)
    mass = mass.to(u.kg)
    masses[i] = mass.value
    mass_unit = mass.unit

    
    area = area.to(u.m**2)
    areas[i] = area.value
    area_unit = area.unit

    radius = np.sqrt(abs(area/np.pi))
    radii[i] = radius.value

#once all the quantities have been calculated, add to database
#eventually = write a function so units are saved in column names 
CLOUDS.insert(len(CLOUDS.columns),'mass',masses)
m_unit = [str(mass_unit)]*len(masses)
CLOUDS.insert(len(CLOUDS.columns),'m_unit',m_unit)

CLOUDS.insert(len(CLOUDS.columns),'area',areas)
a_unit = [str(area_unit)]*len(areas)
CLOUDS.insert(len(CLOUDS.columns),'a_unit',a_unit)

CLOUDS.insert(len(CLOUDS.columns),'velocity_dispersion',sigmas)
sig_unit = [str(sigma_unit)]*len(sigmas)
CLOUDS.insert(len(CLOUDS.columns),'sig_unit',sig_unit)


vir_pars = np.zeros_like(masses)
for i in range(len(vir_pars)):
     vir_pars[i] = (5*radii[i]*((sigmas[i])**2))/(c.G.value*masses[i])
        
        
CLOUDS.insert(len(CLOUDS.columns),'virial_parameter',vir_pars)


CLOUDS.sort_values(by=['d'])


# ### plotting empirical relationships
# 


## 1st Larson's Law

#----data----
radii_pc = radii*u.m.to(u.pc)
sigmas_kms = sigmas*u.m.to(u.km)

plt.scatter(radii_pc,sigmas_kms,s=10,label='data')


#as a check, will also plot the relationship reported in McKee & Ostriker
#plot this in solid lines for ranging from min to max
larson_sig = 0.65*(radii_pc**0.5)
#look for updated version of the sigma-r relatinoship

plt.scatter(radii_pc, larson_sig, color='red',alpha=0.3,label='Larsons Law')

plt.ylabel('Velocity Dispersion [km/s]')
plt.yscale('log')
plt.xlabel('R [pc]')
plt.xscale('log')
plt.title('Linewidth-Size relation')

plt.legend()
plt.show()


## M vs R

#----data----
radii_pc = radii*u.m.to(u.pc)
masses_1000Msun = (masses*u.kg)/(1e3*c.M_sun)

plt.scatter(radii_pc,masses_1000Msun,s=10,label='data')
# plt.scatter(radii_pc, radii_pc/2*c.G*(sigmas*u.km.to(u.pc))**2)

plt.ylabel(r'Masses [$10^3$$M_\odot$]')
plt.yscale('log')
plt.xlabel('R [pc]')
plt.xscale('log')
plt.title('Mass-Radius')
plt.legend()
plt.show()


## 2nd Larson's Law
names = CLOUDS['Cloud']
plt.scatter(vir_pars,names)

plt.vlines(1,'Aquila_S','Ursa_Major',linestyle='--',color='red',alpha=0.5)

plt.show()


orA_cube, orA_vels, orA_lats, orA_lons, orA_units = cld.extract_cloud(17,CLOUDS,galaxy)
orB_cube, orB_vels, orB_lats, orB_lons, orB_units = cld.extract_cloud(18,CLOUDS,galaxy)
per_cube, per_vels, per_lats, per_lons, per_units = cld.extract_cloud(21,CLOUDS,galaxy)
cal_cube, cal_vels, cal_lats, cal_lons, cal_units = cld.extract_cloud(2,CLOUDS,galaxy)


orA_0th = np.nansum(orA_cube*dv, axis=0)
orB_0th = np.nansum(orB_cube*dv, axis=0)
per_0th = np.nansum(per_cube*dv, axis=0)
cal_0th = np.nansum(cal_cube*dv, axis=0)


### 0th moment maps ###


#!! need a better way of determining when lon coords need to be adjusted by +360 or not

fig, axs = plt.subplots(2,2, figsize=(15,15), dpi=200)
kw_imshow = {'cmap':cmap0, 'origin':'lower',  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

#-----orion A------------
im = orA_0th
imx = orA_lons +360 #put lon coords back in 360 form
imy = orA_lats

vmin, vmax = zlims(im)

orA = axs[0,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,0].set_title("Orion A")

wco_thresh = axs[0,0].contour(imx,imy,im,levels=[1.5], colors='r')

orA_cbar = fig.colorbar(orA, ax=axs[0,0], **kw_colorb)
orA_cbar.add_lines(wco_thresh)


#-----orion B------------
im = orB_0th
imx = orB_lons +360 #put lon coords back in 360 form
imy = orB_lats

vmin, vmax = zlims(im)

orB = axs[0,1].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[0,1].set_title("Orion B")

wco_thresh = axs[0,1].contour(imx,imy,im,levels=[1.5], colors='r')


orB_cbar = fig.colorbar(orB, ax=axs[0,1], **kw_colorb)
orB_cbar.add_lines(wco_thresh)


#-----perseus------------
im = per_0th
imx = per_lons
imy = per_lats

vmin, vmax = zlims(im)

per = axs[1,0].imshow(im,\
           extent = [imx[0],imx[-1], imy[0], imy[-1]], **kw_imshow)
axs[1,0].set_title("Perseus")

wco_thresh = axs[1,0].contour(imx,imy,im,levels=[1.5], colors='r')

per_cbar = fig.colorbar(per, ax=axs[1,0], **kw_colorb)
per_cbar.add_lines(wco_thresh)


#-----california------------
im = cal_0th
imx = cal_lons
imy = cal_lats

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


#### GALAXY 0TH MOMENT ####
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
