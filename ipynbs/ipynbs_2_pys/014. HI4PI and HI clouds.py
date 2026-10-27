#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


# %%time

import datacube as dc
import clouds as cld

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c

from astropy.visualization import ZScaleInterval as zscale


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


# %%time

#all-sky map from HI4PI
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI_cube, HI_vels, HI_lats, HI_lons, HI_units, HI_header = HI4PI.values()

HI_dv = ((HI_vels[1] - HI_vels[0])*HI_units['vels']).to(u.km/u.s)
HI_db = (HI_lats[1] - HI_lats[0])*HI_units['lats']
HI_dl = (HI_lons[1] - HI_lons[0])*HI_units['lons']


# %%time

#DHT CO survey
DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
CO_cube, CO_vels, CO_lats, CO_lons, CO_units, CO_header = DHT_CO.values()

#assume arrays have linear spacing
CO_dv = ((CO_vels[1] - CO_vels[0])*CO_units['vels']).to(u.km/u.s)
CO_db = (CO_lats[1] - CO_lats[0])*CO_units['lats']
CO_dl = (CO_lons[1] - CO_lons[0])*CO_units['lons']

CO_0th = np.nansum(CO_cube*CO_dv.value, axis=0)


Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)


# would be interesting/helpful to isolate the different velocity peaks & see how they map onto spatial features


cropped_HI = HI_cube[300:600,900:1200, 2000:2300] #numbers picked pretty arbitrarily


cropped_HI_0th = np.nansum(cropped_HI*HI_dv.value, axis=0) 


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = cropped_HI_0th
imx = HI_lons[2000:2300]
imy = HI_lats[900:1200]


vmin,vmax = zlims(im) 


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')


plt.title('Galaxy in HI: 0th moment')
plt.show()


#surface density maps for solar neighborhood clouds + their HI envelopes
for i in range(len(Zucker_df)):
    cloud, parent = cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, \
                      CO_padding=20, HI_padding=40,  v_padding=10, \
                      show_0th = False, show_spec = False, sigma = 3).values()
    
    cld.mass_and_area(i, Zucker_df, cloud['comp'],1.5, 'CO', show_Sigma=True)
    cld.mass_and_area(i, Zucker_df, parent['comp'],1.5, '21cm', show_Sigma=True)


#0th moment maps & velcoity spectra for solar neighborhood clouds + their HI envelopes

for i in range(len(Zucker_df)):
    cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, \
                      CO_padding=20, HI_padding=40,  v_padding=10, \
                      show_0th = True, show_spec = True, sigma = 3)


# (HI_cube.size * HI_cube.itemsize)/1e9 #the size of HI_cube in gigs = 35.3


HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 
# ^^ this takes too long to run, so instead of making a 0th moment map of the whole
# thing, i will just cut out the regions around the solar neighborhood clouds


# ### GALAXY 0TH MOMENT (in HI) ####
# plt.figure(figsize=(15,3), dpi=200)

# im = HI_0th
# imx = HI_lons
# imy = HI_lats


# vmin,vmax = zlims(im) 


# plt.imshow(im, vmin=vmin, vmax=vmax, \
#             extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')


# plt.title('Galaxy in HI: 0th moment')
# plt.show()


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
