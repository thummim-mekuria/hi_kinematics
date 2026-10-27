#!/usr/bin/env python
# coding: utf-8

from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

import datacube as dc
import clouds as cld

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c

from astropy.visualization import ZScaleInterval as zscale

from scimes import SpectralCloudstering
from astrodendro import Dendrogram
#values from tutorial(https://scimes.readthedocs.io/en/latest/tutorial.html)
#may need to change later
dend_sigma = 0.3 #K, noise level
dend_ppb = 1.3 #pixels/beam
from astrodendro import ppv_catalog

zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}

#all-sky map from HI4PI
LAB = dc.unpack_fits_cube('labh.fits', 'LAB')
HI_cube, HI_vels, HI_lats, HI_lons, HI_units, HI_head = LAB.values()

HI_dv = ((HI_vels[1] - HI_vels[0])*HI_units['vels']).to(u.km/u.s)
HI_db = (HI_lats[1] - HI_lats[0])*HI_units['lats']
HI_dl = (HI_lons[1] - HI_lons[0])*HI_units['lons']

#DHT CO survey
DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
CO_cube, CO_vels, CO_lats, CO_lons, CO_units, CO_head = DHT_CO.values()

CO_dv = ((CO_vels[1] - CO_vels[0])*CO_units['vels']).to(u.km/u.s)
CO_db = (CO_lats[1] - CO_lats[0])*CO_units['lats']
CO_dl = (CO_lons[1] - CO_lons[0])*CO_units['lons']

CO_0th = np.nansum(CO_cube*CO_dv.value, axis=0)

Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)

orAparent = cld.get_parent_HI(17, Zucker_df, DHT_CO, LAB,
                      CO_padding=0, HI_padding=10,  v_padding=0,
                      show_0th = True, show_spec = False, sigma = 3)

orAparent['parent']['comp']['cube'].shape

#plot surface density
cld.surface_density(LAB, 1.5, '21cm', show_Sigma=True)

Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)

orion = cld.extract_cloud(17, Zucker_df, DHT_CO, p_padding=60)

orion_cube = orion['cube']
orion_head = orion['header']
orion_units = orion['units']

plt.imshow(np.nansum(orion_cube, axis=0), origin='lower')

orion_dend = Dendrogram.compute(orion_cube,
                               min_value=dend_sigma,
                               min_delta=2*dend_sigma,
                               min_npix=3*dend_ppb)

orion_mtdt = {}
orion_mtdt['data_unit'] = u.Jy
orion_catalog = ppv_catalog(orion_dend, orion_mtdt)

orion_dclust = SpectralCloudstering(orion_dend, orion_catalog, orion_head)

trunks_flat = np.nanmedian(orion_dclust.trunks_asgn.data, axis=0)

trunks_flat.shape

orion_dclust.trunks_asgn.data.shape

plt.imshow(trunks_flat, origin='lower')
plt.colorbar()

plt.pcolormesh(orion_dclust.trunks_asgn.data[6])
plt.colorbar()
plt.show()

HI_dend = Dendrogram.compute(HI_cube,
                                min_value=dend_sigma,
                                min_delta=2*dend_sigma,
                                min_npix=3*dend_ppb)

HI_mtdt = {}
HI_mtdt['data_unit'] = u.Jy
HI_catalog = ppv_catalog(HI_dend, HI_mtdt)

HI_dclust = SpectralCloudstering(HI_dend, HI_catalog, HI_head)

HI_dclust.clusters_asgn.data.shape

plt.pcolormesh(HI_dclust.clusters_asgn.data[80])

#following along tutorial (https://scimes.readthedocs.io/en/latest/tutorial.html)
#for calculating dendograms
orion_hdu = fits.open('DHT27_Orion_raw.fits')[0]
orion_data = orion_hdu.data
orion_head = orion_hdu.header

orion_dend = Dendrogram.compute(orion_data,
                                min_value=dend_sigma,
                                min_delta=2*dend_sigma,
                                min_npix=3*dend_ppb)

metadata = {}
metadata['data_unit'] = u.Jy
orion_cat = ppv_catalog(orion_dend, metadata)

orion_dclust = SpectralCloudstering(orion_dend, orion_cat, orion_head)

orion_dclust.showdendro()

plt.pcolormesh(orion_dclust.clusters_asgn.data[80])

# would be interesting/helpful to isolate the different velocity peaks & see how they map onto spatial features

HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0)

### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = HI_0th
imx = HI_lons
imy = HI_lats

vmin,vmax = zlims(im)

plt.imshow(im, vmin=vmin, vmax=vmax,
            extent=[imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label=r'$W_{HI}$ $(K\times km/s)$')

plt.title('Galaxy in HI: 0th moment')
plt.show()

#surface density maps for solar neighborhood clouds + their HI envelopes FROM LAB
for i in range(len(Zucker_df)):
    cloud, parent = cld.get_parent_HI(i, Zucker_df, DHT_CO, LAB,
                      CO_padding=20, HI_padding=5,  v_padding=10,
                      show_0th=False, show_spec=False, sigma=3).values()

    cld.surface_density(cloud['comp'],1.5, 'CO', show_Sigma=True)
    cld.surface_density(parent['comp'],1.5, '21cm', show_Sigma=True)

#0th moment maps & velcoity spectra for solar neighborhood clouds + their HI envelopes FROM LAB
for i in range(len(Zucker_df)):
    cld.get_parent_HI(i, Zucker_df, DHT_CO, LAB,
                      CO_padding=20, HI_padding=5,  v_padding=10,
                      show_0th=True, show_spec=True, sigma=3)

# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
