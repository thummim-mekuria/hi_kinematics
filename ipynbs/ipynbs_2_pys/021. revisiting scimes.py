#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


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

from scimes import SpectralCloudstering
from astrodendro import Dendrogram

from astrodendro import ppv_catalog


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


def rms_noise(thing): ##turns out this is just = np.nanstd()
    mean = np.nanmean(thing)
    noise_square = (thing - mean)**2
    
    return np.sqrt(np.nanmean(noise_square))


# * [dendogram documentation](https://dendrograms.readthedocs.io/en/latest/using.html)
# * [scimes documentation](https://scimes.readthedocs.io/en/latest/tutorial.html)


import orion_dendogram as od


dend = od.compute()


len(dend.trunk), len(dend)


# od.launch_viewer(dend)


# orion = dc.unpack_fits_cube('DHT27_Orion_interp.fits', 'DHT')


# im = np.nansum(orion.cube*orion.dv.value, axis = 0)
# vmin, vmax = zlims(im)

# imx, imy = orion.lons, orion.lats
# plt.imshow(im, vmin = vmin, vmax = vmax, \
#            extent = [imx[0], imx[-1], imy[0], imy[-1]],**kw_imshow)

# plt.colorbar()

# plt.show()


# #values from tutorial(https://scimes.readthedocs.io/en/latest/tutorial.html)
# #may need to change later
# dend_sigma = 0.3 #K, noise level
# dend_ppb = 1.3 #pixels/beam


# orion_dend = Dendrogram.compute(orion.cube,\
#                                min_value = 3*dend_sigma,\
#                                min_delta = 2*dend_sigma,\
#                                min_npix = 3*dend_ppb)


# plotter = orion_dend.plotter()


# lines = plotter.get_lines()


# fig, axs = plt.subplots(1, 1, figsize = (5, 5))
# axs.add_collection(lines)
# plt.show()


# fig, axs = plt.subplots(1, 1, figsize = (5, 5))
# plotter.plot_tree(axs)
# plt.show()


dir(dend)


dend.leaves[10].indices()


dend.data.shape


print(dend.structure_at((24,57,157)))


for i in dend.trunk:
    npx = i.get_npix()
    m = i.get_mask()
    plt.imshow(np.nansum(m * dend.data, axis = 0), **kw_imshow)
    plt.show()


orion_mtdt = {}
orion_mtdt['data_unit'] = u.Jy #works for u.Jy (but not u.K ?? )


(1.09/557)*1e3


cogal = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')


np.nanstd(cogal.cube)


dco = Dendrogram.compute(cogal.cube, 3*0.3, 2*0.3,50)


gal_mtdt = {}
gal_mtdt['data_unit'] = u.Jy


galcat


orcat


orhead = dc.unpack_fits_cube('DHT27.fits','DHT').data_header


galhead = cogal.data_header


orclusts = orion_dclust.clusters_asgn.data 


galclusts = gal_clust.clusters_asgn.data 


orion_dclust.clusters_asgn.writeto('orion_dend_clusters.fits',overwrite=True)


gal_clust.clusters_asgn.writeto('gal_dend_clusters.fits',overwrite=True)


for cnum in gal_clust.clusters:
    mskd = (galclusts == cnum)*galclusts
    plt.imshow(np.nansum(mskd,axis=0),**kw_imshow)
    plt.title(str(cnum))
    plt.show()


orcat


galcat


len(orion_dclust.clusters)


orclob = dc.unpack_fits_cube('DHT27.fits','DHT')
orcube = orclob.cube


orcube[30].shape


orclob.lons.shape, orclob.lats.shape


cogal2 = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')


fig = plt.figure(figsize=(20,15))
# ax = fig.add_subplot(projection='3d')
x,y = np.meshgrid(cogal2.lons, cogal2.lats)
z = (np.nansum(cogal2.cube[:,:,:], axis=0))

# ax.plot_wireframe(x,y,z,cmap = 'gray',alpha = 0.9,linewidth=0.3)
# ax.view_init(azim=90,roll=0,elev=60)

# ax2 = fig.add_subplot()
# ax2.imshow(z,origin='lower',cmap=cmap0)
# # fig.colorbar(ax2)
# plt.show()


plt.figure(figsize=(20,15))
z2 = (np.nansum(cogal2.cube[:,:,:],axis=1))
plt.imshow(z2,origin='lower',cmap=cmap0)
plt.show()


plt.figure(figsize=(20,15))
z2 = (np.nansum(cogal2.cube*cogal2.dv.value,axis=0))
plt.imshow(z2,origin='lower',cmap=cmap0)
plt.show()


fig = plt.figure(figsize=(20,15))
# ax = fig.add_subplot(projection='3d')
x,y = np.meshgrid(cogal.lons, cogal.lats)
z = np.sqrt(np.nansum(cogal.cube, axis=0))

# ax.plot_wireframe(x,y,z,cmap = 'gray',alpha = 0.9,linewidth=0.3)
# ax.view_init(azim=90,roll=0,elev=60)

ax2 = fig.add_subplot()
ax2.imshow(z,origin='lower',cmap=cmap0)
# fig.colorbar(ax2)
plt.show()

plt.figure(figsize=(20,15))
z2 = np.nansum(cogal.cube,axis=0)
plt.imshow(z2,origin='lower',cmap=cmap0)
plt.show()


np.nanmax(cogal.cube),np.nanmin(cogal.cube)


np.nanmax(cogal.cube),np.nanmin(cogal.cube)


xs = np.arange(0,10,0.1)


plt.plot(xs,xs)
plt.plot(xs,np.sqrt(xs))
# plt.plot(xs,xs**2)


# dir(orion_dclust)


# trunks_flat = np.nanmedian(orion_dclust.trunks_asgn.data, axis = 0)


# trunks_flat.shape


# orion_dclust.trunks_asgn.data.shape


# plt.imshow(trunks_flat, origin = 'lower')
# plt.colorbar()


# plt.pcolormesh(orion_dclust.trunks_asgn.data[6])
# plt.colorbar()

# plt.show()


# HI_dend = Dendrogram.compute(HI_cube, \
#                                 min_value=dend_sigma,\
#                                 min_delta=2*dend_sigma,\
#                                 min_npix=3*dend_ppb)


# HI_mtdt = {} #metadata
# HI_mtdt['data_unit'] = u.Jy
# HI_catalog = ppv_catalog(HI_dend, HI_mtdt)


# HI_dclust = SpectralCloudstering(HI_dend, HI_catalog, HI_head)


# HI_dclust.clusters_asgn.data.shape


# plt.pcolormesh(HI_dclust.clusters_asgn.data[80])


# would be interesting/helpful to isolate the different velocity peaks & see how they map onto spatial features


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
