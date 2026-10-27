#!/usr/bin/env python
# coding: utf-8


#re-evaluate HI velocity dispersion
#the main point is not to get all the known gmc clouds ever, but just to choose a sample, 
#then base the selection of the non-gmc hi clouds based on ^
#it will be helpful to build in ways to: inspect spectra, show 0th moment maps (in a way where you can 
#overlap the CO & HI data, select different velocity ranges in real time, etc)
#deliverable for next meeting: ^ 
#also for next time: revisit the timelines for the paper, 
#within the next two meetings, have the techniques nailed down for
#selecting the gmc-full & gmc-less hi clouds
#in order to move towards quantitatively analyzing the different populations
#add a column for surface density (in cat)


#look at the literature for why the calculated distances for clouds near l~0 might be 
#larger than the others


cat[cat['accumulation_radius<pc>']<cat['effective_radius<pc>']]


import datacube as dc
import clouds as cld
import interpolate as ip
import bootleg as btlg

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import math

# %matplotlib ipympl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram,periodic_neighbours


1.040664e35*u.kg.to(c.M_sun)


plt.plot(DHT_CO.vels,spec/np.nanmax(spec))
plt.plot(1e-3*HI4PI.vels,spec_hi/np.nanmax(spec_hi))
plt.xlim(-50,50)
plt.show()


big_cloud_co = dnd[167].get_mask()*DHT_CO.cube
spec = dc.spectrum1D(big_cloud_co)


parents = fits.open('HI_parents_from_accR.fits')[0].data


lims_lb = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]


wparents = np.nansum(parents,axis=0)


plt.imshow(wparents,origin='lower',extent=lims_lb)


plt.imshow(wparents,origin='lower',extent=lims_lb)
plt.ylabel(r'galactic latitude$^o$')
plt.xlabel(r'galactic longitude$^o$')
plt.xlim(180,150)
plt.ylim(-30,0)


spec_hi = dc.spectrum1D((parents[:,450:900,1600:2100]).astype('float32'))


plt.show()


cat[cat['accumulation_radius<deg>']>10]


DHT_CO = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')

DHT_CO.cube = DHT_CO.cube[209:285,:,:] #trim to narrower velocity range
DHT_CO.vels = DHT_CO.vels[209:285]

smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole) #interpolate bc 'mom' is staccatoed
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641] #take out |b|<5 for |l|<100
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]

DHT_CO.cube = wo_gal_cent


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]


dnd = Dendrogram.load_from('COGAL_dendrogram.fits')
cat = pd.read_csv('physcl_prprts_catalog.csv')


bad_rs_pc = cat[cat['effective_radius<pc>']>cat['accumulation_radius<pc>']]
bad_rs_deg = cat[np.sqrt(cat['area_sky<deg2>']/np.pi)>cat['accumulation_radius<deg>']]


percent_bad_deg = 100*(np.abs(np.sqrt(bad_rs_deg['area_sky<deg2>']/np.pi)\
                        - bad_rs_deg['accumulation_radius<deg>']\
                            )/bad_rs_deg['accumulation_radius<deg>']).to_numpy()


percent_bad_pc = 100*(np.abs(bad_rs_pc['effective_radius<pc>'] -\
                       bad_rs_pc['accumulation_radius<pc>']\
                           )/bad_rs_pc['accumulation_radius<pc>']).to_numpy()


percent_bad_deg


percent_bad_pc


cat.loc[100]


cat.loc[120]


# aquila: 217(102),248(121),249(122) -- d: x2
# california: 173(75),133(52), 136(53) -- d:x0.5
# perseus: 214(100) -- d:x1
# taurus: 264(134) -- d:x1.5
# ora: 246(120) -- d:x2
# orb: 250(123) -- d:x2
# crossbones:304(171) -- d:1-1.5x
# monr2: 301(168) -- d:1.5x
# scissors: 202(93)
#     
# chamaeleon: 218(103),216(101) -- d:xBIG
# corona: 275(145) -- d:xBIG
#         
# hercules: 257(129) -- d:1x
# polaris: 157(66) -- d:1-1.5x
#         
# ophiuchus: 260(131) -- d:xBIG
# 
# maybe export this or write some sort of routine to find the zucker clouds
# (also find the one that i originally took out bc d>1kpc and add back?)


whi_raw = np.nansum(HI4PI.cube,axis=0)


plt.hist(cat['accumulation_radius<deg>'],bins=50)
plt.xlabel('accumulation radii [deg]')
plt.show()


3/HI4PI.db


4/HI4PI.db


hdu_smooth = fits.PrimaryHDU(HI4PI_smooth)


hdu_smooth.writeto('smoothed_HI4PI_sig40.fits',overwrite=True)


whi_smooth = np.nansum(HI4PI_smooth,axis=0)


x,y = np.meshgrid(HI4PI.lons,HI4PI.lats)
lims_lb = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')
ax.plot_surface(x,y,whi_raw)
ax.view_init(azim=70,elev=10)


x,y = np.meshgrid(HI4PI.lons,HI4PI.lats)


fig = plt.figure(figsize=(10,10))
ax2 = fig.add_subplot(projection='3d')
ax2.plot_surface(x,y,whi_smooth)
ax2.view_init(azim=70,elev=10)


fig = plt.figure(figsize=(10,10))
ax2 = fig.add_subplot(projection='3d')
ax2.plot_surface(x,y,wsmo)
ax2.view_init(azim=70,elev=10)


smo = fits.open('smoothed_HI4PI_sig40.fits')[0].data


wsmo = np.nansum(smo,axis=0)


fig = plt.figure(figsize=(10,10))
ax3 = fig.add_subplot(projection='3d')
ax3.plot_surface(x,y,whi_raw-whi_smooth)
ax3.view_init(azim=70,elev=10)


plt.figure(figsize=(10,10))
plt.imshow(whi_raw-whi_smooth, origin='lower',cmap=dc.cmap_grey0('viridis'),norm='log')#,vmax=0,vmin=-1000)
plt.colorbar(orientation='horizontal')


cont_sub_cube = HI4PI.cube-HI4PI_smooth


plt.figure(figsize=(10,10))
plt.imshow(np.nansum(cont_sub_cube,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),norm='log')


plt.figure(figsize=(10,10))
plt.imshow(np.nansum(cont_sub_cube,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),vmin=0)
plt.colorbar(orientation='horizontal')


plt.figure(figsize=(10,10))
for i in range(0,310,10):
    plt.pcolormesh(cont_sub_cube[i],cmap=dc.cmap_grey0('viridis'))
    plt.show()
#     plt.colorbar(orientation='horizontal')


plt.figure(figsize=(10,10))
plt.imshow(np.nanmedian(cont_sub_cube,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims_lb)
plt.colorbar(orientation='horizontal')


cont_sub_cube.shape


DHT_CO.cube.shape


(76*481*2881/(310*1200*4323))*100


100/7


np.nanstd(DHT_CO.cube)


np.nanmin(DHT_CO.cube)


np.cbrt(500)


DHT_CO.db*8, DHT_CO.dl*8, DHT_CO.dv*8


1/HI4PI.db,1/HI4PI.dl,10/(HI4PI.dv.to(u.km/u.s))


12*12*8


np.nanstd(DHT_CO.cube),np.nanmin(DHT_CO.cube),np.nanmax(DHT_CO.cube)


20*0.3


np.nanstd(cont_sub_cube),np.nanmin(cont_sub_cube),np.nanmax(cont_sub_cube)


-64+18
