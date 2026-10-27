#!/usr/bin/env python
# coding: utf-8


from IPython.core.display import display, HTML
display(HTML("<style>.container { width:85% !important; }</style>"))


import datacube as dc
import clouds as cld
import interpolate as ip
import bootleg as btlg

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage

# %matplotlib ipympl
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]

x,y = np.meshgrid(HI4PI.lons, HI4PI.lats)


zero = np.nansum(HI4PI.cube,axis=0)


minus2 = zero-ndimage.gaussian_filter(zero,sigma=2)
minus5 = zero-ndimage.gaussian_filter(zero,sigma=5)
minus20 = zero-ndimage.gaussian_filter(zero,sigma=20)
minus50 = zero-ndimage.gaussian_filter(zero,sigma=50)
minus100 = zero-ndimage.gaussian_filter(zero,sigma=100)


fig = plt.figure(figsize=(20,20))
x,y = np.meshgrid(HI4PI.lons, HI4PI.lats)


ax1 = fig.add_subplot(2,2,1,projection='3d')
ax1.plot_surface(x,y,ndimage.gaussian_filter(zero,sigma=2),alpha=0.7)
ax1.view_init(2,0,0)


ax2 = fig.add_subplot(2,2,2,projection='3d')
ax2.plot_surface(x,y,ndimage.gaussian_filter(zero,sigma=20),alpha=0.7)
ax2.view_init(2,0,0)


ax3 = fig.add_subplot(2,2,3,projection='3d')
ax3.plot_surface(x,y,ndimage.gaussian_filter(zero,sigma=50),alpha=0.7)
ax3.view_init(2,60,0)

ax4 = fig.add_subplot(2,2,4,projection='3d')
ax4.plot_surface(x,y,ndimage.gaussian_filter(zero,sigma=100),alpha=0.7)
ax4.view_init(2,60,0)

# plt.show()


fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(projection='3d')


x,y = np.meshgrid(HI4PI.lons, HI4PI.lats)
ax.plot_surface(x,y,np.abs(minus20),alpha=0.2)
# ax.contour(x,y,minus100,levels=[150],colors='red',linewidths=0.5)
ax.view_init(15,150,0)
plt.show()


ormaybe = HI4PI.cube[:,:,3900:]


b_spec_ormaybe = np.swapaxes(ormaybe,1,2)


b_spec_ormaybe.shape


plt.figure(figsize=(20,10))
lims = [HI4PI.lats[0],HI4PI.lats[-1],HI4PI.lons[3900],HI4PI.lons[4322]]
plt.imshow(np.nansum(b_spec_ormaybe,axis=0),extent=lims,cmap=dc.cmap_grey0('Paired'),vmax=900)
plt.colorbar()


plt.figure(figsize=(20,10))
lims = [1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1],HI4PI.lats[0],HI4PI.lats[-1]]
plt.imshow(np.nansum(b_spec_ormaybe,axis=1),extent=lims)


plt.figure(figsize=(20,10))
lims = [1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1],HI4PI.lons[3900],HI4PI.lons[-1]]
plt.imshow(np.nansum(b_spec_ormaybe,axis=2),extent=lims)


310*423


bspecs = b_spec_ormaybe.reshape(131130,1200)


cont_sub = np.zeros_like(bspecs)


for i,b in enumerate(bspecs):
    cont_sub[i] = ndimage.gaussian_filter1d(b,sigma=20)


cont_sub_or = cont_sub.reshape(310,423,1200)


plt.figure(figsize=(20,10))
lims = [HI4PI.lats[0],HI4PI.lats[-1],HI4PI.lons[3900],HI4PI.lons[4322]]
plt.imshow(np.nansum(b_spec_ormaybe-cont_sub_or,axis=0),extent=lims,cmap=dc.cmap_grey0('Paired'),vmin=-200,vmax=200)
plt.colorbar()


def b_cont_sub(cube,sigma=20):
    swapped = np.swapaxes(cube,1,2)
    x,y,z = swapped.shape
    flattened = swapped.reshape(x*y,z)
    
    flat_cont_subbed = np.zeros_like(flattened)
    for i,b in enumerate(flattened):
        flat_cont_subbed[i] = ndimage.gaussian_filter1d(b,sigma)
        
    cont_sub = flat_cont_subbed.reshape(x,y,z)
    
    out = np.swapaxes(cont_sub,1,2)
    
    return out


HI4PI_wo_bcont = b_cont_sub(HI4PI.cube)


plt.imshow(np.nansum(HI4PI_wo_bcont,axis=0),origin='lower',cmap=dc.cmap_grey0('Paired'),vmax=600)
# plt.colorbar()
plt.show()


zero.shape


DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')

DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 

smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole)
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641]
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]

DHT_CO.cube = wo_gal_cent


lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


lvs_0 = np.nansum(lvs,axis=(0,1))


DHT_CO.lons


minus200 = zero-ndimage.gaussian_filter(zero,sigma=200)


minus35 = zero-ndimage.gaussian_filter(zero,sigma=35)


plt.figure(figsize=(8,15))
lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
plt.contour(HI4PI.lons,HI4PI.lats,minus200,levels=[-50,-20,20],colors='black',linewidths=[0.1,0.3,0.5])
# plt.contour(HI4PI.lons,HI4PI.lats,minus20,levels=[20,50,150],colors='blue',linewidths=[0.7,0.5,0.3])

# plt.xlim(180,-180)
plt.xlim(-130,-180)

plt.contour(DHT_CO.lons,DHT_CO.lats,lvs_0,levels=[1.5],colors='red',linewidths=[0.7])
# plt.colorbar()
plt.show()


fig,axs = plt.subplots(1,3,figsize=(15,10))
lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]


axs[0].contourf(HI4PI.lons,HI4PI.lats,minus20,levels=[-300,-150,0,150,300],\
                colors=['xkcd:cerulean','xkcd:baby blue','xkcd:sage green','xkcd:light yellow'])#,linewidths=[0.1,0.5,0.1])
axs[0].set_xlim(-130,-180)
axs[0].contour(DHT_CO.lons,DHT_CO.lats,lvs_0,levels=[2],colors='red',linewidths=[0.7])


axs[1].contourf(HI4PI.lons,HI4PI.lats,minus35,levels=[-300,-150,0,150,300],\
                colors=['xkcd:cerulean','xkcd:baby blue','xkcd:sage green','xkcd:light yellow'])
axs[1].set_xlim(-130,-180)
axs[1].contour(DHT_CO.lons,DHT_CO.lats,lvs_0,levels=[2],colors='red',linewidths=[0.7])


axs[2].contourf(HI4PI.lons,HI4PI.lats,np.sqrt(np.abs(zero)),levels=[10,13,16,20,23],\
                colors=['xkcd:cerulean','xkcd:baby blue','xkcd:sage green','xkcd:light yellow'])
axs[2].set_xlim(-130,-180)
axs[2].contour(DHT_CO.lons,DHT_CO.lats,lvs_0,levels=[2],colors='red',linewidths=[0.7])

plt.show()


plt.figure(figsize=(20,5))
lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
# plt.contour(HI4PI.lons,HI4PI.lats,minus50,levels=[20,50,150],colors='black',linewidths=[0.3,0.5,0.7])
plt.contour(HI4PI.lons,HI4PI.lats,minus100,levels=[80],colors='black',linewidths=[0.5])

plt.xlim(180,-180)
plt.contour(DHT_CO.lons,DHT_CO.lats,lvs_0,levels=[2],colors='red',linewidths=[0.3])
# plt.colorbar()
plt.show()


plt.imshow(np.nansum(HI4PI.cube,axis=0),origin='lower',cmap=dc.cmap_grey0('Paired'),vmin=50,vmax=700)
plt.xlim(3900,4323)

# plt.colorbar()
plt.show()


plt.imshow(np.nansum(HI4PI.cube - HI4PI_wo_bcont,axis=0),origin='lower',cmap=dc.cmap_grey0('Paired'),vmin=50,vmax=700)
plt.xlim(3900,4323)

# plt.colorbar()
plt.show()


def l_cont_sub(cube,sigma=20):
    x,y,z = cube.shape
    flattened = cube.reshape(x*y,z)
    
    flat_cont_subbed = np.zeros_like(flattened)
    for i,l in enumerate(flattened):
        flat_cont_subbed[i] = ndimage.gaussian_filter1d(l,sigma)
        
    out = flat_cont_subbed.reshape(x,y,z)
        
    return out


HI4PI_wo_lcont = l_cont_sub(HI4PI.cube)


plt.imshow(np.nansum(HI4PI_wo_lcont,axis=0),origin='lower',cmap=dc.cmap_grey0('Paired'),vmax=600)
# plt.colorbar()
plt.show()


plt.imshow(np.nansum(HI4PI.cube - HI4PI_wo_lcont,axis=0),origin='lower',cmap=dc.cmap_grey0('Paired'),vmin=50,vmax=700)
plt.xlim(3900,4323)
# plt.colorbar()
plt.show()


plt.imshow(minus20,origin='lower',cmap=dc.cmap_grey0('Paired'),vmin=50,vmax=700)
plt.xlim(3900,4323)
# plt.colorbar()
plt.show()


cocat = pd.read_csv('co_catalog.csv',delim_whitespace=True,skiprows=32)


plt.hist(cocat['R_A'],bins=50)


# plt.hist(cocat['major_sigma'],bins=50)
# plt.hist(cocat['minor_sigma'],bins=50)

np.nanmean(np.sqrt(cocat['area_ellipse']/np.pi))


3*0.8/0.125
