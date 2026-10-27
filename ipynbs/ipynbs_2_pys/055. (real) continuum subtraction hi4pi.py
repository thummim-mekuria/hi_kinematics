#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datacube as dc
from astropy.modeling import models,fitting
import astropy.units as u
from specutils.fitting import fit_generic_continuum
from specutils import Spectrum1D


nice_cat = pd.read_csv('lfs_w_simple_spectra_nice_r_acc.csv')
mek = pd.read_csv('Mekuria_catalog.csv',skiprows=2)


hi4pi = fits.open('HI4PI.fits')[0]


dv,db,dl = hi4pi.header['CDELT1'],hi4pi.header['CDELT2'],hi4pi.header['CDELT3']


v0,b0,l0 = hi4pi.header['CRVAL1'],hi4pi.header['CRVAL2'],hi4pi.header['CRVAL3']
l0 = l0*-1


nv, nb, nl = hi4pi.header['NAXIS3'],hi4pi.header['NAXIS2'],hi4pi.header['NAXIS1']


hi_vels = dv*np.arange(nv)+v0
hi_lats = db*np.arange(nb)+b0
hi_lons = dl*np.arange(nl)+l0


bspec = np.nansum(hi4pi.data,axis=(0,2))


smooth_hi = fits.open('smoothed_HI4PI_sig40.fits')[0].data


bspec_sm = np.nansum(smooth_hi,axis=(0,2))


plt.plot(hi_lats,bspec/np.nanmax(bspec))
plt.plot(hi_lats, bspec_sm/np.nanmax(bspec))
plt.plot(hi_lats, bspec/bspec_sm)


for i in range(len(nice_cat)):
    sigv = round(nice_cat['velocity_dispersion_parent<km/s>'][i]/dv)
    ra = round(nice_cat['acc_radius_parent<deg>'][i]/db)
    l,b,v = nice_cat['x_cen_hi4pi<px>'][i],nice_cat['y_cen_hi4pi<px>'][i],nice_cat['v_cen_hi4pi<px>'][i]
    
    sb = hi4pi.data[v-sigv:v+sigv,b-ra:b+ra,l-ra:l+ra]
    sb_smooth = smooth_hi[v-sigv:v+sigv,b-ra:b+ra,l-ra:l+ra]
    
    
    bspc = np.nansum(sb,axis=(0,2))
    bspc_sm = np.nansum(sb_smooth,axis=(0,2))
    
    x,y = np.meshgrid(hi_lons[l-ra:l+ra],hi_lats[b-ra:b+ra])
    w = np.nansum(sb, axis=0)
    www = np.nansum(sb_smooth,axis=0)
    
# #     #--------------3d / side views------------
#     fig = plt.figure(figsize=(7,7))
# #     ax = fig.add_subplot(1,3,1,projection='3d')
#     ax2 = fig.add_subplot(1,2,2,projection='3d')
#     ax3 = fig.add_subplot(1,2,1,projection='3d')#133
    
# #     ax.plot_wireframe(x,y,w)
# #     ax.view_init(azim=10,elev=5)
    
#     ax2.plot_wireframe(x,y,w/np.nanmax(w) - www/np.nanmax(www))
#     ax2.view_init(azim=0,elev=2)
    
# #     ax3.plot_wireframe(x,y,np.nansum(sb-sb_smooth,axis=0))
#     ax3.plot_wireframe(x,y,w/np.nanmax(w),color='orange',alpha=0.5)
#     ax3.plot_wireframe(x,y,www/np.nanmax(www),alpha=0.5)
#     ax3.view_init(azim=0,elev=2)
    
#     plt.show()
# #     #-----------------------------------------
  
    
#     #------------spectra-----------------------------
    plt.plot(hi_lats[b-ra:b+ra],bspc)
    plt.plot(hi_lats[b-ra:b+ra],bspc_sm)
    flat = bspc_sm/np.nanmax(bspc_sm)
    plt.plot(hi_lats[b-ra:b+ra],(bspc-bspc_sm))
    plt.plot(hi_lats[b-ra:b+ra],(bspc/flat))
    plt.show()
#     #-----------------------------------------


whi = np.nansum(hi4pi.data,axis=0)


wsm = np.nansum(smooth_hi,axis=0)


flat_whole = wsm/np.nanmax(wsm)


plt.imshow(whi,origin='lower',cmap=dc.cmap_grey0('viridis'),norm='log')


plt.imshow(flat_whole,origin='lower',cmap=dc.cmap_grey0('viridis'),norm='log')


plt.imshow(whi-wsm,origin='lower',cmap=dc.cmap_grey0('viridis'),vmin=-3e2,vmax=10)#,norm='log')
plt.colorbar(orientation='horizontal')


plt.imshow(whi/flat_whole,origin='lower',cmap=dc.cmap_grey0('viridis'),vmin=-3e2,vmax=4e3)#,norm='log')
plt.colorbar(orientation='horizontal')


data_sm = smooth_hi.flatten()


plt.hist(data_sm,bins=1000)
plt.yscale('log')
plt.show()


neg = hi4pi.data <=0


sigg043 = hi4pi.data >= 0.043 #rms noise from announcement paper


sigg5x043 = hi4pi.data >= 5*0.043


sigg3x043 = hi4pi.data >= 3*0.043


for i in range(0,310,10):
    plt.imshow(neg[i])
    plt.title(str(i))
    plt.show()


for i in range(0,310,10):
    plt.imshow(sigg043[i])
    plt.title(str(i))
    plt.show()


for i in range(0,310,10):
    plt.imshow(sigg3x043[i])
    plt.title(str(i))
    plt.show()


for i in range(0,310,10):
    plt.imshow(sigg5x043[i])
    plt.title(str(i))
    plt.show()


hdu=fits.PrimaryHDU(hi_smooth_043)


hdu.writeto('smoothed_HI4PI_sig40_rms043.fits')


flat043 = hi_smooth_043/(np.nanmax(hi_smooth_043))


adj = np.abs(sigg043*hi4pi.data)/flat043


hdu1 = fits.PrimaryHDU(adj)
hdu1.writeto('HI4PI_flatcalib.fits')


w043_adj = np.nansum(adj,axis=0)


wflat = np.nansum(flat043,axis=0)


plt.imshow(whi,origin='lower',cmap=dc.cmap_grey0('viridis'),vmin=0)#,vmax=4e3,vmin=1e3)
plt.colorbar(orientation='horizontal')


plt.imshow(wflat,origin='lower',cmap=dc.cmap_grey0('viridis'))#,vmax=4e3,vmin=1e3)
plt.colorbar(orientation='horizontal')


x,y = np.meshgrid(hi_lons,hi_lats)


vals = 


vals = np.array((3e3,5e3,9e3))
lims = [hi_lons[0],hi_lons[-1],hi_lats[0],hi_lats[-1]]
plt.figure(figsize=(20,5))
plt.imshow(w043_adj,origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims,vmax=1e4)
plt.colorbar(orientation='horizontal')
plt.contour(x,y,w043_adj,levels=vals,linewidths=vals/1e4,colors=['black','red','white'])#,alpha=0.5)
plt.show()


from astrodendro import Dendrogram


dnd_w043_adj = Dendrogram.compute(w043_adj)


p043 = dnd_w043_adj.plotter()


fig, ax = plt.subplots(1,1,figsize=(7,7))
p043.plot_tree(ax)


plt.figure(figsize=(20,5))
plt.imshow(dnd_w043_adj.index_map,origin='lower',cmap=dc.cmap_grey0('viridis'))
plt.show()


wsm_043 = np.nansum(hi4pi.data-hi_smooth_043,axis=0)


dnd_w043_sub = Dendrogram.compute(wsm_043)


plt.figure(figsize=(20,5))
plt.imshow(dnd_w043_sub.index_map,origin='lower',cmap=dc.cmap_grey0('viridis'))
plt.show()


#when coming back, see if/where there are zeros in the flat, 
#bc that seems to be throwing a divide by 0 error
#email nia asking if it makes more sense (from a data reduction standpoint)
#to divide or to subtract
#also say, i have narrowed down the catalog based on the physical properties
#we are selecting for (one bump in spectrum etc) and (not sure if to phrase this
#as asking for permission, or telling her this is what made more sense to me)
#send email by end of the day today


np.nanmean(adj),np.nanmedian(adj)


np.nanstd(adj)


adj = fits.open('HI4PI_flatcalib.fits')[0].data


dnd_calib_sig = Dendrogram.compute(adj)


plt.hist(nice_cat['acc_radius_parent<deg>'],bins=20)
plt.show()


np.min(nice_cat['acc_radius_parent<deg>'].to_numpy())


np.min(nice_cat['velocity_dispersion_parent<km/s>'].to_numpy())


1.86/db


2.63/dv


2*(22**2)


adj.size


(1000/1608156000)*100
