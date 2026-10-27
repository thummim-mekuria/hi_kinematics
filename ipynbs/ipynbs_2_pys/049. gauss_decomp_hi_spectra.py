#!/usr/bin/env python
# coding: utf-8


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
from astrodendro import Dendrogram,Dendrogram_clobby


DHT_CO = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')
DHT_CO.cube = DHT_CO.cube[209:285,:,:] #trim to narrower velocity range
DHT_CO.vels = DHT_CO.vels[209:285]
#prob also need to take out the center? not sure if that matters yet tho


cat = pd.read_csv('physcl_prprts_catalog.csv')
dnd = Dendrogram.load_from('COGAL_dendrogram.fits')


sigs = np.zeros(len(cat))
specs = np.zeros((len(cat),len(DHT_CO.vels)))

for i,idx in enumerate(cat['_idx']):
    cube = dnd[idx].get_mask() * dnd.data
    spc = dc.spectrum1D(cube)
    sig = dc.dispersion_1sigma(spc)
    
    if len(sig)==0 or (np.argmax(spc) < 5*sig) or (np.argmax(spc) > len(spc) - 5*sig):
        sigs[i] = np.nan
        specs[i] = np.nan
    else:
        sigs[i] = sig[0]*1.3
        specs[i] = spc 
    
cat['sigma<km/s>'] = sigs


good_specs = cat.dropna(subset='sigma<km/s>')


vrms3x = 3*good_specs['sigma<km/s>'].to_numpy()


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]


ls = good_specs['lon<deg>'].to_numpy()#*u.deg
bs = good_specs['lat<deg>'].to_numpy()#*u.deg
vs = good_specs['vel<km/s>'].to_numpy()#*u.km/u.s 
rs = good_specs['accumulation_radius<deg>'].to_numpy()#*u.deg 
vrms = vrms3x#*u.km/u.s

parents = cld.accumulation_prisms(vs,bs,ls,rs,vrms,HI4PI)


def gaussian(xs, peaks, FWHMs, means):
    """Return a Gaussian function."""
    gs = []
    for p,fwhm,mn in zip(peaks,FWHMs,means):
        sigma = fwhm / 2.354820045  # (2 * sqrt( 2 * ln(2)))
        g = p * np.exp(-(xs - mn)**2 / 2. / sigma**2)
        gs.append(g)
        
    return gs


def gauss_fit(spec,thresh=None):
    from scipy.signal import find_peaks,peak_widths
    if thresh==None:
        thresh=np.nanstd(spec)
        
    N = len(spec)
    xs = np.arange(N)
    
    fnd_pk_out = find_peaks(spec,height=thresh)
    pk_xs = fnd_pk_out[0] #peak positions along x
    pk_ys = fnd_pk_out[1]['peak_heights'] #peak positions along y
    
    pk_wdth_out = peak_widths(spec,pk_xs)
    wdths = pk_wdth_out[0] #fwhms
    
    decomposed = gaussian(xs,pk_ys,wdths,pk_xs)
    decomposed = np.array(decomposed)
    
    nump_comps = len(decomposed)
    
    return decomposed


plt.plot(sum_fits)


for i,idx in enumerate(cat['_idx']):
    if i in good_specs.index:
        l_c = good_specs['lon<deg>'].loc[i]
        b_c = good_specs['lat<deg>'].loc[i]
        ra = good_specs['accumulation_radius<deg>'].loc[i]
        
        lat_inds = cld.ind_ang('get index', HI4PI.lats, (b_c-ra,b_c+ra))
        lon_inds = cld.ind_ang('get index', HI4PI.lons, (l_c-ra,l_c+ra))
        
        lat_inds.sort()
        lon_inds.sort()

        b0, b1 = lat_inds #get each element separately
        l0, l1 = lon_inds

        sub_hi = HI4PI.cube[:,b0:b1,l0:l1]
        hi_spec = dc.spectrum1D(sub_hi)
        
        fits = gauss_fit(hi_spec)
        fits = np.array(fits)
        
#         sum_fits = np.nansum(fits,axis=0)
        plt.figure(figsize=(10,10))
        for g in fits:
            plt.plot(HI4PI.vels/1e3,g/np.nanmax(g),linestyle=':')
    
        
        plt.plot(DHT_CO.vels,specs[i]/np.nanmax(specs[i]),color='black',alpha=0.6)
#         plt.plot(HI4PI.vels/1e3, hi_spec/(np.nanmax(hi_spec)),color='black',alpha=0.5)
#         plt.plot(HI4PI.vels/1e3, sum_fits/(np.nanmax(sum_fits)),color='grey',linestyle=':')
        plt.title(str(idx))
        plt.xlim(-50,50)
        plt.show()


for i,idx in enumerate(cat['_idx']):
    if i in good_specs.index:
        l_c = good_specs['lon<deg>'].loc[i]
        b_c = good_specs['lat<deg>'].loc[i]
        ra = good_specs['accumulation_radius<deg>'].loc[i]
        
        lat_inds = cld.ind_ang('get index', HI4PI.lats, (b_c-ra,b_c+ra))
        lon_inds = cld.ind_ang('get index', HI4PI.lons, (l_c-ra,l_c+ra))
        
        lat_inds.sort()
        lon_inds.sort()

        b0, b1 = lat_inds #get each element separately
        l0, l1 = lon_inds

        sub_hi = HI4PI.cube[:,b0:b1,l0:l1]
        hi_spec = dc.spectrum1D(sub_hi)
        
        fits = gauss_fit(hi_spec)
        fits = np.array(fits)
        
#         sum_fits = np.nansum(fits,axis=0)
        plt.figure(figsize=(10,10))
        for g in fits:
            plt.plot(HI4PI.vels/1e3,g/np.nanmax(g),linestyle=':')
    
        
        plt.plot(DHT_CO.vels,specs[i]/np.nanmax(specs[i]),color='black',alpha=0.6)
#         plt.plot(HI4PI.vels/1e3, hi_spec/(np.nanmax(hi_spec)),color='black',alpha=0.5)
#         plt.plot(HI4PI.vels/1e3, sum_fits/(np.nanmax(sum_fits)),color='grey',linestyle=':')
        plt.title(str(idx))
        plt.xlim(-50,50)
        plt.show()


HI4PI.dl


good_specs


for i,idx in enumerate(cat['_idx'].loc[:20]):
    
    if i in good_specs.index:
#         print(i,idx)
        l_c = good_specs['lon<deg>'].loc[i]
        b_c = good_specs['lat<deg>'].loc[i]
        ra = good_specs['accumulation_radius<deg>'].loc[i]
        
        lat_inds = cld.ind_ang('get index', HI4PI.lats, (b_c-ra,b_c+ra))
        lon_inds = cld.ind_ang('get index', HI4PI.lons, (l_c-ra,l_c+ra))
        
        lat_inds.sort()
        lon_inds.sort()

        b0, b1 = lat_inds #get each element separately
        l0, l1 = lon_inds

        sub_hi = HI4PI.cube[:,b0:b1,l0:l1]
        hi_spec = dc.spectrum1D(sub_hi)
        
        lims_ = [(HI4PI.vels/1e3)[0],(HI4PI.vels/1e3)[-1],HI4PI.lons[l0]*100,HI4PI.lons[l1]*100]
#         lims_ = [HI4PI.lons[l0]*100,HI4PI.lons[l1]*100,(HI4PI.vels/1e3)[0],(HI4PI.vels/1e3)[-1]]

#         plt.figure(figsize=(4,20))

        whi_1 = np.nansum(sub_hi,axis=1)
        plt.imshow(np.transpose(whi_1),extent=lims_)
#         plt.xticks(np.arange(HI4PI.lons[l0],HI4PI.lons[l1],0.1))
        
        fitses = gauss_fit(hi_spec)
        fitses = np.array(fitses)
        
        sum_fits = np.nansum(fitses,axis=0)
#         plt.figure(figsize=(10,10))
        y0,ymx = HI4PI.lons[l0]*100,HI4PI.lons[l1]*100
#         for g in fitses:
        y_this = np.nanargmax(sum_fits)
#             plt.hlines((HI4PI.vels/1e3)[y_this],HI4PI.lons[l0]*100,HI4PI.lons[l1]*100,color='red')
        plty = sum_fits/np.nanmax(sum_fits)
            
        
        plt.plot(HI4PI.vels/1e3,y0 + (ymx-y0)*plty,linestyle=':',color='white')
#         plt.show()
#         plt.figure(figsize=(10,10))
#         for f in fitses:
#             plt.plot(HI4PI.vels/1e3,g/np.nanmax(g),linestyle=':')

            
#         plt.plot(DHT_CO.vels,specs[i]/np.nanmax(specs[i]),color='black',alpha=0.6)
        plt.plot(HI4PI.vels/1e3, (ymx-y0)*(hi_spec/(np.nanmax(hi_spec))) + y0,color='orange')#,alpha=0.5)
# # #         plt.plot(HI4PI.vels/1e3, sum_fits/(np.nanmax(sum_fits)),color='grey',linestyle=':')
# #         plt.title(str(idx))
#         plt.xlim(-50,50)
        plt.ylim(y0,ymx)
        plt.show()


help(plt.contour)


dir(ax)


for i,idx in enumerate(cat['_idx'].loc[:20]):
    
    if i in good_specs.index:
#         print(i,idx)
        l_c = good_specs['lon<deg>'].loc[i]
        b_c = good_specs['lat<deg>'].loc[i]
        ra = good_specs['accumulation_radius<deg>'].loc[i]
        
        lat_inds = cld.ind_ang('get index', HI4PI.lats, (b_c-ra,b_c+ra))
        lon_inds = cld.ind_ang('get index', HI4PI.lons, (l_c-ra,l_c+ra))
        
        lat_inds.sort()
        lon_inds.sort()

        b0, b1 = lat_inds #get each element separately
        l0, l1 = lon_inds

        sub_hi = HI4PI.cube[:,b0:b1,l0:l1]
        hi_spec = dc.spectrum1D(sub_hi)
        
        fig = plt.figure(figsize=(10,7))
        ax1 = fig.add_subplot(131,projection='3d')
        ax2 = fig.add_subplot(132,projection='3d')
        ax3 = fig.add_subplot(133,projection='3d')
        
        lims = [HI4PI.lons[l0],HI4PI.lons[l1],HI4PI.lats[b0],HI4PI.lats[b1]]
        whi_0 = np.nansum(sub_hi,axis=0)
        smoof = ndimage.gaussian_filter(whi_0,sigma=5)
        x,y = np.meshgrid(HI4PI.lons[l0:l1],HI4PI.lats[b0:b1])

        mx1 = np.nanmax(np.abs(whi_0))
        ax1.plot_wireframe(x,y,whi_0,alpha=0.5)
        ax1.view_init(azim=0,elev=40)
        ax1.contour(x,y,whi_0,levels=[0.1*mx1,0.7*mx1,mx1],colors=['red'])


        mx2 = np.nanmax(np.abs(whi_0-smoof))
        ax2.plot_wireframe(x,y,whi_0-smoof,alpha=0.5)
        ax2.view_init(azim=0,elev=40)
        ax2.contour(x,y,whi_0-smoof,levels=[0.1*mx2,0.7*mx2,mx2],colors=['red'])
        
        ax3.plot_wireframe(x,y,smoof,alpha=0.3,color='orange')
        ax3.plot_wireframe(x,y,3*(whi_0-smoof),alpha=0.3,color='green')

        ax3.plot_wireframe(x,y,whi_0-(0.3*np.max(whi_0)),alpha=0.2)
        ax3.view_init(azim=0,elev=0)


        fitses = gauss_fit(hi_spec)
        fitses = np.array(fitses)
        
        plt.show()


plt.plot(np.transpose(whi_0))
plt.show()


plt.plot(whi_0)
plt.show()


plt.plot(np.transpose(whi_0),color='green')
plt.plot(whi_0,color='red')
plt.show()


plt.plot(np.nanmean(whi_0,axis=1))
plt.plot(np.nanmean(whi_0,axis=0))


plt.plot(gauss_fit(np.nanmean(whi_0,axis=1))[0])


smoof_3d_sub_hi = ndimage.gaussian_filter(sub_hi,sigma=5)


flat_subhi = sub_hi-smoof_3d_sub_hi


np.nanstd(smoof_3d_sub_hi)


subhid = Dendrogram.compute(flat_subhi,min_value=17,min_delta=8,min_npix=20)


len(subhid.trunk)


colors=[]


from matplotlib import colormaps


prd = colormaps.get_cmap('Paired').colors
vrds = colormaps.get_cmap('viridis').colors


for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()

    w = np.nansum(m,axis=0)
    plt.imshow(w,cmap=dc.cmap_grey0('viridis'),origin='lower')
    
    plt.show()


ws = []
for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()
#     print(np.nanmax(m),np.nanmin(m))
#     print(m.shape)
    w = np.nansum(m,axis=1)
    this = w.astype('bool')* (t+1)
    ws.append(this)
    
for t,w in enumerate(ws):
    plt.contour(x,y,w,colors=[prd[t]],levels=t+1,linewidths=0.5)
plt.show()


ws = []
for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()
#     print(np.nanmax(m),np.nanmin(m))
#     print(m.shape)
    w = np.nansum(m,axis=0)
    this = w.astype('bool')* (t+1)
    ws.append(this)
    
for t,w in enumerate(ws):
    plt.contour(x,y,w,colors=[prd[t]],levels=t+1,linewidths=0.5)
plt.show()


np.nanmin(tr_spec)


tr_specs = [] 
for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()
    tr_spec = dc.spectrum1D(m*sub_hi)
    tr_specs.append(tr_spec)
    plt.plot(np.log10(tr_spec+1))
plt.plot(np.log10(hi_spec+1))
# plt.plot(np.log10(1+np.nansum(tr_specs,axis=0)))
plt.show()


tr_specs = [] 
for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()
    tr_spec = dc.spectrum1D(m*sub_hi)
    tr_specs.append(tr_spec)
#     plt.plot(np.log10(tr_spec+1))
plt.plot(np.log10(hi_spec+1))
plt.plot(np.log10(1+np.nansum(tr_specs,axis=0)))
plt.show()


import cc3d
parents_masked_overlap,N = cc3d.connected_components(parents,connectivity=26,return_N = True)
#connectivity=6 and 18 also give the same N


from numba import njit

@njit
def separate_clouds(labels_map,data,N):
    x,y,z = data.shape
    out = np.zeros((N,x,y,z))
    for n in range(1,N+1):
        m = (labels_map==n)
        mskd = m*data
        out[n-1] = mskd
        
    return out


# #come back and figure out whats wrong with this / the function ^ 
# distinct_hi = separate_clouds(parents_masked_overlap,HI4PI.cube,N)


# lims = [HI4PI.lons[0], HI4PI.lons[-1], HI4PI.lats[0], HI4PI.lats[-1]]
# for i in range(N):
#     plt.imshow(np.nansum(HI4PI.cube*distinct_hi[i],axis=0),\
#                cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)
# #     plt.colorbar()
#     plt.title(str(i+1))
#     plt.show()


# lims = [HI4PI.lons[0], HI4PI.lons[-1], 1e-3*HI4PI.vels[0], 1e-3*HI4PI.vels[-1]]
# for i in range(N):
#     plt.imshow(np.nansum(HI4PI.cube*distinct_hi[i],axis=1),\
#                cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)
# #     plt.colorbar()
#     plt.title(str(i+1))
#     plt.show()


# #once i figure out getting distinct_hi, save it as a fits
# HI_parents = fits.PrimaryHDU(distinct_hi)
# HI_parents.writeto('decomp_HI_parents_frm_RA.fits',overwrite=True)
