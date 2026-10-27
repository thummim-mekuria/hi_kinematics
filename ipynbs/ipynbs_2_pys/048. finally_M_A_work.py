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


cat


cat[np.isnan(cat['velocity_dispersion<km/s>'])]


#can probably add back idx = 13, 18, 34, 37 ... ie decrease the width cutoff from 5sig to 3sig?


# sigs = np.zeros(len(cat))
# specs = []

for i,idx in enumerate(cat['_idx']):
    cube = dnd[idx].get_mask() * dnd.data
    spc = dc.spectrum1D(cube)
    sig = dc.dispersion_1sigma(spc)
    
    
    if len(sig)==0 or (np.argmax(spc) < 4*sig) or (np.argmax(spc) > len(spc) - 4*sig):
        print(np.argmax(spc),4*sig)
        plt.plot(spc)
        plt.title(str(idx))
        plt.show()
#         sigs[i] = np.nan
#     else:
#         sigs[i] = sig[0]*1.3
        
# #     specs[i] = spc        
# cat['sigma<km/s>'] = sigs
# cat['vel_spectrum'] = specs


# sigs = np.zeros(len(cat))
# specs = []

for i,idx in enumerate(cat['_idx']):
    cube = dnd[idx].get_mask() * dnd.data
    spc = dc.spectrum1D(cube)
    sig = dc.dispersion_1sigma(spc)
    
    
    if len(sig)==0 or (np.argmax(spc) < 5*sig) or (np.argmax(spc) > len(spc) - 5*sig):
        plt.plot(spc)
        plt.title(str(idx))
        plt.show()
#         sigs[i] = np.nan
#     else:
#         sigs[i] = sig[0]*1.3
        
# #     specs[i] = spc        
# cat['sigma<km/s>'] = sigs
# cat['vel_spectrum'] = specs


sigs = np.zeros(len(cat))
specs = np.zeros((len(cat),len(DHT_CO.vels)))

for i,idx in enumerate(cat['_idx']):
    cube = dnd[idx].get_mask() * dnd.data
#     print(cube.shape)
    spc = dc.spectrum1D(cube)
    sig = dc.dispersion_1sigma(spc)
    
    if len(sig)==0 or (np.argmax(spc) < 5*sig) or (np.argmax(spc) > len(spc) - 5*sig):
        sigs[i] = np.nan
        specs[i] = np.nan
    else:
        sigs[i] = sig[0]*1.3
        specs[i] = spc 
    
cat['sigma<km/s>'] = sigs
# cat['vel_spectrum'] = specs


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


np.all(parents==False)


whi = np.nansum(HI4PI.cube,axis=0)


good_specs


fig,axs = plt.subplots(2,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
axs[0].imshow(whi,extent=lims,origin='lower')


for l,b,rbb,r_big in \
zip(cat['lon<deg>'],cat['lat<deg>'],\
    cat['effective_radius<deg>'],\
    cat['accumulation_radius<deg>']):    
    axs[0].add_patch(Circle((l,b),rbb,fill=False,edgecolor='white'))
    axs[0].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red'))

axs[0].set_xlim(180,-5)


axs[1].imshow(whi,extent=lims,origin='lower')
for l,b,rbb,r_big in \
zip(cat['lon<deg>'],cat['lat<deg>'],\
    cat['effective_radius<deg>']\
    ,cat['accumulation_radius<deg>']):    
    axs[1].add_patch(Circle((l,b),rbb,fill=False,edgecolor='white'))
    axs[1].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red'))

axs[1].set_xlim(5,-180)


wco = np.nansum(np.abs(DHT_CO.cube), axis=0)


cat['effective_radius<deg>'] = np.sqrt(cat['area_sky<deg2>']/np.pi)


fig,axs = plt.subplots(2,1,figsize=(15,10))

lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]
axs[0].imshow(wco,extent=lims,origin='lower',cmap = dc.cmap_grey0('viridis'),vmax=25)


for l,b,rbb,r_big in \
zip(cat['lon<deg>'],cat['lat<deg>'],\
    cat['effective_radius<deg>']\
    ,cat['accumulation_radius<deg>']):    
    axs[0].add_patch(Circle((l,b),rbb,fill=False,edgecolor='white'))
    axs[0].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red',linewidth=2))

axs[0].set_xlim(180,-5)


axs[1].imshow(wco,extent=lims,origin='lower',cmap = dc.cmap_grey0('viridis'),vmax=25)


for l,b,rbb,r_big in \
zip(cat['lon<deg>'],cat['lat<deg>'],\
    cat['effective_radius<deg>']\
    ,cat['accumulation_radius<deg>']):    
    axs[1].add_patch(Circle((l,b),rbb,fill=False,edgecolor='white',linewidth=2))
    axs[1].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red',linewidth=2))


axs[1].set_xlim(5,-180)


r_co = cat['effective_radius<deg>']
r_hi = cat['accumulation_radius<deg>']


a_co = cat['area_sky<deg2>']
a_hi = np.pi*(cat['accumulation_radius<deg>'])**2


bigger_co = []


# count = 0
for i in range(len(r_co)):
    if a_co.loc[i] > a_hi.loc[i]:
        bigger_co.append(i)
        

bigger_co_cat = cat.loc[bigger_co]


bcc_ms = (bigger_co_cat['mass<kg>'].to_numpy())*u.kg


(bcc_ms/(2*np.pi*(0.025*c.M_sun/(u.pc**3))))**(1/3)


bigger_co_cat['accumulation_radius<pc>'].to_numpy()*u.pc


bigger_co_cat['effective_radius<pc>'].to_numpy()*u.pc


fig,axs = plt.subplots(2,1,figsize=(25,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
axs[0].imshow(whi,extent=lims,origin='lower')


for l,b,rbb,r_big in \
zip(bigger_co_cat['lon<deg>'],bigger_co_cat['lat<deg>'],\
    bigger_co_cat['effective_radius<deg>'],\
    bigger_co_cat['accumulation_radius<deg>']):    
    axs[0].add_patch(Circle((l,b),rbb,fill=False,edgecolor='white'))
    axs[0].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red'))

axs[0].set_xlim(180,-5)


axs[1].imshow(whi,extent=lims,origin='lower')
for l,b,rbb,r_big in \
zip(bigger_co_cat['lon<deg>'],bigger_co_cat['lat<deg>'],\
    bigger_co_cat['effective_radius<deg>']\
    ,bigger_co_cat['accumulation_radius<deg>']):    
    axs[1].add_patch(Circle((l,b),rbb,fill=False,edgecolor='white'))
    axs[1].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red'))

axs[1].set_xlim(5,-180)


fig,axs = plt.subplots(2,1,figsize=(15,10))

lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]
axs[0].imshow(wco,extent=lims,origin='lower',cmap = dc.cmap_grey0('YlGn_r'),vmax=25,alpha=0.6)


for l,b,rbb,r_big in \
zip(bigger_co_cat['lon<deg>'],bigger_co_cat['lat<deg>'],\
    bigger_co_cat['effective_radius<deg>']\
    ,bigger_co_cat['accumulation_radius<deg>']):    
    axs[0].add_patch(Circle((l,b),rbb,fill=False,edgecolor='purple',linewidth=2))
    axs[0].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red',linewidth=2))

axs[0].set_xlim(180,-5)


axs[1].imshow(wco,extent=lims,origin='lower',cmap = dc.cmap_grey0('YlGn_r'),vmax=25,alpha=0.6)


for l,b,rbb,r_big in \
zip(bigger_co_cat['lon<deg>'],bigger_co_cat['lat<deg>'],\
    bigger_co_cat['effective_radius<deg>']\
    ,bigger_co_cat['accumulation_radius<deg>']):    
    axs[1].add_patch(Circle((l,b),rbb,fill=False,edgecolor='purple',linewidth=2))
    axs[1].add_patch(Circle((l,b),r_big,fill=False,edgecolor='red',linewidth=2))


axs[1].set_xlim(5,-180)


count


len(cat)


# fig,ax = plt.subplots(1,1,figsize=(15,10))

# lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]
# ax.imshow(whi,extent=lims)


# for l,v,r,vrms in zip(cocat['lon'],cocat['vel'],ang_R_A,4*vrms_3x):
#     ax.add_patch(Ellipse((l,v),2*r,vrms,fill=False,edgecolor='white',alpha=0.5))
    
    
# ax.set_ylim(-50,50)


# fig,ax = plt.subplots(1,1,figsize=(15,10))

# lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
# ax.imshow(whi,extent=lims)

# for l,b,r in zip(cocat['lon'],cocat['lat'],ang_R_A):
#     ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5))
# #     print(l,b,r)
    
    
# # ax.set_xlim(180,130)


np.all(parents==False)


good_specs


fig,ax = plt.subplots(1,1,figsize=(15,10))
lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]

ax.imshow(np.nansum(parents,axis=0),extent=lims,origin='lower')
for l,b,r in zip(good_specs['lon<deg>'],good_specs['lat<deg>'],good_specs['accumulation_radius<deg>']):
    ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5,linewidth=2))
#     print(l,b,r)


good_specs


fig,ax = plt.subplots(1,1,figsize=(15,10))
lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]

ax.imshow(np.nansum(parents,axis=1),extent=lims,origin='lower')

ax.set_ylim(-50,50)
for l,v,r,vrms in zip(good_specs['lon<deg>'],good_specs['vel<km/s>'],\
                      good_specs['accumulation_radius<deg>'],vrms3x):
    ax.add_patch(Ellipse((l,v),2*r,4*vrms,\
                         fill=False,edgecolor='white',alpha=0.5,linewidth=2))
    
# ax.set_xlim(80,50)
    

np.all(np.isnan(specs))


good_specs


cat


help(HI4PI.create_subcube)


lat_inds = cld.ind_ang('get index', self.lats, (b0,b1))
lon_inds = cld.ind_ang('get index', self.lons, (l0,l1))

#fixes index orders in case data is in decreasing order
lat_inds.sort()
lon_inds.sort()

b0, b1 = lat_inds #get each element separately
l0, l1 = lon_inds


0.00175


70/0.00175


for i,idx in enumerate(cat['_idx'].loc[:3]):
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
        
#         fits = gauss_fit(hi_spec,HI4PI.vels/1e3)
        
#         for g in fits:
#             plt.plot(HI4PI.vels/1e3,g/np.nanmax(g))
    
        
        plt.plot(DHT_CO.vels,specs[i]/np.nanmax(specs[i]))
        plt.plot(HI4PI.vels/1e3, hi_spec/(np.nanmax(hi_spec)))
        plt.title(str(idx))
        plt.show()


from scipy.signal import find_peaks


help(find_peaks)


find_peaks(hi_spec,height=18)[0]


pks = find_peaks(hi_spec,height=18)


np.nanstd(hi_spec)


from scipy.signal import peak_widths


help(peak_widths)


out = peak_widths(hi_spec,pks[0])


pks


out


plt.plot(hi_spec)
plt.hlines(18,0,300)
plt.scatter(165,76.9,color='red')
plt.scatter(183,59.2,color='red')
plt.scatter(203,59.2,color='red')

# plt.hlines(38,0,300,color='orange')
plt.hlines(77-38,0,300,color='orange')

# plt.hlines(41,0,300,color='purple')
plt.hlines(59 - 41,0,300,color='purple')

# plt.hlines(39,0,300,color='green')
plt.hlines(59-39,0,300,color='green')


xs = np.arange(0,300)


def gaussian(xs, peaks, FWHMs, means):
    """Return a Gaussian function."""
    gs = []
    for p,fwhm,mn in zip(peaks,FWHMs,means):
        sigma = fwhm / 2.354820045  # (2 * sqrt( 2 * ln(2)))
        print(p,xs,mn,sigma)
        g = p * np.exp(-(xs - mn)**2 / 2. / sigma**2)
        gs.append(g)
        
    return gs


vx = pks[1]['peak_heights']


vx


gs


gs = np.array(gs)


gs.shape


for g in gs:
    plt.plot(g)
    
plt.plot(hi_spec)
    
plt.show()


gs = gaussian(xs,pks[1]['peak_heights'],out[0],pks[0])
# plt.plot(hi_spec)


def gaussian(xs, peaks, FWHMs, means):
    """Return a Gaussian function."""
    gs = []
    for p,fwhm,mn in zip(peaks,FWHMs,means):
        sigma = fwhm / 2.354820045  # (2 * sqrt( 2 * ln(2)))
        g = p * np.exp(-(xs - mn)**2 / 2. / sigma**2)
        gs.append(g)
    return gs


def gauss_fit(spec,axis,thresh=None):
    from scipy.signal import find_peaks,peak_widths
    
    
    if thresh==None:
        thresh=np.nanstd(spec)
    
    fnd_pk_out = find_peaks(spec,height=thresh)
    pk_xs = fnd_pk_out[0] #peak positions along x
    pk_ys = fnd_pk_out[1]['peak_heights'] #peak positions along y
    
    pk_wdth_out = peak_widths(spec,pk_xs)
    wdths = pk_wdth_out[0] #fwhms
    
    decomposed = gaussian(axis,pk_ys,wdths,pk_xs)
    decomposed = np.array(decomposed)
    
    nump_comps = len(decomposed)
    
    return decomposed


for 


gauss_fit(hi_spec)


import gausspyplus


for i,idx in enumerate()


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
