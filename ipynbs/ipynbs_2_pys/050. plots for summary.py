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


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float16') #changed it from 'float32' in case this makes things quicker
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]


HI4PI.cube.size


HI4PI.cube.itemsize


1608156000*2


3216312000/1e9


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
heavy_cube = trimmed_HI4PI_cube.astype('float32') #changed it from 'float32' in case this makes things quicker


DHT_CO = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')
DHT_CO.cube = DHT_CO.cube[209:285,:,:] #trim to narrower velocity range
DHT_CO.vels = DHT_CO.vels[209:285]
#prob also need to take out the center? not sure if that matters yet tho


cat = pd.read_csv('physcl_prprts_catalog.csv')
dnd = Dendrogram.load_from('COGAL_dendrogram.fits')


dnd2 = Dendrogram.load_from('COGAL_dendrogram_sig20sig.fits')


lvs1 = [dnd[i] for i in range(len(dnd)) if dnd[i].is_leaf]


lvs2 = [dnd2[i] for i in range(len(dnd2)) if dnd2[i].is_leaf]


big1 = np.zeros_like(dnd.data,dtype=bool)
for l in lvs1:
    big1 = big1 | l.get_mask()


plt.imshow(np.nansum(big1,axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower')


big2 = np.zeros_like(dnd2.data,dtype=bool)
for l2 in lvs2:
    big2 = big2 | l2.get_mask()


lims_lb = DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]


from matplotlib.patches import Circle


fig,axs = plt.subplots(2,1,figsize=(20,10))
axs[0].imshow(np.nansum(big2,axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb)
axs[1].imshow(np.nansum(big2,axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb)
axs[0].set_xlim(180,-5)
axs[1].set_xlim(5,-180)


for l,b,rc,rh in zip(bad_rs2l,bad_rs2b,bad_rs2c,bad_rs2h):
    axs[0].add_patch(Circle((l,b),rc,color='orange',fill=False))
    axs[1].add_patch(Circle((l,b),rc,color='orange',fill=False))
    
    axs[0].add_patch(Circle((l,b),rh,color='red',fill=False))
    axs[1].add_patch(Circle((l,b),rh,color='red',fill=False))


cat2 = pd.read_csv('physcl_prprts_catalog_sig20sig.csv')


bad_rs2 = cat2[cat2['effective_radius<pc>']>cat2['accumulation_radius<pc>']]


bad_rs2c = np.sqrt(bad_rs2['area_sky<deg2>'].to_numpy()/np.pi)


bad_rs2h = bad_rs2['accumulation_radius<deg>'].to_numpy()


bads = bad_rs2['_idx'].to_numpy()


fig,axs = plt.subplots(2,1,figsize=(20,10))
bad_msks = [dnd2[i].get_mask() for i in bads]
for bm in bad_msks:
    axs[0].imshow(np.nansum(bm,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims_lb)
    axs[1].imshow(np.nansum(bm,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims_lb)

    axs[0].set_xlim(180,-5)
    axs[1].set_xlim(5,-180)
plt.show()


mekcat = pd.read_csv('Mekuria_catalog.csv',skiprows=5)


18/0.1


8/0.1


mekcat


aquila: 217(102),248(121),249(122) -- d: x2
california: 173(75),133(52), 136(53) -- d:x0.5
perseus: 214(100) -- d:x1
taurus: 264(134) -- d:x1.5
ora: 246(120) -- d:x2
orb: 250(123) -- d:x2
crossbones:304(171) -- d:1-1.5x
monr2: 301(168) -- d:1.5x
scissors: 202(93)
    
chamaeleon: 218(103),216(101) -- d:xBIG
corona: 275(145) -- d:xBIG
        
hercules: 257(129) -- d:1x
polaris: 157(66) -- d:1-1.5x
        
ophiuchus: 260(131) -- d:xBIG


bad_rs2


cat2.loc[[66,131]]


for n,l,b in zip(mekcat['name'].to_numpy(),mekcat['l_c__deg'].to_numpy(),mekcat['b_c__deg'].to_numpy()):
    if l>180:
        l-=360
    plt.figure(figsize=(20,10))
    plt.imshow(np.nansum(big2,axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb)
    plt.xlim(l+10,l-7)
    plt.ylim(b-5,b+7)
    plt.title(n)
    plt.show()
   

len(mekcat)


dnd[248]


dnd2[249]


dnd2[217]


for l,b in zip(mekcat['l_c__deg'].to_numpy(),mekcat['b_c__deg'].to_numpy()):
    print(l,b)


mekcat


bad_msks = [dnd2[i].get_mask() for i in bads]
plt.figure(figsize=(20,10))

for bm in bad_msks:
    plt.imshow(np.nansum(bm,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims_lb)
    
plt.xlim(30,-40)
plt.show()


bad_rs2c


bad_rs2h


bad_rs2b = bad_rs2['lat<deg>'].to_numpy()
bad_rs2l = bad_rs2['lon<deg>'].to_numpy()


sigs = np.zeros(len(cat2))
specs = np.zeros((len(cat2),len(DHT_CO.vels)))

for i,idx in enumerate(cat2['_idx']):
    cube = dnd2[idx].get_mask() * dnd2.data
    spc = dc.spectrum1D(cube)
    sig = dc.dispersion_1sigma(spc)
    
    if len(sig)==0 or (np.argmax(spc) < 5*sig) or (np.argmax(spc) > len(spc) - 5*sig):
        sigs[i] = np.nan
        specs[i] = np.nan
    else:
        sigs[i] = sig[0]*1.3
        specs[i] = spc 
    
cat2['sigma<km/s>'] = sigs


good_specs = cat2.dropna(subset='sigma<km/s>')


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


parents.shape


w_par = np.nansum(parents,axis=0)


w_par_1 = np.nansum(parents,axis=1)


lims_lb_hi = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]


whi = np.nansum(HI4PI.cube,axis=0)


whi_1 = np.nansum(HI4PI.cube,axis=1)


lims_lv_hi = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]


good_specs


HI4PI.dv.to(u.km/u.s)


from matplotlib.patches import Ellipse


fig,axs = plt.subplots(figsize=(10,10))

axs.imshow(np.abs(whi_1)*w_par_1.astype('bool'),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lv_hi)
# axs[1].imshow(w_par_1.astype('bool'),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lv_hi)

# # axs[0].set_xlim(180,-5)
# axs[1].set_xlim(5,-180)
axs.set_ylim(-100,100)
# axs[1].set_ylim(-100,100)


for l,v,rc,rv in zip(good_specs['lon<deg>'],good_specs['vel<km/s>'],\
                     np.sqrt(good_specs['area_sky<deg2>']/np.pi),\
                     (good_specs['sigma<km/s>'])/HI4PI.dv.to(u.km/u.s).value):
    axs.add_patch(Ellipse((l,v),2*rc,4*rv,color='orange',fill=False))
#     axs[1].add_patch(Circle((l,b),rc,color='orange',fill=False))
    
#     axs[0].add_patch(Circle((l,b),rh,color='red',fill=False))
#     axs[1].add_patch(Circle((l,b),rh,color='red',fill=False))


# plt.hist(bad_rs2['dist<kpc>'])
plt.hist(cat2['dist<kpc>'],bins=50)
plt.hist(bad_rs2['dist<kpc>'],bins=50)
plt.yscale('log')


np.nanmean(cat2['dist<kpc>'])


fig,axs = plt.subplots(2,1,figsize=(20,10))

axs[0].imshow(w_par_1.astype('bool'),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb_hi)
axs[1].imshow(w_par_1.astype('bool'),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb_hi)

axs[0].set_xlim(180,-5)
axs[1].set_xlim(5,-180)

# for l,b,rc,rh in zip(cat2['lon<deg>'],cat2['lat<deg>'],np.sqrt(cat2['area_sky<deg2>']/np.pi),cat2['accumulation_radius<deg>']):
#     axs[0].add_patch(Circle((l,b),rc,color='orange',fill=False))
#     axs[1].add_patch(Circle((l,b),rc,color='orange',fill=False))
    
#     axs[0].add_patch(Circle((l,b),rh,color='red',fill=False))
#     axs[1].add_patch(Circle((l,b),rh,color='red',fill=False))


fig,axs = plt.subplots(2,1,figsize=(20,10))

axs[0].imshow(np.abs(whi)*(w_par.astype('bool')),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb_hi)
axs[1].imshow(np.abs(whi)*(w_par.astype('bool')),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims_lb_hi)

axs[0].set_xlim(180,-5)
axs[1].set_xlim(5,-180)

for l,b,rc,rh in zip(cat2['lon<deg>'],cat2['lat<deg>'],np.sqrt(cat2['area_sky<deg2>']/np.pi),cat2['accumulation_radius<deg>']):
    axs[0].add_patch(Circle((l,b),rc,color='orange',fill=False))
    axs[1].add_patch(Circle((l,b),rc,color='orange',fill=False))
    
#     axs[0].add_patch(Circle((l,b),rh,color='red',fill=False))
#     axs[1].add_patch(Circle((l,b),rh,color='red',fill=False))


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


for i,idx in enumerate(cat2['_idx']):
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
        
        sum_fits = np.nansum(fits,axis=0)
        plt.figure(figsize=(7,7))
#         for g in fits:
#             plt.plot(HI4PI.vels/1e3,g/np.nanmax(g),linestyle=':')
    
        
        plt.plot(DHT_CO.vels,specs[i]/np.nanmax(specs[i]),color='red',alpha=0.6,label='CO')
        plt.plot(HI4PI.vels/1e3, hi_spec/(np.nanmax(hi_spec)),color='black',alpha=0.5,label='21cm')
        plt.plot(HI4PI.vels/1e3, sum_fits/(np.nanmax(sum_fits)),color='grey',linestyle=':',label='21cm gaussian fit')
        plt.title('cloud id:'+str(idx))
        plt.xlim(-100,100)
        plt.legend()
        plt.show()


for i,idx in enumerate(cat2['_idx']):
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
        
        sum_fits = np.nansum(fits,axis=0)
        plt.figure(figsize=(7,7))
        for g in fits:
            plt.plot(HI4PI.vels/1e3,g/np.nanmax(g),linestyle=':')
    
        
        plt.plot(DHT_CO.vels,specs[i]/np.nanmax(specs[i]),color='red',alpha=0.6,label='CO')
#         plt.plot(HI4PI.vels/1e3, hi_spec/(np.nanmax(hi_spec)),color='black',alpha=0.5,label='21cm')
#         plt.plot(HI4PI.vels/1e3, sum_fits/(np.nanmax(sum_fits)),color='grey',linestyle=':',label='21cm gaussian fit')
        plt.title('cloud id:'+str(idx))
        plt.xlim(-100,100)
        plt.legend()
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


for i,idx in enumerate(cat2['_idx']):
    
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
        plt.title('cloud id:'+str(idx))
#         plt.yticks(np.arange(HI4PI.lons[l0],HI4PI.lons[l1]))
        plt.show()


help(plt.contour)


dir(ax)


for i,idx in enumerate(cat2['_idx']):
    
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
        
        ax1.set_title('cloud id:'+str(idx))
        
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


from matplotlib import colormaps


prd = colormaps.get_cmap('Paired').colors
vrds = colormaps.get_cmap('viridis').colors


len(vrds)


for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()

    w = np.nansum(m,axis=0)
    plt.imshow(w,cmap=dc.cmap_grey0('viridis'),origin='lower')
    
    plt.show()


x1,y1 = np.meshgrid(HI4PI.lons[l0:l1],1e-3*HI4PI.vels)


ws = []
for t,tr in enumerate(subhid.trunk):
    m = tr.get_mask()
#     print(np.nanmax(m),np.nanmin(m))
#     print(m.shape)
    w = np.nansum(m,axis=0)
    this = w.astype('bool')* (t+1)
    ws.append(this)
    
for t,w in enumerate(ws):
    plt.contour(x,y,w,colors=[vrds[t*10]],levels=t+1,linewidths=0.5)
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
    plt.contour(x1,y1,w,colors=[vrds[t*10]],levels=t+1,linewidths=0.5)
    
plt.ylim(-5,120)
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
