#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datacube as dc
from astropy.modeling import models,fitting
import astropy.units as u
from specutils import Spectrum1D
from specutils.fitting import find_lines_threshold
from scipy.signal import find_peaks,peak_widths


from astrodendro import Dendrogram


dnd = Dendrogram.load_from('COGAL_dendrogram.fits')


oph = dnd[260].get_mask() * dht


c


plt.imshow(np.nansum(oph,axis=0),origin='lower',dc=cmap.)


dht = fits.open('DHT_CO.fits')[0].data


g_init = models.Gaussian1D(amplitude=1.,mean=0.,stddev=1.)
fit_g = fitting.LMLSQFitter()


cat = pd.read_csv('physcl_prprts_catalog.csv')
mek = pd.read_csv('Mekuria_catalog.csv',skiprows=2)


hi4pi = fits.open('HI4PI.fits')[0]


dv,db,dl = hi4pi.header['CDELT1'],hi4pi.header['CDELT2'],hi4pi.header['CDELT3']


v0,b0,l0 = hi4pi.header['CRVAL1'],hi4pi.header['CRVAL2'],hi4pi.header['CRVAL3']
l0 = l0*-1


nv, nb, nl = hi4pi.header['NAXIS3'],hi4pi.header['NAXIS2'],hi4pi.header['NAXIS1']


ls_px = np.array([int((l-l0)/dl) for l in cat['lon<deg>']])
bs_px = np.array([int((b-b0)/db) for b in cat['lat<deg>']])


ras_px = np.array([int(ra/db) for ra in cat['accumulation_radius<deg>']])


hi_sky = hi4pi.data.astype('float32')


hi_vels = dv*np.arange(nv)+v0
hi_lats = db*np.arange(nb)+b0
hi_lons = dl*np.arange(nl)+l0


def decompose(spc):
    #find peaks,calculate their properties
    found_peaks,info = find_peaks(spc,height=np.nanstd(spc))
#     print(found_peaks,info)
    wdths, hghts, lefts, rights = peak_widths(spc,found_peaks)     
    fitss = []
    #for each peak, define a gaussian using its fwhm 
    for h,w,x0 in zip(info['peak_heights'],wdths,found_peaks):
        sigma = w / 2.354820045  # (2 * sqrt( 2 * ln(2)))
#         print(h,w,x0)
         
        gaus = models.Gaussian1D(h,x0,sigma)
        lrnz = models.Lorentz1D(h,x0,w) #probably the wrong fwhm
        voit = models.Voigt1D(x0,h,w,sigma)
        fitss.append(gaus) #this works best for spectra with 1 component
        
    
    return fitss


for i in ([226]):
    l=ls_px[i]
    b=bs_px[i]
    r=ras_px[i]
    
    sub = hi_sky[:,b-r:b+r,l-r:l+r]
    spc = dc.spectrum1D(sub)
    
    fns = decompose(spc)
    plt.plot(spc)
    for f in fns:
        plt.plot(np.arange(310),f(np.arange(310)))
        
    plt.title(cat._idx[i])
        
    plt.show()


# n_comps = np.zeros(len(cat))
for i in range(186):
    l=ls_px[i]
    b=bs_px[i]
    r=ras_px[i]
    
    sub = hi_sky[:,b-r:b+r,l-r:l+r]
    spc = dc.spectrum1D(sub)
    
    fns = decompose(spc)
    plt.plot(spc)
    for f in fns:
        plt.plot(np.arange(310),f(np.arange(310)))
        
    plt.title(cat._idx[i])
        
    plt.show()
#     n_comps[i] = len(fns)


#for taurus go back and make sure that when the ra goes > 180, tell the code to wrap around


comps_df = pd.DataFrame([(cat._idx.iloc[i],int(n_comps[i])) for i in cat.index],columns=['_idx','n_comps_parent'])


comps_df.to_csv('hi_spectrum_components.csv',index=False)


for i in cat.index:
    print(n_comps[i],cat._idx.iloc[i])


for i,idx in enumerate(cat._idx):
    print(i,idx)


len(cat[n_comps==1])


n_comps[129] #hercules


n_comps[52] #california


n_comps[75] #taurus


one_gaus_comp = cat[n_comps==1]


one_gaus_comp


ls_px_1 = np.array([int((l-l0)/dl) for l in one_gaus_comp['lon<deg>']])
bs_px_1 = np.array([int((b-b0)/db) for b in one_gaus_comp['lat<deg>']])


ras_px_1 = np.array([int(ra/db) for ra in one_gaus_comp['accumulation_radius<deg>']])


sigs_his = np.zeros(len(one_gaus_comp))
for i in range(len(one_gaus_comp)):
    l=ls_px_1[i]
    b=bs_px_1[i]
    r=ras_px_1[i]
    
    sub = hi_sky[:,b-r:b+r,l-r:l+r]
    spc = dc.spectrum1D(sub)
    
    sig_hi = dc.dispersion_1sigma(spc)
    sigs_his[i] = sig_hi
    
    fns = decompose(spc)
    fit = fns[0]
    
    plt.plot(np.arange(310),spc)
    plt.plot(np.arange(310),fit(np.arange(310)))
    plt.show()
    
    
two_gaus_comp = cat[n_comps==2]


two_gaus_comp


ls_px_2 = np.array([int((l-l0)/dl) for l in two_gaus_comp['lon<deg>']])
bs_px_2 = np.array([int((b-b0)/db) for b in two_gaus_comp['lat<deg>']])
ras_px_2 = np.array([int(ra/db) for ra in two_gaus_comp['accumulation_radius<deg>']])


#before making the cuts, put the HI spectra in context of the CO spectra (if there are multiple
#peaks in the CO also, then it's okay for the HI to also have multiple peaks)


sigs_his_2 = np.zeros(len(two_gaus_comp))

for i in range(len(two_gaus_comp)):
    l=ls_px_2[i]
    b=bs_px_2[i]
    r=ras_px_2[i]
    
    sub = hi_sky[:,b-r:b+r,l-r:l+r]
    spc = dc.spectrum1D(sub)
    
    sig_hi = dc.dispersion_1sigma(spc)
    sigs_his_2[i] = sig_hi
    
    fns = decompose(spc)
    plt.plot(np.arange(310),spc)
    for f in fns:    
        plt.plot(np.arange(310),f(np.arange(310)))
    plt.show()
    
    
one_gaus_comp['velocity_dispersion_GMC<km/s>'].to_numpy()


sigs_his*1.28


two_gaus_comp.columns


two_gaus_comp['velocity_dispersion<km/s>'].to_numpy()


plt.scatter(one_gaus_comp['velocity_dispersion_GMC<km/s>'].to_numpy(),sigs_his*1.28)
plt.plot(one_gaus_comp['velocity_dispersion_GMC<km/s>'].to_numpy()\
         ,3*one_gaus_comp['velocity_dispersion_GMC<km/s>'].to_numpy())


dht = fits.open('DHT_CO.fits')[0]


dnd = Dendrogram.load_from('COGAL_dendrogram.fits')


one_gaus = [dnd[i] for i in range(len(dnd)) if i in one_gaus_comp['_idx'].to_numpy()]


dv_co, db_co, dl_co = dht.header['CDELT1'],dht.header['CDELT2'],dht.header['CDELT3']


v0_co, b0_co, l0_co = dht.header['CRVAL1'],dht.header['CRVAL2'],dht.header['CRVAL3']


nv_co, nb_co, nl_co = dht.header['NAXIS3'],dht.header['NAXIS2'],dht.header['NAXIS1']


ls_px_1_co = np.array([int((l-l0_co)/dl_co) for l in one_gaus_comp['lon<deg>']])
bs_px_1_co = np.array([int((b-b0_co)/db_co) for b in one_gaus_comp['lat<deg>']])
ras_px_1_co = np.array([int(ra/db_co) for ra in one_gaus_comp['accumulation_radius<deg>']])


co_vels = dv_co*np.arange(nv_co)+v0_co
co_lats = db_co*np.arange(nb_co)+b0_co
co_lons = dl_co*np.arange(nl_co)+l0_co


for i,strut in enumerate(one_gaus):
    m = strut.get_mask()
    w = np.nansum(m*dht.data,axis=0)
    lims_lb = [co_lons[0],co_lons[-1],co_lats[0],co_lats[-1]]
    
    l=ls_px_1_co[i]
    b=bs_px_1_co[i]
    r=ras_px_1_co[i]
    
    
    lmn = max(0,l-r)
    lmx = min(len(co_lons)-1,l+r)
    
#     bmn = max(0,bs_px_1_co[i]-ras_px_1_co[i])
#     bmx = min(len(co_lats)+1,bs_px_1_co[i]-ras_px_1_co[i])
    
#     print(ls_px_1_co[i], ras_px_1_co[i])
#     print(one_gaus_comp.iloc[i])
    
    
    plt.imshow(w,origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims_lb)
# #     plt.title(str(strut.idx))
    plt.xlim(co_lons[lmn],co_lons[lmx])
    plt.ylim(co_lats[b-r],co_lats[b+r])
    plt.show()


#are these salvageable? look into maybe calculating kinematic distances instead of parallax based
one_gaus_comp[one_gaus_comp['effective_radius<pc>']>one_gaus_comp['accumulation_radius<pc>']]


dv


one_gaus_comp['velocity_dispersion_parent<km/s>'] = sigs_his*dv


one_gaus_comp['y_cen_hi4pi<px>'] = bs_px_1


one_gaus_comp['x_cen_hi4pi<px>'] = ls_px_1


one_gaus_comp['x_cen_dht<px>'] = ls_px_1_co


one_gaus_comp['y_cen_dht<px>'] = bs_px_1_co


one_gaus_comp['v_cen_dht<px>'] = vs_px_1_co


one_gaus_comp['v_cen_hi4pi<px>'] = vs_px_1


one_gaus_comp.columns = ['_idx', 'lon<deg>', 'lat<deg>', 'vel<km/s>', 'dist<kpc>',
       'area_sky_GMC<deg2>', 'mass_GMC<kg>', 'area_GMC<pc2>', 'radius_GMC<pc>',
       'radius_GMC<deg>', 'acc_radius_parent<pc>',
       'acc_radius_parent<deg>', 'velocity_dispersion_GMC<km/s>',
       'x_cen_dht<px>', 'y_cen_dht<px>', 'v_cen_dht<px>', 'x_cen_hi4pi<px>',
       'y_cen_hi4pi<px>', 'v_cen_hi4pi<px>','velocity_dispersion_parent<km/s>']


one_gaus_comp


vs_px_1 = np.array([int((v-v0)/dv) for v in one_gaus_comp['vel<km/s>']])


vs_px_1_co = np.array([int((v-v0_co)/dv_co) for v in one_gaus_comp['vel<km/s>']])


hi_lons[ls_px_1]


co_lons[ls_px_1_co]


nice_rs_one_gaus = one_gaus_comp[one_gaus_comp['radius_GMC<pc>']<one_gaus_comp['acc_radius_parent<pc>']]


nice_rs_one_gaus.to_csv('lfs_w_simple_spectra_nice_r_acc.csv',index=False)


sigs_his


from sklearn.mixture import GaussianMixture


#this might be a more robust way of doing gaussian decomposition ? ... idk
help(GaussianMixture)


for i in range(186):
    l=ls_px[i]
    b=bs_px[i]
    r=ras_px[i]
    
    sub = hi_sky[:,b-r:b+r,l-r:l+r]
    spc = dc.spectrum1D(sub)
