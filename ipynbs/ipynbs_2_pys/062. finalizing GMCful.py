#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

from astropy.io import fits
import datacube as dc
import interpolate as ip

import matplotlib.pyplot as plt

import astropy.units as u


hi4pi = fits.open('HI4PI.fits')[0]


wh = np.nansum(hi4pi.data,axis=0)


from astropy.table import QTable


zcat = QTable(fits.open('back2zucker_cat.fits')[1].data)


w = np.nansum(dht.data,axis=0)


ipd_dht_wcent = ip.interpolate_vlb(dht_wcent)


wwcentip = np.nansum(ipd_dht_wcent,axis=0)


zcat.write('final_cat.fits',overwrite=True) #also figure out how to better save header info (aka units)


import clouds as cld


((zcat['b_bounds'][0]+30.125)/0.125).astype('int')


vels = np.arange(494)*1.3-319.8


fits.open('COGAL_all_mom.fits')[0].header


ipd_dht_wcent = ipd_dht_wcent[169:324,:,:]


zcat.remove_column('col8')


zcat['vspectrum_GMC'] = [np.zeros(155)]*len(zcat)


from astropy.modeling import models,fitting
fitter = fitting.LevMarLSQFitter()

def get_g_fit(spec):
    ginit = models.Gaussian1D(np.nanmax(spec),np.argmax(spec),50*np.nanstd(spec))
    #50std bc otherwise, aqr, cam, and lac had fits that were too narrow
#     ginit = models.Gaussian1D(1.,0.,1.)
    gfit = fitter(ginit,np.arange(len(spec)),spec)
    
#     return gfit.parameters
    return gfit#(np.arange(len(spec)))


# aquila r, camelopardis and 
# lacerta have weirdly narrow gaussian fits,
# and lam orionis has two peaks in the spectrum,
# but only the tallest one is getting fitted - so will have to do this by handd?
# 


zcat.rename_column('distance','distance<pc>')


zcat['lon<deg>'][0] = np.array((35.,11.,))


# plt.figure(figsize=(20,7))
for i in range(len(zcat)):
#     plt.imshow(wwcentip,origin='lower',extent=[180,-180,-30,30],cmap=dc.cmap_min0('viridis'),vmax=15,vmin=1)
#     plt.imshow(w,origin='lower',extent=[180,-180,-30,30],cmap=dc.cmap_min0('viridis'),alpha=0.3)

    lo,li = ((zcat['lon<deg>'][i]-180.125)/(-0.125)).astype('int')
    bo,bi = ((zcat['lat<deg>'][i]+30.125)/0.125).astype('int')
    
    sig,ar = cld.surface_density(ipd_dht_wcent[:,bo:bi,lo:li],'CO')
    
    total_area = ar*np.nansum(sig>0)*(zcat['distance<pc>'][i]*u.pc)**2
#     print(total_area)
    r = np.sqrt(total_area/np.pi)
#     plt.imshow(sig.value,origin='lower',cmap=dc.cmap_min0('viridis'))
#     plt.show()
    m = cld.mass(sig,zcat['distance<pc>'][i]*u.pc,ar)
    ra,rad = cld.accumulation_radius(m,zcat['distance<pc>'][i]*u.pc,zbs[i]*u.deg)
    zcat['radius_GMC<pc>'] = r
    zcat['radius_GMC<deg>'] = np.sqrt((u.rad**2).to(u.deg**2)*ar*np.nansum(sig>0)/np.pi)

#     zcat['radius_HI<pc>'][i] = ra.value
#     zcat['radius_HI<deg>'][i] = rad.value
    spc = np.nanmean(ipd_dht_wcent[:,bo:bi,lo:li],axis=(1,2))
#     plt.plot(vels[169:324],spc,marker='x')
    fit = get_g_fit(spc)
#     plt.plot(vels[169:324],fit(np.arange(155)))
    amp, mn, stddev = fit.parameters
#     zcat['vdispersion_GMC<km/s>'][i] = stddev*1.3
#     print(amp,vels[169:324][int(mn)],stddev*1.3*2*np.sqrt(2*np.log(2)))
#     print(dc.dispersion_1sigma(spc)) #need to change this so that the init am/mn/std are taken from the spectrum
#     plt.hlines(0.5*np.nanmax(spc),-25,25)
#     plt.title(zcat['names'][i])
#     plt.xlim(-25,25)
#     plt.show()
#     print(zcat['names'][i],get_g_fit(spc))
    
#     zcat['mass_GMC<kg>'][i] = m.value
#     zcat['vspectrum_GMC<K>'][i] = spc
#     plt.imshow(sig.value,origin='lower',cmap=dc.cmap_min0('viridis'))
#     plt.show()
    
    
#     plt.show()


zcat.write('final_cat.fits',overwrite=True)


for i in range(19):
    plt.plot(vels[169:324],zcat['vspectrum_GMC<K>'][i])
    plt.title(zcat['names'][i])
    plt.show()


gal_ds = pd.read_csv('whole_gal_dists.csv')


zls = np.nanmean(zcat['l_bounds'],axis=1)
zbs = np.nanmean(zcat['b_bounds'],axis=1)


z_idxs = np.zeros(len(zcat),dtype='int')


zls_ext = np.nanstd(zcat['l_bounds'],axis=1)
zbs_ext = np.nanstd(zcat['b_bounds'],axis=1)


gal_ds['name'] = ['']*len(gal_ds)


for i in range(len(zcat)):
    r = np.sqrt((zls[i] - gal_ds['Long.'])**2 + (zbs[i] - gal_ds['Lat.'])**2)
    
    maj_ax = np.nanmax((zls_ext[i],zbs_ext[i]))
    these = gal_ds[r<maj_ax]
    gal_ds['name'][r<maj_ax] = zcat['names'][i]
#     idx=np.argmin(r)
#     print(np.nanmin(r),'\n\n',gal_ds.iloc[idx],'\n\n',zcat[i],'\n\n\n\n')
#     print(r,print(r[np.argmin(r)]))
#     print([gal_ds.iloc[idx] for idx in range(len(gal_ds)) if r[idx] < 0.1])
#     z_idxs[i] = np.argmin(r)


zclouds_reidds = gal_ds[gal_ds['name'] !='']


wwcent = np.nansum(dht_wcent,axis=0)


plt.figure(figsize=(20,7))
plt.scatter(zclouds_reidds['Long.'],zclouds_reidds['Lat.'],\
            c=zclouds_reidds['Dist1'],s=10)
plt.xlim(180,-180)
plt.ylim(-30,30)
plt.imshow(wwcent,origin='lower',extent=[180,-180,-30,30],cmap=dc.cmap_min0('binary_r'))


x = gal_ds['Dist1']*np.sin(gal_ds['Long.']*u.deg.to(u.rad))
y = gal_ds['Dist1']*np.cos(gal_ds['Long.']*u.deg.to(u.rad))


zcat_in_gal = gal_ds.iloc[z_idxs]


xz = zcat_in_gal['Dist1']*np.sin(zcat_in_gal['Long.']*u.deg.to(u.rad))
yz = zcat_in_gal['Dist1']*np.cos(zcat_in_gal['Long.']*u.deg.to(u.rad))


zx = 1e-3*zcat['distance']*np.sin(zls*u.deg.to(u.rad))
zy = 1e-3*zcat['distance']*np.cos(zls*u.deg.to(u.rad))


zcat_in_gal


zcat


for i,idx in enumerate(zcat_in_gal.index):
    plt.figure(figsize=(10,10))
    plt.scatter(x,y,s=5,color='grey')
    plt.scatter(0,8,s=105,marker='x',color='red',linewidth=3)
    plt.scatter(xz[idx],yz[idx],s=55,label='calculated')
    plt.scatter(zx[i],zy[i],s=55,label='zucker')
    plt.title(zcat['names'][i])
    plt.ylim(-20,25)
    plt.legend()
    plt.show()


dht = fits.open('DHT_CO.fits')[0]


dht_wcent = fits.open('COGAL_all_mom.fits')[0]
old_head = dht_wcent.header
dht_wcent = dht_wcent.data.transpose().swapaxes(1,2)
# DHT_CO.cube = DHT_CO.cube[209:285,:,:] #trim to narrower velocity range
# DHT_CO.vels = DHT_CO.vels[209:285]

# smooth_whole = ip.interpolate_vlb(DHT_CO.cube)


dht_wcent.shape


plt.imshow(np.nansum(dht_wcent,axis=0),\
           cmap=dc.cmap_min0('viridis'),origin='lower',extent=[180,-180,-30,30])


zcat['lon'].iloc[2]


for i in range(len(zcat)):
    plt.imshow(np.nansum(dht_wcent,axis=0),\
           cmap=dc.cmap_min0('viridis'),origin='lower',extent=[180,-180,-30,30])
    plt.xlim(zcat['lon'].iloc[i])
    plt.ylim(zcat['lat'].iloc[i])
    plt.show()
