#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
import datacube as dc
from astropy.table import QTable

import fun as fns


hi4pi,dht = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


plt.imshow(np.nansum(dht.data,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))


his = fits.open('hi_clouds.fits')


h2s = fits.open('h2_clouds.fits')


hi_conv = (1.3*u.km/u.s *1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p *u.K).to(u.solMass/u.pc**2) #last term is data

co_conv = (1.3*u.km/u.s *1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p *u.K).to(u.solMass/u.pc**2) #last term is data


siggmc = [np.nansum(h2s[i].data,axis=0)*co_conv for i in range(17)]
sighi = [np.nansum(his[i].data,axis=0)*hi_conv for i in range(17)]


#to further optimize -- shave the top & bottom couple lats from chamaeleon?
#figure out how to set the width of each column (or at least the spacing, so that
#the spae between c2 and c3 isnt gaping like that)
order = np.ravel(np.array(([0,8,10,16,9,5],[7,2,1,15,4,14],[12,11,3,13,6,6])).transpose())


nr,nc,ax_locs = fns.get_ax_locs(17)


fig,axs = plt.subplots(nr,nc,figsize=(20,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    
    this_ax.imshow(siggmc[i].value,origin='lower',cmap=dc.cmap_min0('viridis'),
                  extent=[zcat['l_pix_GMC'][i][0]*-0.125+360,zcat['l_pix_GMC'][i][1]*-0.125+360,
                          zcat['b_pix_GMC'][i][0]*0.125-50,zcat['b_pix_GMC'][i][1]*0.125-50])
    this_ax.set(title=zcat['names'][i])
    
    
for ax in axs.flat:
    ax.set(ylabel=r'galactic latitude [$\degree$]',xlabel=r'galactic longitude [$\degree$]')
#     ax.label_outer()
    
fig.tight_layout()


#figure out how to get the contours that look cleanest
#also lacerta seems to be taller in co than 21cm -> also shave this down bc its empty at the edges


fig,axs = plt.subplots(nr,nc,figsize=(15,30))
for i,ax,whi in zip(range(17),ax_locs,whis):
    this_ax = axs[ax[0],ax[1]]
    
    
    this_ax.imshow(sighi[i].value,origin='lower',cmap=dc.cmap_min0('viridis'),\
                   extent=[zcat['l_pix_HI'][i][0]*-0.125+360,zcat['l_pix_HI'][i][1]*-0.125+360,
                          zcat['b_pix_HI'][i][0]*0.125-50,zcat['b_pix_HI'][i][1]*0.125-50])
    
    this_ax.contour(np.arange(zcat['lon_GMC<deg>'][i][0],zcat['lon_GMC<deg>'][i][1],-0.125),
                   np.arange(zcat['lat_GMC<deg>'][i][0],zcat['lat_GMC<deg>'][i][1],0.125),
                   wcos[i],levels=[2,5],colors='white')
    
    this_ax.set(title=zcat['names'][i])
    
    
for ax in axs.flat:
    ax.set(ylabel=r'galactic latitude [$\degree$]',xlabel=r'galactic longitude [$\degree$]')
#     ax.label_outer()
    
fig.tight_layout()


import pandas as pd


calrsig = pd.DataFrame(np.array((np.ravel(fns.r_arr(siggmc[0].shape[1],siggmc[0].shape[0])),\
                       np.ravel(siggmc[0].value))).transpose(),columns=['distance','surface_density']).sort_values('distance')


def binned_sig_of_r(sigmap,distance):
    df = pd.DataFrame(
            np.array((
                np.ravel(distance*0.125*u.deg.to(u.rad)*
                    fns.r_arr(sigmap.shape[1],sigmap.shape[0])
                ),
                np.ravel(sigmap.value)
            )).transpose(),
        columns=['r','sig'])
    
    drs = np.arange(0,np.nanmax(df['r']),2)
    sigs_of_r = np.zeros(len(drs))
    
    
    for i in range(len(drs)-1):
        sigs_of_r[i] = np.nanmean(df.query(str(drs[i])+'< r <='+str(drs[i+1]))['sig'])
    
    sigs_of_r[-1] = np.nanmean(df.query(str(drs[-1])+'< r'))
    
    return (drs,sigs_of_r)  


binned_hi_profiles = [binned_sig_of_r(sighi[i],zcat['distance<pc>'][i]) for i in range(17)]
binned_co_profiles = [binned_sig_of_r(siggmc[i],zcat['distance<pc>'][i]) for i in range(17)]


fig,axs = plt.subplots(6,3,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    
    this_ax.plot(np.log10(binned_co_profiles[i][0]),np.log10(binned_co_profiles[i][1]),marker='+')
    this_ax.set(title=zcat['names'][i])#,xscale='log',yscale='log')
    
    
for ax in axs.flat:
    ax.set(ylabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',xlabel=r'radial distance [pc]')
#     ax.label_outer()
    
fig.tight_layout()


fig,axs = plt.subplots(6,3,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    
    this_ax.scatter(zcat['distance<pc>'][i]*0.125*u.deg.to(u.rad)*\
                    np.ravel(fns.r_arr(siggmc[i].shape[1],siggmc[i].shape[0])),\
                    np.ravel(siggmc[i]),s=3e3/wcos[i].size)
    
    this_ax.plot(binned_co_profiles[i][0],binned_co_profiles[i][1],marker='+',color='red',linestyle='--')

    this_ax.set(title=zcat['names'][i],xscale='log',yscale='log',ylim=(1e-2,1e3))
    
    
for ax in axs.flat:
    ax.set(ylabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',xlabel=r'radial distance [pc]')
#     ax.label_outer()
    
fig.tight_layout()


#hmmm also the optimal order for showing sigs might be different than the optimal order for showing sig_of_r plots
#but .. would it be better to maintain a consistent cloud order across different types of plots? or group 
#similar features together? (leaning towards the latter)


#try binning in rings of a certain size (parsec)


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    
    this_ax.scatter(zcat['distance<pc>'][i]*0.125*u.deg.to(u.rad)*\
        np.ravel(fns.r_arr(sighi[i].shape[1],sighi[i].shape[0])),\
                    np.ravel(sighi[i]),s=2e3/whis[i].size,alpha=0.7)
    
    this_ax.plot(binned_hi_profiles[i][0],binned_hi_profiles[i][1],marker='+',linestyle='--',color='red')

    
    this_ax.set(title=zcat['names'][i])
   

for ax in axs.flat:
    ax.set(ylabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',xlabel=r'radial distance [pc]')
    
fig.tight_layout()


[(plt.plot(np.log10(binned_hi_profiles[i][0]),
              np.log10(binned_hi_profiles[i][1]),
              marker='+',label=zcat['names'][i]),

 ) for i in range(17)]
# plt.legend()


from astrodendro import Dendrogram


hi_pxbypx_spectra = [his[i].data.transpose().reshape(-1,his[i].data.shape[0]).transpose() for i in range(17)]


his[0].data.shape


caldndspc = Dendrogram.compute(hi_pxbypx_spectra[0],min_value=0.1, min_delta=0.3, min_npix=10)


caldndspc.data.shape


m = np.zeros(caldndspc.data.shape)
for ones in caldndspc.trunk[0].children[0].descendants:
    if ones.is_leaf:
        m += ones.get_mask()
    

def get_momN(cube,vpix,N):
    vels = np.arange(-100 + vpix[0] * 1.3,-100 + vpix[1] * 1.3,1.3)
    vels = vels[:len(cube)] #sometimes vels is 1 longer than the spectral axis of cube
    
    M1 = get_mom1(cube,vpix)
    
    xTx = np.array([cube[i]*(vels[i]-M1)**N for i in range(len(vels))])
    sigxTx = np.nansum(xTx,axis=0)
    
    M0 = 1.3*np.nansum(cube,axis=0)
    
    return (sigxTx/M0)**(1/N)


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    msk = h2s[i].data == 0
    mskd = np.ma.masked_array(data=h2s[i].data, mask=msk, fill_value = np.nan)
    
    im = get_mom1(mskd,zcat['v_pix_GMC'][i])
    
    vmn,vmx = zlims(im)
    mom0_mask = np.nansum(h2s[i].data*1.3,axis=0) > 2
    this_ax.imshow(im*mom0_mask,origin='lower',cmap=dc.cmap_min0('PiYG'),vmin=vmn,vmax=vmx,
                   extent=[zcat['l_pix_GMC'][i][0]*-0.125+360,zcat['l_pix_GMC'][i][1]*-0.125+360,
                          zcat['b_pix_GMC'][i][0]*0.125-50,zcat['b_pix_GMC'][i][1]*0.125-50])

    
    this_ax.set(title=zcat['names'][i])
   

for ax in axs.flat:
    ax.set(ylabel=r'galactic latitude [$\degree$]',xlabel=r'galactic longitude [$\degree$]')
    
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    msk = h2s[i].data == 0
    mskd = np.ma.masked_array(data=h2s[i].data, mask=msk, fill_value = np.nan)
    
    im = get_momN(mskd,zcat['v_pix_GMC'][i],2)
    
    vmn,vmx = zlims(im)
    mom0_mask = np.nansum(h2s[i].data*1.3,axis=0) > 2
    this_ax.imshow(im*mom0_mask,origin='lower',cmap=dc.cmap_min0('viridis'),vmin=vmn,vmax=vmx,
                  extent=[zcat['l_pix_GMC'][i][0]*-0.125+360,zcat['l_pix_GMC'][i][1]*-0.125+360,
                          zcat['b_pix_GMC'][i][0]*0.125-50,zcat['b_pix_GMC'][i][1]*0.125-50])

    
    this_ax.set(title=zcat['names'][i])
   

for ax in axs.flat:
    ax.set(ylabel=r'galactic latitude [$\degree$]',xlabel=r'galactic longitude [$\degree$]')
     
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    msk = h2s[i].data == 0
    mskd = np.ma.masked_array(data=h2s[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_GMC'][i])
    M2 = get_momN(mskd,zcat['v_pix_GMC'][i],2)
    
    
    this_ax.scatter(np.ravel(M1),np.ravel(M2),s=1e2/M1.shape[0])
    
    msk = his[i].data == 0
    mskd = np.ma.masked_array(data=his[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_HI'][i])
    M2 = get_momN(mskd,zcat['v_pix_HI'][i],2)
    
    
    this_ax.scatter(np.ravel(M1),np.ravel(M2),s=10/M1.shape[0],alpha=0.7)
    this_ax.set(title=zcat['names'][i])

for ax in axs.flat:
    ax.set(xlabel=r'$v_{lsr}$ [km/s]',ylabel=r'$\sigma_v$ [km/s]')
    
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    msk = his[i].data == 0
    mskd = np.ma.masked_array(data=his[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_HI'][i])
    M2 = get_momN(mskd,zcat['v_pix_HI'][i],2)
    
    
    this_ax.scatter(np.ravel(sighi[i]), np.ravel(M2),s=1e2/M1.shape[0])
    
    
    this_ax.set(title=zcat['names'][i],xscale='log',yscale='log')

for ax in axs.flat:
    ax.set(xlabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',ylabel=r'$\sigma_v$ [km/s]')
    
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    
    msk = h2s[i].data == 0
    mskd = np.ma.masked_array(data=h2s[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_GMC'][i])
    M2 = get_momN(mskd,zcat['v_pix_GMC'][i],2)
    
    
    this_ax.scatter(np.ravel(siggmc[i]), np.ravel(M2),s=1e2/M1.shape[0])
    
    
    this_ax.set(title=zcat['names'][i],xscale='log',yscale='log')

for ax in axs.flat:
    ax.set(xlabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',ylabel=r'$\sigma_v$ [km/s]')
    
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    msk = his[i].data == 0
    mskd = np.ma.masked_array(data=his[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_HI'][i])
    M2 = get_momN(mskd,zcat['v_pix_HI'][i],2)
    
    
    this_ax.scatter(np.ravel(sighi[i]), np.ravel(M1),s=1e2/M1.shape[0])
    
    
    this_ax.set(title=zcat['names'][i],xscale='log',yscale='log')

for ax in axs.flat:
    ax.set(xlabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',ylabel=r'$v_{lsr}$ [km/s]')
    
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    
    msk = h2s[i].data == 0
    mskd = np.ma.masked_array(data=h2s[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_GMC'][i])
    M2 = get_momN(mskd,zcat['v_pix_GMC'][i],2)
    
    
    this_ax.scatter(np.ravel(siggmc[i]), np.ravel(M1),s=1e2/M1.shape[0])
    
    
    this_ax.set(title=zcat['names'][i],xscale='log',yscale='log')

for ax in axs.flat:
    ax.set(xlabel=r'$\Sigma$ [$M_\odot\ pc^{-2}$]',ylabel=r'$v_{lsr}$ [km/s]')
    
fig.tight_layout()


fig,axs = plt.subplots(nr,nc,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    msk = his[i].data == 0
    mskd = np.ma.masked_array(data=his[i].data, mask=msk, fill_value = np.nan)
    
    M1 = get_mom1(mskd,zcat['v_pix_HI'][i])
    M2 = get_momN(mskd,zcat['v_pix_HI'][i],2)
    
    r = zcat['distance<pc>'][i]*0.125*u.deg.to(u.rad)*\
        np.ravel(fns.r_arr(sighi[i].shape[1],sighi[i].shape[0]))
    
    this_ax.scatter(r, np.ravel(M2),s=10/M1.shape[0])
    
    
    this_ax.set(title=zcat['names'][i],xscale='log',yscale='log')

for ax in axs.flat:
    ax.set(xlabel=r'radial distance [r]',ylabel=r'$\sigma_v$ [km/s]')
    
fig.tight_layout()


# fig,axs = plt.subplots(nr,nc,figsize=(13,20))
# for i,ax in zip(order[:17],ax_locs):
#     this_ax = axs[ax[0],ax[1]]
    
#     r = zcat['distance<pc>'][i]*0.125*u.deg.to(u.rad)*\
#         np.ravel(fns.r_arr(siggmc[i].shape[1],siggmc[i].shape[0]))
    
     
#     msk = h2s[i].data == 0
#     mskd = np.ma.masked_array(data=h2s[i].data, mask=msk, fill_value = np.nan)
    
#     M1 = get_mom1(mskd,zcat['v_pix_GMC'][i])
#     M2 = get_momN(mskd,zcat['v_pix_GMC'][i],2)
    
    
#     alpha = (5*((u.km/u.s * np.ravel(M2))**2)*(r)/\
#               (c.G*zcat['mass_GMC<kg>'][i]*u.kg*(r/zcat['radius_eff_GMC<pc>'][i])**2)).si

#     this_ax.scatter(r,np.ravel(alpha),s=2e3/M1.size,alpha=0.7)
        
#     this_ax.set(title=zcat['names'][i],xscale='log',yscale='log')
   

# for ax in axs.flat:
#     ax.set(xlabel=r'radial distance [pc]',ylabel=r'virial parameter')
    
# fig.tight_layout()


#velocity dispersion as a function of radius
#turbulence on smaller-than-the-cloud scales


plt.plot(np.nanmean(hi_pxbypx_spectra[0],axis=1)) #what would emerge if these were color coded by radius
plt.plot(hi_pxbypx_spectra[0],alpha=0.5,linestyle=':') #what would emerge if these were color coded by radius


#hmmm is trying to correlate the velocity spectra just the same as dendrograming the cube?
plt.imshow(caldndspc.index_map,extent=[0,28561,0,50000])


#regularize bin sizes for all the clouds


fig,axs = plt.subplots(6,3,figsize=(13,20))
for i,ax in zip(order[:17],ax_locs):
    this_ax = axs[ax[0],ax[1]]
    
    this_ax.hist(np.ravel(siggmc[i]),bins=np.arange(0,300,10)*u.solMass/u.pc**2,log=True,color='grey')
    this_ax.hist(np.ravel(sighi[i]),bins=np.arange(0,100,5)*u.solMass/u.pc**2,log=True,fill=False)
    this_ax.set(title=zcat['names'][i])
    
fig.tight_layout()


# zcat['mass_GMC<kg>'] = [(mass(dht.data[
#     zcat['v_pix_GMC'][i][0]:zcat['v_pix_GMC'][i][1],
#     zcat['b_pix_GMC'][i][0]:zcat['b_pix_GMC'][i][1],
#     zcat['l_pix_GMC'][i][0]:zcat['l_pix_GMC'][i][1]]*u.K,
#     zcat['distance<pc>'][i]*u.pc,'CO'))
# for i in range(17)]


# zcat['radius_acc_HI<pc>'],zcat['radius_acc_HI<deg>'] = accumulation_radius(zcat['mass_GMC<kg>'],zcat['distance<pc>']*u.pc)

# zcat['radius_eff_GMC<deg>'], zcat['radius_eff_GMC<pc>'] = np.sqrt(np.array(([(area(dht.data[
#     zcat['v_pix_GMC'][i][0]:zcat['v_pix_GMC'][i][1],
#     zcat['b_pix_GMC'][i][0]:zcat['b_pix_GMC'][i][1],
#     zcat['l_pix_GMC'][i][0]:zcat['l_pix_GMC'][i][1]]*u.K,
#     zcat['distance<pc>'][i]*u.pc))


# for i in range(17)])).transpose()/np.pi)


# zcat['l_pix_HI'] = np.array(((((np.nanmean(zcat['lon_GMC<deg>'],axis=1) +\
#   zcat['radius_acc_HI<deg>'].value)-360)/(-0.125)),\
# (((np.nanmean(zcat['lon_GMC<deg>'],axis=1) -\
#   zcat['radius_acc_HI<deg>'].value)-360)/(-0.125)))).transpose().astype('int')


# zcat['b_pix_HI'] = np.array(((((np.nanmean(zcat['lat_GMC<deg>'],axis=1) -\
#   zcat['radius_acc_HI<deg>'].value)+50)/(0.125)),\
# (((np.nanmean(zcat['lat_GMC<deg>'],axis=1) +\
#   zcat['radius_acc_HI<deg>'].value)+50)/(0.125)))).transpose().astype('int')


# zcat['mass_HI<kg>'] = [(mass(hi4pi.data[
#     zcat['v_pix_HI'][i][0]:zcat['v_pix_HI'][i][1],
#     zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1],
#     zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]]*u.K,
#     zcat['distance<pc>'][i]*u.pc,'21cm'))
# for i in range(17)]


# hdulist = fits.HDUList([fits.ImageHDU(hi4pi.data[
#     zcat['v_pix_HI'][i][0]:zcat['v_pix_HI'][i][1],
#     zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1],
#     zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]])
# for i in range(1,17)])

# i=0
# primary = fits.PrimaryHDU(hi4pi.data[
#     zcat['v_pix_HI'][i][0]:zcat['v_pix_HI'][i][1],
#     zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1],
#     zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]])

# all_hi4pis = [primary] + hdulist

# fits.HDUList(all_hi4pis).writeto('hi_clouds.fits',overwrite=True)


# hdulist = [fits.ImageHDU(dht.data[
#     zcat['v_pix_GMC'][i][0]:zcat['v_pix_GMC'][i][1],
#     zcat['b_pix_GMC'][i][0]:zcat['b_pix_GMC'][i][1],
#     zcat['l_pix_GMC'][i][0]:zcat['l_pix_GMC'][i][1]])
# for i in range(1,17)]

# i=0
# primary = fits.PrimaryHDU(dht.data[
#     zcat['v_pix_GMC'][i][0]:zcat['v_pix_GMC'][i][1],
#     zcat['b_pix_GMC'][i][0]:zcat['b_pix_GMC'][i][1],
#     zcat['l_pix_GMC'][i][0]:zcat['l_pix_GMC'][i][1]])

# all_dhts = [primary] + hdulist

# fits.HDUList(all_dhts).writeto('h2_clouds.fits',overwrite=True)

# zcat.write('final_cat.fits',overwrite=True)


#hmmmm calculating masses according to eq 5 in the carina paper gives masses diff
#from just using an 'x' factor 

# def carmass(cube,d):
#     nhi = 1.82e18*((u.cm)**-2)/(u.K*(u.km/u.s)) * np.nansum((-155*u.K)*(np.log(1 - cube/(155*u.K))),axis=0)*1.3*u.km/u.s
#     sig = nhi*1.36*c.m_p
    
#     npx = np.pi*(0.5*np.sqrt(sig.size))**2
#     dA = npx*(d*0.125*u.deg.to(u.rad))**2
#     m = np.nansum(sig*dA)
    
    
#     return m.to(u.kg).value
# carms = np.array([carmass(u.K*his[i].data,zcat['distance<pc>'][i]*u.pc) for i in range(17)])


# print(100*(carms - zcat['mass_HI<kg>'].value)/zcat['mass_HI<kg>'].value)

#the differences range from 0-14% ... this might be a problem, or might not matter, depending on other errors
#come back if other errors are smaller than this / this would be significant comparatively
