#!/usr/bin/env python
# coding: utf-8


from astrodendro import Dendrogram,periodic_neighbours
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy_healpix import HEALPix as aphp
from astropy.modeling import models

import healpy as hp
from skimage.segmentation import find_boundaries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import bootleg as btlg
import datacube as dc


hi4pi,dht = fits.open('both_same_grid_720d.fits')


dust_lbd = fits.open('mean_and_std_lbd.fits')
mean_dust, std_dust = dust_lbd[1],dust_lbd[2]
long_dust = np.concatenate((mean_dust.data[:,:,786:],mean_dust.data,mean_dust.data[:,:,:786]),axis=2)


hi_conv = 1.828e18#*((u.cm)**-2)/(u.K*(u.km/u.s)) *1.36*c.m_p)#.to(u.solMass/u.pc**2) 
co_conv = 1.97e20 #* ((u.cm)**-2)/(u.K*(u.km/u.s))\
gas2dust = 1.87e21#*(u.cm**-2)/u.mag #from zhu et al: https://arxiv.org/pdf/1706.07109.pdf


hclons = np.arange(5760)*-0.125 + 360
hclats = np.arange(800)*0.125 - 50
hcvels = np.arange(155)*1.3 - 100

dlons = np.arange(3144)*-0.23+360
dlats = np.arange(786)*0.23 - 90
ddist = np.arange(1182) + 68


hi_lb = np.nansum(hi4pi.data,axis=0)
hi_bv = np.nansum(hi4pi.data,axis=2)
hi_lv = np.nansum(hi4pi.data,axis=1)

co_lb = np.nansum(dht.data,axis=0)
co_bv = np.nansum(dht.data,axis=2)
co_lv = np.nansum(dht.data,axis=1)


d_lb = np.nansum(long_dust,axis=0)
d_ld = np.nansum(long_dust,axis=1)
d_bd = np.nansum(long_dust,axis=2)


0.2*gas2dust


fig,axs = plt.subplots(2,1,figsize=(15,15),layout='constrained')

axs[0].imshow(hi_lb*hi_conv,origin='lower',
                   extent=[360,-360,-50,50],norm='log')

axs[0].contour(hclons,hclats,hi_lb*hi_conv,levels=[4e20,1e21],cmap='viridis',
               linewidths=0.7)

axs[0].contour(hclons,hclats,co_lb*co_conv,levels=[1e21],colors='white',
               linewidths=0.7)

axs[0].contour(dlons,dlats,d_lb*gas2dust,levels=[4e20],colors='red',
               linewidths=0.7)


axs[1].imshow(hi_lb*hi_conv,origin='lower',
                   extent=[360,-360,-50,50],norm='log')

axs[1].contour(hclons,hclats,hi_lb*hi_conv,levels=[4e20,1e21],cmap='viridis',
               linewidths=0.7)

axs[1].contour(hclons,hclats,co_lb*co_conv,levels=[1e21],colors='white',
               linewidths=0.7)

axs[1].contour(dlons,dlats,d_lb*gas2dust,levels=[4e20],colors='red',
               linewidths=0.7)
                      
# model = axs[1].imshow(np.log10(hi_conv*nhi_lb),
#            extent=[360,-360,-50,50],origin='lower',vmin=16,vmax=20)

# final = axs[2].imshow(np.log2(np.nanmean(hi4pi.data,axis=0)/nhi_lb_flat),
#                       origin='lower',extent=[360,-360,-50,50],vmin=-5,vmax=5)

# axs[2].contour(lons,lats,np.log10(co_lb*co_conv),levels=[20],colors='white',
#                linewidths=0.7)

axs[0].set_xlim(190,0)
axs[1].set_xlim(0,-190)
for ax in axs.flat:
    ax.set(xlabel=r'galactic longitude [deg]',
          ylabel=r'galactic latitude [deg]',
          ylim=(-50,50))
    ax.label_outer()


# fig.colorbar(og,    ax=axs[0], orientation='vertical',label=r'mean $N_{HI}$ [$cm^{-2}$]')
# fig.colorbar(model, ax=axs[1], orientation='vertical',label=r'mean $N_{HI}$ [$cm^{-2}$]')
# fig.colorbar(final, ax=axs[2], orientation='vertical',label=r'$log_2$ $[\ N_{HI} / \overline{N_{HI}}$ ]')


fig,axs = plt.subplots(2,1,figsize=(15,15),layout='constrained')

axs[0].imshow(hi_lv*hi_conv,origin='lower',
                   extent=[360,-360,-100,100],norm='log')

axs[0].contour(hclons,hcvels,hi_lv*hi_conv,levels=[4e20,1e21,1e22],cmap='viridis',
               linewidths=0.7)

axs[0].contour(hclons,hcvels,co_lv*co_conv,levels=[1e21],colors='white',
               linewidths=0.7)

# axs[0].contour(dlons,dlats,d_lb*gas2dust,levels=[4e20],colors='red',
#                linewidths=0.7)


axs[1].imshow(hi_lv*hi_conv,origin='lower',
                   extent=[360,-360,-100,100],norm='log')

axs[1].contour(hclons,hcvels,hi_lv*hi_conv,levels=[4e20,1e21,1e22],cmap='viridis',
               linewidths=0.7)

axs[1].contour(hclons,hcvels,co_lv*co_conv,levels=[1e21],colors='white',
               linewidths=0.7)

# axs[1].contour(dlons,dlats,d_lb*gas2dust,levels=[4e20],colors='red',
#                linewidths=0.7)
                      
# model = axs[1].imshow(np.log10(hi_conv*nhi_lb),
#            extent=[360,-360,-50,50],origin='lower',vmin=16,vmax=20)

# final = axs[2].imshow(np.log2(np.nanmean(hi4pi.data,axis=0)/nhi_lb_flat),
#                       origin='lower',extent=[360,-360,-50,50],vmin=-5,vmax=5)

# axs[2].contour(lons,lats,np.log10(co_lb*co_conv),levels=[20],colors='white',
#                linewidths=0.7)

axs[0].set_xlim(190,0)
axs[1].set_xlim(0,-190)
for ax in axs.flat:
    ax.set(xlabel=r'galactic longitude [deg]',
          ylabel=r'galactic latitude [deg]')
#           ylim=(-50,50))
    ax.label_outer()


# fig.colorbar(og,    ax=axs[0], orientation='vertical',label=r'mean $N_{HI}$ [$cm^{-2}$]')
# fig.colorbar(model, ax=axs[1], orientation='vertical',label=r'mean $N_{HI}$ [$cm^{-2}$]')
# fig.colorbar(final, ax=axs[2], orientation='vertical',label=r'$log_2$ $[\ N_{HI} / \overline{N_{HI}}$ ]')


nhi_bv = np.outer(nhi_of_v,nhi_of_b)
nhi_lb = np.outer(nhi_of_b,nhi_of_l)
nhi_lv = np.outer(nhi_of_v,nhi_of_l)

nhi_lb_flat = np.broadcast_to(nhi_of_b,(5760,800)).transpose()


nhi_of_b = np.median(hi4pi.data,axis=(0,2))
nhi_of_v = np.median(hi4pi.data,axis=(1,2))
nhi_of_l = np.median(hi4pi.data,axis=(0,1))


((((10*u.solMass/(u.pc**2))/(c.m_p*1.36)).to(u.cm**-2)))/((10*u.pc).to(u.cm))


2*u.deg.to(u.rad)*400*u.pc


def N_H_of_l_notebook(l):
    l = (l- 60)*u.deg.to(u.rad) 
    
    nh = 50
    R = 18*u.kpc.to(u.cm)
    Ro = 8.5*u.kpc.to(u.cm)
    del_l = 0.125*u.deg.to(u.rad)
    integ_b = 50-(-50)
    
    L_l = Ro*np.cos(l) + np.sqrt(R**2 - (Ro*np.sin(l))**2)
    
    N_H_of_l = nh*del_l*integ_b*L_l
    
    return N_H_of_l


fig,axs = plt.subplots(2,1,figsize=(10,10))
axs[0].plot(hclons,np.log10(np.nansum(hi4pi.data,axis=(0,1))*hi_conv))


#plotting N_H(l) that i derived to see how it fits/doesn't
axs[0].plot(hclons,np.log10(N_H_of_l_notebook(hclons)))

axs[1].scatter(hclons,np.log10(np.nansum(hi4pi.data,axis=(0,1))*hi_conv),s=5)
axs[1].scatter(-hclons,np.log10(np.nansum(hi4pi.data,axis=(0,1))*hi_conv),s=5)

axs[1].plot(hclons,np.log10(N_H_of_l_notebook(hclons)))
axs[1].plot(-hclons,np.log10(N_H_of_l_notebook(hclons)))


for ax in axs.flat:
    ax.set(xlabel=r'galactic longitude [$deg$]', ylabel=r'log (integrated $N_{HI}$ [$cm^{-2}$])')
    
axs[1].set(xscale='log')


hi_lb_md = np.nanmedian(hi_lb,axis=0)
hi_bv_md = np.nanmedian(hi_bv,axis=0)
hi_lv_md = np.nanmedian(hi4pi.data,axis=1)

# co_lb_md = np.nanmedian(dht.data,axis=0)
# co_bv_md = np.nanmedian(dht.data,axis=2)
# co_lv_md = np.nanmedian(dht.data,axis=1)


# d_lb_md = np.nanmedian(mean_dust.data,axis=0)
# d_ld_md = np.nanmedian(mean_dust.data,axis=1)
# d_bd_md = np.nanmedian(mean_dust.data,axis=2)


fig,axs = plt.subplots(2,1,figsize=(10,10))
axs[0].plot(hclats,np.log10(np.nansum(hi4pi.data,axis=(0,2))*hi_conv))
axs[0].plot(hclats,np.log10(models.Lorentz1D(4.5e25,0.,5)(hclats)))
            
axs[1].scatter(hclats,np.log10(np.nansum(hi4pi.data,axis=(0,2))*hi_conv),s=10)
axs[1].scatter(-hclats,np.log10(np.nansum(hi4pi.data,axis=(0,2))*hi_conv),s=10)
axs[1].plot(hclats,np.log10(models.Lorentz1D(4.5e25,0.,5)(hclats)))

for ax in axs.flat:
    ax.set(xlabel=r'galactic latitude [$deg$]', ylabel=r'log (integrated $N_{HI}$ [$cm^{-2}$])')
    
axs[1].set(xscale='log')


def N_H_of_b_model(b):
    b = b*u.deg.to(u.rad) 
    
    nh = 50
    R = 18*u.kpc.to(u.cm)
    Ro = 8.5*u.kpc.to(u.cm)
    del_l = 0.125*u.deg.to(u.rad)
    integ_b = 50-(-50)
    
    L_l = Ro*np.cos(l) + np.sqrt(R**2 - (Ro*np.sin(l))**2)
    
    N_H_of_l = nh*del_l*integ_b*L_l
    
    return N_H_of_l


fig,axs = plt.subplots(2,1,figsize=(10,10))
axs[0].plot(vels,np.log10(nhi_of_v*hi_conv))
axs[1].scatter(np.abs(vels),np.log10(nhi_of_v*hi_conv),s=10)

for ax in axs.flat:
    ax.set(xlabel=r'velocity [$km/s$]', ylabel=r'log (median $N_{HI}$ [$cm^{-2}$])')
    
axs[1].set(xscale='log')


def sinu(x):
    return 80*np.sin(0.75*np.pi*x)


def z(x):
    #underlying function is: 1/(1+(sinu(-X) - Y)**2)
    X,Y = np.meshgrid(x*u.deg,np.arange(-100,101,1.3))
    Z1 = 0.9/(1+ 0.5*(np.e**-5)*((np.e**(-np.abs(1e-4*X.value**2))*sinu(-X) - Y)**2))
    
    Z2 = 0.5/(1+ (np.e**-5)*((sinu(0.5*X)) - Y)**2)
    
    Z3 = 1/(1+ 0.2*(np.e**-4)*((0.05*sinu(2*X+90*u.deg)) - Y)**2)
    return Z1+Z2+Z3
    

ideal_lv = z(np.arange(360,-360,-0.125))


ideal_b = models.Lorentz1D(4.5e25,0.,5)(hclats)


np.array([ideal_lv.value]*800).swapaxes(0,1)


np.dot(np.array([ideal_lv.value]*800).swapaxes(0,1),ideal_b)


(np.outer(ideal_lv,ideal_b)).shape


(np.outer(ideal_b,ideal_lv)).shape


ideal_lv.shape,ideal_b.shape


plt.imshow(z(np.arange(180,-180,-0.125)),cmap=dc.cmap_min0('viridis'),
           extent=[180,-180,-100,100])
# plt.colorbar()
# plt.ylim(-3,3)


plt.imshow(hi_lv,extent=[360,-360,-100,100],cmap=dc.cmap_min0('viridis'),
           origin='lower')
plt.xlim(180,-180)


plt.imshow(np.log10(hi_conv*nhi_lb),
           extent=[360,-360,-50,50],origin='lower',vmin=16)
plt.colorbar(orientation='horizontal',
             label=r'mean $N_{HI}$ [$cm^{-2}$]')
plt.xlabel(r'velocity [km/s]')
plt.ylabel(r'galactic latitude [deg]')
plt.xlim(180,-180)
plt.show()


plt.imshow(np.log10(hi_conv*nhi_bv),
           extent=[-100,100,-50,50],origin='lower',vmin=16)
plt.colorbar(orientation='horizontal',
             label=r'mean $N_{HI}$ [$cm^{-2}$]')
plt.xlabel(r'velocity [km/s]')
plt.ylabel(r'galactic latitude [deg]')
plt.show()


plt.imshow(np.log10(hi_conv*nhi_lv),
           extent=[-100,100,-50,50],origin='lower',vmin=16)
plt.colorbar(orientation='horizontal',
             label=r'mean $N_{HI}$ [$cm^{-2}$]')
plt.xlabel(r'velocity [km/s]')
plt.ylabel(r'galactic latitude [deg]')
plt.show()


plt.figure(figsize=(15,15))
plt.imshow(np.log2(np.nanmean(hi4pi.data,axis=0)/nhi_lb_flat),
                      origin='lower',extent=[360,-360,-50,50],vmin=-5,vmax=5)
plt.xlim(180,-180)


plt.figure(figsize=(15,15))
plt.imshow(np.log2(
    np.nanmean(hi4pi.data,axis=0)/np.broadcast_to(
        np.nanmean(hi4pi.data,axis=(0,2)),(5760,800)).transpose()),
          origin='lower',extent=[360,-360,-50,50],vmin=-2,vmax=2)
plt.xlim(180,-180)


fig,axs = plt.subplots(3,1,figsize=(15,15),layout='constrained')

og = axs[0].imshow(np.log10(hi_lb*hi_conv),origin='lower',
                   extent=[360,-360,-50,50],vmin=18)
axs[0].contour(lons,lats,np.log10(co_lb*co_conv),levels=[18.5],colors='white',
               linewidths=0.7)
                      
model = axs[1].imshow(np.log10(hi_conv*nhi_lb),
           extent=[360,-360,-50,50],origin='lower',vmin=16,vmax=20)

final = axs[2].imshow(np.log2(np.nanmean(hi4pi.data,axis=0)/nhi_lb_flat),
                      origin='lower',extent=[360,-360,-50,50],vmin=-5,vmax=5)

axs[2].contour(lons,lats,np.log10(co_lb*co_conv),levels=[20],colors='white',
               linewidths=0.7)


for ax in axs.flat:
    ax.set(xlabel=r'galactic longitude [deg]',
          ylabel=r'galactic latitude [deg]',
          xlim=(180,-180))
    ax.label_outer()


fig.colorbar(og,    ax=axs[0], orientation='vertical',label=r'mean $N_{HI}$ [$cm^{-2}$]')
fig.colorbar(model, ax=axs[1], orientation='vertical',label=r'mean $N_{HI}$ [$cm^{-2}$]')
fig.colorbar(final, ax=axs[2], orientation='vertical',label=r'$log_2$ $[\ N_{HI} / \overline{N_{HI}}$ ]')


#^ add dust here, then show how to decompose?


np.outer(nhi_of_b,nhi_of_l)


plt.figure(figsize=(15,10))
plt.imshow(np.log10(np.nanmedian(hi4pi.data,axis=0)*hi_conv),origin='lower',extent=[360,-360,-50,50],vmin=18,vmax=20)

plt.colorbar()
plt.contour(lons,lats,np.log10(co_lb*co_conv),levels=[18.5],colors='white',linewidths=0.7)
plt.xlim(180,-180)


plt.figure(figsize=(15,10))
plt.imshow(np.log10(hi_lb*hi_conv),origin='lower',extent=[360,-360,-50,50],vmin=18,vmax=20)

plt.colorbar()
plt.contour(lons,lats,np.log10(co_lb*co_conv),levels=[18.5],colors='white',linewidths=0.7)
plt.xlim(180,-180)


np.outer(nhi_of_b,nhi_of_l)


plt.figure(figsize=(15,10))
plt.imshow(np.log10(nhi_lb),
           origin='lower',extent=[360,-360,-50,50])#,vmin=-4,vmax=1)
plt.colorbar()
plt.contour(lons,lats,np.log10(co_lb*co_conv),levels=[18.5],colors='white',linewidths=0.7)

plt.xlim(180,-180)


plt.figure(figsize=(15,10))
plt.imshow(np.log10(np.outer(nhi_of_b,nhi_of_l)),
           origin='lower',extent=[360,-360,-50,50])#,vmin=-4,vmax=1)
plt.colorbar()
plt.contour(lons,lats,np.log10(co_lb*co_conv),levels=[18.5],colors='white',linewidths=0.7)

plt.xlim(180,-180)


fig,axs = plt.subplots(1,3,figsize=(15,15),layout='constrained')

og = axs[0].imshow(np.log10(hi_conv*hi_bv),origin='lower',
                   extent=[-50,50,-100,100],vmin=18,vmax=21)
axs[0].contour(lats,vels,np.log10(co_bv*co_conv),levels=[18.5],colors='white',linewidths=0.7)

                      
model = axs[1].imshow(np.log10(hi_conv*nhi_bv),
           extent=[50,-50,-100,100],origin='lower',vmin=18,vmax=21)

final = axs[2].imshow(np.log2(hi_bv/nhi_bv),origin='lower',extent=[-50,50,-100,100],vmin=-2,vmax=2)
axs[2].contour(lats,vels,np.log10(co_bv*co_conv),levels=[18.5],colors='white',linewidths=0.7)


for ax in axs.flat:
    ax.set(ylabel=r'radial velocity [km/s]',
          xlabel=r'galactic latitude [deg]')
    ax.label_outer()


fig.colorbar(og,    ax=axs[0], orientation='horizontal',label=r'mean $N_{HI}$ [$cm^{-2}$]')
fig.colorbar(model, ax=axs[1], orientation='horizontal',label=r'mean $N_{HI}$ [$cm^{-2}$]')
fig.colorbar(final, ax=axs[2], orientation='horizontal',label=r'$log_2$ $[\ N_{HI} / \overline{N_{HI}}$ ]')


# fig,axs = plt.subplots(1,3,figsize=(15,15),layout='constrained')

# og = axs[0].imshow(np.log10(hi_conv*hi_lv),origin='lower',
#                    extent=[360,-360,-100,100],vmin=18,vmax=21)
# axs[0].contour(lons,vels,np.log10(co_lv*co_conv),levels=[18.5],colors='white',linewidths=0.7)

                      
# # model = axs[1].imshow(np.log10(hi_conv*nhi_bv),
# #            extent=[50,-50,-100,100],origin='lower',vmin=18,vmax=21)

# # final = axs[2].imshow(np.log2(hi_bv/nhi_bv),origin='lower',extent=[-50,50,-100,100],vmin=-2,vmax=2)
# # axs[2].contour(lats,vels,np.log10(co_bv*co_conv),levels=[18.5],colors='white',linewidths=0.7)


# for ax in axs.flat:
#     ax.set(ylabel=r'radial velocity [km/s]',
#           xlabel=r'galactic longitude [deg]',
#           xlim=(180,-180))
#     ax.label_outer()


# # fig.colorbar(og,    ax=axs[0], orientation='horizontal',label=r'mean $N_{HI}$ [$cm^{-2}$]')
# # fig.colorbar(model, ax=axs[1], orientation='horizontal',label=r'mean $N_{HI}$ [$cm^{-2}$]')
# # fig.colorbar(final, ax=axs[2], orientation='horizontal',label=r'$log_2$ $[\ N_{HI} / \overline{N_{HI}}$ ]')
