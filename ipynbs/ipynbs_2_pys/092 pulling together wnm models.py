#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import datacube as dc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
import fun as fns

import astropy.units as u
import astropy.constants as c

from matplotlib.patches import Circle,Rectangle


hi4pi,dht = fits.open('both_same_grid_720d.fits')
zcat = QTable(fits.open('final_cat.fits')[1].data)

lons = np.arange(5760)*-0.125 + 360
lats = np.arange(800)*0.125 - 50
vels = np.arange(155)*1.3 - 100


hi_conv = 1.828e18*1.3#((u.cm)**-2)/(u.K*(u.km/u.s)) \
#            *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = 1.97e20 *1.3# ((u.cm)**-2)/(u.K*(u.km/u.s))\
#            *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


from skimage.filters import gaussian


hi4pi = hi4pi.data*(hi4pi.data > 0)


cmom0 = np.nansum(dht.data,0)*1.3*co_conv


lons = np.arange(5760)*-0.125 + 360
lats = np.arange(800)*0.125 - 50


lbmn = np.nanmean(hi4pi,0)
lvmn = np.nanmean(hi4pi,1)
bvmn = np.nanmean(hi4pi,2)


lbint = np.nansum(hi4pi,0)
lvint = np.nansum(hi4pi,1)
bvint = np.nansum(hi4pi,2)


kernel_smoothed = gaussian(hi4pi,80)


ker_sm_lb = np.nansum(kernel_smoothed,0)


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lbint*hi_conv),vmin=18.5,cmap='viridis',
           extent=[360,-360,-50,50],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10(lbint*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('moment 0 map, raw')
plt.show()


plt.figure(figsize=(20,10))
plt.imshow(np.log10((ker_sm_lb)*hi_conv),
           extent=[360,-360,-50,50],origin='lower',vmin=18.5)

plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
# plt.contour(lons,lats,np.log10((lbint-ker_sm_lb)*hi_conv),
#                                levels=[np.log10(4e20)],
#             colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('gaussian smoothed (kernel size = 80)')
plt.show()


plt.figure(figsize=(20,10))
plt.imshow(np.log10((lbint-ker_sm_lb)*hi_conv),
           extent=[360,-360,-50,50],origin='lower',vmin=18.5)

plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10((lbint-ker_sm_lb)*hi_conv),
                               levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('subtracting gaussian smoothed image')
plt.show()


smooth_lb = np.broadcast_to(np.nanmedian(lb,1),(lb.T.shape)).T


plt.plot(lats[400:],np.nanmean(lb*hi_conv,1)[400:])

plt.xlabel('log ( |galactic latitude [deg]| )')
plt.ylabel(r'log( $\widetilde{N_{HI}}$ [$cm^{-2}$] )')
plt.loglog()
plt.show()


plt.figure(figsize=(20,10))
plt.imshow(np.log10(smooth_lb*hi_conv),cmap='viridis',vmin=18.5,
           extent=[360,-360,-50,50],origin='lower')

plt.colorbar(orientation='horizontal',label=r'$N_{HI}$ /$\widetilde{N_{HI}}$')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)

# plt.contour(lons,lats,lbmn/smooth_lb,levels=[1],
#             colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title(r'median $N_{HI}$ along b')
plt.show()


plt.figure(figsize=(20,10))
plt.imshow(lbmn/smooth_lb,cmap='RdYlBu_r',
           extent=[360,-360,-50,50],origin='lower',vmin=0,vmax=2)

plt.colorbar(orientation='horizontal',label=r'$N_{HI}$ /$\widetilde{N_{HI}}$')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)

# plt.contour(lons,lats,lbmn/smooth_lb,levels=[1],
#             colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title(r'normalizing by median $N_{HI}$ along b')
plt.show()


mom0_fft_i = np.fft.fftn(lbint)
mom0_fft = np.abs(np.fft.fftshift(mom0_fft_i)**2)


plt.plot(np.fft.fftfreq(2880),np.nanmean(mom0_fft,0)[2880:],
         label='spatial power spectrum')
plt.vlines(np.fft.fftfreq(2880)[50],1e6,1e15,color='red',linestyle='--')
# plt.plot(np.fft.fftfreq(400),np.nanmean(mom0_fft,1)[400:],
#         label = 'spatial power spectrum along b')

plt.loglog()
plt.ylabel('power')
plt.xlabel('wavenumber [1/pix]')
plt.legend()
# plt.xlim(np.fft.fftfreq(2880)[50],np.fft.fftfreq(2880)[2000])
plt.show()


fft_msk_l = np.zeros((800,5760))
fft_msk_l[:,2870:2890] = np.ones((800,20))


l_masked = np.abs(np.fft.ifftn(np.fft.ifftshift(fft_msk_l*np.fft.fftshift(mom0_fft_i))))


plt.figure(figsize=(20,10))
plt.imshow(np.log10(l_masked*hi_conv),vmin=18.5,
           extent=[360,-360,-50,50],origin='lower')

plt.colorbar(orientation='horizontal',label=r'log ( $N_{HI}$ [$cm^{-2}$])')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)

# plt.contour(lons,lats,(lbint-b_masked)*hi_conv,levels=[4e20],
#             colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title(r'lowest frequency mode')
plt.show()


plt.figure(figsize=(20,10))
plt.imshow(np.log10((lbint-l_masked)*hi_conv),vmin=18.5,
           extent=[360,-360,-50,50],
           origin='lower')

plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6,zorder=100)

plt.contour(lons,lats,np.log10((lbint-l_masked)*hi_conv),levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title(r'subtracting lowest 10 frequency modes')
plt.show()


lv = np.nanmean(hi4pi.data,1)


lvint_posb = np.nansum(hi4pi.data[:,400:,:],1)
lvint_negb = np.nansum(hi4pi.data[:,:400,:],1)


lvint = np.nansum(hi4pi.data,1)


clv = np.nansum(dht.data,1)*co_conv


clv_posb = np.nansum(dht.data[:,400:,:],1)*co_conv
clv_negb = np.nansum(dht.data[:,:400,:],1)*co_conv


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lvint_negb*hi_conv),vmin=19,cmap='viridis',
           extent=[360,-360,-100,100],origin='lower')

plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')

plt.contour(lons,vels,clv_negb,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,vels,np.log10(lvint_negb*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.hlines([-25,25],180,-180)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('velocity [km/s]')
plt.title('l-v, b < 0')
plt.show()


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lvint_posb*hi_conv),vmin=19,cmap='viridis',
           extent=[360,-360,-100,100],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,vels,clv_posb,levels=[1e21],
            colors='black',linewidths=0.6)

plt.contour(lons,vels,np.log10(lvint_posb*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.hlines([-25,25],180,-180)
plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('velocity [km/s]')
plt.title('l-v, b >= 0')
plt.show()


ker_sm_lv_posb = np.nansum(kernel_smoothed[:,400:,:],1)
ker_sm_lv_negb = np.nansum(kernel_smoothed[:,:400,:],1)


plt.figure(figsize=(20,10))
plt.imshow(np.log10((lvint_posb-ker_sm_lv_posb)*hi_conv),vmin=18.5,
           extent=[360,-360,-100,100],origin='lower')

plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,vels,clv_posb,levels=[1e21],
            colors='black',linewidths=0.6)

plt.contour(lons,vels,np.log10((lvint_posb-ker_sm_lv_posb)*hi_conv),
                               levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('gaussian smoothed (kernel size = 80)')
plt.show()


vels[58]


vels[96]


lb_pos = np.nansum(hi4pi.data[:58,:,:],0)
lb_0 = np.nansum(hi4pi.data[58:96,:,:],0)
lb_neg = np.nansum(hi4pi.data[96:,:,:],0)


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lb_pos*hi_conv),vmin=18.5,cmap='viridis',
           extent=[360,-360,-50,50],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10(lb_pos*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('moment 0 map, raw')
plt.show()


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lb_0*hi_conv),vmin=18.5,cmap='viridis',
           extent=[360,-360,-50,50],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10(lb_0*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)
plt.contour(lons,lats,np.log10(lbint*hi_conv),
            levels=[np.log10(4e20)],
            colors='red',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('moment 0 map, raw')
plt.show()


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lb_neg*hi_conv),vmin=18.5,cmap='viridis',
           extent=[360,-360,-50,50],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10(lb_neg*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('moment 0 map, raw')
plt.show()


plt.plot(vels,np.log10(np.nanmean(lv*hi_conv,1)))

plt.xlabel('velocity [km/s]')
plt.ylabel(r'log( $\widetilde{N_{HI}}$ [$cm^{-2}$] )')
plt.show()


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lbint*hi_conv),vmin=18.5,cmap='viridis',
           extent=[360,-360,-50,50],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10(lbint*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(180,-180)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('moment 0 map, raw')
plt.show()


plt.figure(figsize=(15,7))
plt.imshow(np.log10(lbint*hi_conv),vmin=18.5,cmap='viridis',
           extent=[360,-360,-50,50],origin='lower')
plt.colorbar(orientation='horizontal',label=r'log( $N_{HI}$ [$cm^{-2}$] )')
plt.contour(lons,lats,cmom0,levels=[1e21],
            colors='black',linewidths=0.6)
plt.contour(lons,lats,np.log10(lbint*hi_conv),
            levels=[np.log10(4e20)],
            colors='white',linewidths=0.6)

plt.xlim(360,0)
plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitude [deg]')
plt.title('moment 0 map, raw')
plt.show()
