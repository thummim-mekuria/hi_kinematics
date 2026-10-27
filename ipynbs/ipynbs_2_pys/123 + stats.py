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
import fun as fns


hi4pi,dht = fits.open('both_same_grid_720d.fits')


goodh = hi4pi.data*(hi4pi.data>0)


lvh = np.nansum(goodh,1)


hmom0 = np.nansum(goodh,0)
# hmom1 = fns.get_mom1(goodh,[0,155]) 
# hmom2 = fns.get_momN(goodh,[0,155],2) 
# hmom3 = fns.get_momN(goodh,[0,155],3) 
# hmom4 = fns.get_momN(goodh,[0,155],4) 


# h_bmean = np.broadcast_to(np.nanmean(hmom0,1),(5760,800)).T


hcl = np.arange(5760)*-0.125+360
hcb = np.arange(800)*0.125-50
hcv = np.arange(155)*1.3-100


cmom0 = np.nansum(dht.data,0)


plt.imshow(np.log2(hmom0/h_bmean), origin='lower',extent=[360,-360,-50,50],vmin=-1.5,vmax=1.5)
plt.xlim(180,-180)
# plt.colorbar()


fig,axs = plt.subplots(5,1,figsize=(20,10),layout='constrained')
im0=axs[0].imshow(hmom0,origin='lower',extent=[360,-360,-50,50],vmin=0,vmax=800)
im1=axs[1].imshow(hmom1,origin='lower',extent=[360,-360,-50,50],norm='asinh')
im2=axs[2].imshow(np.sqrt(hmom2),origin='lower',extent=[360,-360,-50,50],vmin=0,vmax=30)
im3=axs[3].imshow(hmom3,origin='lower',extent=[360,-360,-50,50],norm='asinh')
im4=axs[4].imshow(hmom4-3,origin='lower',extent=[360,-360,-50,50],vmax=20)


for ax,im in zip(axs.flat,[im0,im1,im2,im3,im4]):
    ax.set(xlim=(180,-180))
    ax.contour(hcl,hcb,cmom0,levels=[2],colors='white',linewidths=0.5)
    ax.label_outer()
#     fig.colorbar(im,ax=ax,orientation='horizontal')
# plt.colorbar()


lv_fft_i = np.fft.fftn(lvh)
lv_fft = np.abs(lv_fft_i)**2


lv_fft.shape


plt.plot(np.nanmean(np.fft.fftshift(lv_fft),0)[2880:])
# plt.plot(np.nanmean(np.fft.fftshift(lv_fft),1)[77:])
# plt.yscale('log')
plt.loglog()
plt.vlines(10,1e7,1e17,colors='red')


plt.figure(figsize=(20,10))
plt.imshow(lvh,norm='asinh',extent=[360,-360,-100,100],origin='lower')

plt.xlim(180,-180)


lvc = np.nansum(dht.data,1)


lvh_mask = np.ones(lvh.shape)
lvh_mask[:,2870:2890] = np.zeros((155,20))


2800+880


lvh_mask = np.zeros(lvh.shape)
lvh_mask[:,2870:2890] = np.ones((155,20))
# lvh_mask[:,700:1300] = np.ones((155,600))


lvh_new_i = np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(lv_fft_i)*lvh_mask))
lvh_new = np.abs(lvh_new_i)

plt.figure(figsize=(20,10))
plt.imshow(lvh-lvh_new,norm='asinh',extent=[360,-360,-100,100],origin='lower')
plt.contour(hcl,hcv,lvc,levels=[5],colors='white',linewidths=1,alpha=0.7)
# plt.contour(hcl,hcb,back2,levels=[400],colors='red',linewidths=0.4)

plt.xlim(180,-180)
# plt.colorbar()


lb_fft_i = np.fft.fftn(hmom0)
lb_fft = np.abs(lb_fft_i)**2


lbh_mask = np.zeros(hmom0.shape)
lbh_mask[:,2870:2890] = np.ones((800,20))


lbh_new_i = np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(lb_fft_i)*lbh_mask))
lbh_new = np.abs(lbh_new_i)

plt.figure(figsize=(20,10))
plt.imshow(hmom0-lbh_new,norm='asinh',extent=[360,-360,-50,50],origin='lower')
plt.contourf(hcl,hcb,cmom0,levels=[2,5],colors='red',linewidths=0.4,alpha=0.5)
# plt.contour(hcl,hcb,back2,levels=[400],colors='red',linewidths=0.4)

plt.xlim(200,-200)
# plt.colorbar()


hcl[1440]


1440+2880


all_spec = goodh[:,:,1440:4320].T.reshape(-1,155)


al_spc_fft = np.fft.fftn(all_spec)


al_spc_fft_r = np.abs(al_spc_fft)**2


plt.imshow(all_spec,extent=[-100,100,-100,100],norm='asinh')


plt.plot(np.nanmean(np.fft.fftshift(al_spc_fft_r),0)[77:])
# plt.yscale('log')
plt.loglog()


plt.plot(np.nanmean(np.fft.fftshift(al_spc_fft_r),1)[1152000:])
# plt.plot(np.nanmean(np.fft.fftshift(lv_fft),1)[77:])
# plt.yscale('log')
plt.loglog()
# plt.vlines(10,1e7,1e17,colors='red')


allsp_mask.shape


allsp_mask = np.ones(all_spec.shape)
allsp_mask[1151000:1153000,:] = np.zeros((2000,155))


allsp_new_i = np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(al_spc_fft)*allsp_mask))
allsp_new = np.abs(allsp_new_i)


plt.imshow(all_spec,extent=[-100,100,-100,100],norm='asinh')


nooh = allsp_new.reshape(2880,800,155).T


plt.imshow(all_spec-allsp_new,extent=[-100,100,-100,100],norm='asinh')


test_spec = goodh[:,400:401,1509:1510].T[0][0]


plt.plot(test_spec)


eg_spc_fft = np.fft.fft(test_spec)


plt.plot(np.abs(eg_spc_fft)**2)
plt.loglog()


msk = np.zeros(len(eg_spc_fft))
msk[0:10] = np.ones(10)


new_spc = np.abs(np.fft.ifft(msk*eg_spc_fft))


plt.plot(test_spec)
plt.plot(new_spc)

plt.plot(test_spec-new_spc)
# plt.plot(test_spec-new_spc)


test_spec2 = np.nanmean(goodh[:,400:410,2509:2530],(1,2))


plt.plot(test_spec2)


eg2_spc_fft = np.fft.fft(test_spec2)


plt.plot(np.abs(eg2_spc_fft)**2)
plt.loglog()


new_spc2 = np.abs(np.fft.ifft(msk*eg2_spc_fft))


plt.plot(test_spec2)
plt.plot(new_spc2)

plt.plot(test_spec2-new_spc2)
# plt.plot(test_spec-new_spc)


hmom1 = fns.get_mom1(goodh,[0,155]) 


hmom2 = fns.get_momN(goodh,[0,155],2) 


plt.figure(figsize=(20,10))
plt.imshow(lvh,norm='asinh',extent=[360,-360,-100,100],origin='lower')
# plt.contour(hcl,hcv,lvh,levels=np.arange(2e2,2e4,9e2))

[plt.scatter((np.arange(5760)*-0.125 + 360),hmom1[i],s=1e-3*hmom2[i]) 
 for i in range(0,800,37)]
plt.xlim(180,-180)


#okay pivot to: generating l,v plots for a range of Rs & bs


Rs = (np.concatenate((np.arange(1e3,7e3,5e2),
                np.arange(6.5e3,10.5e3,1e2),
                np.arange(1.05e4,13e3,5e2))))*u.pc


lons = np.arange(180,-180,-0.125)*u.deg
vels = (np.arange(155)*1.3 - 100)*u.km/u.s


def vs_of_lnR(l,R):
    return np.nanmean([np.cos(0)*np.sin(l)*
            (R/(8.5*u.kpc) - 1)*220*u.km/u.s for i in range(len(Rs))],axis=0)


v,l,R = np.meshgrid(vels,lons,Rs)


rot_curves_lvr = [1/(1+(v[:,:,i].value - vs_of_lnR(l[:,:,i],R[:,:,i]))**2) for i in range(57)]


plt.imshow(np.nansum(rot_curves_lvr,0).T,norm='asinh',extent=[180,-180,-100,100])
# plt.contour(hcl,hcv,lvh,levels=[1e3,1e4],colors='white')

plt.xlim(180,-180)
plt.show()
# plt.colorbar()


plt.figure(figsize=(20,10))
# plt.imshow(lvh,norm='asinh',extent=[360,-360,-100,100],origin='lower')
plt.contourf(hcl,hcv,lvh,levels=np.log10(np.arange)

# [plt.scatter((np.arange(5760)*-0.125 + 360),hmom1[i],s=1e-3*hmom2[i]) 
#  for i in range(0,800,37)]
plt.xlim(180,-180)


2**(np.arange(10,15))


plt.imshow(np.nansum(goodh[:,:,1440:4320]-nooh,2),origin='lower',norm='asinh',
           extent=[180,-180,-50,50])


plt.imshow(np.nansum(nooh,2),origin='lower',norm='asinh',
           extent=[180,-180,-50,50])


plt.figure(figsize=(20,10))

plt.imshow(hmom0[:,1440:4320]-np.nansum(nooh,0),origin='lower',vmin=0,
           extent=[180,-180,-50,50])
# plt.colorbar()
# plt.contour(hcl,hcb,cmom0,levels=[2],colors='white',linewidths=0.7)
plt.contour(hcl[1440:4320],hcb,hmom0[:,1440:4320]-np.nansum(nooh,0),
            levels=[400],colors='white',linewidths=0.7)

plt.xlim(180,-180)


plt.imshow(np.nansum(nooh,1),origin='lower',vmin=0,
           extent=[180,-180,-100,100])
plt.contour(hcl[1440:4320],hcv,np.nansum(nooh,1),
            levels=[2300],colors='white',linewidths=0.7)
# plt.contour(hcl,hcv,lvh,
#             levels=[1900],colors='orange',linewidths=0.7)
plt.xlim(180,-180)


import matplotlib.image as mpimg
import scipy.stats as scistat

# image = mpimg.imread("clouds.png")
image = hmom0
npix0,npix1 = image.shape

fourier_image = np.fft.fftn(image)
fourier_amplitudes = np.abs(fourier_image)**2


nfa1 = np.nanmean(fourier_amplitudes,axis=1)
fs1 = np.fft.fftfreq(npix0)


nfa0 = np.nanmean(fourier_amplitudes,axis=0)
fs0 = np.fft.fftfreq(npix1)


plt.plot(fs1[:int(npix0/2)],nfa1[:int(npix0/2)])
plt.plot(fs0[:int(npix1/2)],nfa0[:int(npix1/2)],alpha=0.5)
# plt.yscale('log')
plt.loglog()


# back = np.fft.ifftn(fourier_image*np.abs((msk-1)))
# back2 = np.abs(back)
plt.imshow(np.fft.fftshift(fourier_amplitudes),norm='asinh',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.colorbar()


plt.plot(np.nanmean(np.fft.fftshift(fourier_amplitudes),1)[400:])
# plt.yscale('log')
plt.loglog()
# plt.xlim(400,)
plt.vlines([230,350],1e8,1e15)


np.all(fourier_amplitudes*(fourier_amplitudes<1e-2)==0)


msk = np.zeros(hmom0.shape)
msk[630:750,:] = np.ones((120,5760))


back = np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(fourier_image)*np.abs(msk-1)))
back2 = np.abs(back)
plt.figure(figsize=(20,10))
plt.imshow(back2,norm='asinh',extent=[360,-360,-50,50],origin='lower')
plt.contour(hcl,hcb,cmom0,levels=[2],colors='white',linewidths=0.4)
plt.contour(hcl,hcb,back2,levels=[400],colors='red',linewidths=0.4)

plt.xlim(180,-180)
# plt.colorbar()


back = np.fft.ifftn(fourier_image*msk)
back2 = np.abs(back)
plt.imshow(back2,norm='asinh',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.colorbar()


plt.plot(np.nanmean(np.fft.fftshift(fourier_amplitudes),0)[2880:])
# plt.yscale('log')
# plt.xlim(2880,)
plt.loglog()
plt.vlines([10],1e8,1e15,color='red')


msk = np.zeros(hmom0.shape)
msk[:,2870:2890] = np.ones((800,20))


back = np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(fourier_image)*np.abs(msk-1)))
back2 = np.abs(back)
plt.figure(figsize=(20,10))
plt.imshow(back2,norm='asinh',extent=[360,-360,-50,50],origin='lower')
plt.contour(hcl,hcb,cmom0,levels=[2],colors='white',linewidths=0.4)
plt.contour(hcl,hcb,back2,levels=[200],colors='red',linewidths=0.4)

plt.xlim(180,-180)
# plt.colorbar()


plt.plot(np.nanmean(back2,(1))*100)
plt.plot(np.nanmean(hmom0,(1)))


from turbustat.statistics import PDF


hi_conv = 1.3*1.828e18#*((u.cm)**-2)/(u.K*(u.km/u.s)) *1.36*c.m_p)#.to(u.solMass/u.pc**2) 
co_conv = 1.3*1.97e20 #* ((u.cm)**-2)/(u.K*(u.km/u.s))\
# gas2dust = 1.87e21#*(u.cm**-2)/u.mag #from zhu et al: https://arxiv.org/pdf/1706.07109.pdf


PDF(hmom0*hi_conv).run().plot_distrib()


plt.scatter(np.ravel(hmom0/h_bmean)[::39],(np.ravel(np.sqrt(hmom2)))[::39],s=1e-1)
plt.loglog()


plt.scatter((np.ravel(np.sqrt(hmom2)))[::39],np.ravel(hmom3)[::39],s=1e-1)
plt.loglog()


plt.scatter((np.ravel(np.sqrt(hmom2)))[::39],np.ravel(cmom0)[::39],s=1e-1)
plt.loglog()


plt.scatter((np.ravel(hmom4)-3)[::39],np.ravel(cmom0)[::39],s=1e-1)
plt.loglog()


plt.scatter((np.ravel(hmom4)-3)[::39],(np.ravel(np.sqrt(hmom2)))[::39],s=1e-1)
plt.loglog()


from scipy.ndimage import gaussian_filter1d


np.sqrt(hmom2[400,3000])


hmom3[400,3000]


hmom4[400,3000]


np.nanstd(goodh[:,400,3000])


plt.plot(goodh[:,700,1000])


plt.plot(goodh[:,100,3000])


plt.plot(goodh[:,400,3000])
plt.plot(gaussian_filter1d(goodh[:,400,3000],30))
# plt.plot((gaussian_filter1d(goodh[:,400,3000],7)-goodh[:,400,3000]))

# plt.plot(dht.data[:,400,3000]*30)


np.sinh(1e-4*wh)


plt.imshow(wh,norm='asinh')
plt.colorbar()


from fil_finder import FilFinder2D


hifil = FilFinder2D(hmom0)


hifil.preprocess_image(flatten_percent=50)


plt.imshow(hifil.flat_img.value, origin='lower',extent=[360,-360,-50,50])
plt.xlim(180,0)
plt.colorbar()


hifil.create_mask(verbose=True)


h_bmean = np.broadcast_to(np.nanmean(hmom0,1),(5760,800)).T


plt.imshow(np.log2(hmom0/h_bmean), origin='lower',extent=[360,-360,-50,50],vmin=-1.5,vmax=1.5)
plt.xlim(180,-180)
# plt.colorbar()


#still trying to figure out a toy milky way hi pdf


nhb = np.nanmean(hi4pi.data,axis=(0,2))


nhlv = np.nanmean(hi4pi.data,axis=(1))


hi4pi.data.shape


np.array([np.ones((15,576))]*80).shape


plt.imshow(nhlv[::10,::10],extent=[360,-360,-100,100])


bb_hi = np.outer(nhlv[::10,::10],nhb[::10]).reshape(16,576,80).swapaxes(1,2)


plt.plot(np.nanmean(np.nanmean(bb_hi,axis=0),axis=0))


plt.plot(np.nanmean(bb_hi,axis=(1,2)))


plt.plot(np.nanmean(bb_hi,axis=(0,2)))


plt.imshow(np.nansum(bb_hi,axis=2),extent=[-50,50,-100,100],cmap=dc.cmap_min0('viridis'))
# plt.xlim(180,-180)


plt.imshow(np.nansum(bb_hi,axis=0),extent=[360,-360,-50,50],norm='asinh')
plt.xlim(180,-180)
# plt.colorbar()


from skimage.filters import butterworth


from skimage.filters import gaussian


hismo =gaussian(hi4pi.data,80)


plt.imshow(np.nansum(gaussian(hi4pi.data,80),axis=0),
           extent=[360,-360,-50,50],norm='asinh')
plt.xlim(180,-180)
# plt.colorbar()


plt.imshow(np.nansum(bb_hi_model,axis=0),extent=[360,-360,-50,50],norm='asinh')
plt.xlim(180,-180)


-0.05527484/0.10527779


bb_hi_model = bb_hi


bb_hi = hi4pi.data[::10,::10,::10]


np.array([nhlv[::10,::10]]).shape


737280/(576*16*80)


plt.imshow(np.nanmean(hi4pi.data,axis=(1)),extent=[360,-360,-100,100])
