#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

from astropy.io import fits
import astropy.constants as c
import astropy.units as u
from astropy.table import QTable
from astropy.modeling import models,fitting

import datacube as dc
import fun as fns

import matplotlib.pyplot as plt


hi4pi,dht = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


h_3_35 = (hi4pi.data * (hi4pi.data > 0.043*3))[50:105,:,:]
d_3_35 = (dht.data * (dht.data > 0.2*3))[50:105,:,:]


h2 = fns.get_momN(h_3_35,np.arange(55)*1.3 - 35,2)
c2 = fns.get_momN(d_3_35,np.arange(55)*1.3 - 35,2)


h1 = fns.get_mom1(h_3_35,np.arange(55)*1.3 - 35)
c1 = fns.get_mom1(d_3_35,np.arange(55)*1.3 - 35)


h0 = np.nansum(h_3_35,axis=0)*1.3*1.82e18
c0 = np.nansum(d_3_35,axis=0)*1.3*1.97e20


plt.figure(figsize=(10,10))
plt.contourf(h0,levels=[5e20,1e21,2e21,3e21,4e21,5e21,1e22],
             norm='log',extent=[360,-360,-50,50])

plt.xlim(180,150)
plt.ylim(-30,20)
plt.colorbar()

plt.contour(c0,levels=[5e20,1e21,2e21,3e21,4e21,5e21,1e22,5e22,1e23],
             norm='log',extent=[360,-360,-50,50],linewidths=0.7,colors='white')


plt.figure(figsize=(10,10))
plt.contourf(np.sqrt(h2),levels=np.arange(5,15,1),
             extent=[360,-360,-50,50])

plt.xlim(180,150)
plt.ylim(-30,0)
plt.colorbar()

plt.contour(c0,levels=[1e21,2e21],
             norm='log',extent=[360,-360,-50,50],linewidths=0.8,colors='white')


plt.figure(figsize=(10,10))
plt.contourf(np.sqrt(h2),levels=np.arange(5,15,1),
             extent=[360,-360,-50,50])

plt.xlim(230,190)
plt.ylim(-30,0)
plt.colorbar()

plt.contour(c0,levels=[1e21,2e21],
             norm='log',extent=[360,-360,-50,50],linewidths=0.8,colors='white')


plt.figure(figsize=(10,10))
plt.contourf(np.sqrt(h2),levels=np.arange(5,20,2),
             extent=[360,-360,-50,50])

plt.xlim(100,85)
plt.ylim(-20,0)
plt.colorbar()

plt.contour(c0,levels=[5e20,1e21,2e21],
             norm='log',extent=[360,-360,-50,50],linewidths=0.8,colors='white')


plt.figure(figsize=(10,10))
plt.contourf(np.sqrt(h2),levels=np.arange(5,20,2),
             extent=[360,-360,-50,50])

plt.xlim(140,110)
plt.ylim(20,30)
plt.colorbar()

plt.contour(c0,levels=[5e20,1e21,2e21],
             norm='log',extent=[360,-360,-50,50],linewidths=0.6,colors='white')


h3 = fns.get_momN(h_3_35,np.arange(55)*1.3 - 35,3)


plt.hist(np.ravel(h3),bins=50)
plt.plot(np.arange(-2,2,0.5),models.Gaussian1D(5e5,-0.2,0.4)(np.arange(-2,2,0.5)),color='red')
plt.yscale('log')


plt.figure(figsize=(20,8))
plt.contourf(h3,levels=[-2,-1,-0.05,0.05,1,2],cmap='RdYlGn',
             extent=[360,-360,-50,50])

plt.xlim(180,-180)
# plt.ylim(20,30)
plt.colorbar()

plt.contour(h1,levels=np.arange(-30,30,10),
            extent=[360,-360,-50,50],linewidths=0.6,cmap='RdYlGn')
plt.colorbar()


plt.scatter(np.ravel(h1)[::100],np.ravel(h3)[::100],s=1e-2,)
plt.hlines(0,-30,30)
plt.vlines(0,-3,3)


for i in range(17):    
    loh,lih = zcat['l_pix_HI'][i]
    loc,lic = zcat['l_pix_GMC'][i]
    
    boh,bih = zcat['b_pix_HI'][i]
    boc,bic = zcat['b_pix_GMC'][i]
    
    X,Y = np.array(np.meshgrid(np.arange(loh,lih),np.arange(boh,bih)))
    xo,yo = np.nanmean([loh,lih]),np.nanmean([boh,bih])
    Z = np.sqrt((X-xo)**2 + (Y-yo)**2)*0.125*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    rs = np.ravel(Z)
    
    
    fig,axs = plt.subplots(1,3,figsize=[15,5])
    
    axs[0].contourf(h1[boh:bih,loh:lih])
    axs[0].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    axs[1].contourf(h2[boh:bih,loh:lih])
    axs[1].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    axs[2].scatter(rs,np.ravel(np.sqrt(h2[boh:bih,loh:lih])),s=1e-1,alpha=0.5)
    
    
    plt.show()


for i in range(17):    
    loh,lih = zcat['l_pix_HI'][i]
    loc,lic = zcat['l_pix_GMC'][i]
    
    boh,bih = zcat['b_pix_HI'][i]
    boc,bic = zcat['b_pix_GMC'][i]
    
    X,Y = np.array(np.meshgrid(np.arange(loh,lih),np.arange(boh,bih)))
    xo,yo = np.nanmean([loh,lih]),np.nanmean([boh,bih])
    Z = np.sqrt((X-xo)**2 + (Y-yo)**2)*0.125*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    rs = np.ravel(Z)
    
    
    r = int(np.nanstd([loh,lih]))
    X,Y = np.array(np.meshgrid(np.arange(-r+loh,lih+r),
                               np.arange(np.nanmax([-r+boh,0]),np.nanmin([bih+r,800]))))
    
    xo,yo = np.nanmean([loh,lih]),np.nanmean([boh,bih])
    Z = np.sqrt((X-xo)**2 + (Y-yo)**2)*0.125*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    r2s = np.ravel(Z)
    
    
    X,Y = np.array(np.meshgrid(np.arange(-2*r+loh,lih+2*r),
                               np.arange(np.nanmax([-2*r+boh,0]),np.nanmin([bih+2*r,800]))))
    Z = np.sqrt((X-xo)**2 + (Y-yo)**2)*0.125*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    r3s = np.ravel(Z)
    
    
    fig,axs = plt.subplots(1,3,figsize=[15,5])
    
    axs[0].contourf(np.sqrt(h2)[boh:bih,loh:lih])
    axs[0].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    axs[1].contourf(h2[np.nanmax([-r+boh,0]):np.nanmin([bih+r,800]),-r+loh:lih+r])
    axs[1].contour(c0[np.nanmax([-r+boh,0]):np.nanmin([bih+r,800]),-r+loh:lih+r],levels=[5e20])
  
    
    axs[2].scatter(r3s[::200],np.ravel(np.sqrt(h2[np.nanmax([-2*r+boh,0]):np.nanmin([bih+2*r,800]),
                                                  -2*r+loh:lih+2*r]))[::200],s=10,alpha=0.5)
    
    axs[2].scatter(r2s[::200],np.ravel(np.sqrt(h2[np.nanmax([-r+boh,0]):np.nanmin([bih+r,800]),
                                                  -r+loh:lih+r]))[::200],s=10,alpha=0.5)
    axs[2].scatter(rs[::200],np.ravel(np.sqrt(h2[boh:bih,
                                                 loh:lih]))[::200],s=10,alpha=0.5)
    
    
    plt.show()


#from a histogram, can we see that the average sigma for the hi that is associated/within the co contour or ra
#is at a smaller linewidth than outside
#h-alpha, ysos --> for tracing star formation activity, look at what properties of the clouds might be connected
#to this

#


plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(h2),extent=[360,-360,-50,50],
             levels=6,norm='asinh')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],linewidths=0.5,
            levels=[5e20],colors='white')
plt.xlim(360,0)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))
plt.contourf(bav,extent=[360,-360,-50,50],
             levels=6,norm='asinh')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],linewidths=0.5,
            levels=[5e20],colors='white')
plt.xlim(360,0)
# plt.ylim(-30,30)


bav = np.broadcast_to(np.nanmean(np.sqrt(h2),axis=1),(5760,800)).transpose()
plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(h2)-bav,extent=[360,-360,-50,50],
             levels=6,norm='asinh')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],linewidths=0.5,
            levels=[5e20],colors='white')
plt.xlim(360,0)
# plt.ylim(-30,30)


for i in range(17):    
    loh,lih = zcat['l_pix_HI'][i]
    loc,lic = zcat['l_pix_GMC'][i]
    
    boh,bih = zcat['b_pix_HI'][i]
    boc,bic = zcat['b_pix_GMC'][i]
    
    X,Y = np.array(np.meshgrid(np.arange(loh,lih),np.arange(boh,bih)))
    xo,yo = np.nanmean([loh,lih]),np.nanmean([boh,bih])
    Z = np.sqrt((X-xo)**2 + (Y-yo)**2)*0.125*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    rs = np.ravel(Z)
    
    d10 = ((10/zcat['distance<pc>'][i])*u.rad.to(u.deg))/0.125
    #now want i want to do is use this ^ (=how many pixels = 10pc at the cloud distance)
    #to segment the fov & average the spectra within those boxes
    #along the lines of contreras et al
    
    fig,axs = plt.subplots(1,2,figsize=[11,6])
    
    axs[0].contourf(h1[boh:bih,loh:lih])
    axs[0].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    axs[1].contourf(h2[boh:bih,loh:lih])
    axs[1].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    for inc,dec in zip(range(int(-0.5*(lih-loh)//d10),int(0.5*(lih-loh)//d10)),
                 range(int(0.5*(lih-loh)//d10),int(-0.5*(lih-loh)//d10),-1)): 
        #right now this is going down the diagonal, figure out how to do over the whole field
        sp = np.nanmean(h_3_35[:,boh+int(inc*d10):bih+int((inc+1)*d10),
                               loh+int(inc*d10):lih+int((inc+1)*d10)],
                        axis=(1,2))
        sp2 = np.nanmean(h_3_35[:,boh+int(dec*d10):bih+int((dec+1)*d10),
                               loh+int(inc*d10):lih+int((inc+1)*d10)],
                        axis=(1,2))
#     cspec = np.nanmean(d_3_35[:,boc:bic,loc:lic],axis=(1,2))

#         axs[2].plot(sp)
#         axs[2].plot(sp2)
        axs[1].scatter(np.nanmean([int(inc*d10),(lih-loh)+int((inc+1)*d10)]),
                       np.nanmean([int(inc*d10),(bih-boh)+int((inc+1)*d10)]),
                       color='white',alpha=0.6
                       )
        axs[1].scatter(np.nanmean([int(inc*d10),(lih-loh)+int((inc+1)*d10)]),
                       np.nanmean([int(dec*d10),(bih-boh)+int((dec+1)*d10)]),
                       color='white',alpha=0.6
                       )
#     axs[2].plot(cspec/np.nanmax(cspec))

    
    plt.show()


for i in range(17):    
    loh,lih = zcat['l_pix_HI'][i]
    loc,lic = zcat['l_pix_GMC'][i]
    
    boh,bih = zcat['b_pix_HI'][i]
    boc,bic = zcat['b_pix_GMC'][i]
    
    X,Y = np.array(np.meshgrid(np.arange(loh,lih),np.arange(boh,bih)))
    xo,yo = np.nanmean([loh,lih]),np.nanmean([boh,bih])
    Z = np.sqrt((X-xo)**2 + (Y-yo)**2)*0.125*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    rs = np.ravel(Z)
    
    
    fig,axs = plt.subplots(1,3,figsize=[15,5])
    
    axs[0].contourf(h1[boh:bih,loh:lih])
    axs[0].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    axs[1].contourf(h2[boh:bih,loh:lih])
    axs[1].contour(c0[boh:bih,loh:lih],levels=[5e20])
    
    hspec = np.nanmean(h_3_35[:,boh:bih,loh:lih],axis=(1,2))
    cspec = np.nanmean(d_3_35[:,boc:bic,loc:lic],axis=(1,2))

    axs[2].plot(hspec/np.nanmax(hspec))
    axs[2].plot(cspec/np.nanmax(cspec))

    
    plt.show()


plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(h2)/bav,extent=[360,-360,-50,50],
             levels=6,norm='asinh')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],linewidths=0.5,
            levels=[5e20],colors='white')
plt.xlim(360,0)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(np.log10(c0),levels=10,vmin=np.log10(5e21),
            extent=[360,-360,-50,50])
plt.colorbar()

plt.xlim(180,-180)


with_co = h2 * (c0>5e20)
wo_co = h2 - with_co

with_denseco = h2 * (c0>5e21)


plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(wo_co),extent=[360,-360,-50,50],
             levels=6,norm='asinh')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],linewidths=0.5,
            levels=[5e20],colors='white')
plt.xlim(360,0)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(with_co),extent=[360,-360,-50,50],
             levels=6,norm='asinh')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],linewidths=0.5,
            levels=[5e20],colors='white')
plt.xlim(360,0)
# plt.ylim(-30,30)


plt.hist(np.ravel(np.sqrt(with_co)),bins=50,alpha=0.5,density=True)
plt.hist(np.ravel(np.sqrt(wo_co)),bins=50,alpha=0.5,density=True)
plt.hist(np.ravel(np.sqrt(with_denseco)),bins=50,alpha=0.5,density=True)

plt.xlim(0.1,)
plt.yscale('log')
plt.show()


plt.hist(np.ravel(np.sqrt(with_co)-bav),bins=50,alpha=0.5)
plt.hist(np.ravel(np.sqrt(wo_co)-bav),bins=50,alpha=0.5)
plt.yscale('log')
plt.show()


#start with a brief elevator pitch
#pick 2 plots to illustrate the questions i have

#add to calendar
#nearby universe group meeting - 2:30
#share a few slides
#explicitly say questions im working on 


#11 am friday - galaxies meeting
#propose 12pm for meeting with blakeseley

#both on 5th floor
#12pm feticia, 3d dust, 11th floor
