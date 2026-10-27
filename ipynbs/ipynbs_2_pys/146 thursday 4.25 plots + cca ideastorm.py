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


h0 = np.nansum(h_3_35,axis=0)*1.3*1.82e18
c0 = np.nansum(d_3_35,axis=0)*1.3*1.97e20


h1 = fns.get_mom1(h_3_35,np.arange(55)*1.3 - 35)
c1 = fns.get_mom1(d_3_35,np.arange(55)*1.3 - 35)


plt.hist(np.ravel(h1),bins=50)
plt.yscale('log')


h2 = fns.get_momN(h_3_35,np.arange(55)*1.3 - 35,2)
c2 = fns.get_momN(d_3_35,np.arange(55)*1.3 - 35,2)


np.log10([1e20,5e20,
                                                            1e21,5e21,
                                                            1e22,
                                                            ])


np.log10(5e20)


np.log10(4e20)


plt.figure(figsize=(15,7))
plt.contourf(np.log10(h0),extent=[360,-360,-50,50],levels=np.arange(20,22,0.2))

plt.colorbar(label=r'log ($N_{HI}$ [$cm^{-2}$])',orientation='horizontal')


plt.contour(np.log10(h0),extent=[360,-360,-50,50],linewidths=0.7,
            levels=np.log10([4e20]),colors='white')


plt.contourf(np.log10(c0),extent=[360,-360,-50,50],alpha=0.9,
            levels=np.log10([1e21,1e25]),colors='xkcd:pinkish')


plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitiude [deg]')

plt.xlim(180,-180)


plt.figure(figsize=(15,7))
plt.contourf(np.log10(c0),extent=[360,-360,-50,50],levels=np.arange(20.6,23,0.2))

plt.colorbar(label=r'log ($N_{H_2}$ [$cm^{-2}$])',orientation='horizontal')


plt.contour(np.log10(h0),extent=[360,-360,-50,50],linewidths=0.7,
            levels=np.log10([4e20]),colors='black')


# plt.contourf(np.log10(c0),extent=[360,-360,-50,50],alpha=0.9,
#             levels=np.log10([1e21,1e25]),colors='xkcd:pinkish')


plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitiude [deg]')

plt.xlim(180,-180)


make maps of the mach number -> possible for co assuming t=10k
assume that the co temp applies to the cnm
absorption survey - 21 cm sponge

h2 flourescence map

over the weekend -> make a patch by patch map of spectra comparing hi & co for at least one cloud


lombardi alves lada 2010 -> tables of sfr (correlating dense gas fraction with sfr)


narrow down a road map / plan for what to do for *this* paper -> scope scope scope


for i in range(17):
    fig,axs = plt.subplots(2,2,figsize=(15,10),layout='constrained')
    fig.suptitle(zcat['names'][i])
    
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    v0,v1 = zcat['vel_GMC<km/s>'][i]
    
    ho = axs[0,0].contourf(np.log10(h0[boh:bih,loh:lih]),levels=np.arange(20,22,0.2),
                    extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50])
    co = axs[0,0].contourf(np.log10(c0[boh:bih,loh:lih]),levels=15,
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  cmap='binary_r')
    
      
    hi = axs[0,1].contourf(h1[boh:bih,loh:lih],
                    extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                    levels=10,cmap='RdGy',)
    ci = axs[0,1].contourf(c1[boh:bih,loh:lih],
                   extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                   levels=np.arange(v0+1.5*1.3,v1-1.5*1.3,0.6),cmap='PiYG'
                   )
    axs[0,1].contour(np.log10(c0[boh:bih,loh:lih]),levels=[20.6,21],
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  colors='black',linewidths=0.6)
    
    dhidl, dhidb = np.gradient(h1[boh:bih,loh:lih])
    himag = np.sqrt(dhidl**2  + dhidb**2)
    hith = np.arctan(dhidb/dhidl)
    
    dcidl, dcidb = np.gradient(c1[boh:bih,loh:lih])
    cimag = np.sqrt(dcidl**2  + dcidb**2)
#     print(cimag.shape)
    cith = np.arctan(dcidb/dcidl)
    
    x,y = np.meshgrid(np.arange(loh,lih)*-0.125+360,np.arange(boh,bih)*0.125-50)
    hiq = axs[1,1].quiver(x,y,dhidl,dhidb,hith,pivot='mid',scale=0.5*(lih-loh),cmap='RdGy')
#     axs[1,1].quiver(x,y,dcidl,dcidb,cith,pivot='mid',scale=40)
    
    
    axs[1,1].contour(np.log10(c0[boh:bih,loh:lih]),levels=[20.6,21],
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  colors='black',linewidths=0.6)
#     plt.quiver(x,y,dvdx,dvdy,np.arctan(dvdy/dvdx)*u.rad.to(u.deg),angles='xy',scale=30)
# plt.colorbar()

# x,y=np.meshgrid(np.arange(104),np.arange(48))
    
    
    hii = axs[1,0].contourf(np.sqrt(h2[boh:bih,loh:lih]),
                    extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                            levels=5,norm='asinh')
    cii = axs[1,0].contourf(np.sqrt(c2[boh:bih,loh:lih]),
                   extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                            levels=np.arange(0,3,0.2),
                  linewidths=0.3,cmap='binary_r')
    
    axs[1,0].contour(np.log10(c0[boh:bih,loh:lih]),levels=[20.6,21],
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  colors='black',linewidths=0.6)
    
    
    axs[0,0].set_title('0th moment, integrated intensity')
    axs[0,1].set_title('1st moment, mean velocity')
    axs[1,0].set_title('sqrt (2nd moment), standard deviation/linewidth')
    
    for im, aa,lb in zip([ho,co,hi,ci,hii,cii,hiq],
                      [axs[0,0],axs[0,0],axs[0,1],axs[0,1],axs[1,0],axs[1,0],axs[1,1]],
                     [r'log ($N_{HI}$ [$cm^{-2}$])',r'log ($N_{H_2}$ [$cm^{-2}$])',
                      r'$v_{o,HI}$ [$km\ s^{-1}$]',r'$v_{o,H_2}$ [$km\ s^{-1}$]',
                      r'$\sigma_{v,HI}$ [$km\ s^{-1}$]',r'$\sigma_{o,H_2}$ [$km\ s^{-1}$]',
                      r'dv/ds [$km\ s^{-1}$]']):
        fig.colorbar(im, ax=aa,label=lb)
    
    for ax in axs.flat:
        ax.set_xlim(zcat['l_pix_HI'][i]*-0.125 + 360)
        ax.set_ylim(zcat['b_pix_HI'][i]*0.125 - 50)
        ax.set_ylabel('galactic latitude [deg]')
        ax.set_xlabel('galactic longitude [deg]')
#         ax.label_outer()
        
    plt.show()


for i in range(17):
#     fig,axs = plt.subplots(2,2,figsize=(15,10),layout='constrained')
#     fig.suptitle(zcat['names'][i])
    plt.figure(figsize=(10,10))
    
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    v0,v1 = zcat['vel_GMC<km/s>'][i]

    dhidl, dhidb = np.gradient(h1[boh:bih,loh:lih])
    himag = np.sqrt(dhidl**2  + dhidb**2)
    hith = np.arctan(dhidb/dhidl)
    
    dcidl, dcidb = np.gradient(c1[boh:bih,loh:lih])
    cimag = np.sqrt(dcidl**2  + dcidb**2)
    cith = np.arctan(dcidb/dcidl)
    
    x,y = np.meshgrid(np.arange(loh,lih)*-0.125+360,np.arange(boh,bih)*0.125-50)
    plt.quiver(x,y,dhidl,dhidb,hith,pivot='mid',scale=0.3*(lih-loh),cmap='RdGy')
    
    
    plt.contour(np.log10(c0[boh:bih,loh:lih]),levels=[20.6,21],
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  colors='black',linewidths=0.6)

    plt.show()


for i in range(17):
    plt.figure(figsize=(10,10))
    
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    dhidl, dhidb = np.gradient(h1[boh:bih,loh:lih])
    himag = np.sqrt(dhidl**2  + dhidb**2)
    hith = np.arctan(dhidb/dhidl)
    
    dhiidl, dhiidb = np.gradient(h2[boh:bih,loh:lih])
    hiimag = np.sqrt(dhiidl**2  + dhiidb**2)
    hiith = np.arctan(dhiidb/dhiidl)
    
    dciidl, dciidb = np.gradient(c2[boh:bih,loh:lih])
    cimag = np.sqrt(dcidl**2  + dcidb**2)
#     print(cimag.shape)
    ciith = np.arctan(dciidb/dciidl)
    
    x,y = np.meshgrid(np.arange(loh,lih)*-0.125+360,np.arange(boh,bih)*0.125-50)
    plt.quiver(x,y,-dhidl/himag,-dhidb/himag,himag,pivot='mid',cmap=dc.cmap_min0('binary'),
               scale=0.8*(lih-loh),alpha=0.9,norm='asinh')
#     axs[1,1].quiver(x,y,dcidl,dcidb,cith,pivot='mid',scale=40)
#     plt.colorbar()

    
    plt.contour(np.log10(h0[boh:bih,loh:lih]),levels=20,
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50])
    
#     plt.contour(h1[boh:bih,loh:lih],levels=np.arange(-20,20,1),linewidths=0.6,
#                    extent=[loh*-0.125+360,lih*-0.125+360,
#                             boh*0.125-50,bih*0.125-50])
#     plt.colorbar(orientation='horizontal')


    plt.show()


for i in range(17):
    plt.figure(figsize=(10,10))
    
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    dhidl, dhidb = np.gradient(h1[boh:bih,loh:lih])
    himag = np.sqrt(dhidl**2  + dhidb**2)
    hith = np.arctan(dhidb/dhidl)
    
    dcidl, dcidb = np.gradient(c1[boh:bih,loh:lih])
    cimag = np.sqrt(dcidl**2  + dcidb**2)
    cith = np.arctan(dcidb/dcidl)

    
    x,y = np.meshgrid(np.arange(loh,lih)*-0.125+360,np.arange(boh,bih)*0.125-50)
    plt.quiver(x,y,dhidl/himag,dhidb/himag,himag,pivot='mid',cmap=dc.cmap_min0('binary'),
               scale=0.8*(lih-loh),alpha=0.9,norm='asinh')
    
    x,y = np.meshgrid(np.arange(loh,lih)*-0.125+360,np.arange(boh,bih)*0.125-50)
    plt.quiver(x,y,dcidl/cimag,dcidb/cimag,cimag,pivot='mid',cmap=dc.cmap_min0('viridis'),
               scale=0.8*(lih-loh),alpha=0.9,norm='asinh')
#     axs[1,1].quiver(x,y,dcidl,dcidb,cith,pivot='mid',scale=40)
#     plt.colorbar()

    
#     plt.contour(np.log10(h0[boh:bih,loh:lih]),levels=20,
#                    extent=[loh*-0.125+360,lih*-0.125+360,
#                             boh*0.125-50,bih*0.125-50])
    
#     plt.contour(h1[boh:bih,loh:lih],levels=np.arange(-20,20,1),linewidths=0.6,
#                    extent=[loh*-0.125+360,lih*-0.125+360,
#                             boh*0.125-50,bih*0.125-50])
#     plt.colorbar(orientation='horizontal')


    plt.show()


for i in range(17):
    plt.figure(figsize=(10,10))
    
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    dhidl, dhidb = np.gradient(h1[boh:bih,loh:lih])
    himag = np.sqrt(dhidl**2  + dhidb**2)
    hith = np.arctan(dhidb/dhidl)
    
    dcidl, dcidb = np.gradient(c1[boh:bih,loh:lih])
    cimag = np.sqrt(dcidl**2  + dcidb**2)
    cith = np.arctan(dcidb/dcidl)

    dot = dhidl/himag*dcidl/cimag + dhidb/himag*dcidb/cimag
    plt.imshow(dot,origin='lower',extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50])
    plt.colorbar()

    
    plt.show()


dhidlgal, dhidbgal = np.gradient(h1)
dcidlgal, dcidbgal = np.gradient(c1)


himaggal = np.sqrt(dhidlgal**2  + dhidbgal**2)
cimaggal = np.sqrt(dcidlgal**2  + dcidbgal**2)


hithgal = np.arctan(dhidbgal/dhidlgal)
cithgal = np.arctan(dhidbgal/dhidlgal)


galdot = dhidlgal/himaggal*dcidlgal/cimaggal + dhidbgal/himaggal*dcidbgal/cimaggal


plt.hist(np.ravel(galdot),bins=100)
plt.yscale('log')
plt.show()


plt.figure(figsize=(20,10))

plt.imshow(galdot*(galdot==-1),origin='lower',extent=[360,-360,-50,50],cmap=dc.cmap_min0('viridis'))

#            cmap='RdYlGn')
# plt.colorbar()

# plt.xlim(180,150)
# plt.ylim(-20,0)


plt.show()


plt.figure(figsize=(20,10))

plt.imshow(galdot*(galdot>0.9),origin='lower',extent=[360,-360,-50,50],cmap=dc.cmap_min0('viridis'))
#            cmap='RdYlGn')
plt.colorbar()

plt.xlim(180,-180)
# plt.xlim(-150,-180)
# plt.ylim(-20,0)

plt.show()


x,y = np.meshgrid(np.arange(5760)*-0.125+360,np.arange(800)*0.125-50)


for i in range(17):
    fig,axs = plt.subplots(2,2,figsize=(15,10),layout='constrained')
    fig.suptitle(zcat['names'][i])
    
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    v0,v1 = zcat['vel_GMC<km/s>'][i]
    
    ho = axs[0,0].contourf(np.log10(h0[boh:bih,loh:lih]),levels=np.arange(20,22,0.2),
                    extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50])
    co = axs[0,0].contourf(np.log10(c0[boh:bih,loh:lih]),levels=15,
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  cmap='binary_r')
    
      
    hi = axs[0,1].contourf(h1[boh:bih,loh:lih],
                    extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                    levels=10,cmap='RdGy',)
    ci = axs[0,1].contourf(c1[boh:bih,loh:lih],
                   extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                   levels=np.arange(v0+1.5*1.3,v1-1.5*1.3,0.6),cmap='PiYG'
                   )
    axs[0,1].contour(np.log10(c0[boh:bih,loh:lih]),levels=[20.6,21],
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  colors='black',linewidths=0.6)
    
    
    hii = axs[1,0].contourf(np.sqrt(h2[boh:bih,loh:lih]),
                    extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                            levels=5,norm='asinh')
    cii = axs[1,0].contourf(np.sqrt(c2[boh:bih,loh:lih]),
                   extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                            levels=np.arange(0,3,0.2),
                  linewidths=0.3,cmap='binary_r')
    
    axs[1,0].contour(np.log10(c0[boh:bih,loh:lih]),levels=[20.6,21],
                   extent=[loh*-0.125+360,lih*-0.125+360,
                            boh*0.125-50,bih*0.125-50],
                  colors='black',linewidths=0.6)
    
    
    axs[0,0].set_title('0th moment, integrated intensity')
    axs[0,1].set_title('1st moment, mean velocity')
    axs[1,0].set_title('sqrt (2nd moment), standard deviation/linewidth')
    
    for im, aa,lb in zip([ho,co,hi,ci,hii,cii],
                      [axs[0,0],axs[0,0],axs[0,1],axs[0,1],axs[1,0],axs[1,0]],
                     [r'log ($N_{HI}$ [$cm^{-2}$])',r'log ($N_{H_2}$ [$cm^{-2}$])',
                      r'$v_{o,HI}$ [$km\ s^{-1}$]',r'$v_{o,H_2}$ [$km\ s^{-1}$]',
                      r'$\sigma_{v,HI}$ [$km\ s^{-1}$]',r'$\sigma_{o,H_2}$ [$km\ s^{-1}$]']):
        fig.colorbar(im, ax=aa,label=lb)
    
    for ax in axs.flat:
        ax.set_xlim(zcat['l_pix_HI'][i]*-0.125 + 360)
        ax.set_ylim(zcat['b_pix_HI'][i]*0.125 - 50)
        ax.set_ylabel('galactic latitude [deg]')
        ax.set_xlabel('galactic longitude [deg]')
#         ax.label_outer()
        
    plt.show()


bavh0 = np.broadcast_to(np.nanmean(h0,axis=1),(5760,800)).transpose()


plt.figure(figsize=(20,9))

plt.contourf(np.log10(bavh0),
            extent=[360,-360,-50,50],)
plt.colorbar(label=r'log ($N_{HI}$ [$cm^{-2}$])',orientation='horizontal')

plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitiude [deg]')


plt.figure(figsize=(20,9))
plt.contourf(np.log10(h0-bavh0),extent=[360,-360,-50,50],norm='asinh',
             levels=np.arange(20,22,0.2))

plt.colorbar(label=r'log ($N_{HI}$ [$cm^{-2}$])',orientation='horizontal')


plt.contourf(np.log10(c0),extent=[360,-360,-50,50],alpha=0.9,
            levels=np.log10([1e21,1e25]),colors='xkcd:pinkish')

plt.xlim(180,-180)

plt.xlabel('galactic longitude [deg]')
plt.ylabel('galactic latitiude [deg]')

# plt.ylim(-30,30)


plt.figure(figsize=(20,9))
plt.contourf(np.log10(h0),extent=[360,-360,-50,50],norm='asinh',
             levels=np.arange(20,22,0.2))
plt.xlim(180,-180)
# plt.ylim(-30,30)


hercules ysos? -- look in the literature 
send thavisha hercules structure (lb)


the HI around a gmc might be: dissociating H2 or infalling atomic gas
    
calculate virial parameters ^ tells you whether infalling (a<1) or dissociating and going away (a>1)

calculate gravitational potential energy for HI+H2

try a gaussian process to see correlation between co & hi -> then take it to dust
(basically looking at the correlation between spectra across the entire cloud)


if you can separate out the hi that is and isnt traced by the co kinematically, what you'll get
which hi is bound to the co and cloud, and which hi is flowing'


tdharmawardena@flatironinstitute.org


for i in [5]:
    fig,axs = plt.subplots(1,3,figsize=(13,5))

    hlo, hli = zcat['l_pix_GMC'][i]
    hbo, hbi = zcat['b_pix_GMC'][i]
    
#     v0, v1 = zcat['velocity<km/s>'][i]
    
#     h0
    
    s = np.nanstd((np.array([hlo,hli])*u.deg.to(u.rad))*zcat['distance<pc>'][i])

    axs[0].contour(np.nansum(h_3_35[:,hbo:hbi,hlo:hli],axis=(1)),linewidths=0.7,
                   extent=[-s,s,1,15],levels=30,norm='log')

    # axs[0].contour(np.nansum(d_3_35[:,hbo:hbi,hlo:hli],axis=(1)),
    #              levels=10,extent=[-15.85,15.85,1,15],linewidths=0.7,alpha=0.5,
    #              cmap=dc.cmap_min0('binary_r'))

    axs[1].contour(np.nansum(d_3_35[:,hbo:hbi,hlo:hli],axis=(1)),
                 levels=10,extent=[-15.85,15.85,1,15],linewidths=0.7,
                 cmap=dc.cmap_min0('viridis'))


    axs[2].contourf(d_ld[:,572:607],extent=[-15.85,15.85,60,1200],
                   levels=20,linewidths=0.3,
                 cmap=dc.cmap_min0('viridis'))

    # axs[2].contour(np.nansum(d_3_35[:,hbo:hbi,hlo:hli],axis=(1)),
    #              levels=5,extent=[-15.85,15.85,80,290],linewidths=0.9,alpha=0.5,
    #              cmap=dc.cmap_min0('binary_r'))

    for ax in axs.flat:
        ax.set_ylim(6,14)
    #     ax.set_xlim(49,41)

    # # axs[2].set_xlim(49,40)
    axs[2].set_ylim(210,225)


d_bd = np.nansum(dust,axis=2)


plt.contourf(d_bd,extent=[-90,90,60,1200],norm='log',
                   levels=np.arange(0,2,0.05),linewidths=0.3,
                 cmap=dc.cmap_min0('viridis'))
# plt.ylim(300,450)
# plt.xlim(-135,-55)
plt.yscale('log')


plt.contourf(d_ld,extent=[180,-180,60,1200],norm='log',
                   levels=np.arange(0,2,0.05),linewidths=0.3,
                 cmap=dc.cmap_min0('viridis'))
# plt.ylim(300,450)
# plt.xlim(-135,-55)
plt.yscale('log')


plt.contourf(d_ld,extent=[180,-180,60,1200],norm='log',
                   levels=np.arange(0,2,0.05),linewidths=0.3,
                 cmap=dc.cmap_min0('viridis'))
plt.ylim(150,300)
plt.xlim(55,35)
# plt.yscale('log')


generalized morphoogical component analysis -> 
separates out spatially overlpping but spectrally distinct components


picquenot et al 2019


s


s


(np.array([hlo,hli]))


*zcat['distance<pc>'][i]


i=5
hlo, hli = zcat['l_pix_GMC'][i]
hbo, hbi = zcat['b_pix_GMC'][i]

# v0, v1 = zcat['velocity<km/s>'][i]


s = np.nanstd((np.array([hlo,hli])*u.deg.to(u.rad))*zcat['distance<pc>'][i])

plt.contour(np.nansum(h_3_35[:,hbo:hbi,hlo:hli],axis=(1)),linewidths=0.7,
                   extent=[-s,s,1,15],levels=30,norm='log')

plt.ylim(6,14)


i=5
hlo, hli = zcat['l_pix_HI'][i]
hbo, hbi = zcat['b_pix_HI'][i]

# v0, v1 = zcat['velocity<km/s>'][i]


s = np.nanstd((np.array([hlo,hli])*u.deg.to(u.rad))*zcat['distance<pc>'][i])

plt.contour(np.nansum(h_3_35[:,hbo:hbi,hlo:hli],axis=(1)),linewidths=0.7,
                   extent=[-s,s,1,15],levels=30,norm='log')

plt.ylim(6,14)


(np.array([49,41])-180)/-0.229


np.nanstd((np.array([49,41])*u.deg.to(u.rad))*227)


np.nanmean((np.array([49,41])*u.deg.to(u.rad))*227)-


zcat[5]


# plt.figure(figsize=(16,10))
# plt.contour(h0,extent=[360,-360,-50,50],levels=[5e20])
# plt.contour(c0,extent=[360,-360,-50,50],levels=[5e20],colors='xkcd:reddish pink',linewidths=0.9)

plt.contour(np.log10(d_ld*1.87e21),
            extent=[180,-180,60,1200],levels=np.log10([5e19,1e20,2e20,3e20,4e20,5e20,1e21]),
            cmap='binary',linewidths=0.7)


plt.contourf(np.log10(d_ld*1.87e21),extent=[180,-180,60,1200],
             cmap=dc.cmap_min0('viridis'))
plt.colorbar(orientation='horizontal',label=r'log ($N_{H+H_2}$ [$cm^{-2}$])')

# plt.ylim(-50,50)
plt.xlim(49,41)
plt.ylim(150,300)
# plt.yscale('asinh')

plt.ylabel('heliocentric distance [pc]')
plt.xlabel('galactic longitude [deg]')


plt.plot(np.arange(55)*1.3-35,np.nanmean(h_3_35[:,hbo:hbi,hlo:hli],axis=(1,2)))
plt.plot(np.arange(55)*1.3-35,500*np.nanmean(d_3_35[:,hbo:hbi,hlo:hli],axis=(1,2)))


zcat


zcat[1]


lneg60 = int((-60-360)/-0.125)
bneg7 = int((-7+50)/0.125)


lneg65 = int((-65-360)/-0.125)
bneg20 = int((-20+50)/0.125)


lneg70 = int((-70-360)/-0.125)
bneg25 = int((-25+50)/0.125)


np.nanmax(h_3_35[:,bneg25,lneg70])


h1[bneg25,lneg70],np.sqrt(h2[bneg25,lneg70])


np.nanmax(h_3_35[:,bneg20,lneg65])


h1[bneg20,lneg65],np.sqrt(h2[bneg20,lneg65])


np.nanmax(h_3_35[:,bneg7,lneg60])


h1[bneg7,lneg60],np.sqrt(h2[bneg7,lneg60])


fig,axs = plt.subplots(1,3,figsize=(15,5))

axs[0].plot(np.arange(55)*1.3-35,h_3_35[:,bneg25,lneg70],label='real spectrum')
axs[0].plot(np.arange(55)*1.3-35,models.Gaussian1D(43.81,-4.23,7.22)(np.arange(55)*1.3-35),\
        label=r'gaussian spectrum with $v_o$ and $\sigma_v$')

# axs[0].legend(loc='lower left')

axs[0].set_xlabel('v [km/s]')
axs[0].set_ylabel(r'$T_B$ [K]')

axs[0].set_title('l,b = [-70,-25]')


axs[1].plot(np.arange(55)*1.3-35,h_3_35[:,bneg20,lneg65],label='real spectrum')
axs[1].plot(np.arange(55)*1.3-35,models.Gaussian1D(41.13,-0.38,8.45)(np.arange(55)*1.3-35),
        label=r'gaussian spectrum with $v_o$ and $\sigma_v$')

axs[1].legend(loc='lower left')

axs[1].set_xlabel('v [km/s]')
axs[1].set_ylabel(r'$T_B$ [K]')

axs[1].set_title('l,b = [-65,-20]')


axs[2].plot(np.arange(55)*1.3-35,h_3_35[:,bneg7,lneg60],label='real spectrum')
axs[2].plot(np.arange(55)*1.3-35,models.Gaussian1D(48.59,-1.3,16.43)(np.arange(55)*1.3-35),
        label=r'gaussian spectrum with $v_o$ and $\sigma_v$')

# axs[2].legend(loc='lower left')

axs[2].set_xlabel('v [km/s]')
axs[2].set_ylabel(r'$T_B$ [K]')

axs[2].set_title('l,b = [-60,-7]')


#^to these ones, first add in the CO spectra to make sure that theres ambiguity there too
#bc if CO "chooses" one peak then thats that


lneg57 = int((44-360)/-0.125)
bneg16 = int((9+50)/0.125)


plt.plot(np.arange(55)*1.3-35,h_3_35[:,bneg16,lneg57],label='real spectrum')
plt.plot(np.arange(55)*1.3-35,10*d_3_35[:,bneg16,lneg57],label='real spectrum')

# plt.plot(np.arange(55)*1.3-35,models.Gaussian1D(43.81,-4.23,7.22)(np.arange(55)*1.3-35),
#         label=r'gaussian spectrum with $v_o$ and $\sigma_v$')

# plt.legend(loc='lower left')

plt.xlabel('v [km/s]')
plt.ylabel(r'$T_B$ [K]')

plt.title('l,b = [-65,-20]')


plt.plot(np.arange(55)*1.3-35,h_3_35[:,bneg20,lneg65],label='real spectrum')
plt.plot(np.arange(55)*1.3-35,models.Gaussian1D(41.13,-0.38,8.45)(np.arange(55)*1.3-35),
        label=r'gaussian spectrum with $v_o$ and $\sigma_v$')

plt.legend(loc='lower left')

plt.xlabel('v [km/s]')
plt.ylabel(r'$T_B$ [K]')

plt.title('l,b = [-65,-20]')


plt.plot(np.arange(55)*1.3-35,h_3_35[:,bneg7,lneg60],label='real spectrum')
plt.plot(np.arange(55)*1.3-35,models.Gaussian1D(48.59,-1.3,16.43)(np.arange(55)*1.3-35),
        label=r'gaussian spectrum with $v_o$ and $\sigma_v$')

plt.legend(loc='lower left')

plt.xlabel('v [km/s]')
plt.ylabel(r'$T_B$ [K]')

plt.title('l,b = [-60,-7]')


#and make lv for third question, quiver plot to add to first question


dust = fits.open('mean_and_std_lbd.fits')[1].data


d0 = np.nansum(dust,0)


d_ld = np.nansum(dust,1)


plt.figure(figsize=(16,8))
plt.contour(h0,extent=[360,-360,-50,50],levels=[4e20])
plt.contour(c0,extent=[360,-360,-50,50],levels=[4e20],colors='xkcd:reddish pink',linewidths=0.7)
plt.contour(np.log10(d0*1.87e21),
            extent=[180,-180,-90,90],levels=[20.6],colors='white',linewidths=0.7)

plt.contourf(np.log10(d0*1.87e21),extent=[180,-180,-90,90],
             cmap=dc.cmap_min0('viridis'))
plt.colorbar(orientation='horizontal',label=r'log ($N_{H+H_2}$ [$cm^{-2}$])')

plt.ylim(-50,50)
plt.xlim(180,-180)

plt.ylabel('galactic latitude [deg]')
plt.xlabel('galactic longitude [deg]')


plt.figure(figsize=(16,10))
# plt.contour(h0,extent=[360,-360,-50,50],levels=[5e20])
# plt.contour(c0,extent=[360,-360,-50,50],levels=[5e20],colors='xkcd:reddish pink',linewidths=0.9)

plt.contour(np.log10(d_ld*1.87e21),
            extent=[180,-180,60,1200],levels=np.log10([5e19,1e20,2e20,3e20,4e20,5e20,1e21]),
            cmap='binary',linewidths=0.7)


plt.contourf(np.log10(d_ld*1.87e21),extent=[180,-180,60,1200],
             cmap=dc.cmap_min0('viridis'))
plt.colorbar(orientation='horizontal',label=r'log ($N_{H+H_2}$ [$cm^{-2}$])')

# plt.ylim(-50,50)
plt.xlim(180,-180)
plt.yscale('asinh')

plt.ylabel('heliocentric distance [pc]')
plt.xlabel('galactic longitude [deg]')
