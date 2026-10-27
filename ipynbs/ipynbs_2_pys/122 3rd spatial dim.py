#!/usr/bin/env python
# coding: utf-8

# notes for paper:
# in methods, mention that i used the moment masked version, then interpolated it

# the difference caused by degrading the resolution:
# (all in percent)
# mass difference [1.1,0.9,2.2,0.3,2.5,2.0,2.3,0.5,0.1,2.4,0.5,2.0,0.4,0.2,0.0,1.7,0.0,0.8,2.2,0.9,2.8,3.1]
# area difference [0.8,0.1,0.6,0.2,1.6,1.8,0.8,1.5,0.2,0.8,0.0,0.7,1.3,0.2,0.2,2.0,0.8,0.2,0.0,0.9,0.8,0.6]
# 
# 0th moment maps & hists are v similar (sometimes the og has slightly more counts at higher n_h bins, but looks p small
# 
# first & second moment hists are also pretty on top of each other (although for second mom there are some where the bins at smaller linewidths have slightly more counts for the og and the corresponding resampled has slightly higher counts at larger linewidhts, looks p small though
# 
# ^ all of this to be like, the "real" distributions of 0th, 1st, 2nd moment per pixel is preserved in the resampled version
# 
# mention that cam, rcra, aqs have velocity dispersions that is less than the channel widths
# not sure what this introduces in terms of error but yeah
# 
# mention that gem,mad,cma,ros are in the galactic plane and also the furthest away
# resolution stuff, blending stuff
# 
# the velocity channel edges were chosen based on where the first/second derivative of the spectra were at zero
# (second derivative for single peaks / the side of the spectrum where it just goes down, first derivative to find
# the saddle point between multiple peaks)


from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import numpy as np
import pandas as pd
import datacube as dc
import fun as fns

import matplotlib.pyplot as plt


hi_og = fits.open('../not_rn_2/HI4PI_CAR.fits')[0].data


fits.open('../not_rn_2/HI4PI_CAR.fits')[0].header


lns_og = (np.arange(4323)-2162)*-0.083
lts_og = (np.arange(2163)-1082)*0.083
vls_og = (np.arange(945)-472)*1.288


hi4pi,dht = fits.open('both_same_grid_720d.fits')


d0 = np.nansum(dht.data,axis=0)*1.3*1.97e20


cat = QTable(fits.open('post_flash_cat.fits')[1].data)


co_og_mom = fits.open('../not_rn_2/COGAL_all_mom.fits')[0]


co_og_mom.header


hi_og_100 = hi_og[394:549,:,:]


def mass_n_area_og(cube,d,cloud_type):
    w = np.nansum(cube*1.3*u.km/u.s,axis=0)
    
    if cloud_type == '21cm':
        conv = 1.82e18*((u.cm)**-2)/(u.K*(u.km/u.s))
        n = len(np.where(w*conv>4e20/u.cm**2)[0])
        
    elif cloud_type == 'CO':
        conv = 1.97e20*((u.cm)**-2)/(u.K*(u.km/u.s))*2
        n = len(np.where((w*conv)>1e21/u.cm**2)[0]) #for calculating total area
        #but shouldnt be put into dA when using that for mass, because double
        #counting the area

#     print(conv.unit)
    sig = conv*w*1.36*c.m_p
    
    dA = (d*0.083*u.deg.to(u.rad))**2
        
    m = np.nansum(sig)*dA
    
    return m.to(u.kg),(dA*n).to(u.pc**2)


#checking differences between og hi & resampled
for i in range(22):
    logh,ligh = ((cat['lon_GMC<deg>'][i]/-0.083)+2162).astype('int')
    bogh,bigh = ((cat['lat_GMC<deg>'][i]/0.083)+1082).astype('int')
    
    sb = hi_og_100[:,bogh:bigh,logh:ligh]
    sb0 = np.nansum(sb,axis=0)*1.3*1.82e18
    
    sb_rcnt = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    sb0_rcnt = np.nansum(sb_rcnt,axis=0)*1.3*1.82e18
    
#     fig,axs = plt.subplots(1,2)
#     axs[0].contourf(np.log10(sb0),levels=np.arange(20,22.5,0.1),cmap='Paired')
#     axs[1].contourf(np.log10(sb0_rcnt),levels=np.arange(20,22.5,0.1),cmap='Paired')
#     plt.title(cat['name'][i])
#     plt.show()

#     plt.hist(np.ravel(np.log10(sb0)),alpha=0.5,bins=np.arange(20,22.1,0.05),density=True)
#     plt.hist(np.ravel(np.log10(sb0_rcnt)),alpha=0.5,bins=np.arange(20,22.1,0.05),density=True)
#     plt.show()

#     mog,aog = mass_n_area_og(sb*u.K,cat['distance<pc>'][i]*u.pc,'21cm')
#     mnw,anw = fns.mass_n_area(sb_rcnt*u.K,cat['distance<pc>'][i]*u.pc,'21cm')
    
#     print(np.round(np.abs(100-np.round(100*mog.to(u.solMass)/mnw.to(u.solMass),2)),1))
#     print(np.round(np.abs(100-np.round(100*aog/anw,2)),1))


#     ogmom1 = fns.get_mom1(sb,np.arange(155)*1.3-100)
#     nwmom1 = fns.get_mom1(sb_rcnt,np.arange(155)*1.3-100)
    
    
#     plt.hist(np.ravel(ogmom1),alpha=0.5,bins=np.arange(-50,50,0.65/2),density=True)
#     plt.hist(np.ravel(nwmom1),alpha=0.5,bins=np.arange(-50,50,0.65/2),density=True)
#     plt.show()
    
    
#     ogmom2 = fns.get_momN(sb,np.arange(155)*1.3-100,2)
#     nwmom2 = fns.get_momN(sb_rcnt,np.arange(155)*1.3-100,2)
    
    
#     plt.hist(np.ravel(np.sqrt(ogmom2)),alpha=0.5,bins=np.arange(0,40,0.65/2),density=True)
#     plt.hist(np.ravel(np.sqrt(nwmom2)),alpha=0.5,bins=np.arange(0,40,0.65/2),density=True)
#     plt.show()


#checking differences between integrating over all vs limited velocity channels
for i in range(22):
    long = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    short = hi4pi.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    co = dht.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    lmom1 = fns.get_mom1(long,np.arange(155)*1.3-100)
    smom1 = fns.get_mom1(short,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    cmom1 = fns.get_mom1(co,np.arange(155)*1.3-100)
    
    fig,axs = plt.subplots(1,4,figsize=(30,8))
    l = axs[0].contourf(lmom1,levels=10,cmap='Paired')
    fig.colorbar(l,ax=axs[0])
    
    s=axs[1].contourf(smom1,levels=10,cmap='Spectral')
    fig.colorbar(s,ax=axs[1])
    fig.suptitle(cat['name'][i])
    
    
    axs[2].hist(np.ravel(lmom1),alpha=0.5)
    axs[2].hist(np.ravel(smom1),alpha=0.5)
    axs[2].hist(np.ravel(cmom1),alpha=0.5)
    
    
    ysz, xsz = cat['bic'][i]-cat['boc'][i],cat['lic'][i]-cat['loc'][i]
    
    
    [(
      axs[3].plot(np.arange(155)*1.3-100,
              np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))),
              color='black',alpha=0.4,linewidth=0.3),

      axs[3].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
              np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))),
              color='green',alpha=0.4,linewidth=0.3),

      axs[3].plot(np.arange(155)*1.3-100,
              np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))),
              color='red',alpha=0.4,linewidth=0.3),
    #   plt.show()


     ) for y in range(0,24,4) for x in range(0,32,4)]
#     plt.title(cat['name'][i])
#     plt.show()
    
    plt.show()


for i in range(22):
    long = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    short = hi4pi.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    co = dht.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    ysz, xsz = cat['bic'][i]-cat['boc'][i],cat['lic'][i]-cat['loc'][i]
    
    
    [(
      plt.plot(np.arange(155)*1.3-100,
              np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))),
              color='black',alpha=0.4,linewidth=0.3),

      plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
              np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))),
              color='green',alpha=0.4,linewidth=0.3),

      plt.plot(np.arange(155)*1.3-100,
              np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))),
              color='red',alpha=0.4,linewidth=0.3),
    #   plt.show()


     ) for y in range(0,24,4) for x in range(0,32,4)]
    plt.title(cat['name'][i])
    plt.show()
    

for i in range(22):
    long = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    short = hi4pi.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    co = dht.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    ysz, xsz = cat['bic'][i]-cat['boc'][i],cat['lic'][i]-cat['loc'][i]
    lmom1 = fns.get_mom1(long,np.arange(155)*1.3-100)
    smom1 = fns.get_mom1(short,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    [(
      plt.plot(np.arange(155)*1.3-100,
              np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))),
              color='black',alpha=0.4,linewidth=0.3),
        plt.vlines(np.nanmean(lmom1[y:y+4,x:x+4]),0,1,color='orange',
                  alpha=0.4,linewidth=0.3),
        

      plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
              np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))),
              color='green',alpha=0.4,linewidth=0.3),
        
        plt.vlines(np.nanmean(smom1[y:y+4,x:x+4]),0,1,color='blue',
                  alpha=0.4,linewidth=0.3),
        

#       plt.plot(np.arange(155)*1.3-100,
#               np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))/\
#                np.nanmax(np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))),
#               color='red',alpha=0.4,linewidth=0.3),
    #   plt.show()


     ) for y in range(0,24,4) for x in range(0,32,4)]
    plt.title(cat['name'][i])
    plt.show()
    

for i in range(22):
    long = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    short = hi4pi.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    co = dht.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    ysz, xsz = cat['bic'][i]-cat['boc'][i],cat['lic'][i]-cat['loc'][i]
    cmom1 = fns.get_mom1(co,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    [(
#       plt.plot(np.arange(155)*1.3-100,
#               np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))/\
#                np.nanmax(np.nanmean(long[:,y:y+4,x:x+4],axis=(1,2))),
#               color='black',alpha=0.4,linewidth=0.3),
#         plt.vlines(np.nanmean(lmom1[y:y+4,x:x+4]),0,1,color='orange',
#                   alpha=0.4,linewidth=0.3),
        

#       plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
#               np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))/\
#                np.nanmax(np.nanmean(short[:,y:y+4,x:x+4],axis=(1,2))),
#               color='green',alpha=0.4,linewidth=0.3),
        
#         plt.vlines(np.nanmean(smom1[y:y+4,x:x+4]),0,1,color='blue',
#                   alpha=0.4,linewidth=0.3),
        

      plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
              np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))/\
               np.nanmax(np.nanmean(co[:,y:y+4,x:x+4],axis=(1,2))),
              color='red',alpha=0.4,linewidth=0.3),
#       plt.show()


     ) for y in range(0,24,4) for x in range(0,32,4)]
    
    plt.hist(np.ravel(cmom1),density=True,bins=100)
    plt.yscale('log')
    plt.title(cat['name'][i])
    plt.show()
    

for i in range(22):
    long = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    short = hi4pi.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    co = dht.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    ysz, xsz = cat['bic'][i]-cat['boc'][i],cat['lic'][i]-cat['loc'][i]
    cmom1 = fns.get_mom1(co,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    plt.plot(np.arange(155)*1.3-100,np.nansum(co,axis=(1,2)))
    plt.plot(np.arange(155)*1.3-100,100*np.nanmean(long,axis=(1,2)))
    
    plt.title(cat['name'][i])
#     plt.yscale('log')
    plt.hlines([0.2*155,0.2*3*155],-100,100)
    
    chanel = np.around((np.arange(155)*1.3-100)[np.nansum(co,axis=(1,2))>0.2*3*155],2)
    
    plt.vlines([chanel[0],chanel[-1]],0,7e2,color='red')
    
    chanel = np.around((np.arange(155)*1.3-100)[np.nansum(co,axis=(1,2))>0.2*155],2)
    plt.vlines([chanel[0],chanel[-1]],0,7e2,color='orange')
    
    plt.xlim(chanel[0]-5,chanel[-1]+5)
    
    plt.show()


for i in range(22):
    long = hi4pi.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    short = hi4pi.data[cat['vo'][i]:cat['vi'][i],cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    
    co = dht.data[:,cat['boc'][i]:cat['bic'][i],cat['loc'][i]:cat['lic'][i]]
    ysz, xsz = cat['bic'][i]-cat['boc'][i],cat['lic'][i]-cat['loc'][i]
    cmom1 = fns.get_mom1(co,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    chanel = np.around((np.arange(155)*1.3-100)[np.nansum(co,axis=(1,2))>0.2*3*155],2)
    plt.vlines([chanel[0],chanel[-1]],0,1,color='red')
    plt.annotate(str([chanel[0],chanel[-1]]),(-40,0.9))
    
    chanel = np.around((np.arange(155)*1.3-100)[np.nansum(co,axis=(1,2))>0.2*155],2)
    plt.vlines([chanel[0],chanel[-1]],0,1,color='green')
    
    
    plt.vlines(cat['vel<km/s>'][i],0,1,color='black')
    
    plt.annotate(str([chanel[0],chanel[-1]]),(-40,0.8))
    plt.annotate(str(np.arange(chanel[0],chanel[-1],1.3).shape),(-20,0.8))
    plt.annotate(str(cat['vel<km/s>'][i]),(-40,0.7))
    
    
    plt.plot(np.arange(155)*1.3-100,np.nanmean(co,axis=(1,2))/np.nanmax(np.nanmean(co,axis=(1,2))))
    plt.plot(np.arange(155)*1.3-100,np.nanmean(long,axis=(1,2))/np.nanmax(np.nanmean(long,axis=(1,2))))
    plt.title(cat['name'][i])
    plt.xlim(-50,50)
    plt.show()


#even more conservative velocity channels
convs = [
    [0,6],
    [-12,8],
    [-3,2],
    [-2,8],
    [9,20],
    [2,8],
    [8,17],
    [-6,20],
    [1,16],
    [-4,6],
    [-10,14],
    [20,35],
    [-2,12],
    [5,17],
    [-3,8],
    [0,16],
    [-2,16],
    [-6,16],
    [-4,14],
    [-9,6],
    [7,20],
    [-2,15]
]


nvo,nvi = ((np.array(convs).transpose()+100)/1.3).astype('int')


dm_hi = np.zeros(22)
dm_co = np.zeros(22)
da_hi = np.zeros(22)
da_co = np.zeros(22)


mcnew = np.zeros(22)


mcnew*u.kg.to(u.solMass)/cat['mass_GMC<Msun>']


def ra_from_gmmc(m):
    Sigma = 10*u.Msun/u.pc**2
    ra = (np.sqrt(m/(np.pi*Sigma))).to(u.pc)
    return ra


ra_from_gmmc((mcnew*u.kg).to(u.solMass))


ra_from_gmmc((mcold*u.kg).to(u.solMass))/(np.array(cat['radius_HI<pc>'])*u.pc)


ra_from_gmmc((mcnew*u.kg).to(u.solMass))/(np.array(cat['radius_HI<pc>'])*u.pc)


np.array(cat['radius_HI<pc>'])*u.pc


mcnew*u.kg.to(u.solMass)/cat['mass_GMC<Msun>']


mcold*u.kg.to(u.solMass)/cat['mass_GMC<Msun>']


mhold*u.kg.to(u.solMass)


mhnew*u.kg.to(u.solMass)


mhold*u.kg.to(u.solMass)/cat['mass_HI<Msun>']


np.array((cat['name'],mhnew*u.kg.to(u.solMass)/cat['mass_HI<Msun>'])).transpose()


mcold = np.zeros(22)


mhnew = np.zeros(22)
mhold = np.zeros(22)


for i in range(22):
    hlib = hi4pi.data[cat['vo'][i]:cat['vi'][i],
                   cat['boh'][i]:cat['bih'][i],
                   cat['loh'][i]:cat['lih'][i]]
    
    d3 = dht.data*(dht.data>0.2*3)
    clib = d3[cat['vo'][i]:cat['vi'][i],
                   cat['boc'][i]:cat['bic'][i],
                   cat['loc'][i]:cat['lic'][i]]
    
    
    hcsv = hi4pi.data[nvo[i]:nvi[i],
                   cat['boh'][i]:cat['bih'][i],
                   cat['loh'][i]:cat['lih'][i]]
    
    ccsv = d3[nvo[i]:nvi[i],
                   cat['boc'][i]:cat['bic'][i],
                   cat['loc'][i]:cat['lic'][i]]
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,np.nanmean(hlib,axis=(1,2)))
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,100*np.nanmean(clib,axis=(1,2)))
    
    plt.plot(np.arange(nvo[i],nvi[i])*1.3-100,100*np.nanmean(ccsv,axis=(1,2)))
    plt.plot(np.arange(nvo[i],nvi[i])*1.3-100,np.nanmean(hcsv,axis=(1,2)))

    plt.title(cat['name'][i])
    plt.show()
    

def dv_2nd(array):
    dv1 = np.gradient(array)
    return np.gradient(dv1)


for i in [1]:
    hlng = hi4pi.data[:,
                   cat['boh'][i]:cat['bih'][i],
                   cat['loh'][i]:cat['lih'][i]]
    
    d3 = dht.data*(dht.data>0.2*3)
    clng = d3[:,
                   cat['boc'][i]:cat['bic'][i],
                   cat['loc'][i]:cat['lic'][i]]
    
    plt.plot(np.arange(155)*1.3-100,np.nanmean(hlng,axis=(1,2)))
    plt.plot(np.arange(155)*1.3-100,100*np.nanmean(clng,axis=(1,2)))

    plt.plot(np.arange(155)*1.3-100,10*np.gradient(np.nanmean(clng,axis=(1,2))),linestyle='--',color='red')
    plt.plot(np.arange(155)*1.3-100,0.25*np.gradient(np.nanmean(hlng,axis=(1,2))),linestyle='--',color='black')
    
#     plt.plot(np.arange(155)*1.3-100,10*dv_2nd(np.nanmean(clng,axis=(1,2))),color='red')
#     plt.plot(np.arange(155)*1.3-100,0.25*dv_2nd(np.nanmean(hlng,axis=(1,2))),color='black')
    
    
    plt.title(cat['name'][i])
    plt.xlim((nvo[i]-30)*1.3-100,(nvi[i]+30)*1.3-100)
    
    plt.vlines(sighh_final[i],-2,20,color='blue')
#     plt.vlines(cat['vel<km/s>'][i],-2,30,color='red')
#     plt.vlines(dv1d1s[i],-0,20,color='orange')
    
#     plt.ylim(-2,2)
    plt.show()
    

sighh_final = [[-12,10],[-40,20],[-20,20],[-13,20],[6,28],[-15,25],
[-20,30],[-30,40],[-20,50],[-22,30],[-10,25],[20,38],
[-5,15],[-15,30],[-10,20],[-9,23],[-20,30],[-20,25],[-15,25],
[-10,10],[-13,21],[-10,20]]


cat['vel<km/s>'] = sighh_final


cat.write('post_flash_cat.fits',overwrite=True)


#orange
dv1d1s = np.array([[-12,13],[-20,21],[-12,12],[-13,18],[4,32],[-12,23],[-16,25],[-22,35],[-12,45],[-18,20],
[-14,25],[18,40],[-7,18],[-15,32],[-12,20],[-9,23],[-20,30],[-20,25],[-15,25],[-13,10],
[-13,20],[-12,25]])


#blue
dvd0s = np.array([[-10,10],[-20,20],[-20,20],[-15,20],[5,25],[-15,25],[-20,30],[-30,40],[-20,50],[-25,30],
[-10,40],[20,40],[-5,15],[-20,30],[-10,20],[-10,40],[-20,40],[-20,40],[-20,30],[-10,25],
[-10,21],[-10,27]])


np.array(cat['vel<km/s>'])


#checking differences between consv vels & what ive been using
for i in range(22):
    hlib = hi4pi.data[cat['vo'][i]:cat['vi'][i],
                   cat['boh'][i]:cat['bih'][i],
                   cat['loh'][i]:cat['lih'][i]]
    
    d3 = dht.data*(dht.data>0.2*3)
    clib = d3[cat['vo'][i]:cat['vi'][i],
                   cat['boc'][i]:cat['bic'][i],
                   cat['loc'][i]:cat['lic'][i]]
    
    
    hcsv = hi4pi.data[nvo[i]:nvi[i],
                   cat['boh'][i]:cat['bih'][i],
                   cat['loh'][i]:cat['lih'][i]]
    
    ccsv = d3[nvo[i]:nvi[i],
                   cat['boc'][i]:cat['bic'][i],
                   cat['loc'][i]:cat['lic'][i]]
    
#     fig,axs = plt.subplots(1,2)
#     axs[0].contourf(np.log10(np.nansum(hlib,axis=0)*1.3*1.82e18),levels=np.arange(20,22.5,0.1),cmap='Paired')
#     im = axs[1].contourf(np.log10(np.nansum(hcsv,axis=0)*1.3*1.82e18),levels=np.arange(20,22.5,0.1),cmap='Paired')
#     plt.title(cat['name'][i])
#     fig.colorbar(im)
#     plt.show()

#     plt.hist(np.ravel(np.log10(np.nansum(hlib,axis=0)*1.3*1.82e18)),
#              bins=np.arange(19,22,0.05),alpha=0.5,
#              density=True,color='black')
    
#     plt.hist(np.ravel(np.log10(np.nansum(hcsv,axis=0)*1.3*1.82e18)),
#              bins=np.arange(19,22,0.05),alpha=0.5,
#              density=True,edgecolor='black',fill=False)
    
    
#     plt.hist(np.ravel(np.log10(np.nansum(clib,axis=0)*1.3*1.97e20)),
#              bins=np.arange(20,23,0.05),alpha=0.5,
#              density=True,color='red')
    
#     plt.hist(np.ravel(np.log10(np.nansum(ccsv,axis=0)*1.3*1.97e20)),
#              bins=np.arange(20,23,0.05),alpha=0.5,
#              density=True,edgecolor='red',fill=False)
    
    plt.title(cat['name'][i]+' '+str(nvi[i]-nvo[i]))
    
    mlb,alb = fns.mass_n_area(hlib*u.K,cat['distance<pc>'][i]*u.pc,'21cm')
    mcv,acv = fns.mass_n_area(hcsv*u.K,cat['distance<pc>'][i]*u.pc,'21cm')
    
    mhnew[i] = mcv.value
    mhold[i] = mlb.value
    
#     dm_hi[i] = (mlb-mcv)/mlb
#     da_hi[i] = (alb-acv)/alb
    
    
    mlb,alb = fns.mass_n_area(clib*u.K,cat['distance<pc>'][i]*u.pc,'CO')
    mcv,acv = fns.mass_n_area(ccsv*u.K,cat['distance<pc>'][i]*u.pc,'CO')
    
    dm_co[i] = (mlb-mcv)/mlb
    da_co[i] = (alb-acv)/alb
    mcnew[i] = mcv.value
    mcold[i] = mlb.value
#     print(mcv)
#     print(mlb)
    
#     hlb1 = fns.get_mom1(hlib,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
#     hcv1 = fns.get_mom1(hcsv,np.arange(nvo[i],nvi[i])*1.3-100)
    
#     clb1 = fns.get_mom1(clib,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
#     ccv1 = fns.get_mom1(ccsv,np.arange(nvo[i],nvi[i])*1.3-100)
    
#     clb0 = np.nansum(clib,axis=0)*1.3*1.97e20
    
#     plt.hist(np.ravel(hlb1),
#              bins=np.arange(-20,40,0.65/2),alpha=0.5,
#              density=True,color='black')
    
#     plt.hist(np.ravel(hcv1),
#              bins=np.arange(-20,40,0.65/2),alpha=0.5,
#              density=True,edgecolor='black',fill=False)
    
    
#     plt.hist(np.ravel(clb1),
#              bins=np.arange(-20,40,0.65/2),alpha=0.5,
#              density=True,color='red')
    
#     plt.hist(np.ravel(ccv1),
#              bins=np.arange(-20,40,0.65/2),alpha=0.5,
#              density=True,edgecolor='red',fill=False)

    
#     plt.xlim((nvo[i]-5)*1.3-100,(nvi[i]+5)*1.3-100)
    
    
#     hlb2 = fns.get_momN(hlib,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2)
#     hcv2 = fns.get_momN(hcsv,np.arange(nvo[i],nvi[i])*1.3-100,2)
    
#     clb2 = fns.get_momN(clib,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2)
#     ccv2 = fns.get_momN(ccsv,np.arange(nvo[i],nvi[i])*1.3-100,2)
    
    
#     plt.hist(np.ravel(np.sqrt(hlb2)),
#              bins=np.arange(0,20,0.65/2),
#              alpha=0.5,
#              density=True,color='black')
    
#     plt.hist(np.ravel(np.sqrt(hcv2)),
#              bins=np.arange(0,20,0.65/2),
#              alpha=0.5,
#              density=True,edgecolor='black',fill=False)
    
    
#     plt.hist(np.ravel(np.sqrt(clb2)),
#              bins=np.arange(0,20,0.65/2),
#              alpha=0.5,
#              density=True,color='red')
#     plt.show()
    
#     plt.hist(np.ravel(np.sqrt(ccv2)),
#              bins=np.arange(0,20,0.65/2),
#              alpha=0.5,
#              density=True,edgecolor='red',fill=False)
    
    
#     plt.contourf(clb2>5)
#     plt.colorbar()

#     print(cat['name'][i],round(np.nanmean(np.sqrt(clb2[clb0>3e20]))/cat['vel_dispersion_GMC<km/s>'][i],2))

#     plt.show()


for i in range(22):
    thish0 = np.nansum(hi4pi.data[nvo[i]:nvi[i],:,:],axis=0)*1.3*1.82e18
    
    plt.contourf(np.log10(thish0),levels=np.arange(20,22,0.25),extent=[360,-360,-50,50],cmap='Paired')
    plt.colorbar()
    plt.contour(np.log10(thish0),levels=[20.5],extent=[360,-360,-50,50],colors='black',linewidths=0.8)
    
    plt.title(cat['name'][i])
    plt.xlim(cat['lon_GMC<deg>'][i][0]+10,cat['lon_GMC<deg>'][i][0]-10)
    plt.ylim(cat['lat_GMC<deg>'][i][0]-10,cat['lat_GMC<deg>'][i][0]+10)
    
    plt.show()


np.log10(3e20)


cat


np.arange(-2.5,1.4,1.3).shape


0,2,5,  20,11,7,4


22-7


np.log10(cat['mass_GMC<Msun>'][[0,2,5]])


np.log10(np.nanmean(cat['mass_GMC<Msun>']))


np.log10(np.nanmean(cat['mass_GMC<Msun>'][[1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]]))


np.log10(np.nanmean(cat['mass_GMC<Msun>'][[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]]))


np.log10(cat['mass_GMC<Msun>'][[4,7,11,20]])


x  =cat['name']
y = cat['vel_grad_GMC<km/s/pc>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x  =cat['name']
y = cat['mass_GMC<Msun>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x  =cat['name']
y = cat['mass_HI<Msun>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['radius_GMC<pc>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['radius_HI<pc>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x  =cat['name']
y = cat['distance<pc>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['vel_grad_HI<km/s/pc>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['vel_dispersion_GMC<km/s>']*2.35
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.hlines([0.65,1.3,1.3*2,1.3*3],'california','corona_australis')
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['vel_dispersion_HI<km/s>']*2.35
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.hlines([0.65,1.3,1.3*2,1.3*3],'california','corona_australis')
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['rot_ax_GMC<deg>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
# plt.yscale('log')
plt.show()


x = cat['name']
y = cat['rot_ax_HI<deg>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
# plt.yscale('log')
plt.show()


x = cat['name']
y = cat['Sigma_HI<Msun/pc2>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['Sigma_GMC<Msun/pc2>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['a_vir_GMC']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['a_vir_HI']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['b_rot_GMC']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


x = cat['name']
y = cat['b_rot_HI']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


x = cat['name']
y = fcat['SFR<Msun/yr>']
plt.scatter(x[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]],
            y[[1,3,6,8,9,10,12,13,14,15,16,17,18,19,21]])

plt.scatter(x[[4,7,11,20]],
            y[[4,7,11,20]])
    
plt.scatter(x[[0,2,5]],
            y[[0,2,5]])
 
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


cat.colnames
