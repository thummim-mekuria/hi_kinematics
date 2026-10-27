#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting

import fun as fns
import datacube as dc

hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data


fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


new = fcat[['name','lon_GMC<deg>','lat_GMC<deg>','vel<km/s>','distance<pc>']]


loc,lic = ((new['lon_GMC<deg>']-360)/-0.125).astype('int').transpose()
boc,bic = ((new['lat_GMC<deg>']+50)/0.125).astype('int').transpose()


vo,vi = ((new['vel<km/s>']+100)/1.3).transpose()


vo = vo.astype('int')
vi = np.around(vi).astype('int')


h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


vels = [np.arange(new['vel<km/s>'].transpose()[0][i],
          new['vel<km/s>'].transpose()[1][i],1.3) for i in range(22)]


masses = np.zeros(22)
radii = np.zeros(22)
for i in range(22):
    cube = cos[i]*u.K
        
    d = fcat['distance<pc>'][i]*u.pc
    
    m,a = fns.mass_n_area(cube,d,'CO')
    
    masses[i] = m.to(u.kg).value
    radii[i] = np.sqrt(a.to(u.pc**2)/np.pi).value
    

# fcat['mass_HI<kg>'] = masses
# fcat['radius_GMC<pc>'] = radii

# s, ra = fns.accumulation_radius(masses*u.kg,fcat['distance<pc>']*u.pc)
# sh = s.value
# fcat['radius_HI<pc>'] = ra
# fcat.write('re_expanded_cat.fits',overwrite=True)


new['mass_GMC<Msun>'] = masses*u.kg.to(u.solMass)
new['radius_GMC<pc>'] = radii*u.pc


new['radius_HI<pc>'] = fns.accumulation_radius(masses*u.kg,\
                        fcat['distance<pc>']*u.pc)[1]


plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


plt.figure(figsize=(10,10),dpi=200)
plt.hist(new['radius_HI<pc>'].value/new['radius_GMC<pc>'].value,
        bins = np.arange(1,3.5,0.25),color='xkcd:light blue',
        alpha=0.7,edgecolor='xkcd:blue')

plt.xlabel(r'$R_{HI}\ /\ R_{GMC}$')

plt.ylim(0,13)
plt.yticks([0,5,10],[0,5,10])

plt.ylabel(r'$N_{cloud}$')
plt.show()


lcenter = np.nanmean(fcat['lon_GMC<deg>'],axis=1)
bcenter = np.nanmean(fcat['lat_GMC<deg>'],axis=1)

sh = (fcat['radius_HI<pc>']/fcat['distance<pc>'])*u.rad.to(u.deg)


loh_old,lih_old = ((np.array([lcenter+sh,lcenter-sh]).transpose() \
            -360)/-0.125).astype('int').transpose()

boh_old,bih_old = ((np.array([bcenter-sh,bcenter+sh]).transpose()\
            +50)/0.125).astype('int').transpose()

his_old = [h3sig[vo[i]:vi[i],
                 boh_old[i]:bih_old[i],
                 loh_old[i]:lih_old[i]] for i in range(22)]


lcenter = np.nanmean(new['lon_GMC<deg>'],axis=1)
bcenter = np.nanmean(new['lat_GMC<deg>'],axis=1)

sh = (new['radius_HI<pc>'].value/new['distance<pc>'])*u.rad.to(u.deg)


loh,lih = np.around((np.array([lcenter+sh,lcenter-sh]).transpose() \
            -360)/-0.125).astype('int').transpose()

boh,bih = np.around((np.array([bcenter-sh,bcenter+sh]).transpose()\
            +50)/0.125).astype('int').transpose()

his_new = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]


for i in range(22):
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    
    plt.contourf(np.log10(np.nansum(his_old[i],axis=0)*1.3*1.82e18),
                 extent=[l0ho,l1ho,b0ho,b1ho],
                 levels = np.arange(19.5,22.5,0.1),
                cmap='YlGnBu_r')
    
    
    plt.contour(np.log10(np.nansum(his_new[i],axis=0)*1.3*1.82e18),
                 extent=[l0h,l1h,b0h,b1h],
                 levels = 10,
                cmap='magma',linewidths=0.5)
    
    
    plt.contourf(np.log10(np.nansum(his_new[i],axis=0)*1.3*1.82e18),
                 extent=[l0h,l1h,b0h,b1h],
                 levels = 10,
                cmap='magma',alpha=0.5)
    
#     plt.colorbar()
    
    
    plt.show()


h0 = np.nansum(h3sig,axis=0)*1.82e18*1.3


plt.rcParams['font.size']=35


plt.figure(figsize=(35,15),dpi=200)

plt.contourf(np.log10(h0),
                 extent=[360,-360,-50,50],
                 levels = np.arange(19.5,22.5,0.2),
                cmap='binary_r')
    
    
for i in range(22):
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
#     plt.contourf(np.log10(np.nansum(his_old[i],axis=0)*1.3*1.82e18),
#                  extent=[l0ho,l1ho,b0ho,b1ho],
#                  levels = np.arange(19.5,22.5,0.1),
#                 cmap='YlGnBu_r')
    
    
    plt.contourf(np.log10(np.nansum(his_new[i],axis=0)*1.3*1.82e18),
                 extent=[l0h,l1h,b0h,b1h],
                 levels = np.arange(19.5,22.5,0.2),
                cmap='YlGnBu',linewidths=0.5)
    
    
    plt.contourf(np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(19.5,22.5,0.2),
                cmap='magma',alpha=0.5)
    
#     plt.colorbar()
    plt.xlabel(r'galactic longitude [$\degree$]')
    plt.ylabel(r'galactic latitude [$\degree$]')
    
#     plt.


plt.xlim(180,-180)
    
    
# plt.savefig('new0th')
plt.show()


plt.rcParams['font.size']=20


#don't do this now -- cma & maddalena need to have their velocity ranges restrictred
#to just the hi that is co-inciding with the co


for i in range(22):

    fig,ax = plt.subplots(1,1,figsize=(10,10),dpi=200)
    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(h3sig[:,boh[i]:bih[i],loh[i]:lih[i]],
                       axis=(1,2)),color='black')

    ax2 = ax.twinx()

    ax2.plot(np.arange(155)*1.3-100,
            np.nanmean(d3sig[:,boc[i]:bic[i],loc[i]:lic[i]],
                       axis=(1,2)),color='red')

    ax.set_xlabel(r'$v_{lsr}\ [km\ s^{-1}]$')
    ax.set_ylabel(r'$T_{B,HI}\ [K$]')
    ax2.set_ylabel(r'$T_{B,H_2}\ [K$]')

    # ax2.set_ylim(0,0.5)
    # ax.set_ylim(0.80)

    plt.title(fcat['name'][i])

    plt.show()


for i in [4]:

    fig,ax = plt.subplots(1,1,figsize=(10,10),dpi=200)
#     ax.plot(np.arange(155)*1.3-100,
#             np.nanmean(h3sig[:,boh[i]:bih[i],loh[i]:lih[i]],
#                        axis=(1,2)),color='black')

# #     ax2 = ax.twinx()

    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(d3sig[:,boc[i]:bic[i],loc[i]:lic[i]],
                       axis=(1,2)),color='red')

    ax.set_xlabel(r'$v_{lsr}\ [km\ s^{-1}]$')
    ax.set_ylabel(r'$T_{B,H_2}\ [K$]')
    ax2.set_ylabel(r'$T_{B,H_2}\ [K$]')

    # ax2.set_ylim(0,0.5)
    # ax.set_ylim(0.80)

    plt.title(fcat['name'][i])
    plt.xlim(-20,60)

    plt.show()


for i in [0,2]:

    fig,ax = plt.subplots(1,1,figsize=(10,10),dpi=200)
    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(h3sig[:,boh[i]:bih[i],loh[i]:lih[i]],
                       axis=(1,2)),color='black')

    ax2 = ax.twinx()

    ax2.plot(np.arange(155)*1.3-100,
            np.nanmean(d3sig[:,boc[i]:bic[i],loc[i]:lic[i]],
                       axis=(1,2)),color='red')

    ax.set_xlabel(r'$v_{lsr}\ [km\ s^{-1}]$')
    ax.set_ylabel(r'$T_{B,HI}\ [K$]')
    ax2.set_ylabel(r'$T_{B,H_2}\ [K$]')

    # ax2.set_ylim(0,0.5)
    # ax.set_ylim(0.80)

    plt.title(fcat['name'][i])
    plt.xlim(-50,50)

    plt.show()


for i in range(22):

    fig,ax = plt.subplots(1,1,figsize=(10,10),dpi=200)
    ax.plot(np.arange(vo[i],vi[i])*1.3-100,
            np.nanmean(his_new[i],axis=(1,2)),color='black')

    ax2 = ax.twinx()

    ax2.plot(np.arange(vo[i],vi[i])*1.3-100,np.nanmean(cos[i],axis=(1,2)),color='red')

    ax.set_xlabel(r'$v_{lsr}\ [km\ s^{-1}]$')
    ax.set_ylabel(r'$T_{B,HI}\ [K$]')
    ax2.set_ylabel(r'$T_{B,H_2}\ [K$]')

    # ax2.set_ylim(0,0.5)
    # ax.set_ylim(0.80)

    plt.title(fcat['name'][i])

    plt.show()


for i in range(22):
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    
    wold = np.log10(np.nansum(his_old[i],axis=0)*1.3*1.82e18)
                 
    wnew = np.log10(np.nansum(his_new[i],axis=0)*1.3*1.82e18)
    
    wc = np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20)
    
    plt.hist(np.ravel(wold),bins = np.arange(19,23,0.01),alpha=0.7,density=True)
    plt.hist(np.ravel(wnew),bins = np.arange(19,23,0.01),alpha=0.7,density=True)
    plt.hist(np.ravel(wc),bins = np.arange(19,23,0.01),alpha=0.7,density=True)

#     plt.vlines(20.4,1e0,1e2,color='red',linestyle='-')
    
    plt.title(fcat['name'][i])
    plt.yscale('log')
    plt.show()


#from examining this, there are a good number of clouds where the column densities
#of hi are much higher (& with much smaller dispersion) than the co
#the dispersion part is true for all
#but the ones where the hi pdf is peaking when the co is mostly ending are:
#cma?, xbn, gem?, her(ish), lac?, madd, ros .. mostly the ones ~ the galactic plane


#weird dip feature at N ~ 10^(20.5) /cm^2
#looks to be tracing the edges of the clouds
#similar column density for all of them despite, 
#the insides having very different distributions of column density
#lends support to there being some sort of transition happening at N~4e20

for i in range(22):
    
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    
    wold = np.log10(np.nansum(his_old[i],axis=0)*1.3*1.82e18)
                 
    wnew = np.log10(np.nansum(his_new[i],axis=0)*1.3*1.82e18)
    
    wc = np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20)
    
#     plt.hist(np.ravel(wold),bins = np.arange(19,22,0.01),alpha=0.7)
#     plt.hist(np.ravel(wnew),bins = np.arange(19,22,0.01),alpha=0.7)
    plt.hist(np.ravel(wc),bins = np.arange(20,22,0.05),alpha=0.7,density=True)

    plt.vlines(np.log10(3e20),1e0,1e2,color='red',linestyle='-')
    
    plt.title(fcat['name'][i])
    plt.yscale('log')
    plt.show()
    
    
np.log10(3e20)


c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20


plt.figure(figsize=(20,10))

plt.imshow(np.log10(c0)*(np.log10(c0)>20.47),extent=[360,-360,-50,50],origin='lower')
plt.xlim(180,-180)


plt.hist(np.ravel(np.log10(c0)),bins = np.arange(20,22,0.05),alpha=0.7)
plt.vlines(np.log10(3e20),0,1e4)
plt.show()


len(wc[~np.isinf(wc)])


np.nansum(wc[~np.isinf(wc)]<20.5)


np.nansum(wc[~np.isinf(wc)]>20.5)


435+385


wc


np.nansum(wc>20.5)


(np.nansum(10**(less*wc[~np.isinf(wc)]))*c.m_p/u.cm**2).to(u.solMass/u.pc**2)


(np.nansum(10**(more*wc[~np.isinf(wc)]))*c.m_p/u.cm**2).to(u.solMass/u.pc**2)


980/(980+70011)


for i in range(22):
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    wc = np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20)

    tot = len(wc[~np.isinf(wc)])
    less = wc[~np.isinf(wc)]<20.4
    more = wc[~np.isinf(wc)]>20.4
    
    plt.imshow(wc*(wc<20.4),
                 extent=[l0c,l1c,b0c,b1c],origin='lower')
    
    
    print(np.nansum(less)/tot,np.nansum(more)/tot)
#                  levels = np.arange(19.5,22.5,0.1),
#                 cmap='binary_r',alpha=0.5)
    
#     plt.colorbar()
    
    
    plt.show()


for i in range(22):
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    plt.figure(figsize=(10,10))

    plt.contourf(np.log10(np.nansum(his_old[i],axis=0)*1.3*1.82e18),
                 extent=[l0ho,l1ho,b0ho,b1ho],
                 levels = np.arange(19.5,22.5,0.1),
                cmap='binary_r')
    
    plt.contourf(np.log10(np.nansum(his_new[i],axis=0)*1.3*1.82e18),
                 extent=[l0h,l1h,b0h,b1h],
                 levels = np.arange(20,22,0.1),
                cmap='YlGnBu_r',linewidths=0.5)
    
    
    plt.contour(np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(19.5,22.5,0.1),
                cmap='magma',alpha=0.5)
    
#     plt.colorbar()

    xo = np.nanmax([l0ho,l0h,l0c])
    xi = np.nanmin([l1ho,l1h,l1c])
    
    plt.xlim(xo,xi)
    plt.show()


masses_hi = np.zeros(22)
radii_hi = np.zeros(22)
for i in range(22):
    cube = his_new[i]*u.K
        
    d = fcat['distance<pc>'][i]*u.pc
    
    m,a = fns.mass_n_area(cube,d,'21cm')
    
    masses_hi[i] = m.to(u.kg).value
    radii_hi[i] = np.sqrt(a.to(u.pc**2)/np.pi).value
    

# fcat['mass_HI<kg>'] = masses
# fcat['radius_GMC<pc>'] = radii

# s, ra = fns.accumulation_radius(masses*u.kg,fcat['distance<pc>']*u.pc)
# sh = s.value
# fcat['radius_HI<pc>'] = ra
# fcat.write('re_expanded_cat.fits',overwrite=True)


np.array((fcat['name'],\
          masses*u.kg.to(u.solMass),\
          masses_hi*u.kg.to(u.solMass))).transpose()


plt.scatter(fcat['name'],masses*u.kg.to(u.solMass))

plt.xticks(rotation='vertical')

plt.hlines([10,100,1e2,1e3,1e4,1e5,1e6],'aquila_s','taurus')
plt.yscale('log')

plt.show()


plt.rcParams['font.size']=20


plt.hist(masses_hi/masses,color='xkcd:light blue',bins = np.arange(0,8,1),
        alpha=0.7,edgecolor='xkcd:blue')


plt.ylim(0,13)
plt.yticks([0,5,10],[0,5,10])

plt.xlabel(r'$M_{HI}\ /\ M_{GMC}$')

plt.ylabel(r'$N_{cloud}$')
plt.show()


new


plt.scatter(masses*u.kg.to(u.solMass),masses_hi/masses,color='black',s=20)
# plt.xticks(rotation='vertical')

plt.xlabel(r'$M_{GMC}\ [M_\odot]$')
plt.ylabel(r'$M_{HI}\ /\ M_{GMC}$')


plt.loglog()

plt.yticks([1,2,3,4,5,6],[1,2,3,4,5,6])
plt.show()


new['mass_HI<Msun>'] = masses_hi*u.kg.to(u.solMass)


fcat


for i in range(22):    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
    hmom1_old = fns.get_mom1(his_old[i],
                             np.arange(vo[i],vi[i])*1.3-100)
    hmom1_new = fns.get_mom1(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    plt.figure(figsize=(10,10))
    
    plt.contourf(hmom1_old,
                 extent=[l0ho,l1ho,b0ho,b1ho],
                 levels = 10,
                cmap='RdGy',linewidths=0.5)
    
    
    plt.contourf(hmom1_new,
                 extent=[l0h,l1h,b0h,b1h],
                 levels = 10,
                cmap='PiYG',linewidths=0.5)
    
    
#     plt.contour(np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20),)
#                  extent=[l0c,l1c,b0c,b1c],
#                  levels = np.arange(19.5,22.5,0.1),
#                 cmap='binary_r',alpha=0.5)
    
# #     plt.colorbar()
    
    
    plt.show()


fttr = fitting.LinearLSQFitter()
init = models.Linear1D()


OmHs = np.zeros(22)
OmCs = np.zeros(22)

ThHs = np.zeros(22)
ThCs = np.zeros(22)

for i in range(22):
    this_hi = his_new[i]
    this_co = cos[i]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,np.arange(vo[i],vi[i])*1.3-100)
    c1 = fns.get_mom1(this_co,np.arange(vo[i],vi[i])*1.3-100)
    
    
    h2 = fns.get_momN(this_hi,np.arange(vo[i],vi[i])*1.3-100,2)
    c2 = fns.get_momN(this_co,np.arange(vo[i],vi[i])*1.3-100,2)
    
    
    xc,yc = np.meshgrid(np.arange(loc[i],lic[i])*-0.125+360,
                        np.arange(boc[i],bic[i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(loh[i],lih[i])*-0.125+360,
                        np.arange(boh[i],bih[i])*0.125-50)
    
    
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0),
            np.ravel(h1),
            np.ravel(np.sqrt(h2))
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    hdf = hdf_noisy[(hdf_noisy['w']>3e20) &\
                    (hdf_noisy['s']>0) & \
                    ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] 
    #make sure that the linewidths are a valid number 
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0),
            np.ravel(c1),
            np.ravel(np.sqrt(c2))
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>4e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    #converting to pc now makes life wayyy easier
    hxpc = hdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]

    
    dvdxh = fttr(init,hxpc,hdf['v']) 
    dvdyh = fttr(init,hypc,hdf['v'])
    
    dvdxc = fttr(init,cxpc,cdf['v'])
    dvdyc = fttr(init,cypc,cdf['v'])
      
    
    omh = np.sqrt(dvdxh.slope**2 + dvdyh.slope**2)
    omc = np.sqrt(dvdxc.slope**2 + dvdyc.slope**2)

    thh = np.arctan(dvdyh.slope/dvdxh.slope)*u.rad.to(u.deg)
    thc = np.arctan(dvdyc.slope/dvdxc.slope)*u.rad.to(u.deg)
    
#     print(new['name'][i],dvdyc.slope,dvdxc.slope)
    #correcting for arctan only returning |theta|<90
    if dvdxh.slope > 0 : #first & fourth quadrants
        thh = thh
    elif dvdxh.slope < 0: #second quadrant & third
        thh = thh + 180
        
    if dvdxc.slope > 0 : #first & fourth quadrants
        thc = thc  
    elif dvdxc.slope < 0: #second quadrant & third
        thc = thc - 180

    OmHs[i] = omh
    OmCs[i] = omc
    
    ThHs[i] = thh
    ThCs[i] = thc
    

#     plt.scatter(cxpc,cdf['v'],s=1)
#     plt.plot(cxpc,cxpc*dvdxc.slope + dvdxc.intercept,color='red')
#     plt.show()
    
    
#     plt.imshow(h1,origin='lower',cmap='RdGy')
#     plt.show()
#     plt.imshow(c1,origin='lower',cmap='PiYG')
#     plt.show()
    
    
#     print(fcat['name'][i],dvdxh.slope.value,dvdyh.slope.value,
#           dvdxc.slope.value,dvdyc.slope.value,'\n\n')
    
#     dvdxh.slope.value * (0.125*u.deg.to(u.rad)*fcat['distance<pc>'] #km/s/deg * deg/pc
    
#     print(fcat['name'][i],hplane.slope_x,hplane.slope_y,
#           cplane.slope_x,cplane.slope_y,'\n\n\n\n\n')

    
#     plt.scatter(np.ravel(xh),np.ravel(h1),s=1)
#     plt.scatter(np.ravel(xc),np.ravel(c1),s=3)
#     plt.show()

    
new['vel_grad_GMC<km/s/pc>'] = OmCs
new['vel_grad_HI<km/s/pc>'] = OmHs

new['rot_ax_GMC<deg>'] = ThCs-90
new['rot_ax_HI<deg>'] = ThHs-90

for i in range(22):
    if new['rot_ax_GMC<deg>'][i] < -180:
        new['rot_ax_GMC<deg>'][i] +=360
        
    if new['rot_ax_HI<deg>'][i] < -180:
        new['rot_ax_HI<deg>'][i] +=360


for i in range(22):    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
    hmom1_old = fns.get_mom1(his_old[i],
                             np.arange(vo[i],vi[i])*1.3-100)
    hmom1_new = fns.get_mom1(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    cmom1 = fns.get_mom1(cos[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    plt.figure(figsize=(10,10))
    
#     plt.contourf(hmom1_old,
#                  extent=[l0ho,l1ho,b0ho,b1ho],
#                  levels = 10,
#                 cmap='RdGy',linewidths=0.5)
    
    
    plt.contourf(hmom1_new,
                 extent=[l0h,l1h,b0h,b1h],
                 levels = 10,
                cmap='PiYG')
    
    
    plt.contourf(cmom1,extent=[l0c,l1c,b0c,b1c],levels=10,
                cmap='binary_r',alpha=0.5)
    
# #     plt.colorbar()
    
    plt.xlim(max(l0c,l0h),min(l1c,l1h))
    plt.ylim(min(b0c,b0h),max(b1c,b1h))

    
    plt.xlabel(r'galactic longitude [$\degree$]')
    plt.ylabel(r'galactic latitude [$\degree$]')
    plt.show()


fcat['name'][[12,13,15,16,18,21]]


fig,axs = plt.subplots(2,3,figsize=(15,10))
for ax,i in zip(axs.flat,[12,13,15,16,18,21]):
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
    hmom1_new = fns.get_mom1(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    cmom1 = fns.get_mom1(cos[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    
    vmn = np.nanmean(cmom1)
    vst = np.nanstd(cmom1)
    ax.contourf(hmom1_new,
                 extent=[l0h,l1h,b0h,b1h],
                 levels = np.arange(vmn-5*vst,vmn+5*vst,0.5*vst),
                cmap='RdGy')
    
    
    ax.contourf(cmom1,extent=[l0c,l1c,b0c,b1c],
                levels=np.arange(vmn-2.5*vst,vmn+2.5*vst,0.5*vst),
                cmap='PiYG')
    
    ax.contour(np.log10(np.nanmean(cos[i],axis=0)*1.97e20*1.3),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(19,20,0.2),linewidths=0.5,alpha=0.4,
                colors='black')
    
    
    ax.set_xlim(max(l0c,l0h),min(l1c,l1h))
    ax.set_ylim(min(b0c,b0h),max(b1c,b1h))
    
    ax.set_yticks([])
    ax.set_xticks([])
    
plt.subplots_adjust(hspace=0,wspace=0)
# plt.savefig('postagemom1')
# fig.supxlabel('galactic longitude [$\degree$]')

    
#     plt.xlabel(r'galactic longitude [$\degree$]')
#     plt.ylabel(r'galactic latitude [$\degree$]')
#     plt.show()


np.nanmean(cmom1),np.nanstd(cmom1)


for i in [20]:    
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
    hmom1_old = fns.get_mom1(his_old[i],
                             np.arange(vo[i],vi[i])*1.3-100)
    hmom1_new = fns.get_mom1(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    cmom1 = fns.get_mom1(cos[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    plt.figure(figsize=(15,10))

    vmn = 14.1
    vst = 3.85
    
    plt.contourf(hmom1_new,
                 extent=[l0h,l1h,b0h,b1h],
                 levels = np.arange(vmn-0.75*vst,vmn+0.75*vst,0.1*vst) ,
                cmap='RdGy',extend='both')
    
    plt.colorbar(pad=5e-2,label=r'$v_{lsr} [km\ s^{-1}$]',ticks=[12,14,16])
    
    
#     plt.contourf(cmom1,
#                  extent=[l0c,l1c,b0c,b1c],
#                  levels = np.arange(vmn-2*vst,vmn+2*vst,0.5*vst),
#                 cmap='PiYG',extend='both')
    
#     plt.colorbar(pad=5e-2,ticks=[6,10,14,18,22])

    plt.contour(np.log10(np.nanmean(cos[i],axis=0)*1.97e20*1.3),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(19,20.5,0.15),linewidths=0.5,
                colors='black')
    
    
    plt.xlim(max(l0c,l0h),min(l1c,l1h))
    plt.ylim(min(b0c,b0h),max(b1c,b1h))

    
    plt.xlabel(r'galactic longitude [$\degree$]')
    plt.ylabel(r'galactic latitude [$\degree$]')
    
    plt.title('atomic gas [21cm]')
    
    plt.xticks([-151,-153,-155],[-151,-153,-155])
    plt.yticks([-4,-2,-0],[-4,-2,0])
    
    plt.savefig('rosette_mom1_h')
    plt.show()


for i in [13]:    
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
    hmom1_old = fns.get_mom1(his_old[i],
                             np.arange(vo[i],vi[i])*1.3-100)
    hmom1_new = fns.get_mom1(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    cmom1 = fns.get_mom1(cos[i],
                np.arange(vo[i],vi[i])*1.3-100)
    
    plt.figure(figsize=(17,10))

    
    plt.contourf(hmom1_new,
                 extent=[l0h,l1h,b0h,b1h],
                 levels = np.arange(10.86-3*1.53,10.86+3*1.53,0.5*1.53),
                cmap='RdGy',extend='both')
    
    plt.colorbar(pad=5e-2,label=r'$v_{lsr} [km\ s^{-1}$]')
    
    
    plt.contourf(cmom1,
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(10.86-3*1.53,10.86+3*1.53,0.5*1.53),
                cmap='PiYG',extend='both')
    
    plt.colorbar(ticks=[],pad=5e-2)

    plt.contour(np.log10(np.nanmean(cos[i],axis=0)*1.97e20*1.3),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(19,19.5,0.1),linewidths=0.5,
                colors='black')
    
    
    plt.xlim(max(l0c,l0h),min(l1c,l1h))
    plt.ylim(min(b0c,b0h),max(b1c,b1h))

    
    plt.xlabel(r'galactic longitude [$\degree$]')
    plt.ylabel(r'galactic latitude [$\degree$]')
    
    plt.title(fcat['name'][i])
    
    plt.xticks([-142,-145,-148],[-142,-145,-148])
    plt.yticks([-9,-12,-15],[-9,-12,-15])
    plt.show()


np.array((fcat['name'],Lgmc/Lhi > 1)).transpose()


plt.scatter(new['name'],(radii/radii_hi),s=20,color='black')

plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.ylabel(r'$R_{GMC}\ /\ R_{HI}$')
plt.show()


plt.rcParams['font.size']=15


plt.figure(figsize=(10,10))

plt.scatter(new['name'],new['vel_grad_GMC<km/s/pc>'],s=20,color='black')

plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.ylabel(r'$\Omega_{GMC}\ [km\ s^{-1}\ pc^{-1}]$')
# plt.savefig('Omega')
plt.show()


plt.figure(figsize=(10,10))

plt.scatter(new['name'],(radii/radii_hi),s=20,color='black')

plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.ylabel(r'$R_{GMC}\ /\ R_{HI}$')
plt.savefig('RR')
plt.show()


alphavir_gmc = ((((fcat['vel_dispersion_GMC<km/s>']*u.km/u.s)**2)*\
(new['radius_GMC<pc>'])/(c.G*new['mass_GMC<Msun>']*u.solMass).si).value)*((1*u.km**2 *1* u.pc)/(1*u.m**3)).si


betarot_gmc = (((new['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)**2 *\
               (new['radius_GMC<pc>'])**3)/(c.G*new['mass_GMC<Msun>']*u.solMass)).si


plt.figure(figsize=(10,10))
plt.scatter(new['name'],betarot_gmc,s=20,color='black')


plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.ylabel(r'$\beta_{rot,GMC}\ $')
plt.savefig('betarot')
plt.show()


plt.figure(figsize=(10,10))

plt.scatter(new['name'],alphavir_gmc,s=20,color='black')


plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.ylabel(r'$\alpha_{vir,GMC}\ $')
plt.savefig('alphavir')

plt.show()


jgmc = (new['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc*(new['radius_GMC<pc>']**2)).cgs
jhi = (new['vel_grad_HI<km/s/pc>']*u.km/u.s/u.pc*(new['radius_HI<pc>']**2)).cgs


Lgmc = new['mass_GMC<Msun>']*u.solMass*jgmc
Lhi = new['mass_HI<Msun>']*u.solMass*jhi


plt.hist(Lgmc.to(u.solMass*u.pc**2/u.s).value,
         color='xkcd:light blue',
         bins = [1e-13,1e-11,
                 1e-9,1e-7,1e-5,],
        alpha=0.7,edgecolor='xkcd:blue')


# plt.ylim(0,13)
plt.yticks([0,5,10],[0,5,10])

plt.xlabel(r'$L_{GMC}\ [M_\odot\ pc^{2}\ s^{-1}]$')

plt.ylabel(r'$N_{cloud}$')
plt.xscale('log')
plt.show()


plt.rcParams['font.size']=26


Lgmc.to(u.solMass/u.s*u.pc**2)


fttr(init,
     np.log10(fcat['vel_dispersion_GMC<km/s>']),
     np.log10(Lgmc.to(u.solMass/u.s*u.pc**2).value))


plt.figure(figsize=(10,10))

plt.scatter(fcat['vel_dispersion_GMC<km/s>'],
            Lgmc.to(u.solMass/u.s*u.pc**2),s=20,color='black')


# plt.xscale('log')

plt.ylabel(r'$L_{GMC}\ [M_\odot\ pc^{-2}\ s^{-1}]$')
plt.xlabel('$\sigma_{_{LOS}}\ [km\ s^{-1}$]')

plt.plot(np.arange(1e-1,7e0,1e0),
         10**(-9.94)*np.arange(1e-1,7e0,1e0)**(5.81),
        color='red',linestyle='--',label='$L \propto \sigma^{5.8}$')

plt.loglog()

plt.xlim(0.9,6)
plt.ylim(1e-12,1e-5)

plt.legend()
# plt.savefig('Lsigma')

plt.show()


plt.scatter(new['radius_GMC<pc>']/new['radius_HI<pc>'],
            new['vel_grad_GMC<km/s/pc>']/new['vel_grad_HI<km/s/pc>'],
           color='black',s=20)

plt.loglog()

plt.xlabel(r'$R_{GMC}\ /\ R_{HI}$')
plt.ylabel(r'$\Omega_{GMC}\ /\ \Omega_{HI}$')

plt.xticks([0.4,0.5,0.6],[0.4,0.5,0.6])
plt.plot(np.arange(0.35,0.7,0.1),np.arange(0.35,0.7,0.1)**-2,
        color='red',linestyle='--',label=r'$R^2\Omega\ $ = constant')

plt.legend(fontsize='x-small')
plt.show()


plt.figure(figsize=(10,10))

plt.hist(masses_hi/masses,color='xkcd:light blue',bins = np.arange(0,8,1),
        alpha=0.7,edgecolor='xkcd:blue')


plt.ylim(0,13)
plt.yticks([0,5,10],[0,5,10])

plt.xlabel(r'$M_{HI}\ /\ M_{GMC}$')

plt.ylabel(r'$N_{cloud}$')

# plt.savefig('MMhist')
plt.show()


plt.figure(figsize=(10,10))
plt.hist(jgmc/jhi,bins=np.arange(0,2,0.1),color='xkcd:light blue',
        alpha=0.7,edgecolor='xkcd:blue')


plt.xlabel(r'$j_{GMC}\ /\ j_{HI}$')

plt.ylabel(r'$N_{cloud}$')

plt.show()


plt.figure(figsize=(10,10))

plt.hist(Lgmc/Lhi,bins=np.arange(0,2,0.1),color='xkcd:light blue',
        alpha=0.7,edgecolor='xkcd:blue')


plt.xlabel(r'$L_{GMC}\ /\ L_{HI}$')

plt.ylabel(r'$N_{cloud}$')
plt.savefig('LLhist')

plt.show()


delthet = new['rot_ax_GMC<deg>'] - new['rot_ax_HI<deg>']

for i in range(22):
    if delthet[i] > 180:
        delthet[i] -= 360
            
    elif delthet[i] < -180:
        delthet[i] += 360


plt.figure(figsize=(10,10))

plt.hist(delthet,bins=np.arange(-180,210,30),color='xkcd:light blue',
        alpha=0.7,edgecolor='xkcd:blue')


plt.xlabel(r'$\Delta \phi$')

plt.ylabel(r'$N_{cloud}$')
plt.vlines([-90,90],0,6,linestyle='--',color='red')

plt.xticks([-180,-90,0,90,180],[-180,-90,0,90,180])

plt.savefig('delthet')
plt.show()


new[['name','lon_GMC<deg>',
    'lat_GMC<deg>','vel<km/s>','distance<pc>']].pprint()


Sigmah = (new['mass_HI<Msun>']*u.solMass)/(np.pi*(new['radius_HI<pc>'])**2)


Sigmah2 = (new['mass_GMC<Msun>']*u.solMass)/(np.pi*(new['radius_GMC<pc>'])**2)


plt.figure(figsize=(10,10))


# plt.hist(np.log10(jgmc.value),fill=False,bins=np.arange(20,28,0.5),
#        edgecolor='xkcd:black',label=r'$H_2$')

plt.hist(np.log10(jhi.value),color='xkcd:light blue',bins=np.arange(20,28,0.5),
        alpha=0.5,label='$HI$')#,edgecolor='xkcd:blue')


plt.xlabel(r'log $j\ [cm^2\ s^{-1}]$')

plt.ylabel(r'$N_{cloud}$')

# plt.vlines(10,0,6,color='red',linestyle='--')

plt.legend()
plt.xlim(22,27)

plt.savefig('jHIhist')


plt.show()


plt.figure(figsize=(10,10))


plt.hist(np.log10(jgmc.value),fill=False,bins=np.arange(20,28,0.5),
       edgecolor='xkcd:black',label=r'$H_2$')

# plt.hist(np.log10(jhi.value),color='xkcd:light blue',bins=np.arange(20,28,0.5),
#         alpha=0.5,label='$HI$')#,edgecolor='xkcd:blue')


plt.xlabel(r'log $j\ [cm^2\ s^{-1}]$')

plt.ylabel(r'$N_{cloud}$')

# plt.vlines(10,0,6,color='red',linestyle='--')

plt.legend()
plt.xlim(22,27)

plt.savefig('jH2hist')


plt.show()


Lgmc.cgs


plt.figure(figsize=(10,10))


# plt.hist(np.log10(Lgmc.cgs.value),fill=False,bins=np.arange(58,67),
#        edgecolor='xkcd:black',label=r'$H_2$')

plt.hist(np.log10(Lhi.cgs.value),color='xkcd:light blue',bins=np.arange(58,66,0.75),
        alpha=0.5,label='$HI$')#,edgecolor='xkcd:blue')


plt.xlabel(r'log $L\ [erg\ s]$')

plt.ylabel(r'$N_{cloud}$')

# plt.vlines(10,0,6,color='red',linestyle='--')

plt.legend()
# plt.xlim(22,27)

plt.savefig('LHIhist')


plt.show()


plt.figure(figsize=(10,10))


plt.hist(np.log10(Lgmc.cgs.value),fill=False,bins=np.arange(58,66,0.75),
       edgecolor='xkcd:black',label=r'$H_2$')

# plt.hist(np.log10(jhi.value),color='xkcd:light blue',bins=np.arange(20,28,0.5),
#         alpha=0.5,label='$HI$')#,edgecolor='xkcd:blue')


plt.xlabel(r'log $L\ [erg\ s]$')

plt.ylabel(r'$N_{cloud}$')

# plt.vlines(10,0,6,color='red',linestyle='--')

plt.legend()
# plt.xlim(22,27)

plt.savefig('LH2hist')


plt.show()


plt.figure(figsize=(10,10))


plt.hist(Sigmah2.value,fill=False,bins=np.arange(5,80,5),
       edgecolor='xkcd:black',label=r'$H_2$')

plt.hist(Sigmah.value,color='xkcd:light blue',bins=np.arange(5,80,5),
        alpha=0.5,label='$HI$')#,edgecolor='xkcd:blue')


plt.xlabel(r'$\Sigma\ [M_\odot\ pc^{-2}]$')

plt.ylabel(r'$N_{cloud}$')

plt.vlines(10,0,6,color='red',linestyle='--')
# plt.savefig('SigmaHIhist')

plt.legend()

plt.show()


fttr(init,np.log10(new['mass_GMC<Msun>']),np.log10(jgmc.value))


mj_slope = 0.79
mj_int = 21.25


plt.plot(np.arange(1e2,1e6,1e3),10**(mj_int)*np.arange(1e2,1e6,1e3)**mj_slope,label='0.79')

plt.fill_between(np.arange(1e2,1e6,1e3),
                 10**(mj_int)*np.arange(1e2,1e6,1e3)**0.5,
                 10**(mj_int)*np.arange(1e2,1e6,1e3)**0.9,alpha=0.2,
                 label='0.5-0.9')


plt.scatter(new['mass_GMC<Msun>'],jgmc,label='data')


plt.loglog()

plt.xlabel('$M_{GMC}\ [M_\odot]$')
plt.ylabel('$j_{GMC}\ [cm^{2}\ s^{-1}]$')
plt.legend()


plt.show()


#fttr(init,np.log10(new['radius_GMC<pc>'].value),np.log10(new['vel_grad_GMC<km/s/pc>']*new['radius_GMC<pc>'].value**2))


##a fitting relation comparing with theory -- remind myself where this came from / what it is referring to

# plt.fill_between(np.arange(1e0,1e2,1e1),10**(-0.77)*np.arange(1e0,1e2,1e1)**1.6,
#                  10**(-0.77)*np.arange(1e0,1e2,1e1)**2.4,alpha=0.2)


# plt.scatter(new['radius_GMC<pc>'],
#             new['vel_grad_GMC<km/s/pc>']*new['radius_GMC<pc>']**2)

# plt.loglog()

# plt.plot(np.arange(1e0,1e2,1e1),10**(-0.77)*np.arange(1e0,1e2,1e1)**1.74)

# # plt.plot(np.arange(1e0,1e2,1e1),10**(-0.77)*np.arange(1e0,1e2,1e1)**2.4)


# #these are for r not R (ie r starts at 0 in the center of the object and steps outward)
# # plt.plot(np.arange(1e0,1e2,1e1),10**(-0.77)*np.arange(1e0,1e2,1e1)**1.5)
# # plt.plot(np.arange(1e0,1e2,1e1),10**(-0.77)*np.arange(1e0,1e2,1e1)**2)


# plt.show()


new[20]


np.array((fcat['name'],
          new['vel_grad_HI<km/s/pc>']/fcat['vel_grad_HI<km/s/pc>'])).transpose()


np.array((new['rot_ax_HI<deg>'],fcat['rot_ax_HI<deg>'])).transpose()


fcat


for i in range(22):    
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
#     hmom2_old = fns.get_momN(his_old[i],
#                              np.arange(vo[i],vi[i])*1.3-100,2)
    hmom2_new = np.sqrt(fns.get_momN(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100,2))
    
    cmom2 = np.sqrt(fns.get_momN(cos[i],
                np.arange(vo[i],vi[i])*1.3-100,2))
    
    plt.figure(figsize=(17,10))

    
    plt.contourf(hmom2_new,
                 extent=[l0h,l1h,b0h,b1h],
                 levels =10,
                cmap='YlGnBu_r',extend='both')
    
    plt.colorbar(pad=5e-2,label=r'$v_{lsr} [km\ s^{-1}$]')
    
    
    plt.contourf(cmom2,
                 extent=[l0c,l1c,b0c,b1c],
                 levels = 10,
                cmap='magma',extend='both')
    
    plt.colorbar(pad=5e-2)

    plt.contour(np.log10(np.nanmean(cos[i],axis=0)*1.97e20*1.3),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = 10,
                colors='black')
    
    
    plt.xlim(max(l0c,l0h),min(l1c,l1h))
    plt.ylim(min(b0c,b0h),max(b1c,b1h))

    
    plt.xlabel(r'galactic longitude [$\degree$]')
    plt.ylabel(r'galactic latitude [$\degree$]')
    
    plt.title(fcat['name'][i])
    
#     plt.xticks([-142,-145,-148],[-142,-145,-148])
#     plt.yticks([-9,-12,-15],[-9,-12,-15])
    plt.show()


plt.contourf(np.nansum(hcube[90:,:,:],axis=0),extent=[l0h,l1h,b0h,b1h],
             cmap=dc.cmap_min0('viridis'))
plt.colorbar()


plt.contour(np.nansum(ccube[:,:,:],axis=0),extent=[l0c,l1c,b0c,b1c],
             cmap=dc.cmap_min0('binary_r'),levels=10)


plt.contourf(np.nansum(hcube[:70,:,:],axis=0),vmin=20,
             cmap=dc.cmap_min0('viridis'),extent=[l0h,l1h,b0h,b1h])
plt.colorbar()

plt.contour(np.nansum(ccube[:,:,:],axis=0),extent=[l0c,l1c,b0c,b1c],
             cmap=dc.cmap_min0('binary_r'),levels=10)


plt.contourf(np.nansum(hcube[70:90,:,:],axis=0),extent=[l0h,l1h,b0h,b1h],
             vmin=50,cmap=dc.cmap_min0('viridis'))
plt.colorbar()

plt.contour(np.nansum(ccube[:,:,:],axis=0),extent=[l0c,l1c,b0c,b1c],
             cmap=dc.cmap_min0('binary_r'),levels=10)


plt.plot(np.nanmean(hcube,axis=(1,2)))
# plt.vlines([50,70,90],0,60)
# plt.ylim(40,)
plt.xlim(30,130)


plt.plot(np.nanmean(ccube,axis=(1,2)))
plt.vlines([50,70,90],0,0.5)
plt.xlim(30,130)
# plt.ylim(0,0.3)


np.array((92,107))*1.3-100 #maddalena


np.array((82,100))*1.3-100 #cma


np.array((70,115))*1.3-100 


fcat['name'][11]


100*1.3-100


np.array([50,90])*1.3-100


(-40+100)/1.3


for i in [21]:

    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
#     print(boh[i],bih[i],loh[i],lih[i])
    hcube=h3sig[:,boh[i]:bih[i],loh[i]:lih[i]]
    ccube=d3sig[:,boc[i]:bic[i],loc[i]:lic[i]]
    
    plt.contourf(np.log10(np.nansum(hcube,axis=1)*0.125*1.82e18),
                 extent=[l0h,l1h,-100,100],levels=np.arange(20,21.5,0.25),
                cmap='Paired')
    plt.colorbar()
    plt.contour(np.log10(np.nansum(ccube,axis=1)*0.125*1.97e20),
                extent=[l0c,l1c,-100,100],levels=np.arange(20,21.5,0.25))
#     plt.ylim(-50,50)
    plt.ylim(-40,40)
    plt.hlines([-15],l0h,l1h,color='red')
    plt.title(fcat['name'][i])
    plt.show()
    

for i in [19]:
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i]-100,bih[i]+100))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
#     print(boh[i],bih[i],loh[i],lih[i])
    hcube=h3sig[:,boh[i]:bih[i],loh[i]:lih[i]]
    ccube=d3sig[:,boc[i]:bic[i],loc[i]:lic[i]]
    
    plt.contourf(np.log10(np.nansum(hcube,axis=2)*0.125*1.82e18),
                 extent=[b0h,b1h,-100,100],levels=np.arange(19,23,0.25),
                cmap='Paired')
    plt.colorbar()
    plt.contour(np.nansum(ccube,axis=2),extent=[b0c,b1c,-100,100])
#     plt.ylim(-50,50)
    plt.ylim(-60,40)
#     plt.hlines([-10,-15,20,25,30],b0h,b1h,color='red')
    plt.title(fcat['name'][i])
    plt.show()
    

for i in range(22):
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i]-100,bih[i]+100))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
#     print(boh[i],bih[i],loh[i],lih[i])
    hcube=h3sig[:,boh[i]:bih[i],loh[i]:lih[i]]
    ccube=d3sig[:,boc[i]:bic[i],loc[i]:lic[i]]
    
    plt.contourf(np.log10(np.nansum(hcube,axis=2)*0.125*1.82e18),
                 extent=[b0h,b1h,-100,100],levels=np.arange(20,23,0.25),
                cmap='Paired')
    plt.colorbar()
    plt.contour(np.nansum(ccube,axis=2),extent=[b0c,b1c,-100,100])
#     plt.ylim(-50,50)
    plt.ylim(fcat['vel<km/s>'][i][0]-5,fcat['vel<km/s>'][i][1]+5)
#     plt.hlines([-9,-20,-30],l0h,l1h,color='red')
    plt.title(fcat['name'][i])
    plt.show()
    

for i in range(22):
    l0h,l1h = np.array((loh[i]-100,lih[i]+100))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
#     print(boh[i],bih[i],loh[i],lih[i])
    hcube=h3sig[:,boh[i]:bih[i],loh[i]:lih[i]]
    ccube=d3sig[:,boc[i]:bic[i],loc[i]:lic[i]]
    
    plt.contourf(np.log10(np.nansum(hcube,axis=1)*0.125*1.82e18),
                 extent=[l0h,l1h,-100,100],levels=np.arange(20,23,0.25),
                cmap='Paired')
    plt.colorbar()
    plt.contour(np.nansum(ccube,axis=1),extent=[l0c,l1c,-100,100])
#     plt.ylim(-50,50)
    plt.ylim(fcat['vel<km/s>'][i][0]-5,fcat['vel<km/s>'][i][1]+15)
#     plt.hlines([-9,-20,-30],l0h,l1h,color='red')
    plt.title(fcat['name'][i])
    plt.show()
    

new[1]


np.nanmean((257,287))


2614-2589


for i in range(22):    
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0ho,l1ho = np.array((loh_old[i],lih_old[i]))*-0.125+360
    b0ho,b1ho = np.array((boh_old[i],bih_old[i]))*0.125-50
    
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    
    hmom2_old = np.sqrt(fns.get_momN(his_old[i],
                             np.arange(vo[i],vi[i])*1.3-100,2))
    hmom2_new = np.sqrt(fns.get_momN(his_new[i],
                np.arange(vo[i],vi[i])*1.3-100,2))
    
    cmom2 = np.sqrt(fns.get_momN(cos[i],
                np.arange(vo[i],vi[i])*1.3-100,2))
    
#     np.sqrt(fns.get_momN(np.nanmean(his_new[i],axis=0)))
#     spec = np.nanmean(his_old[i],axis=(1,2))
    spec2 = np.nanmean(his_new[i],axis=(1,2))
    spec3 = np.nanmean(h3sig[:,boh[i]:bih[i],loh[i]:lih[i]],axis=(1,2))
    spec4 = np.nanmean(d3sig[:,boc[i]:bic[i],loc[i]:lic[i]],axis=(1,2))
    spec5 = np.nanmean(cos[i],axis=(1,2))
#     spec3 = np.nanmean(cos[i],axis=(1,2))
#     meanmom2 = np.sqrt(fns.get_momN(spec,np.arange(vo[i],vi[i])*1.3-100,2))
    
#     print(fcat['name'][i])
#     print(np.nanmean(hmom2_old))#,meanmom2)#,np.nanmean(cmom2))
#     print(np.nanmean(hmom2_new),'\n')
    
#     plt.plot(spec)
#     plt.plot(np.arange(155)*1.3-100,spec3)
    plt.plot(np.arange(vo[i],vi[i])*1.3-100,spec2)
    plt.plot(np.arange(vo[i],vi[i])*1.3-100,100*spec5)

#     plt.plot(np.arange(155)*1.3-100,spec4*100)
    plt.title(fcat['name'][i])
#     plt.vlines([6,30],0,50)

#     plt.plot(spec3*80)
    plt.show()

    
#     plt.figure(figsize=(17,10))

    
#     plt.contourf(hmom2_new,
#                  extent=[l0h,l1h,b0h,b1h],
#                  levels =10,
#                 cmap='YlGnBu_r',extend='both')
    
#     plt.colorbar(pad=5e-2,label=r'$v_{lsr} [km\ s^{-1}$]')
    
    
#     plt.contourf(cmom2,
#                  extent=[l0c,l1c,b0c,b1c],
#                  levels = 10,
#                 cmap='magma',extend='both')
    
#     plt.colorbar(pad=5e-2)

#     plt.contour(np.log10(np.nanmean(cos[i],axis=0)*1.97e20*1.3),
#                  extent=[l0c,l1c,b0c,b1c],
#                  levels = 10,
#                 colors='black')
    
    
#     plt.xlim(max(l0c,l0h),min(l1c,l1h))
#     plt.ylim(min(b0c,b0h),max(b1c,b1h))

    
#     plt.xlabel(r'galactic longitude [$\degree$]')
#     plt.ylabel(r'galactic latitude [$\degree$]')
    
#     plt.title(fcat['name'][i])
    
# #     plt.xticks([-142,-145,-148],[-142,-145,-148])
# #     plt.yticks([-9,-12,-15],[-9,-12,-15])
#     plt.show()


len(newvels)


np.array((fcat['vel_dispersion_HI<km/s>'],fcat['vel_dispersion_GMC<km/s>'])).transpose()


np.sqrt(fcat['mass_GMC<kg>']*u.kg/(np.pi*10*u.solMass/u.pc**2)).to(u.pc)


# import numpy as np
# import pandas as pd

# import matplotlib.pyplot as plt
# from matplotlib.patches import Rectangle,Circle

# import astropy.units as u
# import astropy.constants as c
# from astropy.io import fits
# from astropy.table import QTable
# from astropy.modeling import models,fitting

# import fun as fns
# import datacube as dc

# #---------------------------------------------------------------------------------------------
# fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


# #-------------------------------------------------------------------------------------------
# hi4pi,dht = fits.open('both_same_grid_720d.fits')
# dust = fits.open('mean_and_std_lbd.fits')[1].data

# h3sig = hi4pi.data*(hi4pi.data>0.043*3)
# d3sig = dht.data*(dht.data>0.2*3)

# h0 = np.nansum(h3sig,axis=0)*1.3*1.82e18
# c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20

# hlv = np.nansum(h3sig,axis=1)*0.125*1.82e18
# clv = np.nansum(d3sig,axis=1)*0.125*1.97e20

# #the 0.5 is to account for the fact that this is 720 deg
# hbv = np.nansum(h3sig,axis=2)*0.125*1.82e18*0.5 
# cbv = np.nansum(d3sig,axis=2)*0.125*1.97e20*0.5


# #---------------------------------------------------------------------------------------------
# lcs_d = np.nanmean(fcat['lon_GMC<deg>'],axis=1)
# bcs_d = np.nanmean(fcat['lat_GMC<deg>'],axis=1)

# ras_d = (fcat['radius_HI<pc>']/fcat['distance<pc>'])*u.rad.to(u.deg)


# loh,lih = ((np.array([lcs_d+ras_d,lcs_d-ras_d]).transpose() -360)/-0.125).astype('int').transpose()
# boh,bih = ((np.array([bcs_d-ras_d,bcs_d+ras_d]).transpose()+50)/0.125).astype('int').transpose()

# loc,lic = ((fcat['lon_GMC<deg>'] -360)/-0.125).astype('int').transpose()
# boc,bic = ((fcat['lat_GMC<deg>'] +50)/0.125).astype('int').transpose()

# vo,vi = np.array([((fcat['vel<km/s>'][i]+100)/1.3).astype('int') for i in range(22)]).transpose()


# vs = [np.arange((vi-vo)[i])*1.3 + fcat['vel<km/s>'][i][0] for i in range(22)]


# his = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
# cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]
# cos_hisz = [d3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
# his_cosz = [h3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


# #---------------------------------------------------------------------------------------------
