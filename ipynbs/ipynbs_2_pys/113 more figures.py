#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting

import fun as fns
import datacube as dc

#---------------------------------------------------------------------------------------------
fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


#-------------------------------------------------------------------------------------------
hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data

h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)

h0 = np.nansum(h3sig,axis=0)*1.3*1.82e18
c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20

hlv = np.nansum(h3sig,axis=1)*0.125*1.82e18
clv = np.nansum(d3sig,axis=1)*0.125*1.97e20

#the 0.5 is to account for the fact that this is 720 deg
hbv = np.nansum(h3sig,axis=2)*0.125*1.82e18*0.5 
cbv = np.nansum(d3sig,axis=2)*0.125*1.97e20*0.5


#---------------------------------------------------------------------------------------------
lcs_d = np.nanmean(fcat['lon_GMC<deg>'],axis=1)
bcs_d = np.nanmean(fcat['lat_GMC<deg>'],axis=1)

ras_d = (fcat['radius_HI<pc>']/fcat['distance<pc>'])*u.rad.to(u.deg)


loh,lih = ((np.array([lcs_d+ras_d,lcs_d-ras_d]).transpose() -360)/-0.125).astype('int').transpose()
boh,bih = ((np.array([bcs_d-ras_d,bcs_d+ras_d]).transpose()+50)/0.125).astype('int').transpose()

loc,lic = ((fcat['lon_GMC<deg>'] -360)/-0.125).astype('int').transpose()
boc,bic = ((fcat['lat_GMC<deg>'] +50)/0.125).astype('int').transpose()

vo,vi = np.array([((fcat['vel<km/s>'][i]+100)/1.3).astype('int') for i in range(22)]).transpose()


vs = [np.arange((vi-vo)[i])*1.3 + fcat['vel<km/s>'][i][0] for i in range(22)]


his = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]
cos_hisz = [d3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
his_cosz = [h3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


#---------------------------------------------------------------------------------------------
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


h0s = []
h1s = []
h2s = []

h0scsz = []
h1scsz = []
h2scsz = []

c0s = []
c1s = []
c2s = []

c0shsz = []
c1shsz = []
c2shsz = []

for i in range(22):
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    thish0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    thish1 = fns.get_mom1(his[i],vs[i])
    thish2 = np.sqrt(fns.get_momN(his[i],vs[i],2))
    
    
    thish0c = np.nansum(his_cosz[i],axis=0)*1.3*1.82e18
    thish1c = fns.get_mom1(his_cosz[i],vs[i])
    thish2c = np.sqrt(fns.get_momN(his_cosz[i],vs[i],2))
    
    
    h0s.append(thish0)
    h1s.append(thish1)
    h2s.append(thish2)
    
    h0scsz.append(thish0c)
    h1scsz.append(thish1c)
    h2scsz.append(thish2c)
    
    
    thisc0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    thisc1 = fns.get_mom1(cos[i],vs[i])
    thisc2 = np.sqrt(fns.get_momN(cos[i],vs[i],2))
    
    thisc0h = np.nansum(cos_hisz[i],axis=0)*1.3*1.97e20
    thisc1h = fns.get_mom1(cos_hisz[i],vs[i])
    thisc2h = np.sqrt(fns.get_momN(cos_hisz[i],vs[i],2))
    
    
    c0s.append(thisc0)
    c1s.append(thisc1)
    c2s.append(thisc2)
    
    c0shsz.append(thisc0h)
    c1shsz.append(thisc1h)
    c2shsz.append(thisc2h)


fcat['radius_GMC<pc>']/fcat['radius_HI<pc>']


import matplotlib as mpl


plt.scatter(fcat['radius_GMC<pc>']/fcat['radius_HI<pc>'],
            fcat['vel_grad_GMC<km/s/pc>']/fcat['vel_grad_HI<km/s/pc>'],
            s=10,color='black')

plt.xlabel(r'$R_{GMC}\ /\ R_{HI}$')
plt.ylabel(r'$\Omega_{GMC}\ /\ \Omega_{HI}$')


plt.plot(np.arange(1e-2,10),np.arange(1e-2,10)**(-2),
         color='crimson',linewidth=0.8,label=r'$R^2 \Omega$ = const.')

plt.plot(np.arange(1e-2,10),np.arange(1e-2,10)**(1),linestyle='--',
         color='grey',linewidth=0.8,label=r'$R /\Omega$ = const.')

plt.loglog()

plt.xticks(ticks=[1e-1,3e-1,1],labels=[r'$10^{-1}$',r'$10^{-0.5}$',r'$10^{0}$'])
plt.ylim(1e-1,1e2)
plt.xlim(0.1,1)
plt.gca().xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
plt.legend(fontsize='x-small')
plt.savefig('omomrr',bbox_inches='tight')
plt.show()


plt.scatter(fcat['mass_GMC<kg>']/fcat['mass_HI<kg>'],
            fcat['spec_ang_mom_GMC<pc.km/s>']/fcat['spec_ang_mom_HI<pc.km/s>'],
            s=10,color='black')

# plt.xlabel(r'$R_{GMC}\ /\ R_{HI}$')
# plt.ylabel(r'$\Omega_{GMC}\ /\ \Omega_{HI}$')


# plt.plot(np.arange(1e-2,10),np.arange(1e-2,10)**(-2),
#          color='crimson',linewidth=0.8,label=r'$R^2 \Omega$ = const.')

# plt.plot(np.arange(1e-2,10),np.arange(1e-2,10)**(1),linestyle='--',
#          color='grey',linewidth=0.8,label=r'$R /\Omega$ = const.')

plt.loglog()

plt.xticks(ticks=[1e-2,1e-1,1],
           labels=[r'$10^{-2}$',r'$10^{-1}$',r'$10^{0}$'])
plt.ylim(1e-2,1e0)
# plt.xlim(0.1,1)
plt.gca().xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
# plt.legend(fontsize='x-small')
# plt.savefig('omomrr',bbox_inches='tight')
plt.show()


plt.scatter(fcat['radius_GMC<pc>']/fcat['radius_HI<pc>'],
            fcat['spec_ang_mom_GMC<pc.km/s>']/fcat['spec_ang_mom_HI<pc.km/s>'],
            s=10,color='black')

# plt.xlabel(r'$R_{GMC}\ /\ R_{HI}$')
# plt.ylabel(r'$\Omega_{GMC}\ /\ \Omega_{HI}$')


# plt.plot(np.arange(1e-2,10),np.arange(1e-2,10)**(-2),
#          color='crimson',linewidth=0.8,label=r'$R^2 \Omega$ = const.')

# plt.plot(np.arange(1e-2,10),np.arange(1e-2,10)**(1),linestyle='--',
#          color='grey',linewidth=0.8,label=r'$R /\Omega$ = const.')

plt.loglog()

plt.xticks(ticks=[1e-1,3e-1,1],labels=[r'$10^{-1}$',r'$10^{-0.5}$',r'$10^{0}$'])
# plt.ylim(1e-1,1e2)
plt.xlim(0.1,1)
plt.gca().xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
# plt.legend(fontsize='x-small')
# plt.savefig('omomrr',bbox_inches='tight')
plt.show()


10**(-0.5)


h0 = np.nansum(h3sig,axis=0)*1.3*1.82e18


c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20


(np.array((-50,50))-360)/-0.125


cbv_l50 = np.nansum(d3sig[:,:,2480:3280],axis=2)*0.125*1.97e20


plt.figure(figsize=(20,17))

plt.contourf(np.log10(cbv_l50).transpose(),extent=[-100,100,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')

# plt.contour(np.log10(c0),extent=[360,-360,-50,50],
#              levels=[20.5,21,21.5],colors=['black','grey','white'])

# plt.xlim(180,-180)
plt.ylabel(r'galactic latitude [$\degree$]')
# plt.xlabel(r'galactic longitude [$\degree$]')
plt.xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')

# plt.xlim(50,-50)
plt.ylim(-10,30)
# plt.savefig('zeroth_h2_zoom',bbox_inches='tight')
plt.colorbar()


(np.array((-20,-10,0,10,20))+100)/1.3


plt.rcParams['font.size'] = 30


np.nanmean(h3sig[:,loh[i]:lih[i],boh[i]:bih[i]],axis=(1,2))


vels = np.arange(155)*1.3-100


plt.rcParams['font.size']=20


# for i in [4]:
i=4
fig,ax = plt.subplots(1,1)

ax.plot(vels,np.nanmean(h3sig[:,boh[i]:bih[i],loh[i]:lih[i]],axis=(1,2)),color='teal')
ax.set_ylim(0,70)
ax2 = ax.twinx()

ax2.plot(vels,np.nanmean(d3sig[:,boc[i]:bic[i],loc[i]:lic[i]],axis=(1,2)),color='crimson')
ax2.set_ylim(0,0.5)

plt.title(fcat['name'][i])

plt.xlim(-30,100)
plt.xticks([-30,0,30,60,90],[30,0,30,60,90])


ax.set_xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')


ax.set_ylabel(r'$T_{B,HI}\ $ [$K$]')
ax2.set_ylabel(r'$T_{B,H_2}\ $ [$K$]')

plt.show()


fig,axs = plt.subplots(2,2,figsize=(25,25),sharex=True,sharey=True)


axs[0,0].contourf(np.log10(np.nansum(d3sig[62:69,:,:],axis=0)*1.97e20*1.3),extent=[360,-360,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')


axs[0,1].contourf(np.log10(np.nansum(d3sig[69:77,:,:],axis=0)*1.97e20*1.3),extent=[360,-360,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')


axs[1,0].contourf(np.log10(np.nansum(d3sig[77:85,:,:],axis=0)*1.97e20*1.3),extent=[360,-360,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')


axs[1,1].contourf(np.log10(np.nansum(d3sig[85:92,:,:],axis=0)*1.97e20*1.3),extent=[360,-360,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')

# plt.contour(np.log10(c0),extent=[360,-360,-50,50],
#              levels=[20.5,21,21.5],colors=['black','grey','white'])

# # plt.xlim(180,-180)
fig.supylabel(r'galactic latitude [$\degree$]')
fig.supxlabel(r'galactic longitude [$\degree$]')
# plt.xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')

for ax in axs.flat:
    ax.set_xlim(50,-50)
    ax.set_ylim(-10,30)
    
    
plt.subplots_adjust(wspace=0,hspace=0)
# plt.ylim(-10,30)
# plt.savefig('zeroth_h2_zoom',bbox_inches='tight')
# plt.colorbar()


plt.figure(figsize=(26,17))

plt.contourf(np.log10(np.nansum(d3sig[60:70,:,:],axis=0)*1.97e20*1.3),extent=[360,-360,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')

# plt.contour(np.log10(c0),extent=[360,-360,-50,50],
#              levels=[20.5,21,21.5],colors=['black','grey','white'])

# plt.xlim(180,-180)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')
# plt.xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')

plt.xlim(50,-50)
plt.ylim(-10,30)
# plt.savefig('zeroth_h2_zoom',bbox_inches='tight')
# plt.colorbar()


plt.figure(figsize=(26,17))

plt.contourf(np.log10(np.nansum(d3sig[80:90,:,:],axis=0)*1.97e20*1.3),extent=[360,-360,-50,50],
             levels=np.arange(20,23,0.2),cmap='magma')

# plt.contour(np.log10(c0),extent=[360,-360,-50,50],
#              levels=[20.5,21,21.5],colors=['black','grey','white'])

# plt.xlim(180,-180)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')
# plt.xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')

plt.xlim(50,-50)
plt.ylim(-10,30)
# plt.savefig('zeroth_h2_zoom',bbox_inches='tight')
# plt.colorbar()


plt.rcParams['font.size'] = 30


plt.figure(figsize=(30,13))

plt.contourf(np.log10(h0),extent=[360,-360,-50,50],
             levels=np.arange(19,22.5,0.2),cmap='YlGnBu_r')

# plt.contour(np.log10(h0),extent=[360,-360,-50,50],
#              levels=[20.5,21,21.5],colors=['black','grey','white'])

plt.xlim(180,-180)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

plt.savefig('zeroth_hi',bbox_inches='tight')
plt.show()


plt.rcParams['font.size'] = 20


plt.figure(figsize=(8,8),dpi = 200)
plt.plot(np.nanmean(cos[19],axis=(1,2)))

plt.xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')
plt.xlim(18,30)
plt.ylim(0,0.35)
plt.ylabel(r'$T_B\ $ [$K$]')
plt.title('polaris')


fcat['name'][19]


plt.figure(figsize=(33,13))

plt.contourf(np.log10(h0),extent=[360,-360,-50,50],
             levels=np.arange(19,22.5,0.2),cmap='binary_r')


for i in range(22):
    plt.contourf(np.log10(h0s[i]),extent=[loh[i]*-0.125+360,lih[i]*-0.125+360,
                                         boh[i]*0.125-50,bih[i]*0.125-50],
                cmap='YlGnBu_r',levels=np.arange(20,22,0.2))

    
plt.contourf(np.log10(c0),extent=[360,-360,-50,50],
             levels=np.arange(20.5,23,0.2),cmap='binary_r')


for i in range(22):
    plt.contourf(np.log10(c0s[i]),extent=[loc[i]*-0.125+360,lic[i]*-0.125+360,
                                         boc[i]*0.125-50,bic[i]*0.125-50],
                cmap='magma',levels=np.arange(20,22.5,0.2))


plt.xlim(180,-180)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

plt.savefig('catalog',bbox_inches='tight')


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()


mjfit = fttr(init,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)),
             np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']))


mjfit


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
            fcat['spec_ang_mom_GMC<pc.km/s>'],s=20,color='black')

plt.loglog()

plt.xlabel(r'$M_{GMC}\ $ [$M_\odot$]')
plt.ylabel(r'$j_{GMC}\ $ [$pc\ km\ s^{-1}$]')


plt.plot(np.arange(1e2,1e6,1e3),
         10**(-2.6306993)*np.arange(1e2,1e6,1e3)**(0.5),
        color='orange',linewidth=0.7,linestyle='--',alpha=0.6,
         label=r'$j \propto M^{0.5}$')


plt.plot(np.arange(1e2,1e6,1e3),
         10**(-2.6306993)*np.arange(1e2,1e6,1e3)**(0.79128551),
        color='red',linewidth=0.7,label=r'$j \propto M^{0.79}$')


plt.plot(np.arange(1e2,1e6,1e3),
         10**(-2.6306993)*np.arange(1e2,1e6,1e3)**(0.9),
        color='purple',linewidth=0.7,linestyle='--',alpha=0.6,
         label=r'$j \propto M^{0.9}$')


plt.legend(fontsize='x-small')

plt.show()


JJfit = fttr(init,np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)*\
                           fcat['spec_ang_mom_HI<pc.km/s>']*u.km.to(u.pc)),
            np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)*\
                     fcat['spec_ang_mom_GMC<pc.km/s>']*u.km.to(u.pc)))


JJfit


JJfit


plt.scatter(fcat['mass_HI<kg>']*u.kg.to(u.solMass)*\
            fcat['spec_ang_mom_HI<pc.km/s>']*u.km.to(u.pc),
            fcat['mass_GMC<kg>']*u.kg.to(u.solMass)*\
            fcat['spec_ang_mom_GMC<pc.km/s>']*u.km.to(u.pc),
            s=20,color='black')

plt.loglog()

plt.xlabel(r'$L_{HI}\ $ [$M_\odot\ pc^2\ s^{-1}$]')
plt.ylabel(r'$L_{GMC}\ $ [$M_\odot\ pc^2\ s^{-1}$]')

plt.plot(np.arange(1e-9,1e-4,1e-8),
         10**(-0.27299537)*np.arange(1e-9,1e-4,1e-8)**(1.2399194),
        color='red',linewidth=0.7,label=r'$L_{GMC} \propto L_{HI}^{1.2}$')


plt.plot(np.arange(1e-9,1e-4,1e-8),
         10**(-3.5299537)*np.arange(1e-9,1e-4,1e-8)**(1),
        color='orange',linewidth=0.7,label=r'$L_{GMC} \propto L_{HI}$',
        linestyle='--',alpha=0.6)


plt.legend(fontsize='x-small')

plt.show()


Lhi = (fcat['spec_ang_mom_HI<pc.km/s>']*u.km.to(u.pc))*fcat['mass_HI<kg>']*u.kg
Lgmc = (fcat['spec_ang_mom_GMC<pc.km/s>']*u.km.to(u.pc))*fcat['mass_GMC<kg>']*u.kg


10**0.5


((10*u.solMass/u.pc**2)/(0.025*c.M_sun/(u.pc**3))).to(u.pc) 
#if we use the rho value given in imara & blitz, it implies that this 
#is the length that we are integrating over to get 10Msun/pc2

#but then we are also saying that the length of the cylinder (the distance over
#which we are integrating) is the same length as the diameter?

#and we're also saying that the cylinder is perpendicular to the galactic plane?
#because in that case we would be seeing a rectangle in projection, with length=diameter
#and width =diameter


#a formula for the accumulation radius that doesn't assume h=200pc 
#(bc this only works for faceone, not-milky way galaxies) 

# (pi*Ra^2 (circle area) * 2*RA (height)) (volume) * rho = M
# rho = Sigma / L = Sigma / 2RA
#M = pi * RA^2 * Sigma -> and here we can assume Sigma=10Msun/pc2 ... but 
#to do it more rigorously, this would need to look at the surface density of HI
#in the vicinity of the gmc

#so RA => sqrt(M/(pi*Sigma))


recalc_RAs = np.sqrt(fcat['mass_GMC<kg>']*u.kg/(np.pi*10*u.solMass/u.pc**2)).to(u.pc)


np.array((fcat['name'],recalc_RAs,fcat['radius_HI<pc>'])).transpose()


plt.scatter(fcat['name'],recalc_RAs/fcat['radius_HI<pc>'])

plt.hlines([1,0.5,0.1],'aquila_s','taurus')

plt.xticks(rotation='vertical')
plt.show()


plt.hist(Lgmc/Lhi,bins=[0.00032,0.0006,0.001,0.0018,0.0032,0.0056,
                        0.01,0.018,0.032,0.056,0.1,0.18,],color='xkcd:light blue',
                    alpha=0.7,edgecolor='xkcd:blue')

plt.xscale('log')

plt.xticks([1e-3,1e-2,1e-1],[r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$'])
plt.show()


plt.hist(Lgmc/Lhi,bins=[0.00032,0.0006,0.001,0.0018,0.0032,0.0056,
                        0.01,0.018,0.032,0.056,0.1,0.18,],color='xkcd:light blue',
                    alpha=0.7,edgecolor='xkcd:blue')

plt.xscale('log')

plt.xticks([1e-3,1e-2,1e-1],[r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$'])
plt.show()


plt.hist(fcat['spec_ang_mom_HI<pc.km/s>']*u.km.to(u.pc)/
            (fcat['spec_ang_mom_GMC<pc.km/s>']*u.km.to(u.pc)),
         bins=10**np.arange(-0.5,4,0.25),alpha=0.5,
        color='xkcd:grey',edgecolor='xkcd:black')


plt.hist(fcat['mass_HI<kg>']*u.kg.to(u.solMass)*\
            fcat['spec_ang_mom_HI<pc.km/s>']*u.km.to(u.pc)/
            (fcat['mass_GMC<kg>']*u.kg.to(u.solMass)*\
            fcat['spec_ang_mom_GMC<pc.km/s>']*u.km.to(u.pc)),
         bins=10**np.arange(-0.5,4,0.25),alpha=0.5,
        color='xkcd:light blue',edgecolor='xkcd:blue')


plt.xscale('log')
plt.show()


fttr(init,np.log10(fcat['spec_ang_mom_HI<pc.km/s>']),
    np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']))


plt.scatter(fcat['spec_ang_mom_HI<pc.km/s>'],fcat['spec_ang_mom_GMC<pc.km/s>'])

plt.loglog()

plt.show()


plt.hist(fcat['mass_GMC<kg>']/fcat['mass_HI<kg>'])


plt.hist((fcat['mass_GMC<kg>']*u.kg.to(u.solMass)*\
 fcat['spec_ang_mom_GMC<pc.km/s>'])/\
(fcat['mass_HI<kg>']*u.kg.to(u.solMass)*fcat['spec_ang_mom_HI<pc.km/s>']))


(u.solMass*u.pc*u.km/u.s).si


jgmchst = np.histogram(fcat['spec_ang_mom_GMC<pc.km/s>'],bins=10**np.arange(0,3,0.25))
jhihst = np.histogram(fcat['spec_ang_mom_HI<pc.km/s>'],bins=10**np.arange(0,3,0.25))


plt.hist(fcat['spec_ang_mom_HI<pc.km/s>'],bins=10**np.arange(0,3.5,0.25),
         edgecolor='black',alpha=0.4,color='grey')

plt.xscale('log')

plt.xlim(0.5,1e3)

plt.xlabel(r'$j_{HI}\ $ [$pc\ km\ s^{-1}$]')

plt.show()


plt.hist(fcat['spec_ang_mom_GMC<pc.km/s>'],bins=10**np.arange(0,3.5,0.25),
         color='xkcd:light blue',alpha=0.7,edgecolor='xkcd:blue')


# plt.hist(fcat['spec_ang_mom_HI<pc.km/s>'],bins=10**np.arange(0,3.5,0.25),
#          edgecolor='black',alpha=0.4,color='grey')

plt.xscale('log')

plt.xlim(0.5,3e2)

plt.xlabel(r'$j_{GMC}\ $ [$pc\ km\ s^{-1}$]')

plt.show()


delthet = fcat['rot_ax_GMC<deg>'] - fcat['rot_ax_HI<deg>']

for i in range(22):
    if delthet[i] > 180:
        delthet[i] -=360


plt.hist(delthet,bins=np.arange(-180,180,30),color='xkcd:light blue',edgecolor='xkcd:blue')

plt.vlines([-90,90],0.,5.,color='red',linestyle='--')
plt.xticks([-180,-90,0,90,180],[-180,-90,0,90,180])

plt.xlabel(r'$\Delta \phi_{rot}\ $ [$\degree$]')
plt.show()


(np.arange(2163)-1082)*0.083


alpha_vir_hi = ((fcat['radius_HI<pc>']*u.pc*(fcat['vel_dispersion_HI<km/s>']*u.km/u.s)**2)/\
(c.G*fcat['mass_HI<kg>']*u.kg)).si


alpha_vir_h2 = ((fcat['radius_GMC<pc>']*u.pc*(fcat['vel_dispersion_GMC<km/s>']*u.km/u.s)**2)/\
(c.G*fcat['mass_GMC<kg>']*u.kg)).si


plt.scatter(alpha_vir_h2,alpha_vir_hi)


# plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.loglog()
plt.show()


plt.scatter(fcat['name'],alpha_vir_hi)
plt.scatter(fcat['name'],alpha_vir_h2)


plt.xticks(rotation='vertical')

plt.hlines([0.1,1,5,10],'aquila_s','taurus',linestyle='--')

# plt.yscale('log')
plt.show()


plt.scatter(fcat['name'],fcat['mass_HI<kg>']/fcat['mass_GMC<kg>'])


plt.xticks(rotation='vertical')

plt.hlines([0.1,1,5,10],'aquila_s','taurus',linestyle='--')

# plt.yscale('log')
plt.show()


# plt.scatter(fcat['name'],alpha_vir_hi)
# plt.scatter(fcat['name'],alpha_vir_h2)

plt.scatter(fcat['name'],alpha_vir_hi/alpha_vir_h2)


plt.xticks(rotation='vertical')

plt.yscale('log')
plt.show()


dust_lb = np.nansum(dust,axis=0)*2.8


plt.figure(figsize=(20,10))

plt.imshow(dust_lb,extent=[180,-180,-90,90],origin='lower',norm='asinh')
plt.contour(dust_lb,extent=[180,-180,-90,90],levels=[0.2],linewidths=0.6,colors='white')
plt.contour(h0,extent=[360,-360,-50,50],levels=[4e20],linewidths=0.6,colors='red')


plt.xlim(180,-180)
plt.ylim(-60,60)


for i in range(22):
    plt.imshow(dust_lb,extent=[180,-180,-90,90],origin='lower',norm='asinh')
    plt.contour(dust_lb,extent=[180,-180,-90,90],levels=[0.2,0.5,1],linewidths=0.6,colors='white')
    plt.contourf(c1s[i],cmap='PiYG',
              extent = [loc[i]*-0.125 + 360, lic[i]*-0.125 + 360, 
                        boc[i]*0.125 - 50, bic[i]*0.125 - 50])
#     plt.colorbar()
    
    plt.xlim((loc[i]-100)*-0.125 + 360, (lic[i]+100)*-0.125 + 360)
    plt.ylim((boc[i]-100)*0.125 - 50, (bic[i]+100)*0.125 - 50)
    plt.show()


plt.rcParams['font.size']=22


plt.figure(figsize=(15,9))

plt.contourf(np.log10(h0),extent=[360,-360,-50,50],
         levels=np.arange(20,21.5,0.1),cmap='YlGnBu_r')

plt.colorbar(label=r'$N_{H}\ $ $[cm^{-2}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])

plt.contourf(np.log10(c0),extent=[360,-360,-50,50],
         levels=np.arange(20.5,22.4,0.2),cmap='magma',linewidths=0.3)

plt.colorbar()#label=r'$N_{H_2}\ $ $[cm^{-2}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.xlim(-50,-73)
plt.ylim(-25,-8)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

# plt.savefig('orionhi',bbox_inches='tight')

plt.title('chamaeleon')


plt.figure(figsize=(16,10))


i = 5
plt.contourf(h2s[i],extent=[loh[i]*-0.125+360,
                             lih[i]*-0.125+360,
                             boh[i]*0.125-50,
                             bih[i]*0.125-50,
                            ],levels=7,cmap='binary_r',alpha=0.7)


plt.colorbar(label=r'$\sigma_{H}\ $ $[km\ s^{-1}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.contourf(c2s[i],extent=[loc[i]*-0.125+360,
                             lic[i]*-0.125+360,
                             boc[i]*0.125-50,
                             bic[i]*0.125-50,
                            ],levels=np.arange(0,1.75,0.25),cmap='viridis',norm='asinh')

plt.colorbar()#,ticks=[20,21],ticklabels=['twenty','twneyon'])

# plt.contourf(np.log10(c0),extent=[360,-360,-50,50],
#              levels=np.arange(20.5,22.4,0.2),cmap='magma',linewidths=0.3)

# plt.colorbar()#label=r'$N_{H_2}\ $ $[cm^{-2}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.xlim(8,-8)
# plt.ylim(-25,-8)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

# # plt.savefig('orionhi',bbox_inches='tight')

plt.title('corona_australis')


fig,ax = plt.subplots(1,1)

i=5

ax.plot(np.nanmean(his[i],axis=(1,2)),color='teal')

ax2 = ax.twinx()

ax2.plot(np.nanmean(cos[i],axis=(1,2)),color='pink')

plt.xlabel(r'$v_{lsr}\ $ [$km\ s^{-1}$]')
plt.xlim(15,50)
ax2.set_ylim(0,0.75)


ax.set_ylabel(r'$T_{B,HI}\ $ [$K$]')
ax2.set_ylabel(r'$T_{B,H_2}\ $ [$K$]')

plt.title('corona_australis')


plt.figure(figsize=(16,10))


i = 12
plt.contourf(h2s[i],extent=[loh[i]*-0.125+360,
                             lih[i]*-0.125+360,
                             boh[i]*0.125-50,
                             bih[i]*0.125-50,
                            ],levels=np.arange(4.5,10,0.5),cmap='binary_r',alpha=0.7)


plt.colorbar(label=r'$\sigma_{H}\ $ $[km\ s^{-1}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.contourf(c2s[i],extent=[loc[i]*-0.125+360,
                             lic[i]*-0.125+360,
                             boc[i]*0.125-50,
                             bic[i]*0.125-50,
                            ],levels=np.arange(0,4,0.5),cmap='viridis',norm='asinh')

plt.colorbar()#,ticks=[20,21],ticklabels=['twenty','twneyon'])

# plt.contourf(np.log10(c0),extent=[360,-360,-50,50],
#              levels=np.arange(20.5,22.4,0.2),cmap='magma',linewidths=0.3)

# plt.colorbar()#label=r'$N_{H_2}\ $ $[cm^{-2}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.xlim(-153,-164)
# plt.ylim(-25,-8)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

# # plt.savefig('orionhi',bbox_inches='tight')

plt.title('mon_ob1')


fcat[12]


fcat[5]


(77.5/745)*u.rad.to(u.deg)/0.125


np.nanmean((-0.5,3.5))


def generate_pars(cat,i):
    d = cat['distance<pc>'][i]
    lo = cat['l_pix_HI'][i][0]
    bo = cat['b_pix_HI'][i][0]
    
    csz = cat['l_pix_HI'][i][1] - cat['l_pix_HI'][i][0]
    
    patchsz = (zcat['radius_acc_HI<pc>'][i]//40 + 1)*10
    
    return {'d':d,'lo':lo,'bo':bo,'cloud_size_px':csz,\
            'patch_size_pc':patchsz}


def generate_patches(cube,params): 
    d,sidelen,lo,bo,ppc = params['d'],params['cloud_size_px'],\
                        params['lo'],params['bo'],params['patch_size_pc']

    pc10 = (ppc/d)*u.rad.to(u.deg)/0.125 
    #npix corresponding to 10pc at cloud distance

    npatch = sidelen//pc10 #how many patches you can make
    extra = sidelen%pc10 #extra pixels
    
    brdrds = [0]+ list(int(extra/2)+np.arange(1,npatch)*int(pc10)) +[sidelen]
    coords = np.array((brdrds[:-1],brdrds[1:])).transpose()

    slices_px = np.array([(i,j) for j in coords for i in coords])

    slices_lb = (slices_px.transpose().reshape(4,-1) + \
                np.broadcast_to(np.array([lo,bo,lo,bo]),
                (int(npatch**2),4)
                ).transpose()
                ).transpose().reshape(-1,2,2).swapaxes(1,2).astype('int')

    axcoords = np.flip(np.array(np.meshgrid(
                        np.arange(npatch),
                        np.arange(npatch))).transpose(),
                axis=0).reshape(-1,2).astype('int')
    
    
    specs = []
    for i in range(int(npatch**2)):
        lo,li = slices_lb[i][0]
        bo,bi = slices_lb[i][1]

        sp = np.nanmean(cube[:,bo:bi,lo:li],axis=(1,2))
        specs.append(sp)
    
    return specs,axcoords,ppc

    
cor_pars = {'d':151,'lo':0.75,'bo':-19.75,'cloud_size_px':65,\
            'patch_size_pc':5}


mob_pars = {'d':745,'lo':-158.5,'bo':1.5,'cloud_size_px':48,\
            'patch_size_pc':20}


specs,axcoords,ppc = fns.generate_patches(his[12],mob_pars)


specs[3]


fig,axs = plt.subplots(4,4,figsize=(30,30),sharex=True,sharey=True)
for j in range(len(axcoords)):
    nr,nc = axcoords[j]
    axs[nr,nc].plot(specs[j])
    
plt.subplots_adjust(wspace=0,hspace=0)


axcoords


plt.figure(figsize=(16,10))


plt.contourf(h2s[18],extent=[loh[18]*-0.125+360,
                             lih[18]*-0.125+360,
                             boh[18]*0.125-50,
                             bih[18]*0.125-50,
                            ],levels=7,cmap='binary_r',alpha=0.7)


plt.colorbar(label=r'$\sigma_{H}\ $ $[km\ s^{-1}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.contourf(c2s[18],extent=[loc[18]*-0.125+360,
                             lic[18]*-0.125+360,
                             boc[18]*0.125-50,
                             bic[18]*0.125-50,
                            ],levels=np.arange(0,4.5,0.5),cmap='viridis',norm='asinh')

plt.colorbar()#,ticks=[20,21],ticklabels=['twenty','twneyon'])

# plt.contourf(np.log10(c0),extent=[360,-360,-50,50],
#              levels=np.arange(20.5,22.4,0.2),cmap='magma',linewidths=0.3)

# plt.colorbar()#label=r'$N_{H_2}\ $ $[cm^{-2}]$')#,ticks=[20,21],ticklabels=['twenty','twneyon'])


plt.xlim(170,148)
# plt.ylim(-25,-8)
plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

# # plt.savefig('orionhi',bbox_inches='tight')

plt.title('perseus')


plt.figure(figsize=(30,15))

plt.contourf(np.log10(hlv),extent=[360,-360,-50,50],
             levels=np.arange(19.8,22,0.2),cmap='YlGnBu_r',alpha=0.5)

plt.contourf(np.log10(clv),extent=[360,-360,-50,50],
             levels=np.arange(20,22,0.2),cmap='magma',linewidths=0.3,alpha=0.5)


plt.plot(np.arange(180,-180,-1),50*np.sin(2*np.arange(180,-180,-1)*u.deg+15*u.deg))
# plt.ylim(-20,40)
plt.xlim(180,-180)

plt.ylabel(r'velocity [$km\ s^{-1}$]')
# plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

# plt.savefig('orionlv',bbox_inches='tight')


plt.figure(figsize=(20,15))

plt.contourf(np.log10(hlv),extent=[360,-360,-50,50],
             levels=np.arange(19.8,22,0.2),cmap='YlGnBu_r')

plt.contourf(np.log10(clv),extent=[360,-360,-50,50],
             levels=np.arange(19.6,22,0.2),cmap='magma',linewidths=0.3)

# plt.ylim(-20,40)
plt.xlim(270,90)

plt.ylabel(r'velocity [$km\ s^{-1}$]')
# plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

# plt.savefig('orionlv',bbox_inches='tight')


plt.figure(figsize=(20,20))

plt.contourf(np.log10(hlv),extent=[360,-360,-50,50],
             levels=np.arange(19.8,22,0.2),cmap='YlGnBu_r')

plt.contourf(np.log10(clv),extent=[360,-360,-50,50],
             levels=np.arange(20,22,0.2),cmap='magma')

plt.ylim(-20,40)
plt.xlim(-130,-170)

plt.ylabel(r'v')
# plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

plt.savefig('orionlv',bbox_inches='tight')


dustAv = dust*2.8


dld = np.nansum(dustAv,axis=1)


dlb = np.nansum(dustAv,axis=0)


plt.figure(figsize=(30,13))
plt.contourf(dlb,extent=[180,-180,-90,90],levels=np.arange(0.15,4,0.5))
plt.contour(dlb,extent=[180,-180,-90,90],levels = [0.2,1],colors='red')
plt.xlim(180,-180)
# plt.xlim(-130,-170)
plt.ylim(-50,50)
# plt.ylim(-40,0)


plt.figure(figsize=(20,10))
plt.contourf(dld,extent=[180,-180,68,1200],levels=np.arange(0.15,2,0.1))

plt.xlim(-130,-155)
plt.ylim(300,800)


i=1
l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50

l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

plt.figure(figsize=(15,10))
plt.contourf(h1s[i],extent=[l0h,l1h,b0h,b1h],cmap='RdGy',
               levels=np.arange(np.nanmedian(c1s[i])-4*np.nanstd(h1s[i]),
                np.nanmedian(c1s[i])+4*np.nanstd(h1s[i]),0.25*np.nanstd(h1s[i])))

plt.colorbar(label=r'$v_{lsr}\ [km\ s^{-1}]$',shrink=0.8,
             ticks=[0.8,2.6,4.4,5.6,6.8,8.6,10.4])


plt.contourf(c1s[i],extent=[l0c,l1c,b0c,b1c],cmap='PiYG',
                               levels=np.arange(np.nanmedian(c1s[i])-3*np.nanstd(c1s[i]),
                               np.nanmedian(c1s[i])+3*np.nanstd(c1s[i]),0.25*np.nanstd(c1s[i])))

plt.colorbar(shrink=0.8,ticks=[4.4,5.0,5.6,6.2,6.8]) #label=r'$v_{lsr,H_2}\ [km\ s^{-1}]$'

plt.contour(c0s[i],extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
                    colors='black',linewidths=0.4)
plt.xlim(l0h,l1h)


plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

plt.title(fcat['name'][i])

# plt.plot(np.arange(170,150,-0.125),np.arange(170,150,-0.125)*\
#          np.tan((-163.1979151135055)*u.deg)-60)


# plt.plot(np.arange(166.5,155.0,-0.125),np.arange(166.5,155.0,-0.125)*\
#          np.tan((-51.95768563566319)*u.deg)+198)

# plt.savefig(fcat['name'][i]+'_mom1')

plt.show()


mjfit.cov_matrix


mjfit


dir(mjfit)


fcat['SFR<Msun/yr>']


fcat


tdepjj = fttr(init,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']),
             np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']/fcat['spec_ang_mom_HI<pc.km/s>']))


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9,
           fcat['spec_ang_mom_GMC<pc.km/s>']/fcat['spec_ang_mom_HI<pc.km/s>'],
           s=20,color='black')

plt.loglog()

# plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9,
#          (fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9)**0.3)


plt.xlabel('$t_{dep}\ $ [$Gyr$]')
plt.ylabel(r'$j_{GMC}\ /\ j_{HI}$' )

plt.show()


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9,
           fcat['mass_GMC<kg>']/fcat['mass_HI<kg>'],
           s=20,color='black')

plt.loglog()

# plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9,
#          (fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9)**0.3)


plt.xlabel('$t_{dep}\ $ [$Gyr$]')
plt.ylabel(r'$M_{GMC}\ /M_{HI}$' )

plt.show()


(fcat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc).si


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>'],
           (2*np.pi/(fcat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)).to(u.yr),
           s=20,color='black')

plt.loglog()

# plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9,
#          (fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9)**0.3)


plt.xlabel('$t_{dep}\ $ [$yr$]')
plt.ylabel(r'$t_{rot}\ $ [$yr$]' )

plt.yticks([1e7,1e8,1e9],[r'$10^7$',r'$10^8$',r'$10^9$'])

plt.show()


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
           (2*np.pi/(fcat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)).to(u.yr),
           s=20,color='black')

plt.loglog()

# plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9,
#          (fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/fcat['SFR<Msun/yr>']*1e-9)**0.3)


plt.xlabel('$M_{GMC}\ [M_\odot]$')
plt.ylabel(r'$t_{rot}\ $ [$Gyr$]' )

plt.yticks([1e7,1e8,1e9],[r'$10^7$',r'$10^8$',r'$10^9$'])

plt.show()


for i in range(22):
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i])
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i])

    
    lc,bc = np.nanmean(fcat['lon_GMC<deg>'][i]),np.nanmean(fcat['lat_GMC<deg>'][i])
    
    hthet = fcat['rot_ax_HI<deg>'][i]*u.deg
    cthet = fcat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
    
#     plt.contourf(bperp.reshape(h1s[i].shape).value,
#                  extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
#     plt.colorbar()
#     plt.scatter(lons,rot_ax_line)
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.show()
    
    
    hv_model = np.nanmean(hmom1) + \
                fcat['vel_grad_HI<km/s/pc>'][i]*np.cos(hthet+90*u.deg)*\
                (hlons-lc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] + \
                fcat['vel_grad_HI<km/s/pc>'][i]*np.sin(hthet+90*u.deg)*\
                (hlats-bc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] 
    
    cv_model = np.nanmean(cmom1) + \
                fcat['vel_grad_GMC<km/s/pc>'][i]*np.cos(cthet+90*u.deg)*\
                (clons-lc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] + \
                fcat['vel_grad_GMC<km/s/pc>'][i]*np.sin(cthet+90*u.deg)*\
                (clats-bc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] 
    
    
    hbperp_pc = hbperp*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    cbperp_pc = cbperp*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    
    fig,axs = plt.subplots(1,2,figsize=(10,5))
    fig.suptitle(fcat['name'][i])
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(c0s[i]),
            cmom1,
            cbperp_pc,
            cv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
    
    axs[0].scatter(cdf['r'],cdf['v'],s=500/cbperp.size,color='black')
    
    axs[0].plot(cdf['r'],cdf['vmod'],color='red')
    
    #not sure if these are the right error bars ... 
    axs[0].plot(cdf['r'],cdf['vmod']-np.nanstd(cdf['vmod']),color='pink') 
    axs[0].plot(cdf['r'],cdf['vmod']+np.nanstd(cdf['vmod']),color='pink')
    
    
    axs[1].scatter(hbperp_pc,hmom1,s=100/hbperp.size,color='black')
    
    axs[1].plot(hbperp_pc,hv_model,color='red')
    
    #not sure if these are the right error bars ... 
    axs[1].plot(hbperp_pc,hv_model-np.nanstd(hv_model),color='pink') 
    axs[1].plot(hbperp_pc,hv_model+np.nanstd(hv_model),color='pink')
    
    
    axs[1].plot(hbperp_pc,hv_model-np.nanstd(hmom1),color='orange') 
    axs[1].plot(hbperp_pc,hv_model+np.nanstd(hmom1),color='orange')
    
    
    for ax in axs.flat:
        ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
        
        
    axs[0].set_ylabel(r'velocity [$km\ s^{-1}$]')
    
#     axs[0].set_ylim(np.nanmean(cdf['vmod'])-5*np.nanstd(cdf['vmod']),
#                     np.nanmean(cdf['vmod'])+5*np.nanstd(cdf['vmod']))
    
#     axs[1].set_ylim(np.nanmean(hv_model)-5*np.nanstd(hv_model),
#                     np.nanmean(hv_model)+5*np.nanstd(hv_model))
    
    
#     plt.savefig(fcat['name'][i]+'vel_grad_fit',bbox_inches='tight')
    
#     print(np.nanstd(hv_model),np.nanstd(hmom1))
    plt.show()
    
    
rsq_co = np.zeros(22)
rsq_hi = np.zeros(22)


for i in range(22):
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i])
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i])

    
    lc,bc = np.nanmean(fcat['lon_GMC<deg>'][i]),np.nanmean(fcat['lat_GMC<deg>'][i])
    
    hthet = fcat['rot_ax_HI<deg>'][i]*u.deg
    cthet = fcat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
        
    hv_model = np.nanmean(hmom1) + \
                fcat['vel_grad_HI<km/s/pc>'][i]*np.cos(hthet+90*u.deg)*\
                (hlons-lc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] + \
                fcat['vel_grad_HI<km/s/pc>'][i]*np.sin(hthet+90*u.deg)*\
                (hlats-bc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] 
    
    cv_model = np.nanmean(cmom1) + \
                fcat['vel_grad_GMC<km/s/pc>'][i]*np.cos(cthet+90*u.deg)*\
                (clons-lc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] + \
                fcat['vel_grad_GMC<km/s/pc>'][i]*np.sin(cthet+90*u.deg)*\
                (clats-bc)*u.deg.to(u.rad)*fcat['distance<pc>'][i] 
    
    
    hbperp_pc = hbperp*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    cbperp_pc = cbperp*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(c0s[i]),
            cmom1,
            cbperp_pc,
            cv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
   
#     plt.scatter(np.nanstd(hv_model),np.nanstd(hmom1))
#     plt.scatter(hbperp_pc,(hv_model-hmom1),s=10)
#     plt.title(fcat['name'][i])
# #     plt.ylim(-5,5)
#     plt.show()
    
#     plt.scatter(np.nanstd((hv_model-hmom1)),np.nanstd((hv_model)-np.nanstd(hmom1)))
    
    rsq_co[i] = (np.corrcoef(cdf['r'],cdf['v'])[0][1])**2
    rsq_hi[i] = (np.corrcoef(hbperp_pc,hmom1)[0][1])**2
    
    #how well the data is fit by a plane
    
    
#     print(fcat['name'][i],np.corrcoef(cdf['r'],cdf['v'])[0][1],
#           np.corrcoef(hbperp_pc,hmom1)[0][1],
#          '\n')
    
#     print()
    
    
dust720 = np.concatenate((dust_lb,dust_lb),axis=1)


180-(-540)


plt.rcParams['font.size']=30


plt.figure(figsize=(31,20))
plt.contourf(np.log10(h0),extent=[360,-360,-50,50],levels=np.arange(20,22,0.1),
            cmap='YlGnBu_r')
plt.contourf(np.log10(c0),extent=[360,-360,-50,50],levels=np.arange(20,23,0.1),
            cmap='magma',linewidths=0.7)
# plt.contour(dust720,extent=[180,-540,-90,90])


plt.xlim(-136,-212)
plt.ylim(-42,5)

plt.xticks([-140,-150,-160,-170,-180,-190,-200,-210],
           np.array([-140,-150,-160,-170,-180,-190,-200,-210])+360)


plt.xlabel(r'galactic longitude [$\degree$]')
plt.ylabel(r'galactic latitude [$\degree$]')


c.m_p/(((4/3)*np.pi*(0.1*1e-6*u.m)**3)*(2.5*u.g*u.cm**(-3))).si


1.3*u.arcsecond.to(u.rad)*20e6*u.pc


0.083*u.deg.to(u.rad)*50000*u.pc


fcat['mass_GMC<kg>']*u.kg.to(u.solMass)*1e-4


0.5*317*0.3076


h1s[21].shape


fcat['distance<pc>']


fcat['radius_HI<pc>']


np.array((fcat['name'],rsq_hi)).transpose()


np.array((fcat['name'],rsq_co)).transpose()


plt.plot(fcat['name'],rsq_co,marker='o')
plt.plot(fcat['name'],rsq_hi,marker='+')


plt.xticks(rotation='vertical')

# plt.yscale('log')

plt.hlines([0.25,0.5],'aquila_s','taurus')

plt.show()


plt.scatter(fcat['vel_grad_HI<km/s/pc>'],rsq_hi)

plt.xticks(rotation='vertical')

# plt.yscale('log')
# plt.loglog()

plt.show()


np.array((fcat['name'],rsq_co,fcat['vel_grad_GMC<km/s/pc>'])).transpose()


plt.scatter(fcat['vel_grad_GMC<km/s/pc>'],rsq_co)

plt.xticks(rotation='vertical')

# plt.yscale('log')
# plt.loglog()

plt.show()


plt.scatter(delthet,rsq_co)


plt.hist(delthet[rsq_co>0.25],bins=np.arange(-180,210,30),alpha=0.5)
plt.hist(delthet[rsq_co<=0.25],bins=np.arange(-180,210,30),alpha=0.5)


plt.show()


plt.hist(fcat['vel_grad_GMC<km/s/pc>'][rsq_co>0.25],alpha=0.5,
        bins=np.arange(0,0.5,0.05))
plt.hist(fcat['vel_grad_GMC<km/s/pc>'][rsq_co<=0.25],alpha=0.5,
        bins=np.arange(0,0.5,0.05))


plt.show()


plt.hist(fcat['spec_ang_mom_GMC<pc.km/s>'][rsq_co>0.25],alpha=0.5,)
#         bins=np.arange(0,0.5,0.05))
plt.hist(fcat['spec_ang_mom_GMC<pc.km/s>'][rsq_co<=0.25],alpha=0.5,)
#         bins=np.arange(0,0.5,0.05))


plt.show()


plt.hist(fcat['vel_grad_HI<km/s/pc>'][rsq_hi>0.5],alpha=0.5,
        bins=np.arange(0,0.3,0.02))
plt.hist(fcat['vel_grad_HI<km/s/pc>'][rsq_hi<=0.5],alpha=0.5,
        bins=np.arange(0,0.3,0.02))


plt.show()


plt.scatter(fcat['name'],fcat['vel_grad_GMC<km/s/pc>'],
           c=rsq_co)

plt.xticks(rotation='vertical')

plt.yscale('log')

plt.show()


np.corrcoef(np.arange(0,3,0.1),-np.arange(0,30,1))
# plt.loglog()


fcat['name'][18]


perh1 = h1s[18]


import scipy.ndimage as ndi


dir(ndi)


plt.imshow(ndi.gaussian_filter(perh1,sigma=10),
                        extent=[loh[18]*-0.125+360,lih[18]*-0.125+360,
                          boh[18]*0.125-50,bih[18]*0.125-50],origin='lower',
            cmap='RdGy')


dvdyp,dvdxp = np.gradient(ndi.gaussian_filter(perh1,sigma=10))


plt.contourf(dvdxp,extent=[loh[18]*-0.125+360,lih[18]*-0.125+360,
                          boh[18]*0.125-50,bih[18]*0.125-50],
             cmap='RdGy',)

plt.colorbar()


plt.contourf(dvdyp,extent=[loh[18]*-0.125+360,lih[18]*-0.125+360,
                          boh[18]*0.125-50,bih[18]*0.125-50],
             cmap='RdGy',)

plt.colorbar()


plt.contourf(np.sqrt(dvdyp**2 + dvdxp**2),
             extent=[loh[18]*-0.125+360,lih[18]*-0.125+360,
                          boh[18]*0.125-50,bih[18]*0.125-50],
             cmap='RdGy',)

plt.colorbar()


plt.contourf(dvdyp + dvdxp,
             extent=[loh[18]*-0.125+360,lih[18]*-0.125+360,
                          boh[18]*0.125-50,bih[18]*0.125-50],
             cmap='RdGy',)

plt.colorbar()


plt.contourf(perh1,extent=[loh[18]*-0.125+360,lih[18]*-0.125+360,
                          boh[18]*0.125-50,bih[18]*0.125-50],levels=20,
            cmap='RdGy')


h1gal = fns.get_mom1(h3sig,np.arange(155)*1.3-100)


plt.figure(figsize=(20,10))
plt.contourf(h1gal,extent=[360,-360,-50,50],levels=20,cmap='RdGy')


plt.xlim(190,-190)


dvdy,dvdx = np.gradient(ndi.gaussian_filter(h1gal,sigma=20))


plt.figure(figsize=(20,10))
plt.contourf(dvdx,extent=[360,-360,-50,50],levels=np.arange(-0.4,0.4,0.1),cmap='RdGy')


plt.xlim(190,-190)

plt.colorbar()


plt.figure(figsize=(20,10))
plt.contourf(dvdy,extent=[360,-360,-50,50],levels=np.arange(-0.6,0.6,0.1),cmap='RdGy')


plt.xlim(190,-190)


plt.figure(figsize=(20,10))
plt.contourf(dvdy+dvdx,extent=[360,-360,-50,50],levels=np.arange(-0.4,0.4,0.1),cmap='RdGy')


plt.xlim(190,-190)


plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(dvdy**2 +dvdx**2),
             extent=[360,-360,-50,50],levels=np.arange(0,0.6,0.05),cmap='RdGy')

plt.contour(c0,extent=[360,-360,-50,50],levels=[1e20],color='blue')

plt.xlim(190,-190)

# plt.colorbar()


plt.contour(c0,extent=[360,-360,-50,50],levels=[1e20],color='blue')
plt.xlim(360,0)


for i in range(22):
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

    plt.contourf(h1s[i],extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
    
    plt.title(fcat['name'][i])
    
    plt.show()


aspect_ratio = np.nanstd(fcat['lon_GMC<deg>'],axis=1)/np.nanstd(fcat['lat_GMC<deg>'],axis=1)
for i in range(22):
    if aspect_ratio[i] < 1:
        aspect_ratio[i] = 1/aspect_ratio[i] 


plt.scatter(aspect_ratio,fcat['vel_grad_GMC<km/s/pc>'])
plt.loglog()
plt.xticks([1,2,3],[1,2,3])
plt.yticks([1e-2,1e-1,1e0],[r'$10^{-2}$',r'$10^{-1}$',r'1'])

plt.xlabel('aspect ratio')
plt.ylabel(r'$\Omega$ [$km\ s^{-1}\ pc^{-1}$]')
plt.show()


plt.scatter(aspect_ratio,fcat['mass_GMC<kg>']*u.kg.to(u.solMass))
plt.loglog()
plt.xticks([1,2,3],[1,2,3])
# plt.yticks([1e-2,1e-1,1e0],[r'$10^{-2}$',r'$10^{-1}$',r'1'])

plt.xlabel('aspect ratio')
plt.ylabel(r'M [$M_\odot$]')
plt.show()


plt.scatter(aspect_ratio,fcat['spec_ang_mom_GMC<pc.km/s>'])
plt.loglog()
plt.xticks([1,2,3],[1,2,3])

plt.xlabel('aspect ratio')
plt.ylabel(r'j [$pc\ km\ s^{-1}$]')
plt.show()


virial_alp = ((((u.km/u.s*fcat['vel_dispersion_GMC<km/s>'])**2)*(u.pc*fcat['radius_GMC<pc>']))/\
                (c.G*fcat['mass_GMC<kg>']*u.kg)).si


fcat


virial_alp


plt.scatter(aspect_ratio,virial_alp)
# plt.loglog()
plt.xticks([1,2,3],[1,2,3])

plt.xlabel('aspect ratio')
# plt.ylabel(r'M [$M_\odot$]')
plt.show()


plt.scatter(virial_alp,
        fcat['spec_ang_mom_HI<pc.km/s>']/fcat['spec_ang_mom_GMC<pc.km/s>'])
plt.loglog()

plt.show()
