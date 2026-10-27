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

hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data


fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


h0 = np.nansum(h3sig,axis=0)*1.3*1.82e18
c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20


h035 = np.nansum(h3sig[50:105],axis=0)*1.3*1.82e18


# plt.rcParams['font.size'] = 20
# plt.rcParams['font.family'] = 'STIXGeneral'
# plt.rcParams['mathtext.fontset'] = 'stix'

plt.figure(figsize=(30,16),dpi=200)
plt.contourf(np.log10(h0),extent=[360,-360,-50,50],cmap='YlGnBu_r',levels=np.arange(19,22.4,0.2))
plt.colorbar(label=r'log ($N_{HI}$ [$cm^{-2}$])')
plt.contour(np.log10(h0),extent=[360,-360,-50,50],colors='white',levels=[np.log10(4e20)],linewidths=0.5)
# plt.contour(np.log10(h035),extent=[360,-360,-50,50],colors='orange',levels=[np.log10(4e20)],linewidths=0.5)

plt.contourf(np.log10(c0),extent=[360,-360,-50,50],cmap='magma',levels=np.arange(20,23.2,0.2))
plt.colorbar(orientation='horizontal',label=r'log ($N_{H_2}$ [$cm^{-2}$])')
plt.xlim(180,-180)

plt.xlabel(r'galactic longitude [$\degree$]')
plt.ylabel(r'galactic latitude [$\degree$]')

# plt.add_patch(Circle((0,0),10))

# plt.savefig('integrated_intensity_gal')

plt.show()
# plt.ylim


#figure out how to line up colorbars nicer


hlv = np.nansum(h3sig,axis=1)*0.125*1.82e18
clv = np.nansum(d3sig,axis=1)*0.125*1.97e20


#can you convert K-arcdeg to cm^-2 the same way as K-km/s
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'

plt.figure(figsize=(24,10))
plt.contourf(np.log10(hlv),extent=[360,-360,-100,100],cmap='YlGnBu_r',levels=np.arange(20,22.2,0.2))
plt.colorbar(label=r'log ($N_{HI}$ [$cm^{-2}$])')

# plt.contour(np.log10(hlv),extent=[360,-360,-100,100],colors='white',levels=[np.log10(4e20)])

plt.contourf(np.log10(clv),extent=[360,-360,-100,100],cmap='magma',linewidths=0.5,
             levels=np.arange(20.2,22.2,0.2))

plt.colorbar(orientation='horizontal',label=r'log ($N_{H_2}$ [$cm^{-2}$])')

plt.xlim(180,-180)

plt.xlabel(r'galactic longitude [$\degree$]')
plt.ylabel(r'velocity [$km\ s^{-1}$]')

# plt.savefig('lv_gal')


plt.show()


hbv = np.nansum(h3sig,axis=2)*0.125*1.82e18*0.5
cbv = np.nansum(d3sig,axis=2)*0.125*1.97e20*0.5


plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'

plt.figure(figsize=(24,16))
plt.contourf(np.log10(hbv),extent=[-50,50,-100,100],cmap='YlGnBu_r',levels=np.arange(19,23.2,0.2))
# plt.contour(np.log10(hlv),extent=[360,-360,-100,100],colors='white',levels=[np.log10(4e20)])
plt.colorbar(label=r'log ($N_{HI}$ [$cm^{-2}$])')

plt.contourf(np.log10(cbv),extent=[-50,50,-100,100],cmap='magma',linewidths=0.5,
             levels=np.arange(20,23.2,0.2))
plt.colorbar(orientation='horizontal',label=r'log ($N_{H_2}$ [$cm^{-2}$])')


# plt.xlim(180,-180)

plt.xlabel(r'galactic latitude [$\degree$]')
plt.ylabel(r'velocity [$km\ s^{-1}$]')

# plt.savefig('bv_gal')

plt.show()


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


# vo,vi = np.array([((newvrange[i]+100)/1.3).astype('int') for i in range(22)]).transpose()


# vs = [np.arange((vi-vo)[i])*1.3 + newvrange[i][0] for i in range(22)]

# his = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
# cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


#hmmmm ^ about this


# fig,ax = plt.subplots(1,1,figsize=(20,7))

for i in range(22):
    fig,axs = plt.subplots(2,2,figsize=(20,20))

    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    h2 = np.sqrt(fns.get_momN(his[i],vs[i],2))
    
    c1 = fns.get_mom1(cos[i],vs[i])
    c2 = np.sqrt(fns.get_momN(cos[i],vs[i],2))

    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    sph = np.nanmean(his[i],axis=(1,2))
    spc = np.nanmean(cos[i],axis=(1,2))
    
    
    h1im = axs[0,0].contourf(h1,extent=[l0h,l1h,b0h,b1h],cmap='RdGy',
                           levels=np.arange(np.nanmedian(c1)-5*np.nanstd(h1),
                                            np.nanmedian(c1)+5*np.nanstd(h1),0.5*np.nanstd(h1)))
    c1im = axs[0,0].contourf(c1,extent=[l0c,l1c,b0c,b1c],cmap='PiYG',
                             levels=np.arange(np.nanmedian(c1)-3*np.nanstd(c1),
                                            np.nanmedian(c1)+3*np.nanstd(c1),0.5*np.nanstd(c1)))
    
    axs[0,0].contour(c0,extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
                    colors='black',linewidths=0.4)
    
    fig.colorbar(h1im,label=r'$\overline{v_{_{HI}}}\ $  [$km\ s^{-1}$]')
    fig.colorbar(c1im,orientation='horizontal',label=r'$\overline{v_{_{H_2}}}\ $ [$km\ s^{-1}$]')
    
    
    h2im = axs[1,0].contourf(h2,extent=[l0h,l1h,b0h,b1h],cmap='YlGnBu_r',levels=6)
    
    c2im = axs[1,0].contourf(c2,extent=[l0c,l1c,b0c,b1c],cmap='magma',levels=6)
    
    axs[1,0].contour(c0,extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
                    colors='black',linewidths=0.4)    
    fig.colorbar(h2im,label=r'$\sigma_{v,_{HI}}\ $  [$km\ s^{-1}$]')
    fig.colorbar(c2im,orientation='horizontal',label=r'$\sigma_{v,_{H_2}}\ $  [$km\ s^{-1}$]')
    
    
    hspec = axs[0,1].plot(vs[i],sph,color='black')
    
    axs01 = axs[0,1].twinx()
    cspec = axs01.plot(vs[i],spc,color='xkcd:kelly green')
    axs[0,1].set_xlabel(r'$v_{_{cloud}}\ $  [$km\ s^{-1}$]')
    axs[0,1].set_ylabel(r'$T_{B,HI}\ $  [$K$]')
    axs01.set_ylabel(r'$T_{B,H_2}\ $  [$K$]')
    
    
#     print(np.nanmean(c2)+np.nanstd(c2),np.nanmax(c2))
    for ax in [axs[0,0],axs[1,0]]:
        ax.set_xlim(l0h,l1h)
        ax.set_ylim(b0h,b1h)
        ax.set_ylabel(r'galactic latitude [$\degree$]')
        ax.set_xlabel(r'galactic longitude [$\degree$]')

    
    fig.savefig(fcat['name'][i]+'_mom1_mom2_meh')
    
    
    plt.show()
    

#^ add patch by patch spectra in axs[1,1]
#^ add in rotation axis / angular momentum vector to axs[0,0]


fcat['clustered'] = [0,1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1]


cluster_colors = ['black','orange']


fcat[fcat['clustered']==0]


cl_names_0 = ['easter purple','purpleish pink','magenta',
            'red','salmon','pumpkin orange','amber','light mustard','lemon']


cl_names_1 = ['cool blue','topaz','blue blue','true green','blue blue',
              'topaz','true green','foam green','foam green','foam green',
              'cool blue','blue blue','cool blue']


all_names = ['easter purple','cool blue','purpleish pink','magenta','red','salmon','topaz',
'blue blue','pumpkin orange','amber','true green','light mustard','blue blue','topaz',
'true green','foam green','foam green','foam green','cool blue','lemon','blue blue','cool blue']


#come back and save all of these vv with fit lines
#and also figure out how to make the marker correspond to the name
#a la larson, so we can track which clouds are doing what

#will have to discuss gemini being so massive -> far away + large on sky area -> maybe its a bunch of stuff
#blended along the line of sight?? esp cos so close to b=0


plt.figure(figsize=(15,7))

for i in range(22):
    plt.scatter(np.nanmean(fcat['lon_GMC<deg>'][i]),
             np.nanmean(fcat['lat_GMC<deg>'][i]),
            c=cluster_colors[fcat['clustered'][i]],
                s=1e3*fcat['vel_grad_GMC<km/s/pc>'][i],alpha=0.5)
    
    plt.scatter(np.nanmean(fcat['lon_GMC<deg>'][i]),
             np.nanmean(fcat['lat_GMC<deg>'][i]),
            edgecolor=cluster_colors[fcat['clustered'][i]],
            s=1e3*fcat['vel_grad_HI<km/s/pc>'][i],
               facecolor='none')
    
plt.xlim(180,-180)

plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(np.nanmean(fcat['lat_GMC<deg>'][i]),
             np.nanmean(fcat['vel<km/s>'][i]),
            c=cluster_colors[fcat['clustered'][i]],
                s=1e3*fcat['vel_grad_GMC<km/s/pc>'][i])
    
    
plt.xlabel(r'galactic latitude[$\degree$]')
plt.ylabel(r'$v_{lsr}$ [$km\ s^{-1}$]')


plt.xlim(30,-30)


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(np.nanmean(fcat['lon_GMC<deg>'][i]),
             np.nanmean(fcat['vel<km/s>'][i]),
            c=cluster_colors[fcat['clustered'][i]],
        s=1e3*fcat['vel_grad_GMC<km/s/pc>'][i])
    
plt.xlabel(r'galactic longitude [$\degree$]')
plt.ylabel(r'$v_{lsr}$ [$km\ s^{-1}$]')


plt.xlim(180,-180)


fcat


plt.figure(figsize=(7,7))

bins = np.arange(0,0.5,5e-2)

plt.hist(fcat[fcat['clustered']==0]['vel_grad_GMC<km/s/pc>'],density=True,bins=bins,
         color='black',alpha=0.6)

plt.hist(fcat[fcat['clustered']==0]['vel_grad_HI<km/s/pc>'],density=True,bins=bins,
         edgecolor='black',fill=False)

plt.hist(fcat[fcat['clustered']==1]['vel_grad_GMC<km/s/pc>'],density=True,bins=bins,
         color='orange',alpha=0.6)

plt.hist(fcat[fcat['clustered']==1]['vel_grad_HI<km/s/pc>'],density=True,bins=bins,
         edgecolor='orange',fill=False)


plt.xlabel(r'$\Omega\ $ [$km\ s^{-1}\ pc^{-1}$]')

plt.savefig('Omegahist')
plt.show()


plt.figure(figsize=(7,7))

bins = np.arange(0,20,0.5)

plt.hist(fcat[fcat['clustered']==0]['vel_dispersion_GMC<km/s>'],density=True,bins=bins,
         color='black',alpha=0.7)

plt.hist(fcat[fcat['clustered']==0]['vel_dispersion_HI<km/s>'],density=True,bins=bins,
         edgecolor='black',alpha=0.7,fill=False)

plt.hist(fcat[fcat['clustered']==1]['vel_dispersion_GMC<km/s>'],density=True,bins=bins,
         color='orange',alpha=0.7)

plt.hist(fcat[fcat['clustered']==1]['vel_dispersion_HI<km/s>'],density=True,bins=bins,
         edgecolor='orange',alpha=0.7,fill=False)


plt.xlabel(r'$\sigma_v\ $ [$km\ s^{-1}$]')

plt.savefig('sigmahist')
plt.show()


plt.figure(figsize=(7,7))

bins = np.arange(0,100,2.5)


plt.hist(fcat[fcat['clustered']==0]['radius_GMC<pc>'],bins=bins,density=True,
         color='black',alpha=0.7)
plt.hist(fcat[fcat['clustered']==1]['radius_GMC<pc>'],bins=bins,density=True,
         color='orange',alpha=0.7)


plt.hist(fcat[fcat['clustered']==0]['radius_HI<pc>'],bins=bins,density=True,
         edgecolor='black',alpha=0.7,fill=False)
plt.hist(fcat[fcat['clustered']==1]['radius_HI<pc>'],bins=bins,density=True,
         edgecolor='orange',alpha=0.7,fill=False)

plt.xlabel(r'$R_{GMC}\ $[$pc$]')
# plt.xlim(0,20)
plt.savefig('rhist')
plt.show()


plt.figure(figsize=(7,7))

bins = np.arange(0,100,2.5)

plt.hist(fcat[fcat['clustered']==0]['mass_GMC<kg>']*u.kg.to(u.solMass)*1e-4,bins=bins,density=True,
         color='black',alpha=0.7)
plt.hist(fcat[fcat['clustered']==1]['mass_GMC<kg>']*u.kg.to(u.solMass)*1e-4,bins=bins,density=True,
         color='orange',alpha=0.7)

plt.hist(fcat[fcat['clustered']==0]['mass_HI<kg>']*u.kg.to(u.solMass)*1e-4,bins=bins,density=True,
         edgecolor='black',alpha=0.7,fill=False)
plt.hist(fcat[fcat['clustered']==1]['mass_HI<kg>']*u.kg.to(u.solMass)*1e-4,bins=bins,density=True,
         edgecolor='orange',alpha=0.7,fill=False)


plt.xlabel(r'$M\ $ [$10^4 \times M_\odot$]')
# plt.xlim(0,20)
plt.savefig('masshist')
# plt.xscale('log')
plt.show()


fcat['vel_dispersion_GMC<km/s>']/fcat['vel_dispersion_HI<km/s>']


fcat['rot_ax_GMC<deg>','rot_ax_HI<deg>']


fcat['rot_ax_GMC<deg>'] - fcat['rot_ax_HI<deg>']


plt.figure(figsize=(7,7))

bins = np.arange(-180,360,30)

theta_diff = fcat['rot_ax_GMC<deg>'] - fcat['rot_ax_HI<deg>']
for i in range(len(theta_diff)):
    if theta_diff[i] > 180:
        theta_diff[i] -=360
        
    elif theta_diff[i] < -180:
        theta_diff[i] +=360


plt.hist(theta_diff[fcat['clustered']==0],bins=bins,density=True,
         color='black',alpha=0.6)
plt.hist(theta_diff[fcat['clustered']==1],bins=bins,density=True,
         color='orange',alpha=0.6)


plt.vlines([-90,90],0,0.015,color='pink',linestyle='--')

plt.xlabel(r'$\Delta \phi_{rot}\ $ [$\degree$]')
plt.xlim(180,-180)
plt.savefig('thetadiffhist')
# plt.xscale('log')
plt.show()


plt.figure(figsize=(7,7))

bins = np.arange(0,80,5)

SigmaG = fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_GMC<pc>']**2)
SigmaH = fcat['mass_HI<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>']**2)


plt.hist(SigmaG[fcat['clustered']==0],bins=bins,density=True,
         color='black',alpha=0.6)
plt.hist(SigmaG[fcat['clustered']==1],bins=bins,density=True,
         color='orange',alpha=0.6)

plt.hist(SigmaH[fcat['clustered']==0],bins=bins,density=True,
         edgecolor='black',fill=False)
plt.hist(SigmaH[fcat['clustered']==1],bins=bins,density=True,
         edgecolor='orange',fill=False)


plt.vlines(10,0,0.1,color='pink',linestyle='--')

plt.xlabel(r'$\Sigma $ [$M_\odot/\ pc^{-2}$]')
# plt.xlim(0,20)
plt.savefig('masshist')
# plt.xscale('log')
plt.show()


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()


sigmarfitG = fttr(init,np.log10(fcat['radius_GMC<pc>']),np.log10(fcat['vel_dispersion_GMC<km/s>']))
sigmarfitH = fttr(init,np.log10(fcat['radius_HI<pc>']),np.log10(fcat['vel_dispersion_HI<km/s>']))


def fitlinelog(x,m,b):
    return (x**m)*(10**b)
#     return 10**(m*np.log10(x) + b)
#     return m*x + b


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['radius_GMC<pc>'][i],
                fcat['vel_dispersion_GMC<km/s>'][i],
                color=cluster_colors[fcat['clustered'][i]],
               alpha=0.7)
    
    
    plt.scatter(fcat['radius_HI<pc>'][i],
                fcat['vel_dispersion_HI<km/s>'][i],
                edgecolor=cluster_colors[fcat['clustered'][i]],
                facecolor='none',
               alpha=0.7)
#     print(fcat['name'][i],1e-3*fcat['mass_GMC<kg>'][i]/(c.M_sun.value))    
# plt.xscale('log')

plt.plot(fcat['radius_GMC<pc>'],
         fitlinelog(fcat['radius_GMC<pc>'],
            sigmarfitG.slope,sigmarfitG.intercept),
        color='red',linewidth=0.3,)

plt.plot(fcat['radius_HI<pc>'],
         fitlinelog(fcat['radius_HI<pc>'],
            sigmarfitH.slope,sigmarfitH.intercept),
        color='pink',linewidth=0.3)


# plt.legend()
plt.loglog()
plt.ylabel(r'$\sigma_v\ $ [$km\ s^{-1}$]')
plt.xlabel(r'$R\ $ [$pc$]')

# plt.savefig('sigmavsr')
plt.show()


sigmamassfitG = fttr(init,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)),
                     np.log10(fcat['vel_dispersion_GMC<km/s>']))
sigmamassfitH = fttr(init,np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)),
                     np.log10(fcat['vel_dispersion_HI<km/s>']))


plt.figure(figsize=(7,7))
for i in range(22):
    plt.scatter(fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                fcat['vel_dispersion_GMC<km/s>'][i],
                color=cluster_colors[fcat['clustered'][i]],
                s=1e3*fcat['vel_grad_GMC<km/s/pc>'][i],
               alpha=0.7)
    
    
    plt.scatter(fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass),
                fcat['vel_dispersion_HI<km/s>'][i],
                facecolor='none',
                edgecolor=cluster_colors[fcat['clustered'][i]],
                s=1e3*fcat['vel_grad_HI<km/s/pc>'][i],
               alpha=0.7)
#     print(fcat['name'][i],1e-3*fcat['mass_GMC<kg>'][i]/(c.M_sun.value))    


plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
         fitlinelog(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
            sigmamassfitG.slope,sigmamassfitG.intercept),
        color='red',linewidth=0.3,)

plt.plot(fcat['mass_HI<kg>']*u.kg.to(u.solMass),
         fitlinelog(fcat['mass_HI<kg>']*u.kg.to(u.solMass),
            sigmamassfitH.slope,sigmamassfitH.intercept),
        color='pink',linewidth=0.3)


plt.loglog()
plt.xlabel(r'$M\ $ [$M_\odot$]')
plt.ylabel(r'$\sigma_v\ $ [$km\ s^{-1}$]')

# plt.savefig('sigmavsmass')
plt.show()


rmassfitG = fttr(init,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)),
                     np.log10(fcat['radius_GMC<pc>']))
rmassfitH = fttr(init,np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)),
                     np.log10(fcat['radius_HI<pc>']))


plt.figure(figsize=(7,7))
for i in range(22):
    plt.scatter(fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                fcat['radius_GMC<pc>'][i],
                color=cluster_colors[fcat['clustered'][i]],
                s=1e3*fcat['vel_grad_GMC<km/s/pc>'][i],
               alpha=0.7)
    
    
    plt.scatter(fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass),
                fcat['radius_HI<pc>'][i],
                facecolor='none',
                edgecolor=cluster_colors[fcat['clustered'][i]],
                s=1e3*fcat['vel_grad_HI<km/s/pc>'][i],
               alpha=0.7)
#     print(fcat['name'][i],1e-3*fcat['mass_GMC<kg>'][i]/(c.M_sun.value))    


plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
         fitlinelog(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
            rmassfitG.slope,rmassfitG.intercept),
        color='red',linewidth=0.3,)

plt.plot(fcat['mass_HI<kg>']*u.kg.to(u.solMass),
         fitlinelog(fcat['mass_HI<kg>']*u.kg.to(u.solMass),
            rmassfitH.slope,rmassfitH.intercept),
        color='pink',linewidth=0.3)


plt.loglog()
plt.xlabel(r'$M\ $ [$M_\odot$]')
plt.ylabel(r'$R\ $ [$pc$]')

# plt.savefig('rvsmass')
plt.show()


sigmaHmassGfit = fttr(init,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)),
                     np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>']**2)))


plt.figure(figsize=(10,7))
for i in range(22):
    plt.scatter(fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>'][i]**2),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e4*fcat['vel_grad_GMC<km/s/pc>'][i],alpha=0.6)
plt.loglog()
plt.xlabel(r'$M_{GMC}\ $ [$M_\odot$]')
plt.ylabel(r'$\Sigma_{HI}\ $ [$M_\odot\ pc^{-2}$]')

plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
         fitlinelog(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
            sigmaHmassGfit.slope,sigmaHmassGfit.intercept),
        color='xkcd:pastel red',linewidth=0.3)

plt.savefig('SigmaHvsmassG')
plt.show()


sigmaHjGfit = fttr(init,np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km)),
                     np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>']**2)))


plt.figure(figsize=(10,7))
for i in range(22):
    plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>'][i]**2),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)


plt.plot(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
         fitlinelog(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
            sigmaHjGfit.slope,sigmaHjGfit.intercept),
        color='xkcd:pastel red',linewidth=0.3)

plt.loglog()
plt.xlabel(r'$j_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$\Sigma_{HI}\ $ [$M_\odot\ pc^{-2}$]')

plt.savefig('SigmaHvsjG')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_grad_GMC<km/s/pc>'][i]*1e2,
                fcat['vel_grad_HI<km/s/pc>'][i]*1e2,
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)
#     print(fcat['name'][i],1e-3*fcat['mass_GMC<kg>'][i]/(c.M_sun.value))    
plt.loglog()
plt.xlabel(r'$\Omega_{GMC}\ $ [$0.01 \times\ km\ s^{-1}\ pc^{-1}$]')
plt.ylabel(r'$\Omega_{HI}\ $ [$0.01 \times\ km\ s^{-1}\ pc^{-1}$]')

plt.savefig('OmHvsOmG')
plt.show()


massHjGfit = fttr(init,np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km)),
                     np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)))


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)

plt.plot(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
         fitlinelog(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
            massHjGfit.slope,massHjGfit.intercept),
        color='xkcd:pastel red',linewidth=0.3)
    
plt.loglog()
plt.xlabel(r'$j_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$M_{HI}\ $ [$M_\odot$]')
# plt.savefig('mHvsjG')
plt.show()


massGjGfit = fttr(init,np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km)),
                     np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)))


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e4*fcat['vel_grad_GMC<km/s/pc>'][i],alpha=0.6)

plt.plot(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
         fitlinelog(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
            massGjGfit.slope,massGjGfit.intercept),
        color='xkcd:pastel red',linewidth=0.3)
    
plt.loglog()
plt.xlabel(r'$j_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$M_{GMC}\ $ [$M_\odot$]')
# plt.savefig('mGvsjG')
plt.show()


jHjGfit = fttr(init,np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km)),
                    np.log10(fcat['spec_ang_mom_HI<pc.km/s>']*u.pc.to(u.km)))


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['spec_ang_mom_HI<pc.km/s>'][i]*u.pc.to(u.km),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)

    
plt.plot(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
         fitlinelog(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
            jHjGfit.slope,jHjGfit.intercept),
        color='xkcd:pastel red',linewidth=0.3)

plt.loglog()
plt.xlabel(r'$j_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$j_{HI}\ $ [$km^{2}\ s^{-1}$]')
# plt.savefig('jHvsjG')
plt.show()


jHperjGjGfit = fttr(init,np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km)),
                    np.log10(fcat['spec_ang_mom_HI<pc.km/s>']*u.pc.to(u.km)/
                (fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km))))


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['spec_ang_mom_HI<pc.km/s>'][i]*u.pc.to(u.km)/
                (fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km)),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)


plt.plot(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
         fitlinelog(fcat['spec_ang_mom_GMC<pc.km/s>']*u.pc.to(u.km),
            jHperjGjGfit.slope,jHperjGjGfit.intercept),
        color='xkcd:pastel red',linewidth=0.3)   
    
plt.loglog()
plt.xlabel(r'$j_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$j_{HI}/j_{GMC}\ $')
plt.savefig('jHperjGvsjG')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_dispersion_GMC<km/s>'][i],
#         fcat['radius_GMC<pc>'][i],
#          fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['vel_dispersion_HI<km/s>'][i]/
                (fcat['vel_dispersion_GMC<km/s>'][i]),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$\sigma_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$\sigma_{HI}/\sigma_{GMC}\ $')
plt.savefig('sigmaHpersigmaGvsjG')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'][i]*u.pc.to(u.km),
                fcat['vel_dispersion_HI<km/s>'][i]/
                (fcat['vel_dispersion_GMC<km/s>'][i]),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$j_{GMC}\ $ [$km^{2}\ s^{-1}$]')
plt.ylabel(r'$\sigma_{HI}/\sigma_{GMC}\ $')
# plt.savefig('sigmaHpersigmaGvsjG')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                fcat['vel_dispersion_HI<km/s>'][i]/
                (fcat['vel_dispersion_GMC<km/s>'][i]),
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$M_{GMC}\ $ [$M_\odot$]')
plt.ylabel(r'$\sigma_{HI}/\sigma_{GMC}\ $')
# plt.savefig('sigmaHpersigmaGvsjG')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                fcat['vel_dispersion_HI<km/s>'][i],
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-2*fcat['mass_GMC<kg>'][i]/(c.M_sun.value),alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$M_{GMC}\ $ [$M_\odot$]')
plt.ylabel(r'$\sigma_{HI}$ [$km\ s^{-1}$]')
# plt.savefig('sigmaHpersigmaGvsjG')
plt.show()


#calculate fitlines for last couple + the goodness-of-fit^ for all
#also more things plotting H vs G, because the whole point is to see how 
#diff atomic enviroments scale with different molecular environments


#then do scatter plots of sigma gmc vs sigma h for clouds that are within the co bounds & without
#+ plots from the list, and more plots from other papers
#and more comparative hi <-> gmc


np.array((fcat['name'],np.zeros(22))).transpose()


fcat


ysoN_names = np.array([['aquila_s', 62],
['california', 770],
['camelopardalis', np.nan],
['chamaeleon', 394],
['cma_ob1', 408],
['corona_australis', 28],
['crossbones', 176],
['gemini_ob1', 64],
['hercules',np.nan],
['lacerta', 525],
['lupus', 1059],
['maddalena', 17],
['mon_ob1', 577],
['mon_r2', 473],
['ophiuchus', 791],
['orion_a', 2072],
['orion_b', 854],
['orion_l', 1075],
['perseus', 435],
['polaris', 15],
['rosette', 78],
['taurus', 432]])


ysos = ysoN_names.transpose()[1].astype('float')


sfrs = 0.25*ysos*1e-6#*u.solMass/u.yr


fcat['SFR<Msun/yr>'] = sfrs


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
                fcat['SFR<Msun/yr>'][i],
                color=cluster_colors[fcat['clustered'][i]],
               s=1e4*fcat['vel_grad_GMC<km/s/pc>'][i],alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$M_{GMC}\ $ [$M_\odot$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.savefig('SFRvsmassG')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass),
                fcat['SFR<Msun/yr>'][i],
                color=cluster_colors[fcat['clustered'][i]],
               s=1e4*fcat['vel_grad_GMC<km/s/pc>'][i],alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$M_{HI}\ $ [$M_\odot$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.savefig('SFRvsmassH')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_dispersion_HI<km/s>'][i]/fcat['vel_dispersion_GMC<km/s>'][i],
                fcat['SFR<Msun/yr>'][i],
                color=cluster_colors[fcat['clustered'][i]],
               s=1e-3*fcat['mass_GMC<kg>'][i],alpha=0.6)
    

plt.loglog()
plt.xlabel(r'$M_{GMC}\ $ [$M_\odot$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.savefig('SFRvssigmaHpersigmaG')
plt.show()


plt.scatter(fcat['vel_dispersion_HI<km/s>']/fcat['vel_dispersion_GMC<km/s>'],fcat['SFR<Msun/yr>'])
plt.loglog()

plt.show()


plt.scatter(fcat['vel_grad_HI<km/s/pc>']/fcat['vel_grad_GMC<km/s/pc>'],fcat['SFR<Msun/yr>'])
plt.loglog()

plt.show()


plt.scatter(fcat['vel_grad_GMC<km/s/pc>'],fcat['SFR<Msun/yr>'])
plt.loglog()

plt.show()


plt.scatter(fcat['vel_grad_HI<km/s/pc>'],fcat['SFR<Msun/yr>'])
plt.loglog()

plt.show()


plt.scatter(fcat['spec_ang_mom_GMC<pc.km/s>'],fcat['SFR<Msun/yr>'])
plt.loglog()

plt.show()


plt.scatter(fcat['spec_ang_mom_HI<pc.km/s>'],fcat['SFR<Msun/yr>'])
plt.loglog()

plt.show()


plt.hist(fcat['SFR<Msun/yr>']*1e6,bins=np.arange(0,2e3,1e2))

plt.show()


plt.scatter(fcat['mass_HI<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>']**2),fcat['SFR<Msun/yr>'])

plt.loglog()
plt.show()


plt.scatter(fcat['clustered'],fcat['SFR<Msun/yr>'])
# plt.loglog()
plt.show()


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
#     fig,axs = plt.subplots(1,3,figsize=(23,7))

    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    h1 = fns.get_mom1(his[i],vs[i])
    h2 = np.sqrt(fns.get_momN(his[i],vs[i],2))
    
    
    h0c = np.nansum(his_cosz[i],axis=0)*1.3*1.82e18
    h1c = fns.get_mom1(his_cosz[i],vs[i])
    h2c = np.sqrt(fns.get_momN(his_cosz[i],vs[i],2))
    
    
    h0s.append(h0)
    h1s.append(h1)
    h2s.append(h2)
    
    h0scsz.append(h0c)
    h1scsz.append(h1c)
    h2scsz.append(h2c)
    
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    c1 = fns.get_mom1(cos[i],vs[i])
    c2 = np.sqrt(fns.get_momN(cos[i],vs[i],2))
    
    c0h = np.nansum(cos_hisz[i],axis=0)*1.3*1.97e20
    c1h = fns.get_mom1(cos_hisz[i],vs[i])
    c2h = np.sqrt(fns.get_momN(cos_hisz[i],vs[i],2))
    
    
    c0s.append(c0)
    c1s.append(c1)
    c2s.append(c2)
    
    c0shsz.append(c0h)
    c1shsz.append(c1h)
    c2shsz.append(c2h)


# fig,ax = plt.subplots(1,1,figsize=(20,7))

for i in range(22):
    fig,axs = plt.subplots(1,3,figsize=(23,7))
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    
    lvs = np.arange(np.nanmean(h2s[i])-5*np.nanstd(h2s[i]),
                    np.nanmean(h2s[i])+5*np.nanstd(h2s[i]),
                    np.nanstd(h2s[i]))
    
    
    mask_sml = c0s[i] > 1e20
    mask_big = c0shsz[i] > 1e20
    
    wmask = axs[0].contourf(h2scsz[i]*mask_sml,cmap=dc.cmap_min0('YlGnBu_r'),
                            extent=[l0c,l1c,b0c,b1c],levels=lvs)
    
    
    womask = axs[1].contourf(h2s[i]*(~mask_big),cmap=dc.cmap_min0('YlGnBu_r'),
                           extent=[l0h,l1h,b0h,b1h],levels=lvs)
    
#     c0im = axs[1].contourf(c2shsz[i]*mask,cmap=dc.cmap_min0('magma'),
#                            extent=[l0h,l1h,b0h,b1h],levels=10,label='hello')
    
    
    #only the dispersion inside the co contours of that cloud factor
    #into the normalization
    nrmh2_big = h2s[i]/np.nanmean(c2s[i])
    nrmh2_sml = h2scsz[i]/np.nanmean(c2s[i])
    
    bns = np.arange(np.nanmin(nrmh2_big),np.nanmax(nrmh2_big),0.5) 
    #the big one has the larger range for most of the clouds
    
    
    sml_dat = np.ravel(nrmh2_sml*mask_sml)
    big_dat = np.ravel(nrmh2_big*(~mask_big))
                 
    #avoid all the 0s from the mask -- and hopefully no real data is affected
    #bc there are no los with linewidth=0
    axs[2].hist(sml_dat[sml_dat>0],bins=bns,density=True,
                alpha=0.5,label='within co bounds')
    
    axs[2].hist(big_dat[big_dat>0],bins=bns,density=True,
                alpha=0.5,label='outside co bounds')
#     axs[2].set_xlim(np.nanmin(nrmh2_big)-0.5,)
    axs[2].legend()
    axs[2].loglog()
    
#     fig.colorbar(wmask,ax=axs[0])
#     fig.colorbar(womask,ax=axs[1])
    fig.colorbar(womask,ax=axs[1])
    
#     axs[1].set_yticks()
    
    for ax in axs[:2]:
        ax.set_xlim(l0h,l1h)
        ax.set_ylim(b0h,b1h)
    
    plt.show()


#annotate this^ 
#look at email to nia to see wht else i said i'd do


#next steps patch by patch 


#also ^ examine sfr plots with c = cluster / individual, as well as s=mass


#for scatter plots that are hi<->gmc maybe useful to color code each cloud?


# fig,ax = plt.subplots(1,1,figsize=(20,7))

for i in range(22):
    fig,axs = plt.subplots(1,1,figsize=(12,12))

    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    h1 = fns.get_mom1(his[i],vs[i])
    h2 = np.sqrt(fns.get_momN(his[i],vs[i],2))
    
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    c1 = fns.get_mom1(cos[i],vs[i])
    c2 = np.sqrt(fns.get_momN(cos[i],vs[i],2))

    
    sph = np.nanmean(his[i],axis=(1,2))
    spc = np.nanmean(cos[i],axis=(1,2))
    
    
#     h1im = axs[0].contourf(h1,extent=[l0h,l1h,b0h,b1h],cmap='RdGy',
#                            levels=np.arange(np.nanmedian(c1)-5*np.nanstd(h1),
#                                             np.nanmedian(c1)+5*np.nanstd(h1),0.5*np.nanstd(h1)))
#     c1im = axs[0].contourf(c1,extent=[l0c,l1c,b0c,b1c],cmap='PiYG',
#                              levels=np.arange(np.nanmedian(c1)-3*np.nanstd(c1),
#                                             np.nanmedian(c1)+3*np.nanstd(c1),0.5*np.nanstd(c1)))
    
#     axs[0].contour(c0,extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
#                     colors='black',linewidths=0.4)
    
#     fig.colorbar(h1im,label=r'$\overline{v_{_{HI}}}\ $  [$km\ s^{-1}$]')
#     fig.colorbar(c1im,orientation='horizontal',label=r'$\overline{v_{_{H_2}}}\ $ [$km\ s^{-1}$]')
    axs.contourf(h2,extent=[l0h,l1h,b0h,b1h],levels=10)
#     axs.contourf(c2,extent=[l0c,l1c,b0c,b1c],levels=10)

    
    x,y = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))

    
    dh1db,dh1dl = np.gradient(h1)
    dh1mag = np.sqrt(dh1db**2+dh1dl**2)
    
#     axs[1].contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],levels=np.arange(20,22,0.1))
    axs.quiver(x,y,dh1dl,dh1db,alpha=0.7)
    
    
    dh0db,dh0dl = np.gradient(h0)
    dh0mag = np.sqrt(dh0db**2+dh0dl**2)
    
#     axs.quiver(x,y,dh0dl,dh0db,dh0mag,headwidth=1.5)
    

#     for ax in axs:
#         ax.set_xlim(l0h,l1h)
#         ax.set_ylim(b0h,b1h)
#         ax.set_ylabel(r'galactic latitude [$\degree$]')
#         ax.set_xlabel(r'galactic longitude [$\degree$]')

    
    plt.show()
    

# newvrange = np.array([[-10,  15], 
#        [-35,  35],
#        [-25,  25],
#        [-18,  25],
#        [-15,  29],
#        [-35,  35],
#        [-15,  29],
#        [-35,  40],
#        [-20,  60],
#        [-25,  35],
#        [-10,  35],
#        [-15,  42],
#        [-25,  15],
#        [-20,  27],
#        [-10,  35],
#        [-25,  40],
#        [-20,  30],
#        [-25,  35],
#        [-25,  35],
#        [-35,  25],
#        [-25,  30],
#        [-12,  30]])

# fcat['vel<km/s>'] = newvrange
# fcat.write('re_expanded_cat.fits',overwrite=True)


np.array((np.arange(22),fcat['name'])).transpose()


fig,ax = plt.subplots(1,1,figsize=(14,10))

for i in [6,13]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
                alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))

for i in [6,13]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
#     h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
#     ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
#                 alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contourf(np.log10(c0),extent=[l0c,l1c,b0c,b1c],cmap='binary_r',levels=np.arange(20,22,0.2))

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))


#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
#     plt.show()
    
# plt.xlim(180,-180)
# plt.ylim(-40,40)
# plt.show()


fig,ax = plt.subplots(1,1,figsize=(14,10))

for i in [10,14]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
                alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 

for i in [10,14]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
#     h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
#     ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
#                 alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contourf(np.log10(c0),extent=[l0c,l1c,b0c,b1c],cmap='binary_r',levels=np.arange(20,22,0.2))

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))


#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
#     plt.show()
    
# plt.xlim(180,-180)
# plt.ylim(-40,40)
# plt.show()


fig,ax = plt.subplots(1,1,figsize=(14,10))

for i in [7,12,20]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
                alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
for i in [7,12,20]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
#     h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
#     ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
#                 alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contourf(np.log10(c0),extent=[l0c,l1c,b0c,b1c],cmap='binary_r',levels=np.arange(20,22,0.2))

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))


#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
#     plt.show()
    
# plt.xlim(180,-180)
# plt.ylim(-40,40)
# plt.show()


fig,ax = plt.subplots(1,1,figsize=(14,10))

for i in [15,16,17]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
                alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
for i in [15,16,17]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
#     h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
#     ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
#                 alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contourf(np.log10(c0),extent=[l0c,l1c,b0c,b1c],cmap='binary_r',levels=np.arange(20,22,0.2))

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))


#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
#     plt.show()
    
# plt.xlim(180,-180)
# plt.ylim(-40,40)
# plt.show()


fig,ax = plt.subplots(1,1,figsize=(10,10))

for i in [1,18,21]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
                alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
for i in [1,18,21]:
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    
#     h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
#     ax.contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],
#                 alpha=0.5,cmap='Paired',levels=np.arange(20,22,0.2))
 
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contourf(np.log10(c0),extent=[l0c,l1c,b0c,b1c],cmap='binary_r',levels=np.arange(20,22,0.2))

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))


#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
#     plt.show()
    
# plt.xlim(180,-180)
# plt.ylim(-40,40)
# plt.show()


from scipy.ndimage import binary_dilation


clvdata = clv>1e20
wider = binary_dilation(clvdata,iterations=10)


plt.figure(figsize=(24,16))

plt.contour(hlv,levels=[1e19],colors='blue')

plt.contour(hlv*(clv>1e20),levels=[0,1],colors='red')
plt.contour(hlv*wider,levels=[0,1],colors='green')

plt.xlim(1440,4320)


hlvposb = np.nansum(h3sig[:60,640:,:],axis=0)*0.125*1.82e18*0.5
# cbv = np.nansum(d3sig,axis=2)*0.125*1.97e20*0.5


plt.figure(figsize=(24,16))

plt.contourf(hlvposb,levels=10,extent=[0,5760,640,800])
plt.xlim(1440,4320)
plt.ylim(0,800)


clv_map = np.broadcast_to(wider,(800,155,5760)).swapaxes(0,1)


h0lvvelrange = np.nansum(clv_map*h3sig,axis=0)*1.3*1.82e18


plt.figure(figsize=(24,16))

plt.contourf(np.log10(h0),cmap='YlGnBu_r',levels=np.arange(19,22.2,0.2))
plt.contour(np.log10(h0),colors='black',linestyles='--',levels=[20.6])
plt.contour(np.log10(h0lvvelrange),colors='white',linestyles='--',levels=[20.6])


plt.contourf(np.log10(c0),cmap='magma',levels=np.arange(20,23.2,0.2))

plt.xlim(1440,4320)


og_hi4pi = fits.open('HI4PI_CAR.fits')[0]


og_hi4pi.header


nv = 945
v0 = 471.92
dv = 1.288


all_hi_v = (np.arange(nv) - v0)*dv


np.argmin(np.abs(all_hi_v-200))


all_hi_v[317]


all_hi_v[627]


shirt200 = og_hi4pi.data[317:628,:,:]


sphi = np.nanmean(shirt200,axis=(1,2))


plt.plot(all_hi_v[317:628],sphi)
plt.hlines(0.043*3,-200,200)
plt.vlines([-110,110],0,1.5)
plt.ylim(0,2)


np.array([fcat['name'],np.nanmean(fcat['lon_GMC<deg>'],axis=1),np.nanmean(fcat['lat_GMC<deg>'],axis=1),
          fcat['radius_HI<pc>'],fcat['distance<pc>']]).transpose()


(np.arange(2163)-1082)*0.083 


(np.arange(4323)-2162)*0.083 


og_hi4pi.header


ora_c = (-149.,-19.)
ora_ra = 81.7*u.pc
ora_d = 432*u.pc


ora_ra_pxog = int((ora_ra/ora_d)*u.rad.to(u.deg)/0.083)


ora_lc_pxog, ora_bc_pxog = ((np.array(ora_c) - (180,-90))/(-0.083,0.083)).astype('int')


ora_hi_og = shirt200[:,
         ora_bc_pxog-ora_ra_pxog:ora_bc_pxog+ora_ra_pxog,
         ora_lc_pxog-ora_ra_pxog:ora_lc_pxog+ora_ra_pxog]


ora0_og = np.nansum(ora_hi_og,axis=0)*1.3*1.82e18


plt.contourf(np.log10(ora0_og),levels=np.arange(20,23,0.2))
plt.colorbar()


ora_ra_px0125 = int((ora_ra/ora_d)*u.rad.to(u.deg)/0.125)
ora_lc_px0125, ora_bc_px0125 = ((np.array(ora_c) - (360,-50))/(-0.125,0.125)).astype('int')


ora_hi_0125 = h3sig[:,
         ora_bc_px0125-ora_ra_px0125:ora_bc_px0125+ora_ra_px0125,
         ora_lc_px0125-ora_ra_px0125:ora_lc_px0125+ora_ra_px0125]


ora0_0125 = np.nansum(ora_hi_0125,axis=0)*1.3*1.82e18


np.nansum(ora0_0125)


np.nansum(ora0_og)


plt.contour(np.log10(ora0_0125),levels=np.arange(20,23,0.5),
            extent=[-149+82,-149-82,-19-82,-19+82],cmap='viridis',linewidths=0.6)
plt.contour(np.log10(ora0_og),levels=np.arange(20,23,0.5),
            extent=[-149+82,-149-82,-19-82,-19+82],cmap='magma',linewidths=0.6)

# plt.colorbar()


mass_n_area(ora_hi_0125*u.K,432*u.pc,'new')


mass_n_area(ora_hi_og[77:232,:,:]*u.K,432*u.pc,'og')


ora_1_og = fns.get_mom1(ora_hi_og[77:232,:,:],np.arange(155)*1.3-100)


ora_1_0125 = fns.get_mom1(ora_hi_0125,np.arange(155)*1.3-100)


plt.contourf(ora_1_0125,levels=np.arange(-10,30,3))
plt.colorbar()


ogl = np.meshgrid(np.arange(260),np.arange(260))[0]*0.083*u.deg.to(u.rad)*432


ogb = np.meshgrid(np.arange(260),np.arange(260))[1]*0.083*u.deg.to(u.rad)*432


o125l = np.meshgrid(np.arange(172),np.arange(172))[0]*0.125*u.deg.to(u.rad)*432


o125b = np.meshgrid(np.arange(172),np.arange(172))[1]*0.125*u.deg.to(u.rad)*432


o125l.shape


fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.,)


fttr(init,np.ravel(o125b),np.ravel(ora_1_0125))


fttr(init,np.ravel(o125l),np.ravel(ora_1_0125))


np.sqrt((0.05943857**2 + (-0.02231205)**2))


np.sqrt((0.05840184**2 + (-0.02178319)**2))


fcat


fttr(init,np.ravel(ogb),np.ravel(ora_1_og))


fttr(init,np.ravel(ogl),np.ravel(ora_1_og))


(-0.02178319+0.02231205)/0.02231205


plt.scatter(ogb,ora_1_og,s=1e-2)
plt.ylim(0,25)


plt.scatter(o125b,ora_1_0125,s=1e-2)


plt.scatter(o125l,ora_1_0125,s=1e-2)
plt.scatter(o125l,o125l*-0.02231205 + 10)


plt.scatter(ogl,ora_1_og,s=1e-2)
plt.ylim(0,25)


plt.contourf(ora_1_og,levels=np.arange(-10,30,4),extent=[-149+82,-149-82,-19-82,-19+82])
# plt.colorbar()
plt.contour(ora_1_0125,levels=np.arange(-10,30,4),cmap='magma',linewidths=0.7,
            extent=[-149+82,-149-82,-19-82,-19+82])

# plt.colorbar()


#there is about a 2% difference for mass calculated with different resolutions 
#just state and move on?


(8.37-8.17)/8.37


def mass_n_area(cube,d,cloud_type):
    w = np.nansum(cube*1.3*u.km/u.s,axis=0)
    conv = 1.82e18*((u.cm)**-2)/(u.K*(u.km/u.s))
    n = len(np.where(w*conv>4e20/u.cm**2)[0])
    if cloud_type == 'og':
        dA = (d*0.083*u.deg.to(u.rad))**2
        
        
    elif cloud_type == 'new':
        dA = (d*0.125*u.deg.to(u.rad))**2

        
    sig = conv*w*1.36*c.m_p
    
    m = np.nansum(sig)*dA
    
    return m.to(u.kg),(dA*n).to(u.pc**2)


masses_og = np.zeros(22)
radii_og = np.zeros(22)
for i in range(22):
    cube = his[i]*u.K
        
    d = fcat['distance<pc>'][i]*u.pc
    
    m,a = fns.mass_n_area(cube,d,'21cm')
    
    masses[i] = m.to(u.kg).value
    radii[i] = np.sqrt(a.to(u.pc**2)/np.pi).value
    
    
# fcat['mass_HI<kg>'] = masses
# fcat['radius_GMC<pc>'] = radii

# s, ra = fns.accumulation_radius(masses*u.kg,fcat['distance<pc>']*u.pc)
# sh = s.value
# fcat['radius_HI<pc>'] = ra
# fcat.write('re_expanded_cat.fits',overwrite=True)
