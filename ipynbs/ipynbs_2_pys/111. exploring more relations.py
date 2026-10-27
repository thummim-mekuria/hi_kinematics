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
    
    axs[0].tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False)
    
    #only the dispersion inside the co contours of that cloud 
    #factor into the normalization
    nrmh2_big = h2s[i]/np.nanmean(c2s[i])
    nrmh2_sml = h2scsz[i]/np.nanmean(c2s[i])
    
    bns = np.arange(np.nanmin(nrmh2_big),np.nanmax(nrmh2_big),0.5) 
    #the big one has the larger range for most of the clouds
    
    
    sml_dat = np.ravel(nrmh2_sml*mask_sml)
    big_dat = np.ravel(nrmh2_big*(~mask_big))
                 
    #avoid all the 0s from the mask -- 
    #and (p likely) no real data is affected bc there are no los with linewidth=0
    axs[2].hist(sml_dat[sml_dat>0],bins=bns,density=True,
                alpha=0.5,label='within co bounds')
    
    axs[2].hist(big_dat[big_dat>0],bins=bns,density=True,
                alpha=0.5,label='outside co bounds')
#     axs[2].set_xlim(np.nanmin(nrmh2_big)-0.5,)
    axs[2].legend()
    axs[2].set_yscale('log')
    axs[2].set_xlabel(r'$\sigma_{HI}\ $/ $\sigma_{H_2}$')
    
    
#     axs[2].loglog()
    
#     fig.colorbar(wmask,ax=axs[0])
#     fig.colorbar(womask,ax=axs[1])
    fig.colorbar(womask,ax=axs[0],location='left')
    
#     axs[1].set_yticks()
    
    for ax in axs[:2]:
        ax.set_xlim(l0h,l1h)
        ax.set_ylim(b0h,b1h)
    
    plt.show()


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
    
    axs[0].tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False)
    
    #only the dispersion inside the co contours of that cloud 
    #factor into the normalization
    nrmh2_big = h2s[i]#/np.nanmean(c2s[i])
    nrmh2_sml = h2scsz[i]#/np.nanmean(c2s[i])
    
    bns = np.arange(np.nanmin(nrmh2_big),np.nanmax(nrmh2_big),0.15) 
    #the big one has the larger range for most of the clouds
    
    
    sml_dat = np.ravel(nrmh2_sml*mask_sml)
    big_dat = np.ravel(nrmh2_big*(~mask_big))
                 
    #avoid all the 0s from the mask -- 
    #and (p likely) no real data is affected bc there are no los with linewidth=0
    axs[2].hist(sml_dat[sml_dat>0],bins=bns,density=True,
                alpha=0.5,label='within co bounds')
    
    axs[2].hist(big_dat[big_dat>0],bins=bns,density=True,
                alpha=0.5,label='outside co bounds')
#     axs[2].set_xlim(np.nanmin(nrmh2_big)-0.5,)
    axs[2].legend()
    axs[2].set_yscale('log')
    axs[2].set_xlabel(r'$\sigma_{HI}\ $/ $\sigma_{H_2}$')
    
    
#     axs[2].loglog()
    
#     fig.colorbar(wmask,ax=axs[0])
#     fig.colorbar(womask,ax=axs[1])
    fig.colorbar(womask,ax=axs[0],location='left')
    
#     axs[1].set_yticks()
    
    for ax in axs[:2]:
        ax.set_xlim(l0h,l1h)
        ax.set_ylim(b0h,b1h)
    
    plt.show()


#change the colorcoding in the histogram to better reflect the difference between
#within the co contour and without


#calculate the means of the two populations (vertical lines in the histograms)
#and see how those two relate with each other

#not a priority but think about reproducing the mach number plot from imara&burkhart


loners = fcat[fcat['clustered']==0]

#maybe also define this slightly differently -> ie how much other co emission is detected
#within the accumulation radius (ie crowded field or not? rather than i can establish that
#these two clouds are connected (could only do for orions, oph/lupus))


clustereds = fcat[fcat['clustered']==1]


nr,nc,axlocs = fns.get_ax_locs(len(loners))
fig,axs = plt.subplots(3,3,figsize=(20,20))
#go back and change so that if the number remainder 3 is 0, don't +1 for nr

for j in range(len(loners)):
    i = np.where(loners['name'][j]==fcat['name'])[0][0]
    
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    
    lvs = np.arange(np.nanmean(h2s[i])-5*np.nanstd(h2s[i]),
                    np.nanmean(h2s[i])+5*np.nanstd(h2s[i]),
                    np.nanstd(h2s[i]))
    
    
    mask_sml = c0s[i] > 1e20
    mask_big = c0shsz[i] > 1e20
    
    #only the dispersion inside the co contours of that cloud 
    #factor into the normalization
    nrmh2_big = h2s[i]/np.nanmean(c2s[i])
    nrmh2_sml = h2scsz[i]/np.nanmean(c2s[i])
    
    bns = np.arange(np.nanmin(nrmh2_big),np.nanmax(nrmh2_big),0.5) 
    #the big one has the larger range for most of the clouds
    
    
    sml_dat = np.ravel(nrmh2_sml*mask_sml)
    big_dat = np.ravel(nrmh2_big*(~mask_big))
                 
    #avoid all the 0s from the mask -- 
    #and (p likely) no real data is affected bc there are no los with linewidth=0
    r,c = axlocs[j]
    
    #come back and turn this into np.histogram then normalize the heights etc 
    #-> basically the same stuff as the npdfs
    
    
    axs[r,c].hist(sml_dat[sml_dat>0],bins=bns,density=True,
                alpha=0.5,label='within co bounds')
    
    axs[r,c].hist(big_dat[big_dat>0],bins=bns,density=True,
                alpha=0.5,label='outside co bounds')
#     axs[2].set_xlim(np.nanmin(nrmh2_big)-0.5,)
    
    axs[r,c].set_yscale('log')
    axs[r,c].set_xlabel(r'$\sigma_{HI}\ $/ $\sigma_{H_2}$')
    axs[r,c].set_title(loners['name'][j])
    
    
axs[0,2].legend(fontsize='small')    
fig.subplots_adjust(hspace=0.3)
for ax in axs.flat:
    ax.set_ylim(1e-3,1e0)
    ax.set_xlim(2,25)
     
# #     axs[2].loglog()
    
# #     fig.colorbar(wmask,ax=axs[0])
# #     fig.colorbar(womask,ax=axs[1])
#     fig.colorbar(womask,ax=axs[0],location='left')
    
# #     axs[1].set_yticks()
    
#     for ax in axs[:2]:
#         ax.set_xlim(l0h,l1h)
#         ax.set_ylim(b0h,b1h)
    
#     plt.show()


nr,nc,axlocs = fns.get_ax_locs(len(clustereds))
fig,axs = plt.subplots(nr,nc,figsize=(20,28))


for j in range(len(clustereds)):
    i = np.where(clustereds['name'][j]==fcat['name'])[0][0]
    
    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50
    
    
    lvs = np.arange(np.nanmean(h2s[i])-5*np.nanstd(h2s[i]),
                    np.nanmean(h2s[i])+5*np.nanstd(h2s[i]),
                    np.nanstd(h2s[i]))
    
    
    mask_sml = c0s[i] > 1e20
    mask_big = c0shsz[i] > 1e20
    
    #only the dispersion inside the co contours of that cloud 
    #factor into the normalization
    nrmh2_big = h2s[i]/np.nanmean(c2s[i])
    nrmh2_sml = h2scsz[i]/np.nanmean(c2s[i])
    
    bns = np.arange(np.nanmin(nrmh2_big),np.nanmax(nrmh2_big),0.5) 
    #the big one has the larger range for most of the clouds
    
    
    sml_dat = np.ravel(nrmh2_sml*mask_sml)
    big_dat = np.ravel(nrmh2_big*(~mask_big))
                 
    #avoid all the 0s from the mask -- 
    #and (p likely) no real data is affected bc there are no los with linewidth=0
    r,c = axlocs[j]
    
    #come back and turn this into np.histogram then normalize the heights etc 
    #-> basically the same stuff as the npdfs
    
    
    axs[r,c].hist(sml_dat[sml_dat>0],bins=bns,density=True,
                alpha=0.5,label='within co bounds')
    
    axs[r,c].hist(big_dat[big_dat>0],bins=bns,density=True,
                alpha=0.5,label='outside co bounds')
#     axs[2].set_xlim(np.nanmin(nrmh2_big)-0.5,)
    
    axs[r,c].set_yscale('log')
    axs[r,c].set_xlabel(r'$\sigma_{HI}\ $/ $\sigma_{H_2}$')
    axs[r,c].set_title(clustereds['name'][j])
    
    
axs[0,2].legend(fontsize='small')    
fig.subplots_adjust(hspace=0.5)
for ax in axs.flat:
    ax.set_ylim(1e-3,1e0)
    ax.set_xlim(2,25)
     
# #     axs[2].loglog()
    
# #     fig.colorbar(wmask,ax=axs[0])
# #     fig.colorbar(womask,ax=axs[1])
#     fig.colorbar(womask,ax=axs[0],location='left')
    
# #     axs[1].set_yticks()
    
#     for ax in axs[:2]:
#         ax.set_xlim(l0h,l1h)
#         ax.set_ylim(b0h,b1h)
    
#     plt.show()


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()

def fitlinelog(x,m,b):
    return (x**m)*(10**b)


# fcat['clustered'] = [0,1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1]

all_names = ['xkcd:easter purple',
 'xkcd:cool blue',
 'xkcd:purpleish pink',
 'xkcd:magenta',
 'xkcd:red',
 'xkcd:salmon',
 'xkcd:topaz',
 'xkcd:blue blue',
 'xkcd:pumpkin orange',
 'xkcd:amber',
 'xkcd:true green',
 'xkcd:light mustard',
 'xkcd:blue blue',
 'xkcd:topaz',
 'xkcd:true green',
 'xkcd:foam green',
 'xkcd:foam green',
 'xkcd:foam green',
 'xkcd:cool blue',
 'xkcd:lemon',
 'xkcd:blue blue',
 'xkcd:cool blue']


d_plane = fcat['distance<pc>']*np.cos(np.nanmean(fcat['lat_GMC<deg>'],axis=1)*u.deg)


x = d_plane*np.sin(np.nanmean(fcat['lon_GMC<deg>'],axis=1)*u.deg)
y = d_plane*np.cos(np.nanmean(fcat['lon_GMC<deg>'],axis=1)*u.deg)


for i in range(22):
    plt.scatter(x[i],y[i],label=fcat['name'][i],
                edgecolor=all_names[i],facecolor='none',
                s=fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass)*1e-3)
plt.xlim(2000,-2000)
plt.ylim(-2000,2000)
# plt.legend()


for i in range(22):
    plt.scatter(x[i],y[i],label=fcat['name'][i],
                edgecolor=all_names[i],facecolor='none',
                s=2*(fcat['vel_dispersion_HI<km/s>'][i]/fcat['vel_dispersion_GMC<km/s>'][i])**2)
plt.xlim(2000,-2000)
plt.ylim(-2000,2000)
plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5),fontsize='x-small')
# plt.legend()


# ysoN_names = np.array([['aquila_s', 62],
# ['california', 770],
# ['camelopardalis', np.nan],
# ['chamaeleon', 394],
# ['cma_ob1', 408],
# ['corona_australis', 28],
# ['crossbones', 176],
# ['gemini_ob1', 64],
# ['hercules',np.nan],
# ['lacerta', 525],
# ['lupus', 1059],
# ['maddalena', 17],
# ['mon_ob1', 577],
# ['mon_r2', 473],
# ['ophiuchus', 791],
# ['orion_a', 2072],
# ['orion_b', 854],
# ['orion_l', 1075],
# ['perseus', 435],
# ['polaris', 15],
# ['rosette', 78],
# ['taurus', 432]])


# ysos = ysoN_names.transpose()[1].astype('float')
# sfrs = 0.25*ysos*1e-6#*u.solMass/u.yr

# fcat['SFR<Msun/yr>'] = sfrs


# fcat.write('re_expanded_cat.fits',overwrite=True)


cluster_colors = ['black','orange']


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_grad_GMC<km/s/pc>'][i],fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
    
plt.loglog()
plt.xlabel(r'$\Omega_{GMC}\ $ [$km\ s^{-1}\ pc^{-1}$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_grad_GMC<km/s/pc>'][i]/fcat['vel_grad_HI<km/s/pc>'][i],fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$\Omega_{GMC}\ $ / $\Omega_{HI}$ ')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.show()


import astropy.constants as c #go back and change how axs[r,c] works to not intefere


virial_par = ((5*((fcat['vel_dispersion_GMC<km/s>']*u.km/u.s)**2)*(u.pc*fcat['radius_GMC<pc>'])\
             )/(c.G*(fcat['mass_GMC<kg>']*u.kg).to(u.solMass))).si


np.array((fcat['name'],virial_par)).transpose()


plt.figure(figsize=(7,7))

bins = np.arange(0,20,1)


plt.hist(virial_par[fcat['clustered']==0],density=True,color='black',alpha=0.7,bins=bins)
plt.hist(virial_par[fcat['clustered']==1],density=True,color='orange',alpha=0.7,bins=bins)


plt.show()


np.array((fcat['name'],fcat['spec_ang_mom_HI<pc.km/s>']/fcat['spec_ang_mom_GMC<pc.km/s>'])).transpose()


np.array((fcat['name'],fcat['SFR<Msun/yr>']*1e6)).transpose()


plt.scatter(fcat['name'],fcat['SFR<Msun/yr>'])
plt.xticks(rotation='vertical')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['name'][i],fcat['vel_dispersion_HI<km/s>'][i]/fcat['vel_dispersion_GMC<km/s>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*u.kg.to(u.solMass)*fcat['mass_GMC<kg>'][i])
    

# plt.loglog()
# plt.yscale('log')
# plt.xlabel(r'$\alpha_{vir,GMC}$')
# plt.ylabel(r'$M_{GMC}$ [$M_\odot$]')
# plt.xlim(0,6)
plt.xticks(rotation='vertical')
plt.show()

#can ish see general trend (ie less massive -> less bound),but that is ignoring
#the three largest alpha clouds (lacerta,hercules,lam-orionis) -> and also the
#spread seems quite large


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_dispersion_HI<km/s>'][i]/fcat['vel_dispersion_GMC<km/s>'][i],
                fcat['spec_ang_mom_HI<pc.km/s>'][i]/fcat['spec_ang_mom_GMC<pc.km/s>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*u.kg.to(u.solMass)*fcat['mass_GMC<kg>'][i])
    

# plt.loglog()
plt.yscale('log')
# plt.xlabel(r'$\alpha_{vir,GMC}$')
# plt.ylabel(r'$M_{GMC}$ [$M_\odot$]')
# plt.xlim(0,6)
# plt.xticks(rotation='vertical')
plt.show()

#can ish see general trend (ie less massive -> less bound),but that is ignoring
#the three largest alpha clouds (lacerta,hercules,lam-orionis) -> and also the
#spread seems quite large


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['name'][i],fcat['spec_ang_mom_HI<pc.km/s>'][i]/fcat['spec_ang_mom_GMC<pc.km/s>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*u.kg.to(u.solMass)*fcat['mass_GMC<kg>'][i])
    
# plt.loglog()
# plt.yscale('log')
# plt.xlabel(r'$\alpha_{vir,GMC}$')
# plt.ylabel(r'$M_{GMC}$ [$M_\odot$]')
# plt.xlim(0,6)
plt.xticks(rotation='vertical')
plt.show()

#can ish see general trend (ie less massive -> less bound),but that is ignoring
#the three largest alpha clouds (lacerta,hercules,lam-orionis) -> and also the
#spread seems quite large


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(virial_par[i],fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=10*fcat['spec_ang_mom_HI<pc.km/s>'][i]/fcat['spec_ang_mom_GMC<pc.km/s>'][i])
plt.loglog()
# plt.yscale('log')
plt.xlabel(r'$\alpha_{vir,GMC}$')
plt.ylabel(r'$M_{GMC}$ [$M_\odot$]')
# plt.xlim(0,6)
plt.show()

#can ish see general trend (ie less massive -> less bound),but that is ignoring
#the three largest alpha clouds (lacerta,hercules,lam-orionis) -> and also the
#spread seems quite large


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(virial_par[i],fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$\alpha_{vir,GMC}$')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['mass_HI<kg>'][i]*u.kg.to(u.solMass),fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$M_{HI}\ $ [$M_\odot$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_grad_HI<km/s/pc>'][i],fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$\Omega_{HI}\ $ [$km\ s^{-1}\ pc^{-1}$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_HI<pc.km/s>'][i],fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$j_{HI}\ $ [$km\ s^{-1}\ pc$]')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['spec_ang_mom_HI<pc.km/s>'][i]/fcat['spec_ang_mom_GMC<pc.km/s>'][i],
                fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$j_{HI}\ $ / $j_{GMC}$')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
# plt.xlim(0,20)
plt.show()


plt.figure(figsize=(7,7))

for i in range(22):
    plt.scatter(fcat['vel_dispersion_HI<km/s>'][i]/fcat['vel_dispersion_GMC<km/s>'][i],
                fcat['SFR<Msun/yr>'][i],
               color=cluster_colors[fcat['clustered'][i]],alpha=0.6,
                s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass))
plt.loglog()
plt.xlabel(r'$\sigma_{v,HI}\ $ / $\sigma_{v,GMC}$')
plt.ylabel(r'SFR [$M_\odot yr^{-1}$]')
# plt.xlim(0,20)
plt.show()


fcat[15]


#do zoomed in 0th moment maps in the catalog
#both for this ^ and the first moment map, reduce the amount
#of things going on
#eg for the first moment maps instead of a filled contour of the co, just do an arrow 
#to show the velocity gradient


#for indicating clustered / not clustered couds -> choose colors / shades closer together
#because with black and orange -> the jump is too big


#we care about quantities /  physical properties that indicate the evolutionary status of the cloud
#the sfr is related to evolution, and so is the ratio of j -> bc the gmc is expected to inherit the
#angular momentum of the parent hi

#another quantity to think about -> sfe -> do this first

#look at https://arxiv.org/pdf/1009.2985 -> how do they narrate expectations vs data


plt.figure(figsize=(7,7))

bins = np.arange(0,0.5,5e-2)

# plt.hist(fcat[fcat['clustered']==0]['vel_grad_GMC<km/s/pc>'],density=True,bins=bins,
#          color='black',alpha=0.6)

plt.hist(fcat[fcat['clustered']==0]['vel_grad_HI<km/s/pc>'],density=True,bins=bins,
         edgecolor='black',fill=False)

# plt.hist(fcat[fcat['clustered']==1]['vel_grad_GMC<km/s/pc>'],density=True,bins=bins,
#          color='orange',alpha=0.6)

plt.hist(fcat[fcat['clustered']==1]['vel_grad_HI<km/s/pc>'],density=True,bins=bins,
         edgecolor='orange',fill=False)


plt.xlabel(r'$\Omega\ $ [$km\ s^{-1}\ pc^{-1}$]')

plt.savefig('Omegahist')
plt.show()


plt.figure(figsize=(7,7))

bins = np.arange(0,0.5,5e-2)

plt.hist(fcat[fcat['clustered']==0]['vel_grad_GMC<km/s/pc>'],density=True,bins=bins,
         color='black',alpha=0.6)

# plt.hist(fcat[fcat['clustered']==0]['vel_grad_HI<km/s/pc>'],density=True,bins=bins,
#          edgecolor='black',fill=False)

plt.hist(fcat[fcat['clustered']==1]['vel_grad_GMC<km/s/pc>'],density=True,bins=bins,
         color='orange',alpha=0.6)

# plt.hist(fcat[fcat['clustered']==1]['vel_grad_HI<km/s/pc>'],density=True,bins=bins,
#          edgecolor='orange',fill=False)


plt.xlabel(r'$\Omega\ $ [$km\ s^{-1}\ pc^{-1}$]')

plt.savefig('Omegahist')
plt.show()


clustereds


clustereds['mass_GMC<kg>']*1e-3*u.kg.to(u.solMass)


plt.figure(figsize=(7,7))

bins = np.arange(0,30,2.5)

plt.hist(fcat[fcat['clustered']==0]['mass_GMC<kg>']*u.kg.to(u.solMass)*1e-3,bins=bins,density=True,
         color='black',alpha=0.7)
plt.hist(fcat[fcat['clustered']==1]['mass_GMC<kg>']*u.kg.to(u.solMass)*1e-3,bins=bins,density=True,
         color='orange',alpha=0.7)

# plt.hist(fcat[fcat['clustered']==0]['mass_HI<kg>']*u.kg.to(u.solMass)*1e-4,bins=bins,density=True,
#          edgecolor='black',alpha=0.7,fill=False)
# plt.hist(fcat[fcat['clustered']==1]['mass_HI<kg>']*u.kg.to(u.solMass)*1e-4,bins=bins,density=True,
#          edgecolor='orange',alpha=0.7,fill=False)


plt.xlabel(r'$M\ $ [$10^3 \times M_\odot$]')
# plt.xlim(0,20)
# plt.savefig('masshist')
# plt.xscale('log')
plt.show()


fcat


#for the plots with j -> change the units to pc km/s
#overplot the trends from m33
#add in nias milky way clouds

#ratio of js -> does the new trend line better or exacerbate the angular momentum problem


#for the histogram of omega -> separate the co and hi
#look at how the negative side of the histogram is defined

#make the plots as similar + clear as possible
#to make comparisons easier

#slow down and describe what's going on in the plots

#also indicate on the plots which clouds do / don't have linear gradients

#look at the bottom left in the four panel figures

#when choosing clouds to zoom in on -> do they fit linear gradients well or not
#how does that compare with other stuff going on in them


#for the gmcs that fit within the same hi cloud, calculate the same properties
#for those associations

#for flash -> give the background, use plots from previous work to give that context


#reproducing npdfs from imara & burkhart


bns = np.arange(1e20,1e22,5e19)
npdfs_h = []
npdfs_hcsz = []
npdfs_c = []


for i in range(22):
#     plt.contourf(np.log10(),levels=)

    thisc0 = np.ravel(c0s[i])
    thish0 = np.ravel(h0s[i])
    thish0csz = np.ravel(h0scsz[i])*np.ravel(c0s[i]>1e20)


    thishi = np.histogram(thish0[thish0>1e20],bins=bns)
    thisco = np.histogram(thisc0[thisc0>1e20],bins=bns)
    thishicsz = np.histogram(thish0csz[thish0csz>1e20],bins=bns)
    
    
    normhi = thishi[0]/np.nanmax(thishi[0])
    normco = thisco[0]/np.nanmax(thisco[0])
    normhicsz = thishicsz[0]/np.nanmax(thishi[0])

    
    npdfs_h.append(normhi)
    npdfs_c.append(normco)
    npdfs_hcsz.append(normhicsz)


for i in range(22):
    fig,axs = plt.subplots()
    plt.plot(bns[:-1],npdfs_h[i],color='pink',label='outside co bounds')
    plt.plot(bns[:-1],npdfs_hcsz[i],color='red',label='within co bounds')
    plt.plot(bns[:-1],npdfs_c[i],color='grey',label='co')
    
    
    plt.title(fcat['name'][i])#+' '+cl_st[fcat['clustered'][i]])
    plt.loglog()
    plt.xlabel(r'$N_{H}\ $ [$cm^{-2}$]')
    plt.ylabel(r'$N/N_{max}$')
    
    plt.ylim(1e-3,)
#     plt.legend(loc='upper left',fontsize='xx-small')
#     plt.savefig(fcat['name'][i]+'_npdf_mrh')
    plt.show()


fttr_new = fitting.LevMarLSQFitter()


def prob_nh (data,No=1,sigma=1):
    coeff = (sigma*np.sqrt(2*np.pi))**-1
    num = (np.log(data) - np.log(No))**2
    den = 2*(sigma**2)
    
    
    return coeff*np.e**(-num/den)


prob_nh_model = models.custom_model(prob_nh)


init_new = prob_nh_model(3,0.1)


for i in range(22):
    no,sig = fttr_new(init_new,bns[:-1]/1e21,npdfs_h[i]).parameters
    
    plt.plot(bns/1e21,prob_nh(bns/1e21,sigma=sig,No=no))
    plt.plot(bns[:-1]/1e21,npdfs_h[i])

    plt.loglog()

    plt.ylim(1e-3,)
    plt.show()


fcat


nr,nc,axsn = fns.get_ax_locs(12)


np.array([[0, 0],
       [0, 1],
       [0, 2],
       [0, 3],
       [1, 0],
       [1, 1],
       [1, 2],
          []
       [2, 0],
       [2, 1],
       [2, 2],
       [3, 0],
       [3, 1],
       [3, 2],
       [4, 0],
       [4, 1],
       [4, 2]])


fig,axs = plt.subplots(3,3,figsize=(35,20))
nr,nc,axsn = fns.get_ax_locs(9)


for j in range(9):
    i = [1,3,8,9,10,13,15,18,19][j]
    
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
    
    rn,cn = axsn[j]
    h0im = axs[rn,cn].contourf(np.log10(h0),extent=[l0h,l1h,b0h,b1h],cmap='YlGnBu_r',
                               levels=np.arange(20,22,0.1))
    
    c0im = axs[rn,cn].contourf(np.log10(c0),extent=[l0c,l1c,b0c,b1c],cmap='magma',
                               levels=np.arange(20,23,0.2))
    
    axs[rn,cn].contour(c0,extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
                    colors='black',linewidths=0.4)
    
#     fig.colorbar(h1im,label=r'$\overline{v_{_{HI}}}\ $  [$km\ s^{-1}$]')
#     fig.colorbar(c1im,orientation='horizontal',label=r'$\overline{v_{_{H_2}}}\ $ [$km\ s^{-1}$]')
    
    
#     h2im = axs[1,0].contourf(h2,extent=[l0h,l1h,b0h,b1h],cmap='YlGnBu_r',levels=6)
    
#     c2im = axs[1,0].contourf(c2,extent=[l0c,l1c,b0c,b1c],cmap='magma',levels=6)
    
#     axs[1,0].contour(c0,extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
#                     colors='black',linewidths=0.4)    
#     fig.colorbar(h2im,label=r'$\sigma_{v,_{HI}}\ $  [$km\ s^{-1}$]')
#     fig.colorbar(c2im,orientation='horizontal',label=r'$\sigma_{v,_{H_2}}\ $  [$km\ s^{-1}$]')
    
    
#     hspec = axs[0,1].plot(vs[i],sph,color='black')
    
#     axs01 = axs[0,1].twinx()
#     cspec = axs01.plot(vs[i],spc,color='xkcd:kelly green')
#     axs[0,1].set_xlabel(r'$v_{_{cloud}}\ $  [$km\ s^{-1}$]')
#     axs[0,1].set_ylabel(r'$T_{B,HI}\ $  [$K$]')
#     axs01.set_ylabel(r'$T_{B,H_2}\ $  [$K$]')
    
    
# #     print(np.nanmean(c2)+np.nanstd(c2),np.nanmax(c2))
#     for ax in axs.flat:
    axs[rn,cn].set_xlim(l0h,l1h)
    axs[rn,cn].set_ylim(b0h,b1h)
    
    
    axs[rn,cn].set_title(fcat['name'][i])
    
for ax in [axs[0,0],axs[1,0],axs[2,0]]:
    ax.set_ylabel(r'galactic latitude [$\degree$]')
    
for ax in [axs[2,0],axs[2,1],axs[2,2]]:
    ax.set_xlabel(r'galactic longitude [$\degree$]')

# fig.colorbar(h1im)

cbar_ax = fig.add_axes([0.99, 0.11, 0.05, 0.38])
cbar_ax2 = fig.add_axes([0.99, 0.5, 0.05, 0.38])

fig.colorbar(h0im, cax=cbar_ax,label=r'log ($N_{HI}$ /$\ cm^{-2}$)')
fig.colorbar(c0im, cax=cbar_ax2,label=r'log ($N_{H_2}$ /$\ cm^{-2}$)')


plt.subplots_adjust(wspace=0.15,hspace=0.2)
fig.tight_layout()

fig.savefig('selected_clouds_0th',bbox_inches="tight")
    
    
plt.show()
    

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
