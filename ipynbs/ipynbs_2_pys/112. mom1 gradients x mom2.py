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


loners_true = [0,2,3,5,8,9,19] #there is nothing nearby
loners_ish = [4,10,11,12,13,20] #there is a little bit of stuff nearby
crowd = [1,6,7,14,15,16,17,18,21] #there is a lot of stuff nearby


plt.scatter


loners = fcat[fcat['clustered']==0]

#maybe also define this slightly differently -> ie how much other co emission is detected
#within the accumulation radius (ie crowded field or not? rather than i can establish that
#these two clouds are connected (could only do for orions, oph/lupus))


clustereds = fcat[fcat['clustered']==1]


sml_dt = [np.ravel((c0s[i]>1e20)*(h2scsz[i]/np.nanmean(c2s[i])))
for i in range(22)]


bg_dt = [np.ravel((c0shsz[i]>1e20)*(h2s[i]/np.nanmean(c2s[i])))
for i in range(22)]


sml_dt_hist = [np.histogram(sd[sd>0],bins=np.arange(0,25,0.25))[0] for sd in sml_dt]
bg_dt_hist = [np.histogram(bd[bd>0],bins=np.arange(0,25,0.25))[0] for bd in bg_dt]


sml_pdf = [sml/np.max(sml) for sml in sml_dt_hist]
bg_pdf = [bg/np.max(bg) for bg in bg_dt_hist]


len(sml_pdf)


for i in range(22):
    wo = np.array(bg_pdf)[i]
    ww = np.array(sml_pdf)[i]
    plt.bar(np.arange(0,24.75,0.25),wo,0.25,alpha=0.5)
    plt.vlines(fns.get_mom1(wo,np.arange(0,24.75,0.25)),0,1,linestyle='--',color='red')
    
    plt.step(np.arange(0,24.75,0.25),ww,0.25)
    plt.vlines(fns.get_mom1(ww,np.arange(0,24.75,0.25)),0,1,linestyle='--',color='orange')


    plt.title(fcat['name'][i])
    plt.show()


for i in range(22):
    plt.scatter(np.ravel(h1s[i]),np.ravel(h2s[i]),s=5)
#     plt.scatter(np.ravel(c1s[i]),np.ravel(c2s[i]),s=5)

    plt.xlabel('lsr velocity')
    plt.ylabel('linewidth')
    plt.title(fcat['name'][i])
    
#     plt.loglog()
#     plt.ylim(1,)
    plt.show()


for i in loners_true:
    wo = np.array(bg_pdf)[i]
    ww = np.array(sml_pdf)[i]
    plt.bar(np.arange(0,24.75,0.25),wo,0.25,alpha=0.5)
    plt.step(np.arange(0,24.75,0.25),ww)

    plt.title(fcat['name'][i])
    plt.show()


loner_true_pdf = np.nansum(np.array(sml_pdf)[loners_true],axis=0)
loner_ish_pdf = np.nansum(np.array(sml_pdf)[loners_ish],axis=0)
crowd_pdf = np.nansum(np.array(sml_pdf)[crowd],axis=0)


bloner_true_pdf = np.nansum(np.array(bg_pdf)[loners_true],axis=0)
bloner_ish_pdf = np.nansum(np.array(bg_pdf)[loners_ish],axis=0)
bcrowd_pdf = np.nansum(np.array(bg_pdf)[crowd],axis=0)


plt.bar(np.arange(0,24.75,0.25),loner_true_pdf,0.25,alpha=0.7,color='black')
# plt.step(np.arange(0,24.5,0.5),bloner_true_pdf,alpha=0.5,linestyle='--')

plt.bar(np.arange(0,24.75,0.25),loner_ish_pdf,0.25,alpha=0.5,color='xkcd:grey')
# plt.step(np.arange(0,24.5,0.5),bloner_ish_pdf,0.5,alpha=0.5,linestyle='--')

plt.bar(np.arange(0,24.75,0.25),crowd_pdf,0.25,alpha=0.5,color='xkcd:light blue')
# plt.step(np.arange(0,24.5,0.5),bcrowd_pdf,0.5,alpha=0.5,linestyle='--')

# plt.loglog()
# plt.xlim(2,)
plt.show()


h0gal = np.nansum(h3sig,axis=0)*1.3*1.82e18


c0gal = np.nansum(d3sig,axis=0)*1.3*1.97e20


#quantify crowdedness as how far out in radius you have to go (in any direction) before running into co that
#is not that cloud (ig maybe a way to do it is to calculate the total co mass at different values of r from
#the center, and when that value starts to change dramatically then thats when you know)


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
fig,axs = plt.subplots(nr,nc,figsize=(20,20))
#go back and change so that if the number remainder 3 is 0, don't +1 for nr

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


loners_ish


loners_true


new_boxes = np.array([0,2,0,0,1,0,2,2,0,0,1,1,1,1,2,2,2,2,2,0,1,2])


loner_color = ['black','grey','xkcd:light blue']


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()


jj_rr_fit = fttr(init,np.log10(fcat['radius_GMC<pc>']/fcat['radius_HI<pc>']),
                 np.log10(fcat['spec_ang_mom_GMC<pc.km/s>']/fcat['spec_ang_mom_HI<pc.km/s>']))


jj_rr_fit


for i in range(22):
    plt.scatter(fcat['radius_GMC<pc>'][i]/fcat['radius_HI<pc>'][i],
        fcat['spec_ang_mom_GMC<pc.km/s>'][i]/fcat['spec_ang_mom_HI<pc.km/s>'][i],
        color='black',s=fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass)*5e-3,
        alpha=0.5)

plt.plot(np.arange(1e-1,5e-1,1e-3),
         10**jj_rr_fit.intercept*np.arange(1e-1,5e-1,1e-3)**jj_rr_fit.slope,
         linestyle='--',color='xkcd:dark blue',
        label='fit to data, n $\sim$'+str(round(jj_rr_fit.slope.value)))

plt.hlines(1,0.1,0.5,label='$j_{HI} = j_{GMC}$',color='xkcd:dark pink')

plt.loglog()
plt.xlim(1e-1,5e-1)
plt.xticks([0.1,0.2,0.3,0.4,0.5],[r'$10^{-1}$',0.2,0.3,0.4,0.5])
plt.xlabel(r'$R_{GMC}$ / $R_{acc,HI}$')
plt.ylabel(r'$j_{GMC}$ / $j_{HI}$')

plt.legend(fontsize='xx-small',loc='upper left')
plt.savefig('jjrr',bbox_inches="tight")
plt.show()


omom_rr_fit = fttr(init,np.log10(fcat['radius_GMC<pc>']/fcat['radius_HI<pc>']),
                 np.log10(fcat['vel_grad_GMC<km/s/pc>']/fcat['vel_grad_HI<km/s/pc>']))


omom_rr_fit


((((375*(c.k_B**3))/(4*np.pi*\
        ((c.m_p*1.36)**4)*c.G**3))**0.5).si)/c.M_sun/1e3*((20*u.K)**1.5)


loner_color = ['xkcd:hunter green','xkcd:seafoam green','xkcd:grass green']


for i in range(22):
    plt.scatter(fcat['radius_GMC<pc>'][i]/fcat['radius_HI<pc>'][i],
        fcat['vel_grad_GMC<km/s/pc>'][i]/fcat['vel_grad_HI<km/s/pc>'][i],
        color='black',s=5e-3*fcat['mass_GMC<kg>'][i]*u.kg.to(u.solMass),
        alpha=0.5)

plt.plot(np.arange(1e-1,5e-1,1e-3),
         10**omom_rr_fit.intercept*np.arange(1e-1,5e-1,1e-3)**omom_rr_fit.slope,
         linestyle='--',color='xkcd:dark blue',
        label='fit to data, n $\sim$'+str(round(omom_rr_fit.slope.value)))

plt.plot(np.arange(1e-1,5e-1,1e-3),
         np.arange(1e-1,5e-1,1e-3)**-2,
         linestyle='--',color='xkcd:dark pink',
        label='angular momentum conserved, $n = -2$')


plt.loglog()
plt.xlim(1e-1,5e-1)
plt.xticks([0.1,0.2,0.3,0.4,0.5],[r'$10^{-1}$',0.2,0.3,0.4,0.5])
plt.xlabel(r'$R_{GMC}$ / $R_{acc,HI}$')
plt.ylabel(r'$\Omega_{GMC}$ / $\Omega_{HI}$')

plt.legend(fontsize='xx-small',loc='upper right')
plt.savefig('omomrr',bbox_inches="tight")
plt.show()


plt.hist(fcat['mass_GMC<kg>'][new_boxes==0]*u.kg.to(u.solMass)*1e-3,density=True)
plt.hist(fcat['mass_GMC<kg>'][new_boxes==1]*u.kg.to(u.solMass)*1e-3,density=True)
plt.hist(fcat['mass_GMC<kg>'][new_boxes==2]*u.kg.to(u.solMass)*1e-3,density=True)

plt.loglog()
plt.show()

#hmmmm is this a real difference between them? like the true loners are systematically less massive than
#the crowded field ones, and those are systematically less massive than the fake loners (mostly b ~0)


vir_pars_GMC = ((5*(fcat['vel_dispersion_GMC<km/s>']*u.km/u.s)**2 *\
            fcat['radius_GMC<pc>']*u.pc)/(c.G*(fcat['mass_GMC<kg>']*u.kg).to(u.solMass))).si

vir_pars_HI = ((5*(fcat['vel_dispersion_HI<km/s>']*u.km/u.s)**2 *\
            fcat['radius_HI<pc>']*u.pc)/(c.G*(fcat['mass_HI<kg>']*u.kg).to(u.solMass))).si


plt.hist(vir_pars_GMC[new_boxes==0],density=True,color='black')
plt.hist(vir_pars_GMC[new_boxes==1],density=True,color='grey',alpha=0.6)
plt.hist(vir_pars_GMC[new_boxes==2],density=True,color='xkcd:light blue',alpha=0.5)


plt.hist(vir_pars_HI[new_boxes==0],density=True,fill=False,edgecolor='black')
plt.hist(vir_pars_HI[new_boxes==1],density=True,fill=False,edgecolor='grey',alpha=0.6)
plt.hist(vir_pars_HI[new_boxes==2],density=True,fill=False,edgecolor='xkcd:light blue',alpha=0.5)


# plt.loglog()
plt.show()

#hmmmm is this a real difference between them? like the true loners are systematically less massive than
#the crowded field ones, and those are systematically less massive than the fake loners (mostly b ~0)


brot_GMC = (((1/3)*((fcat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)**2)
                  *(fcat['radius_GMC<pc>']*u.pc)**3)/(c.G*fcat['mass_GMC<kg>']*u.kg)).si

brot_HI = (((1/3)*((fcat['vel_grad_HI<km/s/pc>']*u.km/u.s/u.pc)**2)
                  *(fcat['radius_HI<pc>']*u.pc)**3)/(c.G*fcat['mass_HI<kg>']*u.kg)).si


brot_GMC.si


plt.hist(brot_HI[new_boxes==0],density=True,fill=False,edgecolor='black')
plt.hist(brot_HI[new_boxes==1],density=True,fill=False,edgecolor='grey',alpha=0.6)
plt.hist(brot_HI[new_boxes==2],density=True,fill=False,edgecolor='xkcd:light blue',alpha=0.5)


# plt.loglog()
plt.show()

#hmmmm is this a real difference between them? like the true loners are systematically less massive than
#the crowded field ones, and those are systematically less massive than the fake loners (mostly b ~0)


fcat['name']


fig,axs = plt.subplots(2,5,figsize=(38,17))


for j in range(10):
    i = [1,3,4,8,9,10,13,15,18,19][j]
    ax = axs.flat[j]
    
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
    
    
    h0im = ax.contourf(h1,extent=[l0h,l1h,b0h,b1h],cmap='RdGy',
                               levels=np.arange(np.nanmedian(c1)-5*np.nanstd(h1),
                                            np.nanmedian(c1)+5*np.nanstd(h1),0.5*np.nanstd(h1)))
    
    c0im = ax.contourf(c1,extent=[l0c,l1c,b0c,b1c],cmap='PiYG',
                               levels=np.arange(np.nanmedian(c1)-3*np.nanstd(c1),
                               np.nanmedian(c1)+3*np.nanstd(c1),0.5*np.nanstd(c1)))
    
    ax.contour(c0,extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
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
    ax.set_xlim(l0h,l1h)
    ax.set_ylim(b0h,b1h)
    
    
    ax.set_title(fcat['name'][i])

    
for ax in [axs[0,0],axs[1,0]]:
    ax.set_ylabel(r'galactic latitude [$\degree$]')
    
for ax in [axs[1,0],axs[1,1],axs[1,2],axs[1,3],axs[1,4]]:
    ax.set_xlabel(r'galactic longitude [$\degree$]')

# fig.colorbar(h1im)

cbar_ax = fig.add_axes([1.01, 0.11, 0.015, 0.38])
cbar_ax2 = fig.add_axes([1.01, 0.53, 0.015, 0.38])

fig.colorbar(h0im, cax=cbar_ax,label=r'$v_{lsr,HI}\ [km\ s^{-1}]$')
fig.colorbar(c0im, cax=cbar_ax2,label=r'$v_{lsr,H_2}\ [km\ s^{-1}]$')


plt.subplots_adjust(wspace=0.1,hspace=0.2)
fig.tight_layout()

fig.savefig('selected_clouds_1th',bbox_inches="tight")
    
    
plt.show()
    

plt.scatter(np.nanmean(fcat['lat_GMC<deg>'],axis=1),
            fcat['rot_ax_GMC<deg>'])
# plt.xlim(180,-180)


plt.rcParams['font.size'] = 32


fcat


plt.rcParams['font.size'] = 20


theta_diff = np.abs(fcat['rot_ax_GMC<deg>']-fcat['rot_ax_HI<deg>'])

theta_nrmd = np.concatenate((theta_diff,np.abs(theta_diff-360)))

plt.hist(theta_nrmd,bins = np.arange(0,360,30),color='xkcd:blue',alpha=0.6,edgecolor='xkcd:blue')

plt.vlines([60,120],0,7,color='black',linestyle='--')
# plt.vlines([45,135,225,315],0,7,color='black',linestyle='--')

# plt.annotate('corotating',(0,7),fontsize=13)
# plt.annotate('misaligned',(80,7),fontsize=13)
# plt.annotate('counter-rotating',(140,7),fontsize=13)

plt.xlabel(r'|$\Delta \phi$| [$\degree$]')
plt.ylabel(r'$N_{clouds}$')

plt.xticks([0,60,120,180,240,300,],[0,60,120,180,240,300])
# plt.xticks([0,45,90,135,180,225,270,315,360],[0,45,90,135,180,225,270,315,360])

plt.xlim(0,180)

plt.savefig('delta_theta',bbox_inches="tight")
plt.show()


i=5
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

plt.savefig(fcat['name'][i]+'_mom1')

plt.show()


2.6-1.8


8.6-6.8


i=1
l0c,l1c = np.array((loc[i],lic[i]))*-0.125 + 360
b0c,b1c = np.array((boc[i],bic[i]))*0.125 -50

l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

plt.figure(figsize=(18,10))
plt.contourf(h1s[1],extent=[l0h,l1h,b0h,b1h],cmap='RdGy',
               levels=np.arange(np.nanmedian(c1s[1])-4*np.nanstd(h1s[1]),
                np.nanmedian(c1s[1])+4*np.nanstd(h1s[1]),0.25*np.nanstd(h1s[1])))

plt.colorbar(label=r'$v_{lsr}\ [km\ s^{-1}]$',shrink=0.8,
             ticks=[-24.5,-13,-2.6,7,18.5])


plt.contourf(c1s[1],extent=[l0c,l1c,b0c,b1c],cmap='PiYG',
                               levels=np.arange(np.nanmedian(c1s[1])-3*np.nanstd(c1s[1]),
                               np.nanmedian(c1s[1])+3*np.nanstd(c1s[1]),0.25*np.nanstd(c1s[1])))

plt.colorbar(shrink=0.8,ticks=[-13,-7.5,-2.6,2.5,7]) #label=r'$v_{lsr,H_2}\ [km\ s^{-1}]$'

plt.contour(c0s[1],extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
                    colors='black',linewidths=0.4)
plt.xlim(l0h,l1h)


plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

plt.title(fcat['name'][1])

# plt.plot(np.arange(170,150,-0.125),np.arange(170,150,-0.125)*\
#          np.tan((-163.1979151135055)*u.deg)-60)


# plt.plot(np.arange(166.5,155.0,-0.125),np.arange(166.5,155.0,-0.125)*\
#          np.tan((-51.95768563566319)*u.deg)+198)

plt.savefig('california_mom1')

plt.show()


plt.figure(figsize=(15,10))
plt.contourf(h1s[1],extent=[l0h,l1h,b0h,b1h],colors='white',
               levels=np.arange(np.nanmedian(c1s[1])-4*np.nanstd(h1s[1]),
                np.nanmedian(c1s[1])+4*np.nanstd(h1s[1]),0.25*np.nanstd(h1s[1])))

# plt.colorbar(label=r'$v_{lsr}\ [km\ s^{-1}]$',shrink=0.8,
#              ticks=[-24.5,-13,-2.6,7,18.5])


plt.contourf(c1s[1],extent=[l0c,l1c,b0c,b1c],cmap='PiYG',
                               levels=np.arange(np.nanmedian(c1s[1])-3*np.nanstd(c1s[1]),
                               np.nanmedian(c1s[1])+3*np.nanstd(c1s[1]),0.25*np.nanstd(c1s[1])))

plt.colorbar(shrink=0.8,ticks=[-13,-7.5,-2.6,2.5,7],label=r'$v_{lsr,H_2}\ [km\ s^{-1}]$')

plt.contour(c0s[1],extent=[l0c,l1c,b0c,b1c],levels=[1e20,5e20,1e21],
                    colors='black',linewidths=0.4)
plt.xlim(l0h,l1h)


plt.ylabel(r'galactic latitude [$\degree$]')
plt.xlabel(r'galactic longitude [$\degree$]')

plt.title(fcat['name'][1])

# plt.plot(np.arange(170,150,-0.125),np.arange(170,150,-0.125)*\
#          np.tan((-163.1979151135055)*u.deg)-60)


# plt.plot(np.arange(166.5,155.0,-0.125),np.arange(166.5,155.0,-0.125)*\
#          np.tan((-51.95768563566319)*u.deg)+198)

plt.savefig('california_mom1_justh2')
# plt.savefig('california_mom1_justh')


plt.show()


plt.quiver(160,-10,np.cos(51.95768*u.deg).value,np.sin(51.95768*u.deg).value)
plt.xlim(170,150)


cal_lons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])


cal_lats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])


cal_mom1 = np.ravel(h1s[1])


fcat['rot_ax_HI<deg>'][1]


plt.plot(np.arange(170,150,-0.125),np.arange(170,150,-0.125)*np.tan(-163.19*u.deg))


b_perp = np.abs(cal_lats - cal_lons*np.tan(-163.19*u.deg))


plt.imshow(b_perp.reshape(151,151))
plt.colorbar()


150*np.tan(-163.19*u.deg)


bperp


l0h-l1h


39//0.5


for i in range(22):
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

    
    lons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    lats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    mom1 = np.ravel(h1s[i])

    
    lc = (l0h+l1h)/2
    bc = (b0h+b1h)/2
    thet = fcat['rot_ax_HI<deg>'][i]*u.deg
    
    plt.scatter(lons,lats,c=mom1,cmap='RdGy')
    plt.colorbar()

    rot_ax_line = bc+(lons-lc)*np.tan(thet)
    
    
    plt.scatter(lons,
            bc+(lons-lc)*np.tan(thet))

    plt.xlim(l0h,l1h)
    plt.ylim(b0h,b1h)
    plt.show()
    
    bperp = np.sqrt((lats - rot_ax_line)**2)
    model_v = np.nanmean(mom1) + fcat['vel_grad_HI<km/s/pc>'][i]*np.cos(thet+90*u.deg)*\
                (lons-np.nanmean((l0h,l1h)))*u.deg.to(u.rad)*fcat['distance<pc>'][i] + \
                fcat['vel_grad_HI<km/s/pc>'][i]*np.sin(thet+90*u.deg)*\
                (lats-np.nanmean((b0h,b1h)))*u.deg.to(u.rad)*fcat['distance<pc>'][i] 
    
    plt.imshow(model_v.reshape(h1s[i].shape),cmap='RdGy',
               extent=[l0h,l1h,b0h,b1h],origin='lower')
    plt.colorbar()

    plt.scatter(lons,
            bc+(lons-lc)*np.tan(thet))
#     plt.imshow(bperp.reshape(h1s[i].shape),extent=[l0h,l1h,b0h,b1h],origin='lower')
    
    plt.xlim(l0h,l1h)
    plt.ylim(b0h,b1h)
    plt.show()


for i in range(22):
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125 -50

    
    lons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    lats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    mom1 = np.ravel(h1s[i])

    
    lc = (l0h+l1h)/2
    bc = (b0h+b1h)/2
    thet = fcat['rot_ax_HI<deg>'][i]*u.deg

    rot_ax_line = bc+(lons-lc)*np.tan(thet)
    bperp = lats - rot_ax_line
#     plt.scatter(lons,lats,c=mom1,cmap='RdGy')
    plt.contourf(bperp.reshape(h1s[i].shape).value,extent=[l0h,l1h,b0h,b1h])
    plt.colorbar()
    plt.scatter(lons,rot_ax_line)
    plt.xlim(l0h,l1h)
    plt.ylim(b0h,b1h)
    
    plt.show()
    
    
    model_v = np.nanmean(mom1) + fcat['vel_grad_HI<km/s/pc>'][i]*np.cos(thet+90*u.deg)*\
                (lons-np.nanmean((l0h,l1h)))*u.deg.to(u.rad)*fcat['distance<pc>'][i] + \
                fcat['vel_grad_HI<km/s/pc>'][i]*np.sin(thet+90*u.deg)*\
                (lats-np.nanmean((b0h,b1h)))*u.deg.to(u.rad)*fcat['distance<pc>'][i] 
    
    bperp_pc = bperp*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    res = model_v - mom1
    plt.scatter(bperp_pc,mom1,s=5e-2)
    plt.scatter(bperp_pc,model_v,s=5e-2)
    plt.scatter(bperp_pc,model_v-np.nanstd(model_v),s=5e-2)
    plt.scatter(bperp_pc,model_v+np.nanstd(model_v),s=5e-2)
    plt.title(fcat['name'][i])
    
    plt.xlabel('offset from rotation axis')
    plt.show()


fcat['radius_HI<pc>']


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
