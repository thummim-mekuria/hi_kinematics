#!/usr/bin/env python
# coding: utf-8

# from convo w pascale (fluid dynamics): anna lisa varri - globular clusters edinburgh looking at the analogy between them and fluid mechanics


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Ellipse
from matplotlib.colors import ListedColormap
from matplotlib import cm

import numpy as np
import pandas as pd

import datacube as dc
import fun as fns

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.modeling import models,fitting

from scipy.optimize import curve_fit

from gaussdecomp import spectrum,fitter

plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


new_lbs = np.array([
       [  0,  64,   5,  95], #0
       [  5,  95,   0, 155],
       [  0,  24,   0,  48], #2
       [  5,  64,   2, 100],
       [  4,  36,   0,  44], #4
       [  0,  55,   3,  42],
       [  3,  47,   0,  32], #6
       [  0,  75,   0,  80],
       [  0,  31,   0,  64], #8
       [  12, 70,   4,  27],
       [  0, 150,   0, 136], #10
       [  0,  33,   3,  50],
       [  0,  32,   0,  56], #12
       [  5,  38,   2,  66],
       [  0, 178,   0, 170], #14
       [  20, 60,  10, 120],
       [  5,  68,   5,  72], #16
       [  10, 100,   0, 116],
       [  0, 100,   0,  56], #18
       [  3,  70,   8,  75],
       [  0,  24,   0,  32], #20
       [  0,  75,   0, 102]])


h2_msks = [fits.open('../rn/h2_0thmom_masks_DISPLAY.fits')[i].data-1 for i in range(22)]


bos,bis,los,lis = new_lbs.transpose()
h2_msks = [h2_msks[i][bos[i]:bis[i],los[i]:lis[i]].astype('bool') for i in range(22)]


hi4pi,dht = fits.open('../rn/both_same_grid_720d.fits')


cat = QTable(fits.open('../rn/sighhplease.fits')[1].data)


bosc = bos + cat['boc']
bisc = bis + cat['boc']

losc = los + cat['loc']
lisc = lis + cat['loc']

coco = [dht.data[:,bosc[i]:bisc[i],losc[i]:lisc[i]] for i in range(22)]


[(
plt.contourf(np.where(~h2_msks[i],np.nan,np.log10(
    np.nansum(coco[i],axis=0)*1.3*1.97e20)),
             levels=20),
    
    
    plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


#the difference that the masks make
np.array([fns.mass_n_area(coco[i]*u.K,cat['distance<pc>'][i]*u.pc,'CO')[0].to(u.solMass).value for i in range(22)])/\
np.array([fns.mass_n_area(coco[i]*u.K,cat['distance<pc>'][i]*u.pc,'CO',h2_msks[i])[0].to(u.solMass).value for i in range(22)])


h2_mass = np.array([fns.mass_n_area(coco[i]*u.K,
                                    cat['distance<pc>'][i]*u.pc,
                                    'CO',h2_msks[i])[0].to(u.solMass).value for i in range(22)])*u.solMass


n_wnm = 0.57/u.cm**3

hi_rad = np.cbrt(0.75*h2_mass/(n_wnm*1.36*c.m_p*np.pi)).to(u.pc) #sphere


hi_rad_ds = np.arctan(hi_rad/(cat['distance<pc>']*u.pc)).to(u.deg).value


hi_rad_px = (hi_rad_ds/0.125).astype('int')


bc = ((np.nanmean(cat['lat_GMC<deg>'],axis=1) + 50)/0.125).astype('int')
lc = ((np.nanmean(cat['lon_GMC<deg>'],axis=1) - 360)/-0.125).astype('int')


bosh = bc - hi_rad_px
bish = bc + hi_rad_px

losh = lc - hi_rad_px
lish = lc + hi_rad_px


hihi = [hi4pi.data[:,bosh[i]:bish[i],losh[i]:lish[i]] for i in range(22)]


[(
plt.hist(np.ravel(
    np.log10(np.nansum(hihi[i],axis=0)*1.3*1.82e18)),
             bins=100),
    
    plt.title(cat['name'][i]),
    plt.yscale('log'),
    plt.vlines([20,21.2],0,1e3,color='red'),
    
#     plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


#filter out weird values for gemini,orion a and b, oph
hi_msks = [
    np.nansum(hihi[i],axis=0)*1.3*1.82e18 > 1e20 for i in range(22)]
hi_msks[7] = np.log10(np.nansum(hihi[7],axis=0)*1.3*1.82e18) > 21.2


[(
#     plt.contourf(hi_msks[i],cmap='magma'),
plt.contourf(np.log10(np.where(~hi_msks[i],np.nan,
    np.nansum(hihi[i],axis=0)*1.3*1.97e20)),
             levels=20,cmap='magma'),
    
    plt.title(cat['name'][i]),
    plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


r_arrs_pix = [np.meshgrid(np.arange(hi_rad_px[i]*-1, hi_rad_px[i]),
             np.arange(hi_rad_px[i]*-1, hi_rad_px[i]))
 
 for i in range(22)]


circles = [
    x**2 + y**2 < hi_rad_px[i]**2
    
    for i in range(22)
    for x,y in [r_arrs_pix[i]]
    
    
]


hi_msks = [hi_msks[i]*circles[i] for i in range(22)]


ppv_circles = [np.broadcast_to(circles[i],(155,circles[i].shape[0],circles[i].shape[1])) for i in range(22)]


sphihi = [hihi[i]*ppv_circles[i] for i in range(22)]


[(
#     plt.contourf(hi_msks[i],cmap='magma'),
plt.contourf(np.log10(np.where(~hi_msks[i],np.nan,
    np.nansum(sphihi[i],axis=0)*1.3*1.97e20)),
             levels=20,cmap='magma'),
    
    plt.title(cat['name'][i]),
    plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


[(
#     plt.contourf(hi_msks[i],cmap='magma'),
plt.contourf(np.log10(
    np.nansum(sphihi[i],axis=1)*0.125*1.82e18),
             levels=20,cmap='magma'),
    
    plt.title(cat['name'][i]),
    plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


for i in range(22):
    fig,axs = plt.subplots(1,3,figsize=(10,5),dpi=200)
    
    axs[0].contourf(np.log10(
            np.nansum(hi4pi.data[:,bosh[i]:bish[i],:],
                    axis=1)*0.125*1.82e18),levels=np.linspace(20,21.5,20),cmap='binary_r'),
    

    axs[0].contourf(np.log10(
            np.nansum(sphihi[i],axis=1)*0.125*1.82e18),
                     levels=np.linspace(20,21.5,20),
                    cmap='magma',extent=[losh[i],lish[i],0,155],
                        alpha=0.7),
    
    
    axs[0].contourf(np.log10(
            np.nansum(coco[i],axis=1)*0.125*1.97e20),
                     levels=np.linspace(20,21.5,20),
                    cmap='YlGnBu_r',extent=[losc[i],lisc[i],0,155],
                        alpha=0.7),
    
    
    axs[0].set_xlim(losh[i]-50,lish[i]+50)
    axs[0].hlines(vs_init[i],losh[i]-50,lish[i]+50)
    axs[0].vlines(lc[i],0,155)
    
    
    axs[1].contourf(np.log10(
                np.nansum(hi4pi.data[:,:,losh[i]:lish[i]],
                         axis=2)*0.125*1.82e18),levels=np.linspace(20,21.5,20),
                    cmap='binary_r'),
    
    axs[1].contourf(np.log10(
            np.nansum(sphihi[i],axis=2)*0.125*1.82e18),
             levels=np.linspace(20,21.5,20),
                    cmap='magma',extent=[bosh[i],bish[i],0,155],
                alpha=0.7),
    
    
    axs[1].contourf(np.log10(
            np.nansum(coco[i],axis=2)*0.125*1.97e20),
             levels=np.linspace(20,21.5,20),
                    cmap='YlGnBu_r',extent=[bosc[i],bisc[i],0,155],
                alpha=0.7),
    
    
    axs[1].set_xlim(bosh[i]-50,bish[i]+50),
    axs[1].hlines(vs_init[i],bosh[i]-50,bish[i]+50)
    axs[1].vlines(bc[i],0,155)
    
    axs[2].plot(np.nanmean(sphihi[i],axis=(1,2)),np.arange(155))
    axs[2].plot(70*np.nanmean(coco[i],axis=(1,2)),np.arange(155))
    axs[2].hlines(vs_init[i],0,np.nanmax(np.nanmean(sphihi[i],axis=(1,2))))
    
    fig.suptitle(cat['name'][i])
    
    plt.show()


vs_init = np.array([
    [70,100], #0
    [30,90],
    [65,85],
    [55,90], #3
    [83,100],
    [75,90],
    [70,95], #6
    [65,100],
    [70,110],
    [65,85], #9
    [55,95],
    [90,107],
    [65,90], #12
    [70,110],
    [60,98],
    [70,100], #15
    [65,110],
    [65,100],
    [65,90], #18
    [55,85],
    [70,100],
    [45,95] #21
])


fttr = fitting.LevMarLSQFitter()
init = models.Gaussian1D(1.,0.,1.,)


hi_sps = [  np.nanmean(sphihi[i][vs_init[i][0]:vs_init[i][1],:,:],
                    axis=(1,2))
        
for i in range(22)]


gausses = [fttr(init,np.arange(vs_init[i][1]-vs_init[i][0]),
                   hi_sps[i])
    
    for i in range(22)
]


co_center = [fns.get_mom1(np.nanmean(coco[i],axis=(1,2)),np.arange(155)) for i in range(22)]


for i in range(22):
    fig,axs = plt.subplots(1,3,figsize=(10,5),dpi=200)
    
    axs[0].contourf(np.log10(
            np.nansum(hi4pi.data[:,bosh[i]:bish[i],:],
                    axis=1)*0.125*1.82e18),levels=np.linspace(20,21.5,20),cmap='binary_r'),
    

    axs[0].contourf(np.log10(
            np.nansum(sphihi[i],axis=1)*0.125*1.82e18),
                     levels=np.linspace(20,21.5,20),
                    cmap='magma',extent=[losh[i],lish[i],0,155],
                        alpha=0.7),
    
    
    axs[0].contourf(np.log10(
            np.nansum(coco[i],axis=1)*0.125*1.97e20),
                     levels=np.linspace(20,21.5,20),
                    cmap='YlGnBu_r',extent=[losc[i],lisc[i],0,155],
                        alpha=0.7),
    
    
    axs[0].contour(np.log10(
            np.nansum(coco[i],axis=1)*0.125*1.97e20),
                     levels=[20],
                    colors='white',extent=[losc[i],lisc[i],0,155],
                        alpha=0.7),
    
    axs[0].set_xlim(losh[i]-50,lish[i]+50)
#     axs[0].hlines(vs_init[i],losh[i]-50,lish[i]+50)
    axs[0].hlines([co_center[i] - 2*gausses[i].stddev,
                   co_center[i] + 2*gausses[i].stddev],
                  losh[i]-50,lish[i]+50,color='red')
    
    axs[0].hlines([co_center[i] - gausses[i].stddev,
                   co_center[i] + gausses[i].stddev],
                  losh[i]-50,lish[i]+50,color='green')


    axs[0].vlines(lc[i],0,155)
    
    
    axs[1].contourf(np.log10(
                np.nansum(hi4pi.data[:,:,losh[i]:lish[i]],
                         axis=2)*0.125*1.82e18),levels=np.linspace(20,21.5,20),
                    cmap='binary_r'),
    
    axs[1].contourf(np.log10(
            np.nansum(sphihi[i],axis=2)*0.125*1.82e18),
             levels=np.linspace(20,21.5,20),
                    cmap='magma',extent=[bosh[i],bish[i],0,155],
                alpha=0.7),
    
    
    axs[1].contourf(np.log10(
            np.nansum(coco[i],axis=2)*0.125*1.97e20),
             levels=np.linspace(20,21.5,20),
                    cmap='YlGnBu_r',extent=[bosc[i],bisc[i],0,155],
                alpha=0.7),
    
    
    axs[1].contour(np.log10(
            np.nansum(coco[i],axis=2)*0.125*1.97e20),
             levels=[20],
                    colors='white',extent=[bosc[i],bisc[i],0,155],
                alpha=0.7),
    
    
    axs[1].set_xlim(bosh[i]-50,bish[i]+50),
#     axs[1].hlines(vs_init[i],bosh[i]-50,bish[i]+50)
    axs[1].vlines(bc[i],0,155)
    axs[1].hlines([co_center[i] - 2*gausses[i].stddev,
                   co_center[i] + 2*gausses[i].stddev],
                  bosh[i]-50,bish[i]+50,color='red')
    
    axs[1].hlines([co_center[i] - gausses[i].stddev,
                   co_center[i] + gausses[i].stddev],
                  bosh[i]-50,bish[i]+50,color='green')
    
    
    axs[2].plot(hi_sps[i],np.arange(vs_init[i][0],vs_init[i][1]))
    axs[2].plot(70*np.nanmean(coco[i],axis=(1,2)),np.arange(155))
#     axs[2].hlines(vs_init[i],0,np.nanmax(np.nanmean(hihi[i],axis=(1,2))))
    axs[2].hlines([co_center[i] - 2*gausses[i].stddev,
                   co_center[i] + 2*gausses[i].stddev],
                  0,np.nanmax(hi_sps[i]),color='red')
    
    axs[2].hlines([co_center[i] - gausses[i].stddev,
                   co_center[i] + gausses[i].stddev],
                  0,np.nanmax(hi_sps[i]),color='green')
    
    
    for ax in axs:
        ax.set_ylim(vs_init[i][0]+gausses[i].mean - 3.5*gausses[i].stddev,
                   vs_init[i][0]+gausses[i].mean + 3.5*gausses[i].stddev)
    
    fig.suptitle(cat['name'][i])
    
    plt.show()


vos,vis = np.array([np.array((co_center[i] - gausses[i].stddev,
                              co_center[i] + gausses[i].stddev)).astype('int')
for i in range(22)]).transpose()


long = np.array([fns.mass_n_area(coco[i]*u.K,cat['distance<pc>'][i]*u.pc,
                 'CO',mask=h2_msks[i])[0].to(u.solMass).value
        for i in range(22)])

short = np.array([fns.mass_n_area(coco[i][vos[i]:vis[i],:,:]*u.K,
                                  cat['distance<pc>'][i]*u.pc,
                 'CO',mask=h2_msks[i])[0].to(u.solMass).value
        for i in range(22)])


cat[(100*short/long)<90] #clouds whose mass changed by >10% 
# error in distance goes up by a factor of 2, and for all of them its 5 + sth else


#change those to be +/- 2sigma
vos[17],vis[17] = np.array([np.array((co_center[i] - 2*gausses[i].stddev,
                              co_center[i] + 2*gausses[i].stddev)).astype('int')
for i in [17]]).transpose()


vels = np.array((vos,vis)).transpose()*1.3-100


#differences for velocity gradient
[(i,100*(-1+
    np.linalg.norm(
    fns.fit_gradient(
    fns.get_mom1(coco[i][vos[i]:vis[i],:,:],np.arange(vos[i],vis[i])*1.3-100),
                  cat['distance<pc>'][i],mask=h2_msks[i],
weights = np.sqrt(fns.get_momN(coco[i][vos[i]:vis[i],:,:],
                               np.arange(vos[i],vis[i])*1.3-100,2)))[:1])/\
 
 np.linalg.norm(
    fns.fit_gradient(
    fns.get_mom1(coco[i],np.arange(155)*1.3-100),
                  cat['distance<pc>'][i],mask=h2_msks[i],
weights = np.sqrt(fns.get_momN(coco[i],np.arange(155)*1.3-100,2)))[:1])
))
 
 for i in range(22)]
#just gonna move past these & say that its not necessarily indicative of a problem
#that omega goes down when delta v goes down


[(
    plt.contourf(
        np.where(~hi_msks[i],np.nan,fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                              np.arange(vos[i],vis[i])*1.3-100)),
                 cmap=fns.segmented_spectral,levels=20),
    plt.title(i),
    plt.gca().set_aspect(1),
    plt.show()


)for i in range(22)]


i = 14

plt.hist(np.ravel(np.nansum(sphihi[i][vos[i]:vis[i],:,:],axis=0)*1.3*1.82e18*(
    np.abs(fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                              np.arange(vos[i],vis[i])*1.3-100))<8)),
         bins=10**np.arange(19,23,0.1))


plt.hist(np.ravel(np.nansum(sphihi[i][vos[i]:vis[i],:,:],axis=0)*1.3*1.82e18*(
    np.abs(fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                              np.arange(vos[i],vis[i])*1.3-100))>8)),
         bins=10**np.arange(19,23,0.1))

plt.vlines(3e20,0,1e4,color='red')

# plt.hist(np.ravel(np.nansum(hihi[i],axis=0)*1.3*1.82e18*(
#     np.abs(fns.get_mom1(hihi[i],
#                               np.arange(155)*1.3-100))>30)),
#          bins=10**np.arange(19,23,0.1))


# plt.hist(np.ravel(np.nansum(hihi[i],axis=0)*1.3*1.82e18*(
#     np.abs(fns.get_mom1(hihi[i],
#                               np.arange(155)*1.3-100))<30)),
#          bins=10**np.arange(19,23,0.1),alpha=0.5)

plt.xscale('log')
plt.yscale('log')


#redoing to account for new bad pixels
hi_msks[14] = hi_msks[14]*np.nansum(sphihi[14][vos[14]:vis[14],:,:],
                                 axis=0)*1.3*1.82e18 > 3e20


[(
    plt.contourf(
        np.where(~hi_msks[i],np.nan,fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                              np.arange(vos[i],vis[i])*1.3-100)),
                 cmap=fns.segmented_spectral,levels=20),
    plt.title(i),
    plt.gca().set_aspect(1),
    plt.show()


)for i in range(22)]


[(
    plt.scatter(np.ravel(np.indices(coco[i].shape[1:])[0]),
    np.ravel(np.where(~h2_msks[i],np.nan,
                      fns.get_mom1(coco[i][vos[i]:vis[i],:,:],
                                   np.arange(vos[i],vis[i])*1.3-100))),
           alpha=0.7,marker='+'),

plt.scatter(np.ravel(np.indices(coco[i].shape[1:])[0]),
    np.ravel(np.where(~h2_msks[i],np.nan,fns.get_mom1(coco[i],np.arange(155)*1.3-100))),
            alpha=0.7,s=5),
plt.title(i),

plt.show()

)for i in range(22)]


[(
    plt.scatter(np.ravel(np.indices(coco[i].shape[1:])[1]),
    np.ravel(np.where(~h2_msks[i],np.nan,
                      fns.get_mom1(coco[i][vos[i]:vis[i],:,:],
                                   np.arange(vos[i],vis[i])*1.3-100))),
           alpha=0.7,marker='+'),

plt.scatter(np.ravel(np.indices(coco[i].shape[1:])[1]),
    np.ravel(np.where(~h2_msks[i],np.nan,fns.get_mom1(coco[i],np.arange(155)*1.3-100))),
            alpha=0.7,s=5),
plt.title(i),

plt.show()

)for i in range(22)]


#not even sure how to interpret these? idk moving on and saying that
#the other emission (showing up in where the centroids are for the blue dots)
#is bc of the influence of the unrelated emission on moving the first moment value around
[(
    plt.scatter(np.ravel(np.indices(sphihi[i].shape[1:])[0]),
    np.ravel(np.where(~hi_msks[i],np.nan,
                      fns.get_mom1(sphihi[i],np.arange(155)*1.3-100),)),
           alpha=0.5,marker='+'),

plt.scatter(np.ravel(np.indices(sphihi[i].shape[1:])[0]),
    np.ravel(np.where(~hi_msks[i],np.nan, 
                      fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                                        np.arange(vos[i],vis[i])*1.3-100),)),
            alpha=0.5,s=5),
plt.title(i),

plt.show()

)for i in range(22)]


[(
    plt.scatter(np.ravel(np.indices(sphihi[i].shape[1:])[1]),
    np.ravel(np.where(~hi_msks[i],np.nan,
                      fns.get_mom1(sphihi[i],np.arange(155)*1.3-100),)),
           alpha=0.7,marker='+'),

plt.scatter(np.ravel(np.indices(sphihi[i].shape[1:])[1]),
    np.ravel(np.where(~hi_msks[i],np.nan, 
                      fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                                        np.arange(vos[i],vis[i])*1.3-100),)),
            alpha=0.7,s=5),
plt.title(i),

plt.show()

)for i in range(22)]


plt.plot(h2_omg_diffs,marker='.')
plt.plot(hi_omg_diffs,marker='.')

plt.hlines(1,0,22)

plt.xticks(np.arange(0,22,3))

plt.yscale('log')
plt.grid()


#okay i think i will just stick with this, esp bc the version with +/-100km/s isnt any more 'real'
#than the ones that are +/- 1 sigma of the fitted gaussian, its all just making choices about whether
#to err on the side of false positives or false negatives


[
    (
        plt.hist(np.ravel(np.where(~hi_msks[i],np.nan,np.nansum(
            sphihi[i][vos[i]:vis[i],:,:],axis=0)*1.3*1.82e18)),
                 bins=10**np.arange(20,22,0.05)),
        plt.xscale('log'),
        plt.title(cat['name'][i]),
        plt.show()
    )for i in range(22)]


hi_0lvls = np.log10(np.array(([2e20,9e20],
[4e20,5e21],
[3e20,8e20],
[2e20,9e20],
[7e20,2e21],
[1e20,6e20],
[8e20,2e21],
[1e21,4e21],
[6e20,3e21],
[5e20,2e21],
[2e20,2e21],
[1e21,4e21],
[5e20,3e21],
[7e20,2e21],
[3e20,3e21],
[2e20,2e21],
[4e20,2e21],
[8e20,3e21],
[2e20,2e21],
[1e20,8e20],
[1e21,4e21],
[4e20,3e21])))


[(
    plt.contourf(np.log10(np.nansum(sphihi[i][vos[i]:vis[i],:,:],axis=0)*1.3*1.82e18),
             extent=[losh[i],lish[i],bosh[i],bish[i]],
             levels=np.linspace(hi_0lvls[i][0],hi_0lvls[i][1],20), 
                 cmap='magma',extend='both'),
    plt.colorbar(),
    
    plt.contour(np.where(~h2_msks[i],np.nan,
            np.nansum(coco[i][vos[i]:vis[i],:,:],axis=0)),
             extent=[losc[i],lisc[i],bosc[i],bisc[i]],
             levels=10, colors='white'),
    
    
  plt.gca().set_aspect(1),
plt.show())

for i in range(22)]


hi_mom1s = [fns.get_mom1(sphihi[i][vos[i]:vis[i],:,:],
                              np.arange(vos[i],vis[i])*1.3-100) for i in range(22)]

hi_mom2s = [np.sqrt(fns.get_momN(sphihi[i][vos[i]:vis[i],:,:],
                              np.arange(vos[i],vis[i])*1.3-100,2)) for i in range(22)]


[(
    plt.contourf(
        np.where(~hi_msks[i],np.nan,hi_mom1s[i]),
                 cmap=fns.segmented_spectral,levels=20),
    plt.title(i),
    plt.gca().set_aspect(1),
    plt.show()


)for i in range(22)]


[(
    plt.contourf(
        np.where(~hi_msks[i],np.nan,hi_mom2s[i]),
                 cmap='YlGnBu_r',levels=20,
        extent=[losh[i],lish[i],bosh[i],bish[i]]),
    
    plt.colorbar(),
    
    plt.contour(np.where(~h2_msks[i],np.nan,np.nansum(coco[i][vos[i]:vis[i],:,:],axis=0)),
             extent=[losc[i],lisc[i],bosc[i],bisc[i]],
             levels=10, colors='red'),
    
    
    plt.title(i),
#     plt.vlines()
    plt.gca().set_aspect(1),
    plt.show()


)for i in range(22)]


bounds = pd.DataFrame(np.array([losc,lisc,bosc,bisc,
        losh,lish,bosh,bish,
        vos,vis
       ]).transpose(),columns=['loc','lic','boc','bic',
               'loh','lih','boh','bih',
               'vo','vi'])


bounds.to_csv('lbv_bounds_both.csv',index=False)


#for later
hih2 = [hihi[i][:,bosc[i]-bosh[i]:bisc[i]-bosh[i],losc[i]-losh[i]:lisc[i]-losh[i]]
    
    for i in range(22)
]


fns.make_fits_from_list([hi_msks[i]+1 for i in range(22)],'hi_masks.fits')
fns.make_fits_from_list([h2_msks[i]+1 for i in range(22)],'h2_masks.fits')

fns.make_fits_from_list(coco,'h2h2.fits')
fns.make_fits_from_list(sphihi,'hihi.fits')
fns.make_fits_from_list(hih2,'hih2.fits')
