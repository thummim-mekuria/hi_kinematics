#!/usr/bin/env python
# coding: utf-8


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

from scipy.optimize import curve_fit


hi4pi,dht = fits.open('../rn/both_same_grid_720d.fits')


from gaussdecomp import spectrum,fitter


anti_center_negb = hi4pi.data[:,:500,1240:1840]


anti_center_negb_co = dht.data[:,:500,1240:1840]


anti_sp = spectrum.Spectrum(np.nanmean(anti_center_negb,axis=(1,2)),np.arange(155)*1.3-100)


comps = fitter.gaussfitter(anti_sp,silent=True)['par'].reshape(-1,3)


from astropy.modeling import models


plt.plot(np.nanmean(anti_center_negb,axis=(1,2)))
[plt.plot(models.Gaussian1D(comps[i][0],
                            comps[i][1],
                            comps[i][2])(np.arange(155)*1.3-100)) 
 for i in range(len(comps))]

plt.plot(np.nansum([models.Gaussian1D(comps[i][0],
                            comps[i][1],
                            comps[i][2])(np.arange(155)*1.3-100)
 for i in range(len(comps))],axis=0))

plt.show()


anti_lv = np.nansum(anti_center_negb,axis=1)*0.125*1.82e18
anti_bv = np.nansum(anti_center_negb,axis=2)*0.125*1.82e18


fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5))

axs[0].contourf(np.log10(anti_lv),levels=np.linspace(20,22.,15))
axs[1].contourf(np.log10(anti_bv),levels=np.linspace(20,22.,15))


# axs[2].plot(np.nanmean(anti_center_negb,axis=(1,2)),np.arange(155)*1.3-100)

[axs[2].plot(models.Gaussian1D(comps[i][0],
                            comps[i][1],
                            comps[i][2])(np.arange(155)*1.3-100),np.arange(155)*1.3-100) 
 for i in range(len(comps))]

axs[2].plot(100*np.nanmean(anti_center_negb_co,axis=(1,2)),np.arange(155)*1.3-100)
# axs[2].plot(np.nansum([models.Gaussian1D(comps[i][0],
#                             comps[i][1],
#                             comps[i][2])(np.arange(155)*1.3-100)
#  for i in range(len(comps))],axis=0),np.arange(155)*1.3-100)

plt.show()


fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5))

axs[0].contourf(np.log10(anti_lv),levels=np.linspace(20,22.,15))
axs[1].contourf(np.log10(anti_bv),levels=np.linspace(20,22.,15))

for ax in axs[:2]:
    ax.set_ylim(70-20,90+20)

# axs[2].plot(np.nanmean(anti_center_negb,axis=(1,2)),np.arange(155)*1.3-100)

[axs[2].plot(models.Gaussian1D(comps[i][0],
                            comps[i][1],
                            comps[i][2])(np.arange(155)*1.3-100),np.arange(155)*1.3-100) 
 for i in range(len(comps))]

axs[2].plot(100*np.nanmean(anti_center_negb_co,axis=(1,2)),np.arange(155)*1.3-100)
axs[2].set_ylim(-4.87-20*1.3,16.16+20*1.3)
# axs[2].plot(np.nansum([models.Gaussian1D(comps[i][0],
#                             comps[i][1],
#                             comps[i][2])(np.arange(155)*1.3-100)
#  for i in range(len(comps))],axis=0),np.arange(155)*1.3-100)

plt.show()


comps[2][1]-2*comps[2][2],comps[2][1]+2*comps[2][2]


(np.array((comps[2][1]-2*comps[2][2],comps[2][1]+2*comps[2][2]))+100)/1.3


anti_smolv = anti_center_negb[73:90,:,:]


anti_smolv_co = anti_center_negb_co[73:90,:,:]


from astrodendro import Dendrogram


np.nanmin(anti_smolv),np.nanstd(anti_smolv),np.nanmean(anti_smolv),np.nanmax(anti_smolv)


for i in range(17):
    plt.contourf(anti_smolv[i],levels=np.linspace(20,70,10))
#     plt.contourf(anti_smolv_co[i],levels=np.linspace(0.5,5,10),cmap='magma')
    plt.title((i+73)*1.3-100)
    plt.show()


plt.contourf(np.nansum(anti_smolv,axis=0),levels=np.linspace(100,1000,15))


def get_mom1(cube,vels,dx=1.3):
    xTx = np.array([cube[i]*np.nanmean(cube[i])*vels[i] for i in range(len(vels))])
    
    sigxTx = np.nansum(xTx*dx,axis=0)

    M0 = np.nansum(cube*dx,axis=0)
    return sigxTx/M0


new_h1 = get_mom1(hi4pi.data,np.arange(155)*1.3-100)


plt.figure(dpi=200,figsize=(10,5))
plt.contourf(new_h1,
             cmap=fns.segmented_spectral,levels=np.linspace(-100,100,30),#np.linspace(-10,2,30),
            extent=[360,-360,-50,50])

plt.xlim(270,90)
plt.colorbar()
# plt.contour(np.nanmedian(anti_smolv,axis=0),
#              levels=np.linspace(10,80,5),
#             antialiased=True,colors='black',linewidths=np.linspace(0,1,5))


plt.figure(dpi=200,figsize=(10,5))
plt.contourf(new_h1,
             cmap=fns.segmented_spectral,levels=np.linspace(-100,100,30),#np.linspace(-10,2,30),
            extent=[360,-360,-50,50])

plt.xlim(90,-90)
plt.colorbar()
# plt.contour(np.nanmedian(anti_smolv,axis=0),
#              levels=np.linspace(10,80,5),
#             antialiased=True,colors='black',linewidths=np.linspace(0,1,5))


hihi = fits.open('../rn/hi_hisized.fits')


#hmm okay maybe the direction to go is in to use the second moment
#to determine which locations are most likely to be have multiple gaussian 
#components ie multiple emission features (then use the linewidth
#of the component we chose from the overall spectra to be 
#some sort of upper limit on the linewidths at individual pixels?)
#then maybe use something like the skewness to apply a correction
#factor (ie if it got pulled too far to the blue side, 
#then we see that it is skewed right, so the value of the first moment
#gets bumped redward)


[(plt.contourf(np.sqrt(fns.get_momN(hihi[i].data[oldcat['vo'][i]:oldcat['vi'][i],:,:],
                        np.arange(oldcat['vo'][i],oldcat['vi'][i])*1.3-100,2)),
              levels=30,cmap='YlGnBu_r'),
plt.colorbar(),
  plt.title(oldcat['name'][i]),
plt.show())
for i in range(22)]


[(plt.contourf(get_mom1(hihi[i].data[oldcat['vo'][i]:oldcat['vi'][i],:,:],
                        np.arange(oldcat['vo'][i],oldcat['vi'][i])*1.3-100),
              levels=30,cmap=fns.segmented_spectral),
plt.colorbar(),
  plt.title(oldcat['name'][i]),
plt.show())
for i in range(22)]


plt.contourf(get_mom1(anti_smolv,np.arange(73,90)*1.3-100),
             cmap=fns.segmented_spectral,levels=np.linspace(0,10,30))
plt.colorbar()
# plt.contour(np.nanmedian(anti_smolv,axis=0),
#              levels=np.linspace(10,80,5),
#             antialiased=True,colors='black',linewidths=np.linspace(0,1,5))


plt.contourf(get_mom1(anti_smolv,np.arange(73,90)*1.3-100),
             cmap=fns.segmented_spectral,levels=30,)
plt.colorbar()
# plt.contour(np.nanmedian(anti_smolv,axis=0),
#              levels=np.linspace(10,80,5),
#             antialiased=True,colors='black',linewidths=np.linspace(0,1,5))


plt.contourf(get_mom1(anti_smolv,np.arange(73,90)*1.3-100),
             cmap=fns.segmented_spectral,levels=30,)
plt.colorbar()
# plt.contour(np.nanmedian(anti_smolv,axis=0),
#              levels=np.linspace(10,80,5),
#             antialiased=True,colors='black',linewidths=np.linspace(0,1,5))


plt.contourf(np.nanmedian(anti_smolv,axis=0),levels=np.linspace(10,100,15))


oldcat = QTable(fits.open('../rn/sighhplease.fits')[1].data)


del_gmc_gal = oldcat['rot_ax_GMC<deg>']-90

for i in range(22):
    if del_gmc_gal[i] > 180:
        del_gmc_gal[i] -=360
        
    if del_gmc_gal[i] <- 180:
        del_gmc_gal[i] +=360


del_hi_gal = oldcat['rot_ax_HI<deg>']-90

for i in range(22):
    if del_hi_gal[i] > 180:
        del_hi_gal[i] -=360
        
    if del_hi_gal[i] <- 180:
        del_hi_gal[i] +=360


plt.hist(del_gmc_gal,bins=np.arange(-180,210,15),alpha=0.5,edgecolor='black')
plt.hist(del_hi_gal,bins=np.arange(-180,210,15),alpha=0.5,edgecolor='black')

plt.vlines([-60,60,-120,120],0,3)


plt.scatter(oldcat['radius_GMC<pc>'],oldcat['j_GMC<cm2/s>'],s=1e-1)
plt.scatter(oldcat['radius_HI<pc>'],oldcat['j_HI<cm2/s>'],s=1e-1)

[plt.plot([oldcat['radius_GMC<pc>'][i],oldcat['radius_HI<pc>'][i]],
         [oldcat['j_GMC<cm2/s>'][i],oldcat['j_HI<cm2/s>'][i]],color='black',linewidth=0.3)

for i in range(22)]
plt.loglog()


plt.scatter(oldcat['radius_GMC<pc>']/oldcat['radius_HI<pc>'],
         oldcat['mass_GMC<Msun>']/oldcat['mass_HI<Msun>'])
plt.loglog()


plt.scatter(oldcat['radius_GMC<pc>']/oldcat['radius_HI<pc>'],
         oldcat['j_GMC<cm2/s>']/oldcat['j_HI<cm2/s>'])
plt.loglog()


#this is the most correct version of the code for making hi pbp spectra
# for i in range(22):
n = anti_smolv.shape[1]

num_ptch = 10

inds = np.linspace(0,n,num_ptch+1).astype('int')
coords = np.array((inds[:-1],inds[1:])).transpose()
slcs = np.array([(i,j) for j in coords for i in coords])


fig,axs = plt.subplots(num_ptch,num_ptch,
                       sharex=True,sharey=True,
                      dpi=200)

map_axs = np.flip(axs,axis=0).flat

        
#     loh = cat['loh'][i]*-0.125+360
#     lih = cat['lih'][i]*-0.125+360
    
#     boh = cat['boh'][i]*0.125-50
#     bih = cat['bih'][i]*0.125-50
    
    
#     loc = cat['loc'][i]*-0.125+360
#     lic = cat['lic'][i]*-0.125+360

#     boc = cat['boc'][i]*0.125-50
#     bic = cat['bic'][i]*0.125-50
    
    
#     background_ax = plt.axes([0, 0, 1, 1]) # create a dummy subplot for the background
# background_ax.set_zorder(-1) # set the background subplot behind the others
# background_ax.imshow(image, aspect='auto') # show the backgroud image
    
    
#     bkg_map = fig.add_subplot(111)#add_axes([0, 0, 1, 1])
#     bkg_map.set_zorder(-10)
#     # bkg_map.set_axis_off()
#     bkg_map.contourf(hi_1s[i].data,cmap=fns.segmented_spectral,extend='both',
#                     extent=[loh,lih,boh,bih],
#                     zorder=-20,alpha=0.4,
#                     levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
#                     antialiased=True)
   
    
# #     this = np.log10(c0s[i].data)*display_msks[i]
# #     smol,big = np.nanmin(this[this>18.]),np.nanmax(this[this>18])


# #     bkg_map.contour(np.log10(c0s[i].data)*display_msks[i],
# #                 levels=[smol],
# #                 extent=[loc,lic,boc,bic],
# #                 colors='black',linewidths=0.6,alpha=0.3)
    
# #     bkg_map.invert_xaxis()
#     bkg_map.set_xlim(loh,lih)
#     bkg_map.set_ylim(boh,bih)
# #     bkg_map.set_aspect('equal')
        
    
for j in range(len(slcs)):
    bo,bi = slcs[j][1]
    lo,li = slcs[j][0]

    sp = np.nanmean(anti_smolv[:,bo:bi,lo:li],axis=(1,2))
    # sp2 = np.nanmean(coshsz[i][:,bo:bi,lo:li],axis=(1,2))

    map_axs[j].plot(np.arange(73,90)*1.3-100,
                    (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.6,
                    linewidth=1.5,color='black')
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
#         map_axs[j].vlines(np.nanmedian(hi_1s[i].data),
#                           0,1.5,color='black',
#                           linestyle='--',linewidth=1.,alpha=0.4)
        
#         map_axs[j].scatter(np.nanmean(fns.get_mom1(hihi[i].data[:,bo:bi,lo:li],
#                                    np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
#                            s=7,color='black',alpha=0.5,linewidth=0)
        
        
#         if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
#             map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
#                             (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')
            
            
#             map_axs[j].scatter(np.nanmean(fns.get_mom1(coshsz[i][:,bo:bi,lo:li],
#                                    np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,s=1,color='black')
#             map_axs[j].vlines(np.nanmedian(np.where(~msks[i],np.nan,c1s[i].data)),
                              # 0,1.5,color='black',linestyle='--',linewidth=0.4)
        
        
    map_axs[j].set_xticks([]),
    map_axs[j].set_yticks([])
    map_axs[j].set_ylim(-0.05,1.2)
    map_axs[j].set_axis_off()


plt.subplots_adjust(wspace=0.07,hspace=0.07)
plt.figaspect(1)

# fig.suptitle(names[i],y=0.93)
fig.supxlabel('Galactic Longitude')
fig.supylabel('Galactic Latitude')

plt.show()


fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5))

axs[0].contourf(np.log10(anti_lv),levels=np.linspace(20,22.,15))
axs[1].contourf(np.log10(anti_bv),levels=np.linspace(20,22.,15))

for ax in axs[:2]:
    ax.set_ylim(90,70)

# axs[2].plot(np.nanmean(anti_center_negb,axis=(1,2)),np.arange(155)*1.3-100)

[axs[2].plot(models.Gaussian1D(comps[i][0],
                            comps[i][1],
                            comps[i][2])(np.arange(155)*1.3-100),np.arange(155)*1.3-100) 
 for i in range(len(comps))]

axs[2].plot(100*np.nanmean(anti_center_negb_co,axis=(1,2)),np.arange(155)*1.3-100)
axs[2].set_ylim(16.16,-4.87,)
# axs[2].plot(np.nansum([models.Gaussian1D(comps[i][0],
#                             comps[i][1],
#                             comps[i][2])(np.arange(155)*1.3-100)
#  for i in range(len(comps))],axis=0),np.arange(155)*1.3-100)

plt.show()
