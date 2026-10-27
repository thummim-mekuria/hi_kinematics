#!/usr/bin/env python
# coding: utf-8


import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.io import fits

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.cm as cm

import fun as fns
import datacube as dc


cat = QTable(fits.open('sighhplease.fits')[1].data)

c0s = fits.open('h2_0s.fits')
c1s = fits.open('h2_1s.fits')
c2s = fits.open('h2_2s.fits')

msks = fits.open('h2_0thmom_masks.fits')
msks = [(msks[i].data - 1).astype('bool') for i in range(22)]


h0s = fits.open('hi_0s.fits')
h1s = fits.open('hi_1s.fits')
h2s = fits.open('hi_2s.fits')


cm.YlGn_r


cm.YlGnBu_r


cm.cmaps_listed


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


names = ['Aquila South', 'California', 'Camelopardalis', 'Chamaeleon',
           'Canis Major OB1', 'Corona Australis', 'Crossbones', 'Gemini OB1',
           'Hercules', 'Lacerta', 'Lupus', 'Maddalena', 'Monoceros OB1', 'Monoceros R2',
           'Ophiuchus', 'Orion A', 'Orion B', r'$\lambda$-Orionis', 'Perseus', 'Polaris',
           'Rosette', 'Taurus']
    
lbls = ['AqS', 'Cal', 'Cam', 'Chm',
           'CMa', 'CrA', 'Xbn', 'Gem',
           'Her', 'Lac', 'Lup', 'Mad', 'MOB', 'MR2',
           'Oph', 'OrA', 'OrB', r'$\lambda$Or', 'Per', 'Pol',
           'Ros', 'Tau']


help(plt.colorbar)


cm.scale


allcmaps = ['Accent',
 'Accent_r',
 'Blues',
 'Blues_r',
 'BrBG',
 'BrBG_r',
 'BuGn',
 'BuGn_r',
 'BuPu',
 'BuPu_r',
 'CMRmap',
 'CMRmap_r',
 'Dark2',
 'Dark2_r',
 'GnBu',
 'GnBu_r',
 'Greens',
 'Greens_r',
 'Greys',
 'Greys_r',
 'OrRd',
 'OrRd_r',
 'Oranges',
 'Oranges_r',
 'PRGn',
 'PRGn_r',
 'Paired',
 'Paired_r',
 'Pastel1',
 'Pastel1_r',
 'Pastel2',
 'Pastel2_r',
 'PiYG',
 'PiYG_r',
 'PuBu',
 'PuBuGn',
 'PuBuGn_r',
 'PuBu_r',
 'PuOr',
 'PuOr_r',
 'PuRd',
 'PuRd_r',
 'Purples',
 'Purples_r',
 'RdBu',
 'RdBu_r',
 'RdGy',
 'RdGy_r',
 'RdPu',
 'RdPu_r',
 'RdYlBu',
 'RdYlBu_r',
 'RdYlGn',
 'RdYlGn_r',
 'Reds',
 'Reds_r',
 'Set1',
 'Set1_r',
 'Set2',
 'Set2_r',
 'Set3',
 'Set3_r',
 'Spectral',
 'Spectral_r',
 'Wistia',
 'Wistia_r',
 'YlGn',
 'YlGnBu',
 'YlGnBu_r',
 'YlGn_r',
 'YlOrBr',
 'YlOrBr_r',
 'YlOrRd',
 'YlOrRd_r',
 'afmhot',
 'afmhot_r',
 'autumn',
 'autumn_r',
 'binary',
 'binary_r',
 'bone',
 'bone_r',
 'brg',
 'brg_r',
 'bwr',
 'bwr_r',
 'cividis',
 'cividis_r',
 'cool',
 'cool_r',
 'coolwarm',
 'coolwarm_r',
 'copper',
 'copper_r',
 'cubehelix',
 'cubehelix_r',
 'gist_earth',
 'gist_earth_r',
 'gist_gray',
 'gist_gray_r',
 'gist_heat',
 'gist_heat_r',
 'gist_ncar',
 'gist_ncar_r',
 'gist_rainbow',
 'gist_rainbow_r',
 'gist_stern',
 'gist_stern_r',
 'gist_yarg',
 'gist_yarg_r',
 'gnuplot',
 'gnuplot2',
 'gnuplot2_r',
 'gnuplot_r',
 'gray',
 'gray_r',
 'hot',
 'hot_r',
 'hsv',
 'hsv_r',
 'inferno',
 'inferno_r',
 'jet',
 'jet_r',
 'magma',
 'magma_r',
 'nipy_spectral',
 'nipy_spectral_r',
 'ocean',
 'ocean_r',
 'pink',
 'pink_r',
 'plasma',
 'plasma_r',
 'prism',
 'prism_r',
 'rainbow',
 'rainbow_r',
 'seismic',
 'seismic_r',
 'spring',
 'spring_r',
 'summer',
 'summer_r',
 'tab10',
 'tab10_r',
 'tab20',
 'tab20_r',
 'tab20b',
 'tab20b_r',
 'tab20c',
 'tab20c_r',
 'terrain',
 'terrain_r',
 'turbo',
 'turbo_r',
 'twilight',
 'twilight_r',
 'twilight_shifted',
 'twilight_shifted_r',
 'viridis',
 'viridis_r',
 'winter',
 'winter_r']


for cmp in  allcmaps:
    i=0
#     plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(np.log10(h0s[i].data),
                cmap=cmp,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,15),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.3)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    plt.title(cmp)
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


for i in range(22):
    plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(np.log10(h0s[i].data),
                cmap=cm.magma,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,10),
                extent=[loc,lic,boc,bic],
                colors='white',linewidths=0.4,alpha=0.5)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


for i in range(22):
    plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(np.log10(h0s[i].data),
                cmap=cm.nipy_spectral_r,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,15),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.5)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


for i in range(22):
    plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(np.log10(h0s[i].data),
                cmap=cm.rainbow_r,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,15),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.5)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


for i in range(22):
    plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(np.log10(h0s[i].data),
                cmap=cm.tab20b_r,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,15),
                extent=[loc,lic,boc,bic],
                colors='white',linewidths=0.5)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


for i in range(22):
    plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(np.log10(h0s[i].data),
                cmap=cm.twilight_shifted_r,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,15),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.3)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


#maybe manually put in the ticks on the colorbars
#get the x and y labels a bit away
#maybe co should be in a diff color?
#save then do mom 1 maps (separately for those)


h0_lvls = np.array((
    [20.4,21.3],
    [20.8,21.8],
    [20.7,21.2],
    [20.7,21.3],
    [21.1,21.5],
    [20.5,21.],
    [21.1,21.5],
    [21.2,21.9],
    [20.7,21.5],
    [20.9,21.4],
    [20.7,22.],
    [20.7,21.6],
    [20.8,21.7],
    [21,21.6],
    [20.6,22.],
    [20,21.5],
    [20.7,21.8],
    [20.8,21.7],
    [20.6,21.4],
    [20.,21.2],
    [21.3,21.7],
    [20.7,21.8]))
    
    
[(plt.hist(np.ravel(np.log10(h0s[i].data)),bins=200,density=True),
plt.show()) for i in range(22)]


from copy import deepcopy


display_msks = deepcopy(msks)


for i in [17]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data,origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        (msky<10)+\
        (msky>93)+\
        ((mskx*0.3 + msky)>112)+\
        ((mskx*3.1 + msky)<75)+\
        ((mskx*0.7 + msky)<43)+\
        ((mskx*-0.7 + msky)<-50)+\
        ((mskx*-4.8+ msky)>30),


                origin='lower',cmap=dc.cmap_min0('binary_r'))#,alpha=0.5)

display_msks[i] = display_msks[i] * ~((msky<10)+\
        (msky>93)+\
        ((mskx*0.3 + msky)>112)+\
        ((mskx*3.1 + msky)<75)+\
        ((mskx*0.7 + msky)<43)+\
        ((mskx*-0.7 + msky)<-50)+\
        ((mskx*-4.8+ msky)>30))


for i in [17]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data,origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        (msky<10)+\
        (msky>93)+\
        ((mskx*0.3 + msky)>112)+\
        ((mskx*3.1 + msky)<75)+\
        ((mskx*0.7 + msky)<43)+\
        ((mskx*-0.7 + msky)<-50)+\
        ((mskx*-4.8+ msky)>30),


                origin='lower',cmap=dc.cmap_min0('binary_r'))#,alpha=0.5)

display_msks[i] = display_msks[i] * ~((msky<10)+\
        (msky>93)+\
        ((mskx*0.3 + msky)>112)+\
        ((mskx*3.1 + msky)<75)+\
        ((mskx*0.7 + msky)<43)+\
        ((mskx*-0.7 + msky)<-50)+\
        ((mskx*-4.8+ msky)>30))


from skimage.morphology import binary_erosion,erosion


for i in [0]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data,origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(erosion(display_msks[i]),origin='lower')
#     plt.imshow(
#         (msky<10)+\
#         (msky>93)+\
#         ((mskx*0.3 + msky)>112)+\
#         ((mskx*3.1 + msky)<75)+\
#         ((mskx*0.7 + msky)<43)+\
#         ((mskx*-0.7 + msky)<-50)+\
#         ((mskx*-4.8+ msky)>30),


#                 origin='lower',cmap=dc.cmap_min0('binary_r'))#,alpha=0.5)

display_msks[i] = erosion(display_msks[i])
#         (msky>93)+\
#         ((mskx*0.3 + msky)>112)+\
#         ((mskx*3.1 + msky)<75)+\
#         ((mskx*0.7 + msky)<43)+\
#         ((mskx*-0.7 + msky)<-50)+\
#         ((mskx*-4.8+ msky)>30))


for i in [19]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data*msks[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*-0.7 + msky)<-37)+\
        ((mskx*-1.7 + msky)<-109)+\
        ((mskx*0.7 + msky)<40)+\
        ((mskx*0.9 + msky)>120),

        
#         ((mskx*-4.8+ msky)>30),


                origin='lower',cmap=dc.cmap_min0('binary_r'))#,alpha=0.5)

display_msks[i] = display_msks[i] * \
~(
    ((mskx*-0.7 + msky)<-37)+\
    ((mskx*-1.7 + msky)<-109)+\
    ((mskx*0.7 + msky)<40)+\
    ((mskx*0.9 + msky)>120)
)


for i in [10]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data*msks[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(

        ((mskx*-0.75 + msky)<-70),

        
#         ((mskx*-4.8+ msky)>30),


#                 origin='lower',cmap=dc.cmap_min0('viridis'),alpha=0.5)
                origin='lower',cmap=dc.cmap_min0('binary_r'))

display_msks[i] = display_msks[i] * \
~(
     ((mskx*-0.75 + msky)<-70)
)


for i in [11]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data*msks[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(

        (mskx<5),

        
#         ((mskx*-4.8+ msky)>30),


                origin='lower',cmap=dc.cmap_min0('viridis'),alpha=0.5)
#                 origin='lower',cmap=dc.cmap_min0('binary_r'))

display_msks[i] = display_msks[i] * \
~(
     (mskx<5)
)


for i in [21]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data*msks[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        
        erosion(erosion(~(((mskx*2+ msky)>237)+\
        ((mskx*-1.1+ msky)<-76)))),

                origin='lower',cmap=dc.cmap_min0('viridis'),alpha=0.5)
#                 origin='lower',cmap=dc.cmap_min0('binary_r'))

display_msks[i] = display_msks[i] * erosion(erosion(~(((mskx*2+ msky)>237)+\
        ((mskx*-1.1+ msky)<-76))))


for i in [11]:
    msky,mskx = np.indices(display_msks[i].shape)
    plt.imshow(c0s[i].data*msks[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(

        (mskx<5),

        
#         ((mskx*-4.8+ msky)>30),


                origin='lower',cmap=dc.cmap_min0('viridis'),alpha=0.5)
#                 origin='lower',cmap=dc.cmap_min0('binary_r'))

display_msks[i] = display_msks[i] * \
~(
     (mskx<5)
)


for i in display_msks:
    plt.imshow(i,origin='lower')
    plt.show()


for n in names:
    print(n)


c0_big = fits.open('big_zeroth_co.fits')[0].data


display_msks[2] = display_msks[2]*(np.log10(c0s[2].data)>19.9)


for i in range(22):
    plt.figure(dpi=250)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    plt.contourf(np.log10(h0s[i].data),
                cmap=cm.magma,extend='both',
                levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
                 shrink=0.75,
                 label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
    
    
#     plt.contour(np.log10(h0s[i].data),
#                 colors='black',alpha=0.3,
#                 linestyle=':',
#                 levels=[20.6],linewidths=0.3,
#                 extent=[loh,lih,boh,bih])
    
    
    this = np.log10(c0s[i].data)*display_msks[i]
    smol,big = np.nanmin(this[this>18]),np.nanmax(this[this>18])


    plt.contour(np.log10(c0_big),
                levels=np.linspace(smol,big,6),
                extent=[360,-360,-50,50],
                colors='white',linewidths=0.2,alpha=0.4)
    
    
    plt.contour(np.log10(c0s[i].data)*display_msks[i],
                levels=np.linspace(smol,big,6),
                extent=[loc,lic,boc,bic],
                colors='white',linewidths=0.6,alpha=0.5)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    plt.title(names[i],loc='left',y=0.05,x=0.05,color='white',fontsize='small')
    
    plt.savefig('hi_density_map'+lbls[i],bbox_inches='tight')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


for cmp in  allcmaps:
    i=3
#     plt.figure(dpi=200)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    hi0th = plt.contourf(h1s[i].data,
                cmap=cmp,extend='both',
                levels=20,
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar()#ticks=np.arange(h0_lvls[i][0],h0_lvls[i][1],0.2),
#                  label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
                  
#     cbar.ax.locator_params(nbins=5)

    
#     plt.contour(np.log10(c0s[i].data),#*msks[i],
#                 levels=np.linspace(20,22.5,15),
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.3,alpha=0.3)
    
    plt.contour(np.log10(c0s[i].data)*msks[i],
                levels=np.linspace(20,22.5,15),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.3)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    plt.title(cmp)
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


cm.tab20b


cm.ocean


cm.gist_ncar_r


cm.gist_rainbow


cm.rainbow_r


cm.cubehelix


cm.Spectral


cm.coolwarm_r


cm.RdGy


cm.RdYlBu


cm.RdBu


from matplotlib.colors import ListedColormap


cm.tab20b


cm.tab20b.colors


#from cm.tab20b, take out a few things, rearrange the darks / lights 
#to be next to each other at the color transitions
#make the green and yellow darks less muddy

new_spectral = ListedColormap([
 (0.2235294117647059, 0.23137254901960785, 0.4745098039215686),
 (0.3215686274509804, 0.32941176470588235, 0.6392156862745098),
 (0.4196078431372549, 0.43137254901960786, 0.8117647058823529),
 (0.611764705882353, 0.6196078431372549, 0.8705882352941177),
 
 (0.807843137254902, 0.8588235294117647, 0.611764705882353),
 (0.7098039215686275, 0.8117647058823529, 0.4196078431372549),
 (0.5490196078431373, 0.7352941176470588, 0.3215686274509804),

 (0.8411764705882353, 0.6196078431372549, 0.2235294117647059),
 (0.9058823529411765, 0.7294117647058823, 0.3215686274509804),
 (0.9058823529411765, 0.796078431372549, 0.5803921568627451),
 
 (0.9058823529411765, 0.5882352941176471, 0.611764705882353),
 (0.8392156862745098, 0.3803921568627451, 0.4196078431372549),
 (0.6784313725490196, 0.28627450980392155, 0.2901960784313726),
 (0.5176470588235295, 0.23529411764705882, 0.2235294117647059),])

new_spectral_r = new_spectral.reversed()


new_spectral


[(plt.hist(np.ravel(h1s[i].data),bins=200,density=True),
plt.show()) for i in range(22)]


for i in [21]:#range(22):
#     plt.figure(dpi=250)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    plt.contourf(h1s[i].data,
                cmap=new_spectral_r,extend='both',
                levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],15),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.linspace(h1_lvls[i][0],h1_lvls[i][1],5),
                 shrink=0.75,)
#                  label=r'$\mathrm{log\ (\ N_{HI}\ /\ cm^{-2}\ )}$')
    
    
# #     plt.contour(np.log10(h0s[i].data),
# #                 colors='black',alpha=0.3,
# #                 linestyle=':',
# #                 levels=[20.6],linewidths=0.3,
# #                 extent=[loh,lih,boh,bih])
    
    
    this = np.log10(c0s[i].data)*display_msks[i]
    smol,big = np.nanmin(this[this>18.]),np.nanmax(this[this>18])


#     plt.contour(np.log10(c0_big),
#                 levels=np.linspace(smol,big,6),
#                 extent=[360,-360,-50,50],
#                 colors='white',linewidths=0.2,alpha=0.4)
    
    
    plt.contour(np.log10(c0s[i].data)*display_msks[i],
                levels=np.linspace(smol,big,6),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.6,alpha=0.5)
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    plt.title(names[i],loc='left',y=0.05,x=0.05,color='black',fontsize='small')
    
#     plt.savefig('zeroth_map'+lbls[i],bbox_inches='tight')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


h1_lvls = np.array((
    [1,13],
    [-27,4],
    [-6,-1],
    [-6,3],
    [14,20],
    [3,8],
    [7,14],
    [-1,14],
    [4,9],
    [-6,0],
    [-25,5],
    [27,31],
    [1,9],
    [6,25],
    [-7,12],
    [5,14],
    [5,25],
    [3,11],
    [-10.,7.],
    [-11,-1.],
    [11,15],
    [-21,9]))
    
    
for j in range(1,16):
    print(cm.Spectral(j/15))


seg_spec = ListedColormap([
(0.7633986928104575, 0.1633986928104575, 0.2928104575163399, 1.0),
    (0.7633986928104575, 0.1633986928104575, 0.2928104575163399, 1.0),

(0.8758169934640523, 0.3045751633986928, 0.29411764705882354, 1.0),
    (0.8758169934640523, 0.3045751633986928, 0.29411764705882354, 1.0),

(0.9568627450980393, 0.42745098039215684, 0.2627450980392157, 1.0),
    (0.9568627450980393, 0.42745098039215684, 0.2627450980392157, 1.0),

(0.9803921568627452, 0.5973856209150326, 0.34117647058823525, 1.0),
    (0.9803921568627452, 0.5973856209150326, 0.34117647058823525, 1.0),

(0.9934640522875817, 0.7477124183006535, 0.4352941176470587, 1.0),
    (0.9934640522875817, 0.7477124183006535, 0.4352941176470587, 1.0),

(0.996078431372549, 0.8784313725490196, 0.5450980392156862, 1.0),
    (0.996078431372549, 0.8784313725490196, 0.5450980392156862, 1.0),

(0.95,0.95,0.95),
(0.7477124183006538, 0.8980392156862746, 0.6274509803921569, 1.0),
    (0.7477124183006538, 0.8980392156862746, 0.6274509803921569, 1.0),

(0.5803921568627453, 0.8313725490196079, 0.6444444444444445, 1.0),
    (0.5803921568627453, 0.8313725490196079, 0.6444444444444445, 1.0),

(0.4, 0.7607843137254902, 0.6470588235294118, 1.0),
    (0.4, 0.7607843137254902, 0.6470588235294118, 1.0),

(0.26405228758169935, 0.6091503267973857, 0.7098039215686275, 1.0),
    (0.26405228758169935, 0.6091503267973857, 0.7098039215686275, 1.0),

(0.25359477124183005, 0.45882352941176474, 0.7058823529411765, 1.0),
    (0.25359477124183005, 0.45882352941176474, 0.7058823529411765, 1.0),

(0.3686274509803922, 0.30980392156862746, 0.6352941176470588, 1.0),
(0.3686274509803922, 0.30980392156862746, 0.6352941176470588, 1.0)])


seg_spec


for i in range(22):
    plt.figure(dpi=250)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    plt.contourf(h1s[i].data,
                cmap=seg_spec,extend='both',
                levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
                extent=[loh,lih,boh,bih])
    
    
    plt.colorbar(ticks=np.around(np.linspace(h1_lvls[i][0],h1_lvls[i][1],7),1),
                 shrink=0.75,
                 label=r'$\mathrm{\overline{v_{LSR}}\ \ [km\ s^{-1}]}$')
    
    
# #     plt.contour(np.log10(h0s[i].data),
# #                 colors='black',alpha=0.3,
# #                 linestyle=':',
# #                 levels=[20.6],linewidths=0.3,
# #                 extent=[loh,lih,boh,bih])
    
    
    this = np.log10(c0s[i].data)*display_msks[i]
    smol,big = np.nanmin(this[this>18.]),np.nanmax(this[this>18])


    plt.contour(np.log10(c0_big),
                levels=np.linspace(smol,big,6),
                extent=[360,-360,-50,50],
                colors='black',linewidths=0.2,alpha=0.4)
    
    
    plt.contour(np.log10(c0s[i].data)*display_msks[i],
                levels=np.linspace(smol,big,6),
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.6,alpha=0.5)
        
    
    plt.xlim(loh,lih)
    plt.ylim(boh,bih)
    
    plt.xticks(np.around(np.linspace(loh-0.5,lih+0.5,5),1))
    plt.yticks(np.around(np.linspace(boh+0.5,bih-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    plt.title(names[i],loc='left',y=0.05,x=0.05,color='black',fontsize='small')
    
    plt.savefig('hi_velocity_map'+lbls[i],bbox_inches='tight')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


[(plt.hist(np.ravel(np.where(~display_msks[i],np.nan,c1s[i].data)),bins=80,density=True),
plt.show()) for i in range(22)]


seg_spec = ListedColormap([
    (0.5633986928104575, 0.0033986928104575, 0.0928104575163399, 1.0),
(0.6633986928104575, 0.0633986928104575, 0.1928104575163399, 1.0),
    (0.7633986928104575, 0.1633986928104575, 0.2928104575163399, 1.0),

(0.8758169934640523, 0.3045751633986928, 0.29411764705882354, 1.0),
    (0.8758169934640523, 0.3045751633986928, 0.29411764705882354, 1.0),

(0.9568627450980393, 0.42745098039215684, 0.2627450980392157, 1.0),
    (0.9568627450980393, 0.42745098039215684, 0.2627450980392157, 1.0),

(0.9803921568627452, 0.5973856209150326, 0.34117647058823525, 1.0),
    (0.9803921568627452, 0.5973856209150326, 0.34117647058823525, 1.0),

(0.9934640522875817, 0.7477124183006535, 0.4352941176470587, 1.0),
    (0.9934640522875817, 0.7477124183006535, 0.4352941176470587, 1.0),

(0.996078431372549, 0.8784313725490196, 0.5450980392156862, 1.0),
    (0.996078431372549, 0.8784313725490196, 0.5450980392156862, 1.0),

(0.95,0.95,0.95),
(0.7477124183006538, 0.8980392156862746, 0.6274509803921569, 1.0),
    (0.7477124183006538, 0.8980392156862746, 0.6274509803921569, 1.0),

(0.5803921568627453, 0.8313725490196079, 0.6444444444444445, 1.0),
    (0.5803921568627453, 0.8313725490196079, 0.6444444444444445, 1.0),

(0.4, 0.7607843137254902, 0.6470588235294118, 1.0),
    (0.4, 0.7607843137254902, 0.6470588235294118, 1.0),

(0.26405228758169935, 0.6091503267973857, 0.7098039215686275, 1.0),
    (0.26405228758169935, 0.6091503267973857, 0.7098039215686275, 1.0),

(0.25359477124183005, 0.45882352941176474, 0.7058823529411765, 1.0),
    (0.25359477124183005, 0.45882352941176474, 0.7058823529411765, 1.0),

(0.3686274509803922, 0.30980392156862746, 0.6352941176470588, 1.0),
(0.2686274509803922, 0.20980392156862746, 0.5352941176470588, 1.0),
(0.1686274509803922, 0.10980392156862746, 0.4352941176470588, 1.0)])


seg_spec


for i in range(22):
#     plt.figure(dpi=250)
    
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360
    
    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
    plt.contourf(np.where(~display_msks[i],np.nan,c1s[i].data),
                cmap=seg_spec,extend='both',
                levels=np.linspace(c1_lvls[i][0],c1_lvls[i][1],28),
                extent=[loc,lic,boc,bic])
    
    
    plt.colorbar(ticks=np.around(np.linspace(c1_lvls[i][0],c1_lvls[i][1],7),1),
                 shrink=0.75,
                 label=r'$\mathrm{\overline{v_{LSR}}\ \ [km\ s^{-1}]}$')
    
    
    this = np.log10(c0s[i].data)*display_msks[i]
    smol,big = np.nanmin(this[this>18.]),np.nanmax(this[this>18])


#     plt.contour(np.log10(c0_big),
#                 levels=np.linspace(smol,big,6),
#                 extent=[360,-360,-50,50],
#                 colors='black',linewidths=0.2,alpha=0.4)
    
    plt.contour(np.log10(c0s[i].data)*display_msks[i],
                levels=[smol],
                extent=[loc,lic,boc,bic],
                colors='black',linewidths=0.6,alpha=0.5)
    
                  
    plt.xlim(loc,lic)
    
    if i == 15:
        plt.ylim(-21.9,bic)
    else:
        plt.ylim(boc,bic)
    
    if i ==2:
        plt.xticks(np.around(np.linspace(loc-0.5,lic+0.5,5),1))
        plt.yticks(np.around(np.linspace(boc+0.5,bic-0.5,3),1))
        
    elif i ==15:
        plt.xticks(np.around(np.linspace(loc-0.5,lic+0.5,5),1))
        plt.yticks(np.around(np.linspace(-21.9,bic-0.5,3),1))
        
    elif i in [6,9,18]:
        plt.xticks(np.around(np.linspace(loc-0.5,lic+0.5,3),1))
        plt.yticks(np.around(np.linspace(boc+0.5,bic-0.5,5),1))
        
    else:
        plt.xticks(np.around(np.linspace(loc-0.5,lic+0.5,5),1))
        plt.yticks(np.around(np.linspace(boc+0.5,bic-0.5,5),1))
    plt.gca().set_aspect('equal')
    
    plt.xlabel(r'Galactic Longitude [$\degree$]')
    plt.ylabel(r'Galactic Latitude [$\degree$]')
    
    
    if i in [0,2,3,5,7,10,14,16,20]:
        plt.title(names[i],loc='left',y=0.87,x=0.05,color='black',fontsize='small')

    else:
        plt.title(names[i],loc='left',y=0.05,x=0.05,color='black',fontsize='small')
    
    plt.savefig('h2_velocity_map'+lbls[i],bbox_inches='tight')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


his = fits.open('hi_hisized.fits')
his = [his[i].data[cat['vo'][i]:cat['vi'][i], :, :] for i in range(22)]


coshsz = [dht[cat['vo'][i]:cat['vi'][i],
        cat['boh'][i]:cat['bih'][i],
        cat['loh'][i]:cat['lih'][i]] for i in range(22)]


for i in range(22):
    n = his[i].data.shape[1]
    
    num_ptch = 10
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    coords = np.array((inds[:-1],inds[1:])).transpose()
    slcs = np.array([(i,j) for j in coords for i in coords])
    
    
    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=300)
    
    map_axs = np.flip(axs,axis=0).transpose().flat

        
    loh = cat['loh'][i]*-0.125+360
    lih = cat['lih'][i]*-0.125+360
    
    boh = cat['boh'][i]*0.125-50
    bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360

    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
#     background_ax = plt.axes([0, 0, 1, 1]) # create a dummy subplot for the background
# background_ax.set_zorder(-1) # set the background subplot behind the others
# background_ax.imshow(image, aspect='auto') # show the backgroud image
    
    
    bkg_map = fig.add_subplot(111)#add_axes([0, 0, 1, 1])
    bkg_map.set_zorder(-10)
    bkg_map.set_axis_off()
    bkg_map.contourf(h1s[i].data,cmap=seg_spec,extend='both',
                    extent=[loh,lih,boh,bih],
                    linewidths=0.4,zorder=-20,alpha=0.1,
                    levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26))
   
    
#     this = np.log10(c0s[i].data)*display_msks[i]
#     smol,big = np.nanmin(this[this>18.]),np.nanmax(this[this>18])


#     bkg_map.contour(np.log10(c0s[i].data)*display_msks[i],
#                 levels=[smol],
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.6,alpha=0.3)
    
#     bkg_map.invert_xaxis()
    bkg_map.set_xlim(loh,lih)
    bkg_map.set_ylim(boh,bih)
#     bkg_map.set_aspect('equal')
        
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(his[i][:,bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(coshsz[i][:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.4,
                        linewidth=2,color='black')
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
        map_axs[j].vlines(np.nanmedian(h1s[i].data),
                          0,1.5,color='black',
                          linestyle='--',linewidth=1.5,alpha=0.2)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(his[i][:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
                           s=7,color='black',alpha=0.3,linewidth=0)
        
        
        if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
            map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                            (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')
            
            
            map_axs[j].scatter(np.nanmean(fns.get_mom1(coshsz[i][:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,s=1,color='black')
            map_axs[j].vlines(np.nanmedian(np.where(~msks[i],np.nan,c1s[i].data)),
                              0,1.5,color='black',linestyle='--',linewidth=0.4)
        
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_ylim(-0.15,1.15)
        map_axs[j].set_axis_off()
    
    
    plt.subplots_adjust(wspace=0.07,hspace=0.07)
    plt.figaspect(1)
        
    fig.suptitle(names[i],y=0.93)
    plt.savefig('patch_spec_'+lbls[i])

    plt.show()


#and for each og these (maybe this is for over the weekend or later), it will be impactful to select
#the spectra along the line of the velocity gradient (ie the line that foes from smallest v to largest v)
#and show how the spectra are changing / how the shifts are a result of additional components
#maybe decompose those? idk


for i in range(22):
    print(cat['name'][i][:4],
          round(np.nanmean(h1s[i].data),4),
          round(np.nanmean(np.where(~msks[i],np.nan,c1s[i].data)),4))


sp2


plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                            (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')


v0s = np.nanmin((h1_lvls.transpose()[0],c1_lvls.transpose()[0]),axis=0)
v1s = np.nanmax((h1_lvls.transpose()[1],c1_lvls.transpose()[1]),axis=0)


v0s = np.array([  0., -30.,  -7.,  -8.,  10.,   2.,   6.,  -10.,   2.,  -8., -32.,
        19.,   -10.,   6.,  -14.,   0.,   3.,  -5., -12., -13.,   5., -25.])


v1s = np.array([15. ,  10. ,  2,  7, 22. ,  10. , 17. , 23. , 17. ,  5. ,  14. ,
       35. , 13. , 27. , 19. , 17. , 30. , 17. , 12. , 4. , 24. ,  10. ])


cm.magma.resampled(20)(0.4)


for i in range(22):
    bns = np.arange(v0s[i],v1s[i],0.26)
    
    hi_hist = np.histogram(np.ravel(h1s[i].data),bins=bns)[0]
    h2_hist = np.histogram(np.ravel(np.where(~display_msks[i],np.nan,c1s[i].data)),bins=bns)[0]
    
    plt.step(bns[:-1],hi_hist/np.nanmax(hi_hist),
             alpha=0.3,color=(0.581819, 0.171596, 0.502777, 1.0),
            linewidth=4),
    
    
    plt.step(bns[:-1],h2_hist/np.nanmax(h2_hist),
            color='black'),


    plt.title(str(i)),

    plt.show()


np.arange(22)


fig,axs = plt.subplots(2,3,dpi=200)

for i,ax in zip([21,3,2,10,9,19,21],axs.flat):
    
    bns = np.arange(v0s[i],v1s[i],0.26)
    
    hi_hist = np.histogram(np.ravel(h1s[i].data),bins=bns)[0]
    h2_hist = np.histogram(np.ravel(np.where(~display_msks[i],np.nan,c1s[i].data)),bins=bns)[0]
    
    ax.step(bns[:-1],hi_hist/np.nanmax(hi_hist),
             alpha=0.3,color=(0.581819, 0.171596, 0.502777, 1.0),
            linewidth=4),
    
    
    ax.step(bns[:-1],h2_hist/np.nanmax(h2_hist),
            color='black'),

    ax.set_yticks([])
    ax.set_title(str(i),y=0.85),
    
axs[0,0].set_yticks([0.,0.5,1.])
axs[1,0].set_yticks([0.,0.5,1.])


fig.tight_layout()
plt.show()


dht = fits.open('both_same_grid_720d.fits')[1].data


for i in range(22):
    c1big = fns.get_mom1(dht[cat['vo'][i]:cat['vi'][i],
              cat['boh'][i]:cat['bih'][i],
              cat['loh'][i]:cat['lih'][i]],
              np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    plt.contourf(c1big)
    plt.show()


c1_lvls = np.array((
    [1,5],
    [-9,7],
    [-3,0.5],
    [1,5.5],
    [12,20],
    [4.5,6.5],
    [10,16],
    [-7,16],
    [5,15],
    [-4,4],
    [3,9],
    [22,33],
    [0,11],
    [8,14],
    [-1,9],
    [3,13],
    [4,12],
    [-4,15],
    [-1.,11.],
    [-6,-1.],
    [7,19],
    [4,8.5]))
    
    
# # for report to nia
# - change the lat range of orion a (bring the bottom up? ughh but then that will change the hi too, so maybe not
# - add in the background faded stuff
# - histograms (also colorcodded like before)
# - scatterplots (at least one?)
# - send
