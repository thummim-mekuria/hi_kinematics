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


from scipy.optimize import curve_fit as fit
def line(x,m,b):
    return m*x + b


def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
    return vo + a*(x-xo) + b*(y-yo)


cat = QTable(fits.open('sighhplease.fits')[1].data)

c0s = fits.open('h2_0s.fits')
c1s = fits.open('h2_1s.fits')
c2s = fits.open('h2_2s.fits')

msks = fits.open('h2_0thmom_masks.fits')
msks = [(msks[i].data - 1).astype('bool') for i in range(22)]


h0s = fits.open('hi_0s.fits')
h1s = fits.open('hi_1s.fits')
h2s = fits.open('hi_2s.fits')


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
    
#     plt.savefig('hi_density_map'+lbls[i],bbox_inches='tight')
    
    #     plt.contourf(c0s[i].data*msks[i])
    plt.show()


from matplotlib.colors import ListedColormap


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
    
#     plt.savefig('hi_velocity_map'+lbls[i],bbox_inches='tight')
    
    #     plt.contourf(c0s[i].data*msks[i])
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
    
#     plt.savefig('h2_velocity_map'+lbls[i],bbox_inches='tight')
    
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


v0s = np.array([  0., -30.,  -7.,  -8.,  10.,   2.,   6.,  -10.,   2.,  -8., -32.,
        19.,   -10.,   6.,  -14.,   0.,   3.,  -5., -12., -13.,   5., -25.])


v1s = np.array([15. ,  10. ,  2,  7, 22. ,  10. , 17. , 23. , 17. ,  5. ,  14. ,
       35. , 13. , 27. , 19. , 17. , 30. , 17. , 12. , 4. , 24. ,  10. ])


for i in range(22):
    bns = np.arange(v0s[i],v1s[i],(v1s[i]-v0s[i])/70)
    
    hi_hist = np.histogram(np.ravel(h1s[i].data),bins=bns)[0]
    h2_hist = np.histogram(np.ravel(np.where(~display_msks[i],np.nan,c1s[i].data)),bins=bns)[0]
    
    plt.step(bns[:-1],hi_hist/np.nanmax(hi_hist),
             alpha=0.3,color=(0.581819, 0.171596, 0.502777, 1.0),
            linewidth=4),
    
    
    plt.step(bns[:-1],h2_hist/np.nanmax(h2_hist),
            color='black'),


    plt.title(cat['name'][i]),
#     plt.yscale('log'),

    plt.show()


[ 0, 
  
 2,  3, 9,10, 18, 19,21,14, 1,
 4,  5,  8,  12, 13, 16,20, 6,  7, 
 
 
 11, 15,17,
 
         ]

#trying to sort them based on similarities .. so far what i have is:
#co to the right, co around hi, and overlapping (0 is the only one with co to the left)


# # for report to nia
# - change the lat range of orion a (bring the bottom up? ughh but then that will change the hi too, so maybe not
# - add in the background faded stuff
# - histograms (also colorcodded like before)
# - scatterplots (at least one?)
# - send


cat


delphi = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']


for i in range(22):
    if delphi[i]>180:
        delphi[i] -=360
        
    if delphi[i]<-180:
        delphi[i] +=360


seg_spec


len(seg_spec.colors)


seg_spec.colors[16]


plt.plot(np.arange(-180,180),np.cos(u.deg*np.arange(-180,180)))
plt.vlines([-120,120,-60,60],-1,1)


plt.figure(dpi=200)

plt.hist(np.abs(delphi),bins=np.arange(0,210,15),
         color=seg_spec.colors[19],alpha=0.8,
         edgecolor='white')

plt.xticks(np.arange(0,210,60),
          [str(dg)+'$\degree$' for dg in np.arange(0,210,60)])

plt.vlines([60,120],0,5.3,colors='black',linestyle='--',linewidth=0.6)

plt.annotate('co-rotating',      (10,5.2))
plt.annotate('un-aligned',       (73,5.2))
plt.annotate('counter-rotating',(127,5.2))

plt.ylim(0,5.7)
plt.xlim(-5,185)
# plt.xlim(0,180)

plt.xlabel('$\mathrm{\Delta} \phi$')

plt.show()


#copied over from the notebook calculating Omega etc
dvdx_gmc = np.array([-0.00123186,  0.07228955, -0.01784782, -0.09577493, -0.01613217,
       -0.06095813,  0.07731865,  0.03445932, -0.24994185,  0.10428457,
       -0.03687923,  0.0389462 ,  0.01370792, -0.03489532,  0.03540596,
       -0.03361667, -0.01717165,  0.05740992,  0.29453495,  0.06868333,
       -0.03902938, -0.06541992])

dvdy_gmc = np.array([ 7.00523379e-02,  1.28886477e-01, -1.20493475e-02,  4.22427693e-02,
        2.90897919e-02, -1.73452040e-04,  4.10811443e-02, -2.27970079e-02,
       -2.85090746e-01,  6.10347562e-02,  2.54940335e-02,  7.38012656e-02,
        1.58945534e-01,  5.25524728e-03, -7.95726341e-02,  2.78899325e-02,
        1.25034889e-02,  1.72436946e-01,  1.75513630e-01,  4.64151366e-02,
       -1.06490156e-01, -1.26040411e-01])


gmc_mags = np.linalg.norm((dvdx_gmc,dvdy_gmc),axis=0)


bigh1 = fits.open('bigh1.fits')[0].data


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-30,30,30),extend='both')
plt.colorbar()
plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle')

plt.xlim(90,-90)
plt.title('gmc')


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-40,40,30),extend='both')
plt.colorbar()
plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle')

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)-360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle')


plt.xlim(-90,-270)
plt.ylim(-30,30)
plt.title('gmc')


dvdx_hi = np.array([-0.05570644,  0.04362601,  0.00269225, -0.02024164,  0.00152525,
         0.00267153, -0.00762758,  0.02532809, -0.00137701, -0.03229255,
         0.10315655, -0.00043652,  0.0088281 ,  0.00716635,  0.08752043,
        -0.00020426,  0.03466059,  0.03927205,  0.03635225,  0.01557406,
         0.00312602,  0.13374581])
                   
dvdy_hi = np.array([ 0.09068724, -0.09058775,  0.04946301, -0.0104879 ,  0.0192844 ,
        -0.03784893,  0.04268686, -0.00606768, -0.04949052,  0.01020186,
         0.12123622,  0.00082334, -0.02239394,  0.07330791, -0.03130099,
         0.01721732,  0.0882288 ,  0.01086052, -0.03815957,  0.00904363,
         0.01043104, -0.17194532])


hi_mags = np.linalg.norm((dvdx_hi,dvdy_hi),axis=0)


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-30,30,30),extend='both')
plt.colorbar()
plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle')

plt.xlim(90,-90)
plt.title('hi')


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-40,40,30),extend='both')
plt.colorbar()
plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle')

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)-360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle')


plt.xlim(-90,-270)
plt.ylim(-30,30)
plt.title('hi')


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-30,30,30),extend='both')
plt.colorbar()
plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle',alpha=0.5,linewidth=3)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle',linewidth=0.6)


plt.xlim(90,-90)


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-40,40,30),extend='both')
plt.colorbar()
plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)-360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle',alpha=0.5,linewidth=3)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle',alpha=0.5,linewidth=3)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)-360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle',linewidth=0.6)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle',linewidth=0.6)

plt.xlim(-90,-270)
plt.ylim(-30,30)


fit(line,np.log10(cat['J_GMC<erg.s>']),np.log10(cat['J_HI<erg.s>']))[0]


np.sqrt(np.diag(fit(line,np.log10(cat['J_GMC<erg.s>']),np.log10(cat['J_HI<erg.s>']))[1]))


sml_lbls = ['aq',
'cl',
'ca',
'ch',
'cm',
'xb',
'cr',
'gm',
'hr',
'lc',
'lp',
'md',
'mb',
'mr',
'op',
'oa',
'ob',
'ol',
'pr',
'pl',
'rs',
'ta']


mrkrs = ['$'+sml_lbls[i]+'$' for i in range(22)]


plt.figure(dpi=200)

# plt.scatter(cat['J_GMC<erg.s>'],cat['J_HI<erg.s>'],marker=mrkrs)


plt.plot(cat['J_GMC<erg.s>'],cat['J_GMC<erg.s>']**0.679 * 10**21.777)
plt.plot(cat['J_GMC<erg.s>'],cat['J_GMC<erg.s>'])

for i in range(22):
    plt.scatter(cat['J_GMC<erg.s>'][i],cat['J_HI<erg.s>'][i],
                marker=mrkrs[i],s=100,color='black')
#     plt.annotate(lbls[i],(cat['J_GMC<erg.s>'][i],cat['J_HI<erg.s>'][i]))

plt.loglog()
plt.show()


plt.figure(dpi=200)

# plt.scatter(cat['J_GMC<erg.s>'],cat['J_HI<erg.s>'],marker=mrkrs)


# plt.plot(cat['j_GMC<cm2/s>'],cat['J_GMC<erg.s>']**0.679 * 10**21.777)
plt.plot(cat['j_GMC<cm2/s>'],cat['j_GMC<cm2/s>'],label='$j_{HI} = j_{GMC}$')

for i in range(22):
    plt.scatter(cat['j_GMC<cm2/s>'][i],cat['j_HI<cm2/s>'][i],
                marker=mrkrs[i],s=100,color='black')
#     plt.annotate(lbls[i],(cat['J_GMC<erg.s>'][i],cat['J_HI<erg.s>'][i]))

plt.loglog()
plt.show()


plt.figure(dpi=200)

plt.scatter(cat['J_GMC<erg.s>'],cat['J_HI<erg.s>'],marker='.')


plt.plot(cat['J_GMC<erg.s>'],cat['J_GMC<erg.s>']**0.679 * 10**21.777,
         label='$J_{HI} \propto J_{GMC}^{0.68}$')

plt.plot(cat['J_GMC<erg.s>'],cat['J_GMC<erg.s>'],label='$J_{HI} = J_{GMC}$')


plt.loglog()
plt.legend()
plt.show()

#dot below the orange line is maddalena


#show ^ this plot in terms of ratio of JHI/JGMC as a function of another GMC property
#show theoretical lines as dashed , same color, set the same xlim/ylim so that the equality line is m=1

#color code ^ according to star formation rate(efficiency) / relate in some way 


#also do ratio for little j 


#look at the delta j (JHI-JGMC) as a function of some other property that can indicate timescales


#find a list of all the timescales (rotation, collapse, depletion,
#turbulent crossing time) that we can calculate from the data so far


#send nia rosolowsky, time since formation paper


fit(line,np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['j_HI<cm2/s>']))[0]


np.sqrt(np.diag(fit(line,np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['j_HI<cm2/s>']))[1]))


plt.figure(dpi=200)

plt.scatter(cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>'],marker='.')


plt.plot(cat['j_GMC<cm2/s>'],cat['j_GMC<cm2/s>']**0.344 * 10**17.061,
         label='$j_{HI} \propto j_{GMC}^{0.34}$')

plt.plot(cat['j_GMC<cm2/s>'],cat['j_GMC<cm2/s>'],label='$j_{HI} = j_{GMC}$')


plt.loglog()
plt.legend()
plt.show()


fit(line,np.log10(cat['radius_HI<pc>']),np.log10(cat['j_HI<cm2/s>']))[0]


fit(line,np.log10(cat['radius_GMC<pc>']),np.log10(cat['j_GMC<cm2/s>']))[0]


np.sqrt(np.diag(fit(line,np.log10(cat['radius_HI<pc>']),np.log10(cat['j_HI<cm2/s>']))[1]))


np.sqrt(np.diag(fit(line,np.log10(cat['radius_GMC<pc>']),np.log10(cat['j_GMC<cm2/s>']))[1]))


0.23/1.88


0.39/1.29


0.16/1.5


fit(line,np.log10(np.concatenate([cat['radius_GMC<pc>'],cat['radius_HI<pc>']])),
np.log10(np.concatenate([cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>']])))[0]


np.sqrt(np.diag(fit(line,np.log10(np.concatenate([cat['radius_GMC<pc>'],cat['radius_HI<pc>']])),
np.log10(np.concatenate([cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>']])))[1]))


#the fit line going through both has 1.5-2 less spread than the fit line going through the gmcs and hi 
#individually 


plt.figure(dpi=200)

plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],marker='.',
           color='black',alpha=0.3,s=200,linewidth=0.)

plt.scatter(cat['radius_GMC<pc>'],cat['j_GMC<cm2/s>'],marker='.',color='black')


plt.plot(np.arange(1e0,1e3,1e1),np.arange(1e0,1e3,1e1)**1.52 * 10**22.97,
         label='$j_ \propto R^{1.52 \pm 0.16}$')

plt.fill_between(np.arange(1e0,1e3,1e1),np.arange(1e0,1e3,1e1)**(1.52-0.16) * 10**22.97,
                np.arange(1e0,1e3,1e1)**(1.52+0.16) * 10**22.97,alpha=0.1)


plt.legend(loc='upper left')

plt.loglog()

plt.xlabel('R [pc]')
plt.ylabel('$j\  [\mathrm{cm^2 s^{-1}}]$')
plt.xlim(3e0,3e2)
plt.ylim(2e23,1e27)
plt.show()

#point way below the fit line is maddalena


#find the cores / smaller structure j vs r/ j vs m and add them to the plot, see if they match the trends
#to explore the story of how the j vs r relation evolves over time, need to 
#show the fits separately for each population, and see if the fit for one vs the other fall within each others
#1 sigma error


fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['j_GMC<cm2/s>']))[0]


np.sqrt(np.diag(fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['j_GMC<cm2/s>']))[1]))


0.10/0.75


fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['j_HI<cm2/s>']))[0]


np.sqrt(np.diag(fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['j_HI<cm2/s>']))[1]))


0.13/0.55


fit(line,np.log10(np.concatenate([cat['mass_GMC<Msun>'],cat['mass_HI<Msun>']])),
np.log10(np.concatenate([cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>']])))[0]


np.sqrt(np.diag(fit(line,np.log10(np.concatenate([cat['mass_GMC<Msun>'],cat['mass_HI<Msun>']])),
np.log10(np.concatenate([cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>']])))[1]))


0.069/0.69


#fit if better again with the combined j vs m 


#check if the hi fits would be better without maddalena


#hmm also check if the vel range restriction on maddalena is more dramatic than for cma?
#ie whether the centroid velocity-ing is artificially restricted


plt.figure(dpi=200)

plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>'],marker='.',
            color='black',alpha=0.3,s=200,linewidth=0.)

plt.scatter(cat['mass_GMC<Msun>'],cat['j_GMC<cm2/s>'],marker='.',color='black')


plt.plot(np.arange(1e2,1e7,1e3),np.arange(1e2,1e7,1e3)**0.69 * 10**21.88,
         label='$j_ \propto M^{0.69 \pm 0.07}$')

plt.fill_between(np.arange(1e2,1e7,1e3),np.arange(1e2,1e7,1e3)**(0.69-0.069) * 10**21.88,
                np.arange(1e2,1e7,1e3)**(0.69+0.069) * 10**21.88,alpha=0.1)


plt.legend(loc='upper left')

plt.loglog()
plt.xlim(5e2,1e7)
plt.ylim(2e23,1e27)
plt.show()

#point way below the fit line is maddalena


plt.figure(dpi=200)

plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],marker='.',
           color='black',alpha=0.3,s=200,linewidth=0.)

plt.scatter(cat['radius_GMC<pc>'],cat['j_GMC<cm2/s>'],marker='.',color='black')


plt.plot(np.arange(1e0,1e3,1e1),np.arange(1e0,1e3,1e1)**1.52 * 10**22.97,
         label='$j_ \propto R^{1.52 \pm 0.16}$')

plt.fill_between(np.arange(1e0,1e3,1e1),np.arange(1e0,1e3,1e1)**(1.52-0.16) * 10**22.97,
                np.arange(1e0,1e3,1e1)**(1.52+0.16) * 10**22.97,alpha=0.1)


plt.legend(loc='upper left')

plt.loglog()

plt.xlabel('R [pc]')
plt.ylabel('$j\  [\mathrm{cm^2 s^{-1}}]$')
plt.xlim(3e0,3e2)
plt.ylim(2e23,1e27)
plt.show()


plt.figure(dpi=200)

plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>'],marker='.',
            color='black',alpha=0.3,s=200,linewidth=0.)

plt.scatter(cat['mass_GMC<Msun>'],cat['j_GMC<cm2/s>'],marker='.',color='black')


plt.plot(np.arange(1e2,1e7,1e3),np.arange(1e2,1e7,1e3)**0.69 * 10**21.88,
         label='$j_ \propto M^{0.69 \pm 0.07}$')

plt.fill_between(np.arange(1e2,1e7,1e3),np.arange(1e2,1e7,1e3)**(0.69-0.069) * 10**21.88,
                np.arange(1e2,1e7,1e3)**(0.69+0.069) * 10**21.88,alpha=0.1)


plt.legend(loc='upper left')

plt.loglog()
plt.xlim(5e2,1e7)
plt.ylim(2e23,1e27)

plt.xlabel('$M\ [\mathrm{M_\odot}$]')
plt.ylabel('$j\  [\mathrm{cm^2 s^{-1}}]$')
plt.show()


plt.figure(dpi=200)

plt.hist(np.abs(delphi),bins=np.arange(0,210,15),
         color=seg_spec.colors[19],alpha=0.8,
         edgecolor='white')


# plt.hist(delphi,bins=np.arange(0,210,15),
#          color=seg_spec.colors[17],alpha=0.8,
#          edgecolor='white')

# plt.hist(-delphi,bins=np.arange(0,210,15),
#          color=seg_spec.colors[21],alpha=0.8,
#          edgecolor='white')

plt.xticks(np.arange(0,210,30),
          [str(dg)+'$\degree$' for dg in np.arange(0,210,30)])

plt.vlines([60,120],0,5.3,colors='black',linestyle='--',linewidth=0.6)

plt.annotate('co-rotating',      (10,5.2))
plt.annotate('un-aligned',       (73,5.2))
plt.annotate('counter-rotating',(127,5.2))

plt.ylim(0,5.7)
plt.xlim(-5,185)
# plt.xlim(0,180)

plt.xlabel('$\mathrm{\Delta} \phi$')

plt.show()


#do a histogram showing the diference between the galactic rotation and each cloud (hi, co)
#change around how the arrows look so they can be more distinct from each other


arr_pars = {'headwidth':4,'headlength':4,'headaxislength':4,'width':3e-3}


hi_arr_pars = {'headwidth':4,'headlength':4,'headaxislength':4,'width':6e-3}


plt.figure(dpi=250)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-30,30,30),extend='both')

plt.colorbar(shrink=0.5,ticks=np.arange(-30,40,10),label='$v\mathrm{_{LSR}\ [km\ s^{-1}]}$')


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle',alpha=0.3,
           **hi_arr_pars,)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle',**arr_pars)


plt.xlim(90,-90)
plt.ylim(-30,30)

plt.title('galactic center')
plt.gca().set_aspect(1.5)

plt.xticks(np.arange(90,-90-30,-30))
plt.show()


plt.figure(dpi=200)
plt.contourf(bigh1,extent=[360,-360,-50,50],
             cmap=seg_spec,levels=np.linspace(-50,50,30),extend='both')


plt.colorbar(shrink=0.5,ticks=np.arange(-50,60,20),label='$v\mathrm{_{LSR}\ [km\ s^{-1}]}$')


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)-360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle',alpha=0.3,**hi_arr_pars)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_hi/hi_mags,dvdy_hi/hi_mags,pivot='middle',alpha=0.3,**hi_arr_pars)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)-360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle',**arr_pars)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -dvdx_gmc/gmc_mags,dvdy_gmc/gmc_mags,pivot='middle',**arr_pars)

plt.gca().set_aspect(1.5)

plt.xlim(-90,-270)
plt.ylim(-30,30)

plt.title('galactic anti-center')
plt.xticks(np.arange(-90,-270-30,-30),labels=[-90,-120,-150,180,150,120,90])

plt.show()


for i in range(22):
    dvdy, dvdx = np.gradient(h1s[i].data)
    plt.contourf(np.log10(h0s[i].data),levels=20,cmap='magma')
    plt.quiver(np.indices(h1s[i].data.shape)[1],np.indices(h1s[i].data.shape)[0],
              dvdx,dvdy)
    plt.show()


def calc_jsum(mom0_maps,mom1_maps,cl_typ,cat,masks=0):
    n = len(mom1_maps)
    
    if type(masks) == int:
        masks = [np.full(mom1_maps[i].shape,1) for i in range(22)]
        

    if n == 1:
        dist = [cat['distance<pc>']]
        if cl_typ=='CO':
            lo = [cat['loc']]
            li = [cat['lic']]
            bo = [cat['boc']]
            bi = [cat['bic']]
            mu = 1.36*2

        elif cl_typ=='21cm':
            lo = [cat['loh']]
            li = [cat['lih']]
            bo = [cat['boh']]
            bi = [cat['bih']]
            mu = 1.36
    
    else:
        dist = cat['distance<pc>']
        if cl_typ=='CO':
            lo = cat['loc']
            li = cat['lic']
            bo = cat['boc']
            bi = cat['bic']
            mu = 1.36*2

        elif cl_typ=='21cm':
            lo = cat['loh']
            li = cat['lih']
            bo = cat['boh']
            bi = cat['bih']
            mu = 1.36
        
        
    com_coords = [(np.array(mom1_maps[i].data.shape)*0.5).astype('int') 
                  for i in range(n)]
    
#     print(com_coords)
    
    coords = [np.array((
                 np.meshgrid(
                 np.arange(lo[i],li[i])*-0.125+360,
                 np.arange(bo[i],bi[i])*0.125-50)
                 ))*u.deg.to(u.rad)*dist[i]
              
                    for i in range(n)]
    
    normed_h1s = [mom1_maps[i].data-\
                  mom1_maps[i].data[com_coords[i][0],com_coords[i][1]]
                     for i in range(n)]
    
    for i in range(22): #if the geometric center is somewhere with no velocity data, set vo = median velocity
        if np.all(np.isnan(normed_h1s[i])):
            normed_h1s[i] = mom1_maps[i].data-\
                      np.nanmedian(mom1_maps[i].data)
                         
    
    mass_maps = [((mom0_maps[i].data*(mu*c.m_p)/u.cm**2\
                  ).to(u.solMass/u.pc**2))*\
                 (0.125*u.deg.to(u.rad)*dist[i]*u.pc)**2
                 for i in range(n)]
    
    mskd_mass = [mass_maps[i]*masks[i] for i in range(n)]
    
    
    r_vector = [np.array((
                 coords[i][0] - np.nanmean(coords[i][0]),
                 coords[i][1] - np.nanmean(coords[i][1]),
                 np.zeros(coords[i][0].shape)
                 ))   for i in range(n)]
    
#     print(r_vector[0].reshape(3,-1).transpose())
    
    v_vector = [np.array((
                np.zeros(normed_h1s[i].shape),
                np.zeros(normed_h1s[i].shape),
                normed_h1s[i]
                ))  for i in range(n)]
    
#     print(v_vector[0].reshape(3,-1).transpose())


    jsum_kmpc_s = [np.nansum(np.cross(
                 r_vector[i].reshape(3,-1).transpose()*u.pc,
                 v_vector[i].reshape(3,-1).transpose()*u.km/u.s)*\
                 mskd_mass[i].ravel()[:,np.newaxis],
                 axis=0)/np.nansum(mskd_mass[i])
                   for i in range(n)]
    
    
#     jsum_map = [np.linalg.norm(np.cross(
#                  r_vector[i].reshape(3,-1).transpose()*u.pc,
#                  v_vector[i].reshape(3,-1).transpose()*u.km/u.s)*\
#                  mskd_mass[i].ravel()[:,np.newaxis],axis=1).reshape(mom0_maps[i].data.shape)
                 
#                    for i in range(n)]
    
    
    return jsum_kmpc_s,mskd_mass,r_vector,v_vector


jmaps_gmc = calc_jsum(c0s,c1s,'CO',cat,msks)[0]


[(
plt.contourf(jmaps_gmc[i].value,levels=np.linspace(1,np.nanmax(jmaps_gmc[i].value),20)),
plt.colorbar(),
plt.show(),
    
plt.contourf(np.where(~msks[i],np.nan,c1s[i].data),levels=20,cmap=seg_spec),
plt.colorbar(),
plt.show(),
    
plt.contourf(np.log10(msks[i]*c0s[i].data),levels=20,cmap='magma'),
plt.colorbar(),
plt.show()
)
for i in range(22)]


jmaps_hi = calc_jsum(h0s,h1s,'21cm',cat)[0]


[(
plt.contourf(jmaps_hi[i].value,levels=np.linspace(1,np.nanmax(jmaps_hi[i].value),20)),
plt.colorbar(),
plt.show(),
    
plt.contourf(h1s[i].data,levels=20,cmap=seg_spec),
plt.colorbar(),
plt.show(),
    
plt.contourf(np.log10(h0s[i].data),levels=20,cmap='magma'),
plt.colorbar(),
plt.show()
)
for i in range(22)]


[print(cat['name'][i],'mean:',round(np.nanmedian(h1s[i].data),4),'\n',
      'com:',round(h1s[i].data[int(h1s[i].data.shape[0]/2),int(h1s[i].data.shape[1]/2)],4),'\n\n')
for i in range(22)]


[print(cat['name'][i],np.nanmedian(h1s[i].data)-\
       h1s[i].data[int(h1s[i].data.shape[0]/2),int(h1s[i].data.shape[1]/2)],'\n',
      np.nanstd(h1s[i].data),'\n\n')
for i in range(22)]


[print(cat['name'][i],np.nanmedian(np.where(~msks[i],np.nan,c1s[i].data))-\
       c1s[i].data[int(c1s[i].data.shape[0]/2),int(c1s[i].data.shape[1]/2)],'\n',
      np.nanstd(np.where(~msks[i],np.nan,c1s[i].data)),'\n\n')
for i in range(22)]


[print(cat['name'][i],'mean:',round(np.nanmedian(np.where(~msks[i],np.nan,c1s[i].data)),4),'\n',
      'com:',round(c1s[i].data[int(c1s[i].data.shape[0]/2),int(c1s[i].data.shape[1]/2)],4),'\n\n')
for i in range(22)]


[(plt.figure(dpi=200),
    plt.contourf(h1s[i].data,levels=27,cmap=seg_spec),
  plt.gca().set_aspect('equal'),
  plt.title(str(round(cat['rot_ax_HI<deg>'][i],3))+'    '+\
            str(round(jsum_hi_theta[i],3))+'   '+str(round(np.nanmedian(h1s[i].data),3))),
  plt.colorbar(),
  plt.scatter(int(0.5*h1s[i].data.shape[1]),int(0.5*h1s[i].data.shape[0]),color='black',marker='x'),
  
plt.show())
for i in range(22)]


jsum_hi


calc_jsum(h0s,h1s,'21cm',cat)[0][21]


jsum_gmc = np.array([np.linalg.norm(calc_jsum(c0s,c1s,'CO',cat,msks)[0][i].cgs.value) for i in range(22)])


jsum_hi = np.array([np.linalg.norm(calc_jsum(h0s,h1s,'21cm',cat)[0][i].cgs.value) for i in range(22)])


plt.hist(jsum_gmc/cat['j_GMC<cm2/s>'],alpha=0.5,bins=np.arange(0,2.1,0.1))

plt.hist(jsum_hi/cat['j_HI<cm2/s>'],alpha=0.5,bins=np.arange(0,2.1,0.1))

plt.vlines([0.1,0.5,1],0,8,linestyle='--')

# plt.xticks(rotation='vertical')
plt.show()


plt.hist(jsum_gmc/jsum_hi,alpha=0.5,bins=np.arange(0,2.1,0.1))
plt.hist(cat['j_GMC<cm2/s>']/cat['j_HI<cm2/s>'],alpha=0.5,bins=np.arange(0,2.1,0.1))

plt.show()


def dxdy_2_theta(dx,dy):
    #arctan recovers things -90< theta < 90
    #angles from second qd (90-180) get sent to 4th
    #3rd to first
    th = np.arctan(dy/dx)*u.rad.to(u.deg)
    
    if dx>0.:
        return th
    
    else:
        th = th + 180
        if th > 180:
            return th-360
        else:
            return th


jsum_gmc_theta = [dxdy_2_theta(calc_jsum(c0s,c1s,'CO',cat,msks)[0][i].value[0],
    calc_jsum(c0s,c1s,'CO',cat,msks)[0][i].value[1])
 
 for i in range(22)]


jsum_hi_theta = [dxdy_2_theta(calc_jsum(h0s,h1s,'21cm',cat)[0][i].value[0],
    calc_jsum(h0s,h1s,'21cm',cat)[0][i].value[1])
 
 for i in range(22)]


plt.scatter(cat['name'],cat['rot_ax_GMC<deg>']-jsum_gmc_theta)

# plt.ylim(-50,50)

plt.xticks(rotation='vertical')
plt.hlines([0,10,-10,15,-15,60,-60],'aquila_s','taurus')
plt.show()


plt.scatter(cat['name'],cat['rot_ax_HI<deg>']-jsum_hi_theta)
plt.scatter(cat['name'],cat['rot_ax_HI<deg>']-jsum_hi_theta-360)

plt.ylim(-30,30)

plt.xticks(rotation='vertical')
plt.hlines([0,10,-10,15,-15,],'aquila_s','taurus')
plt.show()


np.array((jsum_gmc_theta,cat['rot_ax_GMC<deg>'])).transpose()


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)


    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    
    a,b,vo = fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[0]
    

    print(cat['name'][i],'\n',
          np.around((
          np.linalg.norm((cxfit[0][0],cyfit[0][0])),
              np.linalg.norm((a,b)),),4),'\n\n')
              
    
for i in range(22):
   
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)


    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
        np.ravel(h1s[i].data),
            )).transpose(),
            columns=['l','b','v'])
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
    
    
    a,b,vo = fit(plane_rvld,np.array((hxpc,hypc)),hdf['v'])[0]
    
    print(round(np.sqrt(hxfit[1][0,0]),7)/\
          round(np.sqrt(np.diag(fit(plane_rvld,np.array((hxpc,hypc)),hdf['v'])[1])[0]),7),'\n\n',
          round(np.sqrt(hyfit[1][0,0]),7)/\
          round(np.sqrt(np.diag(fit(plane_rvld,np.array((hxpc,hypc)),hdf['v'])[1])[1]),7) ,'\n\n\n')
    

#     print(cat['name'][i],'\n',
#           np.around((
#           np.linalg.norm((hxfit[0][0],hyfit[0][0])),
#               np.linalg.norm((a,b)),),4),'\n\n')
              
    
for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)


    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    
    a,b,vo = fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[0]

    
#     print(cat['name'][i],'\n',
#           np.around((
#           np.linalg.norm((cxfit[0][0],cyfit[0][0])),
#               np.linalg.norm((a,b)),),4),'\n\n')
              
    print(round(np.sqrt(cxfit[1][0,0]),7)/\
          round(np.sqrt(np.diag(fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[1])[0]),7),'\n\n',
          round(np.sqrt(cyfit[1][0,0]),7)/\
          round(np.sqrt(np.diag(fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[1])[1]),7) ,'\n\n\n')
    
    
#hmmmm from this tho it looks like the plane has ~similar number of times when the sigma is smaller than the line
#so theyre kind of on hte same level? 

#aybe other way to check which one is doing better is to calculate the first moment maps with the fit, 
#then subtract from the original, then see which has smaller residuals?


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)


    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    
    ca,cb,cvo = fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[0]
    
    
#     print(round(vo-np.nanmedian(cdf['v']),4)>round(vo-np.nanmean(cdf['v']),4))
    
    
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)


    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
        np.ravel(h1s[i].data),
            )).transpose(),
            columns=['l','b','v'])
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
    
    
    ha,hb,hvo = fit(plane_rvld,np.array((hxpc,hypc)),hdf['v'])[0]
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    axs[0,0].scatter(cxpc-np.nanmean(np.unique(cxpc)),
                   cdf['v']-cvo,s=5e2/len(cdf))
    
    axs[0,0].plot(cxpc-np.nanmean(np.unique(cxpc)),ca*(cxpc-np.nanmean(np.unique(cxpc))),color='red')
    axs[0,0].plot(cxpc-np.nanmean(np.unique(cxpc)),cxfit[0][0]*(cxpc-np.nanmean(np.unique(cxpc))),color='orange')
    
    axs[0,1].scatter(cypc-np.nanmean(np.unique(cypc)),
                   cdf['v']-cvo,s=5e2/len(cdf))
    
    
    axs[0,1].plot(cypc-np.nanmean(np.unique(cypc)),cb*(cypc-np.nanmean(np.unique(cypc))),color='red')
    axs[0,1].plot(cypc-np.nanmean(np.unique(cypc)),cyfit[0][0]*(cypc-np.nanmean(np.unique(cypc))),color='orange')

 
    if i == 10:
        for ax in axs[0]:
            ax.set_ylim(-5,)
            
            
    axs[1,0].scatter(hxpc-np.nanmean(np.unique(hxpc)),
                   hdf['v']-hvo,s=5e2/len(hdf))
    
    axs[1,0].plot(hxpc-np.nanmean(np.unique(hxpc)),
                  ha*(hxpc-np.nanmean(np.unique(hxpc))),color='red')
    
    axs[1,0].plot(hxpc-np.nanmean(np.unique(hxpc)),
                  hxfit[0][0]*(hxpc-np.nanmean(np.unique(hxpc))),color='orange',linestyle='--')
    
    axs[1,1].scatter(hypc-np.nanmean(np.unique(hypc)),
                   hdf['v']-hvo,s=5e2/len(hdf))
    
    
    axs[1,1].plot(hypc-np.nanmean(np.unique(hypc)),
                  hb*(hypc-np.nanmean(np.unique(hypc))),color='red')
    
    axs[1,1].plot(hypc-np.nanmean(np.unique(hypc)),
                  hyfit[0][0]*(hypc-np.nanmean(np.unique(hypc))),color='orange',linestyle='--')

            
    if i == 15:
        for ax in axs[1]:
            ax.set_ylim(-10,10)
            
    if i == [13,16]:
        for ax in axs[1]:
            ax.set_ylim(-20,20)
            
            
    fig.suptitle(cat['name'][i])
            
    plt.show()


def fit_gradient(mom1_map, dist, mask=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
    
    
    xdeg = np.indices(mom1_map.shape)[1]*-0.125+360
    ydeg = np.indices(mom1_map.shape)[0]*0.125-50
    
    
    raveldf = pd.DataFrame(np.array((
            np.ravel(xdeg),np.ravel(ydeg),
            np.ravel(np.where(~mask,np.nan,mom1_map)),
            )).transpose(),
            columns=['xdeg','ydeg','v'])
    
    
    noisy = raveldf
    raveldf = noisy[~(np.isnan(noisy['v']))]
    
    
    xpc = raveldf['xdeg']*u.deg.to(u.rad)*dist
    ypc = raveldf['ydeg']*u.deg.to(u.rad)*dist
    
    
    xfit = fit(line,xpc,raveldf['v'])
    yfit = fit(line,ypc,raveldf['v'])
    
    xslope,xintercept = xfit[0]
    yslope,yintercept = yfit[0]
    
    xmsig, xbsig = np.sqrt(np.diag(xfit[1]))
    ymsig, ybsig = np.sqrt(np.diag(xfit[1]))


#     xslope,yslope,xyintercept = fit(plane_rvld,np.array((xpc,ypc)),raveldf['v'])[0]
    
#     xmsig, ymsig, xybsig = np.sqrt(np.diag(fit(plane_rvld,np.array((xpc,ypc)),raveldf['v'])[1]))


    if not full:
        return xslope,yslope
    else:
        return {'x_m_b_msig_bsig':[xslope, xintercept, xmsig, xbsig],
                'y_m_b_msig_bsig':[yslope, yintercept, ymsig, ybsig]}


[(
    fit_gradient(c1s[i].data,cat['distance<pc>'][i],msks[i])
    
) for i in range(22)]


#output when the fit gradient function has the line lines uncommented
line_dvdxys = np.array([(-0.001231859902104127, 0.07005233791454865),
 (0.07228954940696365, 0.12888647665041997),
 (-0.017847817359207263, -0.012049347482188345),
 (-0.09577492755333053, 0.04224276930490577),
 (-0.016132174961491334, 0.029089791878264337),
 (-0.060958131084062846, -0.00017345203980312363),
 (0.07731864860852007, 0.041081144324292725),
 (0.03445931671090977, -0.022797007943680073),
 (-0.24994184980058964, -0.28509074605649),
 (0.10428457400132217, 0.06103475619611021),
 (-0.0368792268325031, 0.0254940334527748),
 (0.0389462030497509, 0.07380126557905216),
 (0.013707919280511396, 0.1589455338263366),
 (-0.03489532207024431, 0.005255247283522072),
 (0.03540596162384518, -0.07957263406307853),
 (-0.03361666723806247, 0.02788993253267396),
 (-0.01717165479730714, 0.012503488907454698),
 (0.05740992477765461, 0.1724369462380287),
 (0.29453494530328106, 0.17551363038481405),
 (0.06868333159870921, 0.04641513663987504),
 (-0.03902938169424294, -0.10649015558926225),
 (-0.06541992054622685, -0.12604041108972813)]).transpose()


line_velgrad = np.linalg.norm(line_dvdxys,axis=0)


#output when the fit gradient function has the plane rvld lines uncommented
plane_dvdxys = np.array([(-0.02432937640095112, 0.08843305133521862),
 (0.05168468622028937, 0.05622251848617604),
 (-0.021783881492251966, -0.029711282020546737),
 (-0.09565977430315388, 0.04130282413357334),
 (-0.022628191748851556, 0.03629933749084624),
 (-0.10045943718977446, -0.03922657348790204),
 (0.0591693164319321, 0.02267741205849124),
 (0.03160599172751688, -0.011523400763089109),
 (-0.23594382788315207, -0.19115686644286986),
 (0.20360705390556427, 0.08288886912945932),
 (-0.027693166254237456, 0.012308811124517094),
 (0.027451265656149124, 0.045594041250475045),
 (-0.014079665085428239, 0.16851467138948806),
 (-0.038111246576822244, -0.021853882643248064),
 (0.020704308333586163, -0.0776485094583872),
 (-0.03716519700176013, 0.052042156752744706),
 (-0.013150620886133946, 0.004021015935587741),
 (0.03151287787061152, 0.16182145471087273),
 (0.1586337738349951, 0.11357573087242723),
 (0.0626296539594079, 0.013925041774291547),
 (-0.03852000823722122, -0.1059643598600987),
 (-0.04666708445030112, -0.09569220146536324)]).transpose()

plane_velgrad = np.linalg.norm(plane_dvdxys,axis=0)


plt.plot(cat['name'],line_dvdxys[0])
plt.plot(cat['name'],plane_dvdxys[0])


plt.xticks(rotation='vertical')
plt.show()


plt.plot(cat['name'],line_dvdxys[1])
plt.plot(cat['name'],plane_dvdxys[1])


plt.xticks(rotation='vertical')
plt.show()


plt.plot(cat['name'],line_velgrad,marker='x')
plt.plot(cat['name'],plane_velgrad,marker='x')


plt.xticks(rotation='vertical')
plt.show()


# plt.plot(cat['name'],line_dvdxys[0])
# plt.plot(cat['name'],plane_dvdxys[0])

plt.plot(cat['name'],np.abs(line_dvdxys[0]-plane_dvdxys[0]/line_dvdxys[0]),marker='x')
plt.plot(cat['name'],np.abs(line_dvdxys[1]-plane_dvdxys[1]/line_dvdxys[1]),marker='x')


plt.xticks(rotation='vertical')
plt.yscale('log')
plt.show()


j_GMC_plane = ((plane_velgrad*u.km/u.s/u.pc)*(cat['radius_GMC<pc>']*u.pc)**2).cgs


cat['j_GMC<cm2/s>']


c1_lvls_dvds = c1_lvls + np.array([-1,1])


c1_lvls_dvds[1][1] = 8
c1_lvls_dvds[2][1] = 2
c1_lvls_dvds[5] = [4,7]
c1_lvls_dvds[7][1] = 17
c1_lvls_dvds[10] =[1,11]
c1_lvls_dvds[14][1] = 11
c1_lvls_dvds[15][1] = 15
c1_lvls_dvds[18][1] = 11
c1_lvls_dvds[19][0] = -6


h1_lvls_dvds = h1_lvls + np.array([-1,1])
h1_lvls_dvds[3][0] = -8
h1_lvls_dvds[5][1] = 10
h1_lvls_dvds[8][0] = 2
h1_lvls_dvds[9][1] = 0
h1_lvls_dvds[10][0] = -30
h1_lvls_dvds[13][1] = 27
h1_lvls_dvds[14] = [-12,15]
h1_lvls_dvds[16][1] = 29
h1_lvls_dvds[17][1] = 14
h1_lvls_dvds[18][1] = 7
h1_lvls_dvds[21][0] = -24


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)


    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    
    #---------------------------------------------------------------------------------------
    

    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)


    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
        np.ravel(h1s[i].data),
            )).transpose(),
            columns=['l','b','v'])
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
    
    
    c_dots = {'color':'black','s':7e2/len(cdf)}
    h_dots = {'color':'black','s':5e2/len(hdf)}
    
    fit_line = {'color':'red','linewidth':0.8}
    sigma_shade = {'color':'red','alpha':0.2,'zorder':-20}
    #---------------------------------------------------------------------------------------
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    
    
    axs[0,0].scatter(cxpc,cdf['v'],**c_dots)
    axs[0,0].plot(cxpc,cxfit[0][0]*cxpc + cxfit[0][1],**fit_line)
    
    axs[0,0].set_title(r'$\partial_x v = $'+str(round(cxfit[0][0],4)),y=0.9)#+\
#                   ' $\pm$ '+str(round(np.sqrt(cxfit[1][0,0]),4)))
    
    
    ##for some reasong not shading the full range of x
#     axs[0,0].fill_between(cxpc,
#                           (cxfit[0][0]+np.sqrt(cxfit[1][0,0]))*cxpc + cxfit[0][1],# + np.sqrt(cxfit[1][1,1]),
#                           (cxfit[0][0]-np.sqrt(cxfit[1][0,0]))*cxpc + cxfit[0][1],# - np.sqrt(cxfit[1][1,1]),
#                   **sigma_shade)
    
    axs[0,1].scatter(cypc,cdf['v'],**c_dots)
    axs[0,1].plot(cypc,cyfit[0][0]*cypc + cyfit[0][1],**fit_line)
    axs[0,1].set_title(r'$\partial_y v = $'+str(round(cyfit[0][0],4)),y=0.9)#+\
#                   ' $\pm$ '+str(round(np.sqrt(cyfit[1][0,0]),4)))
    
#     axs[0,1].fill_between(cypc,
#                           (cyfit[0][0]+np.sqrt(cyfit[1][0,0]))*cypc + hxfit[0][1],# + np.sqrt(cyfit[1][1,1]),
#                           (cyfit[0][0]-np.sqrt(cyfit[1][0,0]))*cypc + hxfit[0][1],# - np.sqrt(cyfit[1][1,1]),
#                   **sigma_shade)
         
    #---------------------------------------------------------------------------------------

    
    axs[1,0].scatter(hxpc,hdf['v'],**h_dots)
    axs[1,0].plot(hxpc,hxfit[0][0]*hxpc + hxfit[0][1],**fit_line)
    axs[1,0].set_title(r'$\partial_x v = $'+str(round(hxfit[0][0],4)),y=0.9)#+\
#                   ' $\pm$ '+str(round(np.sqrt(hxfit[1][0,0]),4)))
    
    
#     axs[1,0].fill_between(hxpc,
#                           (hxfit[0][0]+np.sqrt(hxfit[1][0,0]))*hxpc + hxfit[0][1] + np.sqrt(hxfit[1][1,1]),
#                           (hxfit[0][0]-np.sqrt(hxfit[1][0,0]))*hxpc + hxfit[0][1] - np.sqrt(hxfit[1][1,1]),
#                   **sigma_shade)
    
    axs[1,1].scatter(hypc,hdf['v'],**h_dots)
    axs[1,1].plot(hypc,hyfit[0][0]*hypc + hyfit[0][1],**fit_line)
    axs[1,1].set_title(label=r'$\partial_y v = $'+str(round(hyfit[0][0],4)),y=0.9)#+\
#                   ' $\pm$ '+str(round(np.sqrt(hyfit[1][0,0]),4)))

#     axs[1,1].fill_between(hypc,
#                           (hyfit[0][0]+np.sqrt(hyfit[1][0,0]))*hypc + hyfit[0][1] + np.sqrt(hyfit[1][1,1]),
#                           (hyfit[0][0]-np.sqrt(hyfit[1][0,0]))*hypc + hyfit[0][1] - np.sqrt(hyfit[1][1,1]),
#                   **sigma_shade)
    #---------------------------------------------------------------------------------------
 
    axs[0,0].set_xticks(np.linspace(np.nanmin(cxpc),np.nanmax(cxpc),5),
                       labels = [round((pcx/cat['distance<pc>'][i])*u.rad.to(u.deg),1)
                                 for pcx in np.linspace(np.nanmin(cxpc),np.nanmax(cxpc),5)])
    
    axs[0,0].set_xlim(np.nanmax(cxpc),np.nanmin(cxpc))
    
    axs[0,1].set_xticks(np.linspace(np.nanmin(cypc),np.nanmax(cypc),5),
                       labels = [round((pcy/cat['distance<pc>'][i])*u.rad.to(u.deg),1)
                                 for pcy in np.linspace(np.nanmin(cypc),np.nanmax(cypc),5)])
    

    axs[0,1].set_xlim(np.nanmin(cypc),np.nanmax(cypc))
    
    
    axs[1,0].set_xticks(np.linspace(np.nanmin(hxpc),np.nanmax(hxpc),5),
                       labels = [round((pcx/cat['distance<pc>'][i])*u.rad.to(u.deg),1)
                                 for pcx in np.linspace(np.nanmin(hxpc),np.nanmax(hxpc),5)])
    
    axs[1,0].set_xlim(np.nanmax(hxpc),np.nanmin(hxpc))
    
    
    axs[1,1].set_xticks(np.linspace(np.nanmin(hypc),np.nanmax(hypc),5),
                       labels = [round((pcy/cat['distance<pc>'][i])*u.rad.to(u.deg),1)
                                 for pcy in np.linspace(np.nanmin(hypc),np.nanmax(hypc),5)])
    
    axs[1,1].set_xlim(np.nanmin(hypc),np.nanmax(hypc))
    
    
    for ax in axs[0]:
        ax.set_ylim(c1_lvls_dvds[i])
        ax.set_yticks(np.around(np.linspace(c1_lvls_dvds[i][0],c1_lvls_dvds[i][1],5),1))
        
    for ax in axs[1]:
        ax.set_ylim(h1_lvls_dvds[i])
        ax.set_yticks(np.around(np.linspace(h1_lvls_dvds[i][0],h1_lvls_dvds[i][1],5),1))
    
    
    plt.subplots_adjust(wspace=0.2)
    
            
    fig.suptitle(names[i],y=0.97)
    axs[1,0].set_xlabel('Galactic Longitude [$\degree$]')
    axs[1,1].set_xlabel('Galactic Latitude [$\degree$]')
    
    fig.supylabel('$v \mathrm{_{LSR}\ [km\ s^{-1}}]$')
            
    plt.show()


for i in range(22):
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50) 
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    
    #---------------------------------------------------------------------------------------
    

    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)


    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
        np.ravel(h1s[i].data),
            )).transpose(),
            columns=['l','b','v'])
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
    
    
    print(cat['name'][i],'\n')
    
    print(cxfit[0][1],cyfit[0][1])
    print(hxfit[0][1],hyfit[0][1])
    
    print('\n\n')
    
    
    cxfit = fit(line,cxpc-np.nanmean(np.unique(cxpc)),cdf['v'])
    cyfit = fit(line,cypc-np.nanmean(np.unique(cypc)),cdf['v'])
    
    hxfit = fit(line,hxpc-np.nanmean(np.unique(hxpc)),hdf['v'])
    hyfit = fit(line,hypc-np.nanmean(np.unique(hypc)),hdf['v'])
    
    
    print(cxfit[0][1],cyfit[0][1])
    print(hxfit[0][1],hyfit[0][1])
    
    print('\n\n\n\n\n\n\n')


plt.plot(np.log2(line_velgrad/plane_velgrad),marker='+')
plt.hlines(0,0,22)


line_dvdxys


line_velgrad


plane_velgrad


for i in range(22):
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50) 
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxl = fit(line,cxpc-np.nanmean(np.unique(cxpc)),cdf['v'])
    cyl = fit(line,cypc-np.nanmean(np.unique(cypc)),cdf['v'])
    
    cpln = fit(plane_rvld,[cxpc,cypc],cdf['v'])
    
    
    vl = (xc-np.nanmean(np.unique(xc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cxl[0][0] + \
            (yc-np.nanmean(np.unique(yc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cyl[0][0] + \
            np.nanmean([cxl[0][1],cyl[0][1]])
    
    vp = (xc-np.nanmean(np.unique(xc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cpln[0][0] + \
            (yc-np.nanmean(np.unique(yc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cpln[0][1] + \
            cpln[0][2]
    
    
    dell = np.where(~msks[i],np.nan,c1s[i].data - vl)
    delp = np.where(~msks[i],np.nan,c1s[i].data - vp)
    
    
    plt.hist(np.ravel(dell),bins=100,histtype='step')
    plt.hist(np.ravel(delp),bins=100,histtype='step')
    plt.title(names[i])
    plt.show()
    
    
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(20,5))
    
    kws = {'levels':27,'cmap':seg_spec}
    
    mdln = axs[0].contourf(np.where(~msks[i],np.nan,vl),**kws)
    mdpn = axs[1].contourf(np.where(~msks[i],np.nan,vp),**kws)
    rl = axs[2].contourf(np.where(~msks[i],np.nan,c1s[i].data),**kws)
    
    for im,ax in zip([mdln,mdpn,rl],axs):
        fig.colorbar(im,ax=ax)
        
        
    plt.subplots_adjust(wspace=0.3)
    plt.show()


#just pick one !! then the whole discussion about the merits of each / comparing them can go in an appendix
#so the next tasks are: pick one, do the r vs v plots, calculate the jsum & plot the ratios of that 
#+ how the slope/angle from that compares to the ones from vel grad
#then look at the listof other plots i told nia i'd do
#also, find the time calculation equation from semadeni-vazquez, then calculate those ages 
#also bring back the n(ysos) and/or the sfr (find it for her,cam), then do the sfe / depletion time
#and move on to doing the timescale plots


for i in range(22):
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50) 
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
        np.ravel(np.where(~msks[i],np.nan,c1s[i].data)),
            )).transpose(),
            columns=['l','b','v'])
    
    cdf_noisy = cdf
    
    cdf = cdf[~np.isnan(cdf['v'])]
    
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxl = fit(line,cxpc-np.nanmean(np.unique(cxpc)),cdf['v'])
    cyl = fit(line,cypc-np.nanmean(np.unique(cypc)),cdf['v'])
    
    cpln = fit(plane_rvld,[cxpc,cypc],cdf['v'])
    
    
    vl = (xc-np.nanmean(np.unique(xc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cxl[0][0] + \
            (yc-np.nanmean(np.unique(yc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cyl[0][0] + \
            np.nanmean([cxl[0][1],cyl[0][1]])
    
    vp = (xc-np.nanmean(np.unique(xc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cpln[0][0] + \
            (yc-np.nanmean(np.unique(yc)))*u.deg.to(u.rad)*cat['distance<pc>'][i]*cpln[0][1] + \
            cpln[0][2]
    
    
    dell = np.where(~msks[i],np.nan,c1s[i].data - vl)
    delp = np.where(~msks[i],np.nan,c1s[i].data - vp)
    
    
    print('line>plane: ',np.abs(np.nansum(dell))>np.abs(np.nansum(delp)))
    #plane is slightly better bc there are 13 times when the abs of the sum of 
    #resiuals from line is bigger aka those are less balance/ have larger numbers of deviation 


#okay gonna go with plane, which also helps that it has just one vo


def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


def fit_gradient(mom1_map, dist, mask=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
    
    
    xdeg = np.indices(mom1_map.shape)[1]*-0.125+360
    ydeg = np.indices(mom1_map.shape)[0]*0.125-50
    
    xpc = xdeg*u.deg.to(u.rad)*dist
    ypc = ydeg*u.deg.to(u.rad)*dist
    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map.data)
    
    raveldf = pd.DataFrame(np.array((
            np.ravel(xpc),np.ravel(ypc),
            np.ravel(mskd_mom1),
            )).transpose(),
            columns=['x','y','v'])
    
    
    noisy = raveldf
    raveldf = noisy[~(np.isnan(noisy['v']))] #get rid of the rows with masked data

    #the plane is pinned down to the geometric center of the map
    #the slope of the fit doesnt change with or without xo
    xin = raveldf['x']-np.nanmean(np.unique(raveldf['x'])) 
    yin = raveldf['y']-np.nanmean(np.unique(raveldf['y']))
    
    
    xslope,yslope,xyintercept = fit(plane_rvld,
                                    np.array((xin,yin)),
                                    raveldf['v'])[0]
    
    xmsig, ymsig, xybsig = np.sqrt(np.diag(fit(plane_rvld,
                                               np.array((xin,yin)),
                                               raveldf['v'])[1]))


    if not full:
        return xslope,yslope,xyintercept
    else:
        return {'fit_results':[xslope,yslope,xyintercept],
                '1_sigma_errors':[xmsig, ymsig, xybsig]}


[(cat['name'][i],round(np.nanmedian(c1s[i].data),5),
  round(np.nanmean(c1s[i].data),5),
  round(fit_gradient(c1s[i],cat['distance<pc>'][i],msks[i])[2],5),
 round(c1s[i].data[int(c1s[i].data.shape[0]*0.5),int(c1s[i].data.shape[1]*0.5)],5), 
 )for i in range(22)]
#the mean vs are diff than the vo from the fit
#the closest thing (ishhh) seems to be the v at the geometric center


gmc_velgrad = np.linalg.norm(np.array(
[(fit_gradient(c1s[i],cat['distance<pc>'][i],msks[i]))for i in range(22)]).transpose()[:2],axis=0)


plt.hist(gmc_velgrad/cat['vel_grad_GMC<km/s/pc>'],bins=np.arange(0.3,2.1,0.05))
#14 are now smaller than before, 8 are bigger
#none of them are diff by more than a factor of 2 in either direction
plt.show()


hi_velgrad = np.linalg.norm(np.array(
[(fit_gradient(h1s[i],cat['distance<pc>'][i]))for i in range(22)]).transpose()[:2],axis=0)


plt.hist(hi_velgrad/cat['vel_grad_HI<km/s/pc>'],bins=np.arange(0.3,2.1,0.05))
#as expected nothing changes for the hi
plt.show()


cat['vel_grad_GMC<km/s/pc>'] = gmc_velgrad


dxdy_gmc_plane = np.array(
[(fit_gradient(c1s[i],cat['distance<pc>'][i],msks[i]))for i in range(22)]).transpose()[:2]


theta_gmc_plane = np.array([dxdy_2_theta(dxdy_gmc_plane[0][i],dxdy_gmc_plane[1][i]) for i in range(22)])


rot_gmc_plane = theta_gmc_plane - 90


for i in range(22):
    if rot_gmc_plane[i] > 180:
        rot_gmc_plane[i] -= 360
        
    if rot_gmc_plane[i] < -180:
        rot_gmc_plane[i] += 360


plt.scatter(gmc_velgrad,rot_gmc_plane - cat['rot_ax_GMC<deg>'])
#hmmm makes sense that the ones with the smaller gradients have the largest swings in theta


cat['rot_ax_GMC<deg>'] = rot_gmc_plane


ca,cb = np.array([(
    fit_gradient(c1s[i],cat['distance<pc>'][i],msks[i],full=True)['fit_results'][:2]
) for i in range(22)]).transpose()

csiga,csigb = np.array([(
    fit_gradient(c1s[i],cat['distance<pc>'][i],msks[i],full=True)['1_sigma_errors'][:2]
) for i in range(22)]).transpose()


ha,hb = np.array([(
    fit_gradient(h1s[i],cat['distance<pc>'][i],full=True)['fit_results'][:2]
) for i in range(22)]).transpose()

hsiga,hsigb = np.array([(
    fit_gradient(h1s[i],cat['distance<pc>'][i],full=True)['1_sigma_errors'][:2]
) for i in range(22)]).transpose()


csig_omega = (np.linalg.norm((ca*csiga,cb*csigb),axis=0)/np.linalg.norm((ca,cb),axis=0))


hsig_omega = (np.linalg.norm((ha*hsiga,hb*hsigb),axis=0)/np.linalg.norm((ha,hb),axis=0))


np.array((cat['name'],gmc_velgrad/csig_omega,cat['vel_grad_HI<km/s/pc>']/hsig_omega)).transpose()


plt.plot(cat['name'],gmc_velgrad/csig_omega,marker='.',color='black')
plt.plot(cat['name'],cat['vel_grad_HI<km/s/pc>']/hsig_omega,marker='o',
         color='black',linewidth=3,alpha=0.2,)


plt.xticks(rotation='vertical')
plt.hlines([1,10,100],'aquila_s','taurus',linestyle='--',color='black')
plt.yscale('log')
plt.show()


cat['Omega_over_error_GMC'] = gmc_velgrad/csig_omega
cat['Omega_over_error_HI'] = cat['vel_grad_HI<km/s/pc>']/hsig_omega


cat['j_GMC<cm2/s>'] = ((gmc_velgrad*u.km/u.s/u.pc)*(cat['radius_GMC<pc>']*u.pc)**2).cgs


cat['J_GMC<erg.s>'] = (cat['j_GMC<cm2/s>']*cat['mass_GMC<Msun>']*u.solMass).cgs


cat['dvdx_GMC'] = ca
cat['sigma_dvdx_GMC'] = csiga

cat['dvdy_GMC'] = cb
cat['sigma_dvdy_GMC'] = csigb


cat['dvdx_HI'] = ha
cat['sigma_dvdx_HI'] = hsiga


cat['dvdy_HI'] = hb
cat['sigma_dvdy_HI'] = hsigb


cat['jsum_GMC<cm2/s>'] = jsum_gmc
cat['jsum_HI<cm2/s>'] = jsum_hi

cat['jsum_ang_GMC<deg>'] = jsum_gmc_theta
cat['jsum_ang_HI<deg>'] = jsum_hi_theta


cat


cat.write('sighhplease.fits',overwrite=True)


plt.contourf(display_msks[i]+1)
plt.colorbar()


fns.make_fits_from_list([display_msks[i]+1 for i in range(22)],'h2_0thmom_masks_DISPLAY.fits')
