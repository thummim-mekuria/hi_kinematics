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


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


big_co = fits.open('big_zeroth_co.fits')[0].data


cat = QTable(fits.open('sighhplease.fits')[1].data)


seg_spec_all_color = ListedColormap([
    (0.5633986928104575, 0.0033986928104575, 0.0928104575163399, 1.0),
(0.6633986928104575, 0.0633986928104575, 0.1928104575163399, 1.0),

    (0.8758169934640523, 0.3045751633986928, 0.29411764705882354, 1.0),

    (0.9568627450980393, 0.42745098039215684, 0.2627450980392157, 1.0),

    (0.9803921568627452, 0.5973856209150326, 0.34117647058823525, 1.0),

    (0.9934640522875817, 0.7477124183006535, 0.4352941176470587, 1.0),

    (0.996078431372549, 0.8784313725490196, 0.5450980392156862, 1.0),

    (0.7477124183006538, 0.8980392156862746, 0.6274509803921569, 1.0),

    (0.5803921568627453, 0.8313725490196079, 0.6444444444444445, 1.0),

    (0.26405228758169935, 0.6091503267973857, 0.7098039215686275, 1.0),

    (0.25359477124183005, 0.45882352941176474, 0.7058823529411765, 1.0),

(0.3686274509803922, 0.30980392156862746, 0.6352941176470588, 1.0),
(0.2686274509803922, 0.20980392156862746, 0.5352941176470588, 1.0),
(0.1686274509803922, 0.10980392156862746, 0.4352941176470588, 1.0)])


h2_0s = fits.open('h2_0s.fits')


h2_msks = [(fits.open('h2_0thmom_masks_DISPLAY.fits')[i].data - 1).astype('bool') for i in range(22)]


plt.figure(dpi=200)
plt.contourf(np.log10(big_co),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50])

# plt.xlim(90,-90)
plt.xlim(180,-180)

[(plt.gca().add_patch(Rectangle((cat['loh'][i]*-0.125+360,cat['boh'][i]*0.125-50),
                                -2*(cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),
                                2*(cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),
                                alpha=0.4,color=cm.tab20(i),linewidth=0.,
                               )),
  
) for i in range(22)]


[(
  plt.contourf(np.log10(h2_0s[i].data)*h2_msks[i],extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                                                          cat['boc'][i]*0.125 - 50,cat['bic'][i]*0.125 - 50],
                                                  levels=np.arange(20,23,0.25))
) for i in range(22)]


plt.gca().set_aspect(1.2)


plt.figure(dpi=200)
plt.contourf(np.log10(big_co),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50])

plt.xlim(270,90)

[(plt.gca().add_patch(Rectangle((cat['loh'][i]*-0.125+360,cat['boh'][i]*0.125-50),
                                -2*(cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),
                                2*(cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),
                                fill=False,linewidth=0.75,linestyle='--',
                               )),
  
  
  plt.gca().add_patch(Rectangle((cat['loh'][i]*-0.125+720,cat['boh'][i]*0.125-50),
                                -2*(cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),
                                2*(cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),
                                fill=False,linewidth=0.75,linestyle='--',
                               )),
  
) for i in range(22)]


[(
    
     plt.contourf(np.log10(h2_0s[i].data)*h2_msks[i],extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                                                          cat['boc'][i]*0.125 - 50,cat['bic'][i]*0.125 - 50],
                                                  levels=np.arange(20,23,0.25)),
    
  plt.contourf(np.log10(h2_0s[i].data)*h2_msks[i],extent=[cat['loc'][i]*-0.125+720,cat['lic'][i]*-0.125+720,
                                                          cat['boc'][i]*0.125 - 50,cat['bic'][i]*0.125 - 50],
                                                  levels=np.arange(20,23,0.25)),
    
) for i in range(22)]


plt.gca().set_aspect(1.)


hihi = fits.open('hi_hisized.fits')
coco = fits.open('h2_h2sized.fits')


fns.segmented_spectral


plt.rcParams['font.size'] = 10


names = [
    'Aquila South',
    'California',
    'Camelopardalis',
    'Chamaeleon',
    'Canis Major OB1',
    'Corona Australis',
    'Crossbones',
    'Gemini OB1',
    'Hercules',
    'Lacerta',
    'Lupus',
    'Maddalena',
    'Monoceros OB1',
    'Monoceros R2',
    'Ophiuchus',
    'Orion A',
    'Orion B',
    r'$\lambda$-Orionis',
    'Perseus',
    'Polaris',
    'Rosette',
    'Taurus'  
]


#mentally prepare to change the scaling factor so it is custom for each cloud, but for now we are going with 100x

fig,axs=plt.subplots(4,3,dpi=200,figsize=(10,10))#,sharey=True)

[(
    axs.flat[i].plot(np.arange(155)*1.3-100,np.nanmean(hihi[i].data,axis=(1,2)),linestyle=':',color='black',alpha=0.8,linewidth=3),
    axs.flat[i].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),color='black',linewidth=4),
    
    axs.flat[i].plot(np.arange(155)*1.3-100,100*np.nanmean(coco[i].data,axis=(1,2)),linestyle=':',color=fns.segmented_spectral(20)),
    axs.flat[i].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,100*np.nanmean(coco[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),color=fns.segmented_spectral(20)),
    
    axs.flat[i].set_xlim((cat['vo'][i]-20)*1.3-100,(cat['vi'][i]+15)*1.3-100),
    axs.flat[i].set_title(names[i],y=0.8,x=0.23,fontsize='small')
    
    
) for i in range(12)]


# for ax in axs.flat[-2:]:
#     ax.set_axis_off()
    
fig.supxlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
fig.supylabel('$T\mathrm{_B\ [K]}$')
    
fig.tight_layout()
plt.savefig('postage_spectra_1',bbox_inches='tight')


#mentally prepare to change the scaling factor so it is custom for each cloud, but for now we are going with 100x


fig,axs=plt.subplots(4,3,dpi=200,figsize=(10,10))#,sharey=True)

[(
    axs.flat[i-12].plot(np.arange(155)*1.3-100,np.nanmean(hihi[i].data,axis=(1,2)),linestyle=':',color='black',alpha=0.8,linewidth=3),
    axs.flat[i-12].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),color='black',linewidth=4),
    
    axs.flat[i-12].plot(np.arange(155)*1.3-100,100*np.nanmean(coco[i].data,axis=(1,2)),linestyle=':',color=fns.segmented_spectral(20)),
    axs.flat[i-12].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,100*np.nanmean(coco[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),color=fns.segmented_spectral(20)),
    
    axs.flat[i-12].set_xlim((cat['vo'][i]-20)*1.3-100,(cat['vi'][i]+15)*1.3-100),
    axs.flat[i-12].set_title(names[i],y=0.8,x=0.23,fontsize='small')
    
    
) for i in range(12,22)]


for ax in axs.flat[-2:]:
    ax.set_axis_off()
    
fig.supxlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
fig.supylabel('$T\mathrm{_B\ [K]}$')
    
fig.tight_layout()
plt.savefig('postage_spectra_2',bbox_inches='tight')


hi_0s = fits.open('hi_0s.fits')


hi_msks = [np.ones(hi_0s[i].data.shape).astype('bool') for i in range(22)]

hi_msks[7] = np.log10(hi_0s[7].data) > 21.1
for i in [13,15,16]:
    hi_msks[i] = np.log10(hi_0s[i].data) > 20
hi_msks[14] = np.log10(hi_0s[14].data) > 20.5


h0_lvls = np.array((
    [20.4,21.3],
    [20.8,21.8],
    [20.7,21.2],
    [20.7,21.3],
    [21.1,21.5],
    [20.5,21.1],
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
    
    
plt.contour(np.log10(big_co))


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    axs.flat[i].contourf(np.where(~hi_msks[i],np.nan,np.log10(hi_0s[i].data)),
                         levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                         extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                                                         cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,],
                         cmap='magma',extend='both',norm='asinh'),
    
    axs.flat[i].contour(np.where(~h2_msks[i],np.nan,np.log10(h2_0s[i].data)),linewidths=7e-1,colors='white',
                        levels = 5,alpha=0.7,
                                 extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                                                         cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,]),
    
    
    axs.flat[i].contour(np.log10(big_co),linewidths=4e-1,colors='white',
                        levels = np.linspace(
                            np.nanmin(np.where(~h2_msks[i],np.nan,np.log10(h2_0s[i].data))),
                            np.nanmax(np.where(~h2_msks[i],np.nan,np.log10(h2_0s[i].data))),5),
                        alpha=0.3,
                                 extent=[360,-360,-50,50]),
    
    axs.flat[i].invert_xaxis(),
    axs.flat[i].set_ylim(cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50),
    # axs.flat[i].set_title(names[i],y=0.1,x=0.3,color='white',fontsize='small')
                        
    # axs.flat[i].set_aspect(1)


    axs.flat[i].set_xlim(cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360),
    axs.flat[i].set_ylim(cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50)
)
for i in range(12)]


[(
    axs.flat[i].set_title(names[i],y=0.1,x=0.2,color='white',fontsize='small')
) for i in range(12) if i not in [2,4,5,7,8,10]]


[(
    axs.flat[i].set_title(names[i],y=0.1,x=0.3,color='white',fontsize='small')
) for i in [5]]


[(
    axs.flat[i].set_title(names[i],y=0.8,x=0.3,color='white',fontsize='small')
) for i in [2,4,7,8,10]]


# for ax in axs.flat[-2:]:
#     ax.set_axis_off()
    
fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()
fig.savefig('postage_0thmom_1',bbox_inches='tight')


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    axs.flat[i-12].contourf(np.where(~hi_msks[i],np.nan,np.log10(hi_0s[i].data)),
                         levels=np.linspace(h0_lvls[i][0],h0_lvls[i][1],20),
                         extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                                                         cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,],
                         cmap='magma',extend='both',norm='asinh'),
    
    axs.flat[i-12].contour(np.where(~h2_msks[i],np.nan,np.log10(h2_0s[i].data)),linewidths=7e-1,colors='white',
                        levels = 5,alpha=0.7,
                                 extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                                                         cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,]),
    
    
    axs.flat[i-12].contour(np.log10(big_co),linewidths=4e-1,colors='white',
                        levels = np.linspace(
                            np.nanmin(np.where(~h2_msks[i],np.nan,np.log10(h2_0s[i].data))),
                            np.nanmax(np.where(~h2_msks[i],np.nan,np.log10(h2_0s[i].data))),5),
                        alpha=0.3,
                                 extent=[360,-360,-50,50]),

    axs.flat[i-12].set_xlim(cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360),
    axs.flat[i-12].set_ylim(cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50),
    # axs.flat[i].set_title(names[i],y=0.1,x=0.3,color='white',fontsize='small')
                        
    # axs.flat[i].set_aspect(1)


)
for i in range(12,22)]


[(
    axs.flat[i-12].set_title(names[i],y=0.1,x=0.2,color='white',fontsize='small')
) for i in range(12,22) if i not in [10,12,13,14,19,20]]


[(
    axs.flat[i-12].set_title(names[i],y=0.1,x=0.3,color='white',fontsize='small')
) for i in [13]]


[(
    axs.flat[i-12].set_title(names[i],y=0.8,x=0.3,color='white',fontsize='small')
) for i in [10,12,14,19,20]]


for ax in axs.flat[-2:]:
    ax.set_axis_off()
    
fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()
fig.savefig('postage_0thmom_2',bbox_inches='tight')


hico_scale = [np.around(np.nanmax(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)))/
 np.nanmax(np.nanmean(coco[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2))))
 
 
 for i in range(22)]


hi_1s = fits.open('hi_1s.fits')


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


i = 21

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    xo, xi =  n*33,(n+1)*33
    yo, yi =  n*33,(n+1)*33
    
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            
#     print(round(vo_clr,3))
    
    plt.plot(np.arange(155)*1.3-100,n*80 + np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*80 + np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*80,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*80,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*80 + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),
             color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(-70,n*80+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-75,45)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,806)


# plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
plt.savefig(names[i][:3]+'_stacked_patch_spectra',bbox_inches='tight')


i = 20

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    xo, xi =  min(66,(10-n-1)*7),(10-n)*7
    yo, yi =  n*7,(n+1)*7
    
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
    
    
#     print(this_vo,round(vo_clr,3),np.nanmean(hi_1s[i].data))
    
    up = 100
    
    plt.plot(np.arange(155)*1.3-100,n*up + np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up + np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(-15,n*up+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-20,80)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,1000,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])

plt.ylim(-10,990)

plt.savefig(names[i][:3]+'_stacked_patch_spectra',bbox_inches='tight')


i = 10

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    xo, xi =  min(270,(10-n-1)*27),(10-n)*27#n*27,(n+1)*27
    yo, yi =  n*27, (n+1)*27#min(270,(10-n-1)*27),(10-n)*27
    
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            
#     print(round(vo_clr,3))
    
    plt.plot(np.arange(155)*1.3-100,n*80 + np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*80 + np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*80,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*80,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*80 + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(-80,n*80+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-90,60)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,806)


# plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
plt.savefig(names[i][:3]+'_stacked_patch_spectra',bbox_inches='tight')


#pick stuff where the slope is easy to pull out a stack of patch by patch spectra


np.array(names)[np.abs(np.around(-cat['dvdy_HI']/cat['dvdx_HI'])) == 1]


np.array(names)[np.abs(np.around(-cat['dvdy_HI']/cat['dvdx_HI'])) > 10]


np.array(names)[np.abs(np.around(-cat['dvdx_HI']/cat['dvdy_HI'])) > 10]


np.array(names)[np.abs(np.around(-cat['dvdx_HI']/cat['dvdy_HI'])) ==1]


plt.figure(dpi=200)
plt.plot(cat['vel_grad_HI<km/s/pc>'],marker='+')
plt.scatter(np.arange(22)[[3,10,18,19,21,2,4,5,8,15,0,11]],cat['vel_grad_HI<km/s/pc>'][[3,10,18,19,21,2,4,5,8,15,0,11]],color='red')
plt.grid()
plt.yscale('log')


plt.figure(dpi=200)
plt.plot(cat['vel_grad_GMC<km/s/pc>'],marker='+')
plt.scatter(np.arange(22)[[3,10,18,19,21,2,4,5,8,15,0,11]],cat['vel_grad_GMC<km/s/pc>'][[3,10,18,19,21,2,4,5,8,15,0,11]],color='red')
plt.scatter(np.arange(22)[[10,18,21,5,8,15,11,0]],cat['vel_grad_GMC<km/s/pc>'][[10,18,21,5,8,15,11,0]],color='orange')

plt.grid()
# plt.yscale('log')


#gonna go with these for now
#made lupus and taurus so far

#gonna do the others now


lupus, aquila, taurus

maddalena, orion a, 

corona australis, hercules, perseus


(-cat['dvdx_HI']/cat['dvdy_HI'])[[0,11,15,5,8,18]]


(cat['dvdy_HI']/-cat['dvdx_HI'])[[0,11,15,5,8,18]]


hihi[11].data.shape[1]


h1_lvls[11]


plt.hist(np.ravel(hi_1s[11].data))


i = 11


plt.figure(dpi=200,figsize=(7,10))
for n in range(1,6):
   
    s = round(hihi[i].data.shape[1]/10)
    
    xo, xi =  12 + n*s, 12 + (n+1)*s
    yo, yi =  2*n*s -6, (2*n+1)*s -6  #min(270,(10-n-1)*27),(10-n)*27
    
    # print(xo,xi,yo,yi)
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            

    up = 230
    plt.plot(np.arange(155)*1.3-100,n*up - up/2 + np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up -up/2 + np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up -up/2,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up-up/2,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up -up/2 + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(-15,n*up-up/2+10),fontsize='small')
    
    
for n in range(0,5):
   
    s = round(hihi[i].data.shape[1]/10)
    
    xo, xi =  12 + n*s, 12 + (n+1)*s   
    yo, yi =  2*n*s , (2*n+1)*s   #min(270,(10-n-1)*27),(10-n)*27
    
    # print(xo,xi,yo,yi)
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            
    
    plt.plot(np.arange(155)*1.3-100,n*up + np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up + np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360,2))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50,2))+'$\degree$',(-15,n*up +10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-20,90)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,1300,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,1130)


plt.savefig(names[i][:3]+'_stacked_patch_spectra',bbox_inches='tight')


hihi[15].data.shape


198/2


99-10


i = 15


plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
   
    s = round(hihi[i].data.shape[1]/10)
    
    xo, xi =  99,110
    yo, yi =  n*s, (n+1)*s  #min(270,(10-n-1)*27),(10-n)*27
    
    # print(xo,xi,yo,yi)
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            

    up = 80
    plt.plot(np.arange(155)*1.3-100,n*up+ np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up+ np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50,1))+'$\degree$',(40,n*up+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-20,50)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,800,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,800)


plt.savefig(names[i][:3]+names[i][-1]+'_stacked_patch_spectra',bbox_inches='tight')


hihi[8].data.shape


13*4


i = 8


plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
   
    s = round(hihi[i].data.shape[1]/10)
    
    xo, xi =  60,75
    yo, yi =  n*s, (n+1)*s  #min(270,(10-n-1)*27),(10-n)*27
    
    # print(xo,xi,yo,yi)
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            

    up = 80
    plt.plot(np.arange(155)*1.3-100,n*up+ np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up+ np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50,1))+'$\degree$',(40,n*up+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-20,50)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,800,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,800)


plt.savefig(names[i][:3]+names[i][-1]+'_stacked_patch_spectra',bbox_inches='tight')


hihi[5].data.shape


148/2


59+15


i = 5


plt.figure(dpi=200,figsize=(7,10))
for n in range(5):
   
    s = round(hihi[i].data.shape[1]/10)
    
    xo, xi =  74,89
    yo, yi =  n*s, (n+1)*s  #min(270,(10-n-1)*27),(10-n)*27
    
    # print(xo,xi,yo,yi)
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            

    up = 80
    plt.plot(np.arange(155)*1.3-100,n*up+ np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up+ np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up + 70*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(40,n*up+10),fontsize='small')
    
    
for n in range(5,10):
   
    s = round(hihi[i].data.shape[1]/10)
    
    xo, xi =  59,74
    yo, yi =  n*s, (n+1)*s  #min(270,(10-n-1)*27),(10-n)*27
    
    # print(xo,xi,yo,yi)
    
    this_vo = fns.get_mom1(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            

    up = 80
    plt.plot(np.arange(155)*1.3-100,n*up+ np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*up+ np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*up,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*up,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*up + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(40,n*up+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-20,50)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,800,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,800)


# plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
plt.savefig(names[i][:3]+names[i][-1]+'_stacked_patch_spectra',bbox_inches='tight')


#when i come back -- go back and change the co scaling for these to be more consistent
#then put in overleaf and start describing them

#and look at other to do list things !


#this is the most correct version of the code for making hi pbp spectra
for i in range(22):
    n = hihi[i].data.shape[1]
    
    num_ptch = 10
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    coords = np.array((inds[:-1],inds[1:])).transpose()
    slcs = np.array([(i,j) for j in coords for i in coords])
    
    
    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=300)
    
    map_axs = np.flip(axs,axis=0).flat

        
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
    # bkg_map.set_axis_off()
    bkg_map.contourf(hi_1s[i].data,cmap=fns.segmented_spectral,extend='both',
                    extent=[loh,lih,boh,bih],
                    zorder=-20,alpha=0.4,
                    levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
                    antialiased=True)
   
    
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
        bo,bi = slcs[j][1]
        lo,li = slcs[j][0]
        
        sp = np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],axis=(1,2))
        # sp2 = np.nanmean(coshsz[i][:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.6,
                        linewidth=1.5,color='black')
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
        map_axs[j].vlines(np.nanmedian(hi_1s[i].data),
                          0,1.5,color='black',
                          linestyle='--',linewidth=1.,alpha=0.4)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(hihi[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
                           s=7,color='black',alpha=0.5,linewidth=0)
        
        
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
    plt.savefig('patch_spec_'+names[i][:3]+names[i][-1])

    plt.show()


res = random.shuffle([
    
    (0.80568627450980394, 0.62921568627450981, 0.9939215686274509), 
    (0.60568627450980394, 0.42921568627450981, 0.7939215686274509), 
    (0.41568627450980394, 0.23921568627450981, 0.6039215686274509),
    
    
 (0.9058823529411765, 0.5411764705882353, 0.8047058823529411),
    (0.8058823529411765, 0.4411764705882353, 0.7047058823529411),
    
    
    (0.954313725490196, 0.4019607843137255, 0.3470588235294118),
    (0.854313725490196, 0.2519607843137255, 0.2070588235294118),
    
 (0.9882352941176471, 0.5529411764705883, 0.3843137254901961),
    (0.8882352941176471, 0.4529411764705883, 0.2843137254901961),
    
    (1.0, 0.8509803921568627, 0.1843137254901961),
    (0.85, 0.7009803921568627, 0.0343137254901961),
    
 
 (0.6509803921568628, 0.8470588235294118, 0.32941176470588235),
    (0.5509803921568628, 0.7470588235294118, 0.22941176470588235),
        
 
 (0.4, 0.7607843137254902, 0.6470588235294118),
    (0.3, 0.6607843137254902, 0.5470588235294118),
    
    
     (0.6519607843137255, 0.8441176470588235, 0.9474509803921568),
     (0.4019607843137255, 0.5941176470588235, 0.7274509803921568),
    

    (0.6029411764705883, 0.6774509803921569, 0.846078431372549),
    (0.4529411764705883, 0.5274509803921569, 0.696078431372549),
    
    
    (0.7019607843137254, 0.7019607843137254, 0.7019607843137254),
    (0.4519607843137254, 0.4519607843137254, 0.4519607843137254),
    (0.30,0.30,0.30),
])


res


distinct = ListedColormap((
    
    (0.80568627450980394, 0.62921568627450981, 0.9939215686274509), 
    (0.60568627450980394, 0.42921568627450981, 0.7939215686274509), 
    (0.41568627450980394, 0.23921568627450981, 0.6039215686274509),
    
    
 (0.9058823529411765, 0.5411764705882353, 0.8047058823529411),
    (0.8058823529411765, 0.4411764705882353, 0.7047058823529411),
    
    
    (0.954313725490196, 0.4019607843137255, 0.3470588235294118),
    (0.854313725490196, 0.2519607843137255, 0.2070588235294118),
    
 (0.9882352941176471, 0.5529411764705883, 0.3843137254901961),
    (0.8882352941176471, 0.4529411764705883, 0.2843137254901961),
    
    (1.0, 0.8509803921568627, 0.1843137254901961),
    (0.85, 0.7009803921568627, 0.0343137254901961),
    
 
 (0.6509803921568628, 0.8470588235294118, 0.32941176470588235),
    (0.5509803921568628, 0.7470588235294118, 0.22941176470588235),
        
 
 (0.4, 0.7607843137254902, 0.6470588235294118),
    (0.3, 0.6607843137254902, 0.5470588235294118),
    
    
     (0.6519607843137255, 0.8441176470588235, 0.9474509803921568),
     (0.4019607843137255, 0.5941176470588235, 0.7274509803921568),
    

    (0.6029411764705883, 0.6774509803921569, 0.846078431372549),
    (0.4529411764705883, 0.5274509803921569, 0.696078431372549),
    
    
    (0.7019607843137254, 0.7019607843137254, 0.7019607843137254),
    (0.4519607843137254, 0.4519607843137254, 0.4519607843137254),
    (0.30,0.30,0.30),


))


distinct


big_h1 = fits.open('bigh1.fits')[0].data


big_h0 = fits.open('bigh0.fits')[0].data


plt.rcParams['font.size'] = 30


plt.figure(dpi=200,figsize=(27,10))

these = [0,3,5,8,10,14]

plt.contourf(big_h1,extent=[360,-360,-50,50],
             levels=np.linspace(-30,30,27),cmap=fns.segmented_spectral,
            extend='both')

# plt.contourf(np.log10(big_h0),extent=[360,-360,-50,50],cmap='magma',levels=np.arange(20,22,0.25))

plt.colorbar(shrink=0.8,ticks=np.arange(-30,40,15),label='$\overline{v_\mathrm{LSR}}\ [\mathrm{km\ s^{-1}}]$')

# plt.contourf(np.log10(big_co),extent=[360,-360,-50,50],
#              levels=np.arange(20.5,23,0.25),cmap='binary_r',alpha=0.5,antialiased=True)


plt.ylim(-40,40)
# plt.gca().set_aspect(1.)

plt.xlim(90,-90)

plt.xticks(np.arange(90,-100,-30),labels = [str(tck)+'$\degree$' for tck in np.arange(90,-100,-30)])
plt.yticks(np.arange(-40,60,20), labels = [str(tck)+'$\degree$' for tck in np.arange(-40,60,20)])

plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')


[(plt.gca().add_patch(
    Rectangle((cat['loh'][i]*-0.125+360,cat['boh'][i]*0.125-50),
             -2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             color='black',linewidth=0.8,fill=False,linestyle='--')
)
)  for i in these
]


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)[these],
           np.nanmean(cat['lat_GMC<deg>'],axis=1)[these],
          (-cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'])[these],
          (cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'])[these],
          pivot='middle',width=10e-3,headwidth=2,headlength=2,headaxislength=1.5,
          scale=15,alpha=0.5)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)[these],
           np.nanmean(cat['lat_GMC<deg>'],axis=1)[these],
          (-cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'])[these],
           (cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'])[these],
          pivot='middle',width=4e-3,headwidth=2,headlength=2,headaxislength=1.5,
          scale=15)


# [(plt.contourf(np.log10(np.where(~h2_msks[i],np.nan,h2_0s[i].data)),levels=np.arange(20.5,22.,0.25),
#                extent = [cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
#                          cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,],cmap='binary_r')
#  )  for i in these
# ]


plt.savefig('catalog_with_velgrad_galcent',bbox_inches='tight')

plt.show()


plt.rcParams['font.size'] = 30


plt.figure(dpi=200,figsize=(27,10))

these = range(22)#[1,9,19,2,18,7,17,11,12,15,4]

plt.contourf(big_h1,extent=[360,-360,-50,50],
             levels=np.linspace(-60,60,25),cmap=fns.segmented_spectral,
            extend='both')

# plt.contourf(np.log10(big_h0),extent=[360,-360,-50,50],cmap='magma',levels=np.arange(20,22,0.25))

plt.colorbar(shrink=0.8,ticks=np.arange(-60,80,20),label='$\overline{v_\mathrm{LSR}}\ [\mathrm{km\ s^{-1}}]$')

# plt.contourf(np.log10(big_co),extent=[360,-360,-50,50],
#              levels=np.arange(20.5,23,0.25),cmap='binary_r',alpha=0.5,antialiased=True)


plt.ylim(-40,40)
plt.gca().set_aspect(1.)

plt.xlim(270,90)

plt.xticks(np.arange(270,60,-30),labels = [str(tck)+'$\degree$' for tck in np.arange(270,60,-30)])
plt.yticks(np.arange(-40,60,20), labels = [str(tck)+'$\degree$' for tck in np.arange(-40,60,20)])

plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')


[(plt.gca().add_patch(
    Rectangle((cat['loh'][i]*-0.125+360,cat['boh'][i]*0.125-50),
             -2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             color='black',linewidth=0.8,fill=False,linestyle='--')
)
)  for i in these
]


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)[these],
           np.nanmean(cat['lat_GMC<deg>'],axis=1)[these],
          (-cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'])[these],
          (cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'])[these],
          pivot='middle',width=10e-3,headwidth=2,headlength=2,headaxislength=1.5,
          scale=25,alpha=0.5)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)[these],
           np.nanmean(cat['lat_GMC<deg>'],axis=1)[these],
          (-cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'])[these],
           (cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'])[these],
          pivot='middle',width=4e-3,headwidth=2,headlength=2,headaxislength=1.5,
          scale=25)


# [(plt.contourf(np.log10(np.where(~h2_msks[i],np.nan,h2_0s[i].data)),levels=20,
#                extent = [cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
#                          cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,],cmap='binary_r')
#  )  for i in these
# ]


[(plt.gca().add_patch(
    Rectangle((cat['loh'][i]*-0.125+720,cat['boh'][i]*0.125-50),
             -2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             color='black',linewidth=0.8,fill=False,linestyle='--')
)
)  for i in these
]


# [(plt.contourf(np.log10(np.where(~h2_msks[i],np.nan,h2_0s[i].data)),levels=20,
#                extent = [cat['loc'][i]*-0.125+720,cat['lic'][i]*-0.125+720,
#                          cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,],cmap='binary_r')
#  )  for i in these
# ]


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)[these]+360,
           np.nanmean(cat['lat_GMC<deg>'],axis=1)[these],
          (-cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'])[these],
          (cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'])[these],
          pivot='middle',width=10e-3,headwidth=2,headlength=2,headaxislength=1.5,
          scale=25,alpha=0.5)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)[these]+360,
           np.nanmean(cat['lat_GMC<deg>'],axis=1)[these],
          (-cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'])[these],
           (cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'])[these],
          pivot='middle',width=4e-3,headwidth=2,headlength=2,headaxislength=1.5,
          scale=25)


plt.savefig('catalog_with_velgrad_galanti',bbox_inches='tight')

plt.show()


plt.rcParams['font.size'] = 30


plt.figure(dpi=200)#figsize=(10,10))


plt.hist(np.abs(cat['dvdy_GMC']/cat['dvdx_GMC']),bins=10**np.arange(-2.,2.,0.2),
        color=fns.segmented_spectral(19),edgecolor='white',)


plt.hist(np.abs(cat['dvdy_HI']/cat['dvdx_HI']),
         bins=10**np.arange(-2.,2.,0.2),
         fill=False,edgecolor='black',linewidth=2)

plt.vlines(1,0,4.5,color=fns.segmented_spectral(4),linestyle='--',linewidth=0.7)

plt.xscale('log')
plt.xlim(1e-1,1e2)
plt.ylim(0,4.5)

plt.ylabel('N (clouds per bin)')
plt.xlabel('$\Omega_y\ /\ \Omega_x$')
plt.show()


delphi = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']
for i in range(22):
    if delphi[i] < -180:
        delphi[i] += 360
        
    if delphi[i] > 180:
        delphi[i] -= 360


help(plt.colorbar)


plt.figure(dpi=200)

plt.contourf(np.log10(big_co),extent=[360,-360,-50,50],levels=np.arange(20.5,23,0.25),cmap='binary_r')

plt.ylim(-40,40)
plt.gca().set_aspect(1.)

plt.xlim(200,90)

# plt.xticks(np.arange(90,-100,-30),labels = [str(tck)+'$\degree$' for tck in np.arange(90,-100,-30)])
plt.yticks(np.arange(-30,45,15), labels = [str(tck)+'$\degree$' for tck in np.arange(-30,45,15)])

plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')


[(plt.gca().add_patch(
    Rectangle((cat['loh'][i]*-0.125+360,cat['boh'][i]*0.125-50),
             -2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             color=distinct.reversed()(i),alpha=0.5,linewidth=0.)
),
  
  
  plt.gca().add_patch(
    Rectangle((cat['loh'][i]*-0.125+720,cat['boh'][i]*0.125-50),
             -2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             2*cat['radius_HI<pc>'][i]/cat['distance<pc>'][i]*u.rad.to(u.deg),
             color=distinct.reversed()(i),alpha=0.5,linewidth=0.)
),
  
  
 )  for i in range(22)
]


[(plt.contourf(np.log10(np.where(~h2_msks[i],np.nan,h2_0s[i].data)),levels=np.arange(20.5,23,0.25),
               extent = [cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                         cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,],cmap='binary_r'),
  
  plt.contourf(np.log10(np.where(~h2_msks[i],np.nan,h2_0s[i].data)),levels=np.arange(20.5,23,0.25),
               extent = [cat['loc'][i]*-0.125+720,cat['lic'][i]*-0.125+720,
                         cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50,],cmap='binary_r')
,
  
  
 )  for i in range(22)
]


plt.show()


phillips99_data = pd.read_csv('J_A+AS_134_241_table1.dat.gz.txt',skiprows=9,delimiter='|')

phillips99_data.columns = ['Name', 'Note', 'Ah Am  RAs  Ed Em  DEs', 'l_Dist', 'Dist',
       'u_Dist', 'RV', 'DV   DV2', 'l_Size', 'Size', 'n_Size', 'Size2',
       'Trans', 'L_OMEGA', 'OMEGA__OMEGA2', 'l_GAMMA', 'GAMMA',
       'l_Mass', 'Mass__Mass2', 'MassEm', 'THETA_THETA2', 'OMEGAo', 'u_OMEGAo', 'Type',
       'r_']


pOmegas = (np.array([(float(n.split()[0])) 
                     for i,n in enumerate(phillips99_data['OMEGA__OMEGA2'][1:-1])])\
           *1e-14/u.s).to(u.km/u.s/u.pc).value


pRs = phillips99_data['Size'][1:-1].astype('float')


pMs = np.zeros(len(phillips99_data['Mass__Mass2'][1:-1]))

for i,n in enumerate(phillips99_data['Mass__Mass2'][1:-1]):
    if len(n.split())==0.:
        pMs[i] = 0
    elif len(n.split())==1.:
        pMs[i] = float(n.split()[0])
        
    elif len(n.split())==3.:
        pMs[i] = np.nanmean([float(n.split()[0]),float(n.split()[2])])
              
    elif n.split()[0][-1]=='>':
        pMs[i] = np.nanmean([float(n.split()[0][:6]),float(n.split()[1])])
        
    elif n.split()[1][:2]=='->':
        pMs[i] = np.nanmean([float(n.split()[0]),float(n.split()[1][2:])])
        
    else:
        pMs[i] = np.nanmean([float(n.split()[0]),float(n.split()[1])])


phillips = pd.DataFrame(np.array([pMs,pRs,pOmegas]).transpose(),columns=['M<Msun>','R<pc>','Omega<km/s/pc>'])

phillips = phillips[~np.isinf(np.log10(pMs))]


phillips['j<cm2/s>'] = (np.array(phillips['Omega<km/s/pc>']*phillips['R<pc>']**2)*u.km/u.s*u.pc).cgs.value


plt.figure(dpi=200,figsize=(10,10))
plt.scatter(phillips['R<pc>'],phillips['j<cm2/s>'],marker='+',s=1e2)
plt.scatter(cat['radius_GMC<pc>'],cat['j_GMC<cm2/s>'],marker='.',s=1e2)
# plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],marker='x',s=1e2)

plt.loglog()


# plt.scatter(cat['radius_GMC<pc>'],cat['sigma_v_GMC<km/s>'])
plt.scatter(cat['radius_GMC<pc>']*cat['mass_GMC<Msun>']/(cat['radius_GMC<pc>'])**2,cat['sigma_v_GMC<km/s>'])

plt.plot(cat['radius_GMC<pc>'],10**(-0.2)*cat['radius_GMC<pc>']**0.5)

plt.loglog()


plt.figure(figsize=(20,20),dpi=200)

plt.scatter(cat['distance<pc>']*np.sin(np.nanmean(cat['lon_GMC<deg>'],axis=1)*u.deg),
cat['distance<pc>']*np.cos(np.nanmean(cat['lon_GMC<deg>'],axis=1)*u.deg),marker='.')
plt.scatter(0,0,color='red',marker='+')

[(
    plt.annotate(cat['name'][i],(cat['distance<pc>'][i]*np.sin(np.nanmean(cat['lon_GMC<deg>'][i])*u.deg),
cat['distance<pc>'][i]*np.cos(np.nanmean(cat['lon_GMC<deg>'][i])*u.deg)),
                fontsize='x-small'),
    
) for i in range(22)]

plt.plot(1e3*np.sin(np.arange(-180,190,10)*u.deg),1e3*np.cos(np.arange(-180,190,10)*u.deg),
        linewidth=0.7,color='pink')


[(plt.plot(np.arange(0,1e4,1e2)*np.sin(ang),np.arange(0,1e4,1e2)*np.cos(ang),
           linewidth=0.5,linestyle='--',color='black',alpha=0.4)


) for ang in np.arange(-180,190,15)*u.deg]


# [(
# plt.gca().add_patch(
#     Circle((cat['distance<pc>'][i]*np.sin(np.nanmean(cat['lon_GMC<deg>'][i])*u.deg),
#             cat['distance<pc>'][i]*np.cos(np.nanmean(cat['lon_GMC<deg>'][i])*u.deg)),
#           cat['radius_HI<pc>'][i],fill=False, )
# )

# ) for i in range(22)]


plt.xlim(-1500,1000)
plt.ylim(-2000,500)

plt.gca().set_aspect(1)


# plt.savefig('birdseye',bbox_inches='tight')
plt.show()


#using orion a as 'this is what a rotating cloud looks like'


h2_1s = fits.open('h2_1s.fits')


for i in [15]:#range(22):
    ny,nx = coco[i].data.shape[1:]
    
    num_ptch_y = 6
    num_ptch_x = 16
    
    inds_x = np.linspace(0,nx,num_ptch_x+1).astype('int')
    coords_x = np.array((inds_x[:-1],inds_x[1:])).transpose()
    
    inds_y = np.linspace(0,ny,num_ptch_y+1).astype('int')
    coords_y = np.array((inds_y[:-1],inds_y[1:])).transpose()
    
    
    slcs = np.array([(i,j) for j in coords_y for i in coords_x])
    
    
    fig,axs = plt.subplots(num_ptch_y,num_ptch_x,
                           sharex=True,sharey=True,
                          dpi=300)
    
    map_axs = np.flip(axs,axis=0).flat 

        
#     loh = cat['loh'][i]*-0.125+360
#     lih = cat['lih'][i]*-0.125+360
    
#     boh = cat['boh'][i]*0.125-50
#     bih = cat['bih'][i]*0.125-50
    
    
    loc = cat['loc'][i]*-0.125+360
    lic = cat['lic'][i]*-0.125+360

    boc = cat['boc'][i]*0.125-50
    bic = cat['bic'][i]*0.125-50
    
    
#     background_ax = plt.axes([0, 0, 1, 1]) # create a dummy subplot for the background
# background_ax.set_zorder(-1) # set the background subplot behind the others
# background_ax.imshow(image, aspect='auto') # show the backgroud image
    
    
    bkg_map = fig.add_subplot(111)#add_axes([0, 0, 1, 1])
    bkg_map.set_zorder(-10)
    # bkg_map.set_axis_off()
    bkg_map.contourf(np.where(~h2_msks[i],np.nan,h2_1s[i].data),
                     cmap=fns.segmented_spectral,extend='both',
                    extent=[loc,lic,boc,bic],
                    zorder=-20,alpha=0.4,
                    levels=26, #np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
                    antialiased=True)
   
    
#     this = np.log10(c0s[i].data)*display_msks[i]
#     smol,big = np.nanmin(this[this>18.]),np.nanmax(this[this>18])


#     bkg_map.contour(np.log10(c0s[i].data)*display_msks[i],
#                 levels=[smol],
#                 extent=[loc,lic,boc,bic],
#                 colors='black',linewidths=0.6,alpha=0.3)
    
#     bkg_map.invert_xaxis()
    bkg_map.set_xlim(loc,lic)
    bkg_map.set_ylim(boc,bic)
#     bkg_map.set_aspect('equal')
        
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][1]
        lo,li = slcs[j][0]
        
        sp = np.nanmean(coco[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],axis=(1,2))
        # sp2 = np.nanmean(coshsz[i][:,bo:bi,lo:li],axis=(1,2))
        
#         map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
#                         (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.6,
#                         linewidth=1.5,color='black')
        
        if np.nansum(np.abs(sp))>3e-1:
            map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,sp/np.nanmax(sp),
                           color='black',linewidth=0.8)
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
        map_axs[j].vlines(np.nanmedian(np.where(~h2_msks[15],np.nan,h2_1s[i].data)),
                          0,1.5,color='black',
                          linestyle='--',linewidth=1.,alpha=0.4)
        
#         map_axs[j].set_title(j,fontsize='xx-small')
        
#         map_axs[j].scatter(np.nanmean(fns.get_mom1(hihi[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],
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
#     plt.savefig('patch_spec_'+names[i][:3]+names[i][-1])

#     bkg_map.set_ylim(-22,)

    plt.show()


plt.figure(dpi=200)

plt.contourf(np.nanmean(coco[15].data,axis=1),levels=np.linspace(0.05,2,10),
             extent = [cat['lon_GMC<deg>'][15][0],cat['lon_GMC<deg>'][15][1],
                       -100,100],cmap='YlGnBu_r')


plt.hlines(np.nanmedian(np.where(~h2_msks[15],np.nan,h2_1s[15].data)),
          cat['lon_GMC<deg>'][15][0],cat['lon_GMC<deg>'][15][1],
          linestyle='--',color='red')

plt.plot((np.arange(cat['loc'][15],cat['lic'][15])*-0.125+360),
         (np.arange(cat['loc'][15],cat['lic'][15])*-0.125+360)*\
         cat['vel_grad_GMC<km/s/pc>'][i]*cat['distance<pc>'][15]*u.rad.to(u.deg),
                
          linestyle='--',color='red')


plt.ylim(-5,20)
plt.xlim(-140,-155)
plt.gca().set_aspect(0.3)
# plt.colorbar()


plt.rcParams['font.size'] = 15


patch_nums = [np.around(np.array(coco[i].data.shape[1:])/8) for i in range(22)]


cos_hisz = fits.open('h2_hisized.fits')


[hihi[i].data.shape for i in range(22)]


for i in range(22):
    ny,nx = cos_hisz[i].data.shape[1:] 
    #ig now its ended up being swuare, but there is a version where i do this for the co sized co maps
    
    num_ptch_y = 15#6
    num_ptch_x = 15#16
    
    inds_x = np.linspace(0,nx,num_ptch_x+1).astype('int')
    coords_x = np.array((inds_x[:-1],inds_x[1:])).transpose()
    
    inds_y = np.linspace(0,ny,num_ptch_y+1).astype('int')
    coords_y = np.array((inds_y[:-1],inds_y[1:])).transpose()
    
    
    slcs = np.array([(i,j) for j in coords_y for i in coords_x])
    
    
    fig,axs = plt.subplots(num_ptch_y,num_ptch_x,
                           sharex=True,sharey=True,
                          dpi=300,figsize=(20,20))
    
    map_axs = np.flip(axs,axis=0).flat 

        
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
    # bkg_map.set_axis_off()
    
#     bkg_map.contourf(cos_hisz[i].data,
#                      cmap=fns.segmented_spectral,extend='both',
#                     extent=[loh,lih,boh,bih],
#                     zorder=-20,alpha=0.,
#                     levels=26, #np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
#                     antialiased=True)
    
    bkg_map.contourf(np.where(~h2_msks[i],np.nan,h2_1s[i].data),
                     cmap=fns.segmented_spectral,extend='both',
                    extent=[loc,lic,boc,bic],
                    zorder=-20,alpha=0.4,
                    levels=26, #np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
                    antialiased=True)
   
    
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
        bo,bi = slcs[j][1]
        lo,li = slcs[j][0]
        
        sp = np.nanmean(cos_hisz[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],axis=(1,2))
        # sp2 = np.nanmean(coshsz[i][:,bo:bi,lo:li],axis=(1,2))
        
#         map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
#                         (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.6,
#                         linewidth=1.5,color='black')
        
        if np.nansum(np.abs(sp))>3e-1:
            map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,sp/np.nanmax(sp),
                           color='black',linewidth=0.8)        
            
            
            map_axs[j].vlines(np.nanmedian(np.where(~h2_msks[i],np.nan,h2_1s[i].data)),
                              0,1.5,color='black',
                              linestyle='--',linewidth=1.,alpha=0.4)
        
#         map_axs[j].set_title(j,fontsize='xx-small')
        
            map_axs[j].scatter(fns.get_mom1(sp,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),0,
                           s=7,color='black',alpha=0.5,linewidth=0)
        
        
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
#     fig.supxlabel('Galactic Longitude')
#     fig.supxlabel('Galactic Latitude')
    plt.savefig('h2_hisz_patch_spec_'+names[i][:3]+names[i][-1])

#     bkg_map.set_ylim(-22,)

    plt.show()


hi_2s = fits.open('hi_2s.fits')


h2_hisz_2s = [np.sqrt(fns.get_momN(cos_hisz[i].data[cat['vo'][i]:cat['vi'][i],:,:],
             np.arange(cat['vo'][i],cat['vi'][i])*1.3,2)) for i in range(22)]


h2_hisz_0s = [np.nansum(cos_hisz[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=0)*1.3*1.97e20
              for i in range(22)]


[(
plt.contourf(h2_hisz_2s[i].data*(np.log10(h2_hisz_0s[i])>20.),
             levels=np.arange(0.5,10)),
    plt.colorbar(),
plt.show()
) for i in range(22)]


[(
plt.hist(np.ravel(h2_hisz_2s[i].data),bins=100,alpha=0.5),
    plt.hist(np.ravel(h2_hisz_2s[i].data*(np.log10(h2_hisz_0s[i])>20.)),bins=100,alpha=0.5),

    plt.yscale('log'),
plt.show()
 )

for i in range(22)]


[(
plt.hist(np.ravel(hi_2s[i].data),bins=100,alpha=0.5),
    plt.hist(np.ravel(hi_2s[i].data*(np.log10(h2_hisz_0s[i])>20.)),
             alpha=0.5,bins=100),
    plt.yscale('log'),
    plt.title(cat['name'][i]),
plt.show()
 )

for i in range(22)]


from astrodendro import Dendrogram


mad_hi_dend = Dendrogram.compute(hihi[11].data,min_npix=50)


len(mad_hi_dend)


len(mad_hi_dend)


len(mad_hi_dend)


for lf in mad_hi_dend.leaves:
    plt.contourf(np.nansum(lf.get_mask(),axis=1))
    plt.contour(np.nansum(cos_hisz[11].data,axis=1))
    plt.title(lf.idx)
    plt.ylim(50,120)
    plt.show()


p = mad_hi_dend.plotter()


mad_hi_dend[352].get_mask()


plt.contourf(np.nansum(mad_hi_dend[352].get_mask(),axis=2),levels=np.arange(0.5,10))


plt.contourf(np.nansum(mad_hi_dend[350].get_mask(),axis=2),levels=np.arange(0.5,10))
plt.contour(np.nansum(cos_hisz[11].data,axis=2))

plt.ylim(50,)


plt.contourf(np.nansum(mad_hi_dend[350].get_mask(),axis=1),levels=np.arange(0.5,10))
plt.contour(np.nansum(cos_hisz[11].data,axis=1))

plt.ylim(50,)


mad_hi_dend[350].parent


plt.contourf(np.nansum(mad_hi_dend[350].get_mask(),axis=0),levels=np.arange(0.5,10))
plt.contour(np.nansum(cos_hisz[11].data,axis=0))


plt.contourf(np.nansum(hihi[11].data,axis=2))
plt.contour(np.nansum(cos_hisz[11].data,axis=2))


plt.contourf(np.nansum(hihi[11].data,axis=1))
plt.contour(np.nansum(cos_hisz[11].data,axis=1))


new_mad_hi = mad_hi_dend[350].get_mask() * hihi[11].data


plt.rcParams['font.size'] = 15


plt.figure(dpi=200)
i = 11
plt.contourf(np.log10(np.nansum(new_mad_hi,axis=0)*1.3*1.82e18),
             cmap='magma',levels=np.linspace(20.6,21.5,20),
            extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                   cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50])

# plt.contourf(np.nansum(mad_hi_dend[342].get_mask(),axis=0),levels=np.arange(0.5,10),
#              extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
#                    cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50])


plt.colorbar(ticks = np.arange(20.6,21.5,0.1),label='log ( $\mathrm{N_{HI}\ /\ [cm^2\ s^{-1}]}$)')

plt.contour(np.log10(h2_hisz_0s[i].data),levels=np.arange(20.5,21.68,0.1),
            extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                   cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],colors='black',
           linewidths=0.2)

plt.contour(np.log10(h2_0s[i].data)*h2_msks[i],levels=np.arange(20.5,21.68,0.1),
            extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                   cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50],colors='black',linewidths=0.9)


plt.gca().invert_xaxis()

plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')

plt.savefig('new_madd_0th')


np.nansum(((np.nansum(new_mad_hi,axis=0)*1.3*1.82e18/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2))*\
(0.125*u.degree.to(u.rad) * cat['distance<pc>'][11])**2*u.pc**2


cat['mass_HI<Msun>'][11]


306107.52/cat['mass_HI<Msun>'][11]


#this is the most correct version of the code for making hi pbp spectra
for i in [11]:#range(22):
    n = hihi[i].data.shape[1]
    
    num_ptch = 10
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    coords = np.array((inds[:-1],inds[1:])).transpose()
    slcs = np.array([(i,j) for j in coords for i in coords])
    
    
    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=300)
    
    map_axs = np.flip(axs,axis=0).flat

        
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
    # bkg_map.set_axis_off()
    bkg_map.contourf(hi_1s[i].data,cmap=fns.segmented_spectral,extend='both',
                    extent=[loh,lih,boh,bih],
                    zorder=-20,alpha=0.4,
                    levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26),
                    antialiased=True)
   
    
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
        bo,bi = slcs[j][1]
        lo,li = slcs[j][0]
        
        sp = np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(new_mad_hi[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],axis=(1,2))
        # sp2 = np.nanmean(coshsz[i][:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        (sp)/np.nanmax(sp),alpha=0.6,
                        linewidth=1.,color='black')
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        (sp2)/np.nanmax(sp2),
                        linewidth=1.,color='purple')
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
        map_axs[j].vlines(np.nanmedian(hi_1s[i].data),
                          0,1.5,color='black',
                          linestyle='--',linewidth=1.,alpha=0.4)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(hihi[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
                           s=7,color='black',alpha=0.5,linewidth=0)
        
        map_axs[j].scatter(fns.get_mom1(sp2,
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),0,
                           s=7,color='purple',alpha=0.5,linewidth=0)
        
        
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
#     plt.savefig('patch_spec_'+names[i][:3]+names[i][-1])

    plt.show()


from scipy.optimize import curve_fit as fit


def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


help(fit)


def fit_gradient(mom1_map, dist, mask=0, weights=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
        
    if type(weights)==int:#all the data is usable and/or a mask wasnt supplied
        weights = np.ones(mom1_map.shape).astype('bool')
    
    
    xdeg = np.indices(mom1_map.shape)[1]*-0.125+360
    ydeg = np.indices(mom1_map.shape)[0]*0.125-50
    
    xpc = xdeg*u.deg.to(u.rad)*dist
    ypc = ydeg*u.deg.to(u.rad)*dist
    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map)
    
    
    raveldf = pd.DataFrame(np.array((
            np.ravel(xpc),np.ravel(ypc),
            np.ravel(mskd_mom1),
            np.ravel(weights),
            )).transpose(),
            columns=['x','y','v','unc'])
    
    
    noisy = raveldf
    raveldf = noisy[~(np.isnan(noisy['v']))] #get rid of the rows with masked data
    
    
    #the plane is pinned down to the geometric center of the map
    #the slope of the fit doesnt change with or without xo
    xin = raveldf['x']-np.nanmean(np.unique(raveldf['x'])) 
    yin = raveldf['y']-np.nanmean(np.unique(raveldf['y']))
    
    
    xslope,yslope,xyintercept = fit(plane_rvld,
                                    np.array((xin,yin)),
                                    raveldf['v'],
                                    sigma=1/raveldf['unc'])[0]
    
    xmsig, ymsig, xybsig = np.sqrt(np.diag(fit(plane_rvld,
                                               np.array((xin,yin)),
                                               raveldf['v'],
                                               sigma=1/raveldf['unc'])[1]))

   
    if not full:
        return xslope,yslope,xyintercept
    else:
        return {'fit_results':[xslope,yslope,xyintercept],
                '1_sigma_errors':[xmsig, ymsig, xybsig]}


plt.scatter(np.indices(new_mad_mom1.shape)[1],np.ravel(new_mad_mom1),alpha=0.5,s=7)
# plt.scatter(np.indices(new_mad_mom1.shape)[1],np.ravel(hi_1s[11].data),alpha=0.5,s=7)


plt.scatter(np.indices(new_mad_mom1.shape)[0],np.ravel(new_mad_mom1),alpha=0.5,s=7)
plt.scatter(np.indices(new_mad_mom1.shape)[0],np.ravel(hi_1s[11].data),alpha=0.5,s=7)


new_mad_mom1 = fns.get_mom1(new_mad_hi[cat['vo'][11]:cat['vi'][11],:,:],
             np.arange(cat['vo'][11],cat['vi'][11])*1.3 - 100)


new_mad_mom2 = np.sqrt(fns.get_momN(new_mad_hi[cat['vo'][11]:cat['vi'][11],:,:],
             np.arange(cat['vo'][11],cat['vi'][11])*1.3 - 100,2))


hi_msks


[(
    plt.scatter(np.ravel(np.log10(hi_0s[i].data)),
           np.ravel(np.where(~hi_msks[i],np.nan,hi_1s[i].data)-\
                   np.nanmedian(np.where(~hi_msks[i],np.nan,hi_1s[i].data))),s=1e-1),
    
    plt.scatter(np.ravel(np.log10(h2_0s[i].data)),
           np.ravel(np.where(~h2_msks[i],np.nan,h2_1s[i].data)-\
                   np.nanmedian(np.where(~h2_msks[i],np.nan,h2_1s[i].data))),s=1),
plt.title(cat['name'][i]),
plt.show()

) for i in range(22)]


plt.contourf(np.log10(hi_0s[21].data))


plt.contourf(hi_1s[21].data)


fit_gradient(new_mad_mom1,
                cat['distance<pc>'][11],weights = new_mad_mom2,full=True)


fit_gradient(new_mad_mom1,
                cat['distance<pc>'][11],full=True)


fit_gradient(hi_1s[11].data,cat['distance<pc>'][11],full=True)


np.linalg.norm([0.011390083877113611, 0.004810981342382735,])


np.linalg.norm([0.017260574292064956, 0.008693657353093687])


np.linalg.norm([-0.000436518399124175, 0.0008233399193423585,])


cat['vel_grad_HI<km/s/pc>'][11]


0.0193263318574657/0.0009318998689643575


np.nansum(np.nansum(new_mad_hi[cat['vo'][11]:cat['vi'][11],:,:],axis=0)!=0)*\
(0.125*u.deg.to(u.rad)*cat['distance<pc>'][11])**2


np.sqrt(31080.17656893926/np.pi)


cat['radius_HI<pc>'][11]


99.46420193331454/134.00406313631336


20*0.23*0.74**2


20*0.74**2


plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'])
plt.scatter(cat['radius_HI<pc>'][11],cat['j_HI<cm2/s>'][11],color='purple')

plt.scatter(99.464,10.952*cat['j_HI<cm2/s>'][11],color='red')
plt.loglog()


np.sqrt(fns.get_momN(np.nanmean(new_mad_hi[cat['vo'][11]:cat['vi'][11],:,:],axis=(1,2)),
            np.arange(cat['vo'][11],cat['vi'][11])*1.3-100,2))


cat['sigma_v_HI<km/s>'][11]


3.2/5.65


plt.plot(np.nanmean(new_mad_hi,axis=(1,2)))
plt.plot(100*np.nanmean(coco[11].data,axis=(1,2)))
plt.plot(np.nanmean(hihi[11].data,axis=(1,2)))
plt.vlines([cat['vo'][11],cat['vi'][11]],0,60)


np.nansum(np.nanmean(new_mad_hi,axis=(1,2)))


np.nansum(np.nanmean(new_mad_hi[cat['vo'][11]:cat['vi'][11],:,:],axis=(1,2)))


#23% of the mass, ^peaks together with maddalena, same width
#vv the first moment map resembles maddalena so much better


plt.contourf(np.where(~h2_msks[11],np.nan,h2_1s[11].data),levels=27,cmap=fns.segmented_spectral)
plt.colorbar()


plt.contourf(fns.get_mom1(new_mad_hi[cat['vo'][11]:cat['vi'][11],:,:],
             np.arange(cat['vo'][11],cat['vi'][11])*1.3 - 100),levels=27,
            cmap=fns.segmented_spectral)
plt.colorbar()


np.nanmean(hi_1s[7].data)


[(plt.contourf(np.nansum(hihi[i].data,axis=2).transpose(),
               extent=[-100,100,cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50]),
  plt.contour(np.nansum(coco[i].data,axis=2).transpose(),
             extent=[-100,100,cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50]),
  plt.vlines(cat['vel<km/s>'][i],cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,color='white'),
  plt.title(cat['name'][i]),
plt.show()) for i in range(22)]


#how far the clouds are off the galactic plane (in the galactic z direction)
np.array((cat['name'],np.cos(u.deg*np.nanmean(cat['lat_GMC<deg>'],axis=1))*cat['distance<pc>'])).transpose()


cat['radius_HI<pc>'][15]


#newnewnew velocity bounds, from looking at ^
noovs = np.array([[ -10. ,  10. ],
       [-10.,  15],
       [-10,  10.],
       [-10,  10],
       [  5.3,  30. ],
       [ -9. ,  19.6],
       [ -6.4,  28.7],
       [-15.5,  32.6],
       [ -9. ,  19.6],
       [-20,  15],
       [-10. ,  13],
       [ 17,  40.4],
       [-12.9,  15. ],
       [ -5.1,  25],
       [-10,  15],
       [ -9. ,  25],
       [-10.3,  25],
       [-12.9,  28.7],
       [-10,  17. ],
       [-10,  13.1],
       [ -5.1,  27.4],
       [-10,  18.3]])


noovos, noovis = ((noovs+100)/1.3).astype('int').transpose()


for i in range(22):

    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
    
    fig.suptitle(names[i])

    axs[0].contourf(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=2)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=2))),15),
                 extent = [cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50, -100,100],
                   cmap='YlGnBu_r'),
    
    axs[0].contour(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=2)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=2))),7),
                 extent = [cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50, -100,100],
                   colors='black',linewidths=0.6),

    axs[0].hlines(cat['vel<km/s>'][i],cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,
              color='red',linestyle='--',linewidth=2)

    axs[0].hlines(noovs[i],cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,
              color='black',linestyle=':',linewidth=2)
    
    axs[0].set_xlabel('latitude')
    axs[0].set_ylabel('velocity [km/s]')
    axs[0].set_xlim(cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50)

    
    axs[1].contourf(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=1)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=1))),15),
                 extent = [cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,-100,100],
                   cmap='YlGnBu_r'       ),
    
    axs[1].contour(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=1)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=2))),7),
                 extent = [cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360, -100,100],
                   colors='black',linewidths=0.6),

    axs[1].hlines(cat['vel<km/s>'][i],cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
              color='red',linestyle='--',linewidth=2)

    axs[1].hlines(noovs[i],cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
              color='black',linestyle=':',linewidth=2)
    axs[1].set_yticks([])
    axs[1].set_xlim(cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360)
    axs[1].set_xlabel('longitude')
    
    
    for ax in axs:
        ax.set_ylim(cat['vel<km/s>'][i]+np.array([-15,15]))
        
        
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.1)

    plt.show()


# plt.colorbar()


hlv,hbv = fits.open('hi4pi_pv_maps.fits')


for i in range(22):

    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
    
    fig.suptitle(names[i])

    
    axs[0].contourf(np.log10(hbv.data),levels=np.linspace(20,22.5,15),extent=[-50,50,-100,100],cmap='binary_r')
    
    
    axs[0].contourf(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=2)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=2))),15),
                 extent = [cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50, -100,100],
                   cmap='YlGnBu_r',alpha=0.8,antialiased=True),
    
    axs[0].contour(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=2)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=2))),7),
                 extent = [cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50, -100,100],
                   colors='black',linewidths=0.6),

    axs[0].hlines(cat['vel<km/s>'][i],cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,
              color='red',linestyle='--',linewidth=2)

    axs[0].hlines(noovs[i],cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50,
              color='black',linestyle=':',linewidth=2)
    
    axs[0].set_xlabel('latitude')
    axs[0].set_ylabel('velocity [km/s]')
    axs[0].set_xlim(cat['boc'][i]*0.125-50-20,cat['bic'][i]*0.125-50+20)
    
    
    axs[1].contourf(np.log10(hlv.data),levels=np.linspace(20,22.5,15),extent=[360,-360,-100,100],cmap='binary_r')

    axs[1].contourf(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=1)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.82e18*np.nansum(hihi[i].data,axis=1))),15),
                 extent = [cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,-100,100],
                   cmap='YlGnBu_r',alpha=0.8,antialiased=True),
    
    axs[1].contour(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=1)),
                 levels=np.linspace(20,np.nanmax(np.log10(0.125*1.97e20*np.nansum(coco[i].data,axis=2))),7),
                 extent = [cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360, -100,100],
                   colors='black',linewidths=0.6),

    axs[1].hlines(cat['vel<km/s>'][i],cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
              color='red',linestyle='--',linewidth=2)

    axs[1].hlines(noovs[i],cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
              color='black',linestyle=':',linewidth=2)
    
    
    axs[1].set_yticks([])
    axs[1].set_xlim(cat['loh'][i]*-0.125+360+30,cat['lih'][i]*-0.125+360-30)
    axs[1].set_xlabel('longitude')
#     plt.xlim(180,-180)
    plt.show()
    

plt.scatter(cat['vel_grad_HI<km/s/pc>'],nooomgs/cat['vel_grad_HI<km/s/pc>'])


np.linalg.norm(fit_gradient(new_mad_mom1,cat['distance<pc>'][11])[:-1])


0.0193263318574657/cat['vel_grad_HI<km/s/pc>'][11]


np.array((np.arange(22),nooomgs/cat['vel_grad_HI<km/s/pc>'])).transpose()


nooms = np.array([fns.mass_n_area(hihi[i].data[noovos[i]:noovis[i],:,:]*u.K,
                 cat['distance<pc>'][i]*u.pc,'21cm')[0].to(u.solMass).value
        for i in range(22)])


noomoo1s = [fns.get_mom1(hihi[i].data[noovos[i]:noovis[i],:,:],
                       np.arange(noovos[i],noovis[i])*1.3-100) for i in range(22)]


np.nanmax(noomoo1s[14])


hi_msks[14] = (noomoo1s[14] > -50)&(noomoo1s[14] < 20)


[(plt.figure(figsize=(3,3,)),
  plt.title(i),
    plt.contourf(np.where(~hi_msks[i],np.nan,hi_1s[i].data),
                 levels=27,cmap=fns.segmented_spectral),
plt.colorbar(),
plt.show()
 ) for i in range(22)]


[(plt.figure(figsize=(3,3,)),
  plt.title(i),
    plt.contourf(np.where(~hi_msks[i],np.nan,noomoo1s[i]),
                 levels=27,cmap=fns.segmented_spectral),
  
plt.colorbar(),
plt.show()
 ) for i in range(22)]


[(
    plt.hist(np.ravel(np.where(~hi_msks[i],np.nan,noomoo1s[i])),bins=100,alpha=0.5,density=True),

plt.hist(np.ravel(np.where(~hi_msks[i],np.nan,hi_1s[i].data)),bins=100,alpha=0.5,density=True),

plt.hist(np.ravel(np.where(~h2_msks[i],np.nan,h2_1s[i].data)),
         bins=100,alpha=0.5,density=True,histtype='step'),
  plt.title(cat['name'][i]),

plt.show()
 ) for i in range(22)]


noomoo2s = [np.sqrt(fns.get_momN(hihi[i].data[noovos[i]:noovis[i],:,:],
                       np.arange(noovos[i],noovis[i])*1.3-100,2)) for i in range(22)]


[(plt.figure(figsize=(3,3,)),
  plt.title(i),
    plt.contourf(noomoo2s[i],
                 levels=27,cmap='magma'),
  
plt.colorbar(),
plt.show()
 ) for i in range(22)]


np.linalg.norm([1.,1.])


def fit_gradient(mom1_map, dist, mask=0, weights=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
        
    if type(weights)==int:#all the data is usable and/or a mask wasnt supplied
        weights = np.ones(mom1_map.shape).astype('bool')
    
    
    xdeg = np.indices(mom1_map.shape)[1]*-0.125+360
    ydeg = np.indices(mom1_map.shape)[0]*0.125-50
    
    xpc = xdeg*u.deg.to(u.rad)*dist
    ypc = ydeg*u.deg.to(u.rad)*dist
    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map)
    
    
    raveldf = pd.DataFrame(np.array((
            np.ravel(xpc),np.ravel(ypc),
            np.ravel(mskd_mom1),
            np.ravel(weights),
            )).transpose(),
            columns=['x','y','v','unc'])
    
    
    noisy = raveldf
    raveldf = noisy[~(np.isnan(noisy['v']))] #get rid of the rows with masked data
    
    
    #the plane is pinned down to the geometric center of the map
    #the slope of the fit doesnt change with or without xo
    xin = raveldf['x']-np.nanmean(np.unique(raveldf['x'])) 
    yin = raveldf['y']-np.nanmean(np.unique(raveldf['y']))
    
    
    xslope,yslope,xyintercept = fit(plane_rvld,
                                    np.array((xin,yin)),
                                    raveldf['v'],
                                    sigma=1/raveldf['unc'])[0]
    
    xmsig, ymsig, xybsig = np.sqrt(np.diag(fit(plane_rvld,
                                               np.array((xin,yin)),
                                               raveldf['v'],
                                               sigma=1/raveldf['unc'])[1]))

#     plt.scatter(yin,raveldf['v'])
# #     plt.plot(xin,xin*xslope + xyintercept,color='red')
#     plt.show()
    
#     print(fit(plane_rvld,
#                                     np.array((xin,yin)),
#                                     raveldf['v'],
#                                     sigma=1/raveldf['unc']))
    
#     return(xin,yin)
    if not full:
        return xslope,yslope,xyintercept
    else:
        return {'fit_results':[xslope,yslope,xyintercept],
                '1_sigma_errors':[xmsig, ymsig, xybsig]}


plt.plot(np.arange(155)*1.3-100,np.nanmean(hihi[14].data,axis=(1,2)))
plt.vlines(cat['vel<km/s>'][14],0,50,color='red',linestyle='--')
plt.vlines(noovs[14],0,50,color='black',linestyle=':')

plt.plot(np.arange(155)*1.3-100,40*np.nanmean(coco[14].data,axis=(1,2)))

plt.xlim(-60,60)

plt.xlabel('velocity [km/s]')
plt.ylabel('brightness temperature [K]')


nooomgs = [np.linalg.norm(fit_gradient(noomoo1s[i],cat['distance<pc>'][i],mask=hi_msks[i],
            weights=noomoo2s[i])[:-1]) for i in range(22)]

#idk why but the fit for ophiuchus is going off the rails when there is the weights

nooomgs[14] = np.linalg.norm([(fit_gradient(noomoo1s[i],
                              cat['distance<pc>'][i],
                              mask=hi_msks[i])[:-1]) 
                for i in [14]])


np.nanmedian(np.around(nooms/cat['mass_HI<Msun>']*100,3))


np.nanmedian(np.around(nooomgs/cat['vel_grad_HI<km/s/pc>']*100,3))


noojs = ((nooomgs*u.km/u.s/u.pc)*(cat['radius_HI<pc>']*u.pc)**2).cgs


plt.scatter(cat['mass_HI<Msun>'], cat['j_HI<cm2/s>'])
plt.scatter(nooms,noojs)

plt.plot([cat['mass_HI<Msun>'],nooms],[cat['j_HI<cm2/s>'],noojs,],linewidth=0.2)

plt.loglog()


plt.scatter(cat['radius_HI<pc>'],noojs)
plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'])

plt.vlines(cat['radius_HI<pc>'],noojs,cat['j_HI<cm2/s>'])

plt.loglog()


plt.hist(noojs.value/cat['j_GMC<cm2/s>'],bins=10**(np.arange(-1,2,0.1)))
plt.xscale('log')

plt.show()


plt.scatter(cat['radius_HI<pc>']/cat['radius_GMC<pc>'],nooomgs/cat['vel_grad_GMC<km/s/pc>'])
plt.loglog()


[(
    plt.plot(np.nanmean(hihi[i].data[noovos[i]:noovis[i],:,:],axis=(1,2))),
    plt.show()


) for i in range(22)]


def line(x,m,b):
    return m*x+b


curve_fit(line,np.log10(nooms),np.log10(noojs.value))


curve_fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['j_HI<cm2/s>']))


np.sqrt(0.01872284)


np.sqrt( 0.01740454)


curve_fit(line,np.log10(cat['radius_HI<pc>']),np.log10(cat['j_HI<cm2/s>']))


np.sqrt(0.15062333)


curve_fit(line,np.log10(cat['radius_HI<pc>']),np.log10(noojs.value))


np.sqrt(0.1495107)


nooms[11],fns.mass_n_area(new_mad_hi[noovos[11]:noovis[11],:,:]*u.K,
                 cat['distance<pc>'][11]*u.pc,'21cm')[0].to(u.solMass).value


306107.53884846694/1494762.417753858


plt.contourf(np.nansum(new_mad_hi[noovos[11]:noovis[11],:,:],axis=1))


plt.contourf(np.nansum(hihi[11].data[noovos[11]:noovis[11],:,:],axis=1))


plt.contourf(np.nansum(new_mad_hi[noovos[11]:noovis[11],:,:],axis=2))


plt.contourf(np.nansum(hihi[11].data[noovos[11]:noovis[11],:,:],axis=2))


np.around(nooms/cat['mass_HI<Msun>']*100,3)


plt.plot(np.around(nooms/cat['mass_HI<Msun>']*100,3),marker='+')


[(plt.contourf(np.nansum(hihi[i].data,axis=2).transpose(),
               extent=[-100,100,cat['boh'][i],cat['bih'][i],]),
  
  plt.contour(np.nansum(coco[i].data,axis=2).transpose(),
             extent=[-100,100,cat['boc'][i],cat['bic'][i],]),
  
  plt.vlines(cat['vel<km/s>'][i],cat['boh'][i],cat['bih'][i],color='white'),
  plt.vlines(noovs[i],cat['boh'][i],cat['bih'][i],color='white',linestyle='--'),
  
  plt.xlim(-50,50),
  plt.title(cat['name'][i]),
plt.show()) for i in range(22)]


[(
plt.plot(np.arange(155)*1.3-100,np.nanmean(hihi[i].data,axis=(1,2))),
    plt.vlines(noovs[i],0,100),
    plt.vlines(cat['vel<km/s>'][i],0,100),
    plt.show()

) for i in range(22)]


[(plt.contourf(np.nansum(hihi[i].data,axis=1),
               extent=[cat['loh'][i],cat['lih'][i],-100,100,]),
  plt.contour(np.nansum(coco[i].data,axis=1),
             extent=[cat['loc'][i],cat['lic'][i],-100,100,]),
  
  plt.hlines(cat['vel<km/s>'][i],cat['loh'][i],cat['lih'][i],color='white'),
  plt.ylim(-50,50),
  plt.title(cat['name'][i]),
plt.show()) for i in range(22)]


# plt.figure(dpi=200)
i = 11
plt.contourf(np.log10(np.nansum(new_mad_hi,axis=0)*1.3*1.82e18),
             cmap='magma',levels=np.linspace(20.6,21.5,20),
            extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                   cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50])

plt.contourf(np.nansum(mad_hi_dend[342].get_mask(),axis=0),levels=np.arange(0.5,10),
             extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                   cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50])


# plt.colorbar(ticks = np.arange(20.6,21.5,0.1))

plt.contour(np.log10(h2_hisz_0s[i].data),levels=np.arange(20.5,21.68,0.1),
            extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                   cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],colors='black',
           linewidths=0.2)

plt.contour(np.log10(h2_0s[i].data)*h2_msks[i],levels=np.arange(20.5,21.68,0.1),
            extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                   cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50],colors='black',linewidths=0.9)


plt.gca().invert_xaxis()


mad_hi_dend[21].children


mad_hi_dend[259].parent


mad_hi_dend[59].descendants


dir(mad_hi_dend[259])


def lbd_to_v_map(l,b,d):
    cosb = np.cos(b*u.deg)
    sinl = np.sin(l*u.deg)
    ro = 8.5*u.kpc
    vo = 220*u.km/u.s
    cosl = np.cos(l*u.deg)
    
    
    neg_term = vo*sinl
    num = ro*sinl*vo
    den1 = (d*u.pc - ro*cosl)**2
    den2 = (ro*sinl)**2
    den = np.sqrt(den1+den2)
    pos_term = num/den
    
    return (pos_term - neg_term)*cosb


np.array(cat['name']).transpose()[:3].transpose()


#these are the velocities that we would expect gas at that distance, towards those viewing angles, 
#to be at if 1) the rotation curve of the galaxy is flat, 
#2) is at rest relative to the rotation curve
#so all the gas that is emitting far way from these locations is gas that is either 1) not at rest relative to
#the rotation curve, or2) is located at a different distance
#(the rotation curve thing would affect both?)
#so if we know how to add back in the gas that really is at that location, but not at rest relative to
#the bulk, then we can decide which velocity channels to integrate over
#(assuming that there isnt gas that is somewhere else, but emitting at these velocities by being in motion 
#relative to the gas near *it*)
#okay so the gas that "is at that location but moving relative to it" still doesnt tell you whether

np.array((cat['name'],
          np.round(lbd_to_v_map(np.nanmean(cat['lon_GMC<deg>'],axis=1),
             np.nanmean(cat['lat_GMC<deg>'],axis=1),
             cat['distance<pc>']).to(u.km/u.s).value,3),
         np.round(lbd_to_v_map(np.nanmean(cat['lon_GMC<deg>'],axis=1),
             np.nanmean(cat['lat_GMC<deg>'],axis=1),
             cat['distance<pc>']-cat['radius_HI<pc>']).to(u.km/u.s).value,3),
          np.round(lbd_to_v_map(np.nanmean(cat['lon_GMC<deg>'],axis=1),
             np.nanmean(cat['lat_GMC<deg>'],axis=1),
             cat['distance<pc>']+cat['radius_HI<pc>']).to(u.km/u.s).value,3),
         )).transpose()


all_ls,all_bs,all_ds = np.array(np.meshgrid(np.arange(-180,190,10),
            np.arange(-31,33,5),np.arange(100,1300,100))).transpose().reshape(-1,3).transpose()


all_ls,all_bs,all_ds = np.array(np.meshgrid(np.arange(-180,190,7),
            np.arange(-31,33,3),np.arange(100,1300,5))).transpose().reshape(-1,3).transpose()


all_vs = lbd_to_v_map(all_ls,all_bs,all_ds).to(u.km/u.s)


all_ls.shape


fig = plt.figure(dpi=200)

ax = fig.add_subplot(projection='3d')

ax.scatter(all_ds,all_ls,all_vs,s=1e-3)
