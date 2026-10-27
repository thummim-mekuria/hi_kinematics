#!/usr/bin/env python
# coding: utf-8

# # when i come back:
# - make a new overleaf, copy over the outlne, then populate the sections with 1) stuff ive already written or 2) stream of consciousness rambling (i really mean that, just write without thinking) based on the bullet point
# - do this for the stuff that won't be affected by whether i choose to go with noovs or not
# - then afterwards, write up the section on choosing hi from pv maps -- maybe decompose the cubes?
# - remake the maps, the scaling relations (with goodman + phillips data), find a timescale that plays well with ji
# - do some analysis with aspect ratios/moments of inertia, and jsum
# - even if i dont have a big conclusion, i have a bunch of tiny things to report, so just focus on those
# - done is better than perfect


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

h2_0s = fits.open('h2_0s.fits')
hi_0s = fits.open('hi_0s.fits')


h2_msks = [(fits.open('h2_0thmom_masks_DISPLAY.fits')[i].data - 1).astype('bool') 
           for i in range(22)]

hihi = fits.open('hi_hisized.fits')
coco = fits.open('h2_h2sized.fits')

hi_1s = fits.open('hi_1s.fits')
h2_1s = fits.open('h2_1s.fits')


hi_2s = fits.open('hi_2s.fits')
cos_hisz = fits.open('h2_hisized.fits')

h2_hisz_2s = [np.sqrt(fns.get_momN(cos_hisz[i].data[cat['vo'][i]:cat['vi'][i],:,:],
             np.arange(cat['vo'][i],cat['vi'][i])*1.3,2)) for i in range(22)]


[(
    plt.contourf(np.linalg.norm(np.gradient(np.log10(
        np.nansum(hihi[i].data,axis=2)*0.125*1.82e18)),axis=0).transpose(),
                 cmap='magma',
                 levels=np.linspace(0,0.1,20)),
    plt.colorbar(),
    plt.plot(np.nanmean(hihi[i].data,axis=(1,2))),

    plt.title(cat['name'][i]),
    plt.show()
) for i in range(22)]


[(
    plt.contourf(np.linalg.norm(np.gradient(np.log10(
        np.nansum(hihi[i].data,axis=1)*0.125*1.82e18)),axis=0).transpose(),
                 cmap='magma',
                 levels=np.linspace(0,0.2,20)),
    plt.colorbar(),
    plt.plot(np.nanmean(hihi[i].data,axis=(1,2))),

    plt.title(cat['name'][i]),
    plt.show()
) for i in range(22)]


[(
    plt.contourf(np.linalg.norm(np.gradient(np.log10(
        np.nansum(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=0
                 )*1.3*1.82e18)),axis=0),
                 cmap='magma'),
    plt.colorbar(),

    plt.title(cat['name'][i]),
    plt.show()
) for i in range(22)]


[(
# plt.hist(np.ravel(np.log10(np.nansum(hihi[i].data,axis=0)*1.3*1.82e18)),bins=np.linspace(10,23,100),alpha=0.5),
plt.hist(np.ravel(np.log10(np.nansum(hihi[i].data,axis=1)*0.125*1.82e18)),bins=np.linspace(10,23,100),alpha=0.5),
# plt.hist(np.ravel(np.log10(np.nansum(hihi[i].data,axis=2)*0.125*1.82e18)),bins=np.linspace(10,23,100),alpha=0.5,),
    plt.title(cat['name'][i]),
plt.show()
) for i in range(22)]


[(
print(np.ravel(np.log10(np.nansum(hihi[i].data,axis=0)*1.3*1.82e18)).size),
print(np.ravel(np.log10(np.nansum(hihi[i].data,axis=1)*0.125*1.82e18)).size),
print(np.ravel(np.log10(np.nansum(hihi[i].data,axis=2)*0.125*1.82e18)).size),
    print('\n\n')
#     plt.title(cat['name'][i]),
# plt.show()
) for i in range(22)]


[(
    plt.contourf(np.log10(
        np.nansum(hihi[i].data,axis=2)*0.125*1.82e18).transpose(),
                 cmap='YlGnBu_r',
                 levels=np.linspace(20.6,21.3,20)),
    plt.colorbar(),
    plt.plot(np.nanmean(hihi[i].data,axis=(1,2))),
    plt.vlines((cat['vel<km/s>'][i]+100)/1.3,0,100,color='red'),

    plt.title(cat['name'][i]),
    plt.show()
) for i in range(22)]


co_pks = [fns.get_mom1(np.nanmean(coco[i].data,axis=(1,2)),
                       np.arange(155)*1.3-100)
    for i in range(22)]


[(
    plt.plot(np.arange(155)*1.3-100,70*np.nanmean(coco[i].data,axis=(1,2))),
    plt.plot(np.arange(155)*1.3-100,np.nanmean(hihi[i].data,axis=(1,2))),
    
    plt.vlines([co_pks[i]-11.4,co_pks[i]+11.4],0,100,linestyle=':'),
    plt.vlines(cat['vel<km/s>'][i],0,100),
    plt.vlines(noovs[i],0,100,linestyle='--'),
    plt.title(cat['name'][i]),

    
plt.show()
) for i in range(22)]


[(
    plt.contourf(np.log10(
        np.nansum(hihi[i].data,axis=1)*0.125*1.82e18).transpose(),
                 cmap='YlGnBu_r',
                 levels=np.linspace(20.6,21.3,20),
                extent = [0,155,cat['loh'][i],cat['lih'][i]]),
    
    plt.contourf(np.nansum(coco[i].data,axis=1).transpose(),
                 levels=np.linspace(1,60,10),cmap='magma',
                extent=[0,155,cat['loc'][i],cat['lic'][i]]),
    
    plt.colorbar(),
    plt.plot(cat['loh'][i]+np.nanmean(hihi[i].data,axis=(1,2))),
    plt.vlines((cat['vel<km/s>'][i]+100)/1.3,cat['loh'][i],cat['loh'][i]+100,color='red'),
    plt.title(cat['name'][i]),
    plt.show()
) for i in range(22)]


zc = np.abs(np.tan(np.nanmean(
              cat['lat_GMC<deg>'],axis=1)*u.deg).value)*cat['distance<pc>']


np.array((np.around(np.tan(np.array((zc-cat['radius_HI<pc>'],
                  zc+cat['radius_HI<pc>']))/cat['distance<pc>'])*\
 u.rad.to(u.deg),3),(np.array((cat['boh'],cat['bih']))*0.125 - 50))).transpose()


(np.array((cat['boh'],cat['bih']))*0.125 - 50)


cat


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


nooms = np.array([fns.mass_n_area(hihi[i].data[noovos[i]:noovis[i],:,:]*u.K,
                 cat['distance<pc>'][i]*u.pc,'21cm')[0].to(u.solMass).value
        for i in range(22)])


noomoo1s = [fns.get_mom1(hihi[i].data[noovos[i]:noovis[i],:,:],
                       np.arange(noovos[i],noovis[i])*1.3-100) for i in range(22)]

noomoo2s = [np.sqrt(fns.get_momN(hihi[i].data[noovos[i]:noovis[i],:,:],
                       np.arange(noovos[i],noovis[i])*1.3-100,2)) for i in range(22)]


from scipy.optimize import curve_fit as fit


def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


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


nooomgs = [np.linalg.norm(fit_gradient(noomoo1s[i],cat['distance<pc>'][i],mask=hi_msks[i],
            weights=noomoo2s[i])[:-1]) for i in range(22)]

#idk why but the fit for ophiuchus is going off the rails when there is the weights

nooomgs[14] = np.linalg.norm([(fit_gradient(noomoo1s[i],
                              cat['distance<pc>'][i],
                              mask=hi_msks[i])[:-1]) 
                for i in [14]])


noojs = ((nooomgs*u.km/u.s/u.pc)*(cat['radius_HI<pc>']*u.pc)**2).cgs


t_ff = np.sqrt((np.pi**2 * (cat['radius_GMC<pc>']*u.pc)**3)/(c.G*cat['mass_GMC<Msun>']*u.solMass*8)).cgs


t_rot = (1/(cat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)).cgs


def line(x,m,b):
    return m*x+b


fit(line,np.log10(cat['j_GMC<cm2/s>']),np.log10(t_ff.value))


fit(line,np.log10(cat['radius_GMC<pc>']),np.log10(cat['mass_GMC<Msun>']))


fit(line,np.log10(cat['radius_GMC<pc>']),np.log10(cat['vel_grad_GMC<km/s/pc>']))


plt.scatter(cat['radius_GMC<pc>'],cat['vel_grad_GMC<km/s/pc>'])
plt.loglog()


plt.scatter(cat['j_GMC<cm2/s>'],t_ff)
plt.loglog()


plt.plot(t_ff/t_rot,marker='+')


plt.plot(nooms/cat['mass_GMC<Msun>'],marker='+')
plt.plot(cat['mass_HI<Msun>']/cat['mass_GMC<Msun>'],marker='+')


j_i = (cat['j_GMC<cm2/s>']*cat['mass_GMC<Msun>'] +\
       cat['j_HI<cm2/s>']*cat['mass_HI<Msun>'])/(cat['mass_GMC<Msun>'] + cat['mass_HI<Msun>'] )

j_i_noo = (cat['j_GMC<cm2/s>']*cat['mass_GMC<Msun>'] +\
       noojs.value*nooms)/(cat['mass_GMC<Msun>'] + nooms)


# plt.plot(j_i,marker='.')
plt.plot(j_i_noo,marker='.')


plt.plot(cat['j_GMC<cm2/s>'],marker='.')
# plt.plot(cat['j_HI<cm2/s>'],marker='.')
# plt.plot(noojs,marker='.')


plt.yscale('log')
plt.show()


plt.plot(cat['j_GMC<cm2/s>']/j_i)
plt.plot(cat['j_GMC<cm2/s>']/j_i_noo)

plt.yscale('log')


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


delphi = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']
for i in range(22):
    if delphi[i] < -180:
        delphi[i] += 360
        
    if delphi[i] > 180:
        delphi[i] -= 360


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


phillips['j<cm2/s>'] = (np.array(phillips['Omega<km/s/pc>']*\
                                 phillips['R<pc>']**2)*u.km/u.s*u.pc).cgs.value
