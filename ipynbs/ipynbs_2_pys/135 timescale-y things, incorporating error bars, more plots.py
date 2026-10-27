#!/usr/bin/env python
# coding: utf-8

# # priorities !
# make the scatter plots
# make the histograms / ratio scatter plots
# find the sfr and one (1) other way of timescaling - can be rotation, crossing time, collapse
# virial alpha and rotational beta


import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.io import fits

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle
import matplotlib.cm as cm

import fun as fns
import datacube as dc


from scipy.optimize import curve_fit as fit

def line(x,m,b):
    return m*x + b


cat = QTable(fits.open('sighhplease.fits')[1].data)

c0s = fits.open('h2_0s.fits')
c1s = fits.open('h2_1s.fits')
# c2s = fits.open('h2_2s.fits')

msks = fits.open('h2_0thmom_masks_DISPLAY.fits')
display_msks = [(msks[i].data - 1).astype('bool') for i in range(22)]


h0s = fits.open('hi_0s.fits')
h1s = fits.open('hi_1s.fits')
# h2s = fits.open('hi_2s.fits')


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


delphi = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']
for i in range(22):
    if delphi[i] > 180:
        delphi[i] -=360
        
    if delphi[i] < -180:
        delphi[i] +=360


seg_spec


fig, ax = plt.subplots(1,1,dpi=200)
ax.add_patch(Rectangle((-60,0,),120,5,alpha=0.2,color=seg_spec(12)))
ax.add_patch(Rectangle((-180,0,),60,5,alpha=0.2,color=seg_spec(25)))
ax.add_patch(Rectangle((120,0),60,5,alpha=0.2,color=seg_spec(25)))


plt.hist(delphi,bins = np.arange(-180,210,30),edgecolor='black',color=[0.20,0.20,0.20,0.95])#seg_spec(18))
plt.vlines([-60,60,-120,120],0,5,
          color='black',linestyle='--',linewidth=0.75)

plt.xticks(np.arange(-180,210,60),labels = [str(tck)+r'$\degree$' for tck in np.arange(-180,210,60)])
plt.xlabel(r'$\mathrm{\Delta}\phi$')
plt.xlabel(r'$\phi_{_{GMC}} - \phi_{_{HI}}$')
# plt.ylabel(r'$\mathrm{N_{clouds}}$')
plt.ylabel(r'N (clouds per bin)')


plt.ylim(0,5)
plt.xlim(-180,180)

# plt.savefig('delphi_hist',bbox_inches='tight')
plt.show()


5760/8*1


(np.arange(5760)*-0.125 + 360)[720:2160]


bigh1 = fits.open('bigh1.fits')[0].data
h1_gc = bigh1[:,2160:3600]
h1_agc = bigh1[:,720:2160]


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
                                    raveldf['v']).transpose()[0]
    
    xmsig, ymsig, xybsig = np.sqrt(np.diag(fit(plane_rvld,
                                               np.array((xin,yin)),
                                               raveldf['v']).transpose()[1]))


    if not full:
        return xslope,yslope,xyintercept
    else:
        return {'fit_results':[xslope,yslope,xyintercept],
                '1_sigma_errors':[xmsig, ymsig, xybsig]}


def fit(func, xdata, ydata):
    res = curve_fit(func,xdata,ydata)
#     print('res0:',res[0],'\n\n\n')

    fit_results = res[0]
    fit_errors = np.sqrt(np.diag(res[1]))

    return np.array([fit_results,fit_errors]).transpose()


def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


from scipy.optimize import curve_fit


[(i,np.linalg.norm(fit_gradient(h1s[i].data,
                             cat['distance<pc>'][i])))
for i in [0]]


[(i,(fit_gradient(c1s[i],
                             cat['distance<pc>'][i],
                             (msks[i].data-1).astype('bool'))[:2]),
 cat['dvdx_GMC'][i],cat['dvdy_GMC'][i])
for i in [2,10]]


np.linalg.norm([-0.023534621121459764, 0.09260196684933897])


cat['vel_grad_GMC<km/s/pc>']


fns.fit_gradient(h1_gc,1000)


plt.contourf(h1_gc,
             levels=np.arange(-50,60,5),extnet='both',cmap=seg_spec)


plt.contourf(h1_agc,
             levels=np.arange(-50,60,5),extnet='both',cmap=seg_spec)


from matplotlib.colors import ListedColormap
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


fns.fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['radius_GMC<pc>']))[0]


fns.fit(line,np.log10(cat['radius_GMC<pc>']),np.log10(cat['mass_GMC<Msun>']))[0]


fns.fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['radius_HI<pc>']))[0]


fns.fit(line,np.log10(cat['radius_HI<pc>']),np.log10(cat['mass_HI<Msun>']))[0]


fns.fit(line,np.log10(cat['radius_GMC<pc>']),np.log10(cat['sigma_v_GMC<km/s>']))[0]


fns.fit(line,np.log10(cat['radius_HI<pc>']),np.log10(cat['sigma_v_HI<km/s>']))[0]


fns.fit(line,np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['j_HI<cm2/s>']))[0]


fns.fit(line,np.log10(cat['J_GMC<erg.s>']),np.log10(cat['J_HI<erg.s>']))[0]


fns.fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['mass_GMC<Msun>']))[0]


fns.fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['mass_HI<Msun>']))[0]


plt.hist(cat['mass_HI<Msun>']/cat['mass_GMC<Msun>'],bins = np.arange(0,20))

plt.vlines(1,0,5,color='red')

plt.show()


plt.hist(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'],bins = np.arange(0,20,))

plt.vlines(1,0,5,color='red')

plt.show()


plt.scatter(cat['radius_GMC<pc>'],cat['sigma_v_GMC<km/s>'])
plt.scatter(cat['radius_HI<pc>'],cat['sigma_v_HI<km/s>'])

plt.loglog()


plt.scatter(cat['mass_GMC<Msun>'],cat['radius_GMC<pc>'])
plt.scatter(cat['mass_HI<Msun>'],cat['radius_HI<pc>'])

plt.loglog()


yso_table = QTable(fits.open('J_A+A_674_A21_kyso.dat.gz.fits')[1].data)


big_c0 = fits.open('big_zeroth_co.fits')[0].data


#the next few cells are just doing sanity checks on the data


yso_table.colnames


np.unique(yso_table['SFR']) #198 diff star forming regions identified
#what i now want to do is match these with the clouds in my catalog


plt.hist(yso_table['angDist']*u.arcsec.to(u.deg),bins=[10.**n for n in np.arange(-15,4)])
plt.yscale('log')
plt.vlines([0.125,0.083],0,1e4)
plt.xscale('log')
plt.show()

#ang dist is the displacement / difference between the gaia location of the yso and the 
#(?previously documented) ref location


yso_table['plx'] #hmm weird that there are parallax values that are <0 and = 0


yso_table_og = yso_table
yso_table = yso_table_og[yso_table_og['plx']>0.]


yso_dists = 1/(yso_table['plx']*1e-3)


plt.figure(dpi=200)
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50])
plt.scatter(yso_table['GLON'],yso_table['GLAT'],s=1e-1,color='red')
plt.scatter(yso_table['GLON']-360,yso_table['GLAT'],s=1e-1,color='red')

plt.xlim(180,0)


plt.figure(dpi=200)
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50])
plt.scatter(yso_table['GLON'],yso_table['GLAT'],s=1e-1,color='red')
plt.scatter(yso_table['GLON']-360,yso_table['GLAT'],s=1e-1,color='red')

plt.xlim(0,-180)


dust = fits.open('mean_and_std_lbd.fits')[1].data


dld = np.nansum(dust*2.8,axis=1)*0.23


plt.figure(dpi=200)
plt.contourf(dld,extent=[180,-180,68,1250],levels=np.arange(5e-2,1.,3e-2))
plt.scatter(yso_table['GLON'],yso_dists,s=1e-1,color='red')
plt.scatter(yso_table['GLON']-360,yso_dists,s=1e-1,color='red')

plt.xlim(180,0)
plt.yscale('log')
plt.ylim(0,1250)


plt.figure(dpi=200)
plt.contourf(dld,extent=[180,-180,68,1250],levels=np.arange(5e-2,1.,3e-2))

plt.scatter(yso_table['GLON'],yso_dists,s=1e-1,color='red')
plt.scatter(yso_table['GLON']-360,yso_dists,s=1e-1,color='red')

plt.xlim(180,-180)
plt.yscale('log')
plt.ylim(0,2e3)


nearby = yso_table[yso_dists<=1250]


plt.figure(dpi=200)
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50])


plt.scatter(yso_table['GLON'],yso_table['GLAT'],s=4e-1,color='pink')
plt.scatter(yso_table['GLON']-360,yso_table['GLAT'],s=4e-1,color='pink')

plt.scatter(nearby['GLON'],nearby['GLAT'],s=1e-1,color='red')
plt.scatter(nearby['GLON']-360,nearby['GLAT'],s=1e-1,color='red')


plt.xlim(180,-180)


yso_ls = yso_table['GLON']
for i in range(len(yso_ls)):
    if yso_ls[i] > 180:
        yso_ls[i] -=360


in_cat = 0
co_kysos = np.zeros(22)
for i in range(22):
#     plt.figure(dpi=200)
    
    plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50])


    plt.scatter(yso_table['GLON'],yso_table['GLAT'],s=1,color='red')
    
    plt.xlim(np.array((cat['loh'][i],cat['lih'][i]))*-0.125 + 360)
    plt.ylim(np.array((cat['boh'][i],cat['bih'][i]))*0.125 - 50)

    plt.gca().set_aspect('equal')
    
    
    in_window = yso_table[(yso_table['GLON'] < cat['loh'][i]*-0.125 + 360) & \
    (yso_table['GLON'] > cat['lih'][i]*-0.125 + 360) & \
    (yso_table['GLAT'] > cat['boh'][i]*0.125 -50) & \
    (yso_table['GLAT'] < cat['bih'][i]*0.125 -50)]
    
    on_co = yso_table[(yso_table['GLON'] < cat['loc'][i]*-0.125 + 360) & \
    (yso_table['GLON'] > cat['lic'][i]*-0.125 + 360) & \
    (yso_table['GLAT'] > cat['boc'][i]*0.125 -50) & \
    (yso_table['GLAT'] < cat['bic'][i]*0.125 -50)]
    
    print(len(in_window),len(on_co))
    in_cat +=len(on_co)
    co_kysos[i] = len(on_co)
    plt.show()
    

in_cat


halpha = fits.open('../not_rn_6/wham-ss-DR1-v161116-170912-grid.fits')[0]


hA0 = np.nansum(halpha.data,axis=0)


(np.arange(1440)-720)*-0.25


(np.arange(720)-360)*0.25


plt.figure(dpi=200)
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),
             extent=[360,-360,-50,50],alpha=0.4)

plt.contour(hA0,levels=np.arange(2,80,2),
            extend='max',extent=[180,-180,-90,90],
            linewidths=0.5,cmap='binary_r',alpha=0.7)

plt.scatter(yso_table['GLON'],yso_table['GLAT'],s=4e-1,color='red')
plt.scatter(yso_table['GLON']-360,yso_table['GLAT'],s=4e-1,color='red')

# plt.colorbar()
plt.xlim(180,-180)

# plt.xlim(-140,-180)
# plt.ylim(-30,0)

# plt.xlim(80,20)
# plt.ylim(-20,40)


hAlv = np.nansum(halpha.data,axis=1)


(np.arange(162)-1)*2.03 - 169


halpha.header


plt.contourf(hAlv,levels=np.arange(2,150,5),extend='max',extent=[180,-180,-170,150])


#hmm okay the halpha stuff was a bit of a side quest, but might be worht htinking about in 
#order to get yso / star formation information 
#also interesting that near hercules there is v little halpha, v little ysos
#compared to eg orion,oph where theres a lot of both,& taupercal where theres a lot of yso


fns.fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>']))


plt.scatter(cat['mass_GMC<Msun>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])

plt.loglog()


plt.xlabel(r'$M_{_{GMC}} \ [\mathrm{M_\odot}$]')
plt.ylabel(r'$j_{_{HI}}\ /\ j_{_{GMC}}$')

plt.hlines(1,1e2,3e6,linestyle='--')
plt.xlim(6e2,2e6)
plt.show()


fns.fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>']))


plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])

plt.loglog()


plt.xlabel(r'$M_{_{HI}} \ [\mathrm{M_\odot}$]')
plt.ylabel(r'$j_{_{HI}}\ /\ j_{_{GMC}}$')

plt.hlines(1,9e3,2e7,linestyle='--')
plt.xlim(9e3,1e7)
plt.show()


fns.fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['J_HI<erg.s>']/cat['J_GMC<erg.s>']))


plt.scatter(cat['mass_GMC<Msun>'],cat['J_HI<erg.s>']/cat['J_GMC<erg.s>'])

plt.loglog()


plt.xlabel(r'$M_{_{GMC}} \ [\mathrm{M_\odot}$]')
plt.ylabel(r'$J_{_{HI}}\ /\ J_{_{GMC}}$')

plt.hlines(1,1e2,3e6,linestyle='--')
plt.xlim(6e2,2e6)
plt.show()


fns.fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['J_HI<erg.s>']/cat['J_GMC<erg.s>']))


plt.scatter(cat['mass_HI<Msun>'],cat['J_HI<erg.s>']/cat['J_GMC<erg.s>'])

plt.loglog()


plt.xlabel(r'$M_{_{HI}} \ [\mathrm{M_\odot}$]')
plt.ylabel(r'$J_{_{HI}}\ /\ J_{_{GMC}}$')

plt.hlines(1,9e3,2e7,linestyle='--')
plt.xlim(9e3,1e7)
plt.show()


fns.fit(line,np.log10(cat['sigma_v_GMC<km/s>']),np.log10(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>']))


plt.scatter(cat['sigma_v_GMC<km/s>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
plt.loglog()


plt.xlabel(r'$\sigma_{_{GMC}} \ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$j_{_{HI}}\ /\ j_{_{GMC}}$')

plt.hlines(1,8e-1,1e1,linestyle='--')
plt.xlim(9e-1,9e0)
plt.show()


fns.fit(line,np.log10(cat['sigma_v_HI<km/s>']),np.log10(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>']))


plt.scatter(cat['sigma_v_HI<km/s>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
plt.loglog()


plt.xlabel(r'$\sigma_{_{HI}} \ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$j_{_{HI}}\ /\ j_{_{GMC}}$')

plt.plot(np.arange(4e0,3e1,1e1),np.arange(4e0,3e1,1e1)**2.96 * 10**(-2.08),
         label='$j_{_{HI}}\ /\ j_{_{GMC}}\ \propto \ \sigma_{_{HI}} ^{2.96 \pm 0.75}$')

plt.hlines(1,4e0,3e1,linestyle='--')
plt.xlim(5e0,2e1)
plt.ylim(5e-1,4e2)

plt.legend(loc='upper left')
plt.show()


fns.fit(line,np.log10(cat['sigma_v_HI<km/s>']),np.log10(cat['j_HI<cm2/s>']))


fns.fit(line,np.log10(cat['sigma_v_HI<km/s>']),np.log10(cat['jsum_HI<cm2/s>']))


fns.fit(line,np.log10(cat['sigma_v_GMC<km/s>']),np.log10(cat['j_GMC<cm2/s>']))


fns.fit(line,np.log10(cat['sigma_v_GMC<km/s>']),np.log10(cat['jsum_GMC<cm2/s>']))


fig, axs = plt.subplots(2,2,dpi=200,figsize=(10,10))


axs[0,0].scatter(cat['sigma_v_GMC<km/s>'],cat['j_GMC<cm2/s>'])
axs[0,0].loglog()


axs[0,0].set_xlabel(r'$\sigma_{_{GMC}} \ [\mathrm{km\ s^{-1}}]$')
axs[0,0].set_ylabel(r'$j_{_{GMC}}\ [\mathrm{cm^2\ s^{-1}}]$')

axs[0,0].plot(np.arange(8e-1,1e1,2e0),np.arange(8e-1,1e1,2e0)**2.61 * 10**(23.75),
         label='$j_{_{GMC}}\  \propto \ \sigma_{_{GMC}} ^{\ 2.61 \pm 0.34}$')

axs[0,0].set_xlim(9e-1,9e0)
axs[0,0].legend(loc='upper left')
axs[0,0].set_xticks(np.arange(1,7),labels=[1,2,3,4,5,6])
axs[0,0].set_xlim(0,7)


axs[0,1].scatter(cat['sigma_v_HI<km/s>'],cat['j_HI<cm2/s>'])
axs[0,1].loglog()


axs[0,1].set_xlabel(r'$\sigma_{_{HI}} \ [\mathrm{km\ s^{-1}}]$')
axs[0,1].set_ylabel(r'$j_{_{HI}}\ [\mathrm{cm^2\ s^{-1}}]$')

axs[0,1].plot(np.arange(4e0,3e1,1e1),np.arange(4e0,3e1,1e1)**2.85 * 10**(22.93),
         label='$j_{_{HI}}\  \propto \ \sigma_{_{HI}} ^{\ 2.85 \pm 0.54}$')

axs[0,1].set_xlim(5e0,2e1)
axs[0,1].legend(loc='upper left')
axs[0,1].set_xticks([6,10,15,20],labels=[6,10,15,20])
axs[0,1].set_ylim(3e24,2e27)


axs[1,0].scatter(cat['sigma_v_GMC<km/s>'],cat['jsum_GMC<cm2/s>'])
axs[1,0].loglog()
axs[1,0].set_xticks(np.arange(1,7),labels=[1,2,3,4,5,6])
axs[1,0].set_xlim(0,7)
axs[1,0].plot(np.arange(8e-1,1e1,2e0),np.arange(8e-1,1e1,2e0)**2.83 * 10**(23.34),
         label='$j_{_{GMC}}\  \propto \ \sigma_{_{GMC}} ^{\ 2.83 \pm 0.37}$')

# axs[1,0].plot(np.arange(8e-1,1e1,2e0),np.arange(8e-1,1e1,2e0)**2.61 * 10**(23.75),
#               linestyle='--')
# #          label='$j_{_{GMC}}\  \propto \ \sigma_{_{GMC}} ^{\ 2.61 \pm 0.34}$')

axs[1,0].legend(loc='upper left')
axs[1,0].set_xlabel(r'$\sigma_{_{GMC}} \ [\mathrm{km\ s^{-1}}]$')
axs[1,0].set_ylabel(r'$j_{_{GMC}}\ [\mathrm{cm^2\ s^{-1}}]$')
axs[1,0].set_ylim(1e23,)


axs[0,0].set_title(r'$j \propto \Omega R^2$',x=1.2,y=1.2)
axs[1,0].set_title(r'$j\ =\ \frac{\sum m r \times v}{\sum m}$',x=1.2,y=1.2)


axs[1,1].scatter(cat['sigma_v_HI<km/s>'],cat['j_HI<cm2/s>'])
axs[1,1].loglog()


axs[1,1].set_xlabel(r'$\sigma_{_{HI}} \ [\mathrm{km\ s^{-1}}]$')
axs[1,1].set_ylabel(r'$j_{_{HI}}\ [\mathrm{cm^2\ s^{-1}}]$')

axs[1,1].plot(np.arange(4e0,3e1,1e1),np.arange(4e0,3e1,1e1)**2.90 * 10**(22.54),
         label='$j_{_{HI}}\  \propto \ \sigma_{_{HI}} ^{\ 2.90 \pm 0.49}$')

axs[1,1].set_xlim(5e0,2e1)
axs[1,1].legend(loc='upper left')
axs[1,1].set_xticks([6,10,15,20],labels=[6,10,15,20])
axs[1,1].set_ylim(3e24,2e27)


# plt.subplots_adjust(wspace=0.7,hspace=1)
fig.tight_layout()
plt.show()


plt.scatter(cat['sigma_v_GMC<km/s>'],cat['vel_grad_GMC<km/s/pc>'])
plt.loglog()
plt.show()


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))


axs[0].hist(0.4*cat['j_GMC<cm2/s>']/cat['jsum_GMC<cm2/s>'],bins=np.arange(0,3,0.1),alpha=0.5,edgecolor='white')
axs[1].hist(0.4*cat['j_HI<cm2/s>']/cat['jsum_HI<cm2/s>'],bins=np.arange(0,3,0.1),alpha=0.5,edgecolor='white')

axs[0].annotate('GMC',(2.2,3.5))
axs[1].annotate('HI',(1.6,3.5))


for ax in axs:
    ax.set_xticks(np.arange(0.5,4.5,0.5))
    ax.vlines(1,0,4,linestyle='--')
    
    
axs[0].set_yticks(np.arange(0,5)) 
axs[1].set_yticks(np.arange(0,5),labels=[]) 


axs[0].set_xlim(0,3)
axs[1].set_xlim(0,2)

fig.supylabel('# clouds per bin')
fig.supxlabel(r'$j_{\Omega R^2}\ /\ j_{r \times v}$')

fig.tight_layout()
# plt.vlines(0.1,0,1)
# plt.hist(cat['jsum_HI<cm2/s>']/cat['j_HI<cm2/s>'],bins=np.arange(0.1,2,0.1),alpha=0.5)

plt.show()


sfrs_old = QTable(fits.open('../not_rn_7/re_expanded_cat.fits')[1].data)['SFR<Msun/yr>']


nysos_marton = np.around(sfrs_old*1e6/0.25) 
#the sfr above was calculated from yso counts reported in
#marton, table 2
#marton 2019 yso catalog is here: https://vizier.cds.unistra.fr/viz-bin/VizieR-4


np.array((cat['name'],nysos_marton,co_kysos)).transpose() 
#the kyso is the training sample (blue dots in marton fig 3)


#theres also table 2 in lada,lombardi,alves


plt.hist((1/(cat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)).to(u.Myr).value,
         bins=10**(np.arange(0,2.3,0.1)),alpha=0.5,edgecolor='white')

plt.xscale('log')

plt.xlabel(r'$\tau_{rot} \ [\mathrm{Myr}]$')

plt.xlim(2e0,1e2)
plt.show()


plt.scatter(sfrs_old/cat['mass_GMC<Msun>'],cat['j_GMC<cm2/s>'])
plt.loglog()

plt.xlabel('star formation efficiency')
plt.ylabel('j GMC')
plt.show()


plt.scatter(sfrs_old/cat['mass_GMC<Msun>'],cat['j_HI<cm2/s>'])
plt.loglog()

plt.xlabel('star formation efficiency')
plt.ylabel('j HI')
plt.show()


plt.scatter(sfrs_old/cat['mass_GMC<Msun>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
plt.loglog()

plt.xlabel('star formation efficiency')
plt.ylabel('j HI / j GMC')
plt.show()


plt.scatter(sfrs_old,cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
plt.loglog()
plt.show()


plt.figure(dpi=200,figsize=(25,10))
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),
             extent=[360,-360,-50,50],alpha=0.6)

# plt.contour(hA0,levels=np.arange(2,80,2),
#             extend='max',extent=[180,-180,-90,90],
#             linewidths=0.5,cmap='binary_r',alpha=0.7)

plt.scatter(yso_table['GLON'],yso_table['GLAT'],s=4e-1,color='red')
plt.scatter(yso_table['GLON']-360,yso_table['GLAT'],s=4e-1,color='red')

# plt.colorbar()
plt.xlim(180,-180)
plt.ylim(-30,30)

# plt.xlim(-140,-180)
# plt.ylim(-30,0)

# plt.xlim(80,20)
# plt.ylim(-20,40)


plt.figure(dpi=200,figsize=(25,10))
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),
             extent=[360,-360,-50,50],alpha=0.6,cmap='binary_r')

plt.contour(hA0,levels=np.arange(2,100,4),
            extend='max',extent=[180,-180,-90,90],
            linewidths=0.5,cmap='magma')

# plt.colorbar()
plt.xlim(180,0)
plt.ylim(-30,30)

# plt.xlim(-140,-180)
# plt.ylim(-30,0)

# plt.xlim(80,20)
# plt.ylim(-20,40)


plt.figure(dpi=200,figsize=(25,10))
plt.contourf(np.log10(big_c0),levels=np.arange(20,23,0.25),
             extent=[360,-360,-50,50],alpha=0.6,cmap='binary_r')

plt.contour(hA0,levels=np.arange(2,100,4),
            extend='max',extent=[180,-180,-90,90],
            linewidths=0.5,cmap='magma')

# plt.colorbar()
plt.xlim(0,-180)
plt.ylim(-30,30)

# plt.xlim(-140,-180)
# plt.ylim(-30,0)

# plt.xlim(80,20)
# plt.ylim(-20,40)


QTable(fits.open('tab_2360_1_II_360_catalog.moc.fits')[1].data)


for i in range(22):
    print(h1s[i].data.size,2e-2*h1s[i].data.size)


for i in range(22):
    plt.hist(np.ravel(h1s[i].data),bins=int(1e-2*h1s[i].data.size),histtype='step',density=True)
    
    plt.hist(np.ravel(np.where(~display_msks[i],np.nan,c1s[i].data)),
             bins=int(1e-2*c1s[i].data.size),histtype='step',density=True)
    
    plt.yscale('log')
    plt.show()


hihi = fits.open('hi_hisized.fits')


kwrds = {'linewidth':0.9,'linestyle':':'}


for i in range(22):
    n = hihi[i].data.shape[1]
    
    plt.plot(np.nanmean(hihi[i].data[:,:int(n/10),:int(n/10)],axis=(1,2)),label='1/10',**kwrds)
#     plt.plot(np.nanmean(hihi[i].data[:,:int(n/9),:int(n/9)],axis=(1,2)),label='1/9',linewidth=0.7)
    plt.plot(np.nanmean(hihi[i].data[:,:int(n/8),:int(n/8)],axis=(1,2)),label='1/8',**kwrds)
#     plt.plot(np.nanmean(hihi[i].data[:,:int(n/7),:int(n/7)],axis=(1,2)),label='1/7',linewidth=0.7)
    plt.plot(np.nanmean(hihi[i].data[:,:int(n/6),:int(n/6)],axis=(1,2)),label='1/6',**kwrds)
#     plt.plot(np.nanmean(hihi[i].data[:,:int(n/5),:int(n/5)],axis=(1,2)),label='1/5',linewidth=0.7)
    plt.plot(np.nanmean(hihi[i].data[:,:int(n/4),:int(n/4)],axis=(1,2)),label='1/4',**kwrds)
#     plt.plot(np.nanmean(hihi[i].data[:,:int(n/3),:int(n/3)],axis=(1,2)),label='1/3',linewidth=0.7)
    plt.plot(np.nanmean(hihi[i].data[:,:int(n/2),:int(n/2)],axis=(1,2)),label='1/2',**kwrds)
    plt.plot(np.nanmean(hihi[i].data,axis=(1,2)),**kwrds,label='1')

    plt.legend()
    plt.show()


for i in range(22):
    n = hihi[i].data.shape[1]
    
    
#     plt.plot(np.nanmean(hihi[i].data[:,:int(n/10),:int(n/10)],axis=(1,2)),**kwrds)
# #     plt.plot(np.nanmean(hihi[i].data[:,int(n/10):2*int(n/10),int(n/10):2*int(n/10)],axis=(1,2)),**kwrds)
#     plt.plot(np.nanmean(hihi[i].data[:,2*int(n/10):3*int(n/10),2*int(n/10):3*int(n/10)],axis=(1,2)),**kwrds)
# #     plt.plot(np.nanmean(hihi[i].data[:,4*int(n/10):5*int(n/10),4*int(n/10):5*int(n/10)],axis=(1,2)),**kwrds)
#     plt.plot(np.nanmean(hihi[i].data[:,9*int(n/10):,9*int(n/10):],axis=(1,2)),**kwrds)

    plt.plot(np.nanmax(hihi[i].data,axis=(1,2)),color='black',linewidth=0.7,label='max')
    plt.plot(np.nanmean(hihi[i].data,axis=(1,2)),color='red',linewidth=0.7,label='mean')

    plt.plot(np.nanmedian(hihi[i].data,axis=(1,2)),color='green',linewidth=0.7,label='median')

    plt.legend()
    plt.show()


def fit_gradient_deg(mom1_map, mask=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
    
    
    xdeg = np.indices(mom1_map.shape)[1]*-0.125+360
    ydeg = np.indices(mom1_map.shape)[0]*0.125-50
    
#     xpc = xdeg*u.deg.to(u.rad)*dist
#     ypc = ydeg*u.deg.to(u.rad)*dist
    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map.data)
    
    raveldf = pd.DataFrame(np.array((
            np.ravel(xdeg),np.ravel(ydeg),
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
                                raveldf['v']).transpose()[0]
    

    xmsig, ymsig, xybsig = np.sqrt(np.diag(fit(plane_rvld,
                                               np.array((xin,yin)),
                                               raveldf['v']).transpose()[1]))


    if not full:
        return xslope,yslope,xyintercept
    else:
        return {'fit_results':[xslope,yslope,xyintercept],
                '1_sigma_errors':[xmsig, ymsig, xybsig]}


ays,bees = np.array([fit_gradient_deg(h1s[i].data) for i in range(22)]).transpose()[:2]


np.array((cat['name'],np.linalg.norm([ays,bees],axis=0))).transpose()


hivdf = pd.DataFrame(np.array((cat['name'],
                               cat['vel_grad_HI<km/s/pc>'],
                               np.linalg.norm([ays,bees],axis=0))).transpose(),
             columns=['n','opc','odeg'])


hivdf['opc'].astype('float64')


hivdf['opc'] = hivdf['opc'].astype('float64')*1e3
hivdf['odeg'] = hivdf['odeg'].astype('float64')*1e3


hivdf.sort_values('odeg')


np.nanstd(cat['vel<km/s>'],axis=1)*2


plt.scatter(cat['vel_grad_HI<km/s/pc>']*1e3,cat['distance<pc>'])
plt.loglog()

for i in range(22):
    plt.annotate(cat['name'][i][:3],
                 (cat['vel_grad_HI<km/s/pc>'][i]*1e3,cat['distance<pc>'][i]))
    

plt.xlim(1e1,)


coco = fits.open('h2_h2sized.fits')


plt.plot(np.nanmean(hihi[20].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[20].data,axis=(1,2))))
plt.plot(np.nanmean(hihi[4].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[4].data,axis=(1,2))))
plt.plot(np.nanmean(hihi[11].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[11].data,axis=(1,2))))
# plt.plot(np.nanmean(hihi[12].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[12].data,axis=(1,2))))


plt.plot(np.nanmean(hihi[10].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[10].data,axis=(1,2))))
plt.plot(np.nanmean(hihi[21].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[21].data,axis=(1,2))))
plt.plot(np.nanmean(hihi[14].data,axis=(1,2))/np.nanmax(np.nanmean(hihi[14].data,axis=(1,2))))


# plt.plot(100*np.nanmean(coco[10].data,axis=(1,2)))

# plt.plot(np.nanmean(hihi[11].data,axis=(1,2)))


plt.figure(dpi=200)
plt.scatter(cat['vel_grad_HI<km/s/pc>'],(cat['radius_HI<pc>']/cat['distance<pc>'])*u.rad.to(u.deg))


plt.figure(dpi=200)

plt.scatter(cat['vel_grad_HI<km/s/pc>'],(cat['radius_HI<pc>']/cat['distance<pc>'])*u.rad.to(u.deg))


for i in range(22):
    plt.annotate(cat['name'][i][:3],
                 (cat['vel_grad_HI<km/s/pc>'][i],
                 (cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg)))


plt.figure(dpi=200)
plt.scatter((cat['radius_HI<pc>']/cat['distance<pc>'])*u.rad.to(u.deg),
            np.nanstd(cat['vel<km/s>'],axis=1)*2)

plt.loglog()


oldcat = QTable(fits.open('../not_rn_8/re_expanded_cat.fits')[1].data)


fit(line,cat['vel_grad_HI<km/s/pc>'],np.nanstd(cat['vel<km/s>'],axis=1)*2)


fit(line,cat['vel_grad_GMC<km/s/pc>'],np.nanstd(cat['vel<km/s>'],axis=1)*2)


# plt.scatter(cat['vel_grad_GMC<km/s/pc>'],np.nanstd(cat['vel<km/s>'],axis=1)*2)
plt.figure(dpi=200)

plt.scatter(cat['vel_grad_HI<km/s/pc>'],np.nanstd(cat['vel<km/s>'],axis=1)*2)
# plt.loglog()
# plt.scatter(oldcat['vel_grad_HI<km/s/pc>'],np.nanstd(oldcat['vel<km/s>'],axis=1)*2)
    
# for i in range(22):
#     plt.annotate(oldcat['name'][i][:3],
#                  (oldcat['vel_grad_HI<km/s/pc>'][i],np.nanstd(oldcat['vel<km/s>'][i])*2))
#     plt.annotate(cat['name'][i][:3],
#                  (cat['vel_grad_HI<km/s/pc>'][i],np.nanstd(cat['vel<km/s>'][i])*2))
# plt.plot(cat['vel_grad_HI<km/s/pc>'],
#         cat['vel_grad_HI<km/s/pc>']*240+46)


plt.figure(dpi=200)
plt.scatter(cat['vel_grad_HI<km/s/pc>'],np.nanstd(cat['vel<km/s>'],axis=1)*2)

for i in range(22):
    plt.annotate(cat['name'][i][:3],
                 (cat['vel_grad_HI<km/s/pc>'][i],np.nanstd(cat['vel<km/s>'][i])*2))


plt.contourf(np.nanmax(hihi[0].data,axis=0))
plt.colorbar()


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


vbnds = np.array([[-30,50],
                [-80,30],
                [-50,30],
                [-60,50],
                [-30,100],
                [-20,30],
                [-30,80],
                [-40,50],
                [-100,80],
                  [-60,40],
                  [-100,50],
                  [-20,90],
                  [-30,70],
                  [-30,80],
                  [-50,50],
                  [-30,50],
                  [-30,70],
                  [-30,40],
                  [-50,30],
                  [-50,30],
                  [-20,80],
                  [-70,40]
                  
                 ])


seg_spec


len(seg_spec.colors)


for i in range(22):
    fig,ax=plt.subplots(1,1,dpi=200)

    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(hihi[i].data,axis=(1,2)),linestyle=':',linewidth=3,
            alpha=0.8,color='black')

    ax.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
            np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),
            linewidth=4,color='black')


    if i in [20,12,11]:
        ax.vlines(np.array([cat['vo'][i],cat['vi'][i]-1])*1.3-100,
                  0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),
                  linestyle='--',linewidth=0.7,color='red')

    elif i in [8,4]:
        ax.vlines(np.array([cat['vo'][i],cat['vi'][i]-1])*1.3-100,
                  0,0.6*np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),
                  linestyle='--',linewidth=0.7,color='red')
    else:
        ax.vlines(np.array([cat['vo'][i],cat['vi'][i]-1])*1.3-100,
                  0,0.4*np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),
                  linestyle='--',linewidth=0.7,color='red')

    ax.set_xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

    ax2 = ax.twinx()

    ax2.plot(np.arange(155)*1.3-100,np.nanmean(coco[i].data,axis=(1,2)),
             color=seg_spec(19),linewidth=1.5)


    if i in [4,5,6,12,13,14,20]:
        ax.set_title(names[i],loc='left',y=0.87,x=0.7,fontsize='small')
    else:
        ax.set_title(names[i],loc='left',y=0.87,x=0.05,fontsize='small')


    ax.set_ylabel('$\mathrm{T_{B,HI}\ [K]}$')
    ax.set_xlim(vbnds[i])
    ax.set_yticks(np.around(np.linspace(0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),5),1),
                   labels = np.around(np.linspace(0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),5),1))
    
    ax2.set_ylabel('$\mathrm{T_{B,H_2}\ [K]}$',color=seg_spec(19))
    ax2.set_yticks(np.around(np.linspace(0,np.nanmax(np.nanmean(coco[i].data,axis=(1,2))),5),3),
                   labels = np.around(np.linspace(0,np.nanmax(np.nanmean(coco[i].data,axis=(1,2))),5),3),
                   color=seg_spec(19))
    plt.savefig('spectrum'+lbls[i],bbox_inches='tight')

    plt.show()


len(goodman_oms),len(goodman_rads)


goodman_oms = np.array([1.7,0.32,0.5,0.61,0.52,2.89,3.9,0.6,2.00,2.01,0.58,
              1.05,1.43,0.88,2.20,2.19,3.32,0.30,1.52,0.51,1.06,2.00,
               0.4,1.31,1.58,0.71,0.12,1.70,0.50,0.40,1.90,0.40,0.80,0.27,
              0.95,0.66,0.87,0.66,0.71,0.17,1.33,3.48,0.40])

goodman_mass = np.array([35,20,140,2.1,0.8,0,4.0,0.,0.,0.,6.4,2.6,1.2,14,1.4,4.1,13,
                2.4,20,0.7,150,48,83,6.1,0.,53,21,33,0.6,4.8,8,39,12,24,
                160,19,170,8.2,49,760,32,0,3.3])

goodman_rads = np.array([[0.42,0.10],
              [0.24,0.10],
              [0.72,0.29],
              [0.09,0.06],
              [0.07,0.05],
              [0.29,0.08],
              [0.23,0.13],
              [0.10,0.08],
              [0.13,0.09],
              [0.14,0.08],
              [0.16,0.08],
              [0.09,0.05],
              [0.10,0.07],
              [0.31,0.14],
              [0.09,0.07],
              [0.18,0.11],
              [0.09,0.07],
              [0.25,0.09],
              [0.07,0.05],
              [0.53,0.25],
              [0.75,0.29],
              [0.10,0.07],
              [0.17,0.04],
              [0.19,0.10],
              [0.19,0.10],
              [0.11,0.05],
              [0.38,0.18],
              [0.09,0.07],
              [0.13,0.10],
              [0.13,0.11],
              [0.17,0.16],
              [0.19,0.13],
              [0.36,0.20],
              [0.29,0.22],
              [0.90,0.38],
              [0.29,0.24],
              [0.35,0.20],
              [0.20,0.18],
              [0.35,0.30],
              [1.05,0.73],
              [0.41,0.12],
              [0.36,0.18],
              [0.15,0.08],])


np.linalg.norm(goodman_rads,axis=1)


fit(line,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['j_GMC<cm2/s>']))


fit(line,np.log10(cat['mass_HI<Msun>']),np.log10(cat['j_HI<cm2/s>']))


fit(line,np.log10(goodman_mass[:23])[~np.isinf(np.log10(goodman_mass[:23]))],
    np.log10(((goodman_oms[:23]*\
              np.linalg.norm(goodman_rads,axis=1)**2)*\
              u.km/u.s*u.pc).cgs.value)[~np.isinf(np.log10(goodman_mass[:23]))])


chenx_om = [6.6,12.5,5.8,10.3,24.2,11.2,8.8,22.9]
chenx_rad = np.array([1703,1232,3530,1440,1272,1329,2762,1774])*u.AU.to(u.pc)
chenx_mass = [0.42,0.30,1.21,0.37,0.65,0.57,0.71,0.72]


chenx_rad


((chenx_rad**2 * chenx_om)*u.km*u.pc/u.s).cgs


tatematsu_mass = [0,22,39,43,10,47,41,0,47,57,9,77,
                  39,21,117,0,0,0,0,0,13,0,59,41,0,0,
                  20,0,0,6,0,67,11,0]
tatematsu_rads = np.array([49,88,83,84,59,99,91,49,98,95,52,104,
                          94,83,111,105,74,57,72,43,62,74,115,
                          83,84,62,74,81,62,52,83,157,57,72])*1e-3

tatemastsu_oms = [2.5,1.8,0.9,0.8,0,0,0.5,1.2,5.2,7.8,2.9,0.9,1.7,
                  1.6,2.5,1.4,1.4,3.3,0,0,3.6,3.3,1.2,5.2,0,1.6,0,1.4,
                 2.5,3.3,3.5,1.2,2.2,0]


tatematsu_rads


plt.figure(dpi=200)
plt.scatter(goodman_mass,
            ((goodman_oms*\
              np.linalg.norm(goodman_rads,axis=1)**2)*u.km/u.s*u.pc).cgs,
            label='cores')


plt.scatter(cat['mass_GMC<Msun>'],cat['j_GMC<cm2/s>'],label='$\mathrm{H_2}$ clouds')
plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>'],label='atomic envelopes')

plt.legend()
plt.loglog()

plt.xlabel('mass [$\mathrm{M_\odot}$]')
plt.ylabel('j [$\mathrm{cm^2\ s^{-1}}$]')


# plt.xlim(1e-1,1e3)

plt.savefig('jM_withcores',bbox_inches='tight')
plt.show()


plt.figure(dpi=200)
plt.scatter(np.linalg.norm(goodman_rads,axis=1),
            ((goodman_oms*\
              np.linalg.norm(goodman_rads,axis=1)**2)*u.km/u.s*u.pc).cgs,label='cores')

plt.scatter(cat['radius_GMC<pc>'],cat['j_GMC<cm2/s>'],label='$\mathrm{H_2}$ clouds')
plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],label='atomic envelopes')

plt.legend()
plt.loglog()

plt.xlabel('radius [pc]')
plt.ylabel('j [$\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('jR_withcores',bbox_inches='tight')


plt.show()


plt.figure(dpi=200)
plt.scatter(0.75*goodman_mass/np.pi/np.linalg.norm(goodman_rads,axis=1)**3,
            ((goodman_oms*\
              np.linalg.norm(goodman_rads,axis=1)**2)*u.km/u.s*u.pc).cgs,
            label='cores')

plt.scatter(0.75*cat['mass_GMC<Msun>']/np.pi/cat['radius_GMC<pc>']**3,
            cat['j_GMC<cm2/s>'],label='$\mathrm{H_2}$ clouds')

plt.scatter(0.75*cat['mass_HI<Msun>']/np.pi/cat['radius_HI<pc>']**3,
            cat['j_HI<cm2/s>'],label='atomic envelopes')


plt.legend()
plt.loglog()

plt.xlabel(r'$\rho\ [\mathrm{pc^3}]$')
plt.ylabel('j [$\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('jR_withcores',bbox_inches='tight')


plt.show()


fit(line,
   np.concatenate((
       np.log10(cat['mass_GMC<Msun>']), np.log10(cat['mass_HI<Msun>']),
       np.log10(goodman_mass[:23])[~np.isinf(np.log10(goodman_mass[:23]))]
   )),
    np.concatenate((
        np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['j_HI<cm2/s>']),
        np.log10(((goodman_oms[:23]*\
              np.linalg.norm(goodman_rads,axis=1)**2)*\
              u.km/u.s*u.pc).cgs.value)[~np.isinf(np.log10(goodman_mass[:23]))]
    ))
   )


fit(line,
   np.concatenate((
       np.log10(cat['radius_GMC<pc>']), np.log10(cat['radius_HI<pc>']),
       np.log10(np.linalg.norm(goodman_rads[:23],axis=1)),
   )),
    np.concatenate((
        np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['j_HI<cm2/s>']),
        np.log10(((goodman_oms[:23]*\
              np.linalg.norm(goodman_rads,axis=1)**2)*\
              u.km/u.s*u.pc).cgs.value)
    ))
   )


lms = np.array([[75, 83], #0
             [63, 86],
             [70, 84], #2
             [74, 85],
             [80, 113], #4
             [76, 88],
             [80, 95], #6
             [65, 105],
             [75, 92], #8
             [70, 85],
             [72, 90], #10
             [80, 110],
             [70, 90], #12
             [78, 93],
             [68, 90], #14
             [75, 90],
             [75, 90], #16
             [70, 90],
             [70, 90], #18
             [68, 83],
             [75, 100], #20
             [75, 90]])


for i in [11]:
    sp = np.nanmean(coco[i].data,axis=(1,2))
    frst = np.gradient(sp)
    scnd = np.gradient(frst)
    
    scndgg0 = np.where(scnd>0)[0]
    plt.plot(sp)
    plt.plot(1e-2*np.nanmean(hihi[i].data,axis=(1,2)))
    plt.vlines([lms[i][0],lms[i][-1]],
               
               0,np.nanmax(sp),color='red')
    
    plt.vlines(cat['vo'][11],
               
               0,np.nanmax(sp),color='green')
    plt.hlines(0.1*np.nanmax(sp),0,155)
    plt.title(cat['name'][i])
#     plt.ylim(-0.5*np.nanmax(),0.5*np.nanmax(frst))
    plt.ylim(-0.1*np.nanmax(sp),0.2*np.nanmax(sp))
    
    plt.show()


new_h1s = [
    fns.get_mom1(hihi[i].data[lms[i][0]:lms[i][1],:,:],
                 np.arange(lms[i][0],lms[i][1])*1.3-100) for i in range(22)
]


dir(plt.colorbar(new))


for i in range(22):
    fig,axs = plt.subplots(1,2,dpi=200,figsize=(20,10))
    
    
    new = axs[1].contourf(new_h1s[i],cmap=seg_spec,levels=15)

    cb = fig.colorbar(new,ax=axs[1])
    mn,mx = cb.vmin,cb.vmax
    
    old = axs[0].contourf(h1s[i].data,
                          levels=np.linspace(mn,mx,15),cmap=seg_spec,
                         extend='both')
    
#     axs[0].contour(np.where(~display_msks[i],np.nan,c1s[i].data),colors='black',levels=10)
#     fig.colorbar(old,ax=axs[0])
    
    fig.tight_layout()
    fig.suptitle(cat['name'][i])
    plt.show()


lms[11]


cat['vo'][11],cat['vi'][11]


for i in range(22):
    plt.hist(np.ravel(h1s[i].data),bins=100,alpha=0.5,density=True)
    plt.hist(np.ravel(new_h1s[i].data),bins=100,alpha=0.5,density=True)
    plt.hist(np.ravel(np.where(~display_msks[i],np.nan,c1s[i].data)),
             bins=100,alpha=0.5,density=True,histtype='step')
    plt.title(cat['name'][i]
             )
    plt.show()


for i in range(22):
    n = h1s[i].data.shape[1]
    
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
    
    
# #     background_ax = plt.axes([0, 0, 1, 1]) # create a dummy subplot for the background
# # background_ax.set_zorder(-1) # set the background subplot behind the others
# # background_ax.imshow(image, aspect='auto') # show the backgroud image
    
    
#     bkg_map = fig.add_subplot(111)#add_axes([0, 0, 1, 1])
#     bkg_map.set_zorder(-10)
#     bkg_map.set_axis_off()
#     bkg_map.contourf(h1s[i].data,cmap=seg_spec,extend='both',
#                     extent=[loh,lih,boh,bih],
#                     linewidths=0.4,zorder=-20,alpha=0.1,
#                     levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26))
   
    
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
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        
        spold = np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],bo:bi,lo:li],axis=(1,2))
        spnew = np.nanmean(hihi[i].data[lms[i][0]:lms[i][1],bo:bi,lo:li],axis=(1,2))
        
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        spold/np.nanmax(spold),linewidth=2,
                        )
        
        map_axs[j].plot(np.arange(lms[i][0],lms[i][1])*1.3-100,
                        spnew/np.nanmax(spnew),linewidth=0.8)
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
#         map_axs[j].vlines(np.nanmedian(h1s[i].data),
#                           0,1.5,color='black',
#                           linestyle='--',linewidth=1.5,alpha=0.2)
        
#         map_axs[j].scatter(np.nanmean(fns.get_mom1(his[i][:,bo:bi,lo:li],
#                                    np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
#                            s=7,color='black',alpha=0.3,linewidth=0)
        
        
#         if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
#             map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
#                             (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')
            
            
#             map_axs[j].scatter(np.nanmean(fns.get_mom1(coshsz[i][:,bo:bi,lo:li],
#                                    np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,s=1,color='black')
#             map_axs[j].vlines(np.nanmedian(np.where(~msks[i],np.nan,c1s[i].data)),
#                               0,1.5,color='black',linestyle='--',linewidth=0.4)
        
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_ylim(-0.15,1.15)
        map_axs[j].set_axis_off()
    
    
    plt.subplots_adjust(wspace=0.07,hspace=0.07)
    plt.figaspect(1)
        
    fig.suptitle(names[i],y=0.93)

    plt.show()


newas,newbs = np.array([fit_gradient(new_h1s[i],cat['distance<pc>'][i]) for i in range(22)]).transpose()[:2]


np.array((cat['name'],np.linalg.norm([newas,newbs],axis=0)/cat['vel_grad_HI<km/s/pc>'])).transpose()


np.linalg.norm([newas,newbs],axis=0)/cat['vel_grad_HI<km/s/pc>']


newomg = np.linalg.norm([newas,newbs],axis=0)


plt.scatter(2*np.nanstd(lms,axis=1)/(cat['vi']-cat['vo']),
   newomg/cat['vel_grad_HI<km/s/pc>'])
            
# plt.loglog()


plt.scatter(cat['radius_HI<pc>'],cat['vel_grad_HI<km/s/pc>'],s=30)
plt.scatter(cat['radius_HI<pc>'],newomg,s=10)

# plt.loglog()


plt.scatter(cat['mass_HI<Msun>'],cat['vel_grad_HI<km/s/pc>'],s=30)
plt.scatter(cat['mass_HI<Msun>'],newomg,s=10)
plt.xscale('log')
# plt.loglog()


plt.hist(np.nanstd(cat['vel<km/s>'],axis=1))
plt.hist(np.nanstd(lms*1.3-100,axis=1),alpha=0.5)

plt.show()


dist_errs = np.array([[2,7],
            [2,24],
            [8,11],
            [3,9],
            [4,60],
            [3,8],
            [4,44],
            [4,89],
            [1,11],
            [5,25],
            [9,9],
            [6,104],
            [3,37],
            [3,39],
            [2,7],
            [2,22],
            [2,21],
            [1,20],
            [2,15],
            [1,18],
            [5,65],
            [2,7]])

norm_d_errs = np.linalg.norm(dist_errs,axis=1)


mass_err_GMC = 2*norm_d_errs/cat['distance<pc>']*cat['mass_GMC<Msun>']
mass_err_HI = 2*norm_d_errs/cat['distance<pc>']*cat['mass_HI<Msun>']


rad_err_gmc = cat['radius_GMC<pc>']*norm_d_errs/cat['distance<pc>']


rad_err_hi = 1/3*mass_err_GMC/cat['mass_GMC<Msun>']*cat['radius_HI<pc>']


omg_err_GMC = np.linalg.norm((cat['dvdx_GMC']*cat['sigma_dvdx_GMC'],
                             cat['dvdy_GMC']*cat['sigma_dvdy_GMC']),axis=0)/cat['vel_grad_GMC<km/s/pc>']


omg_err_HI = np.linalg.norm((cat['dvdx_HI']*cat['sigma_dvdx_HI'],
                             cat['dvdy_HI']*cat['sigma_dvdy_HI']),axis=0)/cat['vel_grad_HI<km/s/pc>']


j_err_GMC = cat['j_GMC<cm2/s>']*np.linalg.norm((omg_err_GMC/cat['vel_grad_GMC<km/s/pc>'], 
                                                2*rad_err_gmc/cat['radius_GMC<pc>']),axis=0)


j_err_HI = cat['j_HI<cm2/s>']*np.linalg.norm((omg_err_GMC/cat['vel_grad_HI<km/s/pc>'], 
                                                2*rad_err_hi/cat['radius_HI<pc>']),axis=0)


help(plt.errorbar)


help(curve_fit)


goodman_rads_eff =  np.linalg.norm(goodman_rads,axis=1)


goodman_oms_err = np.array([0.09,0.07,0.10,0.30,0.56,0.12,0.20,0.30,0.53,0.50,0.20,
                   0.14,0.41,0.31,0.60,0.38,0.62,0.07,0.57,0.25,0.16,0.70,0.70,0.12,0.16,0.38,
                   0.12,0.80,0.40,0.10,0.20,0.50,0.20,0.15,0.09,0.21,0.32,
                   0.18,0.18,0.09,0.07,0.11,0.30
                  ])


goodman_js = ((goodman_oms*goodman_rads_eff**2)*u.km/u.s*u.pc).cgs.value
goodman_js_err = ((goodman_oms_err*goodman_rads_eff**2)*u.km/u.s*u.pc).cgs.value


all_ms = np.concatenate((goodman_mass,cat['mass_GMC<Msun>'],cat['mass_HI<Msun>']),axis=0)
all_rs = np.concatenate((goodman_rads_eff,cat['radius_GMC<pc>'],cat['radius_HI<pc>']),axis=0)
all_js = np.concatenate((goodman_js,cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>']),axis=0)


fit(line,np.log10(all_rs)[~np.isinf(np.log10(all_ms))],
    np.log10(all_js)[~np.isinf(np.log10(all_ms))])


plt.figure(dpi=200)

plt.errorbar(cat['radius_GMC<pc>'],cat['j_GMC<cm2/s>'],
             xerr=rad_err_gmc, yerr=j_err_GMC,fmt='.',elinewidth=0.5,
             label=r'molecular cloud',color=seg_spec(19))

plt.errorbar(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],
             xerr=rad_err_hi, yerr=j_err_HI,fmt='.',elinewidth=0.5,
             label='atomic envelope',color='black')

plt.errorbar(goodman_rads_eff,goodman_js,
            yerr=goodman_js_err,
             label='core [G+93]',fmt='.',elinewidth=0.5,color=seg_spec(6))

#if needing to change size etc in plt.errorbar: ms=markersize, mew=markeredgewidth


plt.loglog()

plt.plot(np.arange(1e-2,1e3,1e1),np.arange(1e-2,1e3,1e1)**1.46 * 10**(23.05),color='grey',
        linewidth=0.7)

# plt.scatter(all_rs,all_js,alpha=0.4)


plt.legend(fontsize='x-small')
plt.loglog()

plt.xlabel('$R$ [pc]')
plt.ylabel('$j$ [$\mathrm{cm^2\ s^{-1}}$]')


plt.xlim(5e-2,3e2)
plt.ylim(3e20,1e27)

plt.savefig('jR_withcores',bbox_inches='tight')
plt.show()


fit(line,np.log10(all_ms)[~np.isinf(np.log10(all_ms))],
    np.log10(all_js)[~np.isinf(np.log10(all_ms))])


all_js_err = np.concatenate((goodman_js_err,j_err_GMC,j_err_HI))


plt.figure(dpi=200)


plt.errorbar(cat['mass_GMC<Msun>'],cat['j_GMC<cm2/s>'],
             xerr=rad_err_gmc, yerr=j_err_GMC,fmt='.',elinewidth=0.5,
             label=r'molecular cloud',color=seg_spec(19))

plt.errorbar(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>'],
             xerr=rad_err_hi, yerr=j_err_HI,fmt='.',elinewidth=0.5,
             label='atomic envelope',color='black')


plt.errorbar(goodman_mass,goodman_js,
            yerr=goodman_js_err,
             label='core [G+93]',fmt='.',elinewidth=0.5,color=seg_spec(6))

plt.loglog()


plt.plot(np.arange(1e-1,1e8,1e1),np.arange(1e-1,1e8,1e1)**0.77 * 10**(21.39),color='grey',
        linewidth=0.7)

plt.legend(fontsize='x-small')
plt.loglog()

plt.xlabel('$M$ [$\mathrm{M_\odot}$]')
plt.ylabel('$j$ [$\mathrm{cm^2\ s^{-1}}$]')

plt.gca().set_xticks(ticks=[1e-1,1e1,1e3,1e5,1e7],
                     minor=True)
# plt.gca().set_xticks([1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7],minor=False)
# 
plt.xlim(3e-1,1e7)
plt.ylim(3e20,1e27)
# plt.gca().set_aspect('equal')

plt.savefig('jM_withcores',bbox_inches='tight')
plt.show()


curve_fit(line,np.log10(all_ms)[~np.isinf(np.log10(all_ms))],
    np.log10(all_js)[~np.isinf(np.log10(all_ms))],sigma=all_js_err[~np.isinf(np.log10(all_ms))])


np.sqrt(0.0007269)


help(curve_fit)


fit(line,np.log10(cat['sigma_v_GMC<km/s>']),np.log10(cat['j_GMC<cm2/s>']))


plt.figure(dpi=200)


plt.errorbar(cat['sigma_v_GMC<km/s>'],cat['j_GMC<cm2/s>'],
             yerr=j_err_GMC,fmt='.',elinewidth=0.5,
             color=seg_spec(19))


plt.loglog()


plt.plot(np.arange(1e0,1e1,1e0),np.arange(1e0,1e1,1e0)**2.61 * 10**(23.75),color='grey',
        linewidth=0.7)

plt.fill_between(np.arange(1e0,1e1,1e0),
                 np.arange(1e0,1e1,1e0)**(2.61-0.36) * 10**(23.75-0.17),
                 np.arange(1e0,1e1,1e0)**(2.61+0.36) * 10**(23.75+0.17),
                 color='grey',linewidth=0.,alpha=0.2,zorder=-20)

# plt.legend(fontsize='x-small')
plt.loglog()

plt.xlabel('$\sigma_{_\mathrm{GMC}}$ [$\mathrm{km\ s^{-1}}$]')
plt.ylabel('$j_{_\mathrm{GMC}}$ [$\mathrm{cm^2\ s^{-1}}$]')

plt.gca().set_xticks(ticks=[1e0,2e0,3e0,4e0,6e0],labels=[1e0,2e0,3e0,4e0,6e0])
 
plt.xlim(1e0,7e0)
plt.ylim(3e23,2e26)
# plt.gca().set_aspect('equal')

plt.show()


plt.scatter(cat['sigma_v_GMC<km/s>'], cat['j_GMC<cm2/s>'])
plt.scatter(cat['sigma_v_HI<km/s>'], cat['j_HI<cm2/s>'])

plt.loglog()


plt.figure(dpi=200)
plt.hist(cat['vel_grad_HI<km/s/pc>']/cat['vel_grad_GMC<km/s/pc>'],
         bins=10**np.arange(-3.,2,0.2),
        edgecolor='black',color=[0.20,0.20,0.20,0.95])

plt.vlines(1,0,5,color=seg_spec(4),linewidth=0.75,linestyle='--')

plt.xscale('log')
plt.ylabel('N (clouds per bin)')

plt.xlabel(r'$\Omega_{_\mathrm{HI}}\ /\ \Omega_{_\mathrm{H_2}}$')

plt.xlim(1e-2,15)
plt.ylim(0,5)
plt.show()


cat[cat['vel_grad_HI<km/s/pc>']>cat['vel_grad_GMC<km/s/pc>']]


plt.figure(dpi=200)
plt.hist(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'],
         bins=10**np.arange(-2.,3,0.2),
        edgecolor='black',color=[0.20,0.20,0.20,0.95])

plt.vlines(1,0,5,color=seg_spec(4),linewidth=0.75,linestyle='--')

plt.xscale('log')
# plt.ylabel('N (clouds per bin)')

# plt.xlabel(r'$\Omega_{_\mathrm{HI}}\ /\ \Omega_{_\mathrm{H_2}}$')

plt.xlim(5e-2,300)
plt.ylim(0,6)
plt.show()


cat[cat['j_GMC<cm2/s>']>cat['j_HI<cm2/s>']]


cat[cat['J_GMC<erg.s>']>cat['J_HI<erg.s>']]


plt.scatter(np.arange(22),np.around(delphi),s=5)

# plt.hlines(np.arange(-180,210,30),0,22)
plt.xticks(np.arange(0,22,2))
plt.yticks(np.arange(-120,180,60))

plt.grid()


#sorting delphi to pick which ones to talk about
algn = [0,3,4,7,8,11,14,15,19]
antalgn = [2,9,12,13,20]
meh = [1,5,6,10,16,17,18,21]


1,5,6    13,2,9   3,14,7


plt.figure(dpi=200)
plt.hist(cat['J_HI<erg.s>']/cat['J_GMC<erg.s>'],
         bins=10**np.arange(-2.,4,0.2),
        edgecolor='black',color=[0.20,0.20,0.20,0.95])

plt.vlines(1,0,5,color=seg_spec(4),linewidth=0.75,linestyle='--')

plt.xscale('log')
# plt.ylabel('N (clouds per bin)')

# plt.xlabel(r'$\Omega_{_\mathrm{HI}}\ /\ \Omega_{_\mathrm{H_2}}$')

# plt.xlim(5e-2,300)
plt.ylim(0,6)
plt.show()


# plt.scatter(cat[ 'sigma_v_GMC<km/s>'], cat['vel_grad_GMC<km/s/pc>'])
plt.scatter(cat['sigma_v_HI<km/s>'], cat['vel_grad_HI<km/s/pc>'])

plt.loglog()


fit(line,np.log10(cat['sigma_v_GMC<km/s>']),np.log10(cat['j_GMC<cm2/s>']))


fit(line,np.log10(cat['sigma_v_HI<km/s>']),np.log10(cat['j_HI<cm2/s>']))


0.54/2.85


0.36/2.61


#making the big catalog figures vv


bigh0 = fits.open('bigh0.fits')[0].data


plt.figure(dpi=200,figsize=(10,5))
plt.contourf(np.log10(bigh0),extent=[360,-360,-50,50],
             levels=np.arange(20.2,22.5,0.2),extend='both',cmap='binary_r')

plt.contour(np.log10(big_c0),extent=[360,-360,-50,50],
            levels=np.arange(19.8,23),colors='black',linewidths=0.3)

for i in [0,3,5,8,14]:
    plt.contourf(np.log10(h0s[i].data),extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                                     cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],
                levels=10,cmap='magma',antialiased=True)
    
for i in [10]:
    plt.contourf(np.log10(h0s[i].data),extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                                     cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],
                levels=10,cmap='magma',alpha=0.75,antialiased=True)
    

for i in range(22):
    plt.contourf(np.where(~display_msks[i],np.nan,np.log10(c0s[i].data)),
                 extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                                     cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50],
                levels=10,cmap='YlGnBu_r',zorder=20)
    

plt.xticks(np.arange(60,-100,-20),[str(tck)+'$\degree$' for tck in np.arange(60,-100,-20)])
plt.yticks(np.arange(-30,50,10),[str(tck)+'$\degree$' for tck in np.arange(-30,50,10)])


plt.gca().set_aspect(1.)
plt.xlim(60,-80)
plt.ylim(-30,40)


plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')


plt.savefig('galactic_center_catalog',bbox_inches='tight')
plt.show()


#grrrr figure out what to do with all the overlap vv


plt.figure(dpi=200,figsize=(10,5))
plt.contourf(np.log10(bigh0),extent=[360,-360,-50,50],levels=np.arange(20,22.2,0.2),cmap='binary_r')

plt.contour(np.log10(big_c0),extent=[360,-360,-50,50],
            levels=np.arange(19.8,23),colors='black',linewidths=0.3)


for i in [4,11,2,9,19,21,12,15]:    
    plt.contourf(np.log10(h0s[i].data),extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
                                     cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],
                levels=10,cmap='magma',antialiased=True,zorder=10)
    
    
    plt.contourf(np.log10(h0s[i].data),extent=[cat['loh'][i]*-0.125+720,cat['lih'][i]*-0.125+720,
                                     cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],
                levels=10,cmap='magma',antialiased=True,zorder=10)

   
for i in [1,18,6,13,20,16,17,7]:
    plt.gca().add_patch(Rectangle((cat['loh'][i]*-0.125+360,cat['boh'][i]*0.125-50),
                                  (cat['lih'][i]-cat['loh'][i])*-0.125,
                                  (cat['bih'][i]-cat['boh'][i])*0.125,
                       fill=False,linestyle='--',linewidth=0.8,edgecolor='white',zorder=15))
    
    
    plt.gca().add_patch(Rectangle((cat['loh'][i]*-0.125+720,cat['boh'][i]*0.125-50),
                                  (cat['lih'][i]-cat['loh'][i])*-0.125,
                                  (cat['bih'][i]-cat['boh'][i])*0.125,
                       fill=False,linestyle='--',linewidth=0.8,edgecolor='white',zorder=15))
    
#     plt.contourf(np.log10(h0s[i].data),extent=[cat['loh'][i]*-0.125+360,cat['lih'][i]*-0.125+360,
#                                      cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],
#                 levels=10,cmap='magma',alpha=0.65,antialiased=True)
    
    
#     plt.contourf(np.log10(h0s[i].data),extent=[cat['loh'][i]*-0.125+720,cat['lih'][i]*-0.125+720,
#                                      cat['boh'][i]*0.125-50,cat['bih'][i]*0.125-50],
#                 levels=10,cmap='magma',alpha=0.65,antialiased=True)

    
for i in range(22):
    plt.contourf(np.where(~display_msks[i],np.nan,np.log10(c0s[i].data)),
                 extent=[cat['loc'][i]*-0.125+360,cat['lic'][i]*-0.125+360,
                                     cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50],
                levels=10,cmap='YlGnBu_r',zorder=20)
    
    
    plt.contourf(np.where(~display_msks[i],np.nan,np.log10(c0s[i].data)),
                 extent=[cat['loc'][i]*-0.125+720,cat['lic'][i]*-0.125+720,
                                     cat['boc'][i]*0.125-50,cat['bic'][i]*0.125-50],
                levels=10,cmap='YlGnBu_r',zorder=20)


plt.gca().set_aspect(1.)
plt.xlim(240,90)


goodman_rads.transpose()


gdmn_df = pd.DataFrame(np.array([goodman_mass,goodman_rads_eff,
                       goodman_rads.transpose()[0],goodman_rads.transpose()[1],
                       goodman_oms,goodman_oms_err,
                       goodman_js,goodman_js_err]).transpose(),
            columns=['M<Msun>','R_eff<pc>','major_ax<pc>','minor_ax<pc>',
                    'Omega<km/s/pc>','err_Omega<km/s/pc>','j<cm2/s>','err_j<cm2/s>'])


gdmn_df = gdmn_df[~np.isinf(np.log10(gdmn_df['M<Msun>']))]


gdmn_df.to_csv('goodman_1993_data.csv',index=False)


cat['err_M_GMC<Msun>'] = mass_err_GMC
cat['err_M_HI<Msun>'] = mass_err_HI

cat['err_R_GMC<pc>'] = rad_err_gmc
cat['err_R_HI<pc>'] = rad_err_hi

cat['err_Omg_GMC<km/s/pc>'] = omg_err_GMC
cat['err_Omg_HI<km/s/pc>'] = omg_err_HI

cat['err_j_GMC<cm2/s>'] = j_err_GMC
cat['err_j_HI<cm2/s>'] = j_err_HI


cat.write('sighhplease.fits',overwrite=True)
