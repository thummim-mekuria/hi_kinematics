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

from gaussdecomp import spectrum,fitter


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
       [  10, 72,   4,  25],
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


bos = bos + cat['boc']
bis = bis + cat['boc']

los = los + cat['loc']
lis = lis + cat['loc']

coco = [dht.data[:,bos[i]:bis[i],los[i]:lis[i]] for i in range(22)]


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


sph_gf = 0.75
cyl_gf = 0.5


n_wnm = 0.57/u.cm**3


sph_hi_rad = np.cbrt(sph_gf/np.pi*(h2_mass/(n_wnm*2*1.36*c.m_p))).to(u.pc)
cyl_hi_rad = np.cbrt(cyl_gf/np.pi*(h2_mass/(n_wnm*2*1.36*c.m_p))).to(u.pc)


def L1 (R,b,d):
    first = d**2
    second = 2*R**2
    third = 2*np.sqrt(2)*R*d
    fourth = np.cos(b+45*u.deg)
    
    return np.sqrt(first + second - third*fourth)


def L2 (R,b,d):
    first = d**2
    second = 2*R**2
    third = 2*np.sqrt(2)*R*d
    fourth = np.cos(135*u.deg-b)
    
    return np.sqrt(first + second - third*fourth)


def L3 (R,b,d):
    first = d**2
    second = 2*R**2
    third = 2*np.sqrt(2)*R*d
    fourth = np.cos(45*u.deg-b)
    
    return np.sqrt(first + second - third*fourth)


def b1 (R,b,d):
    this_L = L1(R,b,d)
    
    first = np.sqrt(2)*R
    second = np.sin(45*u.deg+b)
    
    return np.arcsin(first*second/this_L)


def b2 (R,b,d):
    this_L = L2(R,b,d)
    
    first = np.sqrt(2)*R
    second = np.sin(135*u.deg-b)
    
    return np.arcsin(first*second/this_L)


def b3 (R,b,d):
    this_L = L3(R,b,d)
    
    first = np.sqrt(2)*R
    second = np.sin(45*u.deg-b)
        
    return np.arcsin(first*second/this_L)


db_proj_cyl = (b1(cyl_hi_rad,
                  np.abs(np.nanmean(cat['lat_GMC<deg>'],axis=1))*u.deg,
                  cat['distance<pc>']*u.pc,) + \
               
                b2(cyl_hi_rad,
                   np.abs(np.nanmean(cat['lat_GMC<deg>'],axis=1))*u.deg,
                   cat['distance<pc>']*u.pc,)).to(u.deg)


for i in np.arange(22)[cat['distance<pc>']*u.pc*\
                       np.sin(np.abs(np.nanmean(cat['lat_GMC<deg>'],axis=1))*u.deg) < cyl_hi_rad]:
    db_proj_cyl[i] = (b1(cyl_hi_rad[i],
                         np.abs(np.nanmean(cat['lat_GMC<deg>'],axis=1))[i]*u.deg,
                         cat['distance<pc>'][i]*u.pc,) + \
                      
                        b3(cyl_hi_rad[i],
                           np.abs(np.nanmean(cat['lat_GMC<deg>'],axis=1))[i]*u.deg,
                           cat['distance<pc>'][i]*u.pc,)).to(u.deg)


db_proj_cyl


db_fon_cyl = 2*np.arctan(cyl_hi_rad/(cat['distance<pc>']*u.pc)).to(u.deg)


db_sph = 2*np.arctan(sph_hi_rad/(cat['distance<pc>']*u.pc)).to(u.deg)


plt.figure(dpi=200)
plt.plot(100*(db_proj_cyl-db_fon_cyl)/db_fon_cyl,marker='.',label='cylinder in projection')
# plt.plot(100*(db_fon_cyl*np.cos(np.nanmean(cat['lat_GMC<deg>'],axis=1))-db_fon_cyl)/db_fon_cyl,
#               marker='.',label='cos theta correction')

# plt.plot(100*(db_sph-db_fon_cyl)/db_fon_cyl,marker='.',label='sphere of same volume')

plt.hlines(0,0,22,color='black',label='sphere of same radius (smaller volume)')

plt.ylabel(r"$\Delta \theta$ compared to square aperture")
plt.xticks(np.arange(22),cat['name'],rotation='vertical')
plt.yticks(np.arange(10,50,10),[str(tck)+'%' for tck in np.arange(10,50,10)])

plt.legend(fontsize='x-small')

plt.hlines(np.arange(0,45,5),0,21,linewidth=0.5,color='grey')
plt.vlines(np.arange(0,22,3),40,0,linewidth=0.5,color='grey')
# plt.grid()


rad_one = (0.5*db_fon_cyl.value/0.125).astype('int')
rad_two = (0.5*db_proj_cyl.value/0.125).astype('int')
rad_three = (0.5*db_sph.value/0.125).astype('int')


bc = ((np.nanmean(cat['lat_GMC<deg>'],axis=1) + 50)/0.125).astype('int')
lc = ((np.nanmean(cat['lon_GMC<deg>'],axis=1) - 360)/-0.125).astype('int')


hihi_one = [hi4pi.data[:,
                    bc[i]-rad_one[i]:bc[i]+rad_one[i],
                    lc[i]-rad_one[i]:lc[i]+rad_one[i]] for i in range(22)] 
                                                                            

hihi_two = [hi4pi.data[:,
                    bc[i]-rad_two[i]:bc[i]+rad_two[i],
                    lc[i]-rad_one[i]:lc[i]+rad_one[i]] for i in range(22)] #no stretching along longitude
                                                                          # so its still set by just 2R

hihi_three = [hi4pi.data[:,
                    bc[i]-rad_three[i]:bc[i]+rad_three[i],
                    lc[i]-rad_three[i]:lc[i]+rad_three[i]] for i in range(22)]


hihi_one_0 = [np.nansum(hihi_one[i],axis=0)*1.3*1.82e18 for i in range(22)]
hihi_two_0 = [np.nansum(hihi_two[i],axis=0)*1.3*1.82e18 for i in range(22)]
hihi_three_0 = [np.nansum(hihi_three[i],axis=0)*1.3*1.82e18 for i in range(22)]


#percentage difference between the distance to the center of the cloud
#and the distance to the bottom/top (which makes a difference for the delta l)
#appearing narrower when that side of the cylinder is further away
#aka the rectangular part of the cylinder would appear as a trapezoid
#and the bases will be different by this amount from its width at center


plt.figure(dpi=200)

plt.plot(100*(-1 + np.sqrt(1+\
(cyl_hi_rad.value/cat['distance<pc>'])**2 +\
2*(cyl_hi_rad.value/cat['distance<pc>'])*np.sin(np.abs(np.nanmean(cat['lat_GMC<deg>'],axis=1))[i]*u.deg))),
    marker='.')

plt.yticks(np.arange(0,20,5),[str(tck)+'%' for tck in np.arange(0,20,5)])
plt.xticks(np.arange(22),cat['name'],rotation='vertical')

plt.hlines([0,2.5,5,7.5,10,12.5,15,],0,21,linewidth=0.5,color='grey')
plt.vlines(np.arange(0,22,3),0,15,linewidth=0.5,color='grey')


lc[i]-rad_one[i]:lc[i]+rad_one[i]] for i in range(22)] 
                                                        

hihi_two = [hi4pi.data[:,
bc[i]-rad_two[i]:bc[i]+rad_two[i],
lc[i]-rad_one[i]:lc[i]+rad_one[i]] for i in range(22)] #no stretching along longitude
                                                      # so its still set by just 2R

hihi_three = [hi4pi.data[:,
bc[i]-rad_three[i]:bc[i]+rad_three[i],
lc[i]-rad_three[i]:lc[i]+rad_three[i]


[(
    
#cylinder in projection
plt.contourf(np.log10(hihi_two_0[i]),
             extent=[(lc[i]-rad_one[i])*-0.125+360,(lc[i]+rad_one[i])*-0.125+360,
                 (bc[i]-rad_two[i])*0.125-50,(bc[i]+rad_two[i])*0.125-50],
             levels=20,antialiased=True,cmap='magma'),
    
    
# #sphere of same volume
# plt.contourf(np.log10(hihi_three_0[i]),
#              extent=[(lc[i]-rad_three[i])*-0.125+360,(lc[i]+rad_three[i])*-0.125+360,
#                  (bc[i]-rad_three[i])*0.125-50,(bc[i]+rad_three[i])*0.125-50],
#              levels=20,antialiased=True,alpha=0.5),
    
#how we were doing it before
plt.contourf(np.log10(hihi_one_0[i]),
             extent=[(lc[i]-rad_one[i])*-0.125+360,(lc[i]+rad_one[i])*-0.125+360,
                 (bc[i]-rad_one[i])*0.125-50,(bc[i]+rad_one[i])*0.125-50],
             levels=20,antialiased=True,alpha=0.8,cmap='binary_r'),
    
    
    plt.gca().set_aspect(1),
    plt.gca().invert_xaxis(),
    

plt.title(cat['name'][i]),
plt.show()
 ) for i in range(22)]


vbounds_from_decomp = np.array([(-2.5696620038548375, 6.794540065823298),
 (-8.35431147090679, 12.225397329460073),
 (-5.1880729807415, 5.986613118541698),
 (-7.606591186101394, 9.059522305464654),
 (12.568696890607587, 26.500335197566574),
 (1.7695522050555423, 10.18816299588961),
 (6.4523479182481, 17.413561133438147),
 (-9.035450586925403, 17.884797626508885),
 (-3.240593530004398, 14.859948953825553),
 (-3.516239152098073, 7.472361463362857),
 (-11.700973846547026, 14.238208452606635),
 (21.516062921670724, 37.716429333032934),
 (-3.321259577539468, 17.587129629976182),
 (4.963462296285389, 17.09763801543781),
 (-4.800034971675656, 8.614849509765087),
 (-2.444520292714807, 17.39217308781455),
 (3.3332595821594024, 17.560193659891606),
 (-1.1021775697645495, 16.6796743862138),
 (-8.556163199563162, 12.650429526049372),
 (-11.577160122665761, 8.492378422886313),
 (-1.465290058914693, 18.513101235328694),
 (2.0906136538297186, 14.954705746449859)])


[(plt.plot(np.arange(155)*1.3-100,
           np.nanmean(hihi_one[i],axis=(1,2)),color='black',linewidth=0.7,label='sphere of same radius'),
plt.plot(np.arange(155)*1.3-100,
         np.nanmean(hihi_two[i],axis=(1,2)),color='orange',linewidth=0.7,label='cylinder in projection'),
plt.plot(np.arange(155)*1.3-100,
         np.nanmean(hihi_three[i],axis=(1,2)),color='green',linewidth=0.7,label='sphere of same volume'),
  
plt.plot(np.arange(155)*1.3-100,
         100*np.nanmean(coco[i],axis=(1,2)),color='red',label='co'),
  
  plt.vlines(vbounds_from_decomp[i],0,70),
  plt.xlim(vbounds_from_decomp[i]+20*np.array([-1,1])),
  plt.title(cat['name'][i]),
plt.legend(),
plt.show())
for i in range(22)]


vos_gaus, vis_gaus = ((vbounds_from_decomp.transpose()+100)/1.3).astype('int')


hi_mass_one = np.array([fns.mass_n_area(hihi_one[i][vos_gaus[i]:vis_gaus[i],:,:]*u.K,
                               cat['distance<pc>'][i]*u.pc,
                               '21cm')[0].to(u.solMass).value for i in range(22)])

hi_mass_two = np.array([fns.mass_n_area(hihi_two[i][vos_gaus[i]:vis_gaus[i],:,:]*u.K,
                               cat['distance<pc>'][i]*u.pc,
                               '21cm')[0].to(u.solMass).value for i in range(22)])

hi_mass_three = np.array([fns.mass_n_area(hihi_three[i][vos_gaus[i]:vis_gaus[i],:,:]*u.K,
                               cat['distance<pc>'][i]*u.pc,
                               '21cm')[0].to(u.solMass).value for i in range(22)])


plt.figure(dpi=200)
plt.plot(100*(hi_mass_two-hi_mass_one)/hi_mass_one,marker='.',label='cylinder in projection')
plt.plot(100*(hi_mass_three-hi_mass_one)/hi_mass_one,marker='.',label='sphere of same volume')

plt.hlines(0,0,22,color='black',label='sphere of same radius (smaller volume)')

plt.ylabel(r"$\Delta \mathrm{M}$ compared to square aperture")
plt.xticks(np.arange(22),cat['name'],rotation='vertical')
plt.yticks(np.arange(0,70,10),[str(tck)+'%' for tck in np.arange(0,70,10)])

plt.legend(fontsize='x-small')

plt.hlines(np.arange(0,70,10),0,21,linewidth=0.5,color='grey')
plt.vlines(np.arange(0,22,3),60,0,linewidth=0.5,color='grey')
# plt.grid()


hi_mom1_one = [fns.get_mom1(hihi_one[i][vos_gaus[i]:vis_gaus[i],:,:],
                               np.arange(vos_gaus[i],vis_gaus[i])*1.3-100) for i in range(22)]

hi_mom1_two = [fns.get_mom1(hihi_two[i][vos_gaus[i]:vis_gaus[i],:,:],
                               np.arange(vos_gaus[i],vis_gaus[i])*1.3-100) for i in range(22)]

hi_mom1_three = [fns.get_mom1(hihi_three[i][vos_gaus[i]:vis_gaus[i],:,:],
                               np.arange(vos_gaus[i],vis_gaus[i])*1.3-100) for i in range(22)]


hi_mom2_one = [np.sqrt(fns.get_momN(hihi_one[i][vos_gaus[i]:vis_gaus[i],:,:],
                               np.arange(vos_gaus[i],vis_gaus[i])*1.3-100,2)) for i in range(22)]

hi_mom2_two = [np.sqrt(fns.get_momN(hihi_two[i][vos_gaus[i]:vis_gaus[i],:,:],
                               np.arange(vos_gaus[i],vis_gaus[i])*1.3-100,2)) for i in range(22)]

hi_mom2_three = [np.sqrt(fns.get_momN(hihi_three[i][vos_gaus[i]:vis_gaus[i],:,:],
                               np.arange(vos_gaus[i],vis_gaus[i])*1.3-100,2)) for i in range(22)]


for i in range(22):
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5))
    axs[0].hist(np.ravel(hi_mom1_one[i]),bins=100)
    axs[1].hist(np.ravel(hi_mom1_two[i]),bins=100)
    axs[2].hist(np.ravel(hi_mom1_three[i]),bins=100)

    fig.suptitle(cat['name'][i])
    plt.show()


v_lvsl = np.array(([-0.5,3.5],
[-3,6],
[-2,1],
[-3,3],
[16,20],
[4,6],
[9,12],
[1,7],
[3,7],
[-0.5,1.5],
[-4,5],
[26.5,29.5],
[5,9],
[7.5,11],
[-1,3],
[4,9.5],
[7,11],
[5,9],
[-2,6],
[-6,1],
[7,10],
[4,10]))


for i in range(22):
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5))
    axs[0].contourf(hi_mom1_one[i],cmap=fns.segmented_spectral,
                    levels=np.linspace(v_lvsl[i][0],v_lvsl[i][1],27),
                   extend='both')
    
    
    axs[1].contourf(hi_mom1_two[i],cmap=fns.segmented_spectral,
                   levels=np.linspace(v_lvsl[i][0],v_lvsl[i][1],27),
                   extend='both')

    
    axs[2].contourf(hi_mom1_three[i],cmap=fns.segmented_spectral,
                   levels=np.linspace(v_lvsl[i][0],v_lvsl[i][1],27),
                   extend='both')
    
    axs[0].set_title('sphere with same radius')
    axs[1].set_title('cylinder in projection')
    axs[2].set_title('sphere with same volume')

    fig.suptitle(cat['name'][i])
    plt.show()


for i in range(22):
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5))
    axs[0].contourf(hi_mom2_one[i],cmap='magma',
                    levels=30,
                   extend='both')
    
    axs[1].contourf(hi_mom2_two[i],cmap='magma',
                   levels=30,
                   extend='both')
    
    axs[2].contourf(hi_mom2_three[i],cmap='magma',
                   levels=30,
                   extend='both')
    
    axs[0].set_title('sphere with same radius')
    axs[1].set_title('cylinder in projection')
    axs[2].set_title('sphere with same volume')

    fig.suptitle(cat['name'][i])
    plt.show()


def fit_gradient(mom1_map, dist, mask=0, weights=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
        
    if type(weights)==int:#all the data is usable and/or a mask wasnt supplied
        weights = np.ones(mom1_map.shape).astype('bool')
    
    
    xdeg = np.indices(mom1_map.shape)[1]*-0.125+360
    ydeg = np.indices(mom1_map.shape)[0]*0.125-50
    
    xpc = xdeg*u.deg.to(u.rad)*dist
    ypc = ydeg*u.deg.to(u.rad)*dist
    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map.data)
    
    
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


from scipy.optimize import curve_fit as fit


def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


hi_omg_one = np.array([np.linalg.norm(fit_gradient(hi_mom1_one[i],
                                      cat['distance<pc>'][i],
                                      weights=hi_mom2_one[i])[:1])
          for i in range(22)])

#for some reason the weights mess up the fits here
hi_omg_one[14:17] = [np.linalg.norm(fit_gradient(hi_mom1_one[i],cat['distance<pc>'][i])[:1]) 
                     for i in range(14,17)]


hi_omg_two = np.array([np.linalg.norm(fit_gradient(hi_mom1_two[i],
                                      cat['distance<pc>'][i],
                                      weights=hi_mom2_two[i])[:1])
          for i in range(22)])

hi_omg_two[14:17] = [np.linalg.norm(fit_gradient(hi_mom1_two[i],cat['distance<pc>'][i])[:1]) 
                     for i in range(14,17)]


hi_omg_three = np.array([np.linalg.norm(fit_gradient(hi_mom1_three[i],
                                      cat['distance<pc>'][i],
                                      weights=hi_mom2_three[i])[:1])
          for i in range(22)])

#for some reason the weights mess up the fits here
hi_omg_three[14:17] = [np.linalg.norm(fit_gradient(hi_mom1_three[i],cat['distance<pc>'][i])[:1]) 
                       for i in range(14,17)]


#not sure why the gradient for cra is so small the first time
#it goes up by ~2 orders of magnitude

plt.figure(dpi=200)


plt.plot(100*(hi_omg_two-hi_omg_one)/hi_omg_one,marker='.',label='cylinder in projection')
plt.plot(100*(hi_omg_three-hi_omg_one)/hi_omg_one,marker='.',label='sphere of same volume')

plt.hlines(0,0,22,color='black',label='sphere of same radius (smaller volume)')

plt.ylabel(r"$\Delta \mathrm{\Omega}$ compared to square aperture")
plt.xticks(np.arange(22),cat['name'],rotation='vertical')
plt.yticks(np.arange(-50,70,20),[str(tck)+'%' for tck in np.arange(-50,70,20)])

plt.legend(fontsize='x-small')

plt.hlines(np.arange(-50,70,10),0,21,linewidth=0.5,color='grey')
plt.vlines(np.arange(0,22,3),-50,60,linewidth=0.5,color='grey')

plt.ylim(-60,70)
# plt.grid()


plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>'],label='wider velocity range')
plt.scatter(hi_mass_one,(hi_omg_one*u.km/u.s/u.pc*cyl_hi_rad**2).cgs.value,label='narrower velocity range')

for i in range(22):
    plt.plot((hi_mass_one[i],cat['mass_HI<Msun>'][i]),
             ((hi_omg_one*u.km/u.s/u.pc*cyl_hi_rad**2).cgs.value[i],cat['j_HI<cm2/s>'][i]),
             color='black',linewidth=0.5)
plt.loglog()

plt.xlabel('mass [solar mass]')
plt.legend()
plt.ylabel('j [cm2/s]')


plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],label='wider velocity range')
plt.scatter(cyl_hi_rad,(hi_omg_one*u.km/u.s/u.pc*cyl_hi_rad**2).cgs.value,label='narrower velocity range')

for i in range(22):
    plt.plot((cyl_hi_rad[i].value,cat['radius_HI<pc>'][i]),
             ((hi_omg_one*u.km/u.s/u.pc*cyl_hi_rad**2).cgs.value[i],cat['j_HI<cm2/s>'][i]),
             color='black',linewidth=0.5)
plt.loglog()

plt.xlabel('radius [pc]')
plt.legend()
plt.ylabel('j [cm2/s]')


plt.figure(dpi=200)
plt.hist((hi_omg_one*u.km/u.s/u.pc*cyl_hi_rad**2).cgs.value/cat['j_GMC<cm2/s>'],
         alpha=0.5,bins=10**np.arange(-3,4,1/3),label='narrower velocity range',edgecolor='white')

plt.hist(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'],alpha=0.5,bins=10**np.arange(-3,4,1/3),
         label='wider velocity range',edgecolor='white')
plt.xscale('log')
plt.legend()
plt.xlabel('$j_{HI}/ j_{GMC}$')

plt.vlines(1,0,7,color='black',linestyle='--')

plt.xlim(5e-3,5e2)


plt.figure(dpi=200)

plt.plot(hi_omg_one/cat['vel_grad_HI<km/s/pc>'],marker='.',label='$\Omega$')
plt.plot(hi_mass_one/cat['mass_HI<Msun>'],marker='.',label='M')


plt.yscale('log')
plt.legend()

plt.ylabel(r'changes after choosing narrower velocity range')
plt.hlines(10**np.arange(-3,1,0.5),0,21,linewidth=0.5,color='grey')
plt.vlines(np.arange(0,22,3),1e-3,10**(0.5),linewidth=0.5,color='grey')
plt.xticks(np.arange(22),cat['name'],rotation='vertical')

plt.hlines(1e0,0,21,color='black')

plt.ylim(5e-4,5e0)


# so since we are going to go with the spheres, i will need to make circular masks for the zeroth/first moments
# 
# next steps:
# new notebook
# save the sphere radii
# make new cubes
# make masks


sph_hi_rad
