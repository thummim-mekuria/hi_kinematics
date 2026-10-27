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

coco = fits.open('h2_h2sized.fits')
hihi = fits.open('hi_hisized.fits')

msks = fits.open('h2_0thmom_masks.fits')
msks = [(msks[i].data - 1).astype('bool') for i in range(22)]


subco = [coco[i].data[cat['vo'][i]:cat['vi'][i],:,:] for i in range(22)]
subhi = [hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:] for i in range(22)]


c0s = [np.nansum(subco[i],axis=0)*1.3*1.97e20 for i in range(22)]
h0s = [np.nansum(subhi[i],axis=0)*1.3*1.82e18 for i in range(22)]


c0s_long = [np.nansum(coco[i].data,axis=0)*1.3*1.97e20 for i in range(22)]
h0s_long = [np.nansum(hihi[i].data,axis=0)*1.3*1.82e18 for i in range(22)]


#cfor most of the gmcs, the masses vary by <1%
#cma tho has 1.7, mad has 5.4, mob has 5.2 -> look at these guys


np.array((cat['name'],np.around(100*np.array([np.nansum((c0s[i]*msks[i]/u.cm**2*\
 c.m_p*2*1.36*\
 (0.125*u.deg.to(u.rad)*\
  cat['distance<pc>'][i]*u.pc)**2).to(u.solMass)).value

 
 for i in range(22)])/
np.array(cat['mass_GMC<Msun>']),2))).transpose()


[(
plt.plot(np.nanmean(coco[i].data,axis=(1,2))),
  
plt.vlines([cat['vo'][i],cat['vi'][i]],
           0,np.nanmax(np.nanmean(coco[i].data,axis=(1,2))),color='red'),
plt.show())


for i in [4,11,12]]


#hmm only seeing the difference with the mask off ... but mass was supposed to have been calculated with mask on 
#anyway?


[(
plt.contourf(np.log10(c0s_long[i]),levels=np.arange(20,22.5,0.15),cmap='magma'),
plt.contourf(np.log10(c0s[i]),levels=np.arange(20,22.5,0.15),cmap='binary_r',alpha=0.5),
  

plt.show())


for i in [4,11,12]]


[(
plt.contourf(np.log10(c0s_long[i])*msks[i],levels=np.arange(20,22.5,0.15),cmap='magma'),
plt.contourf(np.log10(c0s[i])*msks[i],levels=np.arange(20,22.5,0.15),cmap='binary_r',alpha=0.5),
  

plt.show())


for i in [4,11,12]]


#okay now checking how big a difference it makes for rh


short_gmc = np.array([np.nansum((c0s[i]*msks[i]/u.cm**2*\
 c.m_p*2*1.36*\
 (0.125*u.deg.to(u.rad)*\
  cat['distance<pc>'][i]*u.pc)**2).to(u.solMass)).value

 
 for i in range(22)])


np.array((cat['radius_HI<pc>'],np.cbrt(short_gmc/(0.04*np.pi)))).transpose()[[4,11,12]]


94.55288641/95.0805506


131.58757843/134.00406314


88.16974062/89.74016082


((np.array((cat['radius_HI<pc>'],
           np.cbrt(short_gmc/(0.04*np.pi)))).transpose()/\
  cat['distance<pc>'][:,np.newaxis])[[4,11,12]])*u.rad.to(u.deg)/0.125


#madd and mob get a ~2% difference, which in the grand scheme of things (ie uncertainty in the distance etc)
#isnt that big, so gonna move on
#also ^ the biggest is a one pixel difference so ... 


hi_mass = [fns.mass_n_area(subhi[i]*u.K,cat['distance<pc>'][i]*u.pc,'21cm')[0].to(u.solMass).value 
           for i in range(22)] #the 1e20 thresh doesnt cut off any pixels


cat['mass_HI<Msun>'] = hi_mass


[
 
   (plt.plot(50*np.nanmean(subco[i],axis=(1,2))),
    plt.plot(np.nanmean(subhi[i],axis=(1,2))),
    plt.show())
    
  for i in range(22)  
]


spco = [np.nanmean(subco[i],axis=(1,2)) for i in range(22)]
sphi = [np.nanmean(subhi[i],axis=(1,2)) for i in range(22)]


spcosm = [np.nansum(subco[i],axis=(1,2)) for i in range(22)]
sphism = [np.nansum(subhi[i],axis=(1,2)) for i in range(22)]


sigmavco = np.sqrt(np.array([fns.get_momN(spco[i],
                                  np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2) 
                     for i in range(22)]))


sigmavco


#check that summing vs meaning isnt changing the spectrum
np.sqrt(np.array([fns.get_momN(spcosm[i],
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2) 
                  for i in range(22)]))


sigmavhi = np.sqrt(np.array([fns.get_momN(sphi[i],
                                  np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2) 
                     for i in range(22)]))


sigmavhi


#check that summing vs meaning isnt changing the spectrum
np.sqrt(np.array([fns.get_momN(sphism[i],
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2) 
                  for i in range(22)]))


cat['sigma_v_GMC<km/s>'] = sigmavco
cat['sigma_v_HI<km/s>'] = sigmavhi


cat


c1s = [fns.get_mom1(subco[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100) for i in range(22)]
h1s = [fns.get_mom1(subhi[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100) for i in range(22)]


[(
# plt.figure(dpi=200),
plt.contourf(np.where(~msks[i],np.nan,c1s[i]),levels=25,cmap='PiYG'),
plt.colorbar(),
plt.gca().set_aspect('equal'),

plt.show())


for i in range(22)]


c2s = [np.sqrt(fns.get_momN(subco[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2)) for i in range(22)]
h2s = [np.sqrt(fns.get_momN(subhi[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2)) for i in range(22)]


oldcat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


for i in [0]:
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    print(cdf)
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    print(cdf)
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
#     print(cat['name'][i],np.around(np.array((np.linalg.norm((cxfit[0][0],
#                                               cyfit[0][0])),oldcat['vel_grad_GMC<km/s/pc>'][i])),4))


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    
    cdfold = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(np.where(~msks[i],np.nan,c1s[i])),
            )).transpose(),
            columns=['l','b','v'])
    
#     print(cdf)
    
    cdf_noisy = cdf
     
    
    cdfnew = cdf_noisy[~(np.isnan(cdf_noisy['v']))]
    
    print(cat['name'][i][:3],len(cdfold),len(cdfnew))
    #looking at how diff ways of selecting datapoints for fitting compare
    
#     cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
#     cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
#     cxfit = fit(line,cxpc,cdf['v'])
#     cyfit = fit(line,cypc,cdf['v'])
    
#     print(cat['name'][i],np.around(np.array((np.linalg.norm((cxfit[0][0],
#                                               cyfit[0][0])),oldcat['vel_grad_GMC<km/s/pc>'][i])),4))


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    
    cdfold = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(np.where(~msks[i],np.nan,c1s[i])),
            )).transpose(),
            columns=['l','b','v'])
    
#     print(cdf)
    
    cdf_noisy = cdf
     
    
    cdfnew = cdf_noisy[~(np.isnan(cdf_noisy['v']))]
    
    plt.scatter(cdfold['l'],cdfold['b'],c=cdfold['v'],cmap='RdGy')
    plt.scatter(cdfnew['l'],cdfnew['b'],c=cdfnew['v'],cmap='PiYG')
    plt.show()
    
    #gonna go with new way
    
#     cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
#     cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
#     cxfit = fit(line,cxpc,cdf['v'])
#     cyfit = fit(line,cypc,cdf['v'])
    
#     print(cat['name'][i],np.around(np.array((np.linalg.norm((cxfit[0][0],
#                                               cyfit[0][0])),oldcat['vel_grad_GMC<km/s/pc>'][i])),4))


from astropy.modeling import fitting,models


plinit = models.Planar2D(1.,1.,1)
fttr = fitting.LevMarLSQFitter(calc_uncertainties=True)


for i in [1,2,4,6,7,11,12,13,16,20,21]: #for some reason these are the only ones astropy is fitting

   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)


    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(np.where(~msks[i],np.nan,c1s[i])),
            )).transpose(),
            columns=['l','b','v'])

#     print(cdf)

    cdf_noisy = cdf


    cdf = cdf_noisy[~(np.isnan(cdf_noisy['v']))]


    #gonna go with new way

    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])


    a,b,vo = fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[0]

#     print(cat['name'][i],'\n',
#           np.around((cxfit[0][0],a),4),'\n',
#           np.around((cyfit[0][0],b),4),'\n',
#           np.around((np.nanmean(cdf['v']),vo),4),'\n\n\n')


    astrofit = fttr(plinit,xc*u.deg.to(u.rad)*cat['distance<pc>'][i],
                 yc*u.deg.to(u.rad)*cat['distance<pc>'][i],
                 c1s[i],weights=msks[i])

    print(cat['name'][i],'\n',
          np.around((
          np.linalg.norm((cxfit[0][0],cyfit[0][0])),
              np.linalg.norm((a,b)),
              np.linalg.norm((astrofit.slope_x.value,
                              astrofit.slope_y.value))),4),'\n')


#     print(cat['name'][i],np.around(np.array((np.linalg.norm((cxfit[0][0],
#                                               cyfit[0][0])),oldcat['vel_grad_GMC<km/s/pc>'][i])),4))


# # for when i come back
# currently trying out diff methods of fitting (ie whether to do x and y separately as linear eqs, or whether to do them as a plane, together)
# doing scipy optimize curve fit works, with the equation of a plane that i defined
# but it gives diff results than the lines separately
# 
# so now doing the astropy plane fitting to see which of the two it corroborates (or gives its own answer?)
# 
# currently leaning towards scipy optimize -- although .. it gives small results for california and perseus
# so maybe should double check taht somehow (eg subtracting the modelled plane from the first moment map and seeing which has smaller residuals


#okay so scipy and astropy are in agreement 100% of the time so gonna go with that
#also check how it does with hi vv


#change this to do hi, see if similar trends as co (whihc clouds astropy can fit or not are here too)


for i in range(22):
   
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)


    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
        np.ravel(h1s[i]),#checked and dont need to do the whole h2 isnt inf thing bc all the h2s are real,>0
            )).transpose(),
            columns=['l','b','v'])
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
    
    
    a,b,vo = fit(plane_rvld,np.array((hxpc,hypc)),hdf['v'])[0]
    
    astrofit = fttr(plinit,xh*u.deg.to(u.rad)*cat['distance<pc>'][i],
                 yh*u.deg.to(u.rad)*cat['distance<pc>'][i],
                 h1s[i])

    print(cat['name'][i],'\n',
          np.around((
          np.linalg.norm((hxfit[0][0],hyfit[0][0])),
              np.linalg.norm((a,b)),
              np.linalg.norm((astrofit.slope_x.value,
                              astrofit.slope_y.value))),4),'\n')

    
#okay so for hi they all give the same result .. makes me think that the co being sparse is causing an issue?
#check visually what works better vv


for i in [1,2,4,6,7,11,12,13,16,20,21]: #for some reason these are the only ones astropy is fitting

   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)


    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(np.where(~msks[i],np.nan,c1s[i])),
            )).transpose(),
            columns=['l','b','v'])

    cdf_noisy = cdf


    cdf = cdf_noisy[~(np.isnan(cdf_noisy['v']))]


    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]


    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])


    a,b,vo = fit(plane_rvld,np.array((cxpc,cypc)),cdf['v'])[0]


    astrofit = fttr(plinit,xc*u.deg.to(u.rad)*cat['distance<pc>'][i],
                 yc*u.deg.to(u.rad)*cat['distance<pc>'][i],
                 c1s[i],weights=msks[i])

    plt.title(cat['name'][i]+' x')

    plt.scatter(cxpc,cdf['v'],s=5)
    plt.plot(cxpc,(cxpc-np.nanmean(np.unique(cxpc)))*a+vo,color='red',label='plane')
    plt.plot(cxpc,cxfit[0][0]*cxpc+cxfit[0][1],color='orange',label='line')
    plt.show()
    
    plt.title(cat['name'][i]+' y')
    
    plt.scatter(cypc,cdf['v'],s=5)
    plt.plot(cypc,(cypc-np.nanmean(np.unique(cypc)))*b+vo,color='red',label='plane')
    plt.plot(cypc,cyfit[0][0]*cypc+cyfit[0][1],color='orange',label='line')
    plt.show()
    
#     plt.show()


#actually i think i will keep using the line for the co, and for hi it doesnt make a difference either way
#prob will also be good to report that (maybe in an appendix, and to show the fits)
#the plane fits end up being shallower (use a couple of the more dramatic examples to illustrate)
#eg california y, maddalena x,y , 


plt.scatter(cdf['l'],cdf['b'],c=(3.19428934-0.02*(cxpc-80.5)+0.08*(cypc+37.3)))
plt.colorbar()
#make this plot (showing the fit profile as one of the figures for the paper? to visually compare and or subtract)


np.indices(c1s[0].shape)[0]*-0.125+360


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)

    
    fakex = np.indices(c1s[i].shape)[1]*-0.125+360
    fakey = np.indices(c1s[i].shape)[0]*0.125-50
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(fakex),np.ravel(fakey),
            np.ravel(np.where(~msks[i],np.nan,c1s[i])),
            )).transpose(),
            columns=['l','b','fl','fb','v'])
    
    cdf_noisy = cdf


    cdf = cdf_noisy[~(np.isnan(cdf_noisy['v']))]

    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    fcxpc = cdf['fl']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    fcypc = cdf['fb']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    fcxfit = fit(line,fcxpc,cdf['v'])
    fcyfit = fit(line,fcypc,cdf['v'])


    print(round(cxfit[0][0],5)==round(fcxfit[0][0],5),'\n',
          round(cyfit[0][0],5)==round(fcyfit[0][0],5),
          '\n\n')
    
    #don't need the actual coordinates, just turn the pixel spacings into pc spacings


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

    
    if not full:
        return xslope,yslope
    else:
        return {'x_m_b_msig_bsig':[xslope, xintercept, xmsig, xbsig],
                'y_m_b_msig_bsig':[yslope, yintercept, ymsig, ybsig]}


dvdx_gmc,dvdy_gmc = np.array([fit_gradient(c1s[i],cat['distance<pc>'][i],mask=msks[i])
                   for i in range(22)]).transpose()


dvdx_gmc


dvdy_gmc


gmc_Omega = [np.linalg.norm([dvdx_gmc[i],dvdy_gmc[i]]) for i in range(22)]


dvdx_hi,dvdy_hi = np.array([fit_gradient(h1s[i],cat['distance<pc>'][i])
            for i in range(22)]).transpose()


dvdx_hi,dvdy_hi


hi_Omega = [np.linalg.norm([dvdx_hi[i],dvdy_hi[i]]) for i in range(22)]


def dxdxy_2_theta(dx,dy):
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


gmc_Theta = np.array([dxdxy_2_theta(dvdx_gmc[i],dvdy_gmc[i]) for i in range(22)])
hi_Theta = np.array([dxdxy_2_theta(dvdx_hi[i],dvdy_hi[i]) for i in range(22)])


gmc_rot = gmc_Theta-90
hi_rot = hi_Theta-90


for i in range(22):
    if gmc_rot[i] > 180:
        gmc_rot[i] -= 360
        
    if gmc_rot[i] < -180:
        gmc_rot[i] += 360


for i in range(22):
    if hi_rot[i] > 180:
        hi_rot[i] -= 360
        
    if hi_rot[i] < -180:
        hi_rot[i] += 360


cat['vel_grad_GMC<km/s/pc>'] = gmc_Omega
cat['vel_grad_HI<km/s/pc>'] = hi_Omega

cat['rot_ax_GMC<deg>'] = gmc_rot
cat['rot_ax_HI<deg>'] = hi_rot


gmc_area = np.array([(np.nansum(
    (c0s[i]*msks[i])>0.)*\
    (0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]*u.pc)**2).to(u.pc**2).value

 for i in range(22)])


cat['j_GMC<cm2/s>'] = ((cat['radius_GMC<pc>']*u.pc)**2 * (cat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)).cgs
cat['j_HI<cm2/s>'] = ((cat['radius_HI<pc>']*u.pc)**2 * (cat['vel_grad_HI<km/s/pc>']*u.km/u.s/u.pc)).cgs


cat['J_GMC<erg.s>'] = (cat['j_GMC<cm2/s>']*cat['mass_GMC<Msun>']*u.solMass).cgs
cat['J_HI<erg.s>'] = (cat['j_HI<cm2/s>']*cat['mass_HI<Msun>']*u.solMass).cgs


cat.write('sighhplease.fits',overwrite=True)


fns.make_fits_from_list(c0s,'h2_0s.fits')
fns.make_fits_from_list(c1s,'h2_1s.fits')
fns.make_fits_from_list(c2s,'h2_2s.fits')


fns.make_fits_from_list(h0s,'hi_0s.fits')
fns.make_fits_from_list(h1s,'hi_1s.fits')
fns.make_fits_from_list(h2s,'hi_2s.fits')
