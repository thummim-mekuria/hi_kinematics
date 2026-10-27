#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import fitting,models
fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.)

import numpy as np
import pandas as pd

import datacube as dc
import fun as fns
import interpolate as ip


c0s = fits.open('h2_mom0.fits')
h0s = fits.open('hi_mom0.fits')

c1s = fits.open('h2_mom1.fits')
h1s = fits.open('hi_mom1.fits')

c2s = fits.open('h2_mom2.fits')
h2s = fits.open('hi_mom2.fits')


cat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


names = ['Aquila South', 'California', 'Camelopardalis', 'Chamaeleon',
       'Canis Major OB1', 'Corona Australis', 'Crossbones', 'Gemini OB1',
       'Hercules', 'Lacerta', 'Lupus', 'Maddalena', 'Monoceros OB1', 'Monoceros R2',
       r'Ophiuchus', 'Orion A', 'Orion B', r'$\lambda$-Orionis', 'Perseus', 'Polaris',
       'Rosette', 'Taurus']


lbls = ['AqS', 'Cal', 'Cam', 'Chm',
       'CMa', 'CrA', 'Xbn', 'Gem',
       'Her', 'Lac', 'Lup', 'Mad', 'MOB', 'MR2',
       r'Oph', 'OrA', 'OrB', r'$\lambda$Or', 'Per', 'Pol',
       'Ros', 'Tau']


hi_dvds = np.zeros(44).reshape(22,2)
co_dvds = np.zeros(44).reshape(22,2)


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)
    
#     xc,yc = np.meshgrid(np.arange(cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],-0.125),
#                         np.arange(cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],0.125))
                        
#     xh,yh = np.meshgrid(np.arange(cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],-0.125),
#                         np.arange(cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],0.125))
    
#     print(xh.shape,yh.shape,h0s[i].shape,h1s[i].shape,h2s[i].shape)
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0s[i].data),
            np.ravel(h1s[i].data),
            np.ravel(h2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    
    hdf = hdf_noisy[
        (hdf_noisy['s']>0) & \
                    ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] 
    #make sure that the linewidths are a valid number 
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    #converting to pc now makes life wayyy easier
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]

    
    dvdxh = fttr(init,hxpc,hdf['v']) 
    dvdyh = fttr(init,hypc,hdf['v'])
    
    hi_dvds[i] = [dvdxh.slope.value,dvdyh.slope.value]
    
    dvdxc = fttr(init,cxpc,cdf['v'])
    dvdyc = fttr(init,cypc,cdf['v'])
    
    co_dvds[i] = [dvdxc.slope.value,dvdyc.slope.value]
      
    
    omh = np.sqrt(dvdxh.slope**2 + dvdyh.slope**2)
    omc = np.sqrt(dvdxc.slope**2 + dvdyc.slope**2)

    thh = np.arctan(dvdyh.slope/dvdxh.slope)*u.rad.to(u.deg)
    thc = np.arctan(dvdyc.slope/dvdxc.slope)*u.rad.to(u.deg)
    
#     print(new['name'][i],dvdyc.slope,dvdxc.slope)
    #correcting for arctan only returning |theta|<90
    if dvdxh.slope > 0 : #first & fourth quadrants
        thh = thh
    elif dvdxh.slope < 0: #second quadrant & third
        thh = thh + 180
        
    if dvdxc.slope > 0 : #first & fourth quadrants
        thc = thc  
    elif dvdxc.slope < 0: #second quadrant & third
        thc = thc - 180

#     print(len(cdf),len(cdf_noisy))#,len(hdf),len(hdf_noisy))
    
    
#     plt.scatter(cdf['l'],cdf['b'],c=cdf['v'])
#     plt.show()
#     plt.scatter(cdf_noisy['l'],cdf_noisy['b'],c=cdf_noisy['v'])
#     plt.show()
    
    
#     plt.scatter(hdf['l'],hdf['b'],c=hdf['v'])
#     plt.show()
#     plt.scatter(hdf_noisy['l'],hdf_noisy['b'],c=hdf_noisy['v'])
#     plt.show()
    
#     plt.scatter(cxpc,cdf['v'],s=1)
#     plt.plot(cxpc,cxpc*dvdxc.slope + dvdxc.intercept,color='red')
#     plt.show()
    
    plt.scatter(hxpc,hdf['v'],s=1)
# #     plt.scatter(np.ravel(xh)*u.deg.to(u.rad)*cat['distance<pc>'][i],np.ravel(h1s[i].data),s=1.5)
    plt.title(cat['name'][i])
    plt.plot(hxpc,hxpc*dvdxh.slope + dvdxh.intercept,color='red')
    plt.show()
    
    
#     plt.imshow(h1,origin='lower',cmap='RdGy')
#     plt.show()
#     plt.imshow(c1,origin='lower',cmap='PiYG')
#     plt.show()
    
    
#     print(fcat['name'][i],dvdxh.slope.value,dvdyh.slope.value,
#           dvdxc.slope.value,dvdyc.slope.value,'\n\n')
    
#     dvdxh.slope.value * (0.125*u.deg.to(u.rad)*fcat['distance<pc>'] #km/s/deg * deg/pc
    
#     print(fcat['name'][i],hplane.slope_x,hplane.slope_y,
#           cplane.slope_x,cplane.slope_y,'\n\n\n\n\n')

    
#     plt.scatter(np.ravel(xh),np.ravel(h1s[i].data),s=1)
#     plt.scatter(np.ravel(xc),np.ravel(c1s[i].data),s=3)
#     plt.show()

    
cat['vel_grad_HI<km/s/pc>'] == np.linalg.norm(hi_dvds,axis=1)


cat['vel_grad_GMC<km/s/pc>'] == np.linalg.norm(co_dvds,axis=1)


help(plt.quiver)


hi4pi,dht = fits.open('both_same_grid_720d.fits')


h1 = fns.get_mom1(hi4pi.data,np.arange(155)*1.3-100)


dh1ds = np.gradient(h1)


plt.contourf(h1,extent=[360,-360,-50,50],extend='both',
             levels=np.arange(-30,30,7),cmap='Spectral')
plt.xlim(200,-200)


x,y = np.meshgrid(np.arange(5760)*-0.125+360,np.arange(800)*0.125-50)


plt.figure(figsize=(20,10))

plt.quiver(x,y,dh1ds[1]/np.linalg.norm(dh1ds,axis=0),
           dh1ds[0]/np.linalg.norm(dh1ds,axis=0))

plt.xlim(180,-180)

# figure out ho to make arrows smaller, and see if the hi arrows match up


plt.figure(figsize=(20,10))

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
          -co_dvds.transpose()[0], #so it looks right for l going +->-
           co_dvds.transpose()[1],
           alpha=0.7,
          pivot='middle',width=5e-3,minshaft=2)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
          -hi_dvds.transpose()[0], #so it looks right for l going +->-
           hi_dvds.transpose()[1],
           alpha=0.7,color='red',
          pivot='middle',width=5e-3,minshaft=2)

plt.xlim(200,-200)


plt.figure(figsize=(20,10))

plt.contourf(h1,extent=[360,-360,-50,50],extend='both',
             levels=np.arange(-30,30,7),cmap='Spectral')


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
          -co_dvds.transpose()[0], #so it looks right for l going +->-
           co_dvds.transpose()[1],
           alpha=0.7,
          pivot='middle',width=5e-3,minshaft=2)

# plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
#           -hi_dvds.transpose()[0],
#            hi_dvds.transpose()[1],
#            color='orange',alpha=0.5,pivot='middle',width=5e-3,minshaft=2)

plt.xlim(200,-200)


plt.figure(figsize=(20,10))

plt.contourf(h1,extent=[360,-360,-50,50],extend='both',
             levels=np.arange(-30,30,7),cmap='Spectral')


# plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
#           -co_dvds.transpose()[0], #so it looks right for l going +->-
#            co_dvds.transpose()[1],
#            alpha=0.7,
#           pivot='middle',width=5e-3,minshaft=2)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
          -hi_dvds.transpose()[0],
           hi_dvds.transpose()[1],
           alpha=0.7,pivot='middle',width=5e-3,minshaft=2)

plt.xlim(200,-200)


np.array((cat['name'],cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>'])).transpose()


hi_dvds[0] = np.array([1,2])


[dvdxc.slope,dvdyc.slope]


for i in range(22):
    l0h,l1h = cat['lon_HI<deg>'][i] #np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = cat['lat_HI<deg>'][i] #np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = cat['lon_GMC<deg>'][i] #np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = cat['lat_GMC<deg>'][i] #np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i].data)
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i].data)

    
    lc,bc = np.nanmean(cat['lon_GMC<deg>'][i]),np.nanmean(cat['lat_GMC<deg>'][i])
    
    hthet = cat['rot_ax_HI<deg>'][i]*u.deg
    cthet = cat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
    
#     plt.contourf(bperp.reshape(h1s[i].shape).value,
#                  extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
#     plt.colorbar()
#     plt.scatter(lons,rot_ax_line)
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.show()
    
    
    hv_model = np.nanmean(hmom1) + \
                hi_dvds[i][0]*(hlons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                hi_dvds[i][1]*(hlats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    cv_model = np.nanmean(cmom1) + \
                co_dvds[i][0]*(clons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                co_dvds[i][1]*(clats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    
    hbperp_pc = hbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cbperp_pc = cbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    fig,axs = plt.subplots(1,2,figsize=(10,5))
    fig.suptitle(cat['name'][i])
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(c0s[i].data),
            cmom1,
            cbperp_pc,
            cv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
    
    axs[0].scatter(cdf['r'],cdf['v'],s=500/cbperp.size,color='black')
    
    axs[0].plot(cdf['r'],cdf['vmod'],color='red')
    
    #not sure if these are the right error bars ... 
#     axs[0].plot(cdf['r'],cdf['vmod']-np.nanstd(cdf['vmod']),color='pink') 
#     axs[0].plot(cdf['r'],cdf['vmod']+np.nanstd(cdf['vmod']),color='pink')
    
    
    axs[1].scatter(hbperp_pc,hmom1,s=500/hbperp.size,color='black')
    
    axs[1].plot(hbperp_pc,hv_model,color='red')
    
    #not sure if these are the right error bars ... 
#     axs[1].plot(hbperp_pc,hv_model-np.nanstd(hv_model),color='pink') 
#     axs[1].plot(hbperp_pc,hv_model+np.nanstd(hv_model),color='pink')
    
    
#     axs[1].plot(hbperp_pc,hv_model-np.nanstd(hmom1),color='orange') 
#     axs[1].plot(hbperp_pc,hv_model+np.nanstd(hmom1),color='orange')
    
    
    for ax in axs.flat:
        ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
        
        
    axs[0].set_ylabel(r'velocity [$km\ s^{-1}$]')
    
#     axs[0].set_ylim(np.nanmean(cdf['vmod'])-5*np.nanstd(cdf['vmod']),
#                     np.nanmean(cdf['vmod'])+5*np.nanstd(cdf['vmod']))
    
#     axs[1].set_ylim(np.nanmean(hv_model)-5*np.nanstd(hv_model),
#                     np.nanmean(hv_model)+5*np.nanstd(hv_model))
    
    
#     plt.savefig(fcat['name'][i]+'vel_grad_fit',bbox_inches='tight')
    
#     print(np.nanstd(hv_model),np.nanstd(hmom1))
    plt.show()
    
    
from scipy.optimize import curve_fit as fit


dvdxh


help(fit)


def line(x,m,b):
    return m*x+b


fit(line,hxpc,hdf['v']) 


np.sqrt(np.diag(fit(line,hxpc,hdf['v'])[1]))


fttr(init,hxpc,hdf['v'])


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)
    
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0s[i].data),
            np.ravel(h1s[i].data),
            np.ravel(h2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    
    hdf = hdf_noisy[
        (hdf_noisy['s']>0) & \
                    ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] 
    #make sure that the linewidths are a valid number 
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    #converting to pc now makes life wayyy easier
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]

    
    dvdxh = fttr(init,hxpc,hdf['v']) 
    dvdyh = fttr(init,hypc,hdf['v'])
    
    #okay so the fftr & fiting a line with curve fit are doing the same thing
    #which means i dont have to recalculate stuff from before, i can judt get the errors from scipy
    
#     print(np.around(dvdxh.parameters,4) == np.around(fit(line,hxpc,hdf['v'])[0],4))
#     print(np.around(dvdyh.parameters,4) == np.around(fit(line,hypc,hdf['v'])[0],4))
#     print('\n')
    
    dvdxc = fttr(init,cxpc,cdf['v'])
    dvdyc = fttr(init,cypc,cdf['v'])
    
    print(np.around(dvdxc.parameters,5) == np.around(fit(line,cxpc,cdf['v'])[0],5))
    print(np.around(dvdyc.parameters,5) == np.around(fit(line,cypc,cdf['v'])[0],5))
    print('\n')      
 

for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)
    
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0s[i].data),
            np.ravel(h1s[i].data),
            np.ravel(h2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    
    hdf = hdf_noisy[
        (hdf_noisy['s']>0) & \
                    ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] 
    #make sure that the linewidths are a valid number 
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    #converting to pc now makes life wayyy easier
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    

    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
    
    hxfit_err = np.sqrt(np.diag(hxfit[1]))
    hyfit_err = np.sqrt(np.diag(hyfit[1]))
    
    
    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    cxfit_err = np.sqrt(np.diag(cxfit[1]))
    cyfit_err = np.sqrt(np.diag(cyfit[1]))
    
    
    print(hxfit_err,hyfit_err,cxfit_err,cyfit_err)
    
    
# hmm so trying ti decide degreeof linearity based on the m/sigma for the r vs v doesnt seem to work..? or at least
# the sigma_m look too small to me
# maybe look at the sigma_b -> aka whether the spread-ed ness of the data is actually being encoded/saved in how much vertical offset there is, rather than how linear
# alsoo look some more in scipy to see if there si a function for calculating goodness of fit or if theres a way to get the errorbars to say something more intuitive
# but! might be useful to show the r vs v plots -> maybe even for a hercules / other ones where there is a strong gradient but its doing something far from linear
# bc i feel like some of this is being influenced (in terms of how big the slope of the fit line is) on the absolute value / difference in velocity between the two ends of r, rather than on the shape of the distribution of the velocities that are inside


fttr(init,np.log10(cat['radius_GMC<pc>']),
     np.log10((cat['ang_mom_GMC<erg.s>']*u.erg*u.s/cat['mass_GMC<Msun>']/u.solMass).cgs.value))


jR_fit = fit(line,np.log10(cat['radius_GMC<pc>']),
   np.log10((cat['ang_mom_GMC<erg.s>']*u.erg*u.s/cat['mass_GMC<Msun>']/u.solMass).cgs.value))


jR_fit = fit(line,np.log10(cat['radius_HI<pc>']),
   np.log10((cat['ang_mom_HI<erg.s>']*u.erg*u.s/cat['mass_HI<Msun>']/u.solMass).cgs.value))


jR_err = np.sqrt(np.diag(jR_fit[1]))


jR_fit[0][0] - jR_err[0],jR_fit[0][0] + jR_err[0]


plt.scatter(cat['radius_GMC<pc>'],(cat['ang_mom_GMC<erg.s>']*u.erg*u.s/cat['mass_GMC<Msun>']/u.solMass).cgs)
plt.scatter(cat['radius_HI<pc>'],(cat['ang_mom_HI<erg.s>']*u.erg*u.s/cat['mass_HI<Msun>']/u.solMass).cgs)


# plt.plot(np.arange(3e0,1e2,1e1),10**(22.64)*np.arange(3e0,1e2,1e1)**(1.81))

# plt.plot(np.arange(3e0,1e2,1e1),10**(jR_fit[0][1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]))

# plt.fill_between(np.arange(3e0,1e2,1e1),
#                  10**(jR_fit[0][1]-jR_err[1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]),
#                 10**(jR_fit[0][1]+jR_err[1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]),
#                 alpha=0.3)


plt.fill_between(np.arange(3e0,1e2,1e1),10**(22.64)*np.arange(3e0,1e2,1e1)**(1.55),
                10**(22.64)*np.arange(3e0,1e2,1e1)**(2.07),alpha=0.1,color='grey')


plt.fill_between(np.arange(3e0,1e2,1e1),
                 10**(jR_fit[0][1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]-jR_err[0]),
                10**(jR_fit[0][1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]+jR_err[0]),
                alpha=0.1,color='grey')

# plt.fill_between(np.arange(3e0,1e2,1e1),
#                  10**(jR_fit[0][1])*np.arange(3e0,1e2,1e1)**(1.2),
#                 10**(jR_fit[0][1])*np.arange(3e0,1e2,1e1)**(1.8),
#                 alpha=0.1,color='grey')


# plt.fill_between(np.arange(3e0,1e2,1e1),
#                  10**(jR_fit[0][1]-jR_err[1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]-jR_err[0]),
#                 10**(jR_fit[0][1]+jR_err[1])*np.arange(3e0,1e2,1e1)**(jR_fit[0][0]+jR_err[0]),
#                 alpha=0.3)


plt.loglog()

plt.show()


for i in range(22):
    l0h,l1h = cat['lon_HI<deg>'][i] #np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = cat['lat_HI<deg>'][i] #np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = cat['lon_GMC<deg>'][i] #np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = cat['lat_GMC<deg>'][i] #np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i].data)
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i].data)

    
    lc,bc = np.nanmean(cat['lon_GMC<deg>'][i]),np.nanmean(cat['lat_GMC<deg>'][i])
    
    hthet = cat['rot_ax_HI<deg>'][i]*u.deg
    cthet = cat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
    
#     plt.contourf(bperp.reshape(h1s[i].shape).value,
#                  extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
#     plt.colorbar()
#     plt.scatter(lons,rot_ax_line)
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.show()
    
    
    hv_model = np.nanmean(hmom1) + \
                hi_dvds[i][0]*(hlons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                hi_dvds[i][1]*(hlats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    cv_model = np.nanmean(cmom1) + \
                co_dvds[i][0]*(clons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                co_dvds[i][1]*(clats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    
    hbperp_pc = hbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cbperp_pc = cbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    fig,axs = plt.subplots(1,2,figsize=(10,5))
    fig.suptitle(cat['name'][i])
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(c0s[i].data),
            cmom1,
            cbperp_pc,
            cv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
     
    axs[0].scatter(cdf['r'],cdf['v'],s=500/cbperp.size,color='black')
    
    axs[0].plot(cdf['r'],cdf['vmod'],color='red')
#     axs[0].fill_between(cdf['r'],cdf['r']*(cm[0]+cs) + cm[1],
#                         cdf['r']*(cm[0]-cs) + cm[1],color='orange',alpha=0.5)
    
    
    axs[1].scatter(hbperp_pc,hmom1,s=500/hbperp.size,color='black')
    
    axs[1].plot(hbperp_pc,hv_model,color='red')
    
    
    for ax in axs.flat:
        ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
        
        
    axs[0].set_ylabel(r'velocity [$km\ s^{-1}$]')
    

    newcdf = cdf.groupby(pd.cut(cdf['r'],
                        bins=np.arange(np.nanmin(cdf['r'])-2,
                        np.nanmax(cdf['r'])+2,
                        0.05*cat['radius_GMC<pc>'][i]))).mean().dropna()
    
    
    hp = fit(line,hbperp_pc,hmom1)[0]
    hs = np.sqrt(np.diag(fit(line,hbperp_pc,hmom1)[1]))
    
    cp = fit(line,cdf['r'],cdf['v'])[0]
    cs = np.sqrt(np.diag(fit(line,cdf['r'],cdf['v'])[1]))
    
    cpn = fit(line,newcdf['r'],newcdf['v'])[0]
    csn = np.sqrt(np.diag(fit(line,newcdf['r'],newcdf['v'])[1]))
    
#     axs[0].plot(newcdf['r'],newcdf['r']*cpn[0]+cpn[1],color='purple')

    axs[0].plot(newcdf['r'],newcdf['v'],color='purple')
#     print(cm/cs,hm/hs,names[i])


#     print(cm/cs,cmn/csn)
#     print(np.sqrt(np.diag(fit(line,hbperp_pc,hmom1)[1])))
    
#     print(round(fit(line,hbperp_pc,hmom1)[0][0],4)==round(fit(line,hbperp_pc,hv_model)[0][0],4))
#     print(fit(line,cdf['r'],cdf['v'])[0][0])
#     ,fit(line,cdf['r'],cdf['vmod'])[0][0],\
#          cat['vel_grad_GMC<km/s/pc>'][i])


    plt.show()


h1s[i].data.shape


for i in range(22):
    l0h,l1h = cat['lon_HI<deg>'][i] #np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = cat['lat_HI<deg>'][i] #np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = cat['lon_GMC<deg>'][i] #np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = cat['lat_GMC<deg>'][i] #np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i].data)
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i].data)

    
    lc,bc = np.nanmean(cat['lon_GMC<deg>'][i]),np.nanmean(cat['lat_GMC<deg>'][i])
    
    hthet = cat['rot_ax_HI<deg>'][i]*u.deg
    cthet = cat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
    
#     plt.contourf(bperp.reshape(h1s[i].shape).value,
#                  extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
#     plt.colorbar()
#     plt.scatter(lons,rot_ax_line)
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.show()
    
    
    hv_model = np.nanmean(hmom1) + \
                hi_dvds[i][0]*(hlons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                hi_dvds[i][1]*(hlats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 

    fig,axs = plt.subplots(1,2,figsize=(10,5))
    fig.suptitle(cat['name'][i])
    
    axs[0].contourf(hv_model.reshape(h1s[i].data.shape))
    axs[1].contourf(h1s[i].data)

    plt.show()
    
    plt.hist(np.ravel(hv_model-hmom1),bins=np.arange(-5,5,0.1))
    plt.show()


171*171


for i in range(22):
    l0h,l1h = cat['lon_HI<deg>'][i] #np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = cat['lat_HI<deg>'][i] #np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = cat['lon_GMC<deg>'][i] #np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = cat['lat_GMC<deg>'][i] #np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i].data)
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i].data)

    
    lc,bc = np.nanmean(cat['lon_GMC<deg>'][i]),np.nanmean(cat['lat_GMC<deg>'][i])
    
    hthet = cat['rot_ax_HI<deg>'][i]*u.deg
    cthet = cat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
    
#     plt.contourf(bperp.reshape(h1s[i].shape).value,
#                  extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
#     plt.colorbar()
#     plt.scatter(lons,rot_ax_line)
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.show()
    
    
    hv_model = np.nanmean(hmom1) + \
                hi_dvds[i][0]*(hlons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                hi_dvds[i][1]*(hlats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    cv_model = np.nanmean(cmom1) + \
                co_dvds[i][0]*(clons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                co_dvds[i][1]*(clats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    
    hbperp_pc = hbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cbperp_pc = cbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    fig,axs = plt.subplots(1,2,figsize=(10,5))
    fig.suptitle(cat['name'][i])
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(c0s[i].data),
            cmom1,
            cbperp_pc,
            cv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
     
    axs[0].scatter(cdf['r'],cdf['v'],s=500/cbperp.size,color='black')
    
    axs[0].plot(cdf['r'],cdf['vmod'],color='red')
#     axs[0].fill_between(cdf['r'],cdf['r']*(cm[0]+cs) + cm[1],
#                         cdf['r']*(cm[0]-cs) + cm[1],color='orange',alpha=0.5)
    
    
    axs[1].scatter(hbperp_pc,hmom1,s=500/hbperp.size,color='black')
    
    axs[1].plot(hbperp_pc,hv_model,color='red')
    
    
    for ax in axs.flat:
        ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
        
        
    axs[0].set_ylabel(r'velocity [$km\ s^{-1}$]')
    

    newcdf = cdf.groupby(pd.cut(cdf['r'],
                        bins=np.arange(np.nanmin(cdf['r'])-2,
                        np.nanmax(cdf['r'])+2,
                        0.05*cat['radius_GMC<pc>'][i]))).mean().dropna()
    
    
    hp = fit(line,hbperp_pc,hmom1)[0]
    hs = np.sqrt(np.diag(fit(line,hbperp_pc,hmom1)[1]))
    
    cp = fit(line,cdf['r'],cdf['v'])[0]
    cs = np.sqrt(np.diag(fit(line,cdf['r'],cdf['v'])[1]))
    
    cpn = fit(line,newcdf['r'],newcdf['v'])[0]
    csn = np.sqrt(np.diag(fit(line,newcdf['r'],newcdf['v'])[1]))
    
#     axs[0].plot(newcdf['r'],newcdf['r']*cpn[0]+cpn[1],color='purple')

    axs[0].plot(newcdf['r'],newcdf['v'],color='purple')
#     print(cm/cs,hm/hs,names[i])


#     print(cm/cs,cmn/csn)
#     print(np.sqrt(np.diag(fit(line,hbperp_pc,hmom1)[1])))
    
#     print(round(fit(line,hbperp_pc,hmom1)[0][0],4)==round(fit(line,hbperp_pc,hv_model)[0][0],4))
#     print(fit(line,cdf['r'],cdf['v'])[0][0])
#     ,fit(line,cdf['r'],cdf['vmod'])[0][0],\
#          cat['vel_grad_GMC<km/s/pc>'][i])


    plt.show()


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


for i in range(22):
    l0h,l1h = cat['lon_HI<deg>'][i] #np.array((loh[i],lih[i]))*-0.125 + 360
    b0h,b1h = cat['lat_HI<deg>'][i] #np.array((boh[i],bih[i]))*0.125 -50

    
    l0c,l1c = cat['lon_GMC<deg>'][i] #np.array((loc[i],lic[i]))*-0.125 + 360
    b0c,b1c = cat['lat_GMC<deg>'][i] #np.array((boc[i],bic[i]))*0.125 -50
    
    
    hlons = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0])
    hlats = np.ravel(np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1])
    hmom1 = np.ravel(h1s[i].data)
    
    
    clons = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0])
    clats = np.ravel(np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1])
    cmom1 = np.ravel(c1s[i].data)

    
    lc,bc = np.nanmean(cat['lon_GMC<deg>'][i]),np.nanmean(cat['lat_GMC<deg>'][i])
    
    hthet = cat['rot_ax_HI<deg>'][i]*u.deg
    cthet = cat['rot_ax_GMC<deg>'][i]*u.deg

    hrot_ax = bc+(hlons-lc)*np.tan(hthet)
    crot_ax = bc+(clons-lc)*np.tan(cthet)
    
    
    hbperp = np.cos(hthet)*(np.tan(hthet)*(hlons-lc)-hlats+bc) #point<->line distance
    cbperp = np.cos(cthet)*(np.tan(cthet)*(clons-lc)-clats+bc) #point<->line distance
    
    
#     plt.contourf(bperp.reshape(h1s[i].shape).value,
#                  extent=[l0h,l1h,b0h,b1h],cmap='RdGy')
#     plt.colorbar()
#     plt.scatter(lons,rot_ax_line)
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.show()
    
    
    hv_model = np.nanmean(hmom1) + \
                hi_dvds[i][0]*(hlons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                hi_dvds[i][1]*(hlats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    cv_model = np.nanmean(cmom1) + \
                co_dvds[i][0]*(clons-lc)*u.deg.to(u.rad)*cat['distance<pc>'][i] + \
                co_dvds[i][1]*(clats-bc)*u.deg.to(u.rad)*cat['distance<pc>'][i] 
    
    
    hbperp_pc = hbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cbperp_pc = cbperp*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(c0s[i].data),
            cmom1,
            cbperp_pc,
            cv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
    
    newcdf = cdf.groupby(pd.cut(cdf['r'],
                        bins=np.arange(np.nanmin(cdf['r'])-2,
                        np.nanmax(cdf['r'])+2,
                        0.05*cat['radius_GMC<pc>'][i]))).mean().dropna()
    
    
    hdf = pd.DataFrame(np.array((
            np.ravel(h0s[i].data),
            hmom1,
            hbperp_pc,
            hv_model
            )).transpose(),
            columns=['w','v','r','vmod'])
    
     
    newhdf = hdf.groupby(pd.cut(hdf['r'],
                        bins=np.arange(np.nanmin(hdf['r'])-2,
                        np.nanmax(hdf['r'])+2,
                        0.05*cat['radius_HI<pc>'][i]))).mean().dropna()
    
    
    #####-------------------#####-------------------#####-------------------#####
    
    fig,axs = plt.subplots(1,4,sharey=True,dpi=200,figsize=(18,4))
    fig.suptitle(names[i])
    
    axs[0].scatter(cdf['r'],cdf['v'],s=400/len(cdf),color='grey')
    axs[1].scatter(hdf['r'],hdf['v'],s=800/len(hdf),color='grey')
    
    
    for ax in axs.flat[:-1]:
        ax.set_xlabel(r'offset from axis [pc]')
        
        
    axs[0].set_ylabel(r'$v\mathrm{_{_{LOS}}\ \ [km\ s^{-1}}$]')
    

    axs[0].plot(newcdf['r'],newcdf['v'],color='black',linewidth=0.9)
    axs[1].plot(newhdf['r'],newhdf['v'],color='black',linewidth=6,alpha=0.3)
    
    axs[0].plot(cdf['r'],cdf['vmod'],color='red',linestyle='--')
    axs[1].plot(hdf['r'],hdf['vmod'],color='red',linestyle='--',linewidth=6,alpha=0.3)

    
    plt.subplots_adjust(wspace=0.)
#     axs[1].set_yticks([])
    axs[0].set_xlim(np.nanmin(cdf['r'])-0.1*cat['radius_GMC<pc>'][i],
                    np.nanmax(cdf['r'])+0.1*cat['radius_GMC<pc>'][i])
    
    axs[1].set_xlim(np.nanmin(hdf['r'])-0.1*cat['radius_HI<pc>'][i],
                    np.nanmax(hdf['r'])+0.1*cat['radius_HI<pc>'][i])
    
    
#     axs[0].set_title('$\mathrm{H_2}$')
#     axs[1].set_title(r'$\mathrm{HI}$')
    
#     axs[1].set_xticks(np.around([np.nanmin(hdf['r']),0.5*np.nanmin(hdf['r']),
#                       0,-np.nanmin(hdf['r']),-0.5*np.nanmin(hdf['r'])]))
    
#     axs[0].set_xticks(np.around([np.nanmax(cdf['r']),0.5*np.nanmax(cdf['r']),
#                       0,-np.nanmax(cdf['r']),-0.5*np.nanmax(cdf['r'])]))

    if i in [15,16]:
        axs[0].set_ylim(2,17)
        
        
    axs[2].plot(newcdf['r'],newcdf['v'],color='black',linewidth=0.9)
    axs[2].plot(newhdf['r'],newhdf['v'],color='black',linewidth=6,alpha=0.3)
    
    axs[2].plot(cdf['r'],cdf['vmod'],color='red',linestyle='--')
    axs[2].plot(hdf['r'],hdf['vmod'],color='red',linestyle='--',linewidth=6,alpha=0.3)
    
    
    vbins = np.arange(min(np.nanmin(hdf['v']),np.nanmin(cdf['v']))-2,
             max(np.nanmax(hdf['v']),np.nanmax(cdf['v']))+2,0.65/2)
    
    
    axs[3].hist(cdf['v'],bins=vbins,color='black',
                orientation='horizontal',density=True,
               histtype='step')
    
    axs[3].hist(hdf['v'],bins=vbins,alpha=0.3,linewidth=6,
                orientation='horizontal',density=True,
               histtype='step',color='black') 
    
    axs[3].set_xticks(np.around(np.linspace(0,np.nanmax(
        np.histogram(hdf['v'],bins=vbins,density=True)[0]),3),1))
    axs[3].set_xlabel('$\mathrm{P\ (v_{_{LOS}})}$')
    
    plt.savefig('rv_mom1'+lbls[i],bbox_inches='tight')
    plt.show()


np.nanmax(np.histogram(hdf['v'],bins=vbins,density=True)[0])


def plane(xsys,vo,a,b):
    xs,ys = xsys
    
    xo = np.nanmean(xs)
    yo = np.nanmean(ys)
    
    return vo + a*(xs-xo) + b*(ys-yo)


r = cat['bih']-cat['boh'] #technically 2r


r3s = np.array((cat['boh'],cat['bih'],cat['loh'],cat['lih'])) + np.array([-r,r,-r,r])


his_3r = [hi4pi.data[cat['vo'][i]:cat['vi'][i],r3s[0][i]:r3s[1][i],r3s[2][i]:r3s[3][i]] for i in range(22)]


cos_3r = [dht.data[cat['vo'][i]:cat['vi'][i],r3s[0][i]:r3s[1][i],r3s[2][i]:r3s[3][i]] for i in range(22)]


big_coords = ((r3s*np.array([0.125,0.125,-0.125,-0.125])[:,np.newaxis])+
              np.array([-50,-50,360,360])[:,np.newaxis]).transpose()


h0s_3r = [np.nansum(his_3r[i],axis=0)*1.3*1.82e18 for i in range(22)]


c0s_3r = [np.nansum(cos_3r[i],axis=0)*1.3*1.97e20 for i in range(22)]


[(
  
plt.contourf(np.log10(h0s_3r[i].data),levels=np.arange(20,22,0.1),
            extent=[big_coords[i][2],big_coords[i][3],
                    big_coords[i][0],big_coords[i][1]],cmap='YlGnBu_r'),
    
plt.contour(np.log10(c0s_3r[i].data),levels=np.arange(20.2,21,0.3),
        extent=[big_coords[i][2],big_coords[i][3],
                big_coords[i][0],big_coords[i][1]],
            colors='white',linewidths=0.4),


# plt.contour(np.log10(h0s[i].data),levels=np.arange(20,22,0.1),
#             extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
#                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
# plt.contourf(np.log10(c0s[i].data),levels=np.arange(20,22,0.25),
#             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                    cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
plt.contour(np.log10(c0s[i].data),levels=np.arange(20.2,21,0.3),linewidths=0.4,
            colors='black',
            extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                   cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
# plt.contour(np.log10(c0s_3r[i].data),levels=[20.5,21],
#         extent=[big_coords[i][2],big_coords[i][3],
#                 big_coords[i][0],big_coords[i][1]],colors='red'),


# plt.xlim(cat['lon_HI<deg>'][i]+np.array([1,-1])),
# plt.ylim(cat['lat_HI<deg>'][i]+np.array([-1,1])),
    
    
plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           colors='black',linewidths=0.6),
plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           colors='black',linewidths=0.6),

plt.show()

)for i in range(22)]


# spend time in the heres how we select the clouds part talking about how we are making a choice to force stuff onto a rectangular grid, esp in how we decide the bounds, and its tricky to decide which emission features go with which cloud, esp when theres two things whose density contours overlap,but we know theyre separate things...


from astrodendro import Dendrogram


dend21 = Dendrogram.compute(cos_3r[21],min_value=0.2,min_delta=0.2,min_npix=5e2)


for i in dend21[0].children:
    print(i)


for i in [6,90,84]:
    plt.contourf(np.nansum(dend21[i].get_mask(),axis=2),levels=np.arange(1,9))
    plt.title(i)
#     plt.show()


for i in [6,90,84]:
    plt.contourf(np.nansum(dend21[i].get_mask(),axis=1),levels=np.arange(1,9))
    plt.title(i)
#     plt.show()


for i in [6,90,84]:
    plt.contourf(np.nansum(dend21[i].get_mask(),axis=0),levels=np.arange(1,9))
    plt.title(i)
#     plt.show()


for t in dend21[3].children:
    plt.contourf(np.nansum(t.get_mask(),axis=0),levels=10)
    plt.title(t.idx)
    plt.show()


for t in dend21[78].children:
    plt.contourf(np.nansum(t.get_mask(),axis=0),levels=10)
    plt.title(t.idx)
    plt.show()


#i cant believe im doing this again sighhhh, okay the thing i want to achieve is some
#way of saying: partitioning the co in this way makes it easier to process while staying
#(somewhat, or at least more) faithful to the physical stuff going on (than a rectangle in l,b)

#might need a diff minnpix for the smaller clouds?

#hmm see if the velocity gradients will end up being too too different
#once i find the clouds in here
#bc maybe there is something to be said about some of the variation in v coming bc 
#we are smooshing tg unrelated stuff thats only in the background on one side of the cloud
dends = [Dendrogram.compute(cos_3r[i],min_value=0.2,min_delta=0.2,min_npix=5e2) for i in range(22)]


plt.contourf(np.log10(c0s[10].data))


for t in dends[10].leaves:
    plt.contourf(np.nansum(t.get_mask(),axis=0),levels=np.arange(1,6))
    plt.title(t.idx)
    plt.show()


#nope never mind this is gonna get too fidgety, esp for the clouds that are made up of a bunch of 
#disparate pieces


#just gonna focus on making mom0, mom1 and mayyybe mom2 maps with the larger context of the hi 


n=12
co_lvls = [
    np.linspace(20.,21.2,n),#0
    np.linspace(20.8,22,n),
    np.linspace(20.5,21.3,n),#2
    np.linspace(20.7,21.7,n),
    np.linspace(20.5,22,n),#4
    np.linspace(20.7,21.8,n),
    np.linspace(20.5,21.5,n), #6
    np.linspace(21.1,22.2,n),
    np.linspace(20.6,21.7,n),#8
    np.linspace(20.5,21.6,n),
    np.linspace(20.6,21.7,n),#10
    np.linspace(20.8,21.6,n),
    np.linspace(20.6,22.2,n),#12
    np.linspace(20.7,22,n),
    np.linspace(20.6,22,n),#14
    np.linspace(20.6,22.4,n),
    np.linspace(20.6,22.4,n),#16
    np.linspace(20.4,21.9,n),
    np.linspace(20.7,22.2,n),#18
    np.linspace(20,21.6,n),
    np.linspace(20.5,21.8,n),#20
    np.linspace(20.8,21.8,n),
]


n=35
hi_lvls = [
    np.linspace(20.5,21.,n),#0
    np.linspace(20.7,21.8,n),
    np.linspace(20.8,21.,n),#2
    np.linspace(20.7,21.1,n),
    np.linspace(21,21.4,n),#4
    np.linspace(20.7,21.,n),
    np.linspace(21.1,21.5,n), #6
    np.linspace(21.1,21.9,n),
    np.linspace(20.9,21.6,n),#8
    np.linspace(20.9,21.3,n),
    np.linspace(20.6,21.8,n),#10
    np.linspace(20.7,21.6,n),
    np.linspace(20.7,21.6,n),#12
    np.linspace(20.7,21.6,n),
    np.linspace(20.7,21.6,n),#14
    np.linspace(20.3,21.6,n),
    np.linspace(20.7,21.8,n),#16
    np.linspace(20.8,21.8,n),
    np.linspace(20.7,21.4,n),#18
    np.linspace(19.9,21.15,n),
    np.linspace(21.2,21.7,n),#20
    np.linspace(20.8,21.6,n),
]


r3s


his_3r = [hi4pi.data[cat['vo'][i]:cat['vi'][i],r3s[0][i]:r3s[1][i],r3s[2][i]:r3s[3][i]] for i in range(22)]


r2s = np.array((cat['boh'],cat['bih'],cat['loh'],cat['lih'])) + np.array([-r/2,r/2,-r/2,r/2]).astype('int')


his_2r = [hi4pi.data[cat['vo'][i]:cat['vi'][i],r2s[0][i]:r2s[1][i],r2s[2][i]:r2s[3][i]] for i in range(22)]


h0s_2r = [np.nansum(his_2r[i],axis=0)*1.3*1.82e18 for i in range(22)]


np.log10(np.nanmean(h0s[0].data)),np.log10(np.nanmax(h0s[0].data))


[(
    plt.hist(np.ravel(np.log10(h0s[i].data)),bins=200,alpha=0.5,density=True),
    plt.hist(np.ravel(np.log10(h0s_2r[i].data)),bins=200,alpha=0.5,density=True),
#     plt.vlines([hi_lvls[i][0],hi_lvls[i][-1]],0,5,colors='red'),
    plt.vlines([np.log10(np.nanmin(h0s[i].data)),np.log10(np.nanmax(h0s[i].data))],
                                                           0,5,colors='red'),
    plt.yscale('log'),
    plt.hlines([3e-1,5e-1],hi_lvls[i][0],hi_lvls[i][-1]),
    plt.show()
    
) for i in range(22)]


r3s = np.array((cat['boh'],cat['bih'],cat['loh'],cat['lih'])) + np.array([-r,r,-r,r])


[(
  
plt.contourf(np.log10(h0s_3r[i].data),
             levels=np.linspace(np.log10(np.nanmin(h0s[i].data)),
                            np.log10(np.nanmax(h0s[i].data)),20),
            extend='both',
            extent=[big_coords[i][2],big_coords[i][3],
                    big_coords[i][0],big_coords[i][1]],cmap='viridis'),
 
plt.colorbar(),
# plt.contour(np.log10(h0s_3r[i].data),levels=hi_lvls[i],extend='both',norm='asinh',
#             extent=[big_coords[i][2],big_coords[i][3],
#                     big_coords[i][0],big_coords[i][1]],cmap='BuPu',alpha=0.3,linewidths=0.2),

    
# plt.contourf(np.log10(c0s_3r[i].data),levels=[20,23],
#         extent=[big_coords[i][2],big_coords[i][3],
#                 big_coords[i][0],big_coords[i][1]],
#             colors='black',alpha=0.15),
    
    
plt.contour(np.log10(c0s_3r[i].data),levels=[21],
        extent=[big_coords[i][2],big_coords[i][3],
                big_coords[i][0],big_coords[i][1]],
            colors='black',linewidths=0.6),

plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           colors='black',linewidths=0.6),
plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           colors='black',linewidths=0.6),

    
# plt.contour(np.log10(h0s[i].data),levels=np.arange(20,22,0.1),
#             extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
#                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
# plt.contourf(np.log10(c0s[i].data),levels=np.arange(20,22,0.25),
#             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                    cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
plt.contourf(np.log10(c0s[i].data),levels=co_lvls[i],
            cmap='magma',zorder=20,
            extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                   cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
# plt.colorbar(),
    
# plt.contour(np.log10(c0s_3r[i].data),levels=[20.5,21],
#         extent=[big_coords[i][2],big_coords[i][3],
#                 big_coords[i][0],big_coords[i][1]],colors='red'),


plt.xlim(cat['lon_HI<deg>'][i]+np.array([rs[i],-rs[i]])),
plt.ylim(cat['lat_HI<deg>'][i]+np.array([-rs[i],rs[i]])),
    
    
plt.show()

)for i in range(22)]


rs = np.around((cat['radius_HI<pc>']/cat['distance<pc>'])*u.rad.to(u.deg))


def make_fits_from_list(array,filename):
    hdulist = fits.HDUList([fits.PrimaryHDU(array[0])]+[fits.ImageHDU(array[i]) for i in range(1,22)])
    hdulist.writeto(filename,overwrite=True)


make_fits_from_list(his_2r,'hi_clouds_2r.fits')


make_fits_from_list(his_3r,'hi_clouds_3r.fits')


make_fits_from_list(cos_3r,'h2_clouds_3rhi.fits')
