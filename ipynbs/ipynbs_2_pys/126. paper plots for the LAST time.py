#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import fitting,models

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


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


mom1_vs = np.array(([1,5,0.65/2],
[-10,5,1.3],
[-3,1,0.65/2],
[1,6,0.65],
[12,20,0.65],
[5,7,0.65/2],
[10,15,0.65],
[-3,11,1.3],
[5,15,1.3],
[-5,5,1.3],
[1,10,0.65],
[20,34,1.3],
[0,11,1.3],
[8,15,0.65],
[-1,7,0.65],
[3,15,1.3],
[6,13,0.65],
[-7,16,2.6],
[-1,10,1.3],
[-6,1,0.65],
[8,20,1.3],
[4,11,0.65]))


linx = np.array((3,3,3,3,3,3,3,4,4,2,4,3,3,3,4,4,3,3,3,3,3,4))
liny = np.array((4,4,2,4,3,4,4,3,3,4,5,4,3,4,4,3,4,3,4,4,3,4))


[( 
    plt.figure(dpi=200),
        
    plt.contourf(c1s[i].data,levels=np.arange(mom1_vs[i][0],mom1_vs[i][1],mom1_vs[i][2]/2),
                                              cmap='RdYlBu',extend='both',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$\overline{v_{_{LOS}}}\ \ [km\ s^{-1}$]',pad=0.05,
                 shrink=0.8,ticks=np.linspace(mom1_vs[i][0],mom1_vs[i][1],3),
                ),
    
    
    plt.contourf(np.log10(c0s[i].data),levels=[18,20.2],colors='white',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linewidths=[0.3,0.3,0.3,0.9],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               ),

 plt.xlim(cat['lon_GMC<deg>'][i]+np.array((0.1,-0.1))),
 plt.ylim(cat['lat_GMC<deg>'][i]+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                        cat['lon_GMC<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                              cat['lon_GMC<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                        cat['lat_GMC<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                              cat['lat_GMC<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom1_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


for ax in axs.flat:
    print(ax)


fig,axs = plt.subplots(2,3,dpi=200,figsize=(10,7))
[( 
        
    ax.contourf(np.log10(c0s[i].data),levels=np.arange(20,22.25,0.1),
                                              cmap='YlGnBu_r',extend='both',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
#     plt.colorbar(label=r'$N_{H_2}\ \ [cm^{-2}]$',pad=0.05,
#                  shrink=0.8,ticks = np.arange(20,23)),
    
    
    ax.contourf(np.log10(c0s[i].data),levels=[18,20.2],colors='white',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    ax.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linewidths=[0.3,0.3,0.3,0.9],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               ),

    ax.set_xlim(cat['lon_GMC<deg>'][i]+np.array((0.1,-0.1))),
     ax.set_ylim(cat['lat_GMC<deg>'][i]+np.array((-0.1,0.1))),

     ax.set_xticks([]),#np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
#                             cat['lon_GMC<deg>'][i][1],linx[i]),1),

#                labels=[str(tck)+r'$\degree$' for tck in
#                       np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
#                                   cat['lon_GMC<deg>'][i][1],linx[i]),1)]),

     ax.set_yticks([]),#np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
#                             cat['lat_GMC<deg>'][i][1],liny[i]),1),
#                labels=[str(tck)+r'$\degree$' for tck in
#                       np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
#                                   cat['lat_GMC<deg>'][i][1],liny[i]),1)]),

     ax.set_aspect('equal'),
    plt.subplots_adjust(wspace=0.05,hspace=0.12),

#      ax.set_xlabel('Galactic Longitude'),
#      plt.ylabel('Galactic Latitude'),

     ax.set_title(lbls[i]),
     

)for ax,i in zip(axs.flat,[20,1,9,17,12,18]) ]

plt.savefig('tile_mom0',bbox_inches='tight')


np.log10(4e20)


fig,axs = plt.subplots(3,3,dpi=200,figsize=(12,9))
[( 
        
    ax.contourf(np.log10(h0s[i].data),levels=np.arange(20.6,21.6,0.05),
                                              cmap='YlGnBu_r',extend='both',
               extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                       cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
    ax.contour(np.log10(h0s[i].data),levels=[21],
                         colors='white',linewidths=0.7,linestyle='--',
               extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                       cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
#     plt.colorbar(label=r'$N_{H_2}\ \ [cm^{-2}]$',pad=0.05,
#                  shrink=0.8,ticks = np.arange(20,23)),
    
    
#     ax.contourf(np.log10(c0s[i].data),levels=[18,20.2],colors='white',
#                extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                        cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    ax.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linewidths=[0.3,0.3,0.3,0.9],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               ),

#     ax.set_xlim(cat['lon_GMC<deg>'][i]+np.array((0.1,-0.1))),
#      ax.set_ylim(cat['lat_GMC<deg>'][i]+np.array((-0.1,0.1))),

     ax.set_xticks([]),#np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
#                             cat['lon_GMC<deg>'][i][1],linx[i]),1),

#                labels=[str(tck)+r'$\degree$' for tck in
#                       np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
#                                   cat['lon_GMC<deg>'][i][1],linx[i]),1)]),

     ax.set_yticks([]),#np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
#                             cat['lat_GMC<deg>'][i][1],liny[i]),1),
#                labels=[str(tck)+r'$\degree$' for tck in
#                       np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
#                                   cat['lat_GMC<deg>'][i][1],liny[i]),1)]),

     ax.set_aspect('equal'),
    plt.subplots_adjust(wspace=0.25,hspace=0.25),

#      ax.set_xlabel('Galactic Longitude'),
#      plt.ylabel('Galactic Latitude'),

     ax.set_title(lbls[i]),
     
    ax.invert_xaxis(),
              
)for ax,i in zip(axs.flat,[20,1,9,17,12,18,3,2]) ]

plt.savefig('tile_mom0_hi',bbox_inches='tight')


[( 
    plt.figure(dpi=200),
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20,22.25,0.1),
                                              cmap='YlGnBu_r',extend='both',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$N_{H_2}\ \ [cm^{-2}]$',pad=0.05,
                 shrink=0.8,ticks = np.arange(20,23)),
    
    
    plt.contourf(np.log10(c0s[i].data),levels=[18,20.2],colors='white',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linewidths=[0.3,0.3,0.3,0.9],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               ),

 plt.xlim(cat['lon_GMC<deg>'][i]+np.array((0.1,-0.1))),
 plt.ylim(cat['lat_GMC<deg>'][i]+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                        cat['lon_GMC<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                              cat['lon_GMC<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                        cat['lat_GMC<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                              cat['lat_GMC<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom0_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


[( 
    plt.figure(dpi=200),
    
    plt.contourf(c2s[i].data,
                 levels=np.linspace(max(0,np.nanmean(c2s[i].data)-1.5*np.nanstd(c2s[i].data)),
                        np.nanmean(c2s[i].data)+2*np.nanstd(c2s[i].data),10),
                 cmap='summer',extend='both',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$\sigma_v\ \ [km\ s^{-1}]$',pad=0.05,shrink=0.8,
                 ticks=np.around(np.linspace(max(0,np.nanmean(c2s[i].data)-1.5*np.nanstd(c2s[i].data)),
                        np.nanmean(c2s[i].data)+2*np.nanstd(c2s[i].data),4),1)),
    
    
    plt.contourf(np.log10(c0s[i].data),levels=[18,20.2],colors='white',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linewidths=[0.3,0.3,0.3,0.9],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               ),

 plt.xlim(cat['lon_GMC<deg>'][i]+np.array((0.1,-0.1))),
 plt.ylim(cat['lat_GMC<deg>'][i]+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                        cat['lon_GMC<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                              cat['lon_GMC<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                        cat['lat_GMC<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                              cat['lat_GMC<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom2_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


cat['loh'],cat['lih'] = ((cat['lon_HI<deg>']-360)/-0.125).astype('int').transpose()


cat['boh'],cat['bih'] = ((cat['lat_HI<deg>']+50)/0.125).astype('int').transpose()


fttr = fitting.LinearLSQFitter()
init = models.Linear1D()


np.array((
            np.ravel(yh).shape,
            np.ravel(h0s[i]).shape,
            np.ravel(h1s[i]).shape,
            
            ))


h0s[i]


OmHs = np.zeros(22)
OmCs = np.zeros(22)

ThHs = np.zeros(22)
ThCs = np.zeros(22)

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
    
    dvdxc = fttr(init,cxpc,cdf['v'])
    dvdyc = fttr(init,cypc,cdf['v'])
      
    
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

    OmHs[i] = omh
    OmCs[i] = omc
    
    ThHs[i] = thh
    ThCs[i] = thc
    

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

    
cat['vel_grad_GMC<km/s/pc>'] = OmCs
cat['vel_grad_HI<km/s/pc>'] = OmHs


cat['rot_ax_GMC<deg>'] = ThCs-90
cat['rot_ax_HI<deg>'] = ThHs-90

for i in range(22):
    if cat['rot_ax_GMC<deg>'][i] < -180:
        cat['rot_ax_GMC<deg>'][i] +=360
        
    if cat['rot_ax_HI<deg>'][i] < -180:
        cat['rot_ax_HI<deg>'][i] +=360


cat.write('lastlastLAST_cat.fits',overwrite=True)


from photutils.centroids import centroid_com as com


np.array(cat[['loc','lic','boc','bic']])


np.transpose(cat[['loc','lic','boc','bic']])


def calc_jsum(mom0_maps,mom1_maps,dist,cl_typ):
    n = len(mom1_maps)
    
    if cl_typ=='CO':
        lo = cat['loc']
        li = cat['lic']
        bo = cat['boc']
        bi = cat['bic']
        conv = 1.36*2
        
    elif cl_typ=='21cm':
        lo = cat['loh']
        li = cat['lih']
        bo = cat['boh']
        bi = cat['bih']
        conv = 1.36
        
        
    com_coords = [com(mom0_maps[i]).astype('int') for i in range(22)]
    
    coords = [np.array((np.meshgrid(
                 np.arange(lo[i],li[i])*-0.125+360,
                 np.arange(bo[i],bi[i])*0.125-50)
                       ))*u.deg.to(u.rad)*dist[i]
              
                    for i in range(n)]
    
    normed_h1s = [mom1_maps[i]-\
                  mom1_maps[i][com_coords[i][1],com_coords[i][0]]
                     for i in range(n)]
    
    
    mass_maps = [((mom0_maps[i]*c.m_p*conv/u.cm**2).to(u.solMass/u.pc**2))*\
            (0.125*u.deg.to(u.rad)*dist[i]*u.pc)**2
             for i in range(n)]
    
    r_vector = [np.array((coords[i][0] - np.nanmean(coords[i][0]),
                 coords[i][1] - np.nanmean(coords[i][1]),
                 np.zeros(coords[i][0].shape)))
                            for i in range(n)]
    
    v_vector = [np.array((np.zeros(normed_h1s[i].shape),
            np.zeros(normed_h1s[i].shape),normed_h1s[i]))
            for i in range(n)]
    
    
#     print(r_vector)
    jsum_kmpc_s = [np.nansum(np.cross(
                                r_vector[i].reshape(3,-1).transpose()*u.pc,
                                v_vector[i].reshape(3,-1).transpose()*u.km/u.s)*\
                             mass_maps[i].ravel()[:,np.newaxis],axis=0)/np.nansum(mass_maps[i])
                   for i in range(n)]
    
    
    return jsum_kmpc_s,mass_maps,r_vector,v_vector


hos = [h0.data for h0 in h0s]
his = [h1.data for h1 in h1s]


coos = [c0.data for c0 in c0s]
cis = [c1.data for c1 in c1s]


mass_maps_hi = calc_jsum(hos,his,cat['distance<pc>'],'21cm')[1] #equal to the ones calculated before
mass_maps_co = calc_jsum(coos,cis,cat['distance<pc>'],'CO')[1] #" " " "


r = calc_jsum(coos,cis,cat['distance<pc>'],'CO')[2]
v = calc_jsum(coos,cis,cat['distance<pc>'],'CO')[3]
m = calc_jsum(coos,cis,cat['distance<pc>'],'CO')[1]


[np.nansum(np.cross(r[i].reshape(3,-1).transpose()*u.pc,
                   v[i].reshape(3,-1).transpose()*u.km/u.s)*m[i].ravel()[:,np.newaxis],axis=0)/np.nansum(m[i])
 for i in range(22)]


calc_jsum(coos,cis,cat['distance<pc>'],'CO')[0]


np.linalg.norm(hj_vec,axis=1)


np.linalg.norm(np.array((calc_jsum(hos,his,cat['distance<pc>'],'21cm')[0])),axis=1)


np.linalg.norm(np.array((calc_jsum(coos,cis,cat['distance<pc>'],'CO')[0])),axis=1)


cj_vec = np.array(calc_jsum(coos,cis,cat['distance<pc>'],'CO')).transpose()[:2].transpose()


oldcat = QTable(fits.open('post_flash_cat.fits')[1].data)


cj_omr = (OmCs*u.km/u.s/u.pc*(cat['radius_GMC<pc>']*u.pc)**2)


hj_omr = (OmHs*u.km/u.s/u.pc*(cat['radius_HI<pc>']*u.pc)**2)


np.array([np.sqrt(x**2 + y**2) for i in range(22) for x,y in [cj_vec[i]]])/cj_omr.value


cat['ang_mom_GMC<erg.s>'] = ((cat['mass_GMC<Msun>']*\
                              cat['vel_grad_GMC<km/s/pc>']*\
                              cat['radius_GMC<pc>']**2)*u.solMass*u.km/u.s*u.pc).cgs


cat['ang_mom_HI<erg.s>'] = ((cat['mass_HI<Msun>']*\
                              cat['vel_grad_HI<km/s/pc>']*\
                              cat['radius_HI<pc>']**2)*u.solMass*u.km/u.s*u.pc).cgs


mom1_vs = np.array((
[0.5,2,0.65/4],
[-15,0,1.3],
[-2.5,-1,0.65/4],
[-0.5,3,0.65/2],
[14,18,0.65/2],
[4,7,0.65/4],
[8.5,13.5,0.65/2],
[4,13,0.65],
[10.5,18.5,0.65],
[-4.5,0,0.65/2],
[-1,7,0.65],
[27,29.5,0.65/2],
[5,8,0.65/2],
[6.5,13.5,0.65],
[0,5,0.65],
[4,10,0.65],
[5,12.5,0.65],
[3,8.5,0.65],
[0,6.5,0.65],
[-4,-1,0.65/2],
[11,14,0.65/2],
[0,7,0.65]))


outside_lats = np.array((
    [-20,-12],
    [-16.2,-2.9],
    [16.15,18.9],
    [-21.1,-11.],
    [-4.2,1.6],
    [-23.5,-16],
    [-12.,-6.5],
    [-3.6,4.9],
    [5.7,12.1],
    [-16.5,-8.],
    [3.5,25.5],
    [-5.5,-0.1],
    [-2.7,5.6],
    [-17.2,-7.9],
    [10.7,26.5],
    [-26.4,-11.8],
    [-20.,-7.1],
    [-18,-7],
    [-27.6,-14],
    [20.,31],
    [-4.6,0.3],
    [-26.7,-5.3]))


outside_lons = np.array((
    [41,28],
    [167,153.6],
    [150,144],
    [-54,-67],
    [-132.5,-138.5],
    [3.8,-2.2],
    [-138,-142.4],
    [-164.8,-173.6],
    [49,41],
    [97.9,94.5],
    [-8,-27],
    [-140.5,-146],
    [-154.2,-162.6],
    [-140.7,-150.1],
    [11.5,-11.5],
    [-141.5,-156.4],
    [-147.8,-160.6],
    [-159.,-172.5],
    [165.4,151.9],
    [127.3,117.7],
    [-151.4,-155.6],
    [180.2,158.8]))


[( 
        plt.figure(dpi=200),
    plt.contourf(h1s[i].data,levels=np.arange(mom1_vs[i][0],mom1_vs[i][1],mom1_vs[i][2]/2),
                                              cmap='RdYlBu',extend='both',
               extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                       cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$\overline{v_{_{LOS}}}\ \ [km\ s^{-1}$]',pad=0.15,
                 shrink=0.8,ticks=np.linspace(mom1_vs[i][0],mom1_vs[i][1],3),
                ),
    
    
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               linewidths=[0.3,0.3,0.3,0.9]),

 plt.xlim(outside_lons[i]),#+np.array((0.1,-0.1))),
 plt.ylim(outside_lats[i]),#+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_HI<deg>'][i][0],
                        cat['lon_HI<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_HI<deg>'][i][0],
                              cat['lon_HI<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_HI<deg>'][i][0],
                        cat['lat_HI<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_HI<deg>'][i][0],
                              cat['lat_HI<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom1_hi_hivs'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


[( 
    plt.figure(dpi=200),
    
    plt.contourf(h2s[i].data,levels=np.linspace(np.nanmean(h2s[i].data)-2*np.nanstd(h2s[i].data),
                                               np.nanmean(h2s[i].data)+2*np.nanstd(h2s[i].data),10),
                 cmap='summer',extend='both',
               extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                       cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$\sigma_v\ \ [km\ s^{-1}]$',pad=0.15,
                 shrink=0.8,ticks=np.around(np.linspace(np.nanmean(h2s[i].data)-2.5*np.nanstd(h2s[i].data),
                                               np.nanmean(h2s[i].data)+2.5*np.nanstd(h2s[i].data),4),1)),
    
    
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linestyles=['solid','solid','solid','solid'],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               linewidths=[0.3,0.3,0.3,0.9]),

 plt.xlim(outside_lons[i]),#+np.array((0.1,-0.1))),
 plt.ylim(outside_lats[i]),#+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_HI<deg>'][i][0],
                        cat['lon_HI<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_HI<deg>'][i][0],
                              cat['lon_HI<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_HI<deg>'][i][0],
                        cat['lat_HI<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_HI<deg>'][i][0],
                              cat['lat_HI<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom2_hi_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


n=25
hi_lvs = [
    np.linspace(20.6,20.8,n),#0
    np.linspace(21,21.7,n),
    np.linspace(20.75,20.95,n),#2
    np.linspace(20.7,21.,n),
    np.linspace(21.1,21.5,n),#4
    np.linspace(20.6,21.,n),
    np.linspace(21.2,21.5,n), #6
    np.linspace(21.5,21.8,n),
    np.linspace(21,21.5,n),#8
    np.linspace(20.9,21.3,n),
    np.linspace(20.6,21.4,n),#10
    np.linspace(20.7,21.5,n),
    np.linspace(20.8,21.4,n),#12
    np.linspace(21.05,21.5,n),
    np.linspace(20.9,21.3,n),#14
    np.linspace(20.7,21.4,n),
    np.linspace(21,21.4,n),#16
    np.linspace(21,21.5,n),
    np.linspace(20.7,21.2,n),#18
    np.linspace(20.3,20.8,n),
    np.linspace(21.3,21.6,n),#20
    np.linspace(20.6,21.4,n),
]


[( 
        plt.figure(dpi=200),
    plt.contourf(np.log10(h0s[i].data),levels=hi_lvs[i],
                                              cmap='YlGnBu_r',extend='both',
               extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                       cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$N_{HI}\ \ [cm^{-2}]$',pad=0.15,
                 shrink=0.8, ticks = np.around(hi_lvs[i][4::8],2)),
    
    
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
                linestyles=['solid','solid','solid','solid'],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               linewidths=[0.3,0.3,0.3,0.9]),

 plt.xlim(outside_lons[i]),#+np.array((0.1,-0.1))),
 plt.ylim(outside_lats[i]),#+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_HI<deg>'][i][0],
                        cat['lon_HI<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_HI<deg>'][i][0],
                              cat['lon_HI<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_HI<deg>'][i][0],
                        cat['lat_HI<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_HI<deg>'][i][0],
                              cat['lat_HI<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom0_hi_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


mom0_bins = np.arange(18.5,22.5,0.1) 

mom0_hi_pdfs = [(
    
    np.histogram(np.ravel(np.log10(h0s[i].data)),bins=mom0_bins,density=True)[0]
    
            ) for i in range(22)]

mom0_h2_pdfs = [(
    
    np.histogram(np.ravel(np.log10(c0s[i].data)),bins=mom0_bins,density=True)[0]
    
            ) for i in range(22)]


[( 
        plt.figure(dpi=200),
    plt.step(mom0_bins[:-1],
        mom0_hi_pdfs[i],#/np.nanmax(mom0_hi_pdfs[i]),
            alpha=0.3,linewidth=6,label='HI',color='xkcd:cerulean'),
    
    
    plt.step(mom0_bins[:-1],
        mom0_h2_pdfs[i],#/np.nanmax(mom0_h2_pdfs[i]),
            label=r'$H_2$',color='xkcd:cerulean'),
    

 plt.gca().set_aspect('equal'),
    
 plt.xlabel(r'$N_{H}\ \ [cm^{-2}]$'),
 plt.ylabel(r'$P\ (\ \log\ N_H\ )$'),
#  plt.legend(),
 plt.yscale('log'),
 plt.ylim(1e-2,),

    
 plt.title(names[i]),
 plt.savefig('mom0_pdf_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


dv = 0.65/3

mom1_pdf_vs = [np.arange(-1,7,dv),
np.arange(-17,7,dv),
np.arange(-5,2,dv),
np.arange(-3,7,dv),
np.arange(10,22,dv),
np.arange(3,8,dv),
np.arange(7,18,dv),
np.arange(-10,15,dv),
np.arange(3,21,dv),
np.arange(-5,5,dv),
np.arange(-2,11,dv),
np.arange(19,35,dv),
np.arange(-5,15,dv),
np.arange(5,15,dv),
np.arange(-2,8,dv),
np.arange(2,14,dv),
np.arange(3,15,dv),
np.arange(-6,16,dv),
np.arange(-3,11,dv),
np.arange(-8,6,dv),
np.arange(4,23,dv),
np.arange(-4,13,dv)]


mom1_hi_pdfs = [(
    
    np.histogram(np.ravel(h1s[i].data),bins=mom1_pdf_vs[i],density=True)[0]
    
            ) for i in range(22)]

mom1_h2_pdfs = [(
    
    np.histogram(np.ravel(c1s[i].data),bins=mom1_pdf_vs[i],density=True)[0]
    
            ) for i in range(22)]


[( 
        plt.figure(dpi=200),
    plt.step(mom1_pdf_vs[i][:-1],
        mom1_hi_pdfs[i],#/np.nanmax(mom0_hi_pdfs[i]),
            alpha=0.3,linewidth=6,label='HI',color='xkcd:scarlet'),
    
    
    plt.step(mom1_pdf_vs[i][:-1],
        mom1_h2_pdfs[i],#/np.nanmax(mom0_h2_pdfs[i]),
            label=r'$H_2$',color='xkcd:scarlet'),
    

#  plt.gca().set_aspect('equal'),
    
 plt.xlabel(r'$\overline{V_{LOS}}\ \ [km\ s^{-1}]$'),
 plt.ylabel(r'$P\ (\ \overline{V_{LOS}}\ )$'),
#  plt.legend(),
 plt.yscale('log'),
 plt.ylim(1e-2,),

    
 plt.title(names[i]),
 plt.savefig('mom1_pdf_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


mom2_bins = np.arange(0,16,0.1) 

mom2_hi_pdfs = [(
    
    np.histogram(np.ravel(h2s[i].data),bins=mom2_bins,density=True)[0]
    
            ) for i in range(22)]

mom2_h2_pdfs = [(
    
    np.histogram(np.ravel(c2s[i].data),bins=mom2_bins,density=True)[0]
    
            ) for i in range(22)]


x2s = np.array(([0,5],
[0,15],
[0,10],
[0,8],
[0,7],
[0,8],
[0,10],
[0,15],
[0,17],
[0,10],
[0,10],
[0,7],
[0,7],
[0,10],
[0,8],
[0,10],
[0,12],
[0,10],
[0,10],
[0,7],
[0,12],
[0,8]))


[( 
      plt.figure(dpi=200),
    plt.step(mom2_bins[:-1],
        mom2_hi_pdfs[i],#/np.nanmax(mom0_hi_pdfs[i]),
            label='HI',color='xkcd:emerald',linewidth=6,alpha=0.3),
    
    
    plt.step(mom2_bins[:-1],
        mom2_h2_pdfs[i],#/np.nanmax(mom0_h2_pdfs[i]),
            label=r'$H_2$',color='xkcd:emerald'),
    

#  plt.gca().set_aspect('equal'),
    
 plt.xlabel(r'$\sigma_{_{V,\ LOS}}\ \ [km\ s^{-1}]$'),
 plt.ylabel(r'$P\ (\ \sigma_{_{V,\ LOS}}\ )$'),
#  plt.legend(),
 plt.yscale('log'),
 plt.ylim(1e-2,),#10**(0.4)),
 plt.xlim(x2s[i]),
    
 plt.title(names[i]),
 plt.savefig('mom2_pdf_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


fttr(init,np.log10(cat['mass_GMC<Msun>']),
            np.log10((cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value))


fttr(init,np.log10(cat['mass_HI<Msun>']),
            np.log10((cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value))


(10**(20)*np.arange(1e2,1e7,1e3)**(0.5)).astype('float')


plt.figure(dpi=200)

plt.fill_between(np.arange(1e2,1e7,1e3),
                (10**(22.1)*np.arange(1e2,1e7,1e3)**(0.5)).astype('float'),
                (10**(22.1)*np.arange(1e2,1e7,1e3)**(0.9)).astype('float'),
                alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['mass_HI<Msun>'],
            (cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=150,alpha=0.3,color='black',linewidth=0)


plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
        color='black')


plt.scatter(cat['mass_GMC<Msun>'],
            (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=10,color='black')


plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

plt.xticks(ticks=[1e3,1e4,1e5,1e6])
plt.xlim(3e2,7e6)
plt.ylim(5e23,4e26)
# plt.yticks(ticks=[1e-1,1e0,1e1])

plt.xlabel(r'$M\ [M_\odot]$')
plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.savefig('jM',bbox_inches='tight')
plt.show()


fttr(init,np.log10(cat['radius_HI<pc>']),
            np.log10((cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value))


fttr(init,np.log10(cat['radius_GMC<pc>']),
            np.log10((cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value))


#it for some reason is getting overwhelmed with 10^20 and up
plt.figure(dpi=200)
plt.fill_between(np.arange(3e0,3e2,3e1),
                10**(22.6)*np.arange(3e0,3e2,3e1)**(1.2),
                10**(23.1)*np.arange(3e0,3e2,3e1)**(1.8),
                alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['radius_HI<pc>'],
            (cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=150,alpha=0.3,color='black',linewidth=0)


plt.plot(np.arange(3e0,3e2,1e1),10**(23.1)*np.arange(3e0,3e2,1e1)**(1.3),linewidth=6,alpha=0.3,
        color='black')


plt.scatter(cat['radius_GMC<pc>'],
            (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=10,color='black')


plt.plot(np.arange(3e0,3e2,1e1),10**(22.6)*np.arange(3e0,3e2,1e1)**(1.8),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
plt.xlim(4e0,2e2)
plt.ylim(3e23,5e26)
# plt.yticks(ticks=[1e-1,1e0,1e1])

plt.xlabel(r'$R\ [pc]$')
plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.savefig('jR',bbox_inches='tight')
plt.show()


(1*u.pc*u.km/u.s).cgs


5e2*u.AU.to(u.pc)


fttr(init,np.log10(cat['mass_GMC<Msun>']/cat['radius_GMC<pc>']**2),
            np.log10((cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value))


fttr(init,np.log10(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2),
            np.log10((cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value))


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
            (cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=150,alpha=0.3,color='black',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


plt.scatter(cat['mass_GMC<Msun>']/cat['radius_GMC<pc>']**2,
            (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
plt.xlim(5e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

plt.xlabel(r'$\Sigma\ [M_\odot\ pc^{-2}]$')
plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


Lgmc = (cat['mass_GMC<Msun>']*u.solMass*cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs


Lhi = (cat['mass_HI<Msun>']*u.solMass*cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs


fttr(init,np.log10((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs.value),
            np.log10(Lgmc.value))


fttr(init,np.log10((cat['mass_HI<Msun>']*u.solMass/cat['radius_HI<pc>']**3/u.pc**3/c.m_p).cgs.value),
            np.log10(Lhi.value))


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)
plt.figure(dpi=200),

plt.scatter((cat['mass_HI<Msun>']*u.solMass/cat['radius_HI<pc>']**3/u.pc**3/c.m_p).cgs,
            (cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=150,alpha=0.3,color='black',linewidth=0)


plt.plot(np.arange(5e0,5e2,1e1),10**(26.9)*np.arange(5e0,5e2,1e1)**(-0.98),linewidth=6,alpha=0.3,
        color='black')


plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
            (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs,
           s=10,color='black')


plt.plot(np.arange(5e0,5e2,1e1),10**(26.5)*np.arange(5e0,5e2,1e1)**(-0.7),color='black')


plt.loglog()

plt.xticks(ticks=[1e1,1e2,1e3])
plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.savefig('jn',bbox_inches='tight')

#find relationship in literature between j & n
#fill between in the back, then save
plt.show()


hi4pi,dht = fits.open('both_same_grid_720d.fits')


hi_specs = [np.nanmean(hi4pi.data[:,bo:bi,lo:li],axis=(1,2))
          
          
          for i in range(22) 
          for bo,bi,lo,li 
              in [list(cat[['boh','bih','loh','lih']][i])] ]


co_specs = [np.nanmean(dht.data[:,bo:bi,lo:li],axis=(1,2))
          
          
          for i in range(22) 
          for bo,bi,lo,li 
              in [list(cat[['boc','bic','loc','lic']][i])] ]


vlims = np.array(([-30,50],
[-70,30],
[-40,30],
[-50,30],
[-20,100],
[-20,40],
[-30,70],
[-40,50],
[-100,80],
[-80,40],
[-60,40],
[-20,100],
[-50,80],
[-30,70],
[-40,40],
[-40,50],
[-50,70],
[-40,50],
[-60,40],
[-100,40],
[-30,80],
[-70,50]))


[(
# plt.figure(dpi=200),
plt.plot(np.arange(155)*1.3-100,hi_specs[i],linewidth=6,alpha=0.3,color='black'),
plt.plot(np.arange(155)*1.3-100,50*co_specs[i],color='black',linewidth=0.9),
plt.vlines(cat['vel<km/s>'][i],0,np.nanmax(hi_specs[i]),linestyle='--',color='red',linewidth=0.7),
plt.xlabel('$V_{LOS}\ [km\ s^{-1}]$'),
plt.ylabel('$T_B\ [K]$'),
plt.title(names[i]),
plt.xlim(vlims[i]),
# plt.xticks([-100,-50,0,50,100]+list(cat['vel<km/s>'][i])),
# plt.savefig('spec_'+lbls[i],bbox_inches='tight'),
plt.show())
for i in range(22)]


delthet = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']


for i in range(22):
    if delthet[i]>180:
        delthet[i]-=360
        

np.arccos(0.5)*u.rad.to(u.deg)


np.arccos(0.25)*u.rad.to(u.deg)


np.cos(75*u.deg)


np.cos(60*u.deg)


plt.scatter(np.arange(0,180,10),np.cos(np.arange(0,180,10)*u.deg))
plt.vlines(60,-1,1)
plt.vlines(60,-1,1)


plt.figure(dpi=200)
plt.hist(delthet,bins=np.arange(-180,210,30),edgecolor='black',color=[0,0,0,0.2])
plt.xticks([-180,-120,-60,0,60,120,180])
plt.xlabel(r'$\phi_{_{GMC}} - \phi_{_{HI}}$')
plt.vlines([-60,60],0,5,color='red',linestyle='--',linewidth=0.7)
plt.ylabel(r'$N_{clouds}$')
plt.ylim(0,5)
# plt.savefig('hist_delphi')
plt.show()


hi_specs_pdf = pd.DataFrame(hi_specs,columns=np.around((np.arange(155)*1.3-100),1).astype('str'))


co_specs_pdf = pd.DataFrame(co_specs,columns=np.around((np.arange(155)*1.3-100),1).astype('str'))


hi_specs_pdf['name'] = cat['name']


co_specs_pdf['name'] = cat['name']


hi_specs_pdf.to_csv('hi_specs.csv',index=False)
co_specs_pdf.to_csv('co_specs.csv',index=False)


cat['ang_mom_GMC<erg.s>'] = Lgmc
cat['ang_mom_HI<erg.s>'] = Lhi


cat.write('lastlastLAST_cat.fits',overwrite=True)


(Lgmc/(cat['mass_GMC<Msun>']*u.solMass)).cgs


fttr(init,np.log10((Lgmc/(cat['mass_GMC<Msun>']*u.solMass)).cgs.value),
    np.log10((Lhi/(cat['mass_HI<Msun>']*u.solMass)).cgs.value))


fttr(init,np.log10(Lgmc.value),
    np.log10(Lhi.value))


(Lgmc/(cat['mass_GMC<Msun>']*u.solMass)).cgs/(Lhi/(cat['mass_HI<Msun>']*u.solMass)).cgs


fig,axs = plt.subplots(1,2,figsize=(10,5),dpi=200)

axs[0].plot(np.arange(1e23,1e27,1e24),
         10**(9.8)*np.arange(1e23,1e27,1e24)**(0.6),
         color='black',linestyle='--',linewidth=0.8)

axs[0].plot(np.arange(1e23,1e27,1e24),
         np.arange(1e23,1e27,1e24),
         color='red',linestyle='--',linewidth=0.8)


axs[0].scatter((Lgmc/(cat['mass_GMC<Msun>']*u.solMass)).cgs,
            (Lhi/(cat['mass_HI<Msun>']*u.solMass)).cgs,
           s=150,color=[0,0,0,0.2],edgecolor='black')


axs[0].set_xlabel(r'$j_{GMC}\ [cm^2\ s^{-1}]$')
axs[0].set_ylabel(r'$j_{HI}\ [cm^2\ s^{-1}]$')

axs[0].set_xlim(3e23,3e26)
axs[0].set_ylim(6e23,3e26)


axs[1].plot(np.arange(7e58,1e66,1e61),
         10**(-4.197)*np.arange(7e58,1e66,1e61)**(1.069),
         color='black',linestyle='--',linewidth=0.8)

axs[1].plot(np.arange(7e58,1e66,1e61),
         np.arange(7e58,1e66,1e61),
         color='red',linestyle='--',linewidth=0.8)


axs[1].scatter(Lgmc,Lhi,
           s=150,color=[0,0,0,0.2],edgecolor='black')


axs[1].set_xlabel(r'$J_{GMC}\ [erg\ s]$')
axs[1].set_ylabel(r'$J_{HI}\ [erg\ s]$')

axs[1].set_xlim(3e59,3e65)
axs[1].set_ylim(4e59,5e66)


plt.subplots_adjust(wspace=0.4,hspace=0.5)


for ax in axs.flat:
#     ax.set_aspect('equal')
    ax.loglog()

    
plt.savefig('jjJJ',bbox_inches='tight')
plt.show()


cat[[0,1,4,2,3,5,6,7,8]]['vel_grad_HI<km/s/pc>']/0.01


cat[[0,1,4,2,3,5,6,7,8]]['mass_HI<Msun>']/1e4


fttr(init, np.log10(cat['vel_dispersion_GMC<km/s>']),np.log10(Lgmc.value))


fttr(init, np.log10(cat['vel_dispersion_HI<km/s>']),np.log10(Lhi.value))


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)
plt.figure(dpi=200),

plt.scatter(cat['vel_dispersion_GMC<km/s>'],
            Lgmc,
           s=150,alpha=0.3,color='black',linewidth=0)


plt.scatter(cat['vel_dispersion_HI<km/s>'],
            Lhi,
           s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(5e0,5e2,1e1),10**(26.9)*np.arange(5e0,5e2,1e1)**(-0.98),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             Lgmc,
#            s=10,color='black')


# plt.plot(np.arange(5e0,5e2,1e1),10**(26.5)*np.arange(5e0,5e2,1e1)**(-0.7),color='black')


plt.loglog()

# plt.xticks(ticks=[1e1,1e2,1e3])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
# plt.savefig('jn',bbox_inches='tight')

#find relationship in literature between j & n
#fill between in the back, then save
plt.show()


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)
plt.figure(dpi=200),

plt.scatter((cat['mass_HI<Msun>']*u.solMass/cat['radius_HI<pc>']**3/u.pc**3/c.m_p).cgs,
            Lhi+Lgmc,
           s=150,alpha=0.3,color='black',linewidth=0)


# plt.plot(np.arange(5e0,5e2,1e1),10**(26.9)*np.arange(5e0,5e2,1e1)**(-0.98),linewidth=6,alpha=0.3,
#         color='black')


plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
            Lgmc,
           s=10,color='black')


# plt.plot(np.arange(5e0,5e2,1e1),10**(26.5)*np.arange(5e0,5e2,1e1)**(-0.7),color='black')


plt.loglog()

plt.xticks(ticks=[1e1,1e2,1e3])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
# plt.savefig('jn',bbox_inches='tight')

#find relationship in literature between j & n
#fill between in the back, then save
plt.show()


plt.hist(Lgmc/Lhi,bins=np.arange(-1,10,0.5))

plt.hist((Lgmc/(cat['mass_GMC<Msun>']*u.solMass)).cgs/\
         (Lhi/(cat['mass_HI<Msun>']*u.solMass)).cgs,bins=np.arange(-1,10,0.5),alpha=0.5)


plt.scatter(Lhi,Lgmc)
plt.plot([1e59,1e66],[1e59,1e66])
plt.loglog()


fttr(init,np.log10(cat['radius_HI<pc>']),np.log10(cat['mass_HI<Msun>']))


fttr(init,np.log10(cat['radius_GMC<pc>']),np.log10(cat['mass_GMC<Msun>']))


fttr(init,np.log10(cat['mass_HI<Msun>']),np.log10(cat['radius_HI<pc>']))


fttr(init,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['radius_GMC<pc>']))


plt.hist(np.log10((cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value),bins=np.arange(23,27,0.25),alpha=0.5)


plt.hist(np.log10((cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value),bins=np.arange(23,27,0.25),alpha=0.5)

# plt.xscale('log')

plt.show()


plt.hist(np.log10((cat['mass_HI<Msun>']*u.solMass*cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2*\
            u.km/u.s/u.pc*u.pc**2).cgs.value),bins=10,alpha=0.5)


# plt.hist(np.log10((cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs.value),bins=np.arange(23,27,0.25),alpha=0.5)

# plt.xscale('log')

plt.show()


plt.scatter(cat['radius_HI<pc>'],
            cat['mass_HI<Msun>'],
           s=150,alpha=0.3,color='black',linewidth=0)


plt.scatter(cat['radius_GMC<pc>'],
            cat['mass_GMC<Msun>'],
           s=10,color='black')


# plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
#             31.6*cat['radius_HI<pc>']**2,
#            s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
            cat['mass_GMC<Msun>'],
           s=150,alpha=0.3,color='black',linewidth=0)

plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
            31.6*cat['radius_HI<pc>']**2,
           s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


plt.hist(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,bins=20,alpha=0.5)
plt.hist(cat['mass_GMC<Msun>']/cat['radius_GMC<pc>']**2,bins=20,alpha=0.5)

plt.show()


plt.hist(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**3,bins=20,alpha=0.5)
plt.hist(cat['mass_GMC<Msun>']/cat['radius_GMC<pc>']**3,bins=20,alpha=0.5)

plt.show()


fttr(init,np.log10(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2),
np.log10(cat['radius_HI<pc>']))


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
            cat['radius_HI<pc>'],
           s=150,alpha=0.3,color='black',linewidth=0)

# plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
#             31.6*cat['radius_HI<pc>']**2,
#            s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['mass_HI<Msun>'],
            cat['mass_GMC<Msun>'],
           s=150,alpha=0.3,color='black',linewidth=0)

# plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
#             31.6*cat['radius_HI<pc>']**2,
#            s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['vel_grad_HI<km/s/pc>'],
            cat['vel_grad_GMC<km/s/pc>'],
           s=150,alpha=0.3,color='black',linewidth=0)

# plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
#             31.6*cat['radius_HI<pc>']**2,
#            s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


fttr(init,np.log10(cat['vel_grad_HI<km/s/pc>']/cat['radius_HI<pc>']**2),
    np.log10(cat['radius_GMC<pc>']**2/cat['vel_grad_GMC<km/s/pc>']))


fttr(init,np.log10(cat['mass_HI<Msun>']*cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2),
    np.log10(cat['mass_GMC<Msun>']*cat['radius_GMC<pc>']**2/cat['vel_grad_GMC<km/s/pc>']))


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['mass_HI<Msun>']*cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2,
            cat['mass_GMC<Msun>']*cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2,
           s=150,alpha=0.3,color='black',linewidth=0)

# plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
#             31.6*cat['radius_HI<pc>']**2,
#            s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()


# plt.fill_between(np.arange(1e2,1e7,1e3),
#                  10**(22.5)*np.arange(1e2,1e7,1e3)**(0.5),
#                 10**(22.5)*np.arange(1e2,1e7,1e3)**(0.9),
#                 alpha=0.1,color='black',linewidth=0)


plt.scatter(cat['vel_dispersion_GMC<km/s>'],
            cat['mass_HI<Msun>']*cat['vel_grad_HI<km/s/pc>']*cat['radius_HI<pc>']**2/(
            cat['mass_GMC<Msun>']*cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2),
           s=150,alpha=0.3,color='black',linewidth=0)

# plt.scatter(cat['mass_HI<Msun>']/cat['radius_HI<pc>']**2,
#             31.6*cat['radius_HI<pc>']**2,
#            s=150,alpha=0.3,color='red',linewidth=0)


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.6)*np.arange(1e2,1e7,1e3)**(0.53),linewidth=6,alpha=0.3,
#         color='black')


# plt.scatter((cat['mass_GMC<Msun>']*u.solMass/cat['radius_GMC<pc>']**3/u.pc**3/c.m_p).cgs,
#             (cat['vel_grad_GMC<km/s/pc>']*cat['radius_GMC<pc>']**2*\
#             u.km/u.s/u.pc*u.pc**2).cgs,
#            s=10,color='black')


# plt.plot(np.arange(1e2,1e7,1e3),10**(22.5)*np.arange(1e2,1e7,1e3)**(0.6),color='black')


plt.loglog()

# plt.xticks(ticks=[1e3,1e4,1e5,1e6])
# plt.xlim(9e0,5e2)
# plt.ylim(5e23,1e27)
# plt.yticks(ticks=[1e-1,1e0,1e1])

# plt.xlabel(r'$n_{_{H}}\ [cm^{-3}]$')
# plt.ylabel(r'$j\ [cm^2\ s^{-1}]$')
plt.show()
