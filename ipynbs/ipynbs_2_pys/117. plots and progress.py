#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle,Ellipse

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting

import fun as fns
import datacube as dc

# his = fits.open('hi_clouds.fits')
# cos = fits.open('h2_clouds.fits')
# dust = fits.open('mean_and_std_lbd.fits')[1].data #for later


pcat = QTable(fits.open('post_flash_cat.fits')[1].data)


h0s,h1s,h2s = fits.open('mom0_hi.fits'),fits.open('mom1_hi.fits'),fits.open('mom2_hi.fits')
c0s,c1s,c2s = fits.open('mom0_co.fits'),fits.open('mom1_co.fits'),fits.open('mom2_co.fits')


coords = [np.array(np.meshgrid(np.arange(pcat['loh'][i],pcat['lih'][i])*-0.125+360,
             np.arange(pcat['boh'][i],pcat['bih'][i])*0.125-50)) for i in range(22)]


plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


from scipy import ndimage as ndi


names = ['Aquila South', 'California', 'Camelopardalis', 'Chamaeleon',
       'Canis Major OB1', 'Corona Australis', 'Crossbones', 'Gemini OB1',
       'Hercules', 'Lacerta', 'Lupus', 'Maddalena', 'Monoceros OB1', 'Monoceros R2',
       r'Ophiuchus', 'Orion A', 'Orion B', r'$\lambda$-Orionis', 'Perseus', 'Polaris',
       'Rosette', 'Taurus']


for i in range(22):
    l0c,l1c = np.array((pcat['loc'][i], pcat['lic'][i]))*-0.125+360
    b0c,b1c = np.array((pcat['boc'][i], pcat['bic'][i]))*0.125-50
    
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    
    r_map = np.sqrt((coords[i][0] - np.nanmean(pcat['lon_GMC<deg>'][i]))**2 + \
                (coords[i][1] - np.nanmean(pcat['lat_GMC<deg>'][i]))**2)
    
    hi_mask = r_map <= h0s[i].shape[0]*0.5*0.125

    
    plt.figure(figsize=(10,10),dpi=200)
    
#     plt.contourf(r_map)
    plt.title(names[i])#+' '+str((pcat['radius_HI<pc>'][i]/pcat['distance<pc>'][i])*u.rad.to(u.deg)))
#     plt.colorbar()
    
#     plt.contourf(coords[i][0])
#     plt.colorbar()
    
    
    plt.contourf(hi_mask*np.log10(h0s[i].data),
                 extent=[l0h,l1h,b0h,b1h],cmap=dc.cmap_min0('binary_r'),
                 levels=np.linspace(np.nanmin(np.log10(h0s[i].data))-0.5,
                                    np.nanmax(np.log10(h0s[i].data)),15),
                 
                extend='both')
    
#     scipy.ndimage.zoom(data, 3)
    plt.contour(np.log10(c0s[i].data),
                 extent=[l0c,l1c,b0c,b1c],
                 levels=np.linspace(20,np.nanmax(np.log10(c0s[i].data))+0.5,30),
                extend='both',cmap='magma',linewidths=0.8)
    
    
    xo,xi = max(l0h,l0c), min(l1h,l1c)
    yo,yi = min(b0h,b0c), max(b1h,b1c)
    
    plt.xlim(xo,xi)
    plt.ylim(yo,yi)
    
    plt.gca().set_aspect('equal')

    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')
    
    plt.xticks(np.linspace(xo,xi,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(xo,xi,4)])
    plt.yticks(np.linspace(yo,yi,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(yo,yi,4)])
    
    
    plt.savefig(pcat['name'][i]+'newest_0th')
    plt.show()


for i in range(22):
    l0c,l1c = np.array((pcat['loc'][i], pcat['lic'][i]))*-0.125+360
    b0c,b1c = np.array((pcat['boc'][i], pcat['bic'][i]))*0.125-50
    
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    
    r_map = np.sqrt((coords[i][0] - np.nanmean(pcat['lon_GMC<deg>'][i]))**2 + \
                (coords[i][1] - np.nanmean(pcat['lat_GMC<deg>'][i]))**2)
    
    hi_mask = r_map <= h0s[i].shape[0]*0.5*0.125

    
    plt.figure(figsize=(10,10),dpi=200)
    
    plt.title(names[i])#+' '+str((pcat['radius_HI<pc>'][i]/pcat['distance<pc>'][i])*u.rad.to(u.deg)))
#     plt.colorbar()
    
#     plt.contourf(coords[i][0])
#     plt.colorbar()
    
    
    plt.contourf(hi_mask*h1s[i].data,
                 extent=[l0h,l1h,b0h,b1h],alpha=0.5,
                 cmap=dc.cmap_min0('Paired'),
                 levels=np.linspace(np.nanmin(h1s[i].data)-0.5,
                                    np.nanmax(h1s[i].data),15),
                 
                extend='both')
    
#     scipy.ndimage.zoom(data, 3)
    plt.contourf(c1s[i].data,
                 extent=[l0c,l1c,b0c,b1c],alpha=0.5,
                 levels=np.linspace(np.nanmin(h1s[i].data)-0.5,
                                    np.nanmax(h1s[i].data),15),
                extend='both',cmap='Paired',linewidths=0.8)
    plt.colorbar()
    
    xo,xi = max(l0h,l0c), min(l1h,l1c)
    yo,yi = min(b0h,b0c), max(b1h,b1c)
    
    plt.xlim(xo,xi)
    plt.ylim(yo,yi)
    
    plt.gca().set_aspect('equal')

    plt.xlabel('Galactic Longitude')
    plt.ylabel('Galactic Latitude')
    
    plt.xticks(np.linspace(xo,xi,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(xo,xi,4)])
    plt.yticks(np.linspace(yo,yi,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(yo,yi,4)])
    
    
#     plt.savefig(pcat['name'][i]+'newest_0th')
    plt.show()


hi4pi,dht = fits.open('both_same_grid_720d.fits')
h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


plt.rcParams['font.size'] = 20


fig,axs = plt.subplots(2,2,figsize=(17,15),dpi=200,)

for ax,i in zip(axs.flat,[4,8,15,18]):
    
    l0c,l1c = np.array((pcat['loc'][i], pcat['lic'][i]))*-0.125+360
    b0c,b1c = np.array((pcat['boc'][i], pcat['bic'][i]))*0.125-50
    
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    
    r_map = np.sqrt((coords[i][0] - np.nanmean(pcat['lon_GMC<deg>'][i]))**2 + \
                (coords[i][1] - np.nanmean(pcat['lat_GMC<deg>'][i]))**2)
    
    hi_mask = r_map <= h0s[i].shape[0]*0.5*0.125

    
    ax.contourf(hi_mask*np.log10(h0s[i].data),
                 extent=[l0h,l1h,b0h,b1h],cmap=dc.cmap_min0('binary_r'),
                 levels=np.linspace(np.nanmin(np.log10(h0s[i].data))-0.5,
                                    np.nanmax(np.log10(h0s[i].data)),15),
                 
                extend='both')
    
#     scipy.ndimage.zoom(data, 3)
    ax.contour(np.log10(c0s[i].data),
                 extent=[l0c,l1c,b0c,b1c],
                 levels=np.linspace(20,np.nanmax(np.log10(c0s[i].data))+0.5,30),
                extend='both',cmap='magma',linewidths=0.8)


    ax.set_title(names[i])

    xo,xi = max(l0h,l0c), min(l1h,l1c)
    yo,yi = min(b0h,b0c), max(b1h,b1c)
    
    ax.set_xlim(xo,xi)
    ax.set_ylim(yo,yi)
    
    plt.gca().set_aspect('equal')
    

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
    
    
    if i in[4,8]:
        sz = 2
    else:
        sz = 4
        
    ax.set_xticks(np.around(np.arange(xo,xi-1,-sz)).astype('int'),
                  labels=[str(num)+r'$\degree$' 
                          for num in np.around(np.arange(xo,xi-1,-sz)).astype('int')])
    
    ax.set_yticks(np.around(np.arange(yo,yi+1,sz)).astype('int'),
                  labels=[str(num)+r'$\degree$' 
                          for num in np.around(np.arange(yo,yi+1,sz)).astype('int')])
    
# fig.supylabel(r'$T_{B,HI}\ [K]$')


fig.supxlabel('Galactic Longitude')
fig.supylabel('Galactic Latitude')
 
plt.savefig('postage_0th.png',bbox_inches='tight')
plt.show()


fig,axs = plt.subplots(2,2,figsize=(15,15),dpi=200,sharex=True,sharey=True)

for ax,i in zip(axs.flat,[4,8,15,18]):
       
    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(h3sig[:,pcat['boh'][i]:pcat['bih'][i],pcat['loh'][i]:pcat['lih'][i]],
                       axis=(1,2)),color='black',linewidth=0.7,
                       linestyle=':',label='HI envelope')

    ax2 = ax.twinx()

    ax2.plot(np.arange(155)*1.3-100,
            np.nanmean(d3sig[:,pcat['boc'][i]:pcat['bic'][i],pcat['loc'][i]:pcat['lic'][i]],
                       axis=(1,2)),color='xkcd:algae',
                       linewidth=0.7,label='GMC')


    ax2.tick_params(axis='y',colors='xkcd:algae')
    ax2.label_outer()#set_yticks([])

    ax.set_title(names[i],loc='right',y=0.9,x=0.9)
    
    ax2.set_ylim(top=1.2)
        
    
#     if i == 15:
#         ax.legend(loc=(0.5,0.2),frameon=False,fontsize='small')
#         ax2.legend(loc=(0.5,0.1),frameon=False,fontsize='small')
        
        
for ax in axs.flat:
#     ax.label_outer()
    ax.set_xlim(-40,80)
    ax.set_xticks([-30,0,30,60],[-30,0,30,60])
    
plt.subplots_adjust(wspace=0,hspace=0)
# fig.supylabel(r'$T_{B,HI}\ [K]$')
ax.set_ylabel(r'$T_{B,HI}\ [K]$',fontsize=25,labelpad=450,y=1.1)
ax2.set_ylabel(r'$T_{B,H_2}\ [K]$',color='xkcd:algae',fontsize=25,labelpad=20,y=1.1)
fig.supxlabel(r'$V_{LSR}\ [km\ s^{-1}]$')

plt.savefig('postage_specs.png',bbox_inches='tight')
plt.show()


fig,axs = plt.subplots(6,4,figsize=(30,30),dpi=200,sharex=True,sharey=True)

for ax,i in zip(axs.flat,range(22)):
       
    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(h3sig[:,pcat['boh'][i]:pcat['bih'][i],pcat['loh'][i]:pcat['lih'][i]],
                       axis=(1,2)),color='black',linewidth=0.7,
                       linestyle=':',label='HI envelope')

    ax2 = ax.twinx()

    ax2.plot(np.arange(155)*1.3-100,
            np.nanmean(d3sig[:,pcat['boc'][i]:pcat['bic'][i],pcat['loc'][i]:pcat['lic'][i]],
                       axis=(1,2)),color='xkcd:algae',
                       linewidth=0.7,label='GMC')


    ax2.tick_params(axis='y',colors='xkcd:algae')
    ax2.label_outer()#set_yticks([])

    ax.set_title(names[i],loc='left',y=0.8,x=0.1)
    
    ax2.set_ylim(top=1.2)
        
    
    if i == 3:
        ax.legend(loc=(0.5,0.2),frameon=False)
        ax2.legend(loc=(0.5,0.1),frameon=False)
        
        
for ax in axs.flat:
#     ax.label_outer()
    ax.set_xlim(-40,80)
    ax.set_xticks([-30,0,30,60],[-30,0,30,60])
    
plt.subplots_adjust(wspace=0,hspace=0)

plt.show()


pcat[['name','vel_grad_GMC<km/s/pc>']]


fig,axs = plt.subplots(2,2,figsize=(20,15),dpi=200,)

for ax,i in zip(axs.flat,[1,13,15,18]):
    
    l0c,l1c = np.array((pcat['loc'][i], pcat['lic'][i]))*-0.125+360
    b0c,b1c = np.array((pcat['boc'][i], pcat['bic'][i]))*0.125-50
    
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    
    r_map = np.sqrt((coords[i][0] - np.nanmean(pcat['lon_GMC<deg>'][i]))**2 + \
                (coords[i][1] - np.nanmean(pcat['lat_GMC<deg>'][i]))**2)
    
    hi_mask = r_map <= h0s[i].shape[0]*0.5*0.125

    
    velmphi = ax.contourf(hi_mask*h1s[i].data,
                 extent=[l0h,l1h,b0h,b1h],cmap=dc.cmap_min0('RdGy'),
                 levels=np.linspace(np.nanmin(h1s[i].data)-0.5,
                                    np.nanmax(h1s[i].data),15),
                extend='both',)
    
    
    fig.colorbar(velmphi,ax=ax,label=r'$v_{los}\ [km\ s^{-1}]$',
                 ticks = np.around(np.linspace(np.nanmin(h1s[i].data),
                                    np.nanmax(h1s[i].data)+0.5,5)).astype('int'))


    ax.contour(c0s[i].data,
                 extent=[l0c,l1c,b0c,b1c],
                 levels=[1e20,5e20,7e20],
                extend='both',colors=['black'],linewidths=0.8)


#     velmph2 = ax.contourf(c1s[i].data,
#                  extent=[l0c,l1c,b0c,b1c],
#                  levels=np.linspace(np.nanmin(c1s[i].data),
#                                     np.nanmax(c1s[i].data)+0.5,30),
#                 extend='both',cmap='PiYG')#,linewidths=0.8)


#     fig.colorbar(velmph2,ax=ax,label=r'$v_{los}\ [km\ s^{-1}]$',
#                  ticks = np.around(np.linspace(np.nanmin(c1s[i].data),
#                                     np.nanmax(c1s[i].data)+0.5,5)).astype('int'))
    
    ax.set_title(names[i])

    xo,xi = max(l0h,l0c), min(l1h,l1c)
    yo,yi = min(b0h,b0c), max(b1h,b1c)
    
    ax.set_xlim(xo,xi)
    ax.set_ylim(yo,yi)

#     ax.set_xlim(l0c,l1c)
#     ax.set_ylim(b0c,b1c)

#     ax.set_xlim(l0h,l1h)
#     ax.set_ylim(b0h,b1h)
    
    plt.gca().set_aspect('equal')
    

    plt.subplots_adjust(wspace=0.3,hspace=0.4)
    
    
    if i in [8,11]:
        sz = 2
    else:
        sz = 4
        
#     ax.set_xticks(np.around(np.arange(xo,xi,-sz)).astype('int'),
#                   labels=[str(num)+r'$\degree$' 
#                           for num in np.around(np.arange(xo,xi,-sz)).astype('int')])
    
#     ax.set_yticks(np.around(np.arange(yo,yi,sz)).astype('int'),
#                   labels=[str(num)+r'$\degree$' 
#                           for num in np.around(np.arange(yo,yi,sz)).astype('int')])
    
    
    ax.set_xticks(np.around(np.arange(l0h,l1h,-sz)).astype('int'),
                  labels=[str(num)+r'$\degree$' 
                          for num in np.around(np.arange(l0h,l1h,-sz)).astype('int')])
    
    ax.set_yticks(np.around(np.arange(b0h,b1h,sz)).astype('int'),
                  labels=[str(num)+r'$\degree$' 
                          for num in np.around(np.arange(b0h,b1h,sz)).astype('int')])
    
# fig.supylabel(r'$T_{B,HI}\ [K]$')


fig.supxlabel('Galactic Longitude')
fig.supylabel('Galactic Latitude')
 
# plt.savefig('postage_1st_h2.png',bbox_inches='tight')
plt.show()


fig,axs = plt.subplots(2,2,figsize=(20,15),dpi=200,)

for ax,i in zip(axs.flat,[13,15,18,20]):
    
    l0c,l1c = np.array((pcat['loc'][i], pcat['lic'][i]))*-0.125+360
    b0c,b1c = np.array((pcat['boc'][i], pcat['bic'][i]))*0.125-50
    
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    
    r_map = np.sqrt((coords[i][0] - np.nanmean(pcat['lon_GMC<deg>'][i]))**2 + \
                (coords[i][1] - np.nanmean(pcat['lat_GMC<deg>'][i]))**2)
    
    hi_mask = r_map <= h0s[i].shape[0]*0.5*0.125

    
    velmphi = ax.contourf(hi_mask*h1s[i].data,
                 extent=[l0h,l1h,b0h,b1h],cmap=dc.cmap_min0('RdGy'),
                 levels=np.linspace(np.nanmin(h1s[i].data)-0.5,
                                    np.nanmax(h1s[i].data),15),
                extend='both',)
    
    
    fig.colorbar(velmphi,ax=ax,label=r'$v_{los}\ [km\ s^{-1}]$',
                 ticks = np.around(np.arange(np.nanmin(h1s[i].data),
                                    np.nanmax(h1s[i].data),2)))
                #np.around(np.linspace(np.nanmin(h1s[i].data),
#                                     np.nanmax(h1s[i].data),5)).astype('int'))


    ax.contour(c0s[i].data,
                 extent=[l0c,l1c,b0c,b1c],
                 levels=[1e20,5e20,7e20],
                extend='both',colors=['black'],linewidths=0.8)


#     velmph2 = ax.contourf(c1s[i].data,
#                  extent=[l0c,l1c,b0c,b1c],
#                  levels=np.linspace(np.nanmin(c1s[i].data),
#                                     np.nanmax(c1s[i].data)+0.5,30),
#                 extend='both',cmap='PiYG')#,linewidths=0.8)


#     fig.colorbar(velmph2,ax=ax,label=r'$v_{los}\ [km\ s^{-1}]$',
#                  ticks = np.around(np.linspace(np.nanmin(c1s[i].data),
#                                     np.nanmax(c1s[i].data)+0.5,5)).astype('int'))
    
    ax.set_title(names[i])

    xo,xi = max(l0h,l0c), min(l1h,l1c)
    yo,yi = min(b0h,b0c), max(b1h,b1c)
    
#     ax.set_xlim(xo,xi)
#     ax.set_ylim(yo,yi)

#     ax.set_xlim(l0c,l1c)
#     ax.set_ylim(b0c,b1c)

    ax.set_xlim(l0h,l1h)
    ax.set_ylim(b0h,b1h)
    
    plt.gca().set_aspect('equal')
    

    plt.subplots_adjust(wspace=0.3,hspace=0.4)
    
    
    if i in [8,12]:
        sz = 2
    else:
        sz = 4
        
#     ax.set_xticks(np.around(np.arange(xo,xi,-sz)).astype('int'),
#                   labels=[str(num)+r'$\degree$' 
#                           for num in np.around(np.arange(xo,xi,-sz)).astype('int')])
    
#     ax.set_yticks(np.around(np.arange(yo,yi,sz)).astype('int'),
#                   labels=[str(num)+r'$\degree$' 
#                           for num in np.around(np.arange(yo,yi,sz)).astype('int')])
    
    
    ax.set_xticks(np.around(np.arange(l0h,l1h,-sz)).astype('int'),
                  labels=[str(num)+r'$\degree$' 
                          for num in np.around(np.arange(l0h,l1h,-sz)).astype('int')])
    
    ax.set_yticks(np.around(np.arange(b0h,b1h,sz)).astype('int'),
                  labels=[str(num)+r'$\degree$' 
                          for num in np.around(np.arange(b0h,b1h,sz)).astype('int')])
    
# fig.supylabel(r'$T_{B,HI}\ [K]$')


fig.supxlabel('Galactic Longitude')
fig.supylabel('Galactic Latitude')
 
# plt.savefig('postage_1st_h2.png',bbox_inches='tight')
plt.show()


im10 = [13,15,18,20]


plt.figure(figsize=(10,10))

plt.hist(pcat['mass_HI<Msun>']/pcat['mass_GMC<Msun>'],
         color='xkcd:light blue',bins = np.arange(0,8,0.5),
         alpha=0.7,edgecolor='xkcd:blue')


plt.ylim(0,7)
plt.yticks([0,3,6,],[0,3,6])

plt.vlines([1],0,6.5,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$M_{HI}\ /\ M_{GMC}$')

plt.ylabel(r'$N_{cloud}$')

plt.savefig('MMhist')
plt.show()


[(np.nansum(h0s[i].data*#circular mask (can't apply post in terms of
            #ratio of circle/square) because the pixels at 
            #the edges won't have the same mass per pixel as the pixels in the middle
 *1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value*\
 (0.125*u.deg.to(u.rad)*pcat['distance<pc>'][i])**2) 
 for i in range(22)] #this is the same as pcat['mass_HI<Msun>']


plt.scatter(pcat['mass_HI<Msun>']/pcat['mass_GMC<Msun>'],pcat['Sigma_HI<Msun/pc2>'])
#this is bc mGMC is used to calculate R_A, and is inversely related to it
#so basically doing MHI/Mgmc is the same as doing Sigma_HI (when calculated as total mass/total area)
# plt.loglog()

plt.show()


plt.figure(figsize=(10,10))

plt.hist(pcat['Sigma_HI<Msun/pc2>'],
         color='xkcd:light blue',bins = np.arange(5,80,5),
         alpha=0.7,edgecolor='xkcd:blue')

# plt.xscale('log')
# plt.ylim(0,7)
plt.yticks([0,3,6,],[0,3,6])

plt.vlines([10],0,6.5,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$\Sigma_{HI}\ [M_\odot\ pc^{-1]}]$')

plt.ylabel(r'$N_{cloud}$')

plt.savefig('SigmaHhist')
plt.show()


avSigH = [(np.nanmean(h0.data)*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value for h0 in h0s]


np.array((avSigH/pcat['Sigma_HI<Msun/pc2>'])).transpose() #hmmm now this is interesting -> whyy

#oo i think its bc -> there is a difference between the area of a circle with radius r
#and the area of a square with size 2r -> the latter is what i summed over to calculate MHI
#but the former is the defined r
#grrrr diff effective radiuseses
#this is a problem for later


#wait no, avSigH isn't including circle-ness vs squareness
#so it's just that the average of the pixels' sigmas is different than the total M/total A


# okay yeah this is the problem -> np.pi/4 ~ 0.785
#avSigH is the correct one because doing total mass/(pi r^2) artificially bumps up sigma
#by taking the mass in the whole square and saying that its in the circle
#my guess is originally did the circle mask and summed up the mass that way, then the avsgima from
#that would be similarish to this
#hmmm actually -> maybe this is a big thing ? bc now the velocity gradient is also being done on
#the square, and im not sure how much the velocities at the edges are contributing to the slope of v


plt.figure(figsize=(10,10))

plt.hist(avSigH,
         color='xkcd:light blue',bins = np.arange(0,60,10/3),
         alpha=0.7,edgecolor='xkcd:blue')

# plt.xscale('log')
# plt.ylim(0,7)
plt.yticks([0,2,4,6],[0,2,4,6])

plt.vlines([10],0,6,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$\Sigma_{HI}\ [M_\odot\ pc^{-1]}]$')

plt.ylabel(r'$N_{cloud}$')

plt.savefig('SigmaHhist')
plt.show()


np.array([pcat['name'],avSigH]).transpose()


plt.scatter(np.nanmean(pcat['lat_GMC<deg>'],axis=1),avSigH)

#nicee this is encouraging & is also telling us that blending is still a problem
#and this is prob because -> the gas at any given lsr velocity isnt guaranteed to only belong
#to one coherent structure -> kinematic distance ambiguity + theres so much stuff in the disk
#so chances go up of encountering something else moving at the same v_los


plt.hist(pcat['vi']-pcat['vo'])
plt.show()


hi_bspec = np.nanmean((np.nansum(h3sig,axis=0)*
                       1.3*1.82e18*\
                       c.m_p/u.cm**2).to(u.solMass/u.pc**2),axis=1)


plt.figure(figsize=(10,10))

plt.scatter(np.nanmean(pcat['lat_GMC<deg>'],axis=1),avSigH,
            s=1e2,marker='+',color='black')
plt.plot(np.arange(-50,50,0.125),hi_bspec,
         color='black')

plt.scatter(-np.nanmean(pcat['lat_GMC<deg>'],axis=1),avSigH,
            s=1e2,marker='_',color='red')
plt.plot(-np.arange(-50,50,0.125),hi_bspec,
         color='red',linestyle='--')


# plt.xticks([0.1,1,10,100],
#            labels=[str(num)+r'$\degree$' for num in [0.1,1,10,100]])

plt.ylabel(r'$\Sigma_{HI}\ [M_\odot\ pc^{-2}]$')
plt.xlabel(r'Galactic Latitude [$\degree$]')

plt.loglog()


plt.savefig('SigmaHvsb')

plt.show()


h2bspec = np.nanmean((np.nansum(d3sig,axis=0)*
                       2*1.3*1.97e20*\
                       c.m_p/u.cm**2).to(u.solMass/u.pc**2),axis=1)


plt.figure(figsize=(10,10))

plt.scatter(np.nanmean(pcat['lat_GMC<deg>'],axis=1),pcat['Sigma_GMC<Msun/pc2>'],
            s=1e2,marker='+',color='black')
plt.plot(np.arange(-50,50,0.125),h2bspec,
         color='black')

plt.scatter(-np.nanmean(pcat['lat_GMC<deg>'],axis=1),pcat['Sigma_GMC<Msun/pc2>'],
            s=1e2,marker='_',color='red')
plt.plot(-np.arange(-50,50,0.125),h2bspec,
         color='red',linestyle='--')


# plt.xticks([0.1,1,10,100],
#            labels=[str(num)+r'$\degree$' for num in [0.1,1,10,100]])

plt.ylabel(r'$\Sigma_{H2}\ [M_\odot\ pc^{-2}]$')
plt.xlabel(r'Galactic Latitude [$\degree$]')

plt.loglog()
plt.ylim(1)

plt.savefig('SigmaH2vsb')

plt.show()


avSigH = [(np.nanmean(h0.data)*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value for h0 in h0s]


plt.figure(figsize=(10,10))

plt.hist(pcat['mass_GMC<Msun>'],
         color='xkcd:light blue',bins=10**np.arange(2,7,0.5),
         alpha=0.7,edgecolor='xkcd:blue')

plt.xscale('log')
# plt.ylim(0,7)
plt.yticks([0,2,4,6],[0,2,4,6])

# plt.vlines([1e4],0,6.5,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$M_{GMC}\ [M_\odot]$')

plt.ylabel(r'$N_{cloud}$')

plt.savefig('Mhist')
plt.show()


[(
    plt.plot(np.arange(0,79),
             np.histogram(np.ravel((h0s[i].data*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value),
               bins=np.arange(0,80))[0]/np.nanmax(np.histogram(np.ravel((h0s[i].data*\
                1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value),
               bins=np.arange(0,80))[0])),
 plt.loglog(),
 plt.xlim(1,100),
    plt.title(pcat['name'][i]),
    )
 for i in range(22) if pcat['mass_GMC<Msun>'][i]<1e4]


5-10 aq,cam
10-50 cal


plt.scatter(np.nanmean(pcat['lon_GMC<deg>'],axis=1),avSigH)
#this makes sense from the perspective of selection -> towards galactic center, we explicitly
#avoid the plane, but then towards the anti-center theres ~5 clouds on the plane
#madd, mob, gem, cma --> maddalena is a great example here bc its at such an extreme velocity 
#relative to the velocity curbe that that
#helps disentangle it (think of the lv plot)


hlv = (np.nansum(h3sig,axis=1)*0.125*1.82e18*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2)


plt.contourf(hlv.value,extent=[360,-360,-50,50],levels=np.arange(10,70,2.5),
            cmap=dc.cmap_min0('Paired'))

plt.scatter(np.nanmean(pcat['lon_GMC<deg>'],axis=1),
            np.nanmean(pcat['vel<km/s>'],axis=1),
           color='black',s=avSigH)

for i in range(22):
    plt.annotate(pcat['name'][i],(np.nanmean(pcat['lon_GMC<deg>'][i]),
                np.nanmean(pcat['vel<km/s>'][i])),fontsize=10)
#            color='black',s=avSigH)


plt.xlim(-50,-180)
plt.show()


pcat


fttr = fitting.LinearLSQFitter()


fttr(init,np.log10(pcat['mass_GMC<Msun>']),np.log10(1e-21*pcat['spec_ang_mom_GMC<cm2/s>']))


fttr(init,np.log10(pcat['mass_GMC<Msun>']),np.log10(1e-21*pcat['spec_ang_mom_GMC<cm2/s>']))


fttr(init,np.log10(pcat['mass_HI<Msun>']),np.log10(1e-21*pcat['spec_ang_mom_HI<cm2/s>']))


plt.scatter(pcat['mass_GMC<Msun>']*1e-4,pcat['spec_ang_mom_HI<cm2/s>']*1e-26)


np.corrcoef(pcat['mass_GMC<Msun>'],pcat['mass_HI<Msun>'])


np.corrcoef(np.log10(pcat['mass_GMC<Msun>']),np.log10(pcat['mass_HI<Msun>']))


np.corrcoef(np.log10(pcat['mass_GMC<Msun>']),np.log10(1e-21*pcat['spec_ang_mom_HI<cm2/s>']))


np.corrcoef(np.log10(pcat['mass_GMC<Msun>']),np.log10(1e-21*pcat['spec_ang_mom_GMC<cm2/s>']))


0.78-0.25


plt.figure(figsize=(10,10))

plt.scatter(pcat['mass_HI<Msun>'],
            pcat['spec_ang_mom_HI<cm2/s>']*1e-21,
            color='red',s=10)

plt.plot(10**np.arange(2.5,8),(10**0.99)*(10**np.arange(2.5,8))**(0.7),
         linestyle='--',alpha=0.7,color='xkcd:blue')


plt.fill_between(10**np.arange(2.5,8),
                 (10**0.99)*(10**np.arange(2.5,8))**(0.5),
                 (10**0.99)*(10**np.arange(2.5,8))**(0.8),
         linestyle='--',alpha=0.15,color='xkcd:light blue')


plt.loglog()
# plt.xscale('log')
# # plt.ylim(0,7)
# # plt.yticks([0,3,6,],[0,3,6])

# plt.vlines([1],0,7.5,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$M_{HI}\ [M_\odot]$')

plt.ylabel(r'$j_{HI}\ \ [10^{21}\ cm^2\ s^{-1}]$')

# plt.savefig('jM')
plt.show()


plt.figure(figsize=(10,10))

plt.scatter(pcat['mass_GMC<Msun>'],
            pcat['spec_ang_mom_GMC<cm2/s>']*1e-21,color='black',s=10)


plt.plot(10**np.arange(2.5,7),(10**0.29)*(10**np.arange(2.5,7))**(0.78),
         linestyle='--',alpha=0.7,color='xkcd:blue')


plt.fill_between(10**np.arange(2.5,7),
                 (10**0.29)*(10**np.arange(2.5,7))**(0.5),
                 (10**0.29)*(10**np.arange(2.5,7))**(0.9),
         linestyle='--',alpha=0.15,color='xkcd:light blue')


plt.loglog()
# plt.xscale('log')
# # plt.ylim(0,7)
# # plt.yticks([0,3,6,],[0,3,6])

# plt.vlines([1],0,7.5,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$M_{GMC}\ [M_\odot]$')

plt.ylabel(r'$j_{GMC}\ \ [10^{21}\ cm^2\ s^{-1}]$')

plt.savefig('jMgmc')
plt.show()


plt.figure(figsize=(10,10))

plt.scatter(np.nanmean(pcat['lon_GMC<deg>'],axis=1),
            pcat['vel_grad_GMC<km/s/pc>'],
            color='black',s=10)

# plt.scatter(pcat['mass_GMC<Msun>'],
#             pcat['spec_ang_mom_HI<cm2/s>']/pcat['spec_ang_mom_GMC<cm2/s>'],
#             color='black',s=10)


# plt.plot(10**np.arange(2.5,7),(10**0.29)*(10**np.arange(2.5,7))**(0.78),
#          linestyle='--',alpha=0.7,color='xkcd:blue')


# plt.fill_between(10**np.arange(2.5,7),
#                  (10**0.29)*(10**np.arange(2.5,7))**(0.5),
#                  (10**0.29)*(10**np.arange(2.5,7))**(0.9),
#          linestyle='--',alpha=0.15,color='xkcd:light blue')


# plt.loglog()
plt.yscale('log')
# # plt.ylim(0,7)
# # plt.yticks([0,3,6,],[0,3,6])

# plt.vlines([1],0,7.5,color='red',linestyle='--',linewidths=0.6)

# plt.xlabel(r'$M_{GMC}\ [M_\odot]$')

# plt.ylabel(r'$j_{GMC}\ \ [10^{21}\ cm^2\ s^{-1}]$')

# plt.savefig('jM')
plt.show()


plt.figure(figsize=(10,10))

plt.hist(pcat['ang_mom_HI<erg.s>']/pcat['ang_mom_GMC<erg.s>'],
         color='xkcd:light blue',bins = 10**np.arange(-1/3,2.5,1/3),
         alpha=0.7,edgecolor='xkcd:blue')

plt.xscale('log')
# plt.ylim(0,7)
# plt.yticks([0,3,6,],[0,3,6])

plt.vlines([1],0,7.5,color='red',linestyle='--',linewidths=0.6)

plt.xlabel(r'$L_{HI}\ /\ L_{GMC}$')

plt.ylabel(r'$N_{cloud}$')

plt.savefig('LLhist')
plt.show()


np.corrcoef(np.log10(pcat['mass_GMC<Msun>']),
    np.log10(pcat['Sigma_HI<Msun/pc2>']))


plt.figure(figsize=(10,10))

plt.scatter(pcat['Sigma_HI<Msun/pc2>'],pcat['mass_GMC<Msun>']) #this is the same as plotting Sigma_H vs R_H
# plt.scatter(pcat['Sigma_HI<Msun/pc2>'],10*(pcat['radius_HI<pc>'])**2)

# plt.scatter(pcat['Sigma_HI<Msun/pc2>'],((1+pcat['Sigma_HI<Msun/pc2>']))**-1)

#          color='xkcd:light blue',bins = 10**np.arange(-0.5,2.5,0.25),
#          alpha=0.7,edgecolor='xkcd:blue')

plt.loglog()

# plt.ylim(0,7)
plt.xticks([1e1,1e2,],[r'$10^1$',r'$10^2$'])

# plt.vlines([1],0,6.5,color='red',linestyle='--',linewidths=0.6)

plt.ylabel(r'$M_{GMC}\ \ [M_\odot]$')
plt.xlabel(r'$\Sigma_{HI}\ \ [M_\odot\ pc^{-2}]$')

# plt.savefig('LLhist')

plt.show()


plt.scatter(pcat['mass_GMC<Msun>'],pcat['spec_ang_mom_GMC<cm2/s>']/pcat['spec_ang_mom_HI<cm2/s>'])

plt.scatter(pcat['mass_GMC<Msun>'][[13,15,18,20]],
            (pcat['ang_mom_GMC<erg.s>']/pcat['ang_mom_HI<erg.s>'])[[13,15,18,20]])


plt.loglog()
plt.show()


np.arange(22)[(pcat['radius_HI<pc>']/pcat['distance<pc>'])*u.rad.to(u.deg) > np.nanstd(pcat['lon_GMC<deg>'],axis=1)]


pcat['name'][(pcat['radius_HI<pc>']/pcat['distance<pc>'])*u.rad.to(u.deg) > np.nanstd(pcat['lon_GMC<deg>'],axis=1)]


np.nanstd(pcat['lon_GMC<deg>'],axis=1)


plt.plot(pcat['name'],np.nanstd(pcat['lon_GMC<deg>'],axis=1),marker='.')

plt.plot(pcat['name'],np.nanstd(pcat['lat_GMC<deg>'],axis=1),marker='.')

plt.xticks(rotation='vertical')
# plt.hlines([0.7,1.,1.3],'aquila_s','taurus')

plt.hlines([2,4,6,8,10,12],'aquila_s','taurus',color='black')


plt.show()


pcat['lon_GMC<deg>']


# # above this is plots for the paper (!! focusing on showing the stuff i've already done!)
# 
# and below is noodling around / moving forward from where i was at flash
# ie looking at linecenters, velocity dispersion on a pixel by pixel basis
# -- for the linecenters one, see if there are single or multiple peaks at each los


for i in range(22):
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    hszh = h3sig[pcat['vo'][i]:pcat['vi'][i],
                 pcat['boh'][i]:pcat['bih'][i],
                 pcat['loh'][i]:pcat['lih'][i]]
    
    hszc = d3sig[pcat['vo'][i]:pcat['vi'][i],
                 pcat['boh'][i]:pcat['bih'][i],
                 pcat['loh'][i]:pcat['lih'][i]]
    
    cszc = d3sig[pcat['vo'][i]:pcat['vi'][i],
                 pcat['boc'][i]:pcat['bic'][i],
                 pcat['loc'][i]:pcat['lic'][i]]
    
    
#     plt.figure(figsize=(7,7),dpi=200)
    plt.title(names[i])
    
    hszh1 = fns.get_mom1(hszh,np.arange(pcat['vel<km/s>'][i][0],pcat['vel<km/s>'][i][1],1.3))
    hszc1 = fns.get_mom1(hszc,np.arange(pcat['vel<km/s>'][i][0],pcat['vel<km/s>'][i][1],1.3))
    cszc1 = fns.get_mom1(cszc,np.arange(pcat['vel<km/s>'][i][0],pcat['vel<km/s>'][i][1],1.3))
    
    
    hi_mask = ~np.isnan(hszc1) #r_map <= hszh[i].shape[0]*0.5*0.125
    plt.hist(hszh1[hi_mask],bins=20,density=True)
    plt.hist(hszh1[~hi_mask],bins=20,alpha=0.5,density=True)
#     plt.hist(np.ravel(hszc1),bins=20,density=True,edgecolor='red',fill=False)
    plt.hist(np.ravel(cszc1),bins=20,density=True,edgecolor='black',fill=False)
#     plt.loglog()
    plt.show()
    
    
#     plt.contourf(hi_mask*hszh1,
#                  extent=[l0h,l1h,b0h,b1h],cmap=dc.cmap_min0('binary_r'),
#                  levels=10,#np.linspace(np.nanmin(hszh1)-0.5,
# #                                     np.nanmax(hszh1),15),
#                 extend='both')
# #     plt.colorbar()
    
# #     plt.contourf(hszc1,
# #                  extent=[l0h,l1h,b0h,b1h],
# #                  levels=10,#p.linspace(np.nanmin(hszh1)-0.5,
# # #                                     np.nanmax(hszh1),15),
# #                 extend='both',cmap='magma')#,linewidths=0.8)
    
    
# #     xo,xi = max(l0h,l0c), min(l1h,l1c)
# #     yo,yi = min(b0h,b0c), max(b1h,b1c)
    
#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.gca().set_aspect('equal')

#     plt.xlabel('Galactic Longitude')
#     plt.ylabel('Galactic Latitude')
    
#     plt.xticks(np.linspace(l0h,l1h,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(xo,xi,4)])
#     plt.yticks(np.linspace(b0h,b1h,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(yo,yi,4)])
    
#     plt.colorbar()
#     plt.show()


av_sig_x = np.zeros(22)


for i in range(22):
    l0h,l1h = np.array((pcat['loh'][i], pcat['lih'][i]))*-0.125+360
    b0h,b1h = np.array((pcat['boh'][i], pcat['bih'][i]))*0.125-50
    
    hszh = h3sig[pcat['vo'][i]:pcat['vi'][i],
                 pcat['boh'][i]:pcat['bih'][i],
                 pcat['loh'][i]:pcat['lih'][i]]
    
    hszc = d3sig[pcat['vo'][i]:pcat['vi'][i],
                 pcat['boh'][i]:pcat['bih'][i],
                 pcat['loh'][i]:pcat['lih'][i]]
    
    cszc = d3sig[pcat['vo'][i]:pcat['vi'][i],
                 pcat['boc'][i]:pcat['bic'][i],
                 pcat['loc'][i]:pcat['lic'][i]]
    
    
#     plt.title(names[i])


    hszh2 = fns.get_momN(hszh,np.arange(pcat['vel<km/s>'][i][0],pcat['vel<km/s>'][i][1],1.3),2)
    hszc2 = fns.get_momN(hszc,np.arange(pcat['vel<km/s>'][i][0],pcat['vel<km/s>'][i][1],1.3),2)
    cszc2 = fns.get_momN(cszc,np.arange(pcat['vel<km/s>'][i][0],pcat['vel<km/s>'][i][1],1.3),2)
    
    
    hi_mask = ~np.isnan(hszc2) #r_map <= hszh[i].shape[0]*0.5*0.125
#     plt.hist(np.sqrt(hszh2)[hi_mask]/np.nanmean(np.sqrt(cszc2)),bins=20,density=True)
    plt.hist(np.sqrt(hszh2)[~hi_mask]/np.nanmean(np.sqrt(cszc2)),bins=20,alpha=0.5,density=True)
    av_sig_x[i] = np.nanmean(np.sqrt(hszh2)[~hi_mask]/np.nanmean(np.sqrt(cszc2)))
#     print(np.nanmedian(np.sqrt(hszh2)[~hi_mask]/np.nanmean(np.sqrt(cszc2))),
#          np.nanmean(np.sqrt(hszh2)[~hi_mask]/np.nanmean(np.sqrt(cszc2))))
#     plt.loglog()
#     plt.hist(np.ravel(hszc1),bins=20,density=True,edgecolor='red',fill=False)
#     plt.hist(np.sqrt(np.ravel(cszc2)),bins=20,density=True,edgecolor='black',fill=False)


#     plt.contourf(hi_mask*hszh2,
#                  extent=[l0h,l1h,b0h,b1h],cmap=dc.cmap_min0('binary_r'),
#                  levels=10,#np.linspace(np.nanmin(hszh1)-0.5,
# #                                     np.nanmax(hszh1),15),
#                 extend='both')
# #     plt.colorbar()
    
#     plt.contourf(hszc2,
#                  extent=[l0h,l1h,b0h,b1h],
#                  levels=10,#p.linspace(np.nanmin(hszh1)-0.5,
# #                                     np.nanmax(hszh1),15),
#                 extend='both',cmap='magma')#,linewidths=0.8)
    

#     plt.xlim(l0h,l1h)
#     plt.ylim(b0h,b1h)
    
#     plt.gca().set_aspect('equal')

#     plt.xlabel('Galactic Longitude')
#     plt.ylabel('Galactic Latitude')
    
#     plt.xticks(np.linspace(l0h,l1h,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(xo,xi,4)])
#     plt.yticks(np.linspace(b0h,b1h,4),labels=[str(round(num,1))+r'$\degree$' for num in np.linspace(yo,yi,4)])
    
#     plt.colorbar()

#     plt.show()


np.corrcoef(np.log10(pcat['mass_GMC<Msun>']),av_sig_x)


fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.)


fttr(init,np.log10(pcat['mass_GMC<Msun>']),av_sig_x)


plt.scatter(np.log10(pcat['mass_GMC<Msun>']),av_sig_x,s=10)

plt.plot(np.arange(2,7),-4*np.arange(2,7)+25)
plt.fill_between(np.arange(2,7),-5*np.arange(2,7)+25,-3*np.arange(2,7)+25,alpha=0.5)

plt.xticks(rotation='vertical')
plt.show()


pcat['vel<km/s>'][[0,1,7,8,10,11,12,14]] 


pcat['vi'][[0,1,7,8,10,11,12,14]]= [85,93,104,108,93,104,89,93]


(np.array((pcat['vo'],pcat['vi']))*1.3-100).transpose()[[0,1,7,8,10,11,12,14]] 
