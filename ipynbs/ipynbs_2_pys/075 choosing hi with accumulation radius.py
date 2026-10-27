#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Ellipse

import datacube as dc

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable


zcat = QTable(fits.open('final_cat.fits')[1].data)


zcat = zcat['names',
 'lon<deg>',
 'lat<deg>',
 'vel<km/s>',
 'distance<pc>',
 'mass_GMC<kg>',
 'radius_GMC<pc>',
 'radius_GMC<deg>',
 'vdispersion_GMC<km/s>',
 'l_pix',
 'b_pix',
 'v_pix']


zcat['radius_HI<pc>'] = ((np.array(zcat['mass_GMC<kg>'])*u.kg)/(2*np.pi*0.025*c.M_sun/(u.pc**3)))**(1/3)


zcat['radius_HI<deg>'] = (zcat['radius_HI<pc>']/(zcat['distance<pc>']*u.pc)) \
* u.rad.to(u.deg) * np.cos(np.nanmean(zcat['lat<deg>'])*u.deg.to(u.rad))


fig,ax = plt.subplots(1,1,figsize=(25,7))

for i in range(17):
    ax.add_patch(Circle((np.nanmean(zcat['lon<deg>'][i]),np.nanmean(zcat['lat<deg>'][i])),
                       zcat['radius_GMC<deg>'][i],fill=None,edgecolor='orange'))
    ax.add_patch(Circle((np.nanmean(zcat['lon<deg>'][i]),np.nanmean(zcat['lat<deg>'][i])),
                       zcat['radius_HI<deg>'][i],fill=None,edgecolor='red'))
    
plt.xlim(180,-180)
plt.ylim(-50,50)
plt.show()


colors=['steel','carolina blue','flat blue','twilight','azure','greenblue','kiwi green','light seafoam',\
        'dark cream','sun yellow','orangeish','grapefruit','lipstick','bubblegum',\
        'deep magenta','light purple','liliac']


fig,ax = plt.subplots(1,1,figsize=(27,7))

for i in range(17):
    ax.add_patch(Circle((np.nanmean(zcat['lon<deg>'][i]),np.nanmean(zcat['lat<deg>'][i])),
                       zcat['radius_GMC<deg>'][i],color='xkcd:'+colors[i]))
    ax.add_patch(Circle((np.nanmean(zcat['lon<deg>'][i]),np.nanmean(zcat['lat<deg>'][i])),
                       zcat['radius_HI<deg>'][i],fill=None,edgecolor='xkcd:'+colors[i]))
    
plt.xlim(200,-180)
plt.ylim(-50,50)
plt.show()


lpxc = np.nanmean(zcat['l_pix'],axis=1).astype('int')
bpxc = np.nanmean(zcat['b_pix'],axis=1).astype('int')

hi_rad = np.array(zcat['radius_HI<deg>']/0.125).astype('int')
co_rad = np.array(zcat['radius_GMC<deg>']/0.125).astype('int')


hi4pi,dht = fits.open('both_same_grid_720d.fits')


vels = np.arange(155)*1.3 - 100.


for i in range(17):
    sbcbhi = hi4pi.data[:,bpxc[i]-hi_rad[i]:bpxc[i]+hi_rad[i],lpxc[i]-hi_rad[i]:lpxc[i]+hi_rad[i]]
    sbcbco = dht.data[:,bpxc[i]-co_rad[i]:bpxc[i]+co_rad[i],lpxc[i]-co_rad[i]:lpxc[i]+co_rad[i]]

    spchi = np.nanmean(sbcbhi,axis=(1,2))
    spcco = np.nanmean(sbcbco,axis=(1,2))
    
    plt.plot(vels,spchi/np.nanmax(spchi),label='hi')
    plt.plot(vels,spcco/np.nanmax(spcco),label='co')
    
    plt.title(str(i)+','+zcat['names'][i])
    plt.legend()
    plt.show()


hi_cubes = [hi4pi.data[:,bpxc[i]-hi_rad[i]:bpxc[i]+hi_rad[i],lpxc[i]-hi_rad[i]:lpxc[i]+hi_rad[i]]\
            for i in range(17)]

co_cubes = [dht.data[:,bpxc[i]-co_rad[i]:bpxc[i]+co_rad[i],lpxc[i]-co_rad[i]:lpxc[i]+co_rad[i]]\
            for i in range(17)]


hilbs = [np.nansum(hi_cubes[i],axis=0) for i in range(17)]
colbs = [np.nansum(co_cubes[i],axis=0) for i in range(17)]


hilvs = [np.nansum(hi_cubes[i],axis=1) for i in range(17)]
colvs = [np.nansum(co_cubes[i],axis=1) for i in range(17)]


hispcs = [np.nansum(hi_cubes[i],axis=(1,2)) for i in range(17)]
cospcs = [np.nansum(co_cubes[i],axis=(1,2)) for i in range(17)]


lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*0.125 - 50


hilbs[0].shape


np.array(zcat['v_pix'])


meanvs = [(vels[np.where(hispcs[i]>(155*0.043))[0][0]],\
  vels[np.where(hispcs[i]>(155*0.043))[0][-1]]) for i in range(17)]


for i in range(17):
    
    fig,axs = plt.subplots(2,1,figsize=(10,7),layout='constrained')

    l0,l1 = np.array((lpxc[i]-hi_rad[i],lpxc[i]+hi_rad[i]))*(-0.125) + 360
    b0,b1 = np.array((bpxc[i]-hi_rad[i],bpxc[i]+hi_rad[i]))*0.125 - 50
    v0,v1 = np.array((zcat['v_pix'][i]))*1.3 - 100
    
    dv = np.nanstd((zcat['v_pix'][i]))*1.3
    
    vs0,vs1 = np.where(hispcs[i]>(155*5*0.043))[0][0],np.where(hispcs[i]>(155*0.043))[0][-1]


    axs[0].imshow(np.sqrt(hilvs[i].transpose()),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,3*l0,3*l1])

    axs[0].contour(vels,\
                   3*lons[lpxc[i]-hi_rad[i]:lpxc[i]+hi_rad[i]],\
                   hilvs[i].transpose(),\
                   levels=[100,200,400,800,1600,2400,3200,4000],colors=['white'])
    
    axs[0].contour(vels,\
                   3*lons[lpxc[i]-co_rad[i]:lpxc[i]+co_rad[i]],\
                   colvs[i].transpose(),levels=[1.5,2.5,3.5,5,7,10],colors=['black'],linewidths=[0.8])


    axs[0].scatter(np.nanmean(zcat['vel<km/s>'][i]),3*np.nanmean(zcat['lon<deg>'][i]),
                marker='x',color='orange')

    axs[0].vlines(widex2p5[i],3*l0,3*l1,linestyle='---',color='pink')
#     axs[0].vlines([v0,v1],3*l0,3*l1,linestyle='---',color='crimson')
#     axs[0].vlines([v0-1.5*dv,v1+1.5*dv],3*l0,3*l1,linestyle='---',color='xkcd:neon purple')
    axs[0].vlines(hivs[i],3*l0,3*l1,linestyle='---',color='xkcd:neon purple')
    

    axs[0].set_title(zcat['names'][i])
    
    
    axs[1].vlines(widex2p5[i],0,1,linestyle='--',color='pink')
    axs[1].plot(vels,hispcs[i]/np.nanmax(hispcs[i]),label='hi')
    axs[1].plot(vels,cospcs[i]/np.nanmax(cospcs[i]),label='co')
#     axs[1].vlines([v0,v1],0,1,linestyle='--',color='crimson')
#     axs[1].vlines([v0-1.5*dv,v1+1.5*dv],0,1,linestyle='--',color='xkcd:neon purple')
    axs[1].vlines(hivs[i],0,1,linestyle='--',color='xkcd:neon purple')
    
    axs[1].legend()


#when fine-tuning the mass and other calculation, choose the hi within the circle defined by the radius,
#rather than the square


zcat


zcat.rename_column('vel_GMC<km/s>','vel_HI<km/s>')
zcat.rename_column('vel<km/s>','vel_GMC<km/s>')


zcat.rename_column('lon<deg>','lon_GMC<deg>')
zcat.rename_column('lat<deg>','lat_GMC<deg>')


zcat['vel_GMC<km/s>'] = widex2p5


zcat['vel_GMC<km/s>'][2] = hivs[2]
zcat['vel_GMC<km/s>'][5] = hivs[5]
zcat['vel_GMC<km/s>'][7] = hivs[7]
zcat['vel_GMC<km/s>'][8] = hivs[8]


# monbo1, maddalena, hercules, cmaob1


hivs = np.array(([-35,25],[-17,25],[7,29],[-15,28],[-25,50],[-10,45],
[-43,25],[17,40],[-20,17],[-7,25],[-15,30],[-15,32],
[-25,40],[-18,25],[-30,20],[-20,30],[-17,30]))


widex2p5 = np.array([[np.nanmean(zcat['vel<km/s>'][i])-2.5*np.nanstd(zcat['vel<km/s>'][i]),\
  np.nanmean(zcat['vel<km/s>'][i])+2.5*np.nanstd(zcat['vel<km/s>'][i])] for i in range(17)])


zcat.rename_column('radius_GMC<pc>','radius_eff_GMC<pc>')
zcat.rename_column('radius_GMC<deg>','radius_eff_GMC<deg>')
zcat.rename_column('radius_HI<deg>','radius_acc_HI<deg>')
zcat.rename_column('radius_HI<pc>','radius_acc_HI<pc>')


zcat.write('final_cat.fits',overwrite=True)
