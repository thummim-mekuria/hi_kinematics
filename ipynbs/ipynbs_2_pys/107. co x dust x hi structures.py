#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable

import fun as fns
import datacube as dc

hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data


fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


l0,l1 = 60,30
b0,b1 = 4,20

lod,lid = ((np.array([l0,l1])-180)/-0.229).astype('int')
bod,bid = ((np.array([b0,b1])+90)/0.229).astype('int')

herd = dust[100-68:400-68,bod:bid,lod:lid]

loch,lich = ((np.array([l0,l1])-360)/-0.125).astype('int')
boch,bich = ((np.array([b0,b1])+50)/0.125).astype('int')

herco = dht.data[:,boch:bich,loch:lich]
herhi = hi4pi.data[:,boch:bich,loch:lich]


fig,axs = plt.subplots(3,3,figsize=(15,10))

axs[0,0].contourf(np.nansum(herco,axis=0),extent=[l0,l1,b0,b1],
               levels=np.arange(0.5,15,3),cmap='binary_r')


axs[0,0].vlines([43,47],b0,b1,linestyles='--',colors=['green','red'])
axs[0,0].set_xlim(50,40)
axs[0,0].set_ylim(6,14)
axs[0,0].add_patch(Rectangle((43,7),-2,4,alpha=0.2,color='blue'))
axs[0,0].add_patch(Rectangle((47,7),-4,4,alpha=0.2,color='orange'))


axs[0,1].contourf(np.nansum(herd,axis=0),levels=10,extent=[l0,l1,b0,b1],cmap='binary_r')

axs[0,1].contour(np.nansum(herco,axis=0),extent=[l0,l1,b0,b1],colors='white',
               levels=np.arange(0.5,1,0.3),linewidths=0.5)
axs[0,1].vlines([43,47],b0,b1,linestyles='--',colors=['green','red'])    
axs[0,1].set_xlim(55,35)
axs[0,1].set_ylim(6,20)
axs[0,1].add_patch(Rectangle((43,7),-2,4,alpha=0.2,color='blue'))
axs[0,1].add_patch(Rectangle((47,7),-4,4,alpha=0.2,color='orange'))


axs[0,2].contourf(np.nansum(herhi,axis=0),levels=10,extent=[l0,l1,b0,b1],cmap='binary_r')
axs[0,2].contour(np.nansum(herco,axis=0),extent=[l0,l1,b0,b1],
               colors='white',levels=np.arange(0.5,1,0.3),linewidths=0.5)

axs[0,2].vlines([43,47],b0,b1,linestyles='--',colors=['green','red'])
axs[0,2].set_xlim(55,35)
axs[0,2].set_ylim(6,20)
axs[0,2].add_patch(Rectangle((43,7),-2,4,alpha=0.2,color='blue'))
axs[0,2].add_patch(Rectangle((47,7),-4,4,alpha=0.2,color='orange'))


for ax in axs[0]:
    ax.set_xlabel('galactic longitude [deg]')
    ax.set_ylabel('galactic latitude [deg]')

#-------------


axs[1,0].contourf(np.nansum(herco,axis=1),extent=[l0,l1,-100,100],cmap='binary_r',
                  levels=np.arange(0.5,50,5))
axs[1,0].vlines([43,47],-10,50,linestyles='--',colors=['green','red'])
axs[1,0].set_ylim(-5,20)
axs[1,0].set_xlim(50,40)
axs[1,0].add_patch(Rectangle((47,8),-6,10,alpha=0.2,color='blue'))
axs[1,0].add_patch(Rectangle((47,3),-4,5,alpha=0.2,color='orange'))
# axs[1,0].hlines([9],60,30,linestyles='--',colors='pink')


axs[1,1].contourf(np.nansum(herd,axis=1),extent=[l0,l1,100,400],cmap='binary_r',
                  levels=np.arange(5e-2,0.5,3e-2))
axs[1,1].vlines([43,47],200,280,linestyles='--',colors=['green','red'])
axs[1,1].set_ylim(210,270)
axs[1,1].set_xlim(50,40)
axs[1,1].add_patch(Rectangle((47.5,240),-7,29,alpha=0.4,fill=False,linewidth=2,color='blue'))
axs[1,1].add_patch(Rectangle((47,220),-4,20,alpha=0.4,fill=False,linewidth=2,color='orange'))


axs[1,2].contourf(np.nansum(herhi,axis=1),extent=[l0,l1,-100,100],cmap='binary_r',
                  levels=np.arange(100,6e3,5e2))
axs[1,2].vlines([43,47],-10,50,linestyles='--',colors=['green','red'])
axs[1,2].add_patch(Rectangle((47,0),-4,10,alpha=0.2,color='orange'))
axs[1,2].add_patch(Rectangle((47,8),-4,10,alpha=0.2,color='blue'))

axs[1,2].add_patch(Rectangle((47,0),-4,10,alpha=0.4,fill=False,linewidth=2,color='orange'))
axs[1,2].add_patch(Rectangle((49,8),-5.5,32,alpha=0.4,fill=False,linewidth=2,color='blue'))

axs[1,2].set_ylim(-10,50)
axs[1,2].set_xlim(50,40)


for ax in axs[1]:
    ax.set_xlabel('galactic longitude [deg]')
    ax.set_ylabel('velocity [km/s]')

axs[1,1].set_ylabel('heliocentric distance [pc]')


#-------------

axs[2,0].contourf(np.nansum(herco,axis=2).transpose(),cmap='binary_r',
                       extent=[-100,100,4,20],levels=np.arange(0.5,100,10))
axs[2,0].set_ylim(6,12)
axs[2,0].set_xlim(-5,20)
# axs[2,0].hlines([7,11],-5,20)
# axs[2,0].vlines([3,8,18],4,14)
axs[2,0].add_patch(Rectangle((8,7),10,4,alpha=0.2,color='blue'))
axs[2,0].add_patch(Rectangle((3,7),5,4,alpha=0.2,color='orange'))
# axs[2,0].vlines([9],4,20,linestyles='--',colors='pink')


axs[2,1].contourf(np.nansum(herd,axis=2).transpose(),cmap='binary_r',
                extent=[100,400,4,20],levels=np.arange(0.1,1.5,0.1))
axs[2,1].set_xlim(210,270)
axs[2,1].set_ylim(6,16)

# axs[2,1].hlines([7,11,14],200,280)
axs[2,1].add_patch(Rectangle((240,7),29,8,alpha=0.4,fill=False,linewidth=2,color='blue'))
axs[2,1].add_patch(Rectangle((220,7),20,4,alpha=0.4,fill=False,linewidth=2,color='orange'))


axs[2,2].contourf(np.nansum(herhi,axis=2).transpose(),cmap='binary_r',
                    extent=[-100,100,4,20],levels=np.arange(100,14e3,800))
axs[2,2].set_xlim(-10,50)
axs[2,2].set_ylim(6,20)
axs[2,2].add_patch(Rectangle((3,7),5,4,alpha=0.2,color='orange'))
axs[2,2].add_patch(Rectangle((8,7),10,4,alpha=0.2,color='blue'))

axs[2,2].add_patch(Rectangle((0,7),8,6,alpha=0.4,fill=False,linewidth=2,color='orange'))
axs[2,2].add_patch(Rectangle((8,10),32,6,alpha=0.4,fill=False,linewidth=2,color='blue'))

# axs[2,2].vlines([9],4,20,linestyles='--',colors='pink')


for ax in axs[2]:
    ax.set_ylabel('galactic latitude [deg]')
    ax.set_xlabel('velocity [km/s]')

axs[2,1].set_xlabel('heliocentric distance [pc]')


axs[0,0].set_title(r'$H_2$')
axs[0,1].set_title(r'dust')
axs[0,2].set_title(r'HI')


fig.tight_layout()


plt.plot(np.arange(155)*1.3-100,np.nanmean(herhi,axis=(1,2)))
plt.plot(np.arange(155)*1.3-100,400*np.nanmean(herco,axis=(1,2)))
plt.xlim(-15,)


plt.plot(np.arange(155)*1.3-100,np.nanmean(calhi,axis=(1,2)))
plt.plot(np.arange(155)*1.3-100,200*np.nanmean(calco,axis=(1,2)))
plt.xlim(-100,25)


l0,l1 = 170,140
b0,b1 = -15,-2

lod,lid = ((np.array([l0,l1])-180)/-0.229).astype('int')
bod,bid = ((np.array([b0,b1])+90)/0.229).astype('int')

cald = dust[300-68:600-68,bod:bid,lod:lid]

loch,lich = ((np.array([l0,l1])-360)/-0.125).astype('int')
boch,bich = ((np.array([b0,b1])+50)/0.125).astype('int')

calco = dht.data[:,boch:bich,loch:lich]
calhi = hi4pi.data[:,boch:bich,loch:lich]


fig,axs = plt.subplots(3,3,figsize=(15,10))

axs[0,0].contourf(np.nansum(calco,axis=0),extent=[l0,l1,b0,b1],
               levels=np.arange(1.5,15,3),cmap='binary_r')


axs[0,0].set_xlim(l0,153)
axs[0,0].set_ylim(-14,-5)
# axs[0,0].add_patch(Rectangle((43,7),-2,4,alpha=0.4))
axs[0,0].add_patch(Rectangle((170,-9),-6,3,alpha=0.2,color='orange'))
axs[0,0].add_patch(Rectangle((160,-13),-5,5,alpha=0.2,color='red'))
axs[0,0].add_patch(Rectangle((164,-10),-4,3,alpha=0.2,color='blue'))


axs[0,1].contourf(np.nansum(cald,axis=0),extent=[l0,l1,b0,b1],levels=np.arange(0.09,0.9,8e-2),
                  cmap='binary_r')

axs[0,1].contour(np.nansum(calco,axis=0),extent=[l0,l1,b0,b1],colors='white',
               levels=np.arange(1.5,2,0.3),linewidths=0.5)
    
axs[0,1].set_xlim(l0,153)
axs[0,2].set_ylim(-14,-5)
# axs[0,2].add_patch(Rectangle((43,7),-2,4,alpha=0.4))
# axs[0,2].add_patch(Rectangle((47,7),-4,4,alpha=0.2,color='orange'))
axs[0,1].add_patch(Rectangle((170,-9),-6,3,alpha=0.2,color='orange'))
axs[0,1].add_patch(Rectangle((160,-13),-5,5,alpha=0.2,color='red'))
axs[0,1].add_patch(Rectangle((164,-10),-4,3,alpha=0.2,color='blue'))


axs[0,2].contourf(np.nansum(calhi,axis=0),extent=[l0,l1,b0,b1],levels=15,
                  cmap='binary_r')
axs[0,2].contour(np.nansum(calco,axis=0),extent=[l0,l1,b0,b1],
               colors='white',levels=np.arange(1.5,2,0.3),linewidths=0.5)

# axs[0,1].hlines([14],l0,l1)
# axs[0,1].vlines([43,47],b0,b1,linestyles='--',colors=['green','red'])
axs[0,2].set_xlim(l0,150)
axs[0,2].set_ylim(-15,-5)

axs[0,2].add_patch(Rectangle((170,-9),-6,3,alpha=0.2,color='orange'))
axs[0,2].add_patch(Rectangle((160,-13),-5,5,alpha=0.2,color='red'))
axs[0,2].add_patch(Rectangle((164,-10),-4,3,alpha=0.2,color='blue'))
# axs[0,2].add_patch(Rectangle((170,-15),-7.7,20,alpha=0.4,fill=False,linewidth=2,edgecolor='orange'))
# axs[0,2].add_patch(Rectangle((143,-35),13.7,20,alpha=0.4,fill=False,linewidth=2,edgecolor='red'))
# axs[0,2].add_patch(Rectangle((162,-27),-5,20,alpha=0.4,fill=False,linewidth=2,edgecolor='blue'))


for ax in axs[0]:
    ax.set_xlabel('galactic longitude [deg]')
    ax.set_ylabel('galactic latitude [deg]')

#-------------


axs[1,0].contourf(np.nansum(calco,axis=1),extent=[l0,l1,-100,100],cmap='binary_r',
                  levels=np.arange(10,100,15))
# axs[1,0].vlines([43,47],-10,50,linestyles='--',colors=['green','red'])
axs[1,0].set_ylim(-15,10)
axs[1,0].set_xlim(l0,150)
axs[1,0].add_patch(Rectangle((170,0),-6,10,alpha=0.2,color='orange'))
axs[1,0].add_patch(Rectangle((160,-12),-5,12,alpha=0.2,color='red'))
axs[1,0].add_patch(Rectangle((164,-5),-4,5,alpha=0.2,color='blue'))


axs[1,1].contourf(np.nansum(cald,axis=1),extent=[l0,l1,300,600],cmap='binary_r',
                  levels=np.arange(3e-2,0.3,3e-2))
# axs[1,1].vlines([43,47],200,280,linestyles='--',colors=['green','red'])
axs[1,1].set_ylim(400,570)
axs[1,1].set_xlim(l0,150)
axs[1,1].add_patch(Rectangle((170,440),-6,100,alpha=0.4,fill=False,linewidth=2,color='orange'))
axs[1,1].add_patch(Rectangle((160,420),-5,60,alpha=0.4,fill=False,linewidth=2,color='red'))
axs[1,1].add_patch(Rectangle((164,460),-4,60,alpha=0.4,fill=False,linewidth=2,color='blue'))


axs[1,2].contourf(np.nansum(calhi,axis=1),extent=[l0,l1,-100,100],cmap='binary_r',
                  levels=np.arange(1000,6e3,5e2))

axs[1,2].add_patch(Rectangle((170,0),-6,10,alpha=0.2,color='orange'))
axs[1,2].add_patch(Rectangle((160,-12),-5,12,alpha=0.2,color='red'))
axs[1,2].add_patch(Rectangle((164,-5),-4,5,alpha=0.2,color='blue'))
axs[1,2].add_patch(Rectangle((170,-15),-7.7,20,alpha=0.4,fill=False,linewidth=2,edgecolor='orange'))
axs[1,2].add_patch(Rectangle((145,-30),11.7,20,alpha=0.4,fill=False,linewidth=2,edgecolor='red'))
axs[1,2].add_patch(Rectangle((162,-20),-5,20,alpha=0.4,fill=False,linewidth=2,edgecolor='blue'))


axs[1,2].set_ylim(-40,10)
axs[1,2].set_xlim(l0,l1)


for ax in axs[1]:
    ax.set_xlabel('galactic longitude [deg]')
    ax.set_ylabel('velocity [km/s]')

axs[1,1].set_ylabel('heliocentric distance [pc]')


# #-------------

axs[2,0].contourf(np.nansum(calco,axis=2).transpose(),cmap='binary_r',
                       extent=[-100,100,b0,b1],levels=np.arange(10,180,20))
axs[2,0].set_ylim(-14,-5)
axs[2,0].set_xlim(-15,10)

axs[2,0].add_patch(Rectangle((0,-9),10,3,alpha=0.2,color='orange'))
axs[2,0].add_patch(Rectangle((-12,-13),12,5,alpha=0.2,color='red'))
axs[2,0].add_patch(Rectangle((-5,-10),5,3,alpha=0.2,color='blue'))


axs[2,1].contourf(np.nansum(cald,axis=2).transpose(),cmap='binary_r',
                extent=[370,570,b0,b1],levels=np.arange(0.04,0.7,6e-2))
axs[2,1].set_xlim(400,570)
axs[2,1].set_ylim(-14,-5)

axs[2,1].add_patch(Rectangle((440,-9),100,3,alpha=0.4,fill=False,linewidth=2,color='orange'))
# axs[2,1].add_patch(Rectangle((440,-10),80,3,alpha=0.4,fill=False,linewidth=2,color='blue'))
axs[2,1].add_patch(Rectangle((420,-13),60,5,alpha=0.4,fill=False,linewidth=2,color='red'))


axs[2,2].contourf(np.nansum(calhi,axis=2).transpose(),cmap='binary_r',
                    extent=[-100,100,b0,b1],levels=np.arange(1e3,14e3,15e2))
axs[2,2].set_xlim(-40,10)
axs[2,2].add_patch(Rectangle((0,-9),10,3,alpha=0.2,color='orange'))
axs[2,2].add_patch(Rectangle((-12,-13),7,5,alpha=0.2,color='red'))
axs[2,2].add_patch(Rectangle((-5,-10),5,3,alpha=0.2,color='blue'))
axs[2,2].set_ylim(-14,-5)
# axs[2,2].add_patch(Rectangle((-15,-10),20,4,alpha=0.4,fill=False,linewidth=2,edgecolor='orange'))
axs[2,2].add_patch(Rectangle((-30,-13),20,7,alpha=0.4,fill=False,linewidth=2,edgecolor='red'))
# axs[2,2].add_patch(Rectangle((-20,-11),20,5,alpha=0.4,fill=False,linewidth=2,edgecolor='blue'))


for ax in axs[2]:
    ax.set_ylabel('galactic latitude [deg]')
    ax.set_xlabel('velocity [km/s]')

axs[2,1].set_xlabel('heliocentric distance [pc]')


axs[0,0].set_title(r'$H_2$')
axs[0,1].set_title(r'dust')
axs[0,2].set_title(r'HI')


fig.tight_layout()
