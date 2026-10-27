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


cat = QTable(fits.open('sighhplease.fits')[1].data)


coco = fits.open('h2_h2sized.fits')
hihi = fits.open('hi_hisized.fits')


msks = fits.open('h2_0thmom_masks.fits')
msks = [msks[i].data - 1 for i in range(22)]


c0s = [np.nansum(coco[i].data,axis=0)*1.3*1.97e20 for i in range(22)]


[(
plt.contourf(np.log10(c0s[i])*msks[i],levels=np.arange(20,22,0.15)),
plt.show())
    
    for i in range(22)]


gmc_mass = [np.nansum((c0s[i]*msks[i]/u.cm**2*\
 c.m_p*2*1.36*\
 (0.125*u.deg.to(u.rad)*\
  cat['distance<pc>'][i]*u.pc)**2).to(u.solMass)).value

 
 for i in range(22)]


oldcat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


cat['mass_GMC<Msun>'] = gmc_mass


new_rhs = np.cbrt(cat['mass_GMC<Msun>']/(0.04*np.pi))


lc = np.nanmean(cat['lon_GMC<deg>'],axis=1)
bc = np.nanmean(cat['lat_GMC<deg>'],axis=1)


np.array((cat['name'],new_rhs/oldcat['radius_HI<pc>'])).transpose()


cat['radius_HI<pc>'] = new_rhs


rh_px = ((new_rhs/cat['distance<pc>'])*u.rad.to(u.deg)/0.125).astype('int')


lc_px = ((lc-360)/-0.125).astype('int')
bc_px = ((bc+50)/0.125).astype('int')


cat['boh'] = (bc_px - rh_px).astype('int')
cat['loh'] = (lc_px - rh_px).astype('int')

cat['bih'] = (bc_px + rh_px).astype('int')
cat['lih'] = (lc_px + rh_px).astype('int')


gmc_area = np.array([np.nansum(msks[i]) for i in range(22)])*(0.125*u.deg.to(u.rad)*cat['distance<pc>']*u.pc)**2


eff_rad = np.sqrt(gmc_area/np.pi)


cat['radius_GMC<pc>'] = eff_rad


vels = np.arange(155)*1.3-100


h0s = [np.nansum(hihi[i].data,axis=0)*1.3*1.82e18 for i in range(22)]


[(
plt.contourf(np.log10(h0s[i]),levels=20),
    plt.colorbar(),
    plt.title(cat['name'][i]),
plt.show())
    
    for i in range(22)]


h1s = [fns.get_mom1(hihi[i].data,vels) for i in range(22)]


[(
    plt.plot(vels,np.nanmean(hihi[i].data,axis=(1,2))),
    plt.plot(vels,80*np.nanmean(coco[i].data,axis=(1,2))),
#     plt.hist(np.ravel(h1s[i]),alpha=0.5,bins=np.arange(-100,100,1.3)),
    plt.title(cat['name'][i]),
    plt.show())
    
    for i in range(22)]


for i in range(22):
    plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    plt.plot(vels,5e-2*hisp,color='black')
    plt.plot(vels,3e-1*np.gradient(hisp))
    plt.plot(vels,5e-1*np.gradient(np.gradient(hisp)))
    plt.hlines(0,-100,100,color='red')
    plt.plot(vels,cosp,color='black',linestyle='--')
    plt.title(cat['name'][i])
    plt.ylim(-3,3)
    plt.show()


#find the peaks


np.nansum(np.arange(-3,3,0.1)),np.nansum(np.abs(np.arange(-3,3,0.1)))


np.nansum(np.arange(3,6,0.1)),np.nansum(np.abs(np.arange(3,6,0.1)))


np.nansum(np.arange(-6,-3,0.1)),np.nansum(np.abs(np.arange(-6,-3,0.1)))


for i in range(22):
#     plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    hidv1 = np.gradient(hisp)
    hidv2 = np.gradient(hidv1)
    
    
#     plt.plot(vels,hidv1,linewidth=0.7)
#     plt.plot(vels,hidv2,linewidth=0.7)
    
#     plt.scatter(vels[np.where((hidv2<0.))],
#                 np.full(len(vels[np.where((hidv2<0.))]),0),
#                 s=5,color='red',zorder=50)
    
#     plt.scatter(vels[np.where((hidv2<0.))],
#                 hidv1[np.where((hidv2<0.))],
#                 s=2,color='red',zorder=50)

    print(np.nansum(hidv1[np.where((hidv2<0.))]),np.nansum(np.abs(hidv1[np.where((hidv2<0.))])))
    
#     plt.hlines(0,-100,100,color='red')
#     plt.plot(vels,cosp,color='black',linestyle='--')
#     plt.title(cat['name'][i])
#     plt.ylim(-3,3)
#     plt.show()


#okay just gonna do it by hand, but describe the ideal (find the locations where)
#the first derivative changes parity and crosses through zero, while simultaneously
#the second derivative is at a negative local minimum
#parity changes when the sum is < the sum of the absolute value


for i in range(22):
    plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    hidv1 = np.gradient(hisp)
    hidv2 = np.gradient(hidv1)
    
    for v in range(154):
        if (hidv1[v] > 0) and (hidv1[v+1] < 0) and (hidv2[v]<0) and (hidv2[v+1]<0):
            print(vels[v])
            plt.vlines(vels[v],-5,5,linewidth=0.5,color='red')
    
    plt.plot(vels,0.1*hisp)
    plt.plot(vels,hidv1,linewidth=0.7)
    plt.plot(vels,hidv2,linewidth=0.7)
    plt.hlines(0,-100,100,color='black',linewidth=0.6)
#     plt.vlines(vels[np.where((hidv2<-0.2))][[0,-1]],-2,2,linewidth=0.5,color='black')
    plt.show()


for i in [3]:
    plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    hidv1 = np.gradient(hisp)
    hidv2 = np.gradient(hidv1)
    
    for v in range(154):
        if (hidv1[v] > 0) and (hidv1[v+1] < 0):
            print(vels[v])
            plt.vlines(vels[v],-5,5,linewidth=0.5,color='red')
    
    plt.plot(vels,0.1*hisp)
    plt.plot(vels,hidv1,linewidth=0.7)
    plt.plot(vels,hidv2,linewidth=0.7)
    plt.hlines(0,-100,100,color='black',linewidth=0.6)
#     plt.vlines(vels[np.where((hidv2<-0.2))][[0,-1]],-2,2,linewidth=0.5,color='black')
    plt.xlim(-50,0)
    plt.ylim(-1,3)
    plt.show()


peaks = [1.4,
1.4,
-1.2,
[-22.,1.4],
[2.7,15.7,37.2,62.5],
5.3,
[10.5,36.5,52.1],
5.3,
[-57.1,4.,23.5],
0.1,
1.4,
[11.8,27.4,52.1],
[7.9,19.6],
9.2,
1.4,
6.6,
7.9,
6.6,
1.4,
-3.8,
[9.2,19.6,37.8,40.4],
4.]


for i in range(22):
    plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    hidv1 = np.gradient(hisp)
    hidv2 = np.gradient(hidv1)
    
    
    plt.plot(vels,0.1*hisp)
    plt.plot(vels,cosp*5,color='black')
#     plt.plot(vels,hidv1,linewidth=0.7)
    plt.plot(vels,hidv2,linewidth=0.7)
#     plt.plot(vels,np.gradient(hidv2),linewidth=0.7)
    plt.hlines(0,-100,100,color='black',linewidth=0.6)
    plt.vlines(peaks[i],0,np.nanmax(hisp*0.1),color='red')
#     plt.vlines(vels[np.where((hidv2<-0.2))][[0,-1]],-2,2,linewidth=0.5,color='black')
    plt.show()


#taking the pixel locations of the old vel ranges and adjusting to where the spectrum goes down to 1e-3
#of the max value
down1e3 = np.array([[ 67+3,  85+15],
       [ 46-10,  93-5],
       [ 61,  93-8],
       [ 66-17,  93-6],
       [ 81,  99+1],
       [ 65+5,  97-5],
       [ 61+11, 101-2],
       [ 53+12, 108-6],
       [ 61+9, 116+11],
       [ 60-5, 101-13],
       [ 69-39 ,  97],#first
       [ 92, 107+1],
       [ 73-6,  89+1],
       [ 65+8, 101+15],
       [ 69-16,  93+10],
       [ 70,  95+7],
       [ 61+8, 101+16],
       [ 61+6,  97+2],
       [ 65-9,  97-7],
       [ 69-16,  85+2],
       [ 66+7,  98],
       [ 69-26,  93-2]]) #first


#get rid of 
#cham peak on negative side
#hercules peak on the positive side


for i in range(22):
    plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    hidv1 = np.gradient(hisp)
    hidv2 = np.gradient(hidv1)
    
    
    plt.plot(vels,0.01*hisp)
    plt.plot(vels,cosp,color='black')
#     plt.plot(vels,5e-2*hidv2,linewidth=0.5)
    plt.hlines(0,-100,100,color='black',linewidth=0.6)
    plt.hlines(1e-3*np.nanmax(hisp),-100,100,color='red',linewidth=0.6)
    plt.fill_between(down1e3[i]*1.3-100,0,np.nanmax(hisp*0.01),alpha=0.3)
#     plt.fill_between(oldcat['vel<km/s>'][i],0,np.nanmax(hisp*0.01),alpha=0.3)
    plt.title(cat['name'][i])
    
#     plt.xlim(oldcat['vel<km/s>'][i]+np.array((-70,25)))
#     plt.ylim(0,1)
    plt.show()


for i in [3,8]:
    plt.figure(dpi=200)
    hisp = np.nanmean(hihi[i].data,axis=(1,2))
    cosp = np.nanmean(coco[i].data,axis=(1,2))
    
    hidv1 = np.gradient(hisp)
    hidv2 = np.gradient(hidv1)
    
    
    plt.plot(vels,0.001*hisp,marker='.')
    plt.plot(vels,cosp,color='black')
    
    plt.plot(vels,5e-2*hidv1,linewidth=0.5)
    plt.plot(vels,5e-2*hidv2,linewidth=0.5)
    
    plt.hlines(0,-100,100,color='black',linewidth=0.6)
#     plt.hlines(1e-3*np.nanmax(hisp),-100,100,color='red',linewidth=0.6)
    
    
    plt.fill_between(down1e3[i]*1.3-100,0,np.nanmax(hisp*0.01),alpha=0.3)
#     plt.fill_between(oldcat['vel<km/s>'][i],0,np.nanmax(hisp*0.01),alpha=0.3)
    plt.title(cat['name'][i])
    
#     plt.xlim(oldcat['vel<km/s>'][i]+np.array((-70,25)))
#     plt.ylim(-0.1,0.1)
    plt.show()


down1e3


new_vs = down1e3
new_vs[3][0] += 15
new_vs[8][1] -= 35


cat['vel<km/s>'] = np.around(new_vs*1.3-100,1)


[coco[i].data[new_vs[i][0]:new_vs[i][1]].shape for i in range(22)]


cat['vo'],cat['vi'] = new_vs.transpose()


cat.write('sighhplease.fits',overwrite=True)
