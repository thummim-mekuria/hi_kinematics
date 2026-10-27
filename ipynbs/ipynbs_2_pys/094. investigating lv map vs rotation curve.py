#!/usr/bin/env python
# coding: utf-8


import astropy.units as u
import astropy.constants as c
from astropy.io import fits

import numpy as np
import matplotlib.pyplot as plt

import datacube as dc


hiconv = 1.828e18#*1.3
coconv = 1.97e20#*1.3


hi4pi,dht = fits.open('both_same_grid_720d.fits')


lvh = np.nansum(hi4pi.data,1)


lons = np.arange(2880)*-0.125 + 180
lats = np.arange(800)*0.125 - 50
vels = np.arange(155)*1.3 - 100


coords


sinL = np.sin(coords.transpose()[0]*u.deg.to(u.rad))


cosB = np.cos(coords.transpose()[1]*u.deg.to(u.rad))


normV = coords.transpose()[2]/220


check = normV/(sinL*cosB)


len(coords)


len(check)


goodR = coords[check <= 0]


import pandas as pd


goodR_df = pd.DataFrame(goodR,columns=['lon','lat','vel'])


len(goodR)/1e6


len(coords)


goodR_df


plt.scatter(goodR_df['lon'][::150967],goodR_df['lat'][::150967],s=1)


plt.scatter(goodR_df['lat'][::150967],goodR_df['vel'][::150967],s=1)


plt.figure(figsize=(20,10))
plt.imshow(lvh,extent=[360,-360,-100,100],origin='lower')
plt.scatter(goodR_df['lon'][::150967],goodR_df['vel'][::150967],s=1,color='white')
plt.xlim(190,-190)
plt.vlines([0,90,180,-90],-100,100)


(np.arange(2880)*-0.125 + 180)[2160]


1440+720


goodrmsk = np.ones((155,2880))
goodrmsk[77:,720:1440] = np.zeros((78,720))
goodrmsk[:78,1440:2160] = np.zeros((78,720))


goodrmsk = np.concatenate((goodrmsk[:,1440:],goodrmsk,goodrmsk[:,:1440]),axis=1)


goodrmsk_cube = np.array([goodrmsk]*800)


goodR_hi = goodrmsk_cube.swapaxes(0,1)*hi4pi.data #galactic radii with no kinematic distance ambiguity


lb_gr = np.nansum(goodR_hi,0)
lv_gr = np.nansum(goodR_hi,1)


all_ds = pd.read_csv('whole_gal_dists.csv')


lons[720:1440]


1440+720


plt.figure(figsize=(20,10))

plt.imshow(lvh,extent=[360,-360,-100,100],origin='lower',
           cmap=dc.cmap_min0('viridis'),norm='log')
plt.xlim(200,-200)
plt.ylim(-100,100)

plt.contour(np.arange(5760)*-0.125 + 360,vels,lv_gr,levels=[400],colors='red',linewidths=0.5)

plt.plot(lons,-np.sin(lons*u.deg)*220,color='white',linestyle='--') 
#velocities more extreme than this ^(further than v = 0) give negative Rgal
plt.plot(lons,-0.5*np.sin(lons*u.deg)*220,color='white')
#velocities more extreme than this ^(further than v = 0) give Rgal > ~15kpc (so likely outside the galaxy (?))

plt.plot(lons,2*np.sin(lons*u.deg),color='black')
#velocities in the same direction as sin(l) (positive when positive etc) are at Rgal < Rsun
#which is okay for |l| < 90 but doesnt make sense for |l| > 180


vc,lc = np.array(np.meshgrid(vels,lons)).transpose().reshape(-1,2).transpose()


(vc > -0.5*np.sin(lc*u.deg)*220) 


vsinlTh_neg =  vc/np.sin(lc*u.deg) < 0    
R_less_20kpc = vc > -0.5*np.sin(lc*u.deg)*220


plt.scatter(lc*(~vsinlTh_neg&R_less_20kpc),
            vc*(vsinlTh_neg),s=1e-2)
plt.xlim(190,-190)
plt.ylim(-100,100)


rot_crmsk = np.zeros((155,2880))


v = 20 l = 50


(8.5*np.sin(170*u.deg)*220)/(5 + 220*np.sin(170*u.deg))


(7.5)**2 - (8.5*np.sin(170*u.deg))**2


plt.figure(figsize=(20,10))

plt.imshow(lb_gr,extent=[360,-360,-50,50],origin='lower',
           cmap=dc.cmap_min0('viridis'),norm='log')

plt.scatter(all_ds['lon<deg>'],all_ds['lat<deg>'],s=1,norm='asinh',cmap='binary',c=1e3*all_ds['d1<kpc>'])
plt.xlim(200,-200)
plt.colorbar()
# plt.contour(np.arange(5760)*-0.125 + 360,lats,np.nansum(hi4pi.data,axis=0),
#             levels=[400],colors='white',linewidths=0.5)
# plt.contour(np.arange(5760)*-0.125 + 360,lats,lb_gr,levels=[400],colors='red',linewidths=0.5)


#resolution at different distances
#will need this to turn l,b,d to x,y,z
#?? l,b,v -> l,b,d should be straightforward if l,b,v -> d?


del_deg = 0.125*u.deg.to(u.rad)


nearest = 100*u.pc #approx
fearest = 28*1e3*u.pc #approx, assuming R_max = 20


74389/120416


clean_ds = all_ds[(all_ds['d1_prob']>=0.9) & (all_ds['d1_unc<kpc>']/all_ds['d1<kpc>'] <=1)]


sol_nb = clean_ds[clean_ds['d1<kpc>']<=1]


plt.figure(figsize=(20,10))
plt.scatter(sol_nb['lon<deg>'],sol_nb['vel<km/s>'],cmap='Paired',marker='x',
            c=sol_nb['d1<kpc>']*1000,s=10)
plt.colorbar()

plt.contour(np.arange(5760)*-0.125 + 360,vels,lvh,levels=[400,800],colors='black',linewidths=0.5)

plt.xlim(180,-180)
plt.ylim(-50,50)


plt.scatter(clean_ds['lon<deg>'],clean_ds['lat<deg>'],c=clean_ds['d1<kpc>']*1000,s=1,norm='log')
plt.colorbar()


[(plt.contourf(lons,vels,np.log10(hiconv*hi4pi.data[:,i,:]),levels=np.arange(17,24,0.5),
           cmap=dc.cmap_min0('viridis')),
# plt.xlim(180,140),
plt.vlines(170,-25,-25+1.3,color='red'),
# plt.ylim(-30,30),
plt.show()) for i in range(0,800,11)]
