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

from photutils.aperture import CircularAnnulus

#------------------------------------------------------------
zcat = QTable(fits.open('final_cat.fits')[1].data)

hi4pi,dht = fits.open('both_same_grid_720d.fits')


#------------------------------------------------------------
hi_conv = (1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = (1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


#------------------------------------------------------------
lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50
vels = np.arange(155)*(1.3) - 100


whi = np.nansum(hi4pi.data,axis=0)*u.K*1.3*u.km/u.s*hi_conv

whole_sig = np.ravel(whi.to(u.solMass/u.pc**2))

wholels, wholebs = np.meshgrid(np.arange(360,-360,-0.125),np.arange(-50,50,0.125))
wholels = np.ravel(wholels)
wholebs = np.ravel(wholebs)

whole_sig_lb = pd.DataFrame(np.array((whole_sig,wholels,wholebs)).transpose(),\
                            columns=['sig<solMass/pc2>', 'lon<deg>', 'lat<deg>'])


all_ds = pd.read_csv('whole_gal_dists.csv')
sig_ds = all_ds.merge(whole_sig_lb,on=['lon<deg>','lat<deg>'])

x = sig_ds['d1<kpc>']*np.sin(sig_ds['lon<deg>']*u.deg.to(u.rad))
y = sig_ds['d1<kpc>']*np.cos(sig_ds['lon<deg>']*u.deg.to(u.rad))-8


plt.figure(figsize=(10,10))

plt.scatter(sig_ds['lon<deg>'],sig_ds['d1<kpc>'],c=sig_ds['vel<km/s>'],s=5,cmap='YlGnBu',norm='asinh')


plt.figure(figsize=(10,10))
plt.scatter(x,y+8,c=sig_ds['vel<km/s>'],s=5,cmap='YlGnBu',norm='asinh')
plt.scatter(0,0,marker='x',color='red',s=20)
plt.scatter(0,8,marker='x',color='red',s=20)

plt.xlim(15,-15)
plt.ylim(-10,25)

plt.show()


def galactocentric_radius(l,b,v):
    Ro = 8.5*u.kpc
    Vo = 229*u.km/u.s
    dvdR = -1.7*u.km/u.s/u.kpc
    
    
    t1 = Ro*np.sin(l)
    t2 = Vo - dvdR*Ro
    t3 = v/np.cos(b) + Vo*np.sin(l) - dvdR*Ro*np.sin(l)
    
    return (t1*t2)/t3


Rgs = galactocentric_radius(np.array(sig_ds['lon<deg>'])*u.deg,\
                      np.array(sig_ds['lat<deg>'])*u.deg,\
                      np.array(sig_ds['vel<km/s>'])*u.km/u.s)


sigs_w_nice_rs = [(Rgs[i].value,sig_ds['sig<solMass/pc2>'][i])
                  for i in range(120416) if Rgs[i]>0 and Rgs[i].value<20]

sigs_w_nice_rs = np.array(sigs_w_nice_rs).transpose()


plt.hist(sigs_w_nice_rs[0],bins=50)
plt.show()


plt.scatter(sigs_w_nice_rs[0],sigs_w_nice_rs[1],s=0.5)
plt.xlim(1,20)
plt.loglog()

# plt.scatter(sig_ds['lon<deg>'],Rgs)
# plt.vlines([-90,90],-40,40,color='red',linestyle='--')
# plt.hlines(8,180,-180,color='red',linestyle='--')
# plt.ylim(-10,10)


plt.scatter(lons,10*np.nanmedian(whi.to(u.solMass/u.pc**2),axis=0),color='orange',s=5)
plt.scatter(lons,np.nanmax(whi.to(u.solMass/u.pc**2),axis=0),color='red',alpha=0.5,s=5)
plt.scatter(lons,10*np.nanmean(whi.to(u.solMass/u.pc**2),axis=0),color='blue',alpha=0.5,s=5)


plt.plot(np.abs(lons),np.nanmedian(whi.to(u.solMass/u.pc**2),axis=0),color='orange',linewidth=0.5)
plt.plot(np.abs(lons),np.nanstd(whi.to(u.solMass/u.pc**2),axis=0),color='red',linewidth=0.5)
plt.xlim(0,190)


#from McKee & Ostriker -- pressure equilibrium of CNM & WNM
def Pth_of_R (R):
    return 1.4e4*np.e**(-R/(5.5*u.kpc))*u.K*u.cm**-3


plt.plot(np.arange(0,20,0.1),Pth_of_R(u.kpc*np.arange(0,20,0.1)))
plt.loglog()
plt.xlabel('galactocentric radius [kpc]')
plt.ylabel(r'thermal pressure [K $cm^{-3}$]')
plt.show()


#when loglog, this looks maybe like the transition from laminar to turbulent flow??
#look into calculating reynolds number etc
#also flaring at higher l -> gas is less gravitationally bound? more turbulent? less pressure?
plt.scatter(wholels,whole_sig,s=1e-3,color='black')
# plt.plot(lons,1e7*np.nanmax(whi,axis=0),color='red',alpha=0.5)
# plt.plot(lons,1e6*np.nanmin(whi,axis=0),color='red',alpha=0.5)

plt.plot(lons,1e7*np.nanmedian(whi,axis=0),color='orange')
plt.plot(lons,1e7*np.nanstd(whi,axis=0),color='red')


plt.xlabel('galactic longitude [deg]')
plt.ylabel(r'surface density [$M_\odot pc^{-2}$]')
# plt.loglog()
plt.show()


plt.scatter(wholebs,whole_sig,s=1e-5,alpha=0.8,color='black')

plt.plot(lats,1e7*np.nanmedian(whi,axis=1),color='orange')

#log normal with powerlaw tail???

plt.xlabel('galactic latitude [deg]')
plt.ylabel(r'surface density [$M_\odot pc^{-2}$]')
plt.loglog()
plt.show()


def calculate_dvtv(array2d,dpx):
    dvtv = np.zeros_like(array2d)
    for i in range(1,len(array2d)-1):
        row = array2d[i]
        for j in range(1,len(row)-1):
            thispx = row[j]
            
            above = array2d[i+1][j]
            below = array2d[i-1][j]
            right = row[j+1]
            left = row[j-1]
            tl = array2d[i+1][j-1]
            tr = array2d[i+1][j+1]
            bl = array2d[i-1][j-1]
            br = array2d[i-1][j+1]
            
            av_ar = np.nanmean((above,below,right,left,tl,tr,bl,br))
            
            dvtv[i][j] = np.abs(thispx-av_ar)/dpx
            
            
    return dvtv
    

whisolm = whi.to(u.solMass/u.pc**2)


whisolm_dvtv = calculate_dvtv(whisolm.value,0.125)


plt.figure(figsize=(20,9))
# plt.imshow(whi.to(u.solMass/u.pc**2).value,extent=[360,-360,-50,50],norm='asinh')
# plt.contourf(lons,lats,whi.to(u.solMass/u.pc**2).value,\
#              levels=[1.25,2.5,5,10,15,20,25,50,100,200],cmap='viridis',norm='asinh')
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(300,60)
# plt.colorbar()
cntrs = plt.contour(lons,lats,whi.to(u.solMass/u.pc**2).value,levels=np.arange(0.5,100,0.5),\
             colors='black',norm='asinh',linewidths=0.3,alpha=0.5)


plt.figure(figsize=(20,9))
plt.imshow(whisolm_dvtv,extent=[360,-360,-50,50],vmax=2)
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(300,60)


plt.figure(figsize=(20,9))
# plt.imshow(whi.to(u.solMass/u.pc**2).value,extent=[360,-360,-50,50],norm='asinh')
plt.contourf(lons,lats,whi.to(u.solMass/u.pc**2).value,\
             levels=[1.25,2.5,5,10,15,20,25,50,100,200],cmap='viridis',norm='asinh')
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(270,90)
plt.colorbar()
plt.contour(lons,lats,whi.to(u.solMass/u.pc**2).value,levels=[10,15],\
             colors='white',norm='asinh',linewidths=0.5)


plt.figure(figsize=(20,9))
# plt.imshow(whi.to(u.solMass/u.pc**2).value,extent=[360,-360,-50,50],norm='asinh')
plt.contourf(lons,lats,whi.to(u.solMass/u.pc**2).value,levels=[1.25,2.5,5,10,15,20,25,50,100,200],\
             cmap='viridis',norm='asinh')
plt.colorbar()
plt.contour(lons,lats,whi.to(u.solMass/u.pc**2).value,levels=[10,15],\
             colors='white',norm='asinh',linewidths=0.5)
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(90,-90)


plt.scatter(np.abs(wholels),whole_sig,s=1e-5,alpha=0.8,color='black')
plt.scatter(lons,np.nanmean(whi.to(u.solMass/u.pc**2),axis=0),color='blue',alpha=0.5,s=5)
plt.scatter(-lons,np.nanmean(whi.to(u.solMass/u.pc**2),axis=0),color='orange',alpha=0.5,s=5)


plt.xlabel('|galactic longitude| [deg]')
plt.ylabel(r'surface density [$M_\odot pc^{-2}$]')
plt.yscale('log')
plt.xlim(0,180)
plt.show()


import fun as fns


gal_mom1_hi = fns.get_mom1(hi4pi.data,[0,155])
gal_mom1_co = fns.get_mom1(dht.data,[0,155])


wco = (np.nansum(dht.data,axis=0)*u.K*1.3*u.km/u.s*co_conv).to(u.solMass/u.pc**2)


plt.figure(figsize=(20,9))
plt.imshow(gal_mom1_co,extent=[360,-360,-50,50],cmap='YlGnBu')
plt.colorbar()

plt.contour(lons,lats,wco.value,levels=[1.25,10],linewidths=0.4,colors='black')
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(180,-180)


plt.figure(figsize=(20,9))
plt.imshow(gal_mom1_hi,extent=[360,-360,-50,50],cmap='YlGnBu')
plt.colorbar()
plt.contour(lons,lats,whisolm.value,levels=[1.25,5,10,50,100],linewidths=0.4,colors='black')

plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(180,-180)


whi_neg = np.nansum(hi4pi.data[:77,:,:],axis=0)


plt.figure(figsize=(20,9))
plt.contourf(lons,lats,gal_mom1_hi,levels=np.arange(-100,0,3),\
             cmap='viridis',norm='asinh')
plt.contour(lons,lats,whi_neg,levels=[200,400,600,800],\
             cmap='rainbow',linewidths=0.7)
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(180,-180)


whi_pos = np.nansum(hi4pi.data[77:,:,:],axis=0)


plt.figure(figsize=(20,9))
plt.contourf(lons,lats,gal_mom1_hi,levels=np.arange(0,100,3),\
             cmap='viridis',norm='asinh')
plt.contour(lons,lats,whi_pos,levels=[200,400,600,800],\
             cmap='rainbow',linewidths=0.7)
plt.vlines([-90,90,270,-270],-50,50,color='black',linestyle='--')
plt.xlim(180,-180)


# okay to help me see how this pattern does / doesnt correspond to just the regular rotation curve, 
# what i want to do is make a model / say heres a grid of points (or sets of rings) that are rotating at a certain
# circular velocity around the galactic center, at a specific galactocentric radius
# what is the transformation for getting the projected line of sight velocity for these rings? (tbh
# eqs 1 in levine blitz and heiles are probably it, but i want to figure it out too to
# convince myself that thats whats happening )
# 
# then after this, use their eq 1 to reproject my data onto a galactocentric cylinder -> will make it 
# easier to then simulate a 'face on' view (by integrating down z) and find high surface density hi that way?


plt.imshow(np.nansum(hi4pi.data,axis=1),extent=[360,-360,-100,100])
plt.xlim(180,-180)


def RzPhiD_outer(l,b,v):
    #using eqs1 in levine/blitz/heiles, + quadratic for d
    #for now just solving the outer galaxy (b+sqrt(stuff))
    SL = np.sin(l)
    CL = np.cos(l)
    SB = np.sin(b)
    CB = np.cos(b)
    
    Ro = 8.5*u.kpc
    Vo = 220*u.km/u.s
    
    R = (Ro*Vo*SL*CB) / (v + Vo*SL*CB)  
    
    d = (1/CB)*(Ro*CL + np.sqrt(R**2 - (SL*Ro)**2))
    
    z = d*SL
    
    phi = np.arcsin(SL*CB*d/R)
    
    
    return R,z,phi


#this takes too long :( find other way to calculate
#-- should be parallelizable since each lbv coordinate is independent of the others when being sent to Rphiz )
maybe=[RzPhiD_outer(l,b,v)
 for l in np.arange(180,-180,-0.125)*u.deg 
 for b in np.arange(-50,50,0.125)*u.deg 
 for v in np.arange(-100,100,1.3)*u.km/u.s]


#alsooo maybe a statistical way of thinking about distances--imagine the milkyway as a circular disk with R=20kpc
#(actually do we know if it is a circle vs major axis / minor axis?)
#then populate that evenly with dots, then calculate some sort of stat correlating Rgal with l??
#kind of like -- maximizing arctan(Rg\Ro) gives |l| = ~68 -> ie its only inside these angles that you're
#getting everything within a specific ring of Rg=something, outside of l>68, youre only getting parts of
#galactocentric annuli, not the whole thing (ie there is a far part of the Rg=something annulus that is at
#smaller l and a near side that is at larger l)


#look at the places where the contour lines bunch up -> ie where the changes in surface density are happening
#quickly -> ie where the local derivative of dsig/dx is steep
#maybe theres a way to calculate this pixel by pixel? then plotting that might reveal "boundaries"??


# from convo with ruth
# 
# vl is prop to k to some power and that power tells you which scale eddies dominate
# to get ^, figure out what is steering the motions
# 
# 
# make a plot of heliocentric distance vs velocity 
# (where the velocity is made up of rotation curve, bulk motion and turbulent motion)
#                                                  
# reynolds number - macro velocity x macro lengthscale / thermal velocity x mean free path
