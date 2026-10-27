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


hi4pi,dht = fits.open('both_same_grid_720d.fits')


cat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


dnoisd = dht.data*(dht.data>3e-2)


cos = [dnoisd[cat['vo'][i]:cat['vi'][i],
         cat['boc'][i]:cat['bic'][i],
         cat['loc'][i]:cat['lic'][i]] for i in range(22)]


c0s = [np.nansum(co,axis=0)*1.3*1.97e20 for co in cos]


[(
 plt.contourf(np.log10(c0s[i]),levels=np.arange(20,22,0.25)),
 plt.contour(np.log10(c0s[i]),levels=[20.5],colors='red'),
 plt.title(cat['name'][i]),
 plt.show()
    
    
    ) for i in range(22)
]


c1s = [fns.get_mom1(cos[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100) for i in range(22)]


c2s = [fns.get_momN(cos[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2) for i in range(22)]


masses_gmc = [fns.mass_n_area(cos[i]*u.K,cat['distance<pc>'][i]*u.pc,'CO')[0].to(u.solMass).value for i in range(22)]


areas_gmc = [fns.mass_n_area(cos[i]*u.K,cat['distance<pc>'][i]*u.pc,'CO')[1].value for i in range(22)]


#i will be using the geometric mean bc angular momentum does care about the distribution of mass
#and doing sqrt(A/pi) affects the calculatd moment of inertia
radii_gmc = np.sqrt(np.nanstd(cat['lon_GMC<deg>'],axis=1)*\
                    np.nanstd(cat['lat_GMC<deg>'],axis=1))*u.deg.to(u.rad)*cat['distance<pc>']


#decide that my hi geometry will be spherical, & explicitly discuss why the accumulation volume 
#being treated as a cylinder makes sense for observations of face-on system (where the height of the
#cylinder is perpendicular to & going through the plane) but because we are dealing with an edge-on view,
#this would mean the area/radius part of cylinder would be along the radial axis, where we dont have spatial
#information, and the height (which for us is projected on to lb) is hand waved to be "arbitrarily/infinitely long" 

#so start here
4/3 pi r3 x rho = mass
#then what we observe of this is a circle, bc we see all the stuff collapsed/integrated along the radial axis 
Sigma = integral(rho x dr)
#which we are assumimg has a constant volume density 
then the integral just evaluates to 2R, bc its going through the diameter


#actually maybe we keep cylinder but use what we know about hi filaments to decide the shape/orientation
#okay how about saying, we assume it is cylinders that we are seeing, and it is still
#a cylinder perpendicular to the plane, and the cylinder's height is 2xradius, so in projection, this 
#gives us a square that is 2R by 2R
#and for this shape, if we integrate radially , it gives us 
Sigma = rho x 2R so (sigma/2R x 4R^2) = rho so (sigma x 2)x4piR^2x2 = M
#blah blah but yeah this is the way to think of it

#okay for the discussion, just say here are the assumptions associated with this method
#and thats it

#so for me now i just gotta derive the expression for R and move on (edited fun.py)


#also one caveat is that we are forcing the hi to be circle/cylindrical 
#-> also affects angular momentum via moment of inertia?


QTable(fits.open('post_flash_cat.fits')[1].data)


ras = [fns.accumulation_radius(masses_gmc[i]*u.solMass,cat['distance<pc>'][i]*u.pc)[1].value for i in range(22)]


ras_deg = [fns.accumulation_radius(masses_gmc[i]*u.solMass,cat['distance<pc>'][i]*u.pc)[0].value for i in range(22)]


cat['lon_HI<deg>'] = np.around(np.array((np.nanmean(cat['lon_GMC<deg>'],axis=1)+ras_deg,
                    np.nanmean(cat['lon_GMC<deg>'],axis=1)-ras_deg)).transpose(),2)

cat['lat_HI<deg>'] = np.around(np.array((np.nanmean(cat['lat_GMC<deg>'],axis=1)-ras_deg,
                    np.nanmean(cat['lat_GMC<deg>'],axis=1)+ras_deg)).transpose(),2)


loh,lih = ((cat['lon_HI<deg>']-360)/-0.125).astype('int').transpose()
boh,bih = ((cat['lat_HI<deg>']+50)/0.125).astype('int').transpose()


his = [hi4pi.data[:,boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]


cat['vel<km/s>'][20][1]= 28


[(plt.plot(np.arange(155)*1.3-100,np.nanmean(his[i],axis=(1,2))),
  plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,100*np.nanmean(cos[i],axis=(1,2))),
  plt.vlines(cat['vel<km/s>'][i],0,80,color='red'),
plt.show())
  for i in range(22) ]


cat['vi'][20] = int((28+100)/1.3)


his = [hi4pi.data[cat['vo'][i]:cat['vi'][i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]


h1s = [fns.get_mom1(his[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100) for i in range(22)]


h2s = [fns.get_momN(his[i],np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2) for i in range(22)]


h2s = [np.sqrt(h2) for h2 in h2s]


c2s = [np.sqrt(c2) for c2 in c2s]


#okay, again just say well the turbulent linewidth can be found in these different ways
#and the one that we will choose is the one that is the second moment of the spectrum
#because this catches all the emission, vs the first & second moments what they do is break the emission
#up into small regions, then say what is the turbulence in those small regions
#so this is helpful maybe when thinking about the power spectrum? as in the relationship between linewidths at 
#small and large scales within the same cloud maybe says something about whether its a true cascade, or theres
#a dominant sizescale for the cloud
#!! and  this is what the patch by patch spectra will do as well !! experiment with different aperture sizes 
#as a way of varying the size scale, then see what those profiles are for the clouds (linewidth vs size scale)


[(
plt.hist(np.ravel(c1s[i]),bins=100,
#          bins=np.arange(cat['vel<km/s>'][i][0],
#                         cat['vel<km/s>'][i][1],0.65/2),
        alpha=0.5),
plt.hist(np.ravel(h1s[i]),bins=100,
#          bins=np.arange(cat['vel<km/s>'][i][0],
#                         cat['vel<km/s>'][i][1],0.65/2),
        alpha=0.5),
    
    plt.yscale('log'),
    plt.title(cat['name'][i]),
plt.show()
  
)for i in range(22)]


vel_disper_12 = np.sqrt(np.array([fns.get_momN(np.histogram(c1s[i],bins=np.arange(cat['vel<km/s>'][i][0],
                                            cat['vel<km/s>'][i][1],0.65))[0],
            np.arange(cat['vel<km/s>'][i][0],
            cat['vel<km/s>'][i][1],0.65)[:-1],2)

            for i in range(22)]))


old_sigmas = QTable(fits.open('post_flash_cat.fits')[1].data)['vel_dispersion_GMC<km/s>']


[co.shape for co in cos]


vel_disper_2spec = np.sqrt(np.array([fns.get_momN(np.nanmean(cos[i],axis=(1,2)),
            np.arange(cat['vel<km/s>'][i][0],
            cat['vel<km/s>'][i][1],1.3),2)

for i in range(22)]))


hivel_disper_2spec = np.sqrt(np.array([fns.get_momN(np.nanmean(his[i],axis=(1,2)),
            np.arange(cat['vel<km/s>'][i][0],
            cat['vel<km/s>'][i][1],1.3),2)

for i in range(22)]))


[(plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,np.nanmean(cos[i],axis=(1,2))),
  plt.title(str(round(vel_disper_2spec[i],2))+' '+\
            str(round(vel_disper_12[i],2))+' '+\
            str(round(old_sigmas[i],2))),
  
  
#   plt.vlines(cat['vel<km/s>'][i],0,80,color='red'),
plt.show())
  for i in range(22) ]


#hmmm i think what the old sigmas (average of the second moment) tell us is actually something closer
#to the thermal broadening, because it is the linewidth of each individual location
#while the turbulent linewidth is necessarily set by how many different values of v have some part of the cloud
#peaking at those velocities
#alternatively, the individual linewidths say something about local teeny turbulence, at the size scale of however
#many pc 0.125 degrees translates to
#check if this (the second moment of the first moment population) is similar to the second moment of the spectrum

#and pay attention the clouds that have more small scale turbulence than whole cloud 


#okay now to looking at om a pixel by pixel level, whether the emissions line up for co & hi
#also notable that the hi emission is on average to the negative side, except for cma(hi in center),
#cra(overlap), gemini(hi to right),her(hi to right),mad(hi in center),mob(h1 in center)
#monr2(hi bigger),oph(hi in center),ora (separateish hi emission at tails),orb(same as a,also hi bigger),
#orl (hi in center),pers(in the center),polaris(in the center),ros(center)
#starging to wonder if the average behavior is hi on the left, or just hi narrower

edges = np.array(([0,2.5,6],[-20,1,10],[-4,-1,3],[-3,3,6],[8,14,19],[4,6,9],[7,14,17],
[-5,1,15],[4,11,20],[-4,0,6],[-5,8,13],
[20,26,30],[-5,5,8],[5,14,16],[-2,5,9],
[0,10,20],[0,10,20],[-20,2,20],[-3,8,12],
[-8,1,6],[5,11,20],[-3,8,18]))


[(
plt.hist(np.ravel(c1s[i]),bins=100,
#          bins=np.arange(cat['vel<km/s>'][i][0],
#                         cat['vel<km/s>'][i][1],0.65/2),
        alpha=0.5),
plt.hist(np.ravel(h1s[i]),bins=100,
#          bins=np.arange(cat['vel<km/s>'][i][0],
#                         cat['vel<km/s>'][i][1],0.65/2),
        alpha=0.5),
    
    plt.yscale('log'),
    plt.title(cat['name'][i]),
plt.show()
  
)for i in range(22)]


i=0
[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(0,3,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(0,3,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 plt.contourf(c1s[i],levels=np.arange(3,6,0.65),cmap='Purples',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),

 
 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [0]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(-10,2,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
 plt.contourf(h1s[i],levels=np.arange(-20,-10,0.65),cmap='Reds',
          extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                 cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(-10,2,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 plt.contourf(c1s[i],levels=np.arange(2,10,0.65),cmap='Purples',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),

 
 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [1]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(-4,0.5,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(-4,0,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 plt.contourf(c1s[i],levels=np.arange(-1,3,0.65),cmap='Purples',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),

 
 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [2]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(-3,1,0.65),cmap='Reds',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
     plt.contourf(h1s[i],levels=np.arange(1,3,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(1,3,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 plt.contourf(c1s[i],levels=np.arange(3,6,0.65),cmap='Purples',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),

 
 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [3]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(14,19,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(14,19,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 plt.contourf(c1s[i],levels=np.arange(8,14,0.65),cmap='Reds',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    plt.contourf(c1s[i],levels=np.arange(19,22,0.65),cmap='Purples',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),

 
 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [4]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(5,7,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.contourf(h1s[i],levels=np.arange(4,5,0.65),cmap='Reds',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    plt.contourf(h1s[i],levels=np.arange(7,8,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(3,8,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [5]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(10,14,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.contourf(h1s[i],levels=np.arange(8,10,0.65),cmap='Reds',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
   
    
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(10,14,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(14,16,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [6]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(0,11,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.contourf(h1s[i],levels=np.arange(11,15,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
   
    
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(-5,0,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(0,11,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [7]
]


[(
  
 plt.contourf(h1s[i],levels=np.arange(10,15,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.contourf(h1s[i],levels=np.arange(15,20,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
   
    
 plt.show(),
    
    
    plt.contourf(h1s[i],levels=np.arange(10,14,0.65),cmap='Reds',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.contourf(h1s[i],levels=np.arange(14,17,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.contourf(h1s[i],levels=np.arange(17,20,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    plt.show(),
    
 plt.contourf(c1s[i],levels=np.arange(5,10,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(10,15,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 plt.show(),
    
    plt.contourf(c1s[i],levels=np.arange(5,7,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    plt.contourf(c1s[i],levels=np.arange(7,15,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    

 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [8]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(-4,0,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(-4,10,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(0,5,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [9]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(-2,8,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
 plt.show(),
 plt.contourf(c1s[i],levels=np.arange(-2,8,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(8,12,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [10]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(26,30,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
 plt.show(),
    
    
 plt.contourf(c1s[i],levels=np.arange(20,26,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(26,30,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

    plt.contourf(c1s[i],levels=np.arange(30,36,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

    plt.show(),
    
    plt.contourf(c1s[i],levels=np.arange(20,28,0.65),cmap='Reds',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    plt.contourf(c1s[i],levels=np.arange(28,36,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
 plt.show(),
 
 plt.title(cat['name'][i]),
#  plt.colorbar(),
 plt.show()
    
    
    ) for i in [11]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(4.5,8,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
    
 plt.show(),
    
    
 plt.contourf(c1s[i],levels=np.arange(-5,5,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(5,8,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

    plt.contourf(c1s[i],levels=np.arange(8,15,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 
    
 plt.title(cat['name'][i]),
#  plt.colorbar(),
    
    plt.show()
    ) for i in [12]
]


[(
    
  
 plt.contourf(h1s[i],levels=np.arange(6,14,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
    
      
 plt.show(),
    
      
 plt.contourf(c1s[i],levels=np.arange(6,14,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(14,16,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

#     plt.contourf(c1s[i],levels=np.arange(8,15,0.65),cmap='Oranges',
#               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                      cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 
    
 plt.title(cat['name'][i]),
#  plt.colorbar(),
    
    plt.show()
    ) for i in [13]
]


[(
plt.hist(np.ravel(c1s[i]),bins=100,
#          bins=np.arange(cat['vel<km/s>'][i][0],
#                         cat['vel<km/s>'][i][1],0.65/2),
        alpha=0.5),
plt.hist(np.ravel(h1s[i]),bins=100,
#          bins=np.arange(cat['vel<km/s>'][i][0],
#                         cat['vel<km/s>'][i][1],0.65/2),
        alpha=0.5),
    plt.yscale('log'),
    plt.title(cat['name'][i]),
plt.show()
  
)for i in range(22)]


#hmmm also a lot of the overlap velocities look like noise on the co? look at when i come back
#also land on sth with sigma clipping dht or no


c1sold = fits.open('../not_rn_7/mom1_co.fits')


[(
    
    plt.contourf(h1s[i],levels=np.arange(-25,0,0.65),cmap='Reds',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 
  
 plt.contourf(h1s[i],levels=np.arange(2,10,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 
    plt.contourf(h1s[i],levels=np.arange(50,60,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 
      
 plt.show(),
    
      
 plt.contourf(c1s[i],levels=np.arange(-2,2,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(2,10,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

    plt.contourf(c1s[i],levels=np.arange(10,20,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 
    
 plt.title(cat['name'][i]),
#  plt.colorbar(),
    
    plt.show()
    ) for i in [15]
]


[(
  
 plt.contourf(h1s[i],levels=np.arange(0,15,0.65),cmap='Oranges',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 
    plt.contourf(h1s[i],levels=np.arange(15,20,0.65),cmap='Purples',
              extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]]),
 
      
 plt.show(),
    
      
 plt.contourf(c1s[i],levels=np.arange(-2,2,0.65),cmap='Reds',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    plt.contourf(c1s[i],levels=np.arange(2,10,0.65),cmap='Oranges',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 

    plt.contourf(c1s[i],levels=np.arange(10,20,0.65),cmap='Purples',
              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                     cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
 
    
 plt.title(cat['name'][i]),
#  plt.colorbar(),
    
    plt.show()
    ) for i in [16]
]


#okay just put things in the cat + save subcubes/maps and move on
#right now looking at the first moment histograms  isnt going anywhere


[(plt.contourf(np.log10(c0s[i]),levels=np.arange(20,22,0.25)),
plt.show())
              
   for i in range(22)           
              
              ]


[(plt.contourf(c1s[i],levels=10),
plt.show(),
  plt.contourf(c1sold[i].data,levels=10),
plt.show()
 
 )
 
              
   for i in range(22)           
              
              ]


[(plt.contourf(np.sqrt(c2s[i]),levels=10),
plt.show(),
 
 )
 
              
   for i in range(22)           
              
              ]


masses_hi = [fns.mass_n_area(his[i]*u.K,
                             cat['distance<pc>'][i]*u.pc,
                             '21cm')[0].to(u.solMass).value for i in range(22)]


cat['mass_GMC<Msun>'] = masses_gmc
cat['mass_HI<Msun>'] = masses_hi


cat['radius_GMC<pc>'] = radii_gmc
cat['radius_HI<pc>'] = ras


cat['vel_dispersion_GMC<km/s>'] = vel_disper_2spec
cat['vel_dispersion_HI<km/s>'] = hivel_disper_2spec


#also for the isolated clouds its worth asking/thinking about what the lowest contour is, like eg for 
#chamaeleon its p close to 1e20 bc not much changes in terms of the how into the cloud the contour lines were
#whenever i was switching around thresholds (1e-2,5e-2,etc for dnoisd) but aquila, cam, etc seem to have a lot of
#low intensity CO emission / low column density gas that is CO bright
#cirrus clouds ??


def make_fits_from_list(array,filename):
    hdulist = fits.HDUList([fits.PrimaryHDU(array[0])]+[fits.ImageHDU(array[i]) for i in range(1,22)])
    hdulist.writeto(filename,overwrite=True)


make_fits_from_list(his,'hi_clouds.fits')


make_fits_from_list(cos,'h2_clouds.fits')


h0s = [np.nansum(hi,axis=0)*1.3*1.82e18 for hi in his]


make_fits_from_list(c0s,'h2_mom0.fits')
make_fits_from_list(c1s,'h2_mom1.fits')
make_fits_from_list(c2s,'h2_mom2.fits')


make_fits_from_list(h0s,'hi_mom0.fits')
make_fits_from_list(h1s,'hi_mom1.fits')
make_fits_from_list(h2s,'hi_mom2.fits')


h0lvs = [np.nansum(hi,axis=1)*0.125*1.82e18 for hi in his]
h0bvs = [np.nansum(hi,axis=2)*0.125*0.5*1.82e18 for hi in his]


c0lvs = [np.nansum(co,axis=1)*0.125*1.97e20 for co in cos]
c0bvs = [np.nansum(co,axis=2)*0.125*0.5*1.97e20 for co in cos]


make_fits_from_list(h0lvs,'hi_mom0lv.fits')
make_fits_from_list(h0bvs,'hi_mom0bv.fits')


make_fits_from_list(c0lvs,'h2_mom0lv.fits')
make_fits_from_list(c0bvs,'h2_mom0bv.fits')


cat.write('lastlastLAST_cat.fits',overwrite=True)
