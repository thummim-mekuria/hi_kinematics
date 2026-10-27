#!/usr/bin/env python
# coding: utf-8


from astrodendro import Dendrogram,periodic_neighbours
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy_healpix import HEALPix as aphp

import healpy as hp
from skimage.segmentation import find_boundaries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import bootleg as btlg
import datacube as dc


galds = pd.read_csv('whole_gal_dists.csv')


Rs = np.arange(0,15e3,1e2)*u.pc


phis = np.arange(0,360,1)*u.deg.to(u.rad)


zs = np.arange(-1e4,1e4,1e2)*u.pc


Ro = 8.5e3*u.pc
Theta_R = 220*u.km/u.s
phio = 180*u.deg.to(u.rad)


R,Phi,Z = np.meshgrid(Rs,phis,zs)


r_primes = np.sqrt(R**2 + Ro**2 + 2*R*Ro*np.cos(Phi))


ls = np.arccos((Ro+R*np.cos(Phi))/r_primes)


bs = np.arctan(Z/r_primes)


vs = np.sin(ls)*np.cos(bs)*(Ro/R - 1)*Theta_R


#okay i think this mostly worked except for l -> look at the condition when phi = 90, sin(l)
#should reduce to R/Ro but it goes to 0 instead


#then the next task is -- if i define a distribution profile in R,Phi,Z .. how do i convert it
#to a profile in l,b,v ?


R.shape


plt.imshow(np.nanmean(vs.value,axis=0))
plt.colorbar()


plt.imshow(np.nanmean(vs.value,axis=1))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=1))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=2))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=0))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=1))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=2))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=2))
plt.colorbar()


np.max(np.nanmean(ls.to(u.deg).value,axis=1))


vs.shape


plt.imshow(np.nanmean(vs.value,axis=2))
plt.colorbar()


plt.imshow(np.nanmean(ls.to(u.deg).value,axis=1))
plt.colorbar()


plt.imshow(np.nanmean(ls,axis=2)*u.rad.to(u.deg))
plt.colorbar()


plt.imshow(ls[:,:,0])


len(zs)
