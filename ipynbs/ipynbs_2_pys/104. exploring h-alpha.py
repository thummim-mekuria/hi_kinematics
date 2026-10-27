#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import QTable
import astropy.units as u
import astropy.constants as c

import datacube as dc
import fun as fns


zcat = QTable(fits.open('final_cat.fits')[1].data)


# yso counts from table 2 of https://arxiv.org/pdf/2206.05796 (gaia dr3)
names = ['aquilas, california, chamaeleon, cma, corona, crossbones,gemini,lupus,maddalena,monob1,monr2,ophiuchus,oriona,orionb,orionl,perseus,polaris,rosette,serpens,taurus,ursa']


names = names[0].split(',')


yso_counts = [62,770,394,408,28,176,64,1059,17,577,473,791,2072,854,
             1075,435,15,78,1136,432,12]


prs = np.array((names,yso_counts)).transpose()


zcat_ysos = [770,394,408,176,64,17,577,473,2072,854,
             1075,435,15,78,432]


len(yso_counts)


zcat


nv,nb,nl = 162,720,1440
dv,db,dl = -0.25,0.25,2.03
v0,b0,l0 = 720 ->0, 360 -> 0, 1-> -169,116


hav = (np.arange(162)-1)*2.03 - 169.116


hav[67]


hav35 = hav[67:102]


hab = (np.arange(720)-360)*0.25


hal = (np.arange(1440)-720)*-0.25


wham = fits.open('wham-ss-DR1-v161116-170912-grid.fits')[0].data


wham35 = wham[67:102,:,:]


hab50 = hab[160:560]


wham3550 = wham35[:,160:560,:]


wham3550720 = np.concatenate((wham35[:,160:560,720:],wham35[:,160:560,:],wham35[:,160:560,:720]),axis=2)


halb = np.nansum(wham3550,axis=0)


hi4pi,dht = fits.open('both_same_grid_720d.fits')


halb_long = np.nansum(wham3550720,axis=0)


h0 = np.nansum((hi4pi.data*(hi4pi.data>0.043*3))[50:105,:,:],axis=0)
c0 = np.nansum((dht.data*(dht.data>0.2*3))[50:105,:,:],axis=0)


wham3550720.shape


ha1 = fns.get_mom1(wham3550720,np.arange(35)*2.03 -35)


h1 = fns.get_mom1((hi4pi.data*(hi4pi.data>0.043*3))[50:105],np.arange(55)*1.3-35)


c1 = fns.get_mom1((dht.data*(dht.data>0.2*3))[50:105],np.arange(55)*1.3-35)


[(plt.plot(np.arange(35)*2.03-35,i),
plt.show())
for i in wham3550720.transpose().swapaxes(0,1).reshape(-1,35)[::1000]]


#look at laura lopez's work : star forming regions with multiple tracers
#look at what she does with h-alpha


plt.figure(figsize=(20,7))
plt.contourf(ha1,extent=[360,-360,-50,50],levels=np.arange(-25,25,2.6),cmap='Set2',alpha=0.4)

# plt.colorbar(orientation='horizontal')
plt.contourf(c1,extent=[360,-360,-50,50],levels=np.arange(-25,25,2.6),cmap='Set2')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],colors='black',linewidths=0.3)


plt.xlim(200,-200)
# plt.xlim(190,130)
# plt.ylim(-20,20)
plt.show()


plt.figure(figsize=(20,7))
plt.contourf(h1,extent=[360,-360,-50,50],levels=np.arange(-25,25,2.6),cmap='Set2',alpha=0.4)
# plt.contour(ha1,extent=[360,-360,-50,50],levels=np.arange(-35,35,5),cmap='Paired',linewidths=0.9)

# plt.colorbar(orientation='horizontal')
plt.contourf(c1,extent=[360,-360,-50,50],levels=np.arange(-25,25,2.6),cmap='Set2')
plt.colorbar()
plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],colors='black',linewidths=0.3)


plt.xlim(200,-200)
# plt.xlim(190,130)
# plt.ylim(-20,20)
plt.show()


plt.figure(figsize=(20,10))
plt.contourf(ha1,extent=[360,-360,-50,50],levels=np.arange(-35,35,5),cmap='Paired',alpha=0.7)
plt.colorbar(orientation='horizontal')
plt.contour(c1,extent=[360,-360,-50,50],levels=np.arange(-35,35,5),linewidths=1.2,cmap='Paired')
plt.colorbar()
# # plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],colors='white',linewidths=0.9,linestyles='--')


plt.xlim(180,-180)
# plt.xlim(-190,-140)
# plt.ylim(-20,20)
plt.show()


plt.figure(figsize=(20,10))
plt.contourf(ha1,extent=[360,-360,-50,50],origin='lower',levels=np.arange(-35,35,5),cmap='Paired',alpha=0.7)
plt.colorbar(orientation='horizontal')
plt.contour(h1,extent=[360,-360,-50,50],levels=np.arange(-35,35,5),linewidths=1.2,cmap='Paired')
plt.colorbar()
# # plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],colors='white',linewidths=0.9,linestyles='--')


plt.xlim(180,-180)
# plt.xlim(-190,-140)
# plt.ylim(-20,20)
plt.show()


plt.figure(figsize=(20,10))
plt.contourf(halb_long,extent=[360,-360,-50,50],origin='lower',levels=np.arange(0,70,5))
plt.colorbar(orientation='horizontal')
plt.contour(h0,extent=[360,-360,-50,50],levels=[400],colors='white',linewidths=0.9)
# plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],colors='white',linewidths=0.9,linestyles='--')


plt.xlim(250,-250)
plt.show()
