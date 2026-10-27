#!/usr/bin/env python
# coding: utf-8


use higher order moment maps to determine if something is a gaussian

do stats on how many diff types of spectra are there (in general, per cloud) -> nice gaussians, two peaks etc

go ahead and do the weighting by second moment & do the velocity field maps

plot for second moment varying as a function of distance from the center
    or histograms of (velocity dispersion of hi within co contours) & 
    (velocity dispersion of hi outside co contours) <- compare these and maybe the place where the
    characteristic location/distance where the linewidth changes is the h1->h2 transition
    (also: set of spectra within, set of spectra outside)
    
    
for the double peaked h2 histogram -> look at the spatial distribution of linewidths / look in more detail
before concluding that these are cnm + wnm


calculate velocity field maps from the lab data, compare to nia's results 
-> compare the two results to say how a higher resolution dataset makes a difference -> define accumulation
radius using masses from nia (try to match iniitial assumptions as much as possible)
-> use the same central velocity & velocity channel width
(look for 13CO data from mark heyer)


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

import fun as fns
from astropy.modeling import models,fitting


hi4pi,dht = fits.open('both_same_grid_720d.fits')


hi4pi35 = (hi4pi.data *( hi4pi.data>0.43*5))[50:105,:,:]
dht35 = (dht.data*(dht.data>0.2*3))[50:105,:,:]


h0 = np.nansum(hi4pi35,0)*1.3*1.82e18
c0 = np.nansum(dht35,0)*1.3*1.97e20


h1 = fns.get_mom1(hi4pi35*(hi4pi35>0.043*5),np.arange(55)*1.3-35)
c1 = fns.get_mom1(dht35  *(dht35>0.2*3),np.arange(55)*1.3-35)


h2 = fns.get_momN(hi4pi35*(hi4pi35>0.043*5),np.arange(55)*1.3-35,2)
c2 = fns.get_momN(dht35  *(dht35>0.2*3),np.arange(55)*1.3-35,2)


# plt.hist(np.ravel(np.sqrt(c2)),bins=100,alpha=0.5)
plt.hist(np.ravel(np.sqrt(h2)),bins=100,alpha=0.5)
plt.plot(np.arange(0,15,0.5),models.Gaussian1D(17e4,6,1.5)(np.arange(0,15,0.5)))
plt.plot(np.arange(0,15,0.5),models.Gaussian1D(11e4,10,1.1)(np.arange(0,15,0.5)))
plt.plot(np.arange(0,15,0.5),models.Gaussian1D(11e4,10,1.1)(np.arange(0,15,0.5))+\
        models.Gaussian1D(17e4,6,1.5)(np.arange(0,15,0.5)))

# plt.plot(np.arange(0,15,0.5),models.Gaussian1D(11e4,11.4,1.1)(np.arange(0,15,0.5))+\
#         models.Gaussian1D(17e4,7.1,1.5)(np.arange(0,15,0.5)))

# plt.loglog()
# plt.yscale('log')
plt.show()


amplitude = np.nanmax(hi4pi35*(hi4pi35>0.043*5),axis=0)
mean = h1
stddv = np.sqrt(h2)


ams = np.ravel(amplitude)
mns = np.ravel(mean)
sds = np.ravel(stddv)


gaus_spec = [model_specs[i](np.arange(-35,36,1.3)) for i in range(len(ams))]


gausmodel_hi4pi = np.array(gaus_spec)


ppvgausmodelh = gausmodel_hi4pi.reshape(800,5760,55).transpose().swapaxes(1,2)


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
plt.show()) for i in range(0,len(gaus_spec))
 if np.nanstd(gaus_spec[i])>20]


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
plt.show()) 
 for i in range(0,len(gaus_spec))
 if 1.5*np.nansum(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]) <
np.nansum(gaus_spec[i])]


gaus0 = np.nansum(ppvgausmodelh,axis=0)


plt.figure(figsize=(20,10))

plt.contourf((gaus0*1.3*1.82e18)/h0,levels=[0.,1.5,10],
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

# plt.contour((gaus0*1.3*1.82e18)/h0,levels=[1],colors='black',
#           origin='lower',extent=[360,-360,-50,50])
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.scatter(np.ravel((gaus0*1.3*1.82e18)/h0),np.ravel(h0),s=1e-2)
plt.loglog()
plt.ylim(1e18,)


plt.scatter(np.ravel((gaus0*1.3*1.82e18)/h0),np.ravel(h2),s=1e-2)
plt.loglog()
plt.ylim(1e-2,)


plt.scatter(np.ravel(h0),np.ravel(h2),s=1e-2,alpha=0.6)
plt.loglog()
plt.ylim(1e-2,)


plt.figure(figsize=(20,10))

plt.contourf(np.log10(h0),levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contour((h0),levels=[4e20],colors='black',
          origin='lower',extent=[360,-360,-50,50],linewidths=0.8)
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(np.log10(gaus0*1.3*1.82e18),levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contour((gaus0*1.3*1.82e18),levels=[4e20],colors='black',linestyles='--',
          origin='lower',extent=[360,-360,-50,50],linewidths=0.8)
plt.contour(h0,levels=[4e20],colors='blue',linestyles='--',
          origin='lower',extent=[360,-360,-50,50],linewidths=0.8)
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(stddv,levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='grey',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(np.nanstd(ppvgausmodelh,axis=0),levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')


plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.hist(np.ravel(mean/(np.argmax(ppvgausmodelh,0)*1.3-35)),bins=100)
plt.yscale('log')
plt.show()


plt.hist(np.ravel(stddv/(np.std(ppvgausmodelh,0)*1.3-35)),bins=100)
plt.yscale('log')
# plt.xlim(-50,50)
plt.show()


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
  plt.plot(gaus_spec[i]/hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
 plt.show())
for i in range(0,len(ams),10000)
if np.ravel(mean/(np.argmax(ppvgausmodelh,0)*1.3-35))[i] < 1.2]


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
  plt.plot(gaus_spec[i]/hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
 plt.show())
for i in range(0,len(ams),1000)
if np.ravel(mean/(np.argmax(ppvgausmodelh,0)*1.3-35))[i] < 1.2]


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
#   plt.plot(gaus_spec[i]/hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
 plt.show())
for i in range(0,len(ams),1000)
if np.ravel(np.abs(mean/(np.argmax(ppvgausmodelh,0)*1.3-35)))[i]>1.2]


plt.figure(figsize=(20,10))

plt.contourf(np.abs(mean-(np.argmax(ppvgausmodelh,0)*1.3-35)),levels=np.arange(0,1,0.1),
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')


plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(np.argmax(ppvgausmodelh,0)*1.3-35,levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')


plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(np.nanstd(hi4pi35,axis=0)/stddv,levels=[0,1,2,5,10],
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')


plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.hist(np.ravel(np.sqrt(h2)),bins=100)
# plt.hist(np.ravel(stddv),bins=100)

plt.show()


plt.hist(np.ravel(stddv),bins=100)

plt.show()


plt.hist(np.ravel(np.nanstd(ppvgausmodelh,axis=0)),bins=100)

plt.show()


plt.hist(np.ravel(np.nanstd(hi4pi35,axis=0)),bins=100)

plt.show()


residuals = hi4pi35 - ppvgausmodelh


plt.imshow(np.nansum(residuals,0))


[(plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[i]),
 plt.plot(gaus_spec[i]),
 plt.show())
for i in range(0,len(gaus_spec)) if (np.ravel(stddv)[i]/np.ravel(np.nanstd(ppvgausmodelh,axis=0))[i]) > 2]


plt.figure(figsize=(20,10))

plt.contourf(),levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contour((gaus0*1.3*1.82e18),levels=[4e20],colors='black',linestyles='--',
          origin='lower',extent=[360,-360,-50,50],linewidths=0.8)
plt.contour(h0,levels=[4e20],colors='blue',linestyles='--',
          origin='lower',extent=[360,-360,-50,50],linewidths=0.8)
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='white',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)
