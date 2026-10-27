#!/usr/bin/env python
# coding: utf-8


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


plt.figure(figsize=(20,10))

plt.contourf(h0,levels=np.array([0.1,0.5,1,5,10])*1e21,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

# plt.contourf((gaus0*1.3*1.82e18)/h0,levels=[3,5],colors=['black','blue'],
#           origin='lower',extent=[360,-360,-50,50],linewidths=0.7)
# plt.contour(smol_neg,extent=[360,-360,-50,50],colors='white',levels=[10])
# plt.contour(smol_pos,extent=[360,-360,-50,50],colors='white',levels=[10])

# plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='red',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


h1 = fns.get_mom1(hi4pi35*(hi4pi35>0.043*5),[50,105])
c1 = fns.get_mom1(dht35  *(dht35>0.2*3),   [50,105])


h2 = fns.get_momN(hi4pi35*(hi4pi35>0.043*5),[50,105],2)
c2 = fns.get_momN(dht35  *(dht35>0.2*3),  [50,105],2)


amplitude = np.nanmax(hi4pi35*(hi4pi35>0.043*5),axis=0)
mean = h1
stddv = np.sqrt(h2)


ams = np.ravel(amplitude)
mns = np.ravel(mean)
sds = np.ravel(stddv)


np.ravel(np.nanstd(ppvgausmodelh,axis=0))[0]


np.nanstd(gaus_spec[0])


np.ravel(np.nanstd(hi4pi35,0))


sds


np.nanstd(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[0])


sds.shape


np.ravel(np.nanstd(gaus_spec,axis=1)).shape


plt.scatter(sds,np.ravel(np.nanstd(gaus_spec,axis=1)))


plt.plot(gaus_spec[0])
plt.plot(hi4pi35.transpose().swapaxes(0,1).reshape(-1,55)[0])


gaus_spec = [model_specs[i](np.arange(-35,36,1.3)) for i in range(len(ams))]


np.nanstd(model_specs[0](np.arange(-35,36.,1.3)))


np.nanmean(model_specs[0](np.arange(-35,36.,1.3)))


model_specs[0]


np.nanstd(gaus_spec[0])*2*np.sqrt(2*np.log(2))


2*np.sqrt(2*np.log(2))


plt.plot(gaus_spec[0])


gaus_spec


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
