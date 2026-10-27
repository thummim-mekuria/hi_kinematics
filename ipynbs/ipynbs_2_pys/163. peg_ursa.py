#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import numpy as np
import pandas as pd
from astropy.table import QTable
import matplotlib.pyplot as plt


img = fits.open('DT22+DHT_masked.fits')[0].data.transpose()


img = img.swapaxes(1,2)


hdr = fits.open('DT22+DHT_masked.fits')[0].header


hdr


np.nanmin(img)


np.nanmax(img) + 12.09968


vls = np.arange(146)*0.65 - 47.1


lts = (np.arange(677)-321)*0.25


lns = np.arange(1441)*0.25-60


plt.figure(figsize=(10,10))
plt.contourf(np.nansum(img,axis=0),norm='log')#,extent=[-60,300,-80,89])
plt.show()


img.shape


plt.figure(figsize=(14,7))
plt.contourf(np.log10(np.nansum(img[:,235:277,885:960],axis=0)*0.65*1.97e20),levels=np.arange(20.7,22.5,0.1))
plt.colorbar()
plt.xlim(72,0)


tau = img[70:95,235:277,885:960]


plt.plot(np.nanmean(tau,axis=(1,2)))
plt.xlim(70,95)


ora = img[60:110,235:252,1060:1110]


plt.figure(figsize=(14,7))
plt.contourf(np.log10(np.nansum(img[:,235:252,1060:1110],axis=0)*0.65*1.97e20),levels=np.arange(20.7,22.5,0.1))
plt.colorbar()
plt.xlim(52,0)


plt.plot(np.nanmean(ora[60:110],axis=(1,2)))


import astropy.units as u
import astropy.constants as c


nhora = np.nansum(ora,axis=0)*0.65*1.97e20


np.nansum(((nhora*(nhora>3e20))*1.36*2*\
c.m_p*u.cm**-2).to(u.solMass/u.pc**2))*\
(0.25*u.deg.to(u.rad)*432*u.pc)**2


np.nansum(((np.nansum(ora,axis=0)*0.65*1.97e20*1.36*2*\
c.m_p*u.cm**-2).to(u.solMass/u.pc**2)))*\
(0.25*u.deg.to(u.rad)*432*u.pc)**2


94990/85883


nhtau = np.nansum(tau,axis=0)*0.65*1.97e20


np.nansum((((nhtau*(nhtau>3e20))*1.36*2*\
c.m_p*u.cm**-2).to(u.solMass/u.pc**2)))*\
(0.25*u.deg.to(u.rad)*141*u.pc)**2


22263.12/18272


22937/18272


#so far with orion a & taurus it looks like using the dht_22 data gives masses 10-25% larger
#than the other dht survey
#but the thing is, the spatial resolution goes down by 1/2 and the velocity resolution goes up by 1/2
#so harder to match with hi4pi 

#it's also kind of weird that this is giving diff mass values, bc all that this dataset did was
#combine the old dht data with the northern sky .. unless orion et al were observed again?


# plt.figure(figsize=(14,7))
plt.contourf(np.log10(np.nansum(img[:,430:510,780:900],axis=0)*0.65*1.97e20),levels=np.arange(19.5,21.5,0.1))
plt.colorbar()
# plt.xlim(52,0)


ursa = img[55:90,430:510,780:900]


plt.plot(np.nanmean(ursa[55:90],axis=(1,2)))


np.nansum(((np.nansum(ursa,axis=0)*0.65*1.97e20*1.36*2*\
c.m_p*u.cm**-2).to(u.solMass/u.pc**2)))*\
(0.25*u.deg.to(u.rad)*371*u.pc)**2


np.nansum((np.nansum(ursa,axis=0)*0.65*1.97e20)>3e20)


np.sqrt((114*(0.25*u.deg.to(u.rad)*371*u.pc)**2)/np.pi)


# plt.figure(figsize=(14,7))
plt.contourf(np.log10(np.nansum(img[:,120:230,570:650],axis=0)*0.65*1.97e20),levels=np.arange(19,21.5,0.1))
plt.colorbar()
# plt.xlim(150,0)


peg = img[45:85,120:230,570:650]


plt.plot(np.nanmean(peg,axis=(1,2))[45:85])


np.nansum(((np.nansum(peg,axis=0)*0.65*1.97e20*1.36*2*\
c.m_p*u.cm**-2).to(u.solMass/u.pc**2)))*\
(0.25*u.deg.to(u.rad)*247*u.pc)**2


np.nansum((np.nansum(peg,axis=0)*0.65*1.97e20)>3e20)


np.sqrt((311*(0.25*u.deg.to(u.rad)*247*u.pc)**2)/np.pi)
