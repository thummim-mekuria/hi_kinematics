#!/usr/bin/env python
# coding: utf-8


import datacube as dc
import clouds as cld

import numpy as np
# %matplotlib

import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram
from astropy.visualization import ZScaleInterval as zscale

import proposal as prp


obsnu2 = 900e9 *u.Hz


nus = np.full(13,obsnu2)*u.Hz


nus


Ts = np.full(13,70)*u.K


pos_gals = pd.read_csv('notion.csv')

sfrs = 10**(pos_gals['SFR_$[M_\\odot$_$yr^{-1}]$'])
mstrs = 10**(pos_gals['M_str_$[M_\\odot]$'])

mus = pos_gals['$\\mu$']
for m, mu in enumerate(mus):
    if not np.isnan(mu):
        mstrs[m] = mstrs[m]/mu
        
mstrs = mstrs.to_numpy() *c.M_sun
sfrs = sfrs.to_numpy()*c.M_sun/u.yr
zs = pos_gals['z']


snus = np.zeros(13)*u.g/u.s**2
for i in range(13):
    snu = prp.S_nu(nus[i], Ts[i], mstrs[i], zs[i], sfrs[i])
    snus[i] = snu
#     print(snu)


snus.to(u.Jansky)


mm11 = c.c/(1.1*1e-3*u.m)


c.c/mm11


(10**(0.09))/10.1


50/4


1.1/10.1


1e-15*(10**(16))/1.7


u.Jy.cgs


prp.T_from_Snu(mm11,0.099*1e-3*u.Jansky,1.97e9*c.M_sun,8.38,20*c.M_sun/u.year)


obsnu2


c.h*obsnu2/(c.k_B*50*u.K)


plt.plot(snus.to(u.Jansky).value,marker='+')


np.abs(53.15500-53.15528)


0.00028 is 0.06


0.06 seconds


obsnu2


mstrs


snus


secondgo = pd.read_csv('notion2.csv')


zs2 = secondgo['z']


mstrs2 = c.M_sun*(10**(secondgo['M_str'])).to_numpy()


sfrs2 = c.M_sun/u.yr*(10**(secondgo['SFR'])).to_numpy()


obsnu2 = 300*1e9*u.Hz


nus = np.full(13,obsnu2)*u.Hz


snus2 = np.zeros(13)*u.g/u.s**2
for i in range(13):
    snu = prp.S_nu(nus[i], Ts[i], mstrs2[i], zs2[i], sfrs2[i])
    snus2[i] = snu
#     print(snu)


snus2.to(u.Jansky)


201734*u.day.to(u.year)
