#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


all_ds = pd.read_csv('../new_reid/summary.prtmary.prt',skiprows=4,header=None,delim_whitespace=True)


all_ds.columns = ['lon<deg>','lat<deg>','vel<km/s>','v_unc<km/s>','d1<kpc>','d1_unc<kpc>','prob1','arm1',\
                 'd2<kpc>','d2_unc<kpc>','prob2','arm2']


good_snr = all_ds[all_ds['d1<kpc>']/all_ds['d1_unc<kpc>'] > 10]


good_snr5 = all_ds[all_ds['d1<kpc>']/all_ds['d1_unc<kpc>'] > 5]


local = all_ds[all_ds['d1<kpc>']<2]


plt.figure(figsize=(20,9))
plt.scatter(local['lon<deg>'][::37],local['vel<km/s>'][::37],
            s=15,c=local['d1<kpc>'][::37]*1000,
            cmap='Paired')
plt.colorbar()
plt.xlim(180,-180)
# plt.ylim(-25,25)


plt.figure(figsize=(20,9))
plt.scatter(all_ds['lon<deg>'][::37],all_ds['vel<km/s>'][::37],
            s=15,c=all_ds['d1<kpc>'][::37],
            cmap='Paired',vmax=10)
plt.colorbar()
plt.xlim(180,-180)


plt.figure(figsize=(20,9))
plt.scatter(local['lon<deg>'][::37],local['lat<deg>'][::37],
            s=15,c=local['d1<kpc>'][::37]/local['d1_unc<kpc>'][::37],
            cmap='Paired',vmin=1,vmax=10)
plt.colorbar()
plt.xlim(180,-180)


plt.figure(figsize=(20,9))
plt.scatter(all_ds['lon<deg>'][::37],all_ds['vel<km/s>'][::37],
            s=15,c=all_ds['d1<kpc>'][::37]/all_ds['d1_unc<kpc>'][::37],
            cmap='Paired',vmax=10)
plt.colorbar()
plt.xlim(180,-180)


plt.figure(figsize=(10,5))
plt.scatter(good_snr5['lon<deg>'],good_snr5['vel<km/s>'],s=10)
plt.xlim(180,-180)


plt.figure(figsize=(10,5))
plt.scatter(good_snr['lon<deg>'],good_snr['vel<km/s>'],s=10)
plt.xlim(180,-180)


1901901/1e6
