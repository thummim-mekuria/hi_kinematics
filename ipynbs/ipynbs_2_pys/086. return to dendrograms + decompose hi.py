#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import datacube as dc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
import fun as fns
import astropy.units as u
import astropy.constants as c
from matplotlib.patches import Circle,Rectangle


hi4pi,dht = fits.open('both_same_grid_720d.fits')
zcat = QTable(fits.open('final_cat.fits')[1].data)


lons = np.arange(5760)*-0.125 + 360
lats = np.arange(800)*0.125 - 50
vels = np.arange(155)*1.3 - 100


hi_conv = (1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = (1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


whi = (hi_conv*u.K*1.3*u.km/u.s*np.nansum(hi4pi.data,axis=0)).to(u.solMass/u.pc**2).value
wco = (co_conv*u.K*1.3*u.km/u.s*np.nansum(dht.data,axis=0)).to(u.solMass/u.pc**2).value


whole_gal_ds = pd.read_csv('whole_gal_dists.csv')
zcat = QTable(fits.open('final_cat.fits')[1].data)
gless_cat = QTable(fits.open('gmcless_cat.fits')[1].data)


h2s = fits.open('h2_clouds.fits')
his = fits.open('hi_clouds.fits')
gls = fits.open('gmcless_tiles.fits')


from astrodendro import Dendrogram


np.nanstd(hi4pi.data)


hiqm = Dendrogram.compute(hi4pi.data[:,:,1360:4400],min_value=10,min_delta=50,min_npix=1000)#l=190->-190


hiqm30 = Dendrogram.compute(hi4pi.data[:,:,1360:4400],min_value=10,min_delta=30,min_npix=1000) 


hiqm10 = Dendrogram.compute(hi4pi.data[:,:,1360:4400],min_value=10,min_delta=10,min_npix=1000) #best one?


mom2_gal_shrt = fns.get_momN(hi4pi.data[:,:,1360:4400],[0,155],2)


[(plt.imshow(np.nansum(t.get_mask()*hi4pi.data[:,:,1360:4400],axis=0),origin='lower',
             cmap=dc.cmap_min0('viridis')),
) for t in hiqm10.leaves]


[(plt.imshow(np.nansum(t.get_mask()*hi4pi.data[:,:,1360:4400],axis=0),origin='lower',
             cmap=dc.cmap_min0('viridis')),
) for t in hiqm30.leaves]


[print(t.get_npix(),t.get_peak())
 for t in hiqm30.leaves]


[(np.nanmean(np.nansum(his[i].data,axis=0)*hi_conv*1.3*u.K*u.km/u.s).to(u.solMass/u.pc**2)) for i in range(17)]


dens_enf_5 = [t.idx for t in hiqm10.leaves if (np.nanmax(np.nansum(t.get_mask()*hi4pi.data[:,:,1360:4400],
                            axis=0))*hi_conv*1.3*u.K*u.km/u.s).to(u.solMass/u.pc**2) > 5*u.solMass/u.pc**2 ]


#choose things that dont fall off edges


goodvs = [hiqm10[i].idx for i in dens_enf_5 
               if np.nanmean(hiqm10[i].indices()[0])*1.3-100 - 2*np.nanstd(hiqm10[i].indices()[0])*1.3 > -100.
               and np.nanmean(hiqm10[i].indices()[0])*1.3-100 + 2*np.nanstd(hiqm10[i].indices()[0])*1.3 < 100.]


good_vs_bs = [hiqm10[i].idx for i in rvsd_goodvs 
               if np.nanmean(hiqm10[i].indices()[1])*0.125-50 - 2*np.nanstd(hiqm10[i].indices()[0])*0.125 > -50.
               and np.nanmean(hiqm10[i].indices()[1])*0.125-50 + 2*np.nanstd(hiqm10[i].indices()[0])*0.125 < 50.]


good_lbvs = [hiqm10[i].idx for i in good_vs_bs 
               if np.nanmean(hiqm10[i].indices()[2]+1360)*-0.125+360 + 2*np.nanstd(hiqm10[i].indices()[0])*0.125 < 360.
               and np.nanmean(hiqm10[i].indices()[2]+1360)*-0.125+360 - 2*np.nanstd(hiqm10[i].indices()[0])*0.125 > -360.]


good_lf_coords = [(np.nanmean(np.array(hiqm10[i].indices()),axis=1)+\
                       [0,0,1360])*[1.3,0.125,-0.125] + [-100,-50,360]
 
 for i in good_lbvs]


dfgood = pd.DataFrame(good_lf_coords)


dfgood['idx'] = good_lbvs


dfgood.columns=['vel<km/s','lat<deg>','lon<deg>','idx']


dfgood = dfgood[['idx','lon<deg>','lat<deg>','vel<km/s']]


import bootleg as btlg


btlg.dend_to_fits(hiqm10,'decomposed_HI4PI.fits')


dfgood.to_csv('hi4pi_dendleaves_5solMass.csv',index=False)


dfgood
