#!/usr/bin/env python
# coding: utf-8


import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.io import fits

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import fun as fns
import datacube as dc


((0.57*u.cm**-3)*c.m_p*1.4).to(u.solMass/u.pc**3)


cat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


new_rhs = np.cbrt(cat['mass_GMC<Msun>']/(0.02*2*np.pi))


#came back and realized that i made a mistake, ie i used 1.43 as the mu before, and when correcting it (1.36) whats below is the difference in size
#in pixels that the accumulation radius ends up being ... but im not gonna redo anything until later (ughh is this shoowing myself in the foot?
#maybe i can get away with not using the california ae so much

#or find a way to argue that 1.4


270/296


np.array((np.around((new_rhs/cat['distance<pc>'])*u.rad.to(u.deg)/0.125),
          np.around((QTable(fits.open('sighhplease.fits')[1].data)['radius_HI<pc>']/cat['distance<pc>'])*u.rad.to(u.deg)/0.125))).transpose()


lc = np.nanmean(cat['lon_GMC<deg>'],axis=1)
bc = np.nanmean(cat['lat_GMC<deg>'],axis=1)


newr,oldr=((np.array((new_rhs,cat['radius_HI<pc>'])).transpose()/\
cat['distance<pc>'][:,np.newaxis])*u.rad.to(u.deg)).transpose()


another_one = cat[['name',
 'lon_GMC<deg>',
 'lat_GMC<deg>',
 'distance<pc>',
 'loc',
 'lic',
 'boc',
 'bic']]


another_one = QTable(fits.open('sighhplease.fits')[1].data) #the lons & lats are diff here


another_one[1]['lon_GMC<deg>'] = [174,154]
another_one[1]['lat_GMC<deg>'] = [-14,-1]


another_one[21]['lat_GMC<deg>'] = [-21,-10]


another_one[20]['lon_GMC<deg>'] = [-151,-155]


another_one['lon_GMC<deg>'][1] = [20.,30.]


another_one['lon_GMC<deg>'][6] = [-138,-142]
another_one['lat_GMC<deg>'][6] = [-12.5,-6.5]


another_one['lon_GMC<deg>'][18] = [162,154]


another_one['lon_GMC<deg>'][7] = [-164,-174]
another_one['lat_GMC<deg>'][7] = [-4,6]


hi4pi,dht = fits.open('both_same_grid_720d.fits')


c0 = np.nansum(dht.data,axis=0)*1.3*1.97e20


h0 = np.nansum(hi4pi.data,axis=0)*1.3*1.82e18


dust = fits.open('mean_and_std_lbd.fits')[1].data


#for the distances just say that not all of the cloud has los through it from zucker, and also there are 
#significant (~R_GMC) variations across the face, and we are flattening that 3d detail when using a single value


#esp for the clouds with background stuff (that doesnt get removed by the spectra)
#say this is the n_H thershold we're using to define the cloud (eg gemini--21,madd--20.5,mob1--21)


#also cite papers for each cloud (ie discovery, description)


lo,li = ((another_one['lon_GMC<deg>']-360)/-0.125).astype('int').transpose()
bo,bi = ((another_one['lat_GMC<deg>']+50)/0.125).astype('int').transpose()


another_one['loc'] = lo
another_one['lic'] = li
another_one['boc'] = bo
another_one['bic'] = bi


bc = np.nanmean((bo,bi),axis=0).astype('int')


lc = np.nanmean((lo,li),axis=0).astype('int')


newr_px = (newr/0.125).astype('int')


[(
 plt.contourf(np.log10(h0),levels=np.arange(20,22.0,0.1)),#,extent=[360,-360,-50,50]),
 plt.contour(np.log10(c0),levels=[20.5,21],colors='black',linewidths=0.6),
 plt.xlim(lo[i]-newr_px[i],li[i]+newr_px[i]),
 plt.ylim(bo[i]-newr_px[i],bi[i]+newr_px[i]),
 plt.gca().set_aspect('equal'),                           
 plt.title(cat['name'][i]),
 plt.show())
 for i in range(22)]


vels = np.arange(155)*1.3-100


coco


coco = [dht.data[:,bo[i]:bi[i],lo[i]:li[i]] for i in range(22)]


hico = [hi4pi.data[:,bo[i]:bi[i],lo[i]:li[i]] for i in range(22)]


hihi = [hi4pi.data[:,
                   bc[i]-newr_px[i]:bc[i]+newr_px[i],
                   lc[i]-newr_px[i]:lc[i]+newr_px[i]] 
        for i in range(22)]

cohi = [dht.data[:,
                   bc[i]-newr_px[i]:bc[i]+newr_px[i],
                   lc[i]-newr_px[i]:lc[i]+newr_px[i]] 
        for i in range(22)]


hi2hi = [hi4pi.data[:,
                   bc[i]-2*newr_px[i]:bc[i]+2*newr_px[i],
                   lc[i]-2*newr_px[i]:lc[i]+2*newr_px[i]] 
        for i in range(22)]

co2hi = [dht.data[:,
                   bc[i]-2*newr_px[i]:bc[i]+2*newr_px[i],
                   lc[i]-2*newr_px[i]:lc[i]+2*newr_px[i]] 
        for i in range(22)]


def make_fits_from_list(array,filename):
    hdulist = fits.HDUList([fits.PrimaryHDU(array[0])]+[fits.ImageHDU(array[i]) for i in range(1,22)])
    hdulist.writeto(filename,overwrite=True)


make_fits_from_list(coco,'h2_h2sized.fits')
make_fits_from_list(hico,'hi_h2sized.fits')
make_fits_from_list(hihi,'hi_hisized.fits')
make_fits_from_list(cohi,'co_hisized.fits')
make_fits_from_list(hi2hi,'hi_2xhisized.fits')
make_fits_from_list(co2hi,'co_2xhisized.fits')


another_one['radius_HI<pc>'] = new_rhs


another_one['boh'] = bc-newr_px
another_one['bih'] = bc+newr_px
another_one['loh'] = lc-newr_px
another_one['lih'] = lc+newr_px


another_one.write('sighhplease.fits',overwrite=True)


[(
 plt.plot(100*np.nanmean(dht.data[:,
                                  bo[i]:bi[i],
                                  lo[i]:li[i]],
                         axis=(1,2))),
 plt.plot(np.nanmean(hi4pi.data[:,
                                  bo[i]:bi[i],
                                  lo[i]:li[i]],
                         axis=(1,2))),
 plt.plot(np.nanmean(hi4pi.data[:,
                                bc[i]-newr_px[i]:bc[i]+newr_px[i],
                                lc[i]-newr_px[i]:lc[i]+newr_px[i]],
                     axis=(1,2))),
 plt.title(cat['name'][i]),
 plt.show())
 for i in range(22)]


i=1


cal=fns.get_mom1(hi4pi.data[:,bo[i]-50:bi[i]+50,lo[i]-50:li[i]+50], vels)


plt.contourf(aq1*(aq1<12))


plt.hist(np.ravel(aq1),bins=300)
plt.vlines(12,0,700,color='red')
plt.show()


plt.hist(np.ravel(cal),bins=100)
plt.vlines([-27,-19,-8,-4,1],0,600,color='red')
plt.show()


plt.contourf(cal*(cal<-4)*(cal>-8))
