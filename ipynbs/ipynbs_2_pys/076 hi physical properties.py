#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Ellipse

import datacube as dc

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable


zcat = QTable(fits.open('final_cat.fits')[1].data)


hi4pi,dht = fits.open('both_same_grid_720d.fits')


v0s,v1s = np.transpose(np.array((zcat['vel_HI<km/s>'] + 100.)/1.3).astype('int'))


l0s = np.array(np.nanmean(zcat['l_pix'],axis=1) - np.array(zcat['radius_acc_HI<deg>']/0.125)).astype('int')
l1s = np.array(np.nanmean(zcat['l_pix'],axis=1) + np.array(zcat['radius_acc_HI<deg>']/0.125)).astype('int')

b0s = np.array(np.nanmean(zcat['b_pix'],axis=1) - np.array(zcat['radius_acc_HI<deg>']/0.125)).astype('int')
b1s = np.array(np.nanmean(zcat['b_pix'],axis=1) + np.array(zcat['radius_acc_HI<deg>']/0.125)).astype('int')


hicubes = [hi4pi.data[v0s[i]:v1s[i],b0s[i]:b1s[i],l0s[i]:l1s[i]] for i in range(17)]


hi_conv = (1.3*u.km/u.s *1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p *u.K).to(u.solMass/u.pc**2) #last term is data


px_size = 0.125*u.deg.to(u.rad)


whi   = [np.nansum(hicubes[i],axis=0) for i in range(17)]
sighi = [(whi[i] * hi_conv).to(u.solMass/u.pc**2) for i in range(17)]
mhi   = [(np.nansum(sighi[i]* (px_size*zcat['distance<pc>'][i]*u.pc)**2)).to(u.kg) for i in range(17)]


zcat['mass_HI<kg>'] = mhi


np.array(zcat['mass_GMC<kg>'])/zcat['mass_HI<kg>'].value


zcat.rename_columns(['l_pix','b_pix','v_pix'],['l_pix_GMC','b_pix_GMC','v_pix_GMC'])


zcat['l_pix_HI'] = np.array((l0s,l1s)).transpose()
zcat['b_pix_HI'] = np.array((b0s,b1s)).transpose()
zcat['v_pix_HI'] = np.array((v0s,v1s)).transpose()


specs = [np.nansum(hicubes[i],axis=(1,2)) for i in range(17)]


hicubes_long = [hi4pi.data[:,b0s[i]:b1s[i],l0s[i]:l1s[i]] for i in range(17)]


vos


hicubes_long_o = [hi4pi.data[vos[i]:vis[i],b0s[i]:b1s[i],l0s[i]:l1s[i]] for i in range(17)]


specs_long = [np.nansum(hicubes_long_o[i],axis=(1,2)) for i in range(17)]


zcat['vdispersion_HI<km/s>'] = np.array([(get_dispers(specs[i]).stddev.value)
  for i in range(17)])*1.3


# [(plt.plot(specs_long[i]),\
#  plt.plot(get_dispers(specs_long[i])(np.arange(len(specs_long[i])))),\
#   plt.title(zcat['names'][i]),
#   plt.show(),
[(print(zcat['names'][i]),
    print(100*(get_dispers(specs_long[i]).stddev.value - \
           get_dispers(specs[i]).stddev.value)/\
          get_dispers(specs[i]).stddev.value),
    print(np.nanmean(zcat['lat_GMC<deg>'][i])),
  print('\n\n')

  
 ) for i in range(17)]


zcat['v_pix_HI'] = np.array([[50,  100],
[ 65,  95],
[ 82,  102],
[ 60,  100],
[ 40, 120],
[ 60, 117],
[ 50,  100],
[ 92, 107],
[ 58,  102],
[ 50,  120],
[ 50,  126],
[ 50,  120],
[ 40, 120],
[ 55,  104],
[ 50,  100],
[ 60, 104],
[ 55,  95]])


zcat['vel_HI<km/s>'] = zcat['v_pix_HI']*1.3-100.


vos,vis = np.array([[50,  100],
[ 65,  95],
[ 82,  102],
[ 60,  100],
[ 40, 120],
[ 60, 117],
[ 50,  100],
[ 92, 107],
[ 58,  102],
[ 50,  120],
[ 50,  126],
[ 50,  120],
[ 40, 120],
[ 55,  104],
[ 50,  100],
[ 60, 104],
[ 55,  95]]).transpose()


[(plt.plot(specs_long[i]),\
 plt.plot(get_dispers(specs_long[i])(np.arange(len(specs_long[i])))),\
  plt.title(zcat['names'][i]),
  plt.show(),
  print(get_dispers(specs_long[i])),
  
 ) for i in range(17)]


[(plt.plot(specs[i]),\
 plt.plot(get_dispers(specs[i])(np.arange(len(specs[i])))),\
  plt.title(zcat['names'][i]),
  plt.show(),
  print(get_dispers(specs[i])),
  
 ) for i in range(17)]


from astropy.modeling import models, fitting
def get_dispers (spec):
    fitter = fitting.LevMarLSQFitter()
#     print(np.nanmax(spec),np.nanargmax(spec),5)
    
#     ginit = models.Gaussian1D(1.,0.,1.)
    ginit = models.Gaussian1D(np.nanmax(spec),np.nanargmax(spec),5.)

    gfit = fitter(ginit,np.arange(len(spec)),spec)
#     print(gfit.parameters)
    return gfit


zcat.colnames


zcat = zcat[['names',
 'lon_GMC<deg>',
 'lat_GMC<deg>',
 'vel_GMC<km/s>',
 'vel_HI<km/s>',
 'distance<pc>',
 'mass_GMC<kg>',
 'mass_HI<kg>',
 'radius_eff_GMC<pc>',
 'radius_eff_GMC<deg>',
 'radius_acc_HI<pc>',
 'radius_acc_HI<deg>',
 'vdispersion_GMC<km/s>',
 'vdispersion_HI<km/s>',
 'l_pix_GMC',
 'b_pix_GMC',
 'v_pix_GMC',
 'l_pix_HI',
 'b_pix_HI',
 'v_pix_HI'
 ]]


zcat.write('final_cat.fits',overwrite=True)
