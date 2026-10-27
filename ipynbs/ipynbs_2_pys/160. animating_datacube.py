#!/usr/bin/env python
# coding: utf-8


import os, sys, string, matplotlib, aplpy

import matplotlib.pyplot as plt
import numpy as np
import pylab as py
from astropy.io import fits
from matplotlib import colors, cm

import datacube as dc


#Colormap I like to use
cdict_bkcmap={
    'red'  :  ((0., 0., 0.), (0.25,0.,0.), (0.5,1.,1.), (0.75,1.0,1.0),  (1., 1., 1.)),
    'green':  ((0., 0., 0.), (0.25,0.,0.), (0.5,0.,0.), (0.75,1.0,1.0),  (1., 1., 1.)),
    'blue' :  ((0., 0., 0.), (0.25,1.,1.), (0.5,0.,0.), (0.75,0.0,0.0),  (1., 1., 1.))
    }

bkcmap = matplotlib.colors.LinearSegmentedColormap('bkcmap', cdict_bkcmap,1024)


#LAB fits -- has been superseded but this is the one that works so will use it for now
LAB_HI = dc.unpack_fits_cube('labh.fits')
HI_cube = LAB_HI[0]
HI_vels = LAB_HI[1]
HI_lats = LAB_HI[2]
HI_lons = LAB_HI[3]
HI_units = LAB_HI[4] 
HI_units['vels'] = u.m/u.s
#had to go in and change units of vels to m/s, so need to edit dc.unpack_fits_cube
#to make sure it reads the header for units instea of hard-coding 

HI_dv = ((HI_vels[1] - HI_vels[0])*HI_units['vels']).to(u.km/u.s)
HI_db = (HI_lats[1] - HI_lats[0])*HI_units['lats']
HI_dl = (HI_lons[1] - HI_lons[0])*HI_units['lons']

HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0)
