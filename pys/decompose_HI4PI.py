import datacube as dc
import clouds as cld
import progress as pgs
import interpolate as ip
import bootleg as btlg

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram

#=============================================================================================

HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')

rms_noise = 43*1e-3
HI4PI.data_header['T_rms'] = rms_noise #u.K, median value from Table 1 of Dame paper: 

d = Dendrogram.compute(HI4PI.cube,rms_noise,rms_noise,50)

btlg.dend_to_fits(d, 'HI4PI_decomp.fits')

