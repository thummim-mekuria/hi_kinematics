#===================================================================================
import datacube as dc
import clouds as cld
import interpolate as ip
import bootleg as btlg

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import math

# %matplotlib ipympl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram,periodic_neighbours
#===================================================================================

#===================================================================================
DHT_CO = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')

DHT_CO.cube = DHT_CO.cube[209:285,:,:] #trim to narrower velocity range
DHT_CO.vels = DHT_CO.vels[209:285]

smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole) #interpolate bc 'mom' is staccatoed
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641] #take out |b|<5 for |l|<100
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]

DHT_CO.cube = wo_gal_cent


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]


def display(thing_3d,order=0,cmp = dc.cmap_grey0('viridis'),**kwargs):
    if type(thing_3d) == u.Quantity:
        thing_3d = thing_3d.value
    plt.imshow(np.nansum(thing_3d,axis=order),origin='lower',cmap=cmp,**kwargs)
    plt.show()