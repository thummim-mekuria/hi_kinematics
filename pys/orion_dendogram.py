# %%time

import datacube as dc
import clouds as cld

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c

from astropy.visualization import ZScaleInterval as zscale

from scimes import SpectralCloudstering
from astrodendro import Dendrogram

from astrodendro import ppv_catalog





zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


def rms_noise(thing):
    mean = np.nanmean(thing)
    noise_square = (thing - mean)**2
    
    return np.sqrt(np.nanmean(noise_square))



orion_clob = dc.unpack_fits_cube('DHT27.fits', 'DHT')

# cube, vels, lats, lons, units, header = orion.values()

#values from tutorial(https://scimes.readthedocs.io/en/latest/tutorial.html)
#may need to change later
dend_sigma = 0.3 #K, noise level
dend_ppb = 1.3 #pixels/beam

def compute():
    orion_dend = Dendrogram.compute(orion_clob.cube, min_value = dend_sigma,\
        min_delta = 2*dend_sigma, min_npix = 3*dend_ppb)

    return orion_dend

def launch_viewer(d):
    v = d.viewer()
    v.show()


dend = compute()
# launch_viewer(dend)


# mtdt = {}
# mtdt['data_unit'] = u.Jy #works for u.Jy (but not u.K ?? )

# catalog = ppv_catalog(dend, mtdt)

# cluster = SpectralCloudstering(dend, catalog, header)

# cluster.showdendro()














