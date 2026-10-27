import datacube as dc
import clouds as cld
import progress as pgs
import interpolate as ip

from numba import njit

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram

#=============================================================================================

DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')
DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 

smoothed_cube = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smoothed_cube)
wo_gal_cent[:,:,:641]        = smoothed_cube[:,:,:641]
wo_gal_cent[:,:,2241:]       = smoothed_cube[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smoothed_cube[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smoothed_cube[:,281:,641:2241]

#=============================================================================================

sig = np.nanstd(wo_gal_cent)
d = Dendrogram.compute(wo_gal_cent,sig,sig,50)

# d_wocent_wo_npx.save_to('d_wocent_wo_npx.fits')

#=============================================================================================

colnames=['ID','l(deg)','b(deg)','v_lsr(km/s)','pfar']
sources = pd.DataFrame(columns=colnames)

N = len(d.trunk)
IDs = np.zeros(N)
ls = np.zeros(N)
bs = np.zeros(N)
vs = np.zeros(N)
far_ps = np.full(N, 0.5)
# clobs = np.zeros(N)

sources['pfar'] = far_ps

#=============================================================================================
prev = 0
for t,tr in enumerate(d.trunk):
    m = tr.get_mask()
    cb = m*wo_gal_cent
    clob = DHT_CO.clone()
    clob.cube = cb
    v,b,l = dc.weighted_vbl(clob)
    IDs[t] = str(t)
    ls[t] = l
    bs[t] = b
    vs[t] = v
    # clobs[t] = clob


    prog = pgs.progress(t,len(d.trunk),prev)
    # print(t)
    if prog:
        print(t,prog,'%')
        prev +=10


for l in range(len(ls)):
    if ls[l] < 0:
        ls[l] +=360



sources['ID'] = IDs
sources['l(deg)'] = ls
sources['b(deg)'] = bs
sources['vlsr(km/s)'] = vs


sources.to_csv('sources_info.inp',sep=' ',index = False) #then go in & add ! in front of the line with column names


print('done')











