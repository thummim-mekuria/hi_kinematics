#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))


# moving forward:
# start by reading in more detail the Carina Arm paper to look at their methodoloy
# look into clumpfind & decomposing the HI4PI dataset before settling on a specific approach for identifying the controld HI group
# 
# & in parallel start examining the GMC-ful HI clouds and see if i can reproduce plots / results that are already knonwn
# take a closer look at the methodology of papers that use the Dame CO data & use that as a way to determine a better way to fdefine the contours of the GMCs
# compare accumulation radius bounds with nia's / other published values as a sanity check


# %%time

import datacube as dc
import clouds as cld

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram
from astropy.visualization import ZScaleInterval as zscale


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf

DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 


#----------------------
cat = pd.read_csv('Mekuria_catalog.csv',sep=',',skiprows=5)
names = cat['name']


l_unit = u.deg
b_unit = u.deg
d_unit = u.pc


#----------------------


def read_DHT(ns='all'):
    output={}
    if ns == 'all':
        ns = ['03','06','09','12','16','19','21','26','27','37']
        #for now, figure out DHT01 (bad data) 
        #& DHT18 (two ceph ds, beeg) 
        #& DHT08 (beeg) 
        #& DHT24, DHT34 (lam or & chamaeleon are disjointed down to the really low sig levels)
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output


#---------------------
surveys = read_DHT()
for s in surveys:
    k = np.where(cat['survey']==s)[0]

    these = cat['name'][k]
    
    nums = np.zeros(len(these))
    for i, this in enumerate(these):
        nums[i] = np.where(cat['name']==this)[0]
        
    nums = nums.astype(np.int64)
    cube = surveys[s]
    surveys[s] = [cube, nums]


cat


lis = [0,1,2,4,5,6,7,8,9,10,11,12,14,15,17]


#^this are the clouds that i have selected so far, for tomorrow: plot the data from the park et al tables, 
#also decompose the dataset and see if i get a similar population of MCs and HIs


for k in lis:
    clob = cat['clob'][k]
    lim = [clob.lons[0],clob.lons[-1], clob.lats[0], clob.lats[-1]]

    plt.imshow(np.nansum(clob.cube, axis = 0),**kw_imshow, extent = lim)
    plt.title(cat['name'][k])
    plt.show()


cat


DHT06 = surveys['DHT06'][0]


s = 0.25


plt.imshow(np.nansum(Dendrogram.compute(DHT34.cube,s,0,400).trunk[0].get_mask()*DHT34.cube, axis = 0))


#baseline: 3sig, sig, 400 --> do position cut here if len(d.trunk) > 1
#06, 37, 34: sig, sig, 500
#21 6sig, sig, 1000


for k in d21.trunk:
    print(k.smallest_index)


print(DHT21.lats[152], DHT21.lons[277],'\n',
DHT21.lats[123], DHT21.lons[169],'\n',
DHT21.lats[40], DHT21.lons[274])


cat.columns
