#===================================================================================
import datacube as dc
import clouds as cld
import interpolate as ip
import bootleg as btlg

from numba import njit

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c

from astrodendro import ppv_catalog, Dendrogram,periodic_neighbours 

from scimes import SpectralCloudstering
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
#===================================================================================


#===================================================================================
mtdt = {}
mtdt['data_unit'] = u.K
mtdt['beam_major'], mtdt['beam_minor'] = 8.4*u.arcmin, 8.4*u.arcmin
mtdt['spatial_scale'] = DHT_CO.db #same as DHT_CO.dl
mtdt['velocity scale'] = DHT_CO.dv
mtdt['wavelength'] = c.c/(115e9*u.Hz)
#===================================================================================


#===================================================================================
def compute_dend(min_val, min_del, min_px, cube=DHT_CO.cube,neighbors=periodic_neighbours(2)):
    d = Dendrogram.compute(cube, min_val, min_del, min_px)
    return d
#------------------------------------------------------------------------
def make_cat(dend,metadata=mtdt,head=DHT_CO.data_header,**kwargs):
    catalog = ppv_catalog(dend,metadata,**kwargs)
    return catalog
#------------------------------------------------------------------------
def compute_cluster(dend,metadata=mtdt,head=DHT_CO.data_header,**kwargs):
    catalog = ppv_catalog(dend,metadata)
    clust = SpectralCloudstering(dend,catalog,head,**kwargs)

    return clust
#===================================================================================


#===================================================================================
sig = DHT_CO.calculate_rms_noise()

d = compute_dend(sig,20*sig,500)
# kw_bool = {'save_all_leaves':True}
# clus = compute_cluster(d,**kw_bool)

lvs = [d[i] for i in range(len(d)) if d[i].is_leaf]
lf_ids = [d[i].idx for i in range(len(d)) if d[i].is_leaf]

lf_catalog = make_cat(lvs,verbose=False)
#===================================================================================


#===================================================================================
lvs_mask = np.zeros_like(DHT_CO.cube)
for leaf in lvs:
    m = leaf.get_mask()
    m_labeled = m*leaf.idx
    lvs_mask += m_labeled
#===================================================================================


#===================================================================================
ids = lf_ids

ls_px = round(pd.Series(lf_catalog['x_cen'].data)).astype('int')
bs_px = round(pd.Series(lf_catalog['y_cen'].data)).astype('int')
vs_px = round(pd.Series(lf_catalog['v_cen'].data)).astype('int')

ls = cld.ind_ang('get angular', DHT_CO.lons, ls_px)
bs = cld.ind_ang('get angular', DHT_CO.lats, bs_px)
vs = cld.ind_ang('get angular', DHT_CO.vels, vs_px)

sources_w_catalog_px = dc.generate_sources(['!id','lon','lat','vel'],[ids,ls,bs,vs])
#===================================================================================


#===================================================================================
lf_catalog.write('lvs_cat.csv',format='pandas.csv',overwrite=True)
btlg.dend_to_fits(d,'COGAL_dendrogram.fits')

hdu = fits.PrimaryHDU(lvs_mask)
hdu.writeto('lvs_ind_map.fits',overwrite=True)
#===================================================================================






