import datacube as dc
import clouds as cld
import transform as tr
import interpolate as ip

from numba import njit

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astropy.visualization import ZScaleInterval as zscale

from astrodendro import ppv_catalog
from astrodendro import Dendrogram


from scimes import SpectralCloudstering


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}

print("imports done, 0th functions defined")

#--------------------------------------------------------------------------------------------------------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')

#----------------------------

DHT_CO = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')
DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 
smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole)
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641]
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]

DHT_CO.cube = wo_gal_cent


print("clobs set up")
#--------------------------------------------------------------------------------------------------------------------

co_sig = np.nanstd(DHT_CO.cube)
dco = Dendrogram.compute(DHT_CO.cube, co_sig, co_sig, 50)

co_mtdt = {}
co_mtdt['data_unit'] = u.Jy

co_cat = ppv_catalog(dco, co_mtdt)
co_cluster_ob = SpectralCloudstering(dco, co_cat, DHT_CO.data_header)

co_cluster_nums = np.array(co_cluster_ob.clusters)
cluster_label_cube = co_cluster_ob.clusters_asgn.data

print("dendrogram, clusters done")
#--------------------------------------------------------------------------------------------------------------------
@njit
def cluster_masked_list (cnums_list, data):
    x,y,z = data.shape
    print('xyz\n\n')
    N = len(cnums_list)
    print('N')
    output = np.zeros((N,x,y,z))
    print('output initiated\n\n')
    for n, cnum in enumerate(cnums_list):
        print('in for loop', n, cnum)
        msk = (data == cnum)
        print('made mask')
        mskd = data * msk
        print('timesed mask & data')
        output[n] = mskd
        print('sent to output')

    return output


@njit
def weighted_coords(data, vbl_axes, ax_ress): #prob very imprecise, since the clusters are spread out in l, b, v space
    #might need to do multiple distance calculations???

    if len(data.shape)==3: #then it's just one cube
        cb = data
        v,b,l = dc.weighted_vbl(cb,vbl_axes,ax_ress)

        output = np.array(v,b,l)

    else: #a list of cubes
        N = len(data)

        IDs, ls, bs, vs = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
        
        for i in range(N):
            cb = mskd_list[i]
            
            v,b,l = dc.weighted_vbl(cb,vbl_axes,ax_ress)
            IDs[i] = i
            ls[i] = l
            bs[i] = b
            vs[i] = v


        output = np.array(IDs, vs, bs, ls)
        
    return output


def generate_sources(colnames, data):
    #assume data is a list of arrays, where each array corresponds to a column from colnames
    
    if colnames[0][0] != '!':
        colnames[0] = '! '+ colnames[0]
        
    if len(colnames) < 6: #assume that only id, lon, lat, and vlsr info is given
        N = len(data[0])
        pfar = np.full(N, 0.5)
        v_unc = np.full(N, 0.7)
        
        data = data+[v_unc, pfar]
        colnames = colnames + ['v_unc', 'pfar']
    
    df = pd.DataFrame(np.transpose(data),columns = colnames)
    df.columns = [df.columns.to_list(),['!'] + [' ']*(len(df.columns)-1)] #add second row thats commented out
    
    df.to_csv('sources_info.inp',sep=' ',index=False)
    
    return df



@njit
def apply_clust_mask(masks,cube):
    output = np.zeros_like(masks)

    for i,m in enumerate(masks):
        msk = m/m #normalize because currently the values are the id number assigned to the cluster
        output[i] = cube*msk
        
    return output

print("1th functions defined\n\n")
#--------------------------------------------------------------------------------------------------------------------
ind_labels = cluster_masked_list(co_cluster_nums,cluster_label_cube) #isolate each cluster
print('cluster masked list generated \n\n')
co_clusters = apply_clust_mask(ind_labels, DHT_CO.cube) #apply the cluster masks to the data
print('apply cluster mask ... applied\n\n')
ids, vs, bs, ls = weighted_coords(co_clusters, \
    [DHT_CO.vels, DHT_CO.lats, DHT_CO.lons], [DHT_CO.dv.value, DHT_CO.db.value, DHT_CO.dl.value])

print('weighted coords calculated\n\n')

columns = ['! ID', 'l(deg)', 'b(deg)', 'v_lsr(km/s)']

clusters_df = generate_sources(columns, [ids, ls, bs, vs])
print("sources_info made")
#cleaning up the df so it is ready to be merged with the output from BDC

#get rid of the second header that was added 
clusters_df.columns = [clusters_df.columns.to_list()[i][0] for i in range(len(clusters_df.columns))]

#recast data as floats bc for some reason, generate_sources casts the data in the columns as str
clusters_df[clusters_df.columns] = clusters_df[clusters_df.columns].to_numpy().astype('float64')

#round stuff to the same number of decimal pts as summary.prt
clusters_df[['lon','lat']] = round(clusters_df[['lon','lat']],2)
clusters_df['vlsr'] = round(clusters_df['vlsr'],1)

clusters_df.to_csv('DHT_decomp.csv',index=False, sep=' ')

##and then come back to calculate weighted l b & v
## then calculate distance from here (try  both bayesian method & velocity gradient defined as y=mx+b method)
#then from there, calc masses, ares, compare with known values, do accumulation radii 

print("all done")
#--------------------------------------------------------------------------------------------------------------------




