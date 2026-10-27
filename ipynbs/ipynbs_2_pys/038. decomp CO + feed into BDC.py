#!/usr/bin/env python
# coding: utf-8


# %%time

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
from astrodendro import Dendrogram
from astropy.visualization import ZScaleInterval as zscale


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')


DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')
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
        ns = ['03','06','09','12','16','19','21','24','26','27','34','37','18']
        #for now, figure out DHT01 (bad data) 
        #& DHT18 (two ceph ds, beeg) 
        #& DHT08 (beeg) 
        #& DHT24, DHT34 (lam or & chamaeleon are disjointed down to the really low sig levels)
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output

#----------------------

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


wo_gal_cent = np.zeros_like(smooth_whole)
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641]
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]


np.nanstd(smooth_whole),np.nanstd(wo_gal_cent)


# d = Dendrogram.compute(smooth_whole, 0.2,0.2, 50)


npx_wo_cent2 = [d2.trunk[i].get_npix() for i in range(len(d2.trunk))]
peaks2 = [d2.trunk[i].get_peak()[1] for i in range(len(d2.trunk))]


# npx_wcent = [d.trunk[i].get_npix() for i in range(len(d.trunk))]


# peaks = [d.trunk[i].get_peak()[1] for i in range(len(d.trunk))]


# plt.hist(np.log(peaks),bins=50,label='peaks val,wcent',density=True,alpha=0.7)

plt.hist(np.log(peaks2),bins=50,fill=0,label='peaks val,wocent',density=True)

plt.xlabel('log[N_pxs] / log[T_peak]')
plt.legend()
plt.show()


# plt.hist(np.log10(npx_wcent),bins=50,label='pixel count,wcent',density=True,alpha=0.7)

plt.hist(np.log10(npx_wcent2),bins=50,label='pixel count,wocent',fill=False,density=True)


plt.xlabel('log[N_pxs] / log[T_peak]')
plt.legend()
plt.show()


vels,lats,lons = DHT_CO.vels,DHT_CO.lats,DHT_CO.lons
dv,db,dl = DHT_CO.dv.value,DHT_CO.db.value,DHT_CO.dl.value


masks = [d2.trunk[i].get_mask() for i in range(len(d2.trunk))]


@njit
def calc_trunk_coords(cube, masks_list, axes, ress):
    N = len(masks_list)
    IDs, ls, bs, vs = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    
    for i,m in enumerate(masks_list):
        cb = m*cube
        
        v,b,l = dc.weighted_vbl(cb,axes,ress)
        IDs[i] = i
        ls[i] = l
        bs[i] = b
        vs[i] = v
        
    return IDs, vs, bs, ls


ids = ids.astype(int).astype('str')


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


def fix_colnames(filename,num):
    df = pd.read_csv(filename,skiprows=num,delim_whitespace=True,header=None)
    head = pd.read_csv(filename,nrows=num-1).loc[1].to_list()[0].split()[1:] 
    #asssume the format of the file is summary.prt from BDC
    
    df.columns = head
    
    return df


sm = fix_colnames('summary.prt',4)


src[src.columns] = src[src.columns].to_numpy().astype('float64')


src[['lon','lat']] = round(src[['lon','lat']],2)


src['vlsr'] = round(src['vlsr'],1)


src.columns = [src.columns.to_list()[i][0] for i in range(len(src.columns))]


src.columns.to_list()[1:4] = sm.columns[0:3]


new_cols = src.columns[1:4].to_list()+sm.columns[3:8].to_list()+\
           (' '.join(sm.columns[8:].to_list()).upper()).split() #want to differentiate column names between d1&d2


sm.columns = new_cols


mg = sm.merge(src,on=['lon','lat','vlsr'])


ids_w_ds = mg['! id'].to_numpy().astype('int')


trunks_w_ds = np.array(d2.trunk)[ids_w_ds]


masks_w_ds = [trunks_w_ds[i].get_mask() for i in range(len(trunks_w_ds))]
masks_w_ds = np.array(masks_w_ds)


(1/60)*8 + (1/3600)*8


@njit
def apply_mask(masks,cube):
    output = np.zeros_like(masks)
    for i,m in enumerate(masks):
        output[i] = cube*m
        
    return output


plt.hist(np.log10(np.sqrt(As/np.pi)),bins=50)
plt.show()


plt.hist(np.log10(Ms*units['M'].to(c.M_sun)),bins=50)
plt.show()


plt.hist(np.log10(1e3*mg['Dist.']),bins=50)
# plt.xscale('log')
plt.show()


ndesc = [len(trunks_w_ds[i].descendants) for i in range(len(trunks_w_ds))]


plt.hist(ndesc,bins=100)
plt.vlines(10,0,1e3,color='red',linestyle='--')
plt.yscale('log')
# plt.xscale('log')
# plt.xlim(1e-5,1e3)
plt.show()


peaks_w_ds = np.array(peaks2)[ids_w_ds]


npx_w_ds = np.array(npx_wo_cent2)[ids_w_ds]


ndesc_long = [len(d2.trunk[i].descendants) for i in range(len(d2.trunk))]


leaves_w_ds = [trunks_w_ds[i].sorted_leaves() for i in range(len(trunks_w_ds))]
leaves_w_ds = np.array([leaf for row in leaves_w_ds for leaf in row])


leaf_masks_ds = [leaves_w_ds[i].get_mask() for i in range(len(leaves_w_ds))]
leaf_npix_ds = [leaves_w_ds[i].get_npix() for i in range(len(leaves_w_ds))]


leaf_peaks_ds = [leaves_w_ds[i].get_peak()[1] for i in range(len(leaves_w_ds))]


10**(0.25)


np.nanstd(wo_gal_cent)


np.log10(5*0.1168)


trunks_w_ds[1].newick


plt.scatter(np.log10(leaf_peaks_ds),np.log10(leaf_npix_ds),marker='.')
plt.hlines(2,-0.5,1.5,color='red',linestyle='--')
plt.vlines(-0.233587,1.75,3.5,color='red',linestyle='--')

plt.xlabel('peak T')
plt.ylabel('npix')
plt.show()


plt.hist(np.log10(leaf_peaks_ds),bins=50,alpha=0.7,edgecolor='grey')
plt.vlines(-0.233587,0,80,color='red',linestyle='--')
plt.show()


plt.hist(np.log10(leaf_npix_ds),bins=50,alpha=0.7,edgecolor='grey')
plt.show()


trunks_w_ds[5].sorted_leaves()


plt.scatter(npx_w_ds,ndesc,c=np.log10(peaks_w_ds),marker='x')# np.log10(Ms*units['M'].to(c.M_sun)))
# plt.scatter(npx_wo_cent2,ndesc_long,c=np.log10(peaks2),cmap='YlGn' ,marker='+')# np.log10(Ms*units['M'].to(c.M_sun)))

# plt.scatter(np.log10(Ms*units['M'].to(c.M_sun)),ndesc)
plt.hlines(10,1,1e6)
plt.vlines(1e2,0,1e3)

plt.yscale('log')
plt.xscale('log')
plt.show()


plt.scatter(peaks_w_ds,ndesc,c=np.log10(Ms*units['M'].to(c.M_sun)))
# plt.scatter(np.log10(Ms*units['M'].to(c.M_sun)),ndesc)
plt.hlines(10,1,8)
plt.yscale('log')
plt.xscale('log')
plt.show()


Ms[20:30]*units['M'].to(c.M_sun)/1e3


msk = dSig.value != 0.


plt.imshow(msk*dA.value,**kw_imshow)
plt.colorbar()


np.nansum(dA)


mg['Dist.'][20:30].to_numpy()*u.pc*1e3


u.pc**2*As[20]/(mg['Dist.'][20]*u.pc*1e3)**2


u.pc**2*As[21]/(mg['Dist.'][21]*u.pc*1e3)**2


As[21]


(mg['Dist.'][20]*u.pc*1e3)**2


(mg['Dist.'][21]*u.pc*1e3)**2


for mc in masked_cubes[20:30]:
    plt.imshow(np.nansum(mc,axis=0),**kw_imshow)
    plt.show()


id(masked_clobs[0])


id(masked_clobs[1])


id(masked_clobs[1])


k,l = np.array(([DHT_CO.clone()]*2))


id(k)==id(l)


k,l
