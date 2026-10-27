# %%time

import datacube as dc
import clouds as cld
import kd as kd

import numpy as np
# %matplotlib

import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram
from astropy.visualization import ZScaleInterval as zscale

from matplotlib.patches import Rectangle

zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf

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
        ns = ['03','06','09','12','16','19','21','24','26','27','34','37']
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






#=============================================================================================

# DHT_CO.lats[201],DHT_CO.lats[281]
# DHT_CO.lons[641],DHT_CO.lons[2241]
wo_gal_cent = np.zeros_like(DHT_CO.cube)
wo_gal_cent[:,:,:641] = DHT_CO.cube[:,:,:641]
wo_gal_cent[:,:,2241:] = DHT_CO.cube[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = DHT_CO.cube[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = DHT_CO.cube[:,281:,641:2241]


#=============================================================================================

d_wocent = Dendrogram.compute(wo_gal_cent,5*0.2,5*0.2,100)

#=============================================================================================

num_tings = len(d_wocent)
Ms = np.zeros(num_tings)
ds = np.zeros(num_tings)
As = np.zeros(num_tings)
zs = np.zeros(num_tings)

#=============================================================================================

def mass_and_area_kd(surface_density_map, area_map, d):
    """
    the first two help with retrieving the distance to the cloud (this is)
    """
    d=d[0]
    dAphys = (d**2 * area_map.value * area_map.unit.to(u.rad**2)) #to get rid of 'rad' as a unit
    
    A_total =  np.nansum(dAphys)
    M_total = np.nansum(surface_density_map * dAphys)
    M_total_kg = M_total.to(u.kg)


    return M_total_kg, A_total

#=============================================================================================



for t, tr in enumerate(d_wocent.trunk[:10]):
    print(t)
    m = tr.get_mask()
    this_cube = m*wo_gal_cent

    this_clob = DHT_CO.clone()
    this_clob.cube = this_cube


    dSig, dA = cld.surface_density(this_clob, 'CO',show_Sigma=False)
    v,b,l = tr.smallest_index
    d = kd.calculate_kd(this_clob.lons[l],this_clob.lats[b],this_clob.vels[v])

    M, A = mass_and_area_kd(dSig, dA,d)
    #the first time this goes through, save the units
    if t == 0:
        Ms = Ms*M.unit
        As = As*A.unit
        ds = ds*d[0].unit
    
    Ms[t] = M
    As[t] = A
    ds[t] = d

print(ds)
#     z = np.array(d.value) * this_clob.lats[b]*u.deg.to(u.rad) #height above the midplane
#     zs[t] = z

# colnames = ['M ['+Ms.unit.to_string()+']', 'A ['+As.unit.to_string()+']', 'd ['+ds.unit.to_string()+']']
# df = pd.DataFrame({colnames[0]:Ms, colnames[1]:As,colnames[2]:ds})


# df.to_csv('wo_cent_MsAsds_from_dend.csv',index = False)



