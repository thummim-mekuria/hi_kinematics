#!/usr/bin/env python
# coding: utf-8


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
from astrodendro import Dendrogram


lfcat = pd.read_csv('COGAL_decomposed_catalog_lvs_only.csv')
dend = Dendrogram.load_from('COGAL_dendrogram.fits')
src = pd.read_csv('sources_info.inp',delim_whitespace=True,skiprows=1)
sm = pd.read_csv('summary.prt',delim_whitespace=True,skiprows=2)


sm = sm.loc[1:]
sm = sm[['!','Long.','Lat.','+/-']]
sm.columns = ['lon','lat','vel','dist']


src.columns = ['_idx','lon','lat','vel','pfar','vunc']
src = src[['_idx','lon','lat','vel']]


sm['lon'] = sm['lon'].astype('float64')
sm['lat'] = sm['lat'].astype('float64')
sm['vel'] = sm['vel'].astype('float64')
sm['dist'] = sm['dist'].astype('float64')


src['_idx'] = src['_idx'].astype('int')
src['lon'] = round(src['lon'],2)
src['lat'] = round(src['lat'],2)
src['vel'] = round(src['vel'],1)


lvs = [dend[i] for i in range(len(dend)) if dend[i].is_leaf]


mg = src.merge(sm,on=['lon','lat','vel'])
have_ds= [i for i in range(179) if src['_idx'].loc[i] in mg['_idx'].to_list()]
lvs_w_ds = np.array(lvs)[have_ds]


DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')

DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 

smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole)
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641]
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

x,y = np.meshgrid(HI4PI.lons, HI4PI.lats)


big_leafy_cat = lfcat.merge(mg,on='_idx')


rsnbl_leafy_cat = big_leafy_cat[['_idx','lon','lat','vel','dist','area_exact']]


rsnbl_leafy_cat.columns = ['_idx','lon<deg>','lat<deg>','vel<km/s>','dist<kpc>','area_sky<deg2>']


rsnbl_leafy_cat


Ms = np.zeros(len(rsnbl_leafy_cat)) * u.kg


conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
mmw = 2 * 1.36 * c.m_p


rsnbl_leafy_cat.loc[132]


rsnbl_leafy_cat.loc[132]['effective_radius<pc>']


rsnbl_leafy_cat.loc[132]['lon<deg>']


lims = [DHT_CO.lons[0], DHT_CO.lons[-1], DHT_CO.lats[0], DHT_CO.lats[-1]]


lims


rsnbl_leafy_cat.loc[132]['_idx']


rsnbl_leafy_cat.loc[132]['effective_radius<deg>']


rsnbl_leafy_cat.loc[132]['accumulation_radius<deg>']


rsnbl_leafy_cat.loc[132]['total_mass<kg>']


(4.85e41 * u.kg/(2*np.pi*rho_HI))**(1/3)


oldRAs[132]


rsnbl_leafy_cat.loc[132]['accumulation_radius<pc>']


Ms[132]


oldMs[132]


plt.imshow(np.nansum(dend[264].get_mask(),axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)

plt.vlines(-19.1-3.7, 8.12-3.7, 8.12+3.7,color='red')
plt.vlines(-19.1+3.7, 8.12-3.7, 8.12+3.7,color='red')

plt.hlines(8.12-3.7, -19.1+3.7,-19.1-3.7,color='red')
plt.hlines(8.12+3.7, -19.1+3.7,-19.1-3.7,color='red')

plt.vlines(-19.1-40.2, 8.12-40.2, 8.12+40.2,color='orange')
plt.vlines(-19.1+40.2, 8.12-40.2, 8.12+40.2,color='orange')

plt.hlines(8.12-40.2, -19.1+40.2,-19.1-40.2,color='orange')
plt.hlines(8.12+40.2, -19.1+40.2,-19.1-40.2,color='orange')

# plt.xlim(0,-30)
# plt.ylim(0,15)


M.to(c.M_sun)


d/(423*u.pc)


a


np.nanmax(W)


pd.read_csv('summary.prt').loc[2]


conversion_factor


mold.to(c.M_sun)


M.to(c.M_sun)/1e3


M.to(c.M_sun)


M.to(c.M_sun)


M.to(c.M_sun)


M.to(c.M_sun)


-159+360


rsnbl_leafy_cat['effective_radius<pc>'].loc[115]


d


DHT_CO_unip = dc.unpack_fits_cube('COGAL_local_raw.fits','DHT')


DHT_CO_unip = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


plt.imshow(N.value,cmap=dc.cmap_grey0('viridis'))
plt.xlim(2600,2800)
plt.ylim(50,255)
plt.colorbar()


aangold*(930*u.pc)**2


np.sqrt(5775.62/np.pi)


dir(dend[244])


conversion_factor


a.to(u.rad**2).value*(450*u.pc)**2


a/(0.125**2)


W.shape


481*2881


ar = (np.nan_to_num(W))==0


np.nansum(ar)


1384358+1403==1385761


M.to(c.M_sun)


0.25*(np.nansum(area_sig)*((DHT_CO.db*DHT_CO.dl).to(u.rad**2)).value*d**2).to(u.pc**2)


(0.8e5*c.M_sun)/(1500*u.pc**2)


np.nansum(dSig).to(u.kg/u.pc**2)


for i in [115]:#range(len(rsnbl_leafy_cat)):
    idx = rsnbl_leafy_cat['_idx'][i]
    lf_mask = dend[idx].get_mask()
    dat = dend.data
    
    a = rsnbl_leafy_cat['area_sky<deg2>'].loc[i] *u.deg**2
    d = rsnbl_leafy_cat['dist<kpc>'].loc[i] *u.kpc
    
    masked_leaf = dat*lf_mask
    
    clob = DHT_CO.clone()
    clob.cube = masked_leaf
    
    dsigold, aangold = cld.surface_density(clob,'CO')
    mold,aold = cld.mass_and_area(dsigold,aangold,d)
    
    
    W = np.nansum(masked_leaf*u.K*DHT_CO.dv,axis=0)
    N = W*conversion_factor
    dSig = N*mmw
    
    area_sig = (W!=0)
    A = np.nansum(area_sig) * ((DHT_CO.db*DHT_CO.dl).to(u.rad**2)).value * (d**2)
    
    M = np.nansum(dSig)*A
    
#     Sig = np.nansum(dSig)
#     M = Sig * (a.to(u.rad**2)).value * (d)**2
    
#     Ms[i] = M
    plt.imshow(np.log10(W.value),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims)
#     plt.imshow(np.nansum(masked_leaf,axis=1),origin='lower',cmap=dc.cmap_grey0('viridis'))#,extent=lims)

    plt.xlim(-145,-160)
#     plt.xlim(2600,2800)
    plt.ylim(-20,-5)
#     plt.ylim(20,45)
    plt.colorbar()
    plt.show()
# #     print(np.nanmax(dSig.to(c.M_sun/u.pc**2)))
#     print(Sig.to(c.M_sun/u.pc**2) * (a.to(u.rad**2)).value * (d)**2)
    

for i in [115]:#range(len(rsnbl_leafy_cat)):
    idx = rsnbl_leafy_cat['_idx'][i]
    lf_mask = dend[idx].get_mask()
#     dat = dend.data
    dat=DHT_CO_unip.cube
    
    a = rsnbl_leafy_cat['area_sky<deg2>'].loc[i] *u.deg**2
    d = rsnbl_leafy_cat['dist<1000pc>'].loc[i] *1e3*u.pc
    
    masked_leaf = dat*lf_mask
    
    clob = DHT_CO.clone()
    clob.cube = masked_leaf
    
    dsigold, aangold = cld.surface_density(clob,'CO')
    mold,aold = cld.mass_and_area(dsigold,aangold,d)
    
    
    W = np.nansum(masked_leaf*u.K*DHT_CO.dv,axis=0)
    N = W*conversion_factor
    dSig = N*mmw
    
    Sig = np.nansum(dSig)
    M = Sig * (a.to(u.rad**2)).value * (d)**2
    
#     Ms[i] = M
    plt.imshow(dSig.value,origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims)
    plt.xlim(-145,-160)
    plt.ylim(-20,-5)
    plt.colorbar()
    plt.show()
# #     print(np.nanmax(dSig.to(c.M_sun/u.pc**2)))
#     print(Sig.to(c.M_sun/u.pc**2) * (a.to(u.rad**2)).value * (d)**2)
    

rsnbl_leafy_cat


cocat_old = pd.read_csv('co_catalog.csv',skiprows=26,delim_whitespace=True)


oldRAs = (cocat_old['R_A'].to_numpy()*u.pc)


oldMs = (cocat_old['M'].to_numpy()*u.kg)


plt.hist(np.log10(Ms/oldMs),bins=50)
# plt.loglog()


rsnbl_leafy_cat['total_mass<kg>'] = Ms.value


Rs = np.zeros(len(rsnbl_leafy_cat)) * u.pc


for i, ar in enumerate(rsnbl_leafy_cat['area_sky<deg2>']):
    d = rsnbl_leafy_cat['dist<1000pc>'].loc[i]  *1e3*u.pc

    ar_phys = ar*(u.deg**2).to(u.rad**2) * (d**2)
    r_phys = np.sqrt(ar_phys/np.pi)
    
    Rs[i] = r_phys


rsnbl_leafy_cat['effective_radius<pc>'] = Rs.value


rho_HI = 0.025*c.M_sun/(u.pc**3)


RAs = np.zeros(len(rsnbl_leafy_cat)) * u.pc


plt.hist(np.log10(Ms.to(c.M_sun).value),bins=50)
plt.yscale('log')


RAs


for i, m in enumerate(rsnbl_leafy_cat['total_mass<kg>']):
    r_a = (m*u.kg/(2*np.pi*rho_HI))**(1/3)
    
    RAs[i] = r_a


rsnbl_leafy_cat['accumulation_radius<pc>'] = RAs.value


ds = rsnbl_leafy_cat['dist<1000pc>'].to_numpy()*1e3*u.pc


bs = rsnbl_leafy_cat['lat<deg>'].to_numpy()*u.deg


acc_rad_ang = (RAs/ds)*u.rad.to(u.deg) * np.cos(bs.to(u.rad))


eff_rad_ang = (Rs/ds)*u.rad.to(u.deg) * np.cos(bs.to(u.rad))


rsnbl_leafy_cat['effective_radius<deg>'] = eff_rad_ang


rsnbl_leafy_cat['accumulation_radius<deg>'] = acc_rad_ang


parents = np.zeros_like(HI4PI.cube,dtype=bool)

for i in range(len(rsnbl_leafy_cat)):
    l = rsnbl_leafy_cat['lon'][i]*u.deg
    b = rsnbl_leafy_cat['lat'][i]*u.deg
    v = rsnbl_leafy_cat['vel'][i]*u.km/u.s 
    r = rsnbl_leafy_cat['accumulation_radius<deg>'][i]*u.deg
    
    vrms = vrms_3x[i]*DHT_CO.dv #DHT px -> km/s
    parents = cld.add_accumulation_prism(v,b,l,r,vrms,parents)

# parents = HI4PI.cube * parents


whi = np.nansum(HI4PI.cube,axis=0)


rsnbl_leafy_cat.columns


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
ax.imshow(whi,extent=lims,origin='lower')


for l,b,rbb,r_big in \
zip(rsnbl_leafy_cat['lon<deg>'],rsnbl_leafy_cat['lat<deg>'],\
    rsnbl_leafy_cat['effective_radius<deg>']\
    ,rsnbl_leafy_cat['accumulation_radius<deg>']):    
    
    ax.add_patch(Circle((l,b),rbb,fill=False,edgecolor='white'))
#     ax.add_patch(Circle((l,b),r_big,fill=False,edgecolor='red'))

#         print(ind)
# ax.set_xlim(180,130)


# fig,ax = plt.subplots(1,1,figsize=(15,10))

# lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]
# ax.imshow(whi,extent=lims)


# for l,v,r,vrms in zip(cocat['lon'],cocat['vel'],ang_R_A,4*vrms_3x):
#     ax.add_patch(Ellipse((l,v),2*r,vrms,fill=False,edgecolor='white',alpha=0.5))
    
    
# ax.set_ylim(-50,50)


# fig,ax = plt.subplots(1,1,figsize=(15,10))

# lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
# ax.imshow(whi,extent=lims)

# for l,b,r in zip(cocat['lon'],cocat['lat'],ang_R_A):
#     ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5))
# #     print(l,b,r)
    
    
# # ax.set_xlim(180,130)


# fig,ax = plt.subplots(1,1,figsize=(15,10))
# lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]

# ax.imshow(np.nansum(parents,axis=0),extent=lims,origin='lower')
# for l,b,r in zip(cocat['lon'],cocat['lat'],ang_R_A):
#     ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5,linewidth=2))
# #     print(l,b,r)


# fig,ax = plt.subplots(1,1,figsize=(15,10))
# lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]

# ax.imshow(np.nansum(parents,axis=1),extent=lims,origin='lower')

# ax.set_ylim(-50,50)
# for l,v,r,vrms in zip(cocat['lon'],cocat['vel'],ang_R_A,vrms_3x):
#     #everything else is in physical units so convert vrms from px -> km/s
#     vrms = vrms*1e-3*HI4PI.dv.value
#     ax.add_patch(Ellipse((l,v),2*r,4*vrms,\
#                          fill=False,edgecolor='white',alpha=0.5,linewidth=2))
    
# # ax.set_xlim(80,50)
    

import cc3d
parents_masked_overlap,N = cc3d.connected_components(parents,connectivity=26,return_N = True)
#connectivity=6 and 18 also give the same N


from numba import njit

@njit
def separate_clouds(labels_map,data,N):
    x,y,z = data.shape
    out = np.zeros((N,x,y,z))
    for n in range(1,N+1):
        m = (labels_map==n)
        mskd = m*data
        out[n-1] = mskd
        
    return out


# #come back and figure out whats wrong with this / the function ^ 
# distinct_hi = separate_clouds(parents_masked_overlap,HI4PI.cube,N)


# lims = [HI4PI.lons[0], HI4PI.lons[-1], HI4PI.lats[0], HI4PI.lats[-1]]
# for i in range(N):
#     plt.imshow(np.nansum(HI4PI.cube*distinct_hi[i],axis=0),\
#                cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)
# #     plt.colorbar()
#     plt.title(str(i+1))
#     plt.show()


# lims = [HI4PI.lons[0], HI4PI.lons[-1], 1e-3*HI4PI.vels[0], 1e-3*HI4PI.vels[-1]]
# for i in range(N):
#     plt.imshow(np.nansum(HI4PI.cube*distinct_hi[i],axis=1),\
#                cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)
# #     plt.colorbar()
#     plt.title(str(i+1))
#     plt.show()


# #once i figure out getting distinct_hi, save it as a fits
# HI_parents = fits.PrimaryHDU(distinct_hi)
# HI_parents.writeto('decomp_HI_parents_frm_RA.fits',overwrite=True)
