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


cocat = pd.read_csv('co_catalog.csv',delim_whitespace=True,skiprows=32)


pd.read_csv('co_catalog.csv',skiprows=4,nrows=28)


ang_R_A = (cocat['R_A']/(cocat['Dist.']*1e3) * np.cos(cocat['lat'] * u.deg.to(u.rad))) *\
            u.rad.to(u.deg) 


vrms_3x = 3*cocat['v_rms']


whi = np.nansum(HI4PI.cube,axis=0)


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


parents = np.zeros_like(HI4PI.cube,dtype=bool)

for i in range(len(cocat)):#[130,69,121,118,64]:
    l = cocat['lon'][i]*u.deg #add unit
    b = cocat['lat'][i]*u.deg #add unit
    v = cocat['vel'][i]*u.km/u.s #add unit
    r = ang_R_A[i]*u.deg #add unit
    vrms = vrms_3x[i]*DHT_CO.dv #DHT px -> km/s
    parents = cld.add_accumulation_prism(v,b,l,r,vrms,parents)

# parents = HI4PI.cube * parents


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


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
ax.imshow(whi,extent=lims,origin='lower')


for l,b,a_bb,r_A,ind in zip(cocat['lon'],cocat['lat'],cocat['area_exact'],ang_R_A,cocat.index):
    r_bb = np.sqrt(a_bb/np.pi)
#     r_A = r #1 degree -> 8pixels bc 0.125deg/pix
    if r_bb > r_A:
        ax.add_patch(Circle((l,b),r_A,fill=False,edgecolor='white'))
        ax.add_patch(Circle((l,b),r_bb,fill=False,edgecolor='red'))

#         print(ind)
# ax.set_xlim(180,130)


big_rbb =[1,5,31,83,88,94,98,99,110,111,123,124,125,132,150,151,152]


#accumulation radii
ras = (cocat['M'].to_numpy()[big_rbb]*u.kg /(2*np.pi*0.025*c.M_sun/(u.pc**3)))**(1/3)


#effective radii
rss = np.sqrt(cocat['area_exact'].to_numpy()[big_rbb]/np.pi)*u.deg.to(u.rad)\
* cocat['Dist.'].to_numpy()[big_rbb] * 1e3 * u.pc


ras < rss


ds_big_rbb = u.pc*1e3*cocat['Dist.'].to_numpy()[big_rbb]


rs_all = np.sqrt(cocat['area_exact'].to_numpy()/np.pi)*u.deg.to(u.rad)\
* cocat['Dist.'].to_numpy() * 1e3 * u.pc


rho_all = cocat['M'].to_numpy()*u.kg/(4/3*np.pi*rs_all**3)
rho_big_rbb = cocat['M'].to_numpy()[big_rbb]*u.kg/(4/3*np.pi*rss**3)


all_rhos_myA = cocat['M'].to_numpy()*u.kg/(4/3*np.pi*all_rs_with_myA**3)
big_rhos_myA = cocat['M'].to_numpy()[big_rbb]*u.kg/(4/3*np.pi*big_rs_with_myA**3)


for i in range(len(cocat)):
    if all_rs_with_myA[i] > cocat['R_A'][i]*u.pc:
        print(i)


big_rbb


plt.hist(np.log10((rho_all/c.M_sun.value).value),bins=50,color='xkcd:crimson',alpha=0.7,density=True)
plt.hist(np.log10((rho_big_rbb/c.M_sun.value).value),bins=5,color='xkcd:sage green',alpha=0.7,density=True)

plt.hist(np.log10(all_rhos_myA.value),density=True,bins=50,edgecolor='xkcd:crimson',fill=False)
plt.hist(np.log10(big_rhos_myA.value),bins=8,density=True,edgecolor='xkcd:sage green',fill=False)


plt.vlines(0.025,0,2)
# plt.yscale('log')
# plt.xscale('log')
plt.title(r'$\rho$ [$M_\odot/pc^3$]')
# plt.ylim(1e-2,1e2)
# plt.xlim(-5,5)
# plt.xlim(25,33)
plt.show()


plt.hist(1e-6*cocat['M'].to_numpy()*u.kg.to(c.M_sun),density=True,bins=90,color='xkcd:crimson',alpha=0.7)
plt.hist(1e-6*cocat['M'].to_numpy()[big_rbb]*u.kg.to(c.M_sun),bins=8,density=True,color='xkcd:sage green',alpha=0.7)
plt.yscale('log')
plt.title(r'masses [$10^6\ M_\odot$]')
plt.show()


plt.hist(rs_all,density=True,bins=50,color='xkcd:crimson',alpha=0.7)
plt.hist(rss,bins=8,density=True,color='xkcd:sage green',alpha=0.7)

plt.hist(all_rs_with_myA,density=True,bins=50,edgecolor='xkcd:crimson',fill=False)
plt.hist(big_rs_with_myA,bins=8,density=True,edgecolor='xkcd:sage green',fill=False)

plt.yscale('log')
plt.title('effective radii [pc]')
plt.show()


plt.hist(u.pc*cocat['R_A'].to_numpy(),density=True,bins=50,color='xkcd:crimson',alpha=0.7)
plt.hist(ras,bins=8,density=True,color='xkcd:sage green',alpha=0.7)
plt.title('accumulation radii [pc]')
plt.show()


plt.hist(u.pc*cocat['Dist.'].to_numpy(),density=True,bins=30,color='xkcd:crimson',alpha=0.7)
plt.hist(ds_big_rbb/1e3,bins=10,density=True,color='xkcd:sage green',alpha=0.7)
plt.title('distances [kpc]')
plt.show()


((1e3*u.pc*cocat['Dist.'].to_numpy())**2)*(cocat['area_exact'].to_numpy()*(u.deg**2).to(u.rad**2))


cocat['A'].to_numpy() * u.pc**2


cocat.index.to_numpy()


myAs = cocat['A'].to_numpy() * u.pc**2


all_rs_with_myA = np.sqrt((cocat['A'].to_numpy() * u.pc**2)/np.pi)


big_rs_with_myA = np.sqrt((cocat['A'].to_numpy()[big_rbb] * u.pc**2)/np.pi)


ppvAs = ((u.pc*1e3*cocat['Dist.'].to_numpy())**2)*\
            (cocat['area_exact'].to_numpy()*(u.deg**2).to(u.rad**2))


plt.scatter(cocat.index.to_numpy(), 100*(myAs-ppvAs)/myAs,marker='+')
plt.scatter(big_rbb, ((myAs-ppvAs)/myAs)[big_rbb],marker='+',color='red')

plt.yscale('log')
plt.title('diff in calculated area (%)')
plt.show()


big_rbb


lb_pairs_big_rbb = cocat[['lon','lat']].to_numpy()[big_rbb]


lims


for pair in lb_pairs_big_rbb:
    l,b=pair
    plt.imshow(whi,extent=lims)
    plt.ylim(b-5,b+5)
    plt.xlim(l-5,l+5)
    plt.scatter(l,b,marker='+',color='red')
    plt.show()


dend = Dendrogram.load_from('COGAL_dendrogram.fits')


dend


help(dend[0])


onelf = dend._structures_dict[308]


indices = onelf.indices(subtree=True)


indices2 = zip(*tuple(indices[i] for i in range(3) if i != 0))


arr = np.transpose(np.array(indices[1:]))


arr.sort()


len(set(indices2))*DHT_CO.db**2


ind_lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


lvs_w_ds = [ind_lvs[i] for i in range(172) if i in cocat['!id']]


onelf_in_cocat = cocat.iloc[-1]


(1e3*u.pc*onelf_in_cocat['Dist.'])**2*(onelf_in_cocat['area_exact']*(u.deg**2).to(u.rad**2))


onelf_in_cocat['A']


msk = (np.nansum(onelf.get_mask(),axis=0) != 0.)


np.sum(msk)


264*DHT_CO.db**2


cocat.iloc[-1]


msk2 = (np.nansum(ind_lvs[-1]*u.K*1.3*u.km/u.s,axis=0)*1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))) != 0.


np.nansum(np.ones_like(msk2) * msk2) * (DHT_CO.db**2).to(u.rad**2)


np.nansum(np.ones_like(msk2) * msk2) *(u.deg**2).to(u.rad**2)


DHT_CO.dl.to(u.rad) * DHT_CO.db.to(u.rad)


(DHT_CO.db**2).to(u.rad**2)


np.nansum((np.ones_like(msk2)*(DHT_CO.db**2)*msk2).to(u.rad**2) * (1e3*u.pc*onelf_in_cocat['Dist.']**2))


np.nansum((np.ones_like(msk2)*(DHT_CO.db**2)*msk2).to(u.rad**2)) * (1e3*u.pc*onelf_in_cocat['Dist.']**2)
