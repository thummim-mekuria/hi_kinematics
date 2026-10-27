#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astrodendro import Dendrogram
from astropy.io import fits


import astropy.units as u
import astropy.constants as c


nice_catalog = pd.read_csv('lfs_w_simple_spectra_nice_r_acc.csv')


#add volume/surface density as one of the physical properties that are a part of the table
#build a way to see multiple thumbnails of spectra together, to more easily compare them


rho = (big_catalog['mass<kg>']*u.kg.to(c.M_sun) / (4/3 * np.pi * (big_catalog['effective_radius<pc>'])**3)).to_numpy()


plt.hist(rho,bins=100)
plt.hist(rho_w,density=True)
plt.xlim(0,0.2)
plt.vlines(0.025,0,70,color='red')
plt.xlabel(r'$M_\odot / pc^3$')
plt.show()


plt.scatter(rho,big_catalog['effective_radius<pc>'])


np.where(big_catalog._idx==260)


1e-6*1.594765e38*u.kg.to(c.M_sun)


(((40.5*1e3*c.M_sun)/(2*np.pi*0.025*c.M_sun/u.pc**3))**(1/3)).to(u.pc)


(((1.2*1e3*c.M_sun)/(2*np.pi*0.025*c.M_sun/u.pc**3))**(1/3)).to(u.pc)


#keep track of how many clouds there are at each stage of elimination


((m/(2*np.pi*rho_HI))**(1/3)).to(u.pc)


np.sqrt((23*(144)**2)/np.pi)


big_catalog.iloc[131]


rho_w = (weird_r['mass<kg>']*u.kg.to(c.M_sun) / (4/3 * np.pi * (weird_r['effective_radius<pc>'])**3)).to_numpy()


weird_r = big_catalog[big_catalog['accumulation_radius<deg>']<big_catalog['effective_radius<deg>']]


big_catalog = pd.read_csv('physcl_prprts_catalog.csv')


mekcat = pd.read_csv('Mekuria_catalog.csv',skiprows=2)


n_comps = pd.read_csv('hi_spectrum_components.csv')


nice_catalog


big_cat = big_catalog.merge(n_comps,on='_idx')


big_cat.to_csv('all_lvs_all_info.csv',index=False)


big_cat[big_cat['n_comps_parent']==1]


# Aquila_S [248, 249]
# California [173]
# Cam [182]
# Chamaeleon [218]
# CMa_OB1 [306]
# Corona_Australis [275]
# Crossbones [304]
# Gem_OB1 [208]
# Hercules [257]
# Lacerta [152]
# Lupus [271]
# Maddalena [311, 313]
# Mon_OB1 [189]
# Mon_R2 [301]
# Ophiuchus [242, 260]
# Orion_A [246]
# Orion_B [250]
# Orion_Lam []
# Perseus [214]
# Polaris [157]
# Rosette [226]
# Serpens_AqR []
# Taurus [264]


nice_solar = mekcat.iloc[[0,2,3,7,9,13,15,16,18,19]] 
#these are the ones that are in the catalog with: 1 gaussian component + r_gmc < r_hi


# maaaaybe bring back:
# -california
# -taurus
# -hercules
# -


nice_solar


all_coords = big_cat[['lon<deg>','lat<deg>']].to_numpy()


def minimize_ds (catalog,options):
    matches = []
    for cloud in catalog.values:
        l,b = cloud[1],cloud[2]
        ds = [np.sqrt((l-c[0])**2 + (b-c[1])**2) for c in options]
#         print(ds,'\n\n\n')
        c_in_cat = np.argmin(ds)
        matches.append((cloud[0],c_in_cat))
        
    return matches
 

solar_matches = np.array(minimize_ds(mekcat,all_coords))


these_idx = np.transpose(solar_matches)


solar_neighbors = big_cat.iloc[these_idx[1]]


solar_neighbors['name'] = these_idx[0]


solar_neighbors.columns = ['_idx', 'lon<deg>', 'lat<deg>', 'vel<km/s>', 'dist<kpc>',
       'area_sky_GMC<deg2>', 'mass_GMC<kg>', 'area_GMC<pc2>', 'radius_GMC<pc>',
       'radius_GMC<deg>', 'radius_parent<pc>',
       'radius_parent<deg>', 'vel_dispersion_GMC<km/s>',
       'n_comps_parent', 'name']


solar_neighbors = solar_neighbors[['name','_idx', 'lon<deg>', 'lat<deg>', 'vel<km/s>', 'dist<kpc>',
       'radius_GMC<deg>', 'radius_GMC<pc>', 'area_sky_GMC<deg2>',
       'area_GMC<pc2>', 'mass_GMC<kg>', 'vel_dispersion_GMC<km/s>',
       'radius_parent<deg>', 'radius_parent<pc>','n_comps_parent']]


solar_neighbors[['lon_z<deg>','lat_z<deg>']] = mekcat[['l_c<deg>','b_c<deg>']].to_numpy()


solar_neighbors


solar_neighbors.to_csv('solar_neighbors.csv',index=False)
