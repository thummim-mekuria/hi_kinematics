#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astrodendro import Dendrogram
import astropy.units as u
import astropy.constants as c

from astropy.io import fits

import datacube as dc
import math


hi4pi = fits.open('HI4PI_flatcalib.fits')[0]
dht = fits.open('DHT_CO.fits')[0]
struts = fits.open('structures_w_ds-Copy1.fits')[0]


nice_cat = pd.read_csv('lfs_w_simple_spectra_nice_r_acc.csv')


struts.shape


plt.imshow(struts.data[35],cmap=dc.cmap_grey0('viridis'))


nice_cat


def decomp_hi(ind,catalog=nice_cat,survey=hi4pi):
    this = catalog.iloc[ind]
#     name = this['name']
    idx = this['_idx']
    
    ds = 0.0833
    dv = 1.29
    
    
    l = math.ceil(this['l_cen_hi4pi<px>'])
    b = math.ceil(this['b_cen_hi4pi<px>'])
    v = math.ceil(this['v_cen_hi4pi<px>'])
    
    vsig = this['radius_parent<deg>']
    ra = this['radius_parent<deg>']
    
    vsig = math.ceil(5*vsig/dv)
    r = math.ceil(2.5*ra/ds) 
    
    vlen,blen,llen = hi4pi.data.shape
    
    lmn = max(0,l-r)
    lmx = min(llen-1,l+r)
    
    bmn = max(0,b-r)
    bmx = min(blen-1,b+r)
    
    vmn = max(0,v-vsig)
    vmx = min(vlen-1,v+vsig)
    
    this_cube = survey.data[vmn:vmx,bmn:bmx,lmn:lmx]
    rms = np.nanstd(this_cube) 
    dend = Dendrogram.compute(this_cube,rms,3*rms,(v*r**2)/250)
    lvs = [dend[i] for i in range(len(dend)) if dend[i].is_leaf]
    return lvs,[v,vsig,b,l,r]


nice_cat.columns


nice_cat[10*nice_cat['vel_dispersion_parent<km/s>']/1.29 > nice_cat['v_cen_hi4pi<px>']]


nice_cat[5*nice_cat['radius_parent<deg>']/0.0833 > nice_cat['l_cen_hi4pi<px>']]


nice_cat[5*nice_cat['radius_parent<deg>']/0.0833 > nice_cat['b_cen_hi4pi<px>']]


class HIC(cube,bounds,dend,survey=hi4pi):
    def __init__(self,cube,survey,bounds,dend)
        self.data = cube
        self.survey = survey
        vc, ch, bc, lc, r = bounds
        v0_ind = vc - ch
        b0_ind = bc - r
        l0_ind = lc - r
        
        sdv,sdb,sdl = survey.header['CDELT1'],survey.header['CDELT2'],survey.header['CDELT3']
        sv0,sb0,sl0 = survey.header['CRVAL1'],survey.header['CRVAL2'],survey.header['CRVAL3']
        sl0 = sl0*-1
        
        self.v0 = v0_ind*dv + sv0 
        self.b0 = b0_ind*db + sb0
        self.l0 = l0_ind*dl + sl0
        
        self.vc = vc
        self.bc = bc
        self.lc = lc
        
        lf_ms = [dend[i].get_mask() for i in range(len(dend)) if dend[i].is_leaf]
        self.leaf_masks = np.array(lf_ms)
        
        
    def __distance(self,coords1,coords2):
        x1,y1,z1 = coords1
        x2,y2,z2 = coords2
        
        r_sqd = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**3
        return np.sqrt(r_sqda)
        
    def match_GMC2leaf(self,gmc_coords,gmc_survey=dht):
        v_gmc, b_gmc, l_gmc = gmc_coords
        
        dv,db,dl = gmc_survey.header['CDELT1'],gmc_survey.header['CDELT2'],gmc_survey.header['CDELT3']
        
        
        #do something with the dendrogram?


(l-l0)/dl = x


len(nice_cat)


ora_parent


for l in ora_parent[0]:
    m = l.get_mask()
    v,vsig,b,l,r = ora_parent[1]
    
    vlen,blen,llen = hi4pi.data.shape
    
    lmn = max(0,l-r)
    lmx = min(llen-1,l+r)
    
    bmn = max(0,b-r)
    bmx = min(lblen-1,b+r)
    
    vmn = max(0,v-vsig)
    vmx = min(vlen-1,v+vsig)
    
    cb = m*hi4pi.data[vmn:vmx,bmn:bmx,lmn:lmx]
    w = np.nansum(cb,axis=0)
    plt.imshow(w,origin='lower',cmap=dc.cmap_grey0('viridis'))
    plt.show()


# #rewrite better so it is more of a plotty function
# def cloud_0th_moment(cloud,cat,data,co_or_hi):
#     if type(cloud)==int: #id
#         ind = np.where(cat['_idx']==cloud)[0]
        
#     elif type(cloud)==str:
#         ind = np.where(cat['name']==cloud)[0]
        
    
#     this = cat.iloc[ind]
# #     name = this['name']
#     idx = this['_idx']
    
#     if co_or_hi == 'CO':
#         cloud_type = 'GMC'
#         survey = 'dht'
#         ds = 0.125
#         dv = 1.3
        
#         l0 = 180.125
#     elif co_or_hi == 'HI':
#         cloud_type = 'parent'
#         survey = 'hi4pi'
#         ds = 0.0833
#         dv = 1.29
        
#     lcent = this['x_cen_'+survey+'<px>'].to_numpy()[0]
#     bcent = this['y_cen_'+survey+'<px>'].to_numpy()[0]
#     vcent = this['v_cen_'+survey+'<px>'].to_numpy()[0]
    
#     vsig = this['radius_'+cloud_type+'<deg>'].to_numpy()[0]
#     ra = this['radius_'+cloud_type+'<deg>'].to_numpy()[0]
    
#     v = math.ceil(5*vsig/dv) #5x sigma velocity range 
#     r = math.ceil(2.5*ra/ds)
    

#     this_cube = data[vcent-v:vcent+v,bcent-r:bcent+r,lcent-r:lcent+r]
    
#     print(vcent,v,bcent,r,lcent)
    
# #     if co_or_hi == 'HI':
# #         d = Dendrogram.compute(this_cube,np.nanstd(this_cube),3*np.nanstd(this_cube),(v*r**2)/250)
# #     return d
# #     w = np.nansum(this_cube*dv,axis=0)
    
# #     plt.imshow(w,origin='lower',cmap=dc.cmap_grey0('viridis'))
    
# #     plt.show()


# okay from here what i want to do is 
# --the cloud object will have: the cube, a pandas.series containing the info, and option to add more info (eg name if its known)
#     for the parent hi, there needs to be a step that dendrogra decomposes a cube obtained from the ra & vsig values
#     and then match the most relevant leaf to the GMC (can first do by eye to check that there will always/usually be 
#     a leaf to straightforwardly match (and if not, i think just kick out instead of fighting it) thenn use that mask to create the cube
#         
# --also for the cloud, contain coords in some way: db,dl,dv from survey 
#         &  l0,b0,v0 based on the indexing (might also move the x_cen etc info here?)
# --then there will aslo be functions, such as quick plotting (zeroth moment etc)
# --also important !! figure out how to save all of this info (the clobs, complexes) so i can just read them
#         instead of computing dendrograms blah blah all over again each time
#         (i think the way to save it would be to save a primary hdu with data=cube, header=strings that say what
#         files to import (dendrogram file name, info file name) 
#         
#      ^this might need be 2 or 3 functions used tg instead of one that tries to do everything (ie the rabbit hole for .load_from
#      
# oooo i think it would be better (instead of cloud & complex) to have GMC and HI objects (the GMC will be basically what i said ^ + a pointer to its parent) and the HI will also be similar, but will have an attribute 'has children' that affects how it is initialized(?) and a pointer to its child(ren)
# 
# and then maybe, one more big function that when called, loads all the objects (or you can tell it which objecs to load)


plt.hist(np.log10(mass_msun.value),bins=20,fill=False,edgecolor='blue')
plt.xlabel(r'log(mass/$M_\odot$)')
plt.show()


rs_co = u.pc*nice_cat['effective_radius<pc>']
rs_hi = u.pc*nice_cat['accumulation_radius<pc>']


plt.hist(rs_co,bins=20,fill=False,edgecolor='blue')
plt.hist(rs_hi,bins=20,fill=False,edgecolor='orange')

# plt.xlabel(r'log(mass/$M_\odot$)')
plt.show()


len(nice_cat)
