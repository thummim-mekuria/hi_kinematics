#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


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
from astropy.visualization import ZScaleInterval as zscale

from astrodendro import ppv_catalog
from astrodendro import Dendrogram


from scimes import SpectralCloudstering


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


# * [dendogram documentation](https://dendrograms.readthedocs.io/en/latest/using.html)
# * [scimes documentation](https://scimes.readthedocs.io/en/latest/tutorial.html)


# import orion_dendogram as od


# dend = od.compute()


# len(dend.trunk), len(dend)


# od.launch_viewer(dend)


# orion = dc.unpack_fits_cube('DHT27_Orion_interp.fits', 'DHT')


# im = np.nansum(orion.cube*orion.dv.value, axis = 0)
# vmin, vmax = zlims(im)

# imx, imy = orion.lons, orion.lats
# plt.imshow(im, vmin = vmin, vmax = vmax, \
#            extent = [imx[0], imx[-1], imy[0], imy[-1]],**kw_imshow)

# plt.colorbar()

# plt.show()


# #values from tutorial(https://scimes.readthedocs.io/en/latest/tutorial.html)
# #may need to change later
# dend_sigma = 0.3 #K, noise level
# dend_ppb = 1.3 #pixels/beam


# orion_dend = Dendrogram.compute(orion.cube,\
#                                min_value = 3*dend_sigma,\
#                                min_delta = 2*dend_sigma,\
#                                min_npix = 3*dend_ppb)


# plotter = orion_dend.plotter()


# lines = plotter.get_lines()


# fig, axs = plt.subplots(1, 1, figsize = (5, 5))
# axs.add_collection(lines)
# plt.show()


# fig, axs = plt.subplots(1, 1, figsize = (5, 5))
# plotter.plot_tree(axs)
# plt.show()


# dir(dend)


# dend.leaves[10].indices()


# dend.data.shape


# print(dend.structure_at((24,57,157)))


# for i in dend.trunk:
#     npx = i.get_npix()
#     m = i.get_mask()
# #     plt.imshow(np.nansum(m * dend.data, axis = 0), **kw_imshow)
#     plt.show()


# (1.09/557)*1e3


cogal = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')


# np.nanstd(cogal.cube)


len(dco), len(dco.trunk)


xs = np.arange(1.7,17,0.1)


ys = (a1*(xs/8.5)**a2 + a3)*240


ys2 = (xs*0.2 + (240-0.2*8.5))


plt.plot(xs,ys)
# plt.plot(xs,ys2)
plt.plot(xs,xs*1.2+233)


def d_Russeil17(l, b, vlsr, ThetaSun=240*u.km/u.s, Rsun = 8.34*1e3*u.pc):
    cosb = np.cos(b*u.deg.to(u.rad))
    cosl = np.cos(l*u.deg.to(u.rad))
    sinl = np.sin(l*u.deg.to(u.rad))
    vlsr = vlsr*u.km/u.s
    a2 = 0.0803
    pwr = 2*(1/(1-a2))
    coeff = 1.022
    
    factor = (coeff*ThetaSun*cosb*sinl/(vlsr+ThetaSun*cosb*sinl))
    
    R = Rsun * factor**pwr
    
    c = Rsun**2 - R**2
    a = cosb**2
    b = -2*Rsun*cosb*cosl
    
    term1 = -b
    term2 = b**2
    term3 = 4*a*c
    term4 = 2*a
    
#     print(coeff*ThetaSun*cosb*sinl,(vlsr+ThetaSun*cosb*sinl))#, term2, term3, term4)
    
    dfar = (term1 + np.sqrt(term2-term3))/term4
    dnear = (term1 - np.sqrt(term2-term3))/term4
    
    if R < Rsun:
        return (dnear, dfar)

    else:
        return dfar
        

def d_Mroz18(l, b, vlsr, ThetaSun=233.6*u.km/u.s, Rsun = 8.122*1e3*u.pc, M = -1.34*u.km/(u.s*1e3*u.pc)):
    cosb = np.cos(b*u.deg.to(u.rad))
    cosl = np.cos(l*u.deg.to(u.rad))
    sinl = np.sin(l*u.deg.to(u.rad))
    vlsr = vlsr*u.km/u.s
    
    B = ThetaSun-(M*Rsun)
        
    num = Rsun*sinl*cosb*B
    den = vlsr+ThetaSun*sinl*cosb
    const = num/den
        
    R = const/(1-M/B*const)
    
    c = Rsun**2 - R**2
    a = cosb**2
    b = -2*Rsun*cosb*cosl
    
#     term1 = -b
#     term2 = b**2
#     term3 = 4*a*c
#     term4 = 2*a
    
# #     print(coeff*ThetaSun*cosb*sinl,(vlsr+ThetaSun*cosb*sinl))#, term2, term3, term4)
    
#     dfar = (term1 + np.sqrt(term2-term3))/term4
#     dnear = (term1 - np.sqrt(term2-term3))/term4
    
    term1 = Rsun*cosl
    term2 = R**2
    term3 = (Rsun*sinl)**2
    
    dnear = term1 - np.sqrt(term2-term3)
    dfar = term1 + np.sqrt(term2-term3)
    
#     print(M/B)
    if R < Rsun:
        return dnear
#         return (dnear, dfar)

    else:
        return dfar
        

d_Mroz18(160,-19,-2)


xs = np.arange(-180,180,0.5)


y0pos = np.zeros_like(xs)
yneg10pos = np.zeros_like(xs)
ypos10pos = np.zeros_like(xs)
yneg10neg = np.zeros_like(xs)
ypos10neg = np.zeros_like(xs)


for x,xx in enumerate(xs):
    y1 = d_Mroz18(xx,0,5)
    if y1.value < 0:
        y0pos[x] = 1
    
    y2 = d_Mroz18(xx,10,5)
    if y2.value < 0 :
        ypos10pos[x] = 1
    
    y3 = d_Mroz18(xx,-10,5)
    if y3.value < 0 :
        yneg10pos[x] = 1
    
    y4 = d_Mroz18(xx,-10,-5)
    if y4.value < 0 :
        yneg10neg[x] = 1
    
    y5 = d_Mroz18(xx,10,-5)
    if y5.value < 0 :
        ypos10neg[x] = 1
    
    
# plt.plot(xs, y0pos)
plt.plot(xs,ypos10neg)
plt.plot(xs,ypos10pos)
# plt.plot(xs,yneg10neg)
# plt.plot(xs,yneg10pos)


lims = [dht21.lons[0], dht21.lons[-1], dht21.vels[0], dht21.vels[-1]]
plt.imshow(np.nansum(np.abs(dht21.cube),axis=1),**kw_imshow,extent=lims)


plt.figure(figsize=(20,6))
lims = [cogal.lons[0],cogal.lons[-1],cogal.vels[0],cogal.vels[-1]]
plt.imshow(np.nansum(np.abs(cogal.cube),axis=1),**kw_imshow,extent=lims)
plt.plot(xs,-150+ypos10neg*300)
plt.plot(xs,-150+ypos10pos*300)
plt.show()


plt.plot(xs*u.deg.to(u.rad)/np.pi, y0pos)
# # plt.plot(xs,ypos10neg)
# plt.plot(xs,ypos10pos)
# # plt.plot(xs,yneg10neg)
# plt.plot(xs,yneg10pos)


cogal = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')


xs = np.arange(-2*np.pi,2*np.pi,0.01)


ys = np.sin(xs)


plt.plot(xs,ys)
plt.vlines(np.pi/2,-1,1)
plt.vlines(-np.pi/2,-1,1)

plt.hlines(0,-6.28,6.28)


dht21 = dc.unpack_fits_cube('DHT21.fits','DHT')


dht21.cube.shape


pers = dc.spectrum1D(dht21.cube[:,25:100,200:300])


dc.weighted_v(np.nan_to_num(dht21.cube[:,25:100,200:300])*u.K,dht21.vels*u.km/u.s,dht21.dv)


dc.spectrum1D(np.nan_to_num(dht21.cube[:,25:100,200:300]))


wco = dc.numba_nansum_ax0(dht21.cube[:,25:100,200:300]*dht21.dv) #include data units at somepoint
wco_tot = np.nansum(wco)

num = 0
for v,vv in enumerate(dht21.vels):
    this_slice = dht21.cube[:,25:100,200:300][v] 
    wcoi = np.nansum(this_slice*dht21.dv*vv) #include data units at somepoint
    num += wcoi

vc = num/wco_tot


vc


dht21.vels[44]


plt.plot(dht21.vels,pers)
plt.scatter(dht21.vels[44],pers[44])


lims = [dht21.lons[0],dht21.lons[-1], dht21.lats[0],dht21.lats[-1]]
# plt.imshow(dht21.cube[40,25:100,200:300],**kw_imshow)#,cmap=dc.cmap_grey0('viridis'))
plt.imshow(np.nansum(dht21.cube,axis=0),extent = lims,**kw_imshow)#,cmap=dc.cmap_grey0('viridis'))


cat


v_unc = np.full(18,2)


pfar = np.full(18,0.5)


cat['v_unc'] = v_unc


cat['pfar'] = pfar


cat


(-2.07)**2.


(2)**0.9


(8.34*u.pc*1e3*(231.211/-111.519)**2)**0.1746


(-3)**1.5


3*3**0.5


8.34*u.pc*1e3*(231.211/-111.519)**(2/(1-0.083))


nak = pd.read_csv('untitled.txt',skiprows=2,delim_whitespace=True,header=None)


ls, bs, vs = nak[1], nak[2], nak[3]


for i in range(20):
    d = d_Russeil17(ls[i],bs[i],vs[i])
    print(d)


nak


(1/(1-0.0803))*2


theta/thetasun = a1(R/Rsun)^a2 + a3


1/(1-0.0803)


1.0873*2


a1=1.00767
a2 = 0.0394
a3 = 0.00712


len(dco), len(dco.trunk)


dco = Dendrogram.compute(cogal.cube, 3*0.3, 5*0.3,50)


dco = Dendrogram.compute(cogal.cube, 3*0.3, 2*0.3,50)


co_mtdt = {}
co_mtdt['data_unit'] = u.K
co_mtdt['beam_major'], co_mtdt['beam_minor'] = 8.4*u.arcmin, 8.4*u.arcmin
co_mtdt['spatial_scale'] = DHT_CO.db #same as DHT_CO.dl
co_mtdt['velocity scale'] = DHT_CO.dv
co_mtdt['wavelength'] = c.c/(115e9*u.Hz)


# orion_mtdt = {}
# orion_mtdt['data_unit'] = u.Jy #works for u.Jy (but not u.K ?? )
# orion_mtdt['spatial_scale'] = (1/8)*u.deg
# # orion_mtdt['beam_major'], orion_mtdt['beam_minor'] = (1/4)*u.deg, (1/4)*u.deg


# %%time
# orcat = ppv_catalog(dend,orion_mtdt)

# #the catalog creates a table of the 'structures'
# #(ie individual leaves) and their properties


galcat


# orcat


# orhead = dc.unpack_fits_cube('DHT27.fits','DHT').data_header


galhead = cogal.data_header


# %%time
# orion_dclust = SpectralCloudstering(dend, orcat,orhead)


# orclusts = orion_dclust.clusters_asgn.data 


galclusts = gal_clust.clusters_asgn.data 


# orion_dclust.clusters_asgn.writeto('orion_dend_clusters.fits',overwrite=True)


gal_clust.clusters_asgn.writeto('gal_dend_clusters.fits',overwrite=True)


# dir(orion_dclust)


# for cnum in gal_clust.clusters:
#     mskd = (galclusts == cnum)*galclusts
#     plt.imshow(np.nansum(mskd,axis=0),**kw_imshow)
#     plt.title(str(cnum))
#     plt.show()


# orcat


trs = galcat[galcat['structure_type']=='T']


trs.sort('n_leaves')


trs.sort('_idx')


own_a = galcat[galcat['ancestor']==galcat['_idx']]


own_a.sort('_idx')


plt.hist(1*u.deg.to(u.rad)*1e3*own_a['radius'],bins=100)
plt.yscale('log')
# plt.show()


len(trs)


for c in own_a['_idx']:
    print(c)


plt.figure(figsize=(14,5))
plt.imshow(np.nanmean(galclusts,axis=1),cmap=dc.cmap_grey0('Paired'))#[150]*(galclusts[150]==880),cmap=dc.cmap_grey0('Paired'))
plt.colorbar()


plt.figure(figsize=(14,5))
plt.pcolormesh(galclusts[250])#*(galclusts[150]==1200))#,cmap=dc.cmap_grey0('Paired'))
plt.colorbar()


for c in trs['_idx']:
    plt.imshow(np.nansum((galclusts==c)*galclusts,axis=0),**kw_imshow)
    plt.show()


galclusts.shape


galcat


# orclusts.shape


# from numba import njit


# @njit
# def cluster_masked_list (cnums_list, data):
#     x,y,z = data.shape
#     output = np.zeros((len(cnums_list),x,y,z))
#     for n, cnum in enumerate(cnums_list):
#         msk = (data == cnum)
#         mskd = data * msk
#         output[n] = mskd

#     return output


# %%time
# gal_cl_list = cluster_masked_list(gal_clust.clusters, galclusts)


# or_cl_list.shape


# orion_dclust.clusters.index(18)


# np.nanmin(or_cl_list[7]),np.nanmax(or_cl_list[7])


# orb = or_cl_list[7]/or_cl_list[7]


# plt.pcolormesh((orb*orcube)[32])
# plt.colorbar()


# plt.pcolormesh(or_cl_list[7][30],cmap =dc.cmap_grey0('Paired'))
# plt.colorbar()


# plt.pcolormesh(orclusts[30],cmap =dc.cmap_grey0('Paired'))
# plt.colorbar()


# len(orion_dclust.clusters)


# orclob = dc.unpack_fits_cube('DHT27.fits','DHT')
# orcube = orclob.cube


# orcube[30].shape


# orclob.lons.shape, orclob.lats.shape


# cogal2 = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')


# fig = plt.figure(figsize=(20,15))
# # ax = fig.add_subplot(projection='3d')
# x,y = np.meshgrid(cogal2.lons, cogal2.lats)
# z = (np.nansum(cogal2.cube[:,:,:], axis=0))

# # ax.plot_wireframe(x,y,z,cmap = 'gray',alpha = 0.9,linewidth=0.3)
# # ax.view_init(azim=90,roll=0,elev=60)

# # ax2 = fig.add_subplot()
# # ax2.imshow(z,origin='lower',cmap=cmap0)
# # # fig.colorbar(ax2)
# # plt.show()


# plt.figure(figsize=(20,15))
# z2 = (np.nansum(cogal2.cube[:,:,:],axis=1))**0.5
# plt.imshow(z2,origin='lower',cmap=cmap0)
# plt.show()


# plt.figure(figsize=(20,15))
# z2 = (np.nansum(cogal2.cube*cogal2.dv.value,axis=0))
# plt.imshow(z2,origin='lower',cmap=cmap0)
# plt.show()


# fig = plt.figure(figsize=(20,15))
# # ax = fig.add_subplot(projection='3d')
# x,y = np.meshgrid(cogal.lons, cogal.lats)
# z = np.sqrt(np.nansum(cogal.cube, axis=0))

# # ax.plot_wireframe(x,y,z,cmap = 'gray',alpha = 0.9,linewidth=0.3)
# # ax.view_init(azim=90,roll=0,elev=60)

# ax2 = fig.add_subplot()
# ax2.imshow(z,origin='lower',cmap=cmap0)
# # fig.colorbar(ax2)
# plt.show()

# plt.figure(figsize=(20,15))
# z2 = np.nansum(cogal.cube,axis=0)
# plt.imshow(z2,origin='lower',cmap=cmap0)
# plt.show()


# np.nanmax(cogal.cube),np.nanmin(cogal.cube)


# np.nanmax(cogal.cube),np.nanmin(cogal.cube)


# xs = np.arange(0,10,0.1)


# plt.plot(xs,xs)
# plt.plot(xs,np.sqrt(xs))
# # plt.plot(xs,xs**2)


# dir(orion_dclust)


# trunks_flat = np.nanmedian(orion_dclust.trunks_asgn.data, axis = 0)


# trunks_flat.shape


# orion_dclust.trunks_asgn.data.shape


# plt.imshow(trunks_flat, origin = 'lower')
# plt.colorbar()


# plt.pcolormesh(orion_dclust.trunks_asgn.data[6])
# plt.colorbar()

# plt.show()


# HI_dend = Dendrogram.compute(HI_cube, \
#                                 min_value=dend_sigma,\
#                                 min_delta=2*dend_sigma,\
#                                 min_npix=3*dend_ppb)


# HI_mtdt = {} #metadata
# HI_mtdt['data_unit'] = u.Jy
# HI_catalog = ppv_catalog(HI_dend, HI_mtdt)


# HI_dclust = SpectralCloudstering(HI_dend, HI_catalog, HI_head)


# HI_dclust.clusters_asgn.data.shape


# plt.pcolormesh(HI_dclust.clusters_asgn.data[80])


# would be interesting/helpful to isolate the different velocity peaks & see how they map onto spatial features


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
