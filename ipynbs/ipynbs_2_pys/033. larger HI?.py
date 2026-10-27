#!/usr/bin/env python
# coding: utf-8


# %%time

import datacube as dc
import clouds as cld

import numpy as np
# %matplotlib

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


for s in surveys.values():
    if len(s[1])==1: #only one cloud in this survey
        clob = s[0]
        sig = clob.calculate_rms_noise()
        d = Dendrogram.compute(clob.cube, sig, sig, 100)
        bigm = d.trunk[0].get_mask()
        for t in d.trunk:
            bigm += t.get_mask()
            
        new_cube = clob.cube*bigm
        
        new_clob = clob.clone()
        new_clob.cube = new_cube
        
        dSig, dA = cld.surface_density(new_clob, 'CO',show_Sigma=False)
        M, A = cld.mass_and_area(s[1][0],cat,dSig,dA)
        print(cat['name'][s[1][0]],1e-3*M.to(c.M_sun).value, np.sqrt(A/np.pi).value)
#         print(type(M))
#         l0,l1 = clob.lons[0],clob.lons[-1]
#         b0,b1 = clob.lats[0],clob.lats[-1]
#         v0,v1 = clob.vels[0],clob.vels[-1]
        
#         print(l0,l1,b0,b1)


surveys['DHT37'][0].calculate_rms_noise()


doph = Dendrogram.compute(surveys['DHT37'][0].cube,0.59,0.59,100)


bigm_plrs = dplrs.trunk[0].get_mask()
for t in dplrs.trunk:
    bigm_plrs += t.get_mask()


dplrs = Dendrogram.compute(surveys['DHT16'][0].cube,0.25,0.25,200)


bigm_plrs = dplrs.trunk[0].get_mask()
for t in dplrs.trunk:
    bigm_plrs += t.get_mask()


plt.imshow(np.nansum(bigm_plrs*surveys['DHT16'][0].cube,axis=1),**kw_imshow)


plt.imshow(np.nansum(bigm_plrs*surveys['DHT16'][0].cube,axis=0),**kw_imshow)
plt.colorbar()


surveys['DHT16'][0].calculate_rms_noise()


cat['d__pc'][15]


from scipy.signal import find_peaks,peak_widths
for s in surveys.values():
#     print(len(s[1]))
    if len(s[1])==1: #only one cloud in this survey
        clob = s[0]
        sig = clob.calculate_rms_noise()
        d1,d2,d3 = clob.cube.shape
        
        p = find_peaks(dc.spectrum1D(clob.cube), height=np.nanstd(dc.spectrum1D(clob.cube)))
#         if len(p[0])>1:
        plt.plot(clob.vels,dc.spectrum1D(clob.cube))
        plt.title(cat['name'][s[1][0]])
        plt.show()
        print(peak_widths(dc.spectrum1D(clob.cube),p[0]))
        print(np.nanstd(dc.spectrum1D(clob.cube)))

#             for pk in p[0]:
#                 print(peak_widths(dc.spectrum1D(clob.cube),[pk]))
# #         print(d1*d2*d3,'pix')
#         d = Dendrogram.compute(clob.cube,sig,sig,100) #min_npix is 1% of the cube
#         bigm = d.trunk[0].get_mask()
#         for t in d.trunk:
#             bigm += t.get_mask()
            
#         plt.imshow(np.nansum(clob.cube*bigm,axis=1),**kw_imshow)
#         plt.title(cat['name'][s[1][0]])
#         plt.show()


#when i come back, mess around with different sigma values?? other parameters?? 
#to decompose HI (maybe look at pycupid/clumpfind stuff too?
#read carina paper again if they have a reasoning for choosing the thresholds they did
#figure out accumulation radius stuffies (do they overlap etc)
#and also HI spectra?? maybe that will say which cloud things are aligned with


cocube = DHT_CO.cube[20:50,50:220,2500:2750]


cos = np.nanstd(cocube)


hicube = HI4PI.cube[400:600,600:1090,3600:4300]


hos = np.nanstd(hicube)


hod = Dendrogram.compute(hicube,6*hos,hos,1000)


fig,ax = plt.subplots(1,1,figsize=(8,8))
hod.plotter().plot_tree(ax)


hod.trunk


hod.trunk[2].descendants


fig,ax = plt.subplots(1,1,figsize=(10,10))
hi0 = np.nansum(hicube*HI4PI.dv.value,axis=0)
limshi = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
vnhi,vxhi = zlims(hi0)
plt.imshow(hi0,**kw_imshow,vmin=vnhi,vmax=vxhi,extent=limshi)
plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],co0,levels=[4],colors=['white'])
ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2,label='crossbones'))
ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2,label='orion a'))
ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2,label='orion b'))
ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2,label='mon r2'))
plt.legend()
plt.show()


for t in hod.trunk[2].descendants:
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    limshi = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
#     if t.is_leaf:
    ax.imshow(np.nansum(t.get_mask()*hicube,axis=0),**kw_imshow,extent=limshi)
    ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
    ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
    ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
    ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))

    plt.show()


# calculated accumulation radii & got:
# crossbones = 4 degrees
# mon_r2 = 7 degrees
# orion_a = 11 degrees
# orion_b = 11 degrees


from matplotlib.patches import Circle
fig,ax = plt.subplots(1,1,figsize=(10,10))
limsco = [DHT_CO.lons[2500], DHT_CO.lons[2750],DHT_CO.lats[50],DHT_CO.lats[220]]
co0 = np.nansum(cocube*DHT_CO.dv.value,axis=0)
vnco,vxco=zlims(co0)
ax.imshow(co0,**kw_imshow,extent=limsco,vmin=vnco,vmax=vxco)
# plt.colorbar()
ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
ax.add_patch(Circle((-149,-19),11,edgecolor='white',fill=False,linewidth=2))
ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))

plt.show()


206-360


cod = Dendrogram.compute(cocube,cos,cos,100)


fig,ax = plt.subplots(1,1,figsize=(8,8))
cod.plotter().plot_tree(ax)


hod = Dendrogram.compute(hicube,5*hos,hos,1000)


fig,ax = plt.subplots(1,1,figsize=(8,8))
hod.plotter().plot_tree(ax)


cat


oraclob = DHT_CO.create_subcube([-21,-17,-143,-155])


dSigora, dAora = cld.surface_density(oraclob,'CO')


Mora, Aora = cld.mass_and_area(11,cat,dSigora,dAora)


Mora, Aora


Morb, Aorb


cld.calculate_accumulation_radius(Mora, cat['d__pc'][11]*u.pc)


plt.imshow(oraclob.calculate_0th_moment().value,**kw_imshow)


cat['l_c__deg'][11]-360-6


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[7],colors='red')
plt.scatter(cat['l_c__deg'][11]-360,cat['b_c__deg'][11],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][11]-360+11,cat['l_c__deg'][11]-360-11)
plt.ylim(cat['b_c__deg'][11]-11,cat['b_c__deg'][11]+11)


orbclob = DHT_CO.create_subcube([-18,-10,-150,-158])


dSigorb, dAorb = cld.surface_density(orbclob,'CO')


Morb, Aorb = cld.mass_and_area(12,cat,dSigorb,dAorb)


Morb


cld.calculate_accumulation_radius(Morb,cat['d__pc'][12]*u.pc)


cat


mr2clob = DHT_CO.create_subcube([-15, -10,-143,-149])


dSigmr2, dAmr2 = cld.surface_density(mr2clob,'CO')


Mmr2, Amr2 = cld.mass_and_area(9, cat, dSigmr2,dAmr2)


cld.calculate_accumulation_radius(Mmr2,cat['d__pc'][9]*u.pc)


cat['b_c__deg'][9]-3


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[4.5],colors='red')
plt.scatter(cat['l_c__deg'][9]-360,cat['b_c__deg'][9],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][9]-360+7,cat['l_c__deg'][9]-360-7)
plt.ylim(cat['b_c__deg'][9]-7,cat['b_c__deg'][9]+7)


219-360


crsclob = DHT_CO.create_subcube([-13,-6,-138,-142])


dSigcrs, dAcrs = cld.surface_density(crsclob,'CO')


Mcrs, Acrs = cld.mass_and_area(5, cat, dSigcrs, dAcrs)


cld.calculate_accumulation_radius(Mcrs,cat['d__pc'][5]*u.pc)


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[5],colors='red')
plt.xlim(-135,-160)


hicube.shape


cocube.shape


DHT_CO.lats[50],DHT_CO.lats[220]


for i in range(0,200,5):
    plt.pcolormesh(hicube[i],cmap=cmap0)
#     plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
#                 np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[7],colors='red')
    plt.show()


#crossbones - 886pc
#mpn_r2 - 778pc
#orion_a - 432pc
#orion_b - 423pc


hicube_T = hicube.transpose().swapaxes(0,1)


from scipy.signal import find_peaks

find_peaks(dc.spectrum1D(hicube[:,200:210,300:310]),height=25)


plt.plot(dc.spectrum1D(hicube[:,200:210,300:310]))


specs = np.zeros((490,700))
for b,latslc in enumerate(hicube_T):
    for l,lonslc in enumerate(latslc):
        H = 10
        p = find_peaks(lonslc,height=H)
        nps = len(p[0])
        specs[b][l] = nps


HI4PI.lats[600],HI4PI.lats[1050]


1050-600


plt.figure(figsize=(10,10))
lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
plt.imshow(np.nansum(hicube_T,axis=2),origin='lower',cmap=cmap0,extent=lims)
plt.colorbar()
# plt.ylim(0,450)
plt.show()


hicube = HI4PI.cube[400:600,600:1090,3600:4300]


HI4PI.lons[4300]


cut2 = hicube[:,370:430,180:300].transpose().swapaxes(0,1)


cut2.shape


plt.plot(HI4PI.vels[400:600]/1e3,dc.spectrum1D(hicube[:,370:400,180:190]))
plt.plot(oraclob.vels,10*dc.spectrum1D(oraclob.cube),label='ora')
plt.plot(orbclob.vels,10*dc.spectrum1D(orbclob.cube),label='orb')
plt.plot(mr2clob.vels,10*dc.spectrum1D(mr2clob.cube),label='mr2')
plt.plot(crsclob.vels,10*dc.spectrum1D(crsclob.cube),label='crs')
plt.xlim(-20,20)


for k in range(0,60,5):
    for m in range(0,120,10):
        plt.plot(HI4PI.vels[400:600]/1e3,cut2[k][m],label='this px')
        plt.plot(oraclob.vels,10*dc.spectrum1D(oraclob.cube),label='ora')
        plt.plot(orbclob.vels,10*dc.spectrum1D(orbclob.cube),label='orb')
        plt.plot(mr2clob.vels,10*dc.spectrum1D(mr2clob.cube),label='mr2')
        plt.plot(crsclob.vels,10*dc.spectrum1D(crsclob.cube),label='crs')
        plt.legend()
        plt.xlim(-10,25)
        plt.show()


fig,ax = plt.subplots(figsize=(10,10))
# lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
plt.imshow(specs[370:430,180:300],origin='lower',cmap='Paired')#,extent=lims)
plt.colorbar()
# ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
# ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
# ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
# ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))
plt.show()


fig,ax = plt.subplots(figsize=(10,10))
lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
ax.imshow(specs,origin='lower',cmap='Paired',extent=lims)
# plt.colorbar()
ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))
# plt.xlim(-135,-145)
# plt.ylim(-10,-6)
plt.show()


cutout = hicube[:,100:301,300:502]


cutout_T = cutout.transpose().swapaxes(0,1)


plt.plot(dc.spectrum1D(oraclob.cube))


# when coming back -> look more at promise of HISA? esp where multiple peaks are reported
# also, look at taurus/california/perseus dht survey to see if similar pattern (CO peaks at the valleys of HI)
# is replicated
# look at HISA papers some more, look into gaussian decomposition
# also -- look at physical sizes of accumulation radii and if it can be justified to smoosh orions / mon / xbns tgt?
# ^ or if some are actually separate, do axis=1 np.nansum on the HI dendrogram to see how the leafs spread out in v
# ^ or also, gaussian decomposition??


for k in range(0,200,20):
    for m in range(0,200,20):
        plt.plot(HI4PI.vels[400:600]/1e3,cutout_T[k][m],label='this px')
        plt.plot(oraclob.vels,10*dc.spectrum1D(oraclob.cube),label='ora')
        plt.plot(orbclob.vels,10*dc.spectrum1D(orbclob.cube),label='orb')
        plt.plot(mr2clob.vels,10*dc.spectrum1D(mr2clob.cube),label='mr2')
        plt.plot(crsclob.vels,10*dc.spectrum1D(crsclob.cube),label='crs')
        plt.legend()
        plt.xlim(-10,25)
        plt.show()


cutout_T.shape


hicube.shape


fig,ax = plt.subplots(figsize=(10,10))
# lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
plt.imshow(specs[100:301,300:502],origin='lower',cmap='Paired')#,extent=lims)
plt.colorbar()
# ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
# ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
# ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
# ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))
plt.show()


plt.imshow(np.nansum(hicube,axis=0),origin='lower')


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[2.5],colors='red')
plt.scatter(cat['l_c__deg'][5]-360,cat['b_c__deg'][5],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][5]-360+4,cat['l_c__deg'][5]-360-4)
plt.ylim(cat['b_c__deg'][5]-4,cat['b_c__deg'][5]+4)


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[7],colors='red')
plt.scatter(cat['l_c__deg'][12]-360,cat['b_c__deg'][12],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][12]-360+11,cat['l_c__deg'][12]-360-11)
plt.ylim(cat['b_c__deg'][12]-11,cat['b_c__deg'][12]+11)
