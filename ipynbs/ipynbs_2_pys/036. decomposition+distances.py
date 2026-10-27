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


import transform as tr
LAB = dc.unpack_fits_cube('labh.fits','LAB')


len(HI4PI.vels), len(HI4PI.lats), len(HI4PI.lons)


HI4PI.dv, HI4PI.db, HI4PI.dl


len(LAB.vels), len(LAB.lats), len(LAB.lons)


LAB.dv, LAB.db, LAB.dl


1e-3*HI4PI.vels[0],1e-3*LAB.vels[0]


1e-3*HI4PI.vels[-1],1e-3*LAB.vels[-1]


HI4PI.lats[0],LAB.lats[0]


HI4PI.lons[0],LAB.lons[0]


HI4PI.lats[-1],LAB.lats[-1]


dDHT = Dendrogram.compute(DHT_CO.cube,3*0.4,0.2,1e3)


for t in dDHT.trunk:
    m = t.get_mask()
    this = m*DHT_CO.cube
    this0th = np.nansum(this, axis=0)
    vs,bs,ls = t.indices()
    v0 = np.min(vs)
    v1 = np.max(vs)
    
    b0 = np.min(bs)
    b1 = np.max(bs)
    
    l0 = np.min(ls)
    l1 = np.max(ls)
    
    lims = [DHT_CO.lons[l0], DHT_CO.lons[l1], DHT_CO.lats[b0], DHT_CO.lats[b1]]
    
#     plt.imshow(this0th,**kw_imshow)#,extent = lims)
#     plt.show()


    fig, axs = plt.subplots(1,3,figsize=(15,3))
    axs[0].plot(DHT_CO.vels,np.nanmean(this, axis=(1,2)))
    axs[1].plot(DHT_CO.lons,np.nanmean(this, axis=(0,1)))
    axs[2].plot(DHT_CO.lats,np.nanmean(this, axis=(0,2)))
#     plt.title('vels are '+str(DHT_CO.vels[v0])+' and '+str(DHT_CO.vels[v1]))
#     plt.colorbar()
    plt.show()
    

# DHT_CO.lats[201],DHT_CO.lats[281]
# DHT_CO.lons[641],DHT_CO.lons[2241]
wo_gal_cent = np.zeros_like(DHT_CO.cube)
wo_gal_cent[:,:,:641] = DHT_CO.cube[:,:,:641]
wo_gal_cent[:,:,2241:] = DHT_CO.cube[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = DHT_CO.cube[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = DHT_CO.cube[:,281:,641:2241]


np.nanstd(wo_gal_cent)


len(d_wocent.trunk)


d_wocent.trunk[80].smallest_index


v,b,l = d_wocent.trunk[9].smallest_index


tr = d_wocent.trunk[9]
m = tr.get_mask()
this_cube = m*wo_gal_cent

this_clob = DHT_CO.clone()
this_clob.cube = this_cube


import kd as kd


this_clob.lons[l]


this_clob.lats[b]


this_clob.vels[v]


possible_Rs


dht24 = dc.unpack_fits_cube('DHT21.fits','DHT')


plt.plot(dht24.vels,dc.spectrum1D(dht24.cube[:,25:100,230:300]))


plt.imshow(np.nansum(dht24.cube[:,25:100,230:300],axis=0))


plt.imshow(dht24.calculate_0th_moment().value)


cat


R_poly_abys2(possible_Rs,160.,-19.,5*u.km/u.s)[60:65]


possible_Rs[62]


R_poly_abys(7980*u.pc,160.,-19.5,5*u.km/u.s)


R_poly(7500*u.pc,160.,-19.5,5*u.km/u.s)


R_poly(8956*u.pc,211.,-19.,5*u.km/u.s)


dist(7500*u.pc,160.,-19.)


pd.read_csv('summary-Copy1.prt',delim_whitespace=True,skiprows=2)


R_sun


dist(8956*u.pc,211.,-19.)


dist_yarn(7500*u.pc,R_sun,160.,-19.5)


A1,A2,A3


def dist_yarn(R, r_sun, l, b):
    term1 = r_sun*np.cos(b*u.deg.to(u.rad))*np.cos(l*u.deg.to(u.rad))
    
    term2_1 = term1**2
    term2_2 = (r_sun**2)
    term2_3 = R**2
    
    term2 = np.sqrt(term2_1 - term2_2 + term2_3)
    
    final = term1 + term2
    
    return final


orRs[72],orRs[73]


possible_Rs[62],possible_Rs[63]


plt.plot(perss[60:65],marker='+')


#velocity curve from https://doi.org/10.1051/0004-6361/202346474
#slope -- 0.4±0.6km/s/kpc, V_sun = 233±7km/s 


(1/60)*8 + (1/3600)*4


1/8


def R_poly_abys2(R, l, b, v_lsr, theta_sun=Theta_sun, r_sun=R_sun, a1=A1, a2=A2, a3=A3):
    sinl = np.sin(l*u.deg.to(u.rad))
    cosb = np.cos(b*u.deg.to(u.rad))
    
    coeff1 = a1*Theta_sun/(r_sun**a2)
    term1 = coeff1*(R**a2)
    
    coeff2 = (v_lsr/(sinl * cosb) + Theta_sun)/r_sun
    term2 = coeff2*R
    
    const = a3*Theta_sun
    
    final = term1 - term2 + const
    
    return final


def R_poly_abys(R, l, b, v_lsr, theta_sun=Theta_sun, r_sun=R_sun, a1=A1, a2=A2, a3=A3):
    sinl = np.sin(l*u.deg.to(u.rad))
    cosb = np.cos(b*u.deg.to(u.rad))
    
    coeff1 = a1*sinl*Theta_sun*r_sun
    term1 = coeff1*(R**a2)
    
    coeff2 = (v_lsr/cosb + Theta_sun*sinl)*r_sun**(a2)
    term2 = coeff2*R
    
    const = a3*sinl*Theta_sun*r_sun**(a2+1)
    
    final = term1 - term2 + const
    
    return final


def R_poly(R, l, b, v_lsr, theta_sun=Theta_sun, r_sun=R_sun, a1=A1, a2=A2, a3=A3):

    term1 = R/r_sun
    term2 = (R/r_sun)**a2
    
    coeff1 = (v_lsr/np.cos(b))+ (theta_sun*np.sin(l))
    coeff2 = a1*theta_sun*np.sin(l)
    coeff3 = a3*theta_sun*np.sin(l)
             
    final = coeff1*term1 - coeff2*term2 - coeff3
    
    print(final)
    
    return final


def dist(R, l, b, r_sun=R_sun):
    term1 = r_sun*np.cos(l*u.deg.to(u.rad))/np.cos(b*u.deg.to(u.rad))
    
    term2_1 = R**2
    term2_2 = (r_sun*np.sin(l*u.deg.to(u.rad)))**2
    term2 = np.sqrt(term2_1 - term2_2)/np.cos(b*u.deg.to(u.rad))
    
    print(term2_1, term2_2)
        
    if R<r_sun:
        near = term1 - term2
        far = term1 + term2
        
        return [near,far]
    else:
        d = term1 + term2
        return [d]


#===========================================
possible_Rs = np.arange(1.7, 17, 0.1)*1000*u.pc

poly_out = np.zeros(len(possible_Rs))
def calculate_kd(l,b,v_lsr):
    l = l*u.deg
    b = b*u.deg
    v_lsr = v_lsr * u.km/u.s
    
    for i, posR in enumerate(possible_Rs):
        p = R_poly(posR,l,b,v_lsr)

        poly_out[i] = p.value
        
#     print(poly_out)

    this = np.argmin(np.abs(0-poly_out))
    
    R_gal = possible_Rs[this]
#     print(poly_out)

    d = dist(R_gal,l,b)

    return d


this_clob.lons[l],this_clob.vels[v]


np.sin(359*u.deg.to(u.rad))


for t, tr in enumerate(d_wocent.trunk):
    m = tr.get_mask()
    this_cube = m*wo_gal_cent

    this_clob = DHT_CO.clone()
    this_clob.cube = this_cube


    dSig, dA = cld.surface_density(this_clob, 'CO',show_Sigma=False)
    v,b,l = tr.smallest_index
    d = calculate_kd(this_clob.lons[l],this_clob.lats[b],this_clob.vels[v])
    
#     if d[0] < 0*u.pc:
#         print(t)


plt.plot(R_poly_abys2(possible_Rs,160.,-19.5,5*u.km/u.s)[60:65],marker='+')


#values from Russeil et al 2017
Theta_sun = 220*u.km/u.s
R_sun = 8.5*1000*u.pc
A1 = 1.00767
A2 = 0.0394
A3 = 0.00712


cat.to


cat.to_csv('cat',sep=' ')


pd.read_csv('wo_cent_MsAsds_from_dend.csv')


from_wocent = pd.read_csv('wo_cent_MsAsds_from_dend.csv')


from_wocent['M [kg]']


np.array((8,9))*u.kg


52850*u.lightyear.to(u.pc)


possible_Rs_finer = np.arange(1,17e3,1)*u.pc


possible_Rs_finer


lst = R_poly_abys2(possible_Rs_finer,190,-19.5,5*u.km/u.s)


np.argmin(np.abs(lst-0))


lst


possible_Rs_finer[7961]


(R_sun*np.sin(300*u.deg.to(u.rad)))**2


R_sun*R_sun*(np.sin(310*u.deg.to(u.rad)))**2


6000*u.pc/R_sun


np.sin(310)


ds = np.arange(0,360)*u.deg


plt.plot(ds,(np.sin(ds.to(u.rad)))**2)


DHT_CO.vels


plt.plot(ds,(np.sin(ds.to(u.rad))))


0.01*R_sun


A1*0.01**A2


HI4PI.vels/1e3


R_sun*0.87


7000*u.pc/R_sun


np.arange(1,15)*u.pc*1e3


dist_yarn


R_sun


plt.plot(np.arange(1,15)*u.pc*1e3,dist_yarn(np.arange(1,15)*u.pc*1e3,R_sun,70,-19.5))
plt.hlines(0,4000,14000)
plt.vlines(8500,-6000,6000)


possible_Rs[57:59]


np.sin(160*u.deg.to(u.rad))


np.cos(285*u.deg.to(u.rad))


90-180 cos is < 0


d_wocent = Dendrogram.compute(wo_gal_cent,5*0.2,5*0.2,100)


bigm = d_wocent.trunk[1].get_mask()
smolm = d_wocent.trunk[0].get_mask()

thresh = 1e3
for t,tr in enumerate(d_wocent.trunk):
    npx = tr.get_npix()
    if npx>=thresh:
        m = tr.get_mask()
        bigm +=m
        
    else:
        m = tr.get_mask()
        smolm +=m
        

fig,ax = plt.subplots(2,1,figsize=(10,5))


lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0], DHT_CO.lats[-1]]
ax[0].imshow(np.nansum((smolm)*np.abs(wo_gal_cent),axis=0),**kw_imshow,extent=lims)
#solar neighbors
ax[0].add_patch(Rectangle((115,20),15,10,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((-70,-20),20,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((150,-25),30,20,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-160,-25),20,20,color='cyan',fill=False))
ax[0].add_patch(Rectangle((40,5),10,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-10,10),10,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-5,-25),10,10,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((-175,-20),15,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((90,-20),10,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-165,-5),15,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((25,-25),15,15,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((145,15),5,5,color='cyan',fill=False,linestyle='--'))

ax[1].imshow(np.nansum((smolm)*np.abs(wo_gal_cent),axis=1),**kw_imshow,extent=lims)
plt.show()


len(d_wocent.trunk)


fig,ax = plt.subplots(2,1,figsize=(10,5))


lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0], DHT_CO.lats[-1]]
ax[0].imshow(np.nansum((bigm)*np.abs(wo_gal_cent),axis=0),**kw_imshow,extent=lims)
#solar neighbors
ax[0].add_patch(Rectangle((115,20),15,10,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((-70,-20),20,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((150,-25),30,20,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-160,-25),20,20,color='cyan',fill=False))
ax[0].add_patch(Rectangle((40,5),10,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-10,10),10,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-5,-25),10,10,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((-175,-20),15,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((90,-20),10,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-165,-5),15,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((25,-25),15,15,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((145,15),5,5,color='cyan',fill=False,linestyle='--'))

ax[1].imshow(np.nansum((bigm)*np.abs(wo_gal_cent),axis=1),**kw_imshow,extent=lims)
plt.show()


plt.imshow(np.nansum(np.abs(wo_gal_cent),axis=0),**kw_imshow)


fig,ax = plt.subplots(2,1,figsize=(10,5))

ax[0].imshow(np.nansum((bigm)*np.abs(wo_gal_cent),axis=2),**kw_imshow,extent=lims)
ax[1].imshow(np.nansum((smolm)*np.abs(wo_gal_cent),axis=2),**kw_imshow,extent=lims)
ax[0].set_xlim(-130,-160)


d_wocent.trunk[1].get_npix()


smolm = d_wocent.trunk[0].get_mask()
for t,tr in enumerate(d_wocent.trunk):
    npx = tr.get_npix()
    count_px[t] = npx
    if npx<800:
#         print(npx)
        m = tr.get_mask()
        smolm +=m
        

fig,ax = plt.subplots(figsize=(10,10))

#solar neighbors
ax.add_patch(Rectangle((115,20),15,10,color='cyan',fill=False,linestyle='--'))
ax.add_patch(Rectangle((-70,-20),20,10,color='cyan',fill=False))
ax.add_patch(Rectangle((150,-25),30,20,color='cyan',fill=False))
ax.add_patch(Rectangle((-160,-25),20,20,color='cyan',fill=False))
ax.add_patch(Rectangle((40,5),10,10,color='cyan',fill=False))
ax.add_patch(Rectangle((-10,10),10,10,color='cyan',fill=False))
ax.add_patch(Rectangle((-5,-25),10,10,color='cyan',fill=False,linestyle='--'))
ax.add_patch(Rectangle((-175,-20),15,15,color='cyan',fill=False))
ax.add_patch(Rectangle((90,-20),10,15,color='cyan',fill=False))
ax.add_patch(Rectangle((-165,-5),15,15,color='cyan',fill=False))
ax.add_patch(Rectangle((25,-25),15,15,color='cyan',fill=False,linestyle='--'))
ax.add_patch(Rectangle((145,15),5,5,color='cyan',fill=False,linestyle='--'))


lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0], DHT_CO.lats[-1]]
ax.imshow(np.nansum(smolm*wo_gal_cent,axis=0),**kw_imshow,extent=lims)
plt.show()


cat


plt.plot(count_px, marker='+')
# plt.yscale('log')
plt.ylim(0,800)


plt.hist(count_px,bins=10)
# plt.xscale('log')
# plt.yscale('log')
plt.show()


DHT_CO.lats[200],DHT_CO.lats[282]


# 0.37526998 one sigma od DHT_CO.cube


DHT_CO.cube.shape


55*481*2881


76216855/1e2


h = plt.hist((DHT_CO.cube.flat),bins=100)
plt.xlim(0,7)
plt.vlines(0.37, 0, 7)
plt.yscale('log')
plt.show()


plt.plot(DHT_CO.vels,np.nanmean(DHT_CO.cube, axis=(1,2)))


len(DHT_CO.vels)


mask = DHT_CO


np.zeros((481,2881))


DHT_CO.cube.shape


a = np.array((3, 5, 6, 1, 7, 3, 4,  0, 9))


a * (a>6) + a * (a==3)


a * (a==3)


a[:3]


wo_gal_cent = np.zeros_like(DHT_CO.cube)


# DHT_CO.lats[201],DHT_CO.lats[281]
# DHT_CO.lons[641],DHT_CO.lons[2241]
wo_gal_cent = np.zeros_like(DHT_CO.cube)
wo_gal_cent[:,:,:641] = DHT_CO.cube[:,:,:641]
wo_gal_cent[:,:,2241:] = DHT_CO.cube[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = DHT_CO.cube[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = DHT_CO.cube[:,281:,641:2241]


DHT_CO.lons[641],DHT_CO.lons[2241]


plt.imshow(np.nansum(np.abs(wo_gal_cent),axis=0),**kw_imshow)
