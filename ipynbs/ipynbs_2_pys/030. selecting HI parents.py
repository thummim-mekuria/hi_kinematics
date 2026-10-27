#!/usr/bin/env python
# coding: utf-8


# %%time

import datacube as dc
import clouds as cld

import numpy as np

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
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf

DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
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
        ns = ['03','06','09','12','16','19','21','26','27','37']
        #for now, figure out DHT01 (bad data) 
        #& DHT18 (two ceph ds, beeg) 
        #& DHT08 (beeg) 
        #& DHT24, DHT34 (lam or & chamaeleon are disjointed down to the really low sig levels)
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output


#---------------------
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


for s in surveys:
    clob = surveys[s][0]
    sig = np.nanstd(clob.cube)

    clouds = surveys[s][1]
    
    if (s not in ['DHT03','DHT06', 'DHT34','DHT37']) and (s not in ['DHT21']) and (s != 'DHT08'):
        d = Dendrogram.compute(clob.cube, min_value=3*sig,min_delta=sig, min_npix=400)
        
    elif s in ['DHT03','DHT06', 'DHT34','DHT37']:
        d = Dendrogram.compute(clob.cube, min_value=sig,min_delta=sig, min_npix=500)
        
    elif s == 'DHT21':
        d = Dendrogram.compute(clob.cube, min_value=6*sig,min_delta=sig, min_npix=1000)

        
    for cl in clouds:
        l_c = cld.ind_ang('get index', clob.lons, [cat['l_c__deg'][cl]])
        b_c = cld.ind_ang('get index', clob.lats, [cat['b_c__deg'][cl]])


        for tr in d.trunk:
            inds = np.array(tr.indices())
            coords = np.transpose(inds[1:]).tolist()            
            if [b_c, l_c] in coords:
#                 print(s, cl,'success')
                m = tr.get_mask()

                new_cube = clob.cube * m
                new_clob = clob.clone()
                new_clob.cube = new_cube

                cat['clob'][cl] = new_clob


locmom = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


locmom0th = locmom.calculate_0th_moment()


(10*c.M_sun*u.pc**-2/(1.36*c.m_p*1.82e18*u.s*(u.cm**2*u.K*u.km)**-1)).to(u.K*u.km/u.s)


((1.5*u.K*u.km/u.s)*(2*1.36*c.m_p)*(1.97e20*u.s/(u.K*u.km*u.cm**2))).to(c.M_sun/u.pc**2)


cat


loHI = 0
liHI = 550

boHI = 750
biHI = 1150

imHI_tau = np.nansum(HI4PI.cube[:,boHI:biHI,loHI:liHI],axis=0)
xHI_tau = HI4PI.lons[loHI:liHI]
yHI_tau = HI4PI.lats[boHI:biHI]


# HI4PI.lons[0],HI4PI.lons[50]
HI4PI.lats[750],HI4PI.lats[800]


plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,750:800,0:50]),linewidth=0.7,label='l=180-175, b=-27.5- -23.5')
plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,750:800,500:550]),linewidth=0.7,label='l=138-134, b=-27.5- -23.5')

# plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,100:150,850:900]),linewidth=0.7,label='100-150,850-900')

# plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,10:20,760:770]),linewidth=0.7)
# plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,20:30,770:780]),linewidth=0.7)
# plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,30:40,780:790]),linewidth=0.7)
# plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,40:50,790:800]),linewidth=0.7)
# plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,400:450,900:950]),linewidth=0.7,label='400-450,900-950')
plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,1100:1150,0:50]),linewidth=0.7,label='l=138-134, b=-27.5- -23.5')

plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,1100:1150,500:550]),linewidth=0.7,label='l=138-134, b=1.5-5.6')


plt.xlim(-50000,50000)
plt.legend()
plt.show()


plt.imshow(imHI_tau,**kw_imshow,extent=[xHI_tau[0],xHI_tau[-1],yHI_tau[0],yHI_tau[-1]])


# plt.contour(xCO,yCO,imCO,levels = [1,5],colors=['black','grey'],linewidths=1)
plt.contour(xHI_tau,yHI_tau,imHI_tau,levels = [500,1000,1500],\
            colors=['white','red','yellow'],linewidths=0.5)
# #levels = [500,750,1210]
#,1210,1260,1320,1430,1700,1810,1850]
#,'red','orange','yellow','white','red','orange','yellow','white']
plt.xlim(180,150)
plt.ylim(-30,0)
plt.show()


imHI2 = np.nansum(HI4PI.cube[:,b0HI:b1HI,l0HI:l1HI],axis=1)
vHI = HI4PI.vels


plt.imshow(imHI2,**kw_imshow,extent=[xHI[0],xHI[-1],vHI[0]/1e3,vHI[-1]/1e3])


plt.contour(xHI,vHI/1e3,imHI2,colors=['red','white','yellow'],levels=[1e4,2e4,3e4])
plt.contour(locmom.lons[2500:2750],locmom.vels,np.nansum(ors,axis=1),colors=['white','red','yellow'],\
            levels=[10,20,50],linewidths=[0.5])
plt.xlim(-130,-180)
plt.ylim(-20,20)
# plt.colorbar()
plt.show()


l0HI = 3800
l1HI = 4100

b0HI = 750
b1HI = 1050

imHI = np.nansum(HI4PI.cube[:,b0HI:b1HI,l0HI:l1HI],axis=0)
xHI = HI4PI.lons[l0HI:l1HI]
yHI = HI4PI.lats[b0HI:b1HI]


xCO = locmom.lons
yCO = locmom.lats
imCO = locmom0th.value


dHI = Dendrogram.compute(imHI,500,100,100)


# fig, ax = plt.subplots(1,1,figsize=(10,10))
dHI.viewer().show()
# k.show()
# fig.show()


len(dHI.trunk)


m


for i in dHI.trunk:
    m = i.get_mask()
    plt.imshow(imHI*m,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])


#     plt.contour(xCO,yCO,imCO,levels = [4],colors=['white'],linewidths=1.5)
#     plt.contour(xHI,yHI,imHI,levels = [500,1210,1260,1320,1430,1700,1810,1850],\
#         colors=['white','red','orange','yellow','white','red','orange','yellow','white'],linewidths=0.5)
    #levels = [500,750,1210]
#     plt.colorbar()
    plt.xlim(-135,-160)
    plt.ylim(-25,-5)
    plt.show()
    

30*6.4


np.nanstd(HI4PI.cube[:,b0HI:b1HI,l0HI:l1HI])


dHI.trunk


mHI = dHI.trunk[0].get_mask()


plt.imshow(imHI*mHI,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])


plt.contour(xCO,yCO,imCO,levels = [4],colors=['white'],linewidths=1.5)
plt.contour(xHI,yHI,imHI,levels = [500,1210,1260,1320,1430,1700,1810,1850],\
            colors=['white','red','orange','yellow','white','red','orange','yellow','white'],linewidths=0.5)
#levels = [500,750,1210]
plt.xlim(-130,-180)
plt.ylim(-25,10)
plt.show()


cat


d2 = Dendrogram.compute(np.nansum(ors,axis=0),3*s2,3*s2,100)


3*s2


s2 = np.nanstd(ors)


d2.trunk


for i in d2.trunk:
    plt.imshow(i.get_mask() * d2.data,**kw_imshow)
    plt.show()


fig, axs = plt.subplots(1,2,figsize=(10,10))
y = locmom.lats
x = locmom.lons
vmn, vmx = zlims(locmom0th)
axs[0].imshow(locmom0th.value,**kw_imshow,extent=[x[0],x[-1],y[0],y[-1]],vmin=vmn,vmax=vmx)
axs[0].contour(x,y,locmom0th.value,levels = [4],colors=['red'],linewidths=0.3)
axs[0].set_xlim(-130,-160)
axs[0].set_ylim(-25,-1)

m = locmom0th.value > 1.5
axs[1].imshow(locmom0th.value*m,**kw_imshow,extent=[x[0],x[-1],y[0],y[-1]],vmin=vmn,vmax=vmx)
axs[1].contour(x,y,locmom0th.value,levels = [5],colors=['red'],linewidths=0.3)
axs[1].set_xlim(180,150)
# axs[1].set_ylim(-25,-10)

plt.show()


ors = locmom.cube[:,0:200,2500:2750]


trs = locmom.cube[:,0:250,:250]


plt.imshow(np.nansum(trs,axis=1),**kw_imshow)


cat


plt.imshow(np.nansum(ors,axis = 1),**kw_imshow)
plt.colorbar()


fig, axs = plt.subplots(1,2,figsize=(10,10))
y = locmom.lats
x = locmom.lons
vmn, vmx = zlims(locmom0th)
axs[0].imshow(locmom0th.value,**kw_imshow,vmin=vmn,vmax=vmx)#,extent=[x[0],x[-1],y[0],y[-1]]
# axs[0].contour(x,y,locmom0th.value,levels = [4],colors=['red'],linewidths=0.3)
# axs[0].set_xlim(-130,-160)
# axs[0].set_ylim(-25,-1)
axs[0].set_xlim(2500,2750)
axs[0].set_ylim(0,200)

m = locmom0th.value > 1.5
axs[1].imshow(locmom0th.value*m,**kw_imshow,vmin=vmn,vmax=vmx)#,extent=[x[0],x[-1],y[0],y[-1]]
# axs[1].contour(x,y,locmom0th.value,levels = [5],colors=['red'],linewidths=0.3)
axs[1].set_xlim(0,250)
axs[1].set_ylim(0,250)

plt.show()


plt.plot(dc.spectrum1D(locmom.cube))
plt.plot(dc.spectrum1D(locraw.cube))
plt.show()


cat


o = np.array((False, False, True))


t = np.array((True, False, False))


from copy import deepcopy
sss = deepcopy(surveys)


DHT21 = surveys['DHT21']


np.nanstd(DHT21[0].cube)


z.shape


y.shape


x.shape


DHT21[0].calculate_rms_noise()


z = DHT21[0].calculate_0th_moment() 
x,y = np.meshgrid(DHT21[0].lons, DHT21[0].lats)

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
ax.plot_surface(x,y,z,cmap = 'viridis',alpha = 0.7)
plt.show()


z = DHT21[0].calculate_0th_moment() * (DHT21[0].calculate_0th_moment() > 10*DHT21[0].calculate_rms_noise())
x,y = np.meshgrid(DHT21[0].lons, DHT21[0].lats)

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
ax.plot_surface(x,y,z,cmap = 'viridis',alpha = 0.7)
plt.show()


locraw = dc.unpack_fits_cube('COGAL_local_raw.fits','DHT')


locmom = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


im = locraw.calculate_0th_moment()
vmin,vmax = zlims(im)
plt.imshow(locraw.calculate_0th_moment(),vmin=vmin, vmax=vmax,**kw_imshow)
plt.show()


im = locmom.calculate_0th_moment()
vmin,vmax = zlims(im)
plt.imshow(locmom.calculate_0th_moment(),vmin=vmin, vmax=vmax,**kw_imshow)
plt.show()


cat


for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] - 6
    b1 = cat['b_c__deg'][n] + 6
    
    l0 = cat['l_c__deg'][n] - 6
    l1 = cat['l_c__deg'][n] + 6
    
    s = locraw.create_subcube([b0,b1,l0,l1])
    plt.imshow(s.calculate_0th_moment(),**kw_imshow)
    plt.title(cat['name'][n])
    plt.show()


for s in surveys:
    clobsss = surveys[s][0]
    sigsss = np.nanstd(clob.cube)

#     clouds = surveys[s][1]
    
#     if (s not in ['DHT03','DHT06', 'DHT34','DHT37']) and (s not in ['DHT21']) and (s != 'DHT08'):
#         d = Dendrogram.compute(clob.cube, min_value=3*sig,min_delta=sig, min_npix=400)
        
#     elif s in ['DHT03','DHT06', 'DHT34','DHT37']:
#         d = Dendrogram.compute(clob.cube, min_value=sig,min_delta=sig, min_npix=500)
        
#     elif s == 'DHT21':
#         d = Dendrogram.compute(clob.cube, min_value=6*sig,min_delta=sig, min_npix=1000)

        
    dsss = Dendrogram.compute(clob.cube, sigsss, sigsss, )
    for cl in clouds:
        l_c = cld.ind_ang('get index', clob.lons, [cat['l_c__deg'][cl]])
        b_c = cld.ind_ang('get index', clob.lats, [cat['b_c__deg'][cl]])


        for tr in d.trunk:
            inds = np.array(tr.indices())
            coords = np.transpose(inds[1:]).tolist()            
            if [b_c, l_c] in coords:
#                 print(s, cl,'success')
                m = tr.get_mask()

                new_cube = clob.cube * m
                new_clob = clob.clone()
                new_clob.cube = new_cube

                cat['clob'][cl] = new_clob
