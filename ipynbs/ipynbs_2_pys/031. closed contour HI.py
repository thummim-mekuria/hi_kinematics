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


locmom = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')
locmom0th = locmom.calculate_0th_moment()


#calculating the W_HI equivalent to Sigma_HI = 10Msun/pc**2
(10*c.M_sun*u.pc**-2/(1.36*c.m_p*1.82e18*u.s*(u.cm**2*u.K*u.km)**-1)).to(u.K*u.km/u.s)


#calculating the W_CO equivalent to Sigma_CO = 10Msun/pc**2
(10*c.M_sun*u.pc**-2/(2*1.36*c.m_p*1.97e20*u.s/(u.K*u.km*u.cm**2))).to(u.K*u.km/u.s)


((1210*u.K*u.km/u.s)*(1.36*c.m_p)*(1.82e18*u.s/(u.K*u.km*u.cm**2))).to(c.M_sun/u.pc**2)


((4*u.K*u.km/u.s)*(2*1.36*c.m_p)*(1.97e20*u.s/(u.K*u.km*u.cm**2))).to(c.M_sun/u.pc**2)


#calculating the Sigma_CO equivalent to W_CO = 1.5 K km/s (from Lewis) 
((1.5*u.K*u.km/u.s)*(2*1.36*c.m_p)*(1.97e20*u.s/(u.K*u.km*u.cm**2))).to(c.M_sun/u.pc**2)


#compare the masses when decomposing using the zeroth moment map vs the ppv cube


cb = HI4PI.cube[433:511,b0HI:b1HI,l0HI:l1HI]
sp = dc.spectrum1D(cb)


np.nanstd(HI4PI.cube[390:400,b0HI:b1HI,l0HI:l1HI])


d = Dendrogram.compute(cb,3*13,3*13,100)


d.trunk


bigm = d.trunk[0].get_mask()
for k in d.trunk:
    bigm +=k.get_mask()


plt.imshow(np.nansum(cb,axis=1))


v = HI4PI.vels[450:500]
cb = HI4PI.cube[450:500,b0HI:b1HI,l0HI:l1HI]


for vch in cb:
    

np.nansum(cb,axis=2).shape


yHI.shape


xHI.shape


y.shape,x.shape,z.shape


cb.shape


len(cb)


cbt = cb.transpose()


var=var/1e3


var


dvd = Dendrogram.compute(var,min_npix=500)


pdvd = dvd.plotter()


dvd.trunk


dvd.trunk[0].sorted_leaves()[0]


cat


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.imshow(var,extent=[yCO[0],yCO[-1],xCO[0],xCO[-1]])
for l in dvd.trunk[0].sorted_leaves():
    pdvd.plot_contour(ax,structure=l,colors='white',linewidths=0.5)
    
plt.contour(yCO,xCO,CO_0.transpose(),colors='red',linewidths=0.7,levels=3)
plt.ylim(-130,-160)
plt.xlim(-25,-5)
plt.show()


t = t[0]


len(t.sorted_leaves())


len(t.descendants)


len(t.children)#[0].children[0].children[1].children)


t.children[0].children[0]


for ch in t.descendants:
    if ch.is_leaf:
        mch = ch.get_mask()
        plt.imshow(mch*var,**kw_imshow)
        plt.show()


fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(1,2,2,projection='3d')
z = (var/1e3)*(var/1e3 < 0)
x,y = np.meshgrid(yHI,xHI)
ax.plot_surface(x,y,z,cmap = 'gray',alpha = 0.7,linewidth=0.3)
ax.view_init(azim=30,roll=0,elev=20)


var = np.zeros((500,300))
for l,ll in enumerate(cbt):
    for b,bb in enumerate(ll):
        print(bb.shape)
        print(len(bb))
        print(np.nanargmax(bb))
        vp = np.nanargmax(bb)
        vo = v[vp]
        var[l][b] = vo


plt.plot(cbt[0][0])


fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(1,2,2,projection='3d')
z = np.nansum(cb,axis=2)
x,y = np.meshgrid(yHI, v)

ax.plot_wireframe(x,y,z,cmap = 'gray',alpha = 0.7,linewidth=0.3)
ax.view_init(azim=180,roll=0,elev=20)

ax = fig.add_subplot(1,2,1,projection='3d')
z = np.nansum(cb,axis=1)
x,y = np.meshgrid(xHI, v)

ax.plot_wireframe(x,y,z,cmap = 'gray',alpha = 0.7,linewidth=0.3)
ax.view_init(azim=-0,roll=0,elev=8)


plt.show()


for i in d.trunk:
    imHI = np.nansum(cb*i.get_mask(),axis=0)
    plt.imshow(imHI,**kw_imshow,\
               extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
    plt.contour(xHI,yHI,imHI,levels=[100,500],colors=['white','red'])
    plt.contour(xCO,yCO,CO_0,levels=[3],colors=['orange'],linewidths=[0.9])
    # plt.colorbar()
    plt.ylim(-25,-5)
    plt.show()


imHI = np.nansum(cb*bigm,axis=1)
plt.imshow(imHI,**kw_imshow,\
           extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
# plt.contour(xHI,yHI,imHI,levels=[100,500],colors=['white','red'])
# plt.contour(xCO,yCO,CO_0,levels=[3],colors=['orange'],linewidths=[0.9])
# plt.colorbar()
# plt.ylim(-25,-5)
plt.show()


imHI = np.nansum(cb*bigm,axis=0)
plt.imshow(imHI,**kw_imshow,\
           extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
plt.contour(xHI,yHI,imHI,levels=[100,500],colors=['white','red'])
plt.contour(xCO,yCO,CO_0,levels=[3],colors=['orange'],linewidths=[0.9])
# plt.colorbar()
plt.ylim(-25,-5)
plt.show()


sg = dc.dispersion_1sigma(sp)


ct = np.nanargmax(sp)


vo = HI4PI.vels[ct]/1e3


np.nanargmin(np.abs(HI4PI.vels - 5e4))


np.nanargmin(np.abs(HI4PI.vels + 5e4))


plt.plot(HI4PI.vels/1e3,sp)
plt.hlines(0.5*np.nanmax(sp),vo-7*sg,vo+7*sg)
plt.hlines(0.07*np.nanmax(sp),vo-7*sg,vo+7*sg)

plt.xlim(-100,100)


l0HI = 3700
l1HI = 4200

b0HI = 700
b1HI = 1000

# HI_0 = np.nansum(HI4PI.cube[:,b0HI:b1HI,l0HI:l1HI],axis=0)
xHI = HI4PI.lons[l0HI:l1HI]
yHI = HI4PI.lats[b0HI:b1HI]

l0CO = 2500
l1CO = 2750
b0CO = 60
b1CO = 200

xCO = locmom.lons[l0CO:l1CO]
yCO = locmom.lats[b0CO:b1CO]
CO_0 = np.nansum(locmom.cube[:,b0CO:b1CO,l0CO:l1CO],axis=0)


plt.figure(figsize=(10,10))
plt.imshow(HI_0,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
plt.colorbar(label=r'$W_{HI}\ [K\ km\ s^{-1}$]')

plt.contour(xCO,yCO,CO_0,levels=[1],colors=['white'],linewidths=[0.9])
plt.ylabel(r'$b\ [\degree$]')
plt.xlabel(r'$l\ [\degree$]')
plt.show()


d = Dendrogram.compute(HI_0,min_value=500,min_delta=100)


m = d.trunk[0].get_mask()


dc = Dendrogram.compute(CO_0,min_value=3)#,min_delta=2,min_npix=50)


bigm = dc.trunk[0].get_mask()
for n in dc.trunk[1:]:
    bigm +=n.get_mask()


fig,axs = plt.subplots(2,2,figsize=(10,10))

axs[0,0].imshow(CO_0*bigm,**kw_imshow,extent=[xCO[0],xCO[-1],yCO[0],yCO[-1]])

axs[0,0].set_ylabel(r'$b\ [\degree$]')
axs[0,0].set_xlabel(r'$l\ [\degree$]')
# axs[0,0].xlim(-135,-160)
plt.show()


plt.figure(figsize=(10,10))
# plt.imshow(HI_0*m,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
# plt.contour(xHI,yHI,HI_0,levels=[500,750,1000],\
#             colors=['cyan','orange','red'],linewidths=[0.9])
# plt.colorbar(label=r'$W_{HI}\ [K\ km\ s^{-1}$]')

# plt.contour(xCO,yCO,CO_0,levels=[1.5,2,4,8],colors=['black','red','orange','green'],linewidths=[0.9])
plt.contour(xCO,yCO,CO_0,levels=[2,3,4],colors=['red','orange','cyan'],linewidths=[0.9])

plt.ylabel(r'$b\ [\degree$]')
plt.xlabel(r'$l\ [\degree$]')
plt.xlim(-135,-160)
plt.show()z


plt.figure(figsize=(10,10))
# plt.imshow(HI_0*m,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
# plt.contour(xHI,yHI,HI_0,levels=[500,750,1000],\
#             colors=['cyan','orange','red'],linewidths=[0.9])
# plt.colorbar(label=r'$W_{HI}\ [K\ km\ s^{-1}$]')

# plt.contour(xCO,yCO,CO_0,levels=[1.5,2,4,8],colors=['black','red','orange','green'],linewidths=[0.9])
plt.contour(xCO,yCO,CO_0,levels=[2,3,4],colors=['red','orange','cyan'],linewidths=[0.9])

plt.ylabel(r'$b\ [\degree$]')
plt.xlabel(r'$l\ [\degree$]')
plt.xlim(-135,-160)
plt.show()


d3d = Dendrogram.compute(HI4PI.cube[:,b0HI:b1HI,l0HI:l1HI])


d3d.trunk


plt.figure(figsize=(10,10))
plt.imshow(HI_0*m,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
plt.contour(xHI,yHI,HI_0,levels=[500,750,1000],\
            colors=['cyan','orange','red'],linewidths=[0.9])
plt.colorbar(label=r'$W_{HI}\ [K\ km\ s^{-1}$]')

# plt.contour(xCO,yCO,CO_0,levels=[2],colors=['white'],linewidths=[0.9])
plt.ylabel(r'$b\ [\degree$]')
plt.xlabel(r'$l\ [\degree$]')
plt.show()


plt.figure(figsize=(10,10))
plt.imshow(HI_0*m,**kw_imshow,extent=[xHI[0],xHI[-1],yHI[0],yHI[-1]])
plt.contour(xHI,yHI,HI_0,levels=[500,750,1000],\
            colors=['cyan','orange','red'],linewidths=[0.9])
plt.colorbar(label=r'$W_{HI}\ [K\ km\ s^{-1}$]')

# plt.contour(xCO,yCO,CO_0,levels=[2],colors=['white'],linewidths=[0.9])
plt.ylabel(r'$b\ [\degree$]')
plt.xlabel(r'$l\ [\degree$]')
plt.show()


fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(1,2,2,projection='3d')
z = HI_0
x,y = np.meshgrid(xHI, yHI)

# for a in axs.flat:
#     a = plt.axes(projection='3d')

ax.plot_surface(x,y,z,cmap = 'binary',alpha = 0.6)
ax.view_init(azim=-90,roll=0,elev=45)
ax.contour3D(x,y,z,levels=[500,750,1000],colors=['xkcd:cyan','orange','red'])

xco,yco = np.meshgrid(xCO,yCO)
# ax.contour3D(xco,yco,CO_0,levels=[1.5],colors=['red'])


ax = fig.add_subplot(1,2,1,projection='3d')
ax.plot_surface(x,y,z,cmap = 'binary',alpha = 0.6)
ax.view_init(azim=-40,roll=0,elev=8)
ax.contour3D(x,y,z,levels=[500,750,1000],colors=['xkcd:cyan','orange','red'])


plt.show()


for i in d.trunk:
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
    

plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,750:800,0:50]),linewidth=0.7,label='l=180-175, b=-27.5- -23.5')
plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,750:800,500:550]),linewidth=0.7,label='l=138-134, b=-27.5- -23.5')
plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,1100:1150,0:50]),linewidth=0.7,label='l=138-134, b=-27.5- -23.5')
plt.plot(HI4PI.vels,dc.spectrum1D(HI4PI.cube[:,1100:1150,500:550]),linewidth=0.7,label='l=138-134, b=1.5-5.6')


plt.xlim(-50000,50000)
plt.legend()
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


plt.imshow(np.nansum(trs,axis=1),**kw_imshow)


plt.imshow(np.nansum(ors,axis = 1),**kw_imshow)
plt.colorbar()
