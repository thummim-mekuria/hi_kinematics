#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
import datacube as dc
import fun as fns


cat = QTable(fits.open('post_flash_cat.fits')[1].data)


hi4pi,dht = fits.open('both_same_grid_720d.fits')


h3sig = hi4pi.data*(hi4pi.data > 0.043*3)
d3sig = dht.data*(dht.data > 0.2*3)


dust = fits.open('mean_and_std_lbd.fits')[1].data


dlb = np.nansum(dust,axis=0)
dld = np.nansum(dust,axis=1)
dbd = np.nansum(dust,axis=2)


hlb = np.nansum(h3sig,axis=0)*1.3*1.82e18
hlv = np.nansum(h3sig,axis=1)*0.125*1.82e18
hbv = np.nansum(h3sig,axis=2)*0.125*0.5*1.82e18


clb = np.nansum(d3sig,axis=0)*1.3*1.97e20
clv = np.nansum(d3sig,axis=1)*0.125*1.97e20
cbv = np.nansum(d3sig,axis=2)*0.125*0.5*1.97e20


plt.contour(clb)
plt.colorbar()


plt.figure(figsize=(20,15))
plt.contourf(dlb*2.8,levels=np.arange(5e-2,5,5e-2),extent=[180,-180,-90,90])
plt.contour(dlb*2.8,levels=[0.2],colors='red',extent=[180,-180,-90,90])

plt.contour(np.log10(clb),levels=[20.5,21],colors=['orange','white'],extent=[360,-360,-50,50])

plt.xlim(180,-180)


plt.figure(figsize=(20,15))
plt.contourf(dld*2.8,levels=np.arange(5e-2,5,5e-2),extent=[180,-180,68,1250])
plt.contour(dld*2.8,levels=[0.2],colors='red',extent=[180,-180,68,1250])

# plt.contour(np.log10(clb),levels=[20.5,21],colors=['orange','white'],extent=[360,-360,-50,50])

plt.xlim(180,-180)


plt.figure(figsize=(20,15))
plt.contourf(dbd*2.8,levels=np.arange(5e-2,5,5e-2),extent=[-90,90,68,1250])
plt.contour(dbd*2.8,levels=[0.2],colors='red',extent=[-90,90,68,1250])


dld_posb = np.nansum(dust[:,393:,:],axis=1)
dld_negb = np.nansum(dust[:,:393,:],axis=1)


plt.figure(figsize=(20,15))
plt.contourf(dld_posb*2.8,levels=np.arange(5e-2,5,5e-2),extent=[180,-180,68,1250])
plt.contour(dld_posb*2.8,levels=[0.2],colors='red',extent=[180,-180,68,1250])

# plt.contour(np.log10(clb),levels=[20.5,21],colors=['orange','white'],extent=[360,-360,-50,50])

plt.xlim(180,-180)


plt.figure(figsize=(20,15))
plt.contourf(dld_negb*2.8,levels=np.arange(5e-2,5,5e-2),extent=[180,-180,68,1250])
plt.contour(dld_negb*2.8,levels=[0.2],colors='red',extent=[180,-180,68,1250])

# plt.contour(np.log10(clb),levels=[20.5,21],colors=['orange','white'],extent=[360,-360,-50,50])

plt.xlim(180,-180)


hmom1 = fns.get_mom1(h3sig,np.arange(155)*1.3-100)
hmom2 = np.sqrt(fns.get_momN(h3sig,np.arange(155)*1.3-100,2))
hmom3 = fns.get_momN(h3sig,np.arange(155)*1.3-100,3)
hmom4 = fns.get_momN(h3sig,np.arange(155)*1.3-100,4)


plt.figure(figsize=(20,10))
plt.contourf(hmom1,extent=[360,-360,-50,50],cmap='YlGnBu',levels=np.linspace(-30,30,10))
plt.colorbar()

plt.xlim(180,-180)


plt.figure(figsize=(20,10))
plt.contourf(hmom3,extent=[360,-360,-50,50],cmap='YlGnBu',levels=np.linspace(-2,2,10))
plt.colorbar()

plt.xlim(180,-180)
plt.vlines([-90,-60,-45,45,60,90],-50,50,color=['red','orange','green','green','orange','red'])


plt.figure(figsize=(20,10))
plt.contourf(hmom4-3,extent=[360,-360,-50,50],cmap='YlGnBu',levels=np.linspace(-3,0,10))
plt.contourf(hmom4-3,extent=[360,-360,-50,50],cmap='magma',levels=np.linspace(-1,1,10))

plt.colorbar()

plt.xlim(180,-180)
# plt.vlines([-90,-60,-45,45,60,90],-50,50,color=['red','orange','green','green','orange','red'])


plt.figure(figsize=(20,10))
plt.contourf(hmom4-3,extent=[360,-360,-50,50],cmap='YlGnBu',levels=np.linspace(-1,1,10))
plt.colorbar()

plt.xlim(180,-180)
# plt.vlines([-90,-60,-45,45,60,90],-50,50,color=['red','orange','green','green','orange','red'])


plt.figure(figsize=(20,10))
plt.contourf(hmom4-3,extent=[360,-360,-50,50],cmap='YlGnBu',levels=np.linspace(-0,3,10))
plt.colorbar()

plt.xlim(180,-180)
# plt.vlines([-90,-60,-45,45,60,90],-50,50,color=['red','orange','green','green','orange','red'])


plt.figure(figsize=(20,10))
plt.contourf(hmom4-3,extent=[360,-360,-50,50],cmap='YlGnBu',levels=np.linspace(3,10,10))
plt.colorbar()

plt.xlim(180,-180)
# plt.vlines([-90,-60,-45,45,60,90],-50,50,color=['red','orange','green','green','orange','red'])


h3sig.shape


hmom1.shape


np.ravel(hmom1).shape


vls = np.arange(155)*1.3-100


mns = np.ravel(hmom1)
sgs = np.ravel(hmom2)
skw = np.ravel(hmom3)
krt = np.ravel(hmom4)


sp_list = h3sig.transpose().swapaxes(0,1).reshape(-1,155)


lns, lts = np.meshgrid(np.arange(360,-360,-0.125),np.arange(-50,50,0.125))


lns = np.ravel(lns)
lts = np.ravel(lts)


hmom5 = fns.get_momN(h3sig,vls,5)


fth = np.ravel(hmom5)


plt.figure(figsize=(20,10))
plt.contourf(hmom5,extent=[360,-360,-50,50],levels=np.linspace(-50,50,10))
plt.xlim(180,-180)
plt.colorbar()


for i,sp in zip(range(4608000)[::100],sp_list[::100]):
    plt.plot(vls,sp)
    plt.vlines(mns[i],0,np.nanmax(sp),color='red')
    plt.annotate('('+str(lns[i])+','+str(lts[i])+')'+
                         '\nsig: '+str(sgs[i])+
                     '\nskw: '+str(skw[i])+
                     '\nkrt: '+str(krt[i])+
                    '\n5th: '+str(fth[i]),
                        (-100,0.8*np.nanmax(sp)))
    plt.show()


46080/1e3


from astropy.modeling import models


models.Gaussian1D(1.,0.,1.)


for i,sp in zip(range(4608000)[::10000],sp_list[::10000]):
    plt.plot(vls,sp)
    
    mdl = models.Gaussian1D(np.nanmax(sp),mns[i],sgs[i])
    mdl2 = models.Lorentz1D(np.nanmax(sp),mns[i],sgs[i])
    plt.plot(vls,mdl(vls))
    plt.plot(vls,sp-mdl(vls))
    
    pos_df = (sp-mdl(vls)>0)*(sp-mdl(vls))
    plt.annotate(round(100*np.nansum(pos_df)/np.nansum(sp),2),(-100,0.8*np.nanmax(sp)))
    plt.annotate(round(100*np.nansum(sp-mdl(vls))/np.nansum(sp),2),(-100,0.7*np.nanmax(sp)))
    plt.annotate(round(100*np.nansum(sp-mdl2(vls))/np.nansum(sp),2),(-100,0.7*np.nanmax(sp)))
#     plt.vlines(mns[i],0,np.nanmax(sp),color='red')
    
    
#     plt.annotate('('+str(lns[i])+','+str(lts[i])+')'+
#                          '\nsig: '+str(sgs[i])+
#                      '\nskw: '+str(skw[i])+
#                      '\nkrt: '+str(krt[i])+
#                     '\n5th: '+str(fth[i]),
#                         (-100,0.8*np.nanmax(sp)))
    plt.show()


for i,sp in zip(range(4608000)[::10000],sp_list[::10000]):
    plt.plot(vls,sp)
    
    mdl = models.Gaussian1D(np.nanmax(sp),mns[i],sgs[i])
    plt.plot(vls,mdl(vls))
    plt.plot(vls,sp-mdl(vls))
    
    plt.annotate(round(100*np.nansum(sp-mdl(vls))/np.nansum(sp),2),(-100,0.8*np.nanmax(sp)))
#     plt.vlines(mns[i],0,np.nanmax(sp),color='red')
    
    
#     plt.annotate('('+str(lns[i])+','+str(lts[i])+')'+
#                          '\nsig: '+str(sgs[i])+
#                      '\nskw: '+str(skw[i])+
#                      '\nkrt: '+str(krt[i])+
#                     '\n5th: '+str(fth[i]),
#                         (-100,0.8*np.nanmax(sp)))
    plt.show()


fk_cube = np.arange(1000).reshape(20,5,10)


fk0 = np.nansum(fk_cube,axis=0)


plt.contourf(fk0)
plt.colorbar()


rshpd = fk_cube.transpose().reshape(-1,20)


fk02 = np.nansum(rshpd,axis=1).reshape(5,10)


plt.contourf(fk02)
plt.colorbar()


rshpd = fk_cube.transpose().swapaxes(0,1).reshape(-1,20)


fk03 = np.nansum(rshpd,axis=1).reshape(5,10)


plt.contourf(fk03)
plt.colorbar()


fk03==fk0


np.array(np.meshgrid(np.arange(20),np.arange(5),np.arange(10),indexing='ij'))


ahelp(np.meshgrid)


786/2


dust.shape


dust[:,:393,:].shape


dust[:,393:,:].shape


cat['radius_GMC<pc>']


1182+68
