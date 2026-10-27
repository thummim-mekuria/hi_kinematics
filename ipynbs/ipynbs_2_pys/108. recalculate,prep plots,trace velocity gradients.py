#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable

import fun as fns
import datacube as dc

hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data


fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


# #okay now the masses are just a factor of two off from lewis

# masses = np.zeros(22)
# radii = np.zeros(22)
# for i in range(22):
# #     fig,axs = plt.subplots(1,2,figsize=(15,8))
#     vo,vi = ((fcat['vel<km/s>'][i]+100)/1.3).astype('int')
#     lo,li = ((fcat['lon_GMC<deg>'][i]-360)/-0.125).astype('int')
#     bo,bi = ((fcat['lat_GMC<deg>'][i]+50)/0.125).astype('int')
    
    
#     cube = (dht.data*(dht.data>0.2*3))[vo:vi,bo:bi,lo:li]*u.K
    
#     d = fcat['distance<pc>'][i]*u.pc
    
#     m,a = mass_n_area(cube,d,'CO')
    
#     masses[i] = m.to(u.kg).value
#     radii[i] = np.sqrt(a.to(u.pc**2)/np.pi).value
    
    
# fcat['mass_GMC<kg>'] = masses
# fcat['radius_GMC<pc>'] = radii

# s, ra = fns.accumulation_radius(masses*u.kg,fcat['distance<pc>']*u.pc)
# sh = s.value
# fcat['radius_HI<pc>'] = ra
# fcat.write('re_expanded_cat.fits',overwrite=True)


lcenter = np.nanmean(fcat['lon_GMC<deg>'],axis=1)
bcenter = np.nanmean(fcat['lat_GMC<deg>'],axis=1)


sh = (fcat['radius_HI<pc>']/fcat['distance<pc>'])*u.rad.to(u.deg)


lpixhi = ((np.array([lcenter+sh,lcenter-sh]).transpose() -360)/-0.125).astype('int')


bpixhi = ((np.array([bcenter-sh,bcenter+sh]).transpose()+50)/0.125).astype('int')


lpixgmc = ((fcat['lon_GMC<deg>'] -360)/-0.125).astype('int')
bpixgmc = ((fcat['lat_GMC<deg>'] +50)/0.125).astype('int')


vpix = ((fcat['vel<km/s>']+100)/1.3).astype('int')


h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


l0g,l1g = lpixgmc.transpose()
b0g,b1g = bpixgmc.transpose()

v0,v1 = vpix.transpose()

l0h,l1h = lpixhi.transpose()
b0h,b1h = bpixhi.transpose()


vs = [np.arange(v0[i],v1[i])*1.3-100 for i in range(22)]


h0s = [np.nansum(h3sig[v0[i]:v1[i],
                       b0h[i]:b1h[i],
                       l0h[i]:l1h[i]],
                 axis=0)*1.82e18*1.3 for i in range(22)]


c0s = [np.nansum(d3sig[v0[i]:v1[i],
                       b0g[i]:b1g[i],
                       l0g[i]:l1g[i]],
                 axis=0)*1.97e20*1.3 for i in range(22)]


h1s = [fns.get_mom1(h3sig[v0[i]:v1[i],
                       b0h[i]:b1h[i],
                       l0h[i]:l1h[i]],
                 vs[i])for i in range(22)]


c1s = [fns.get_mom1(d3sig[v0[i]:v1[i],
                       b0g[i]:b1g[i],
                       l0g[i]:l1g[i]],
                 vs[i])for i in range(22)]


delvmean = np.zeros(22)


v0[i]


delvmean2 = np.zeros(22)


plt.scatter(delvmean,delvmean2)
plt.plot(delvmean,delvmean,marker='+')


for i in range(22):
    hsp = np.nanmean(h3sig[v0[i]:v1[i],b0h[i]:b1h[i],l0h[i]:l1h[i]],axis=(1,2))
    print(np.nanmean(h1s[i]),
          fns.get_mom1(hsp,np.arange(v0[i],v1[i])*1.3-100))
    
    csp = np.nanmean(d3sig[v0[i]:v1[i],b0h[i]:b1h[i],l0h[i]:l1h[i]],axis=(1,2))

    print(np.nanmean(c1s[i]),
          fns.get_mom1(csp,np.arange(v0[i],v1[i])*1.3-100))
    
    delvmean2[i]=fns.get_mom1(csp,np.arange(v0[i],v1[i])*1.3-100)-\
    fns.get_mom1(hsp,np.arange(v0[i],v1[i])*1.3-100)
    print('\n\n')


for i in range(22):
    plt.hist(np.ravel(h1s[i]),alpha=0.5,density=True)
    plt.hist(np.ravel(c1s[i]),alpha=0.5,density=True)
    plt.vlines(np.nanmean(h1s[i]),0,1)
    plt.vlines(np.nanmean(c1s[i]),0,1,color='orange')
    plt.title(str(i)+': '+fcat['name'][i])
    
    print()
    
    delvmean[i] = np.nanmean(c1s[i]) - np.nanmean(h1s[i])
#     plt.yscale('log')
    plt.show()


plt.hist(delvmean,bins=10)
plt.show()


plt.scatter(delvmean,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)),s=5)
plt.scatter(delvmean,np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)),s=5)
plt.scatter(delvmean,np.log10((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass)),s=5)


plt.scatter(delvmean,fcat['spec_ang_mom_GMC<pc.km/s>'],s=5)
plt.loglog()


plt.scatter(np.abs(delvmean),
            fcat['spec_ang_mom_GMC<pc.km/s>']/fcat['spec_ang_mom_HI<pc.km/s>'],
           s=5)
plt.loglog()
plt.show()


plt.scatter(delvmean,fcat['vel_dispersion_GMC<km/s>'],s=5)
plt.scatter(delvmean,fcat['vel_dispersion_HI<km/s>'],s=5)


plt.scatter(fcat['vel_dispersion_GMC<km/s>'],fcat['vel_dispersion_HI<km/s>'])
# plt.loglog()
# plt.show()


fcat


for i in range(22):
    print(fcat['name'][i],np.nanstd(h1s[i]))


(1.3/atd)*


1.3/5


1.3/1.6


aligned velocity gradient arrows can trace the path of a coherent flow (eg hercules)
- if successive pixels have velocities that are monotonically increasing/decreasing, then
the velocity vectors at those pixels will be lined up head-to-tail

this breaks down for a largescale velocity gradient (eg california)


^decribe in more words, starting with how i calculate the vectors
give all the context neceesarry


for i in [8]: #just hercules

    vo,vi = ((fcat['vel<km/s>'][i]+100)/1.3).astype('int')
    loc,lic = lpixgmc[i]
    boc,bic = bpixgmc[i]
    
    loh,lih = lpixhi[i]
    boh,bih = bpixhi[i]

    vs = np.arange(vi-vo)*1.3+fcat['vel<km/s>'][i][0]
    
    this_hi = (hi4pi.data*(hi4pi.data>0.043*3))[vo:vi,boh:bih,loh:lih]
    this_co = (dht.data*(dht.data>0.2*3))[vo:vi,boc:bic,loc:lic]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,vs)
    c1 = fns.get_mom1(this_co,vs)
        

clv = np.nansum(this_co,axis=1)*1.82e18*1.3


dvdl,dvdb = np.gradient(c1)
dvmag = np.sqrt(dvdl**2 + dvdb**2)

fig,axs = plt.subplots(1,2,figsize=(20,8))

x,y = np.meshgrid(np.arange(lic-loc)*-0.125 + fcat['lon_GMC<deg>'][i][0],\
                  np.arange(bic-boc)*0.125 + fcat['lat_GMC<deg>'][i][0])


axs[0].contourf(c1,extent=[loc*-0.125+360,lic*-0.125+360,boc*0.125-50,bic*0.125-50],
                levels=10,cmap='PiYG')
axs[0].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.5*(lic-loc),cmap='binary_r')


# axs[0].set_xlim(fcat['lon_GMC<deg>'][i])

axs[1].contourf(clv,extent=[loc*-0.125+360,lic*-0.125+360,vo*1.3-100,vi*1.3-100],
                levels=10,cmap='YlGnBu_r')

axs[1].set_ylim(-5,20)
# axs[1].set_xlim(fcat['lon_GMC<deg>'][i])


for i in [1]: #just california

    vo,vi = ((fcat['vel<km/s>'][i]+100)/1.3).astype('int')
    loc,lic = lpixgmc[i]
    boc,bic = bpixgmc[i]
    
    loh,lih = lpixhi[i]
    boh,bih = bpixhi[i]

    vs = np.arange(vi-vo)*1.3+fcat['vel<km/s>'][i][0]
    
    this_hi = (hi4pi.data*(hi4pi.data>0.043*3))[vo:vi,boh:bih,loh:lih]
    this_co = (dht.data*(dht.data>0.2*3))[vo:vi,boc:bic,loc:lic]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,vs)
    c1 = fns.get_mom1(this_co,vs)
        

clv = np.nansum(this_co,axis=1)*1.82e18*1.3


dvdl,dvdb = np.gradient(c1)
dvmag = np.sqrt(dvdl**2 + dvdb**2)

x,y = np.meshgrid(np.arange(lic-loc)*-0.125 + fcat['lon_GMC<deg>'][i][0],\
                  np.arange(bic-boc)*0.125 + fcat['lat_GMC<deg>'][i][0])

fig,axs = plt.subplots(1,2,figsize=(20,8))

axs[0].contourf(c1,extent=[loc*-0.125+360,lic*-0.125+360,boc*0.125-50,bic*0.125-50],
                levels=10,cmap='PiYG')
axs[0].quiver(x[::2,::2],y[::2,::2],(dvdl/dvmag)[::2,::2],
              (dvdb/dvmag)[::2,::2],dvmag[::2,::2],scale=0.3*(lic-loc),cmap='binary_r')

# axs[0].set_xlim(fcat['lon_GMC<deg>'][i])


axs[1].contourf(clv,extent=[loc*-0.125+360,lic*-0.125+360,vo*1.3-100,vi*1.3-100],
                levels=10,cmap='YlGnBu_r')

# axs[1].set_xlim(fcat['lon_GMC<deg>'][i])

axs[1].set_ylim(-10,5)


np.nanmean(fcat['lon_GMC<deg>'][4]),np.nanmean(fcat['lat_GMC<deg>'][4])


(fcat['radius_HI<pc>'][4]/fcat['distance<pc>'][4])*u.rad.to(u.deg)


plt.scatter(np.log10(np.ravel(h0)),np.ravel(ddmag),s=4)
# plt.loglog()
plt.yscale('log')


x.shape


for i in range(22): #just cma
    
    
    vo,vi = ((fcat['vel<km/s>'][i]+100)/1.3).astype('int')
    loc,lic = lpixgmc[i]
    boc,bic = bpixgmc[i]
    
    loh,lih = lpixhi[i]
    boh,bih = bpixhi[i]

    vs = np.arange(vi-vo)*1.3+fcat['vel<km/s>'][i][0]
    
    this_hi = (hi4pi.data*(hi4pi.data>0.043*3))[vo:vi,boh:bih,loh:lih]
    this_co = (dht.data*(dht.data>0.2*3))[vo:vi,boc:bic,loc:lic]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,vs)
    c1 = fns.get_mom1(this_co,vs)

    fig,axs = plt.subplots(2,2,figsize=(20,20),sharey=True)

    h0im = axs[0,0].contourf(np.log10(h0),extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                    levels=np.arange(20,22,0.05),cmap='YlGnBu_r')

    h1im = axs[0,1].contourf(h1,extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                    levels=10,cmap='RdGy')


    dvdl,dvdb = np.gradient(h1)
    dvmag = np.sqrt(dvdl**2 + dvdb**2)

    dvdl2,dvdb2 = np.gradient(np.abs(h1-13))
    dvmag2 = np.sqrt(dvdl2**2 + dvdb2**2)


    lc = lcenter[i] 
    bc = bcenter[i] 
    
    shi = sh[i]
    
    x,y = np.meshgrid(np.arange(loh,lih)*-0.125 + 360,
                      np.arange(boh,bih)*0.125 -50)


    axs[0,1].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.7*(lih-loh),cmap='binary')
    # axs[0,0].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.7*(lih-loh),cmap='binary')


    dddl,dddb = np.gradient(h0/1e21)
    ddmag = np.sqrt(dddl**2 + dddb**2)   


    dh0im = axs[1,0].contourf(ddmag,extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
                    levels=10,cmap='YlGnBu_r')
#     axs[0,0].contour(ddmag,extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
#                     levels=10,cmap='YlGnBu_r',linewidths=0.7)
    
    axs[1,1].quiver(x,y,dvdl2/dvmag2,dvdb2/dvmag2,dvmag2,scale=0.7*(lih-loh),cmap='binary')
    
    plt.show()


h1gal = fns.get_mom1(hi4pi.data*(hi4pi.data>0.043*3),np.arange(155)*1.3-100)


h3sig = (hi4pi.data*(hi4pi.data>0.043*3))
d3sig = (dht.data*(dht.data>0.2*3))

h0 = np.nansum(h3sig,axis=0)*1.82e18*1.3
c0 = np.nansum(d3sig,axis=0)*1.97e20*1.3

    
vs = np.arange(155)*1.3-100
# h1 = fns.get_mom1(h3sig,vs)
# c1 = fns.get_mom1(d3sig,vs)

fig,axs = plt.subplots(2,1,figsize=(20,20),sharey=True)

h0im = axs[0].contourf(np.log10(h0),extent=[360,-360,-50,50],
                levels=np.arange(20,22,0.05),cmap='YlGnBu_r')

# h1im = axs[0,1].contourf(h1,extent=[360,-360,-50,50],
#                     levels=10,cmap='RdGy')


# dvdl,dvdb = np.gradient(h1)
# dvmag = np.sqrt(dvdl**2 + dvdb**2)

# dvdl2,dvdb2 = np.gradient(np.abs(h1-13))
# dvmag2 = np.sqrt(dvdl2**2 + dvdb2**2)


x,y = np.meshgrid(np.arange(5760)*-0.125 + 360,
                  np.arange(800)*0.125 -50)


# axs[0,1].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.7*(lih-loh),cmap='binary')
    # axs[0,0].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.7*(lih-loh),cmap='binary')


dddl,dddb = np.gradient(h0/1e21)
ddmag = np.sqrt(dddl**2 + dddb**2)   


dh0im = axs[1].contourf(np.sqrt(ddmag),extent=[360,-360,-50,50],norm='asinh',
                levels=np.arange(0,1,0.05),cmap='YlGnBu_r')

axs[1].contour(np.log10(c0),extent=[360,-360,-50,50],colors='white',linewidths=0.6,
                levels=np.arange(20,22,0.5))
#     axs[0,0].contour(ddmag,extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
#                     levels=10,cmap='YlGnBu_r',linewidths=0.7)
    
# axs[1,1].quiver(x,y,dvdl2/dvmag2,dvdb2/dvmag2,dvmag2,scale=0.7*(lih-loh),cmap='binary')

fig.colorbar(dh0im)

for ax in axs.flat:
    ax.set_xlim(180,-180)
    ax.set_ylim(-50,50)

plt.show()


h3sig = (hi4pi.data*(hi4pi.data>0.043*3))
d3sig = (dht.data*(dht.data>0.2*3))

h0 = np.nansum(h3sig,axis=0)*1.82e18*1.3
c0 = np.nansum(d3sig,axis=0)*1.97e20*1.3

    
vs = np.arange(155)*1.3-100
# h1 = fns.get_mom1(h3sig,vs)
# c1 = fns.get_mom1(d3sig,vs)

fig,axs = plt.subplots(2,1,figsize=(20,20),sharey=True)

h0im = axs[0].contourf(np.log10(h0),extent=[360,-360,-50,50],
                levels=np.arange(20,22,0.05),cmap='YlGnBu_r')

# h1im = axs[0,1].contourf(h1,extent=[360,-360,-50,50],
#                     levels=10,cmap='RdGy')


# dvdl,dvdb = np.gradient(h1)
# dvmag = np.sqrt(dvdl**2 + dvdb**2)

# dvdl2,dvdb2 = np.gradient(np.abs(h1-13))
# dvmag2 = np.sqrt(dvdl2**2 + dvdb2**2)


x,y = np.meshgrid(np.arange(5760)*-0.125 + 360,
                  np.arange(800)*0.125 -50)


# axs[0,1].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.7*(lih-loh),cmap='binary')
    # axs[0,0].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.7*(lih-loh),cmap='binary')


dddl,dddb = np.gradient(h0/1e21)
ddmag = np.sqrt(dddl**2 + dddb**2)   


dh0im = axs[1].contourf(dddl+dddb,extent=[360,-360,-50,50],norm='asinh',
                levels=np.arange(0,1,0.05),cmap='YlGnBu_r')

axs[1].contour(np.log10(c0),extent=[360,-360,-50,50],colors='white',linewidths=0.6,
                levels=np.arange(20,22,0.5))
#     axs[0,0].contour(ddmag,extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
#                     levels=10,cmap='YlGnBu_r',linewidths=0.7)
    
# axs[1,1].quiver(x,y,dvdl2/dvmag2,dvdb2/dvmag2,dvmag2,scale=0.7*(lih-loh),cmap='binary')

fig.colorbar(dh0im)

for ax in axs.flat:
    ax.set_xlim(180,-180)
    ax.set_ylim(-50,50)

plt.show()


for i in range(22):
    fig,axs = plt.subplots(1,2,figsize=(20,8),sharey=True)
    vo,vi = ((fcat['vel<km/s>'][i]+100)/1.3).astype('int')
    loc,lic = lpixgmc[i]
    boc,bic = bpixgmc[i]
    
    loh,lih = lpixhi[i]
    boh,bih = bpixhi[i]

    vs = np.arange(vi-vo)*1.3+fcat['vel<km/s>'][i][0]
    
    this_hi = (hi4pi.data*(hi4pi.data>0.043*3))[vo:vi,boh:bih,loh:lih]
    this_co = (dht.data*(dht.data>0.2*3))[vo:vi,boc:bic,loc:lic]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,vs)
    c1 = fns.get_mom1(this_co,vs)
    
    
    dvdl,dvdb = np.gradient(h1)
    dvmag = np.sqrt(dvdl**2 + dvdb**2)
    
    x,y = np.meshgrid(np.arange(lih-loh)*-0.125 + fcat['lon_GMC<deg>'][i][0],
                      np.arange(bih-boh)*0.125 +fcat['lat_GMC<deg>'][i][0])
    axs[1].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.3*(li-lo),cmap='binary')
    
    
    dvdl,dvdb = np.gradient(h1)
    dvmag = np.sqrt(dvdl**2 + dvdb**2)
    
    x,y = np.meshgrid(np.arange(lih-loh)*-0.125 + fcat['lon_GMC<deg>'][i][0],
                      np.arange(bih-boh)*0.125 +fcat['lat_GMC<deg>'][i][0])
    axs[1].quiver(x,y,dvdl/dvmag,dvdb/dvmag,dvmag,scale=0.3*(li-lo),cmap='binary')
    
    
    fig.suptitle(fcat['name'][i])
    
    fig.colorbar(ax0,ax=axs[0])
    
    
    plt.show()
