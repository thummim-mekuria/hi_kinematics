#!/usr/bin/env python
# coding: utf-8


#main takeaway of this notebook -> narrow down to +-30km/s (or 40 being generous)
#this should do a lot to avoid multiple distinct clouds along the line of sight (<- assumes that a single
#cloud/system wouldnt have two widely separated peaks in v ... but havent proved that)
#^ which is a separate problem to blending (ie when velocity peaks are very close & overlapping)
#recalcualte first & second moment maps + do the model each spectrum as a gaussian thing after ^
#then iterate again to see which lines of sight have the largest residuals -> indicative of hisa maybe? 
#complex stuff? ... at the very least, this way we can be sure that the velocities from the first moment 
#map are the actual peak velocities along that line of sight


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import datacube as dc
import fun as fns

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting


hi4pi,dht = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


h0 = np.nansum(hi4pi.data*(hi4pi.data>0.043*5), 0)*1.3*1.82e18
c0 = np.nansum(dht.data  *(dht.data  >0.2*3),   0)*1.3*1.97e20


h1 = fns.get_mom1(hi4pi.data*(hi4pi.data>0.043*5 ),[0,155])
c1 = fns.get_mom1(dht.data  *(dht.data  >0.2*3),   [0,155])


h2 = fns.get_momN(hi4pi.data*(hi4pi.data>0.043*5),[0,155],2)
c2 = fns.get_momN(dht.data*(dht.data    >0.2*3),  [0,155],2)


amplitude = np.nanmax(hi4pi.data*(hi4pi.data>0.043*5),axis=0)
mean = h1
stddv = np.sqrt(h2)


ams = np.ravel(amplitude)
mns = np.ravel(mean)
sds = np.ravel(stddv)


gaus_spec = [model_specs[i](np.arange(-100,101,1.3)) for i in range(len(ams))]


gausmodel_hi4pi = np.array(gaus_spec)


ppvgausmodelh = gausmodel_hi4pi.reshape(800,5760,155).transpose().swapaxes(1,2)


gaus0 = np.nansum(ppvgausmodelh,axis=0)


plt.hist(np.ravel((gaus0*1.3*1.82e18)/h0),bins=100)
# plt.loglog()
plt.yscale('log')
plt.vlines(1,1,1e6,color='red')
plt.show()


plt.scatter(h0,np.ravel(h0-(gaus0*1.3*1.82e18)),s=1e-3,alpha=0.2)
# plt.scatter(h0,-1*np.ravel(h0-(gaus0*1.3*1.82e18)),s=1e-3,alpha=0.2)
plt.vlines([6e19,4e20],1e15,1e21)
plt.plot(10**(np.arange(19,22,0.5)),0.3*10**(np.arange(19,22,0.5)))
plt.loglog()


zcat['b_pix_HI'][10]


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi.data.transpose().swapaxes(0,1).reshape(-1,155)[i]),
plt.show()) for i in range(0,len(gaus_spec))
 if np.nanargmax(gaus_spec[i])<40]


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi.data.transpose().swapaxes(0,1).reshape(-1,155)[i]),
plt.show()) for i in range(0,len(gaus_spec))
 if np.nanstd(gaus_spec[i])>20]


[(plt.plot(gaus_spec[i]),
plt.plot(hi4pi.data.transpose().swapaxes(0,1).reshape(-1,155)[i]),
plt.show()) for i in range(0,len(gaus_spec))
 if 3*np.nansum(hi4pi.data.transpose().swapaxes(0,1).reshape(-1,155)[i]) <
np.nansum(gaus_spec[i])]


smol_neg = np.nansum((hi4pi.data*(hi4pi.data>0.43*5))[:50,:,:],axis=0)


smol_pos = np.nansum((hi4pi.data*(hi4pi.data>0.43*5))[120:,:,:],axis=0)


(35+100)/1.3


#narrow the velocity channles down to -30 to +30 ... or something along those lines, in order
#to cut out the high velocity peaks stacked along the line of sight (look at spectra for when gaus0 > 3xh0)
#^ whats happening here is that the second moment / variance / standard deviation is being
#over-estimated for the spectrum, since there are peaks at well separated velocities


plt.figure(figsize=(20,10))

# plt.contourf(h1,extent=[360,-360,-50,50],levels=np.arange(-100,-40,5))
# plt.contourf(h1,extent=[360,-360,-50,50],levels=np.arange(40,100,5))
plt.contour(c1,extent=[360,-360,-50,50],levels=np.arange(-100,-35,5))
plt.contour(c1,extent=[360,-360,-50,50],levels=np.arange(35,100,5))

plt.xlim(180,-180)


50*1.3-100


120*1.3-100


plt.figure(figsize=(20,10))

plt.imshow(smol_neg,extent=[360,-360,-50,50],origin='lower',norm='log')
plt.xlim(180,-180)
plt.colorbar()


plt.figure(figsize=(20,10))

plt.imshow(smol_pos,extent=[360,-360,-50,50],origin='lower',norm='log')
plt.xlim(180,-180)
plt.colorbar()


plt.figure(figsize=(20,10))

plt.imshow(smol_neg,extent=[360,-360,-50,50],origin='lower',norm='asinh',cmap=dc.cmap_min0('Oranges'))
plt.imshow(smol_pos,extent=[360,-360,-50,50],origin='lower',norm='asinh',cmap=dc.cmap_min0('Blues'))

plt.xlim(210,-210)
# plt.colorbar()


plt.figure(figsize=(20,10))

plt.contourf((gaus0*1.3*1.82e18)/h0,levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contourf((gaus0*1.3*1.82e18)/h0,levels=[3,5],colors=['black','blue'],
          origin='lower',extent=[360,-360,-50,50],linewidths=0.7)
plt.contour(smol_neg,extent=[360,-360,-50,50],colors='white',levels=[10])
plt.contour(smol_pos,extent=[360,-360,-50,50],colors='white',levels=[10])

# plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='red',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(h0/(gaus0*1.3*1.82e18),levels=10,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contour(h0/(gaus0*1.3*1.82e18),levels=[1],colors='black',
          origin='lower',extent=[360,-360,-50,50])
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='red',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(h0-(gaus0*1.3*1.82e18),levels=np.array([-1,-0.5,-0.05,0.05,0.5,1])*1e21,
          origin='lower',extent=[360,-360,-50,50],cmap='RdYlGn_r')
plt.colorbar(orientation='horizontal')

plt.contour(h0-gaus0*1.3*1.82e18,levels=[0],colors='black',
          origin='lower',extent=[360,-360,-50,50])
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='red',extent=[360,-360,-50,50])
plt.xlim(180,-180)
# plt.ylim(-30,30)


plt.figure(figsize=(20,10))

plt.contourf(h0,levels=[4e20,1e21,5e21,1e22,1e23,1e24],
          origin='lower',extent=[360,-360,-50,50],norm='log')
plt.colorbar(orientation='horizontal')

# plt.contour(h0-gaus0*1.3*1.82e18,levels=[0],colors='orange',
#           origin='lower',extent=[360,-360,-50,50])
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='red',extent=[360,-360,-50,50])
plt.xlim(250,-180)


plt.figure(figsize=(20,10))

plt.contourf(gaus0*1.3*1.82e18,levels=[4e20,1e21,5e21,1e22,1e23,1e24],
          origin='lower',extent=[360,-360,-50,50],norm='log')
plt.colorbar(orientation='horizontal')

# plt.contour(h0-gaus0*1.3*1.82e18,levels=[0],colors='orange',
#           origin='lower',extent=[360,-360,-50,50])
plt.contour(c0,levels=np.array([0.5,1])*1e21,linewidths=0.5,colors='red',extent=[360,-360,-50,50])
plt.xlim(250,-180)


plt.figure(figsize=(20,10))
plt.contourf(np.sqrt(h2),levels=10)
#            origin='lower',extent=[360,-360,-50,50])
# plt.xlim(180,-180)
plt.colorbar(orientation='horizontal')


plt.imshow(np.ravel(amplitude).reshape(800,5760),
           origin='lower',extent=[360,-360,-50,50])
plt.xlim(180,-180)


plt.figure(figsize=(20,10))
plt.imshow(np.log10(h0),origin='lower',extent=[360,-360,-50,50],cmap=dc.cmap_min0('viridis'))
plt.contourf(np.arange(5760)*-0.125+360,np.arange(800)*0.125-50,c0,levels=[4e20,1e21],colors='black')
plt.xlim(180,-180)


plt.hist(np.ravel(np.sqrt(h2)),bins=100)
plt.show()


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()


plt.figure(figsize=(20,10))
plt.imshow(np.log10(h0),origin='lower',extent=[360,-360,-50,50],cmap=dc.cmap_min0('viridis'))
plt.contourf(np.arange(5760)*-0.125+360,np.arange(800)*0.125-50,c0,levels=[4e20,1e21],colors='black')
plt.xlim(180,-180)


from astropy.visualization import ZScaleInterval
zlim = ZScaleInterval().get_limits


#gonna go with fitting l & b separately
Omegas = np.zeros(17)
thetas = np.zeros(17)
for i in range(17):
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    
    loc,lic = zcat['l_pix_GMC'][i]
    boc,bic = zcat['b_pix_GMC'][i]
    

    hl,hb = np.meshgrid(np.arange(loh,lih),np.arange(boh,bih))
    
    
    hdf = pd.DataFrame(np.array((
            np.ravel(hl),np.ravel(hb),
            np.ravel(h0[boh:bih,loh:lih]),
            np.ravel(h1[boh:bih,loh:lih]),
            np.ravel(h2[boh:bih,loh:lih])
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]
    
    
    lpc = ((hdf['l']-np.nanmean(np.unique(hdf['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((hdf['b']-np.nanmean(np.unique(hdf['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    
    delvdell = fttr(init,lpc,hdf['v'])
    delvdelb = fttr(init,bpc,hdf['v'])
    
    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    #account for similar angles, based on the signs of dx and dy
    if delvdell.slope > 0 : #first & fourth quadrants
        theta = theta
    elif delvdell.slope < 0: #second quadrant & third
        theta = theta + 180
        
    if theta > 180:
        theta -=360
        
    Omegas[i] = Omega
    thetas[i] = theta
        
        
    fig,axs = plt.subplots(1,2,figsize=(10,5),dpi=200)  
    fig.suptitle(zcat['names'][i])
    vmn,vmx = zlim(h1[boh:bih,loh:lih]*(h0[boh:bih,loh:lih]>4e20))
    imh = axs[0].imshow(h1[boh:bih,loh:lih]*(h0[boh:bih,loh:lih]>4e20),
#                vmin=vmn,vmax = vmx,
                vmin=1000,vmax=10001,
               extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],
               origin='lower',cmap=dc.cmap_min0('RdYlGn'))
    
    imh = axs[0].contourf(np.arange(5760)[loh:lih]*-0.125+360, np.arange(800)[boh:bih]*0.125-50,
               h1[boh:bih,loh:lih]*(h0[boh:bih,loh:lih]>4e20),levels=10,
#                vmin=vmn,vmax = vmx,
#                extent=[loh*-0.125+360,lih*-0.125+360,boh*0.125-50,bih*0.125-50],origin='lower',
                cmap=dc.cmap_min0('RdYlGn'))
    fig.colorbar(imh,ax=axs[0])
    axs[0].contour(np.arange(loc-3,lic+3)*-0.125+360,np.arange(boc-3,bic+3)*0.125-50,
            c1[boc-3:bic+3,loc-3:lic+3],linewidths=0.5,colors='black',levels=5)
    axs[0].contour(np.arange(loc-3,lic+3)*-0.125+360,np.arange(boc-3,bic+3)*0.125-50,
            c0[boc-3:bic+3,loc-3:lic+3],linewidths=0.5,colors='black',levels=[4e20])
    
    
#     vmn,vmx = zlim(c1[boc:bic,loc:lic]*(c0[boc:bic,loc:lic]>4e20))
    axs[1].imshow(c1[boc:bic,loc:lic]*(c0[boc:bic,loc:lic]>4e20),
#                vmin=zcat['vel_GMC<km/s>'][i][0],vmax=zcat['vel_GMC<km/s>'][i][1],
                vmin=1000,vmax=10001,
               extent=[loc*-0.125+360,lic*-0.125+360,boc*0.125-50,bic*0.125-50],
               origin='lower',cmap=dc.cmap_min0('RdYlGn'))
    
    imc = axs[1].contourf(np.arange(5760)[loc:lic]*-0.125+360, np.arange(800)[boc:bic]*0.125-50,
                c1[boc:bic,loc:lic]*(c0[boc:bic,loc:lic]>4e20),levels=10,
#                vmin=zcat['vel_GMC<km/s>'][i][0],vmax=zcat['vel_GMC<km/s>'][i][1],
#                extent=[loc*-0.125+360,lic*-0.125+360,boc*0.125-50,bic*0.125-50],origin='lower',
                 cmap='RdYlGn')
    
    
    fig.colorbar(imc,ax=axs[1],orientation='horizontal')

#     plt.scatter(((dffl['l']-np.nanmean(np.unique(dffl['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i], 
#             ((dffl['b']-np.nanmean(np.unique(dffl['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i],
#             c=dffl['v'],alpha=0.6,cmap='RdYlGn')
    
#     axs[0].plot(np.arange(5760)[loh:lih]*-0.125+360,
#          (np.arange(5760)[loh:lih]*-0.125+360)*np.sin(-90*u.deg+theta*u.deg),
#              color='black',linestyle='--') #rotation axis, 90 degree away from the velocity gradient vector 
    
    plt.show()
  

zcat['vel_grad_HI<km/s/pc>'] = Omegas
zcat['vel_grad_HI<deg>'] = thetas


#binning the data might make a difference
#propagate uncertainties -> use them as weights to calculate dvdx and dvdy
#start by fitting a gaussian to each line of sight -> get amplitude, fwhm, mean 
#sanity check -> calculate Omega and theta on the LAB data
#plot velocity vs offset from rotation axis
#look into large scale flows


zcat.write('final_cat.fits',overwrite=True)
