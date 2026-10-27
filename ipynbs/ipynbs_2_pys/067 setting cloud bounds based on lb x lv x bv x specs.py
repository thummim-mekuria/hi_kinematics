#!/usr/bin/env python
# coding: utf-8


import numpy as np
from astropy.table import QTable
from astropy.io import fits

import datacube as dc

import matplotlib.pyplot as plt


interm = QTable(fits.open('interm_cat.fits')[1].data)


zcat =  QTable(fits.open('final_cat.fits')[1].data)


hivels = np.arange(310)*1.28 - 200.
covels = np.arange(155)*1.3 - 100.

hilons = np.arange(4323)*(-0.083) + 180.166
colons = np.arange(2881)*(-0.125) + 180.125

hilats = np.arange(1200)*0.083 - 49.99
colats = np.arange(481)*0.125 - 30.125


dht = fits.open('DHT_CO_longer.fits')[0]


lvdht = np.nansum(dht.data,axis=1)


bvdht = np.nansum(dht.data,axis=2)


lbdht = np.nansum(dht.data,axis=0)


dht.data[:bc_i-bext_i:bc_i+bext_i,:].shape


zcat['lon<deg>'] = [np.array([35., 15.]),
 np.array([168., 155.]),
 np.array([152., 141.]),
 np.array([-54., -68.]),
 np.array([-132., -136.]),
 np.array([-138., -142.]),
 np.array([-165., -173.]),
 np.array([49., 41.]),
 np.array([98., 94.]),
 np.array([-139., -146.]),
 np.array([-155., -162.]),
 np.array([-143., -149.]),
 np.array([-143., -155.]),
 np.array([-151., -158.]),
 np.array([-160., -172.]),
 np.array([162., 157.]),
 np.array([126., 118.]),
 np.array([-151., -154.]),
 np.array([177., 164.])]

zcat['vel<km/s>']=[[0.,17.],[-20.,10.],[-7.,7.],[-3.,10.],[5.,22.],[7.,20.],\
[-15.,23.],[0.,20.],[-10.,10.],[20.,40.],[0.,15.],[0.,20.],\
[-3.,20.],[-3.,20.],[-5.,20.],[0.,15.],[-7.,0.],[10.,18.],[3.,15.] ]       

zcat['lat<deg>']=[np.array([ 2., 13.]),
 np.array([-12.,  -6.]),
 np.array([15., 26.]),
 np.array([-20., -12.]),
 np.array([-3.,  1.]),
 np.array([-12.,  -7.]),
 np.array([-3.,  3.]),
 np.array([ 7., 11.]),
 np.array([-17.,  -8.]),
 np.array([-4.,  0.]),
 np.array([0.,  3.]),
 np.array([-14., -11.]),
 np.array([-21., -17.]),
 np.array([-17.,  -10.]),
 np.array([-17.,  -7.]),
 np.array([-25., -17.]),
 np.array([20., 28.]),
 np.array([-3.,  -1.]),
 np.array([-18., -12.])]


cospecs = []
colbs = []
colvs = []
cobvs = []
co_cubes = []


for i in range(19):
    lc = np.nanmean(zcat['lon<deg>'][i])
    lext = np.nanstd(zcat['lon<deg>'][i])
    
    bc = np.nanmean(zcat['lat<deg>'][i])
    bext = np.nanstd(zcat['lat<deg>'][i])
    
    vc = np.nanmean(zcat['vel<km/s>'][i])
    vext = np.nanstd(zcat['vel<km/s>'][i])
    
    lc_i = int((lc-180.125)/(-0.125))
    lext_i = int(lext/0.125)
    
    bc_i = int((bc+30.125)/(0.125))
    bext_i = int(bext/0.125)
    
    vc_i = int((vc+100)/1.3)
    vext_i = int(vext/1.3)
    
    lv = np.nansum(dht.data[:,bc_i-bext_i:bc_i+bext_i,:],axis=1)
    bv = np.nansum(dht.data[:,:,lc_i-lext_i:lc_i+lext_i],axis=2)
    lb = np.nansum(dht.data[vc_i-vext_i:vc_i+vext_i,:,:],axis=0)
    cospec = np.nanmean(dht.data[:,bc_i-bext_i:bc_i+bext_i,lc_i-lext_i:lc_i+lext_i],axis=(1,2))

    cube = dht.data[vc_i-vext_i:vc_i+vext_i,bc_i-bext_i:bc_i+bext_i,lc_i-lext_i:lc_i+lext_i]
    co_cubes.append(cube)
    colbs.append(lb)
    colvs.append(lv)
    cobvs.append(bv)
    cospecs.append(cospec)
    
    fig,axs = plt.subplots(4,1,figsize=(8,10),layout='constrained')
    
    axs[0].plot(hivels,zcat['vspectrum_HI<K>'][i]/np.nanmax(zcat['vspectrum_HI<K>'][i]))
    
    
    axs[0].plot(covels,cospec/np.nanmax(cospec))
    axs[0].set_xlim(-50,70)
    axs[0].set_title(zcat['names'][i])
#     plt.show()
    
    axs[1].imshow(lv.transpose(),extent=[-100,100,180,-180],origin='lower',cmap=dc.cmap_min0('viridis'))
    axs[1].contour(covels,colons,np.sqrt(lv.transpose()),colors=['white'],\
                   linewidths=0.9,linestyles=['solid','--',],levels=[1,5,7,10])
    axs[1].set_ylim(lc+1.5*lext,lc-1.5*lext)
    axs[1].set_xlim(vc-1.5*vext,vc+1.5*vext)
    
    
    axs[2].imshow(bv.transpose(),extent=[-100,100,-30,30],origin='lower',cmap=dc.cmap_min0('viridis'))
    axs[2].contour(covels,colats,np.sqrt(bv.transpose()),colors=['white'],\
                   linewidths=0.9,linestyles=['solid','--',],levels=[1,5,7,10])
    axs[2].set_ylim(bc-1.5*bext,bc+1.5*bext)
    axs[2].set_xlim(vc-1.5*vext,vc+1.5*vext)
    
    
    axs[3].imshow(lb,extent=[180,-180,-30,30],origin='lower',cmap=dc.cmap_min0('viridis'),vmax=50)
    axs[3].contour(colons,colats,lb,colors=['white'],\
                   linewidths=0.9,linestyles=['solid','--'],levels=[2,5])
    axs[3].set_ylim(bc-1.5*bext,bc+1.5*bext)
    axs[3].set_xlim(lc+1.5*lext,lc-1.5*lext)
    
    
#     fig.tight_layout()
    plt.show()


import astropy.units as u


# zcat['lb_GMC<K.km/s>'] = colbs
# zcat['lv_GMC<K.deg>'] = colvs
# zcat['bv_GMC<K.deg>'] = cobvs


hi4pi = fits.open('HI4PI.fits')[0]


zcat.write('final_cat.fits',overwrite=True)


for i in range(19):
    lc = np.nanmean(zcat['lon<deg>'][i])
    lext = np.nanstd(zcat['lon<deg>'][i])
    
    bc = np.nanmean(zcat['lat<deg>'][i])
    bext = np.nanstd(zcat['lat<deg>'][i])
    
    vc = np.nanmean(zcat['vel<km/s>'][i])
    vext = np.nanstd(zcat['vel<km/s>'][i])
    
    
    lc_i = int((lc-180.166)/(-0.083))
    lext_i = int(lext/0.083)
    
    bc_i = int((bc+49.999)/(0.083))
    bext_i = int(bext/0.083)
    
    vc_i = int((vc+200)/1.288)
    vext_i = int(vext/1.288)
    
    
    lb = np.nansum(hi4pi.data[vc_i-vext_i:vc_i+vext_i,:,:],axis=0)
    lv = np.nansum(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,:],axis=1)
    bv = np.nansum(hi4pi.data[:,:,lc_i-lext_i:lc_i+lext_i],axis=2)
    spec = np.nanmean(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,lc_i-lext_i:lc_i+lext_i],axis=(1,2))

    
    fig,axs = plt.subplots(4,1,figsize=(8,10),layout='constrained')
    
    axs[0].set_title(zcat['names'][i])
    axs[0].imshow(lb,extent=[180,-180,-50,50],origin='lower',cmap=dc.cmap_min0('viridis'),vmax=600)
    axs[0].contour(hilons,hilats,lb,colors=['red'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[50,100,200,400,800])
#     axs[0].contour(colons,colats,np.nansum(cocubes[i],axis=0),colors=['red'],\
#                    linewidths=0.4,linestyles=['solid','--'],levels=[2,5])
    axs[0].set_ylim(bc-5*bext,bc+5*bext)
    axs[0].set_xlim(lc+5*lext,lc-5*lext)
    
    
    axs[1].imshow(np.sqrt(lv.transpose()),extent=[-200,200,180,-180],origin='lower',\
                  cmap=dc.cmap_min0('viridis'))
    axs[1].contour(hivels,hilons,lv.transpose(),colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
#     axs[1].contour(covels,colons,np.nansum(cocubes[i],axis=1).transpose(),colors=['red'],\
#                    linewidths=0.4,linestyles=['solid','--'],levels=[1,5,7,10])
    axs[1].set_xlim(vc-5*vext,vc+5*vext)
    axs[1].set_ylim(lc+5*lext,lc-5*lext)
    
    
    axs[2].imshow(np.sqrt(bv.transpose()),extent=[-200,200,-50,50],origin='lower',\
                  cmap=dc.cmap_min0('viridis'))
    axs[2].contour(hivels,hilats,bv.transpose(),colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
#     axs[2].contour(covels,colats,np.nansum(cocubes[i],axis=2).transpose(),colors=['red'],\
#                    linewidths=0.4,linestyles=['solid','--'],levels=[1,5,7,10])
    axs[2].set_xlim(vc-5*vext,vc+5*vext)
    axs[2].set_ylim(bc+5*bext,bc-5*bext)
    
    axs[3].plot(hivels,spec/np.nanmax(spec))
#     axs[3].plot(covels,np.nanmean(cocubes[i],axis=(1,2))/np.nanmax(np.nanmean(cocubes[i],axis=(1,2))))
    axs[3].set_xlim(vc-5*vext,vc+5*vext)
    
    
    plt.show()
    

import astropy.constants as c


400*(u.K*u.km/u.s * 1.288*1.36*1.67e-27*u.kg*1.828e18*u.cm**-2/(u.K*u.km/u.s)).to(c.M_sun/u.pc**2)


2*(u.K*u.km/u.s * 1.3*2*1.36*1.67e-27*u.kg*1.9e20*u.cm**-2/(u.K*u.km/u.s)).to(c.M_sun/u.pc**2)


for i in range(19):
    lc = np.nanmean(zcat['lon<deg>'][i])
    lext = np.nanstd(zcat['lon<deg>'][i])
    
    bc = np.nanmean(zcat['lat<deg>'][i])
    bext = np.nanstd(zcat['lat<deg>'][i])
    
    vc = np.nanmean(zcat['vel<km/s>'][i])
    vext = np.nanstd(zcat['vel<km/s>'][i])
    
    
    lc_i = int((lc-180.166)/(-0.083))
    lext_i = int(lext/0.083)
    
    bc_i = int((bc+49.999)/(0.083))
    bext_i = int(bext/0.083)
    
    vc_i = int((vc+200)/1.288)
    vext_i = int(vext/1.288)
    
    
    lb = np.nansum(hi4pi.data[vc_i-vext_i:vc_i+vext_i,:,:],axis=0)
    lv = np.nansum(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,:],axis=1)
    bv = np.nansum(hi4pi.data[:,:,lc_i-lext_i:lc_i+lext_i],axis=2)
    spec = np.nanmean(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,lc_i-lext_i:lc_i+lext_i],axis=(1,2))

    
    fig,axs = plt.subplots(2,1,figsize=(8,10),layout='constrained')
    
    axs[0].set_title(zcat['names'][i])
    axs[0].imshow(lb,extent=[180,-180,-50,50],origin='lower',cmap=dc.cmap_min0('viridis'),vmax=1200)
    axs[0].contour(hilons,hilats,lb,colors=['red'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[400,800])
    axs[0].contour(colons,colats,colbs[i],colors=['white'],\
                   linewidths=0.4,linestyles=['solid','--'],levels=[2,5])
    axs[0].set_ylim(bc-5*bext,bc+5*bext)
    axs[0].set_xlim(lc+5*lext,lc-5*lext)
    
    
#     axs[1].imshow(np.sqrt(lv.transpose()),extent=[-200,200,180,-180],origin='lower',\
#                   cmap=dc.cmap_min0('viridis'))
#     axs[1].contour(hivels,hilons,lv.transpose(),colors=['white'],\
#                    linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
# #     axs[1].contour(covels,colons,np.nansum(cocubes[i],axis=1).transpose(),colors=['red'],\
# #                    linewidths=0.4,linestyles=['solid','--'],levels=[1,5,7,10])
#     axs[1].set_xlim(vc-5*vext,vc+5*vext)
#     axs[1].set_ylim(lc+5*lext,lc-5*lext)
    
    
#     axs[2].imshow(np.sqrt(bv.transpose()),extent=[-200,200,-50,50],origin='lower',\
#                   cmap=dc.cmap_min0('viridis'))
#     axs[2].contour(hivels,hilats,bv.transpose(),colors=['white'],\
#                    linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
# #     axs[2].contour(covels,colats,np.nansum(cocubes[i],axis=2).transpose(),colors=['red'],\
# #                    linewidths=0.4,linestyles=['solid','--'],levels=[1,5,7,10])
#     axs[2].set_xlim(vc-5*vext,vc+5*vext)
#     axs[2].set_ylim(bc+5*bext,bc-5*bext)
    
#     axs[3].plot(hivels,spec/np.nanmax(spec))
# #     axs[3].plot(covels,np.nanmean(cocubes[i],axis=(1,2))/np.nanmax(np.nanmean(cocubes[i],axis=(1,2))))
#     axs[3].set_xlim(vc-5*vext,vc+5*vext)
    
    
    plt.show()
    

from sklearn.decomposition import PCA


hmm = np.nanmean(np.array((zcat['lon<deg>'],zcat['lat<deg>'],zcat['vel<km/s>'])),axis=1)


hmm_0 = np.nanmean(np.array((zcat['lon<deg>'],zcat['lat<deg>'],zcat['vel<km/s>'])),axis=2)


pca = PCA().fit(hmm_0)


pca.components_


#come back and figure out how to plot this? / see how they are clustered
#basically what i wanna know is which clouds are close enough spatially and kinematically that it makes
#more sense to look for a common HI for them rather than a separate parent for each
#also there is a class of these that are: off on their own (chamaeleon,lacerta maybe) and /or the HI near
#them is <10 Msun/pc^2 (check how low you can (if you can) go to draw a closedish contour (acros lb,lv,bv))
#also, look at the hi spectrum with just an lb cut, if the hi peak is offset, sigma is being underestimated bc the
#integration window i used was for the co
#ad then there are the ones @ b~0 (gemini, maddalena) where the closed contour around them is >=20Mun/pc^2


for i in [2,3,7,8,15,16,18]:
    print(zcat['names'][i])


cs_and_exts_q


cs_and_exts_q = np.array([(np.nanmean(zcat[column][i]),np.nanstd(zcat[column][i])) \
                         for column in ['lon<deg>','lat<deg>','vel<km/s>']\
                         for i in [2,3,7,8,15,16,18]]).reshape(3,-1,2).transpose().reshape()
# cs = np.array([cs_and_exts_q[i].transpose()[0] for i in range(len(cs_and_exts_q))])

# l0s,b0s,v0s = inits.transpose()


for i in [2,3,7,8,15,16,18]:
    lc = np.nanmean(zcat['lon<deg>'][i])
    lext = np.nanstd(zcat['lon<deg>'][i])
    
    bc = np.nanmean(zcat['lat<deg>'][i])
    bext = np.nanstd(zcat['lat<deg>'][i])
    
    vc = np.nanmean(zcat['vel<km/s>'][i])
    vext = np.nanstd(zcat['vel<km/s>'][i])
    
    
    lc_i = int((lc-180.166)/(-0.083))
    lext_i = int(lext/0.083)
    
    bc_i = int((bc+49.999)/(0.083))
    bext_i = int(bext/0.083)
    
    vc_i = int((vc+200)/1.288)
    vext_i = int(vext/1.288)
    
    
    lb = np.nansum(hi4pi.data[vc_i-vext_i:vc_i+vext_i,:,:],axis=0)
    lv = np.nansum(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,:],axis=1)
    bv = np.nansum(hi4pi.data[:,:,lc_i-lext_i:lc_i+lext_i],axis=2)
    spec = np.nanmean(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,lc_i-lext_i:lc_i+lext_i],axis=(1,2))

    
    thiscolons = np.arange(co_cubes[i].shape[2])*(-0.125) + (lc+lext)
    thiscolats = np.arange(co_cubes[i].shape[1])*0.125 + (bc-bext)
    thiscovels = np.arange(co_cubes[i].shape[0])*0.125 + (vc-vext)
    
#     print(zcat['lon<deg>'][i],thiscolons)
#     print(zcat['lat<deg>'][i],thiscolats)
#     print(zcat['vel<km/s>'][i],thiscovels)
    
    fig,axs = plt.subplots(4,1,figsize=(8,10),layout='constrained')
    axs[0].set_title(zcat['names'][i])
    axs[0].imshow(lb,extent=[180,-180,-50,50],origin='lower',cmap=dc.cmap_min0('viridis'),vmax=1200)
    axs[0].contour(hilons,hilats,lb,colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[100,200,400])
    axs[0].contour(thiscolons,thiscolats,np.nansum(co_cubes[i],axis=0),colors=['purple'],\
                   linewidths=0.4,linestyles=['solid','--'],levels=[0.01,0.05,0.1,0.3])
    axs[0].set_ylim(bc-5*bext,bc+5*bext)
    axs[0].set_xlim(lc+5*lext,lc-5*lext)
    
    
    axs[1].imshow(np.sqrt(lv.transpose()),extent=[-200,200,180,-180],origin='lower',\
                  cmap=dc.cmap_min0('viridis'))
    axs[1].contour(hivels,hilons,lv.transpose(),colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
    axs[1].contour(thiscovels,thiscolons,np.nansum(co_cubes[i],axis=1).transpose(),colors=['purple'],\
                   linewidths=0.4,linestyles=['solid','--'],levels=[0.01,0.05,0.1,0.3])
    axs[1].set_xlim(vc-5*vext,vc+5*vext)
    axs[1].set_ylim(lc+5*lext,lc-5*lext)
    
    
    axs[2].imshow(np.sqrt(bv.transpose()),extent=[-200,200,-50,50],origin='lower',\
                  cmap=dc.cmap_min0('viridis'))
    axs[2].contour(hivels,hilats,bv.transpose(),colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
    axs[2].contour(thiscovels,thiscolats,np.nansum(co_cubes[i],axis=2).transpose(),colors=['purple'],\
                   linewidths=0.4,linestyles=['solid','--'],levels=[0.01,0.05,0.1,0.3])
    axs[2].set_xlim(vc-5*vext,vc+5*vext)
    axs[2].set_ylim(bc+5*bext,bc-5*bext)
    
    axs[3].plot(hivels,spec/np.nanmax(spec))
    axs[3].plot(thiscovels,np.nanmean(co_cubes[i],axis=(1,2))/np.nanmax(np.nanmean(co_cubes[i],axis=(1,2))))
    axs[3].set_xlim(vc-5*vext,vc+5*vext)
    
    
    plt.show()
    

for i in [2,3,7,8,15,16,18]:
    lc = np.nanmean(zcat['lon<deg>'][i])
    lext = np.nanstd(zcat['lon<deg>'][i])
    
    bc = np.nanmean(zcat['lat<deg>'][i])
    bext = np.nanstd(zcat['lat<deg>'][i])
    
    vc = np.nanmean(zcat['vel<km/s>'][i])
    vext = np.nanstd(zcat['vel<km/s>'][i])
    
    
    lc_i = int((lc-180.166)/(-0.083))
    lext_i = int(lext/0.083)
    
    bc_i = int((bc+49.999)/(0.083))
    bext_i = int(bext/0.083)
    
    vc_i = int((vc+200)/1.288)
    vext_i = int(vext/1.288)
    
    
    lb = np.nansum(hi4pi.data[vc_i-vext_i:vc_i+vext_i,:,:],axis=0)
    lv = np.nansum(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,:],axis=1)
    bv = np.nansum(hi4pi.data[:,:,lc_i-lext_i:lc_i+lext_i],axis=2)
    spec = np.nanmean(hi4pi.data[:,bc_i-bext_i:bc_i+bext_i,lc_i-lext_i:lc_i+lext_i],axis=(1,2))

    
    thiscolons = np.arange(co_cubes[i].shape[2])*(-0.125) + (lc+lext)
    thiscolats = np.arange(co_cubes[i].shape[1])*0.125 + (bc-bext)
    thiscovels = np.arange(co_cubes[i].shape[0])*1.3 + (vc-vext)
    
#     print(zcat['lon<deg>'][i],thiscolons)
#     print(zcat['lat<deg>'][i],thiscolats)
#     print(zcat['vel<km/s>'][i],thiscovels)
    
    fig,axs = plt.subplots(4,1,figsize=(8,10),layout='constrained')
    axs[0].set_title(zcat['names'][i])
    axs[0].imshow(lb,extent=[180,-180,-50,50],origin='lower',cmap=dc.cmap_min0('viridis'),vmax=1200)
    axs[0].contour(hilons,hilats,lb,colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[100,200,400])
    axs[0].contour(thiscolons,thiscolats,np.nansum(co_cubes[i],axis=0),colors=['purple'],\
                   linewidths=0.4,linestyles=['solid','--'])
    axs[0].set_ylim(bc-3.5*bext,bc+3.5*bext)
    axs[0].set_xlim(lc+3.5*lext,lc-3.5*lext)
    
    
    axs[1].imshow(np.sqrt(lv.transpose()),extent=[-200,200,180,-180],origin='lower',\
                  cmap=dc.cmap_min0('viridis'))
    axs[1].contour(hivels,hilons,lv.transpose(),colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
    axs[1].contour(thiscovels,thiscolons,np.nansum(co_cubes[i],axis=1).transpose(),colors=['purple'],\
                   linewidths=0.4,linestyles=['solid','--'])
    axs[1].set_xlim(vc-3.5*vext,vc+3.5*vext)
    axs[1].set_ylim(lc+3.5*lext,lc-3.5*lext)
    
    
    axs[2].imshow(np.sqrt(bv.transpose()),extent=[-200,200,-50,50],origin='lower',\
                  cmap=dc.cmap_min0('viridis'))
    axs[2].contour(hivels,hilats,bv.transpose(),colors=['white'],\
                   linewidths=0.6,linestyles=['solid','--'],levels=[1e3,2e3,3e3,4e3,6e3,8e3])
    axs[2].contour(thiscovels,thiscolats,np.nansum(co_cubes[i],axis=2).transpose(),colors=['purple'],\
                   linewidths=0.4,linestyles=['solid','--'])
    axs[2].set_xlim(vc-3.5*vext,vc+3.5*vext)
    axs[2].set_ylim(bc+3.5*bext,bc-3.5*bext)
    
    axs[3].plot(hivels,spec/np.nanmax(spec))
    axs[3].plot(thiscovels,np.nanmean(co_cubes[i],axis=(1,2))/np.nanmax(np.nanmean(co_cubes[i],axis=(1,2))))
#     axs[3].plot(covels,cospecs[i]/np.nanmax(cospecs[i]))

    axs[3].set_xlim(vc-3.5*vext,vc+3.5*vext)
    
    
    plt.show()
    

# sort into groups based on lbv proximity / hi overlap possibility:
# 
# group 1 = aqr
# group 2 = california
# group 3 = camelopardalis
# group 4 = chamaeleon
# group 5 = cma ob1
# group 6 = crossbones + oriona+b+lam+monr2
# group 7 = gemini
# group 8 = hercules
# group 8 = lacerta
# group 9 = maddalena
# group 10 = mon ob1
# group 11 = perseus + taurus
# group 12 = polaris
# group 13 = rosette


hi_table = QTable()
hi_table['id'] = np.arange(14)


hi_table['lon<deg>'] = np.array(([35.,15.],[168.,155.],[152.,141.],[-54.,-68.],\
                                 [-132.,-136.],[-138.,-172.],[-165.,-173.],[49.,41.],[98.,94.],\
                                 [-139.,-146.],[-155.,-162.],[177.,157.],[126.,118.],[-151.,-154.]))

hi_table['lat<deg>'] = np.array(([2.,13.],[-12.,-6.],[15.,26.],[-20.,-12.],\
                                 [-3.,1.],[-21.,-7.],[-3.,3.],[7.,11.],[-17.,-8.],\
                                 [-4.,0.],[0.,3.],[-25.,-12.],[20.,28.],[-3.,-1.]))

hi_table['vel<km/s>'] = np.array(([0.,17.],[-20.,10.],[-7.,7.],[-3.,10.],\
                                 [5.,22.],[-3.,20.],[-15.,23.],[0.,20.],[-10.,10.],\
                                 [20.,40.],[0.,15.],[0.,15.],[-7.,0.],[10.,18.]))


hi_table['names'] = ['aquarius rift','california','camelopardalis','chamaeleon','cma ob1',\
                    'crossbones, orion complex', 'gemini ob1', 'hercules','lacerta','maddalena','mon ob1',\
                    'perseus, taurus','polaris','rosette']


hi_table.write('hi_table.fits',overwrite=True)


for i in range(14):
    lc = np.nanmean(hi_table['lon<deg>'][i])
    lext = np.nanstd(hi_table['lon<deg>'][i])
    
    bc = np.nanmean(hi_table['lat<deg>'][i])
    bext = np.nanstd(hi_table['lat<deg>'][i])
    
    vc = np.nanmean(hi_table['vel<km/s>'][i])
    vext = np.nanstd(hi_table['vel<km/s>'][i])
    
    
    lc_co = int((lc-180.125)/(-0.125))
    lext_co = int(lext/0.125)
    
    bc_co = int((bc+30.125)/(0.125))
    bext_co = int(bext/0.125)
    
    vc_co = int((vc+100)/1.3)
    vext_co = int(vext/1.3)
    
    
    lc_hi = int((lc-180.166)/(-0.083))
    lext_hi = int(lext/0.083)
    
    bc_hi = int((bc+49.999)/(0.083))
    bext_hi = int(bext/0.083)
    
    vc_hi = int((vc+200)/1.288)
    vext_hi = int(vext/1.288)
    
    
    lb_co = np.nansum(dht.data[vc_co-vext_co:vc_co+vext_co,:,:],axis=0)
    lv_co = np.nansum(dht.data[:,bc_co-bext_co:bc_co+bext_co,:],axis=1)
    bv_co = np.nansum(dht.data[:,:,lc_co-lext_co:lc_co+lext_co],axis=2)
#     spec = np.nanmean(dht.data[:,bc_co-bext_co:bc_i+bext_co,lc_co-lext_co:lc_co+lext_co],axis=(1,2))

    lb_hi = np.nansum(hi4pi.data[vc_hi-vext_hi:vc_hi+vext_hi,:,:],axis=0)
    lv_hi = np.nansum(hi4pi.data[:,bc_hi-bext_hi:bc_hi+bext_hi,:],axis=1)
    bv_hi = np.nansum(hi4pi.data[:,:,lc_hi-lext_hi:lc_hi+lext_hi],axis=2)

    
    fig,axs=plt.subplots(2,2,figsize=(12,12),layout='constrained')
    
    axs[0,0].imshow(lb_hi,extent=[180,-180,-50,50],origin='lower',\
                  cmap=dc.cmap_min0('viridis'),vmax = 1200)
    axs[0,0].contour(hilons,hilats,lb_hi,levels=[200,400,600,800],colors=['white','pink'],alpha=0.7)

    axs[0,0].contour(colons,colats,lb_co,levels=[0.1,2,5],colors=['purple','green','blue'],\
                   alpha=0.7,linewidths=[0.5])
    axs[0,0].set_ylim(bc-2*bext,bc+2*bext)
    axs[0,0].set_xlim(lc+2*lext,lc-2*lext)
#     print(0.5*np.nanmax(lb_hi))
    
    
    axs[1,0].imshow(lv_hi,origin='lower',cmap=dc.cmap_min0('viridis'),\
                  vmax = 0.5*np.nanmax(lv_hi),extent=[180,-180,-200,200])
    axs[1,0].contour(colons,covels,lv_co,levels=[0.1,2,5],\
                   colors=['purple','green','blue'],alpha=0.7,linewidths=[0.5])
    axs[1,0].set_ylim(vc-2*vext,vc+2*vext)
    axs[1,0].set_xlim(lc+2*lext,lc-2*lext)
    
    
    axs[0,1].imshow(bv_hi,origin='lower',cmap=dc.cmap_min0('viridis'),\
                  vmax = 0.5*np.nanmax(bv_hi),extent=[180,-180,-200,200])
    axs[0,1].contour(colats,covels,bv_co,levels=[0.1,2,5],\
                   colors=['purple','green','blue'],alpha=0.7,linewidths=[0.5])
    axs[0,1].set_ylim(vc-2*vext,vc+2*vext)
    axs[0,1].set_xlim(bc+2*bext,bc-2*bext)
    
    
    axs[0,0].set_title(hi_table['names'][i])
    
    plt.show()


i = [0,7]
zcat['names'][i],zcat['vel<km/s>'][i]


zcat[['names','lon<deg>','lat<deg>','vel<km/s>']]
