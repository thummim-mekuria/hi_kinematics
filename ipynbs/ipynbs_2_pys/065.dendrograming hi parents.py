#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import datacube as dc
import clouds as cld

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import matplotlib.pyplot as plt


zcat = QTable(fits.open('final_cat.fits')[1].data)


zcat.write('final_cat.fits',overwrite=True)


zcat['names']=['aquila_r',
 'california',
 'camelopardis',
 'cepheus_f',
 'cepheus_n',
 'chamaeleon',
 'cma_ob1',
 'crossbones',
 'gemini_ob1',
 'hercules',
 'lacerta',
 'maddalena',
 'mon_ob1',
 'mon_r2',
 'orion_a',
 'orion_b',
 'orion_l',
 'perseus',
 'polaris',
 'rosette',
 'taurus']


hi4pi_r = fits.open('HI4PI.fits')[0]


dht = fits.open('DHT_CO_longer.fits')[0]


# for i in range(21):
#     lo,li = ((zcat['lon<deg>'][i] - 180.125)/(-0.125)).astype('int')
#     bo,bi = ((zcat['lat<deg>'][i] + 30.125)/(0.125)).astype('int')
    
#     sb = dht.data[:,bo:bi,lo:li]
#     sig,ar = cld.surface_density(sb,'CO')
    
#     ang_ar = (u.rad**2).to(u.deg**2) *ar*len(np.where(sig>0)[0])
#     rgad = np.sqrt(ang_ar/np.pi)
#     pc2_ar = ar*len(np.where(sig>0)[0])*(zcat['distance<pc>'][i]*u.pc)**2
#     rg = np.sqrt(pc2_ar/np.pi)
# #     zcat['radius_GMC<deg>'][i] = rgad
#     zcat['radius_GMC<pc>'][i] = rg.value
    
# #     m = cld.mass(sig,zcat['distance<pc>'][i]*u.pc,ar)
# #     ra,s = cld.accumulation_radius(m,zcat['distance<pc>'][i]*u.pc,u.deg*np.nanmean(zcat['lat<deg>'][i]))
    
# #     zcat['mass_GMC<kg>'][i] = m.value
# #     zcat['radius_HI<pc>'][i] = ra.value
# #     zcat['radius_HI<deg>'][i] = s.value
    
    
# #     spc = np.nanmean(sb,axis=(1,2))
    
# #     g = gft(spc)
# #     zcat['vdispersion_GMC<km/s>'][i] = 1.3*g.stddev.value
    
# # #     plt.plot(spc)
# # #     plt.plot(g(np.arange(155)))
# # #     plt.show()
    
# # #     zcat['vspectrum_GMC<K>'][i] = spc
    
    
# # #     w = np.nansum(sb,axis=0)
# #     plt.imshow(sig.value,origin='lower',cmap=dc.cmap_min0('viridis'))
# #     plt.show()


# for i in range(21):
#     hicube = hi4pi_r.data[:,:,1:-1]
#     lc,bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
#     rhi = zcat['radius_HI<deg>'][i]

#     lc = round((lc - 180.1666)/(-0.08333))
#     bc = round((bc + 49.999)/(0.08333))
#     rhi = round(rhi/0.0833)
    
    
#     l0,l1 = int(lc-rhi),int(lc+rhi)
#     b0,b1 = int(bc-rhi),int(bc+rhi)
    
#     if l0 < 0:
#         sb = np.concatenate((hicube[:,b0:b1,l0:],hicube[:,b0:b1,:l1]),axis=2)
        
#     elif l1 > hicube.shape[2]:
#         l1_n = hicube.shape[2] - l1
#         sb = np.concatenate((hicube[:,b0:b1,l0:],hicube[:,b0:b1,:l1_n]),axis=2)
        
#     else:
#         sb = hicube[:,b0:b1,l0:l1]
        
       
# #     sig,ar = cld.surface_density(sb,'21cm')
        
# # #     w = np.nansum(sb,axis=0)
# # #     w1 = np.nansum(sb,axis=1)
# # #     plt.imshow(w1,origin='lower',cmap=dc.cmap_min0('viridis'))
# # #     plt.title(zcat['names'][i])
# # #     plt.show()
    
#     hivels = np.arange(310)*1.288 - 200.
#     covels = np.arange(155)*1.3 - 100.
    
#     spc = np.nanmean(sb,axis=(1,2))
#     plt.plot(hivels,spc)
#     plt.plot(covels,50*zcat['vspectrum_GMC<K>'][i])
# #     plt.vlines(covels[np.argmax(zcat['vspectrum_GMC<K>'][i])]\
# #                - 5*zcat['vdispersion_GMC<km/s>'][i],0,np.nanmax(spc))
# #     plt.vlines(covels[np.argmax(zcat['vspectrum_GMC<K>'][i])]\
# #                + 5*zcat['vdispersion_GMC<km/s>'][i],0,np.nanmax(spc))
# #     plt.title(zcat['names'][i])
#     plt.xlim(-80,80)
#     plt.show()
#     pxs = sb.transpose().reshape(-1,310)
#     vps = np.array([hivels[np.nanargmax(pxs[i])] for i in range(len(pxs))])
#     vps_sq = (vps.reshape(sb.transpose().shape[:2])).transpose()
    
#     plt.imshow(vps_sq,origin='lower')
#     plt.title(zcat['names'][i])
#     plt.colorbar()
#     plt.show()
    

cmplq = zcat[[6, 11, 19]]


def get_vps(cube,vels):
     
    pxs = cube.transpose().reshape(-1,cube.shape[0])
    vps = np.array([vels[np.nanargmax(pxs[i])] for i in range(len(pxs))])
    vps_sq = (vps.reshape(cube.transpose().shape[:2])).transpose()
 
    return vps_sq


from astrodendro import Dendrogram


cmplq_dnds = []


cmplq_dnds_flt = []


hi4pi_f = fits.open('HI4PI_flatcalib_80.fits')[0]


for i in range(3):
    hicube = hi4pi_r.data[:,:,1:-1]
    lc,bc = np.nanmean(cmplq['lon<deg>'][i]), np.nanmean(cmplq['lat<deg>'][i])
    rhi = cmplq['radius_HI<deg>'][i]
    rco = cmplq['radius_GMC<deg>'][i]

    lchi = round((lc - 180.1666)/(-0.08333))
    bchi = round((bc + 49.999)/(0.08333))
    rhi = round(rhi/0.0833)
    
    lcco = round((lc - 180.125)/(-0.125))
    bcco = round((bc + 30.125)/(0.125))
    rco = round(rco/0.125)
    
    hivels = np.arange(310)*1.288 - 200.
    covels = np.arange(155)*1.3 - 100.
   
    sbhi = hicube[:,bchi-rhi:bchi+rhi,lchi-rhi:lchi+rhi]
    sbhi_flt = hi4pi_f.data[:,bchi-rhi:bchi+rhi,lchi-rhi:lchi+rhi]

    sbco = dht.data[:,bcco-rco:bcco+rco,lcco-rco:lcco+rco]

#     vps_sq_hi = get_vps(sbhi,hivels)
#     vps_sq_co = get_vps(sbco,covels)
    
    dnd = Dendrogram.compute(sbhi,min_value=np.nanstd(sbhi),\
                             min_delta=np.nanstd(sbhi),min_npix=1e-3*sbhi.size)
#     print(time.time())
    dnd2 = Dendrogram.compute(sbhi_flt,min_value=np.nanstd(sbhi),\
                              min_delta=np.nanstd(sbhi),min_npix=1e-3*sbhi.size)
#     print(time.time())
    cmplq_dnds.append(dnd)
    cmplq_dnds_flt.append(dnd2)
    
#     plt.imshow(vps_sq_hi,origin='lower', cmap=dc.cmap_min0('Paired'),vmin=-30,vmax=100)

#     plt.title(cmplq['names'][i])
#     plt.colorbar()
#     plt.show()
    
    
#     plt.imshow(vps_sq_co,origin='lower',cmap=dc.cmap_min0('Paired'),vmin=-30,vmax=100)
#     plt.title(cmplq['names'][i])
#     plt.colorbar()
#     plt.show()
    
#     spc = np.nanmean(sbhi,axis=(1,2))
#     plt.plot(hivels,spc)
#     plt.plot(covels,50*cmplq['vspectrum_GMC<K>'][i])
#     plt.xlim(-30,100)
#     plt.show()
    
#     print('\n\n\n')
    

cma_d,mad_d,ros_d = cmplq_dnds[3:]


cma_d_f,mad_d_f,ros_d_f = cmplq_dnds_flt[3:]


for i in cma_d:
    if i.is_leaf:
        plt.plot(hivels,np.nanmean(i.get_mask(),\
                                   axis=(1,2))/np.nanmax(np.nanmean(i.get_mask(),axis=(1,2))),\
                label=str(i.idx))

plt.plot(covels,cmplq['vspectrum_GMC<K>'][0]/np.nanmax(cmplq['vspectrum_GMC<K>'][0]))
plt.xlim(-30,100)
plt.legend()
plt.show()
#     plt.imshow(np.nansum(i.get_mask(),axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))
#     plt.show()
#     plt.imshow(np.nansum(i.get_mask(),axis=1),origin='lower',cmap=dc.cmap_min0('viridis'))
#     plt.show()
#     print('\n\n')


for i in cma_d_f:
    if i.is_leaf:
        plt.plot(hivels,np.nanmean(i.get_mask(),\
                                   axis=(1,2))/np.nanmax(np.nanmean(i.get_mask(),axis=(1,2))),\
                label=str(i.idx))

plt.plot(covels,cmplq['vspectrum_GMC<K>'][0]/np.nanmax(cmplq['vspectrum_GMC<K>'][0]))
plt.xlim(-30,100)
plt.legend()
plt.show()
#     plt.imshow(np.nansum(i.get_mask(),axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))
#     plt.show()
#     plt.imshow(np.nansum(i.get_mask(),axis=1),origin='lower',cmap=dc.cmap_min0('viridis'))
#     plt.show()
#     print('\n\n')


plt.imshow(np.nansum(cma_d[2].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()
plt.imshow(np.nansum(cma_d_f[4].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()


for i in mad_d:
    if i.is_leaf:
        plt.plot(hivels,np.nanmean(i.get_mask(),\
                                   axis=(1,2))/np.nanmax(np.nanmean(i.get_mask(),axis=(1,2))),\
                label=str(i.idx))

plt.plot(covels,cmplq['vspectrum_GMC<K>'][1]/np.nanmax(cmplq['vspectrum_GMC<K>'][1]))
plt.xlim(-30,100)
plt.legend()
plt.show()


for i in mad_d_f:
    if i.is_leaf:
        plt.plot(hivels,np.nanmean(i.get_mask(),\
                                   axis=(1,2))/np.nanmax(np.nanmean(i.get_mask(),axis=(1,2))),\
                label=str(i.idx))

plt.plot(covels,cmplq['vspectrum_GMC<K>'][1]/np.nanmax(cmplq['vspectrum_GMC<K>'][1]))
plt.xlim(-30,100)
plt.legend()
plt.show()


plt.imshow(np.nansum(mad_d[3].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()
plt.imshow(np.nansum(mad_d_f[3].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()
#     print('\n\n')


for i in ros_d:
    if i.is_leaf:
        plt.plot(hivels,np.nanmean(i.get_mask(),\
                                   axis=(1,2))/np.nanmax(np.nanmean(i.get_mask(),axis=(1,2))),\
                label=str(i.idx))

plt.plot(covels,cmplq['vspectrum_GMC<K>'][2]/np.nanmax(cmplq['vspectrum_GMC<K>'][2]))
plt.xlim(-30,100)
plt.legend()
plt.show()


for i in ros_d_f:
    if i.is_leaf:
        plt.plot(hivels,np.nanmean(i.get_mask(),\
                                   axis=(1,2))/np.nanmax(np.nanmean(i.get_mask(),axis=(1,2))),\
                label=str(i.idx))

plt.plot(covels,cmplq['vspectrum_GMC<K>'][2]/np.nanmax(cmplq['vspectrum_GMC<K>'][2]))
plt.xlim(-30,100)
plt.legend()
plt.show()


np.concatenate((cmplq['lon<deg>'][0],cmplq['lat<deg>'][0]))


plt.imshow(np.nansum(ros_d[1].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()
plt.imshow(np.nansum(ros_d_f[2].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()
plt.imshow(np.nansum(ros_d_f[4].get_mask(),axis=0),origin='lower',\
           cmap=dc.cmap_min0('viridis'))
plt.show()
#     print('\n\n')


sbs = []
sbs_f = []


all_dnds = []
all_dnds_f = []


import bootleg as btlg


for i in range(21):
    btlg.dend_to_fits(all_dnds[i],zcat['names'][i]+'_dend.fits')
    btlg.dend_to_fits(all_dnds_f[i],zcat['names'][i]+'_dend_flat.fits')


diff_nums = [len(all_dnds[i]) != len(all_dnds_f[i]) for i in range(21)]


zcat['names'][diff_nums]


for i in range(21):
    if diff_nums[i]:
        print(zcat['names'][i])
        print('raws')     
        for tr in all_dnds[i]:
            if tr.is_leaf:
                plot_struct(tr)
            
        print('flats')
        for tr in all_dnds_f[i]:
            if tr.is_leaf:
                plot_struct(tr)
            
        print('\n\n\n')


len(np.nanmean(sbs[i],axis=(1,2))),len(hivels)


for i in range(21):
    plt.plot(hivels,np.nanmean(sbs[i],axis=(1,2))/np.nanmax(np.nanmean(sbs[i],axis=(1,2))),\
             color='black',linestyle='--')
    plt.plot(covels,zcat['vspectrum_GMC<K>'][i]/np.nanmax(zcat['vspectrum_GMC<K>'][i]),color='black')
    for lf in all_dnds_f[i]:
        if lf.is_leaf:
            this = np.nanmean(lf.get_mask()*all_dnds[i].data,axis=(1,2))
            plt.plot(hivels,this/np.nanmax(this),label=lf.idx)
    plt.legend()
    plt.title(zcat['names'][i])
    plt.xlim(-20,20)
    plt.show()


for i in range(21):
    plt.plot(hivels,np.nanmean(sbs[i],axis=(1,2))/np.nanmax(np.nanmean(sbs[i],axis=(1,2))),\
             color='black',linestyle='--')
    plt.plot(covels,zcat['vspectrum_GMC<K>'][i]/np.nanmax(zcat['vspectrum_GMC<K>'][i]),color='black')
    for lf in all_dnds[i]:
        if lf.is_leaf:
            this = np.nanmean(lf.get_mask()*all_dnds[i].data,axis=(1,2))
            plt.plot(hivels,this/np.nanmax(this),label=lf.idx)
    plt.legend()
    plt.title(zcat['names'][i])
    plt.xlim(-20,20)
    plt.show()


all_dnds[5].trunk


[len(dnd) for dnd in all_dnds_f]


plot_struct(all_dnds_f[5].trunk[0])


plot_struct(all_dnds[5].trunk[0])


plot_struct(all_dnds[5].trunk[1])


sbsco = []
for i in range(21):
    l0,l1 = ((zcat['lon<deg>'][i] - 180.125)/(-0.125)).astype('int')
    b0,b1 = ((zcat['lat<deg>'][i] + 30.125)/0.125).astype('int')
     
    sbco = dht.data[:,b0:b1,l0:l1]
    sbsco.append(sbco)


def plot_struct(structure,kind=0):
    m = structure.get_mask()
    w = np.nansum(m,axis=kind)
    
    plt.imshow(w,origin='lower',cmap=dc.cmap_min0('viridis'))
    plt.colorbar()
    plt.show()


def plot_hi_w_cocontour(i,hdend):
    l0,l1 = zcat['lon<deg>'][i]
    b0,b1 = zcat['lat<deg>'][i]
    rco = zcat['radius_GMC<deg>'][i]
    rhi = zcat['radius_HI<deg>'][i]
    
    co = sbsco[i]
    hi = sbs[i]
    
    
    colons = np.arange(co.shape[2])*(-0.125) + l0
    colats = np.arange(co.shape[1])*(0.125) + b0
    
    hilons = np.arange(hi.shape[2])*(-0.0833) + np.nanmean((l0,l1))+rhi
    hilats = np.arange(hi.shape[1])*(0.0833) + np.nanmean((b0,b1))-rhi
    
    for t in hdend[i]:
        if t.is_leaf:
            plt.imshow(np.nansum(t.get_mask()*hi,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilons[0],hilons[-1],hilats[0],hilats[-1]])
            plt.contour(colons,colats,np.nansum(co,axis=0),cmap='binary_r')
            plt.title(zcat['names'][i]+str(t.idx))
            plt.show()


for i in range(21):
    plot_hi_w_cocontour(i,all_dnds)


for i in range(21):
    plot_hi_w_cocontour(i,all_dnds_f)


zcat['names'] = ['aquila_r',
 'california',
 'camelopardis',
 'cepheus_f',
 'cepheus_n',
 'chamaeleon',
 'cma_ob1',
 'crossbones',
 'gemini_ob1',
 'hercules',
 'lacerta',
 'maddalena',
 'mon_ob1',
 'mon_r2',
 'orion_a',
 'orion_b',
 'orion_l',
 'perseus',
 'polaris',
 'rosette',
 'taurus']
