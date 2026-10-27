#!/usr/bin/env python
# coding: utf-8


import datacube as dc

from astropy.io import fits
from astropy.table import QTable
from astrodendro import Dendrogram

import matplotlib.pyplot as plt
import numpy as np


zcat = QTable(fits.open('final_cat.fits')[1].data)


dht = fits.open('DHT_CO_longer.fits')[0]
hi4pi = fits.open('HI4PI.fits')[0]


dnds = []
for cloud in zcat['names']:
    d = Dendrogram.load_from(cloud+'_dend.fits')
    dnds.append(d)


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
    
    hivels = np.arange(310)*1.288 - 200.
    covels = np.arange(155)*1.3 - 100.
    
    for t in hdend[i]:
        if t.is_leaf:
            plt.imshow(np.nansum(t.get_mask()*hi,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilons[0],hilons[-1],hilats[0],hilats[-1]])
            plt.contour(colons,colats,np.nansum(co,axis=0),cmap='binary_r')
            plt.title(zcat['names'][i]+str(t.idx))
            plt.show()


sbsco = []
for i in range(21):
    l0,l1 = ((zcat['lon<deg>'][i] - 180.125)/(-0.125)).astype('int')
    b0,b1 = ((zcat['lat<deg>'][i] + 30.125)/0.125).astype('int')
     
    sbco = dht.data[:,b0:b1,l0:l1]
    sbsco.append(sbco)


sbs = []
for i in range(21):
    hicube = hi4pi.data[:,:,1:-1]
#     hicubef = hi4pi_f.data[:,:,1:-1]
    
    lc,bc = np.nanmean(zcat['lon<deg>'][i]), np.nanmean(zcat['lat<deg>'][i])
    rhi = zcat['radius_HI<deg>'][i]

    lc = round((lc - 180.1666)/(-0.08333))
    bc = round((bc + 49.999)/(0.08333))
    rhi = round(rhi/0.0833)
    
    
    l0,l1 = int(lc-rhi),int(lc+rhi)
    b0,b1 = int(bc-rhi),int(bc+rhi)
    
    if l0 < 0:
        sb = np.concatenate((hicube[:,b0:b1,l0:],hicube[:,b0:b1,:l1]),axis=2)
#         sbf = np.concatenate((hicubef[:,b0:b1,l0:],hicubef[:,b0:b1,:l1]),axis=2)
        
    elif l1 > hicube.shape[2]:
        l1_n = hicube.shape[2] - l1
        sb = np.concatenate((hicube[:,b0:b1,l0:],hicube[:,b0:b1,:l1_n]),axis=2)
#         sbf = np.concatenate((hicubef[:,b0:b1,l0:],hicubef[:,b0:b1,:l1_n]),axis=2)
        
    else:
        sb = hicube[:,b0:b1,l0:l1]
#         sbf = hicubef[:,b0:b1,l0:l1]
    
    sbs.append(sb)    


hivels = np.arange(310)*1.288 - 200.
covels = np.arange(155)*1.3 - 100.


for i in range(21):
    plt.plot(hivels,np.nanmean(sbs[i],axis=(1,2))/np.nanmax(np.nanmean(sbs[i],axis=(1,2))),\
             color='black',linestyle='--')
    plt.plot(covels,zcat['vspectrum_GMC<K>'][i]/np.nanmax(zcat['vspectrum_GMC<K>'][i]),color='black')
    for lf in dnds[i]:
        if lf.is_leaf:
            this = np.nanmean(lf.get_mask()*dnds[i].data,axis=(1,2))
            plt.plot(hivels,this/np.nanmax(this),label=lf.idx)
    plt.legend()
    plt.title(zcat['names'][i])
    plt.xlim(-30,100)
    plt.show()


for i in range(21):
    plt.plot(hivels,np.nanmean(sbs[i],axis=(1,2))/np.nanmax(np.nanmean(sbs[i],axis=(1,2))),\
             color='black',linestyle='--')
    plt.plot(covels,zcat['vspectrum_GMC<K>'][i]/np.nanmax(zcat['vspectrum_GMC<K>'][i]),color='black')
    for lf in dnds[i]:
        if lf.is_leaf:
            this = np.nanmean(lf.get_mask()*dnds[i].data,axis=(1,2))
            plt.plot(hivels,this/np.nanmax(this),label=lf.idx)
    plt.legend()
    plt.title(zcat['names'][i])
    plt.xlim(-20,20)
    plt.show()


spec_match = [-1,-1,2,2,2,0,2,1,0,2,-1,3,2,0,-1,-1,0,0,-1,1,0]


for i in range(21):
    plot_hi_w_cocontour(i,dnds)


lb_match =[-1,-1,2,2,2,0,2?,1?,0,2,2?,3?,2?,0,-1,-1,0,0,-1?,1?,0]


def plot_hi_w_cocontour_lv(i,hdend):
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
    
    hivels = np.arange(310)*1.288 - 200.
    covels = np.arange(155)*1.3 - 100.
    
    for t in hdend[i]:
        if t.is_leaf:
            plt.imshow(np.nansum(t.get_mask()*hi,axis=1),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilons[0]*10,hilons[-1]*10,hivels[0],hivels[-1]])
            plt.contour(colons*10,covels,np.nansum(co,axis=1),cmap='binary_r')
            plt.title(zcat['names'][i]+str(t.idx))
            plt.ylim(-30,100)
            plt.show()


for i in range(21):
    plot_hi_w_cocontour_lv(i,dnds)


lv_match = [-1,-1,2,2,2,0,2,1,0,2,2?,3?,1?,0,1?,2?,0,0,3?,1,0]


def plot_hi_w_cocontour_bv(i,hdend):
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
    
    hivels = np.arange(310)*1.288 - 200.
    covels = np.arange(155)*1.3 - 100.
    
    for t in hdend[i]:
        if t.is_leaf:
            plt.imshow(np.nansum(t.get_mask()*hi,axis=2),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilats[0]*10,hilats[-1]*10,hivels[0],hivels[-1]])
            plt.contour(colats*10,covels,np.nansum(co,axis=2),cmap='binary_r')
            plt.title(zcat['names'][i]+str(t.idx))
            plt.ylim(-30,100)
            plt.show()


for i in range(21):
    plot_hi_w_cocontour_bv(i,dnds)


bv_match = [-1,-1,2,2,2,0,2,1,0,2,2,3,1?,0,1,-1,0,0,3?,1?,0]


def plot_spc_lb_lv_bv(i,hdend):
    l0,l1 = zcat['lon<deg>'][i]
    b0,b1 = zcat['lat<deg>'][i]
    rco = zcat['radius_GMC<deg>'][i]
    rhi = zcat['radius_HI<deg>'][i]
    
    co = sbsco[i]
    hi = sbs[i]
    
    
    colons = np.arange(co.shape[2])*(-0.125) + l0
    colats = np.arange(co.shape[1])*(0.125) + b0
    covels = np.arange(155)*1.3 - 100.
    
    hilons = np.arange(hi.shape[2])*(-0.0833) + np.nanmean((l0,l1))+rhi
    hilats = np.arange(hi.shape[1])*(0.0833) + np.nanmean((b0,b1))-rhi
    hivels = np.arange(310)*1.288 - 200.

    
    for t in hdend[i]:
        if t.is_leaf:
            fig,axs = plt.subplots(2,2,figsize=(8,8))
            
            this = np.nanmean(t.get_mask()*dnds[i].data,axis=(1,2))
            axs[0,0].plot(hivels,this/np.nanmax(this),label=lf.idx)
            axs[0,0].plot(covels,\
                          zcat['vspectrum_GMC<K>'][i]/np.nanmax(zcat['vspectrum_GMC<K>'][i]),color='black')
            axs[0,0].set_xlabel('v [km/s]')
            axs[0,0].set_xlim(-30,50)
            
            
            axs[0,1].imshow(np.nansum(t.get_mask()*hi,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilons[0],hilons[-1],hilats[0],hilats[-1]])
            axs[0,1].contour(colons,colats,np.nansum(co,axis=0),cmap='binary_r')
            axs[0,1].set_ylabel('b [deg]')
            axs[0,1].set_xlabel('l [deg]')
            
            
            axs[1,0].imshow(np.nansum(t.get_mask()*hi,axis=1),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilons[0]*10,hilons[-1]*10,hivels[0],hivels[-1]])
            axs[1,0].contour(colons*10,covels,np.nansum(co,axis=1),cmap='binary_r')
            axs[1,0].set_ylim(-30,50)
            axs[1,0].set_xlabel('lx10 [deg]')
            axs[1,0].set_ylabel('v [km/s]')
            
            axs[1,1].imshow(np.nansum(t.get_mask()*hi,axis=2),origin='lower',cmap=dc.cmap_min0('viridis'),\
                      extent=[hilats[0]*10,hilats[-1]*10,hivels[0],hivels[-1]])
            axs[1,1].contour(colats*10,covels,np.nansum(co,axis=2),cmap='binary_r')
            axs[1,1].set_ylim(-30,50)
            axs[1,1].set_ylabel('v [km/s]')
            axs[1,1].set_xlabel('bx10 [deg]')
            
            fig.suptitle(zcat['names'][i]+str(t.idx))
            
            plt.show()


for i in range(21):
    plot_spc_lb_lv_bv(i,dnds)


spec_match = [-1,-1, 2,2,2,0,1,1,0,2,2,1,  2,  0,0,0,0,0,2,1,0]
lb_match =   [-1,-1, 2,2,2,0,1,1,0,2,2,1,  2?, 0,0,0,0,0,2,1,0]
lv_match =   [-1,-1, 2,2,2,0,1,1,0,2,2,1,  1?, 0,0,0,0,0,2,1,0]
lv_match =   [-1,-1, 2,2,2,0,1,1,0,2,2,1,  1?, 0,0,0,0,0,2,1,0]


aq,cal  cam,cef,cef,chm,cma,crs,gem,her,lac,mad  mob mr2,ora,orb,orl,per,pol,ros,tau 


0,1     2,3,4,5,6,7,8,9,10,11               12           13,14,15,16,17,18,19,20


def plot_all_spc_lb_lv_bv(i,hdend):
    l0,l1 = zcat['lon<deg>'][i]
    b0,b1 = zcat['lat<deg>'][i]
    rco = zcat['radius_GMC<deg>'][i]
    rhi = zcat['radius_HI<deg>'][i]
    
    co = sbsco[i]
    hi = sbs[i]
    
    
    colons = np.arange(co.shape[2])*(-0.125) + l0
    colats = np.arange(co.shape[1])*(0.125) + b0
    covels = np.arange(155)*1.3 - 100.
    
    hilons = np.arange(hi.shape[2])*(-0.0833) + np.nanmean((l0,l1))+rhi
    hilats = np.arange(hi.shape[1])*(0.0833) + np.nanmean((b0,b1))-rhi
    hivels = np.arange(310)*1.288 - 200.

    
    for t in hdend[i]:
        fig,axs = plt.subplots(2,2,figsize=(8,8))

        this = np.nanmean(t.get_mask()*dnds[i].data,axis=(1,2))
        axs[0,0].plot(hivels,this/np.nanmax(this),label=lf.idx)
        axs[0,0].plot(covels,\
                      zcat['vspectrum_GMC<K>'][i]/np.nanmax(zcat['vspectrum_GMC<K>'][i]),color='black')
        axs[0,0].set_xlabel('v [km/s]')
        axs[0,0].set_xlim(-70,100)


        axs[0,1].imshow(np.nansum(t.get_mask()*hi,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'),\
                  extent=[hilons[0],hilons[-1],hilats[0],hilats[-1]])
        axs[0,1].contour(colons,colats,np.nansum(co,axis=0),cmap='binary_r')
        axs[0,1].set_ylabel('b [deg]')
        axs[0,1].set_xlabel('l [deg]')


        axs[1,0].imshow(np.nansum(t.get_mask()*hi,axis=1),origin='lower',cmap=dc.cmap_min0('viridis'),\
                  extent=[hilons[0]*10,hilons[-1]*10,hivels[0],hivels[-1]])
        axs[1,0].contour(colons*10,covels,np.nansum(co,axis=1),cmap='binary_r')
        axs[1,0].set_ylim(-70,100)
        axs[1,0].set_xlabel('lx10 [deg]')
        axs[1,0].set_ylabel('v [km/s]')

        axs[1,1].imshow(np.nansum(t.get_mask()*hi,axis=2),origin='lower',cmap=dc.cmap_min0('viridis'),\
                  extent=[hilats[0]*10,hilats[-1]*10,hivels[0],hivels[-1]])
        axs[1,1].contour(colats*10,covels,np.nansum(co,axis=2),cmap='binary_r')
        axs[1,1].set_ylim(-70,100)
        axs[1,1].set_ylabel('v [km/s]')
        axs[1,1].set_xlabel('bx10 [deg]')

        fig.suptitle(zcat['names'][i]+str(t.idx))

        plt.show()


dnds_f = []
for cloud in zcat['names']:
    d = Dendrogram.load_from(cloud+'_dend_flat.fits')
    dnds_f.append(d)           


sbs[12].shape


plt.plot(hivels,np.nanmean(sbs[12][:,50:,:],axis=(1,2)))
plt.plot(hivels,np.nanmean(sbs[12][:,70:,:],axis=(1,2)))
plt.plot(hivels,np.nanmean(sbs[12][:,:100,:],axis=(1,2)))

plt.plot(covels,70*zcat['vspectrum_GMC<K>'][12],color='black')
plt.xlim(-70,100)


#for auila, just do a cut at b~2, for california, cut at b~-5, for mob keep everything


for i in [0,1,12]:
    plot_all_spc_lb_lv_bv(i,dnds)


for i in [0,1,12]:
    plot_all_spc_lb_lv_bv(i,dnds_f)


for i,idx in enumerate([-1,-1, 2,2,2,0,1,1,0,2,2,1,-1,0,0,0,0,0,2,1,0]):
    if idx != -1:
        plt.imshow(np.nansum(dnds[i][idx].get_mask(),axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))
        print(i,idx)
        plt.show()
    

plt.imshow(np.nansum(dnds[18][1].get_mask()*sbs[18],axis=1),origin='lower',cmap=dc.cmap_min0('viridis'))


plt.plot(hivels,np.nanmean(dnds[18][1].get_mask()*sbs[18],axis=(1,2)))
plt.plot(hivels,np.nanmean(sbs[18],axis=(1,2)))


for i in range(21):
    plot_all_spc_lb_lv_bv(i,dnds)


#get rid of cepheus


idxs = [0,0, 2,2,2,0,1,1,0,2,2,1,0,0,0,0,0,0,2,1,0] #do 0 for aqr, cal, and mob but come back to change


sbs_mod = [sbs[i]*dnds[i][idxs[i]].get_mask() for i in range(21)]
sbs_mod[0] = sbs_mod[0][:,105:,:]
sbs_mod[1] = sbs_mod[1][:,:-103,:]


zcat['names'][6]


dnds[6].trunk[0].children


dir(dnds[19])


sbs_mod[19] = sbs[19]


3*0.043


vs2 = vs


zcat['vel<km/s>'] = [np.array(vs[i]) for i in range(21)]


zcat['vspectrum_HI<K>']=[np.nanmean(sbs_mod[i],axis=(1,2)) for i in range(21)]


np.array(zcat['vel<km/s>'])


[(covels[coveledges[i][0]],covels[coveledges[i][1]]) for i in range(21)]


coveledges


coveledges = [((zcat['vel<km/s>'][i]+100)/1.3).astype('int') for i in range(21)]


zcat['vel<km/s>'][0] = np.array([-68.624, 70])


zcat['vel<km/s>'][19] = np.array([-10 , 60])
# zcat['vel<km/s>'][11] = np.array([ 2.216, 60.   ])


zcat['vel<km/s>'][11]


zcat['names'][11]


[zcat['vspectrum_GMC<K>'][i][coveledges[i][0]:coveledges[i][1]] for i in range(21)]


zcat['vspectrum_GMC<K>'] = [zcat['vspectrum_GMC<K>'][i][coveledges[i][0]:coveledges[i][1]] for i in range(21)]


[len(zcat['vspectrum_HI<K>'][i]) for i in range(21)]


[len(zcat['vspectrum_GMC<K>'][i]) for i in range(21)]


zcat.remove_column('vspectrum_HI<K>')


for i in range(21):
    plt.plot(np.nanmean(sbs_mod[i]*(0.43< sbs_mod[i]),axis=(1,2)))
    plt.show()


for i in range(21):
    plt.plot(np.nanmean(sbs_mod[i],axis=(1,2)))
    plt.show()


import clouds as cld


import astropy.units as u


zcat['vspectrum_HI<K>']=[np.nanmean(sbs[i],axis=(1,2)) for i in range(21)]


zcat.remove_column('vel<km/s>')


zcat.write('final_cat.fits',overwrite=True)


id(sbs)


id(sbs_mod)


sbs[i].


for i in range(21):
    sig,ar = cld.surface_density(sbs_mod[i],'21cm')
    m = cld.mass(sig,zcat['distance<pc>'][i]*u.pc,ar)
    print(zcat['names'][i])
#     print('hi mod',m)
    sig2,ar = cld.surface_density(sbs[i],'21cm')
    m2 = cld.mass(sig2,zcat['distance<pc>'][i]*u.pc,ar)
    print('diff',100*(m-m2)/m2,'%')
    print('co',zcat['mass_GMC<kg>'][i])
    print('\n\n')


for i in range(21):
    plt.plot(zcat['vspectrum_GMC<K>'][i])
#     plt.plot(zcat['vspectrum_GMC<K>'][i][coveledges[i][0]:coveledges[i][1]])
    plt.title(zcat['names'][i])
    plt.show()


zcat['vspectrum_GMC<K>'].shape


-10 + (n)*1.3 


covels


vs


vs = []
for i,sb in enumerate(sbs_mod):
    psc = np.nanmean(sb,axis=(1,2))
#     print(psc)
#     here = np.where(psc>np.nanstd(psc))[0]
    here = np.where(psc>0)[0]
#     print(here)
    vs.append((hivels[here[0]],hivels[here[-1]]))
#     print(here)
#     plt.plot(here)
#     plt.title(zcat['names'][i])
# plt.show()
#     here2 = np.where(psc>0.13)[0]
#     sbs_mod[i] = sb[here[0]:here[-1],:,:]
#     plt.plot(hivels,psc)
#     plt.vlines((hivels[here[0]],hivels[here[-1]]),0,np.nanmax(psc),color='red')
# # #     plt.vlines((hivels[here2[0]],hivels[here2[-1]]),0,np.nanmax(psc),color='purple')

#     plt.plot(covels,100*zcat['vspectrum_GMC<K>'][i])
#     plt.title(zcat['names'][i])
#     plt.show()


pscs = [np.nanmean(sbs_mod[i],axis=(1,2)) for i in range(21)]


for p in pscs:
    plt.plot(p)
    plt.show()


for sb in sbs_mod:
    plt.imshow(np.nansum(sb,axis=0),origin='lower',cmap=dc.cmap_min0('viridis'))
    plt.show()


sbs_mod[1] = sbs_mod[1][:,:-103,:]


l0,l1 = zcat['lon<deg>'][i]
b0,b1 = zcat['lat<deg>'][i]
rco = zcat['radius_GMC<deg>'][i]
rhi = zcat['radius_HI<deg>'][i]

co = sbsco[i]
hi = sbs[i]


colons = np.arange(co.shape[2])*(-0.125) + l0
colats = np.arange(co.shape[1])*(0.125) + b0
covels = np.arange(155)*1.3 - 100.

hilons = np.arange(hi.shape[2])*(-0.0833) + np.nanmean((l0,l1))+rhi
hilats = np.arange(hi.shape[1])*(0.0833) + np.nanmean((b0,b1))-rhi
hivels = np.arange(310)*1.288 - 200.
