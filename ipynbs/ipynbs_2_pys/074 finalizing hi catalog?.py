#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
from astropy.table import QTable

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Ellipse

import datacube as dc


hi4pi,dht = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


xbn_mr2, per_tau, orions, others = QTable(fits.open('segmented_cats.fits')[1].data),\
                                   QTable(fits.open('segmented_cats.fits')[2].data),\
                                   QTable(fits.open('segmented_cats.fits')[3].data),\
                                   QTable(fits.open('segmented_cats.fits')[4].data)


xbn_mr2_ls = [3984,4072]
xbn_mr2_bs = [288,344]
xbn_mr2_vs = [81,90]


xbm2_lb_hi = np.nansum(hi4pi.data[xbn_mr2_vs[0]:xbn_mr2_vs[1],:,:],axis=0)
xbm2_lv_hi = np.nansum(hi4pi.data[:,xbn_mr2_bs[0]:xbn_mr2_bs[1],:],axis=1)
xbm2_bv_hi = np.nansum(hi4pi.data[:,:,xbn_mr2_ls[0]:xbn_mr2_ls[1]],axis=2)


xbm2_lb_co = np.nansum(dht.data[xbn_mr2_vs[0]:xbn_mr2_vs[1],:,:],axis=0)
xbm2_lv_co = np.nansum(dht.data[:,xbn_mr2_bs[0]:xbn_mr2_bs[1],:],axis=1)
xbm2_bv_co = np.nansum(dht.data[:,:,xbn_mr2_ls[0]:xbn_mr2_ls[1]],axis=2)


lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50
vels = np.arange(155)*1.3 - 100


fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

# lbounds, bbounds, vbounds = [(360,-360), (-50,50), (-75,50)]
lbounds, bbounds, vbounds = [(-130,-160), (-20,0), (-10,25)]


axs[0,0].imshow(xbm2_lb_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])
axs[0,0].contour(lons,lats,xbm2_lb_hi,levels=[400,700],colors=['white'])
axs[0,0].contour(lons,lats,xbm2_lb_co,levels=[3.5],colors=['black'],linewidths=[0.8])

# axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='x',color='red')
axs[0,0].scatter(np.nanmean(xbn_mr2['lon<deg>'],axis=1),np.nanmean(xbn_mr2['lat<deg>'],axis=1),
            marker='+',color='orange')

axs[0,0].set(xlim=lbounds,ylim=bbounds)


axs[0,1].imshow(xbm2_bv_hi.transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
axs[0,1].contour(vels,lats,xbm2_bv_hi.transpose(),levels=[2300,4900,6000],colors=['white'])
axs[0,1].contour(vels,lats,xbm2_bv_co.transpose(),levels=[12],colors=['black'],linewidths=[0.8])


# axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='x',color='red')

axs[0,1].scatter(np.nanmean(xbn_mr2['vel<km/s>'],axis=1),np.nanmean(xbn_mr2['lat<deg>'],axis=1),
            marker='+',color='orange')

axs[0,1].set(xlim=vbounds, ylim=bbounds)


axs[1,0].imshow(xbm2_lv_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
axs[1,0].contour(lons,vels,xbm2_lv_hi,levels=[2400,3000,3500],colors=['white'])
axs[1,0].contour(lons,vels,xbm2_lv_co,levels=[5],colors=['black'],linewidths=[0.8])


# axs[1,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['vel<km/s>'],axis=1),
#             marker='x',color='red')
axs[1,0].scatter(np.nanmean(xbn_mr2['lon<deg>'],axis=1),np.nanmean(xbn_mr2['vel<km/s>'],axis=1),
            marker='+',color='orange')

axs[1,0].set(xlim=lbounds, ylim=vbounds)


li,bi,vi = [-136,-15,0]
lr,br,vr = [-15,9,21]
axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


per_tau_ls = [1464,1624]
per_tau_bs = [200,304]
per_tau_vs = [75,87]


peta_lb_hi = np.nansum(hi4pi.data[per_tau_vs[0]:per_tau_vs[1],:,:],axis=0)
peta_lv_hi = np.nansum(hi4pi.data[:,per_tau_bs[0]:per_tau_bs[1],:],axis=1)
peta_bv_hi = np.nansum(hi4pi.data[:,:,per_tau_ls[0]:per_tau_ls[1]],axis=2)


peta_lb_co = np.nansum(dht.data[xbn_mr2_vs[0]:xbn_mr2_vs[1],:,:],axis=0)
peta_lv_co = np.nansum(dht.data[:,per_tau_bs[0]:per_tau_bs[1],:],axis=1)
peta_bv_co = np.nansum(dht.data[:,:,per_tau_ls[0]:per_tau_ls[1]],axis=2)


fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

lbounds, bbounds, vbounds = [(200,130), (-50,0), (-20,30)]
# lbounds, bbounds, vbounds = [(-110,-210), (-45,10), (-10,40)]


axs[0,0].imshow(peta_lb_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])
axs[0,0].contour(lons,lats,peta_lb_hi,levels=[200,400,700],colors=['white'])
axs[0,0].contour(lons,lats,peta_lb_co,levels=[3.5],colors=['black'],linewidths=[0.8])

axs[0,0].scatter(np.nanmean(per_tau['lon<deg>'],axis=1),np.nanmean(per_tau['lat<deg>'],axis=1),
            marker='x',color='orange')
# axs[0,0].scatter(np.nanmean(zcat['lon<deg>']+360,axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='x',color='orange')

axs[0,0].set(xlim=lbounds,ylim=bbounds)


axs[0,1].imshow(peta_bv_hi.transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
axs[0,1].contour(vels,lats,peta_bv_hi.transpose(),levels=[3000,5500,8500],colors=['white'])
axs[0,1].contour(vels,lats,peta_bv_co.transpose(),levels=[90],colors=['black'],linewidths=[0.8])

axs[0,1].scatter(np.nanmean(per_tau['vel<km/s>'],axis=1),np.nanmean(per_tau['lat<deg>'],axis=1),
            marker='x',color='orange')

axs[0,1].set(xlim=vbounds, ylim=bbounds)


axs[1,0].imshow(peta_lv_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
axs[1,0].contour(lons,vels,peta_lv_hi,levels=[2200,3800,5000],colors=['white'])
axs[1,0].contour(lons,vels,peta_lv_co,levels=[30],colors=['black'],linewidths=[0.8])

axs[1,0].scatter(np.nanmean(per_tau['lon<deg>'],axis=1),np.nanmean(per_tau['vel<km/s>'],axis=1),
            marker='x',color='orange')

axs[1,0].set(xlim=lbounds, ylim=vbounds)


li,bi,vi = [199,-47,-7]
lr,br,vr = [-55,43,23]
axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


orab = orions[[0,1]]
orl = orions[2]


orl


ors_ab_ls = [4024,4144]
ors_ab_bs = [232,320]
ors_ab_vs = [76,89]


orl_ls = [4160,4256]
orl_bs = [264,344]
orl_vs = [73,89]


orab_lb_hi = np.nansum(hi4pi.data[ors_ab_vs[0]:ors_ab_vs[1],:,:],axis=0)
orab_lv_hi = np.nansum(hi4pi.data[:,ors_ab_bs[0]:ors_ab_bs[1],:],axis=1)
orab_bv_hi = np.nansum(hi4pi.data[:,:,ors_ab_ls[0]:ors_ab_ls[1]],axis=2)


orab_lb_co = np.nansum(dht.data[ors_ab_vs[0]:ors_ab_vs[1],:,:],axis=0)
orab_lv_co = np.nansum(dht.data[:,ors_ab_bs[0]:ors_ab_bs[1],:],axis=1)
orab_bv_co = np.nansum(dht.data[:,:,ors_ab_ls[0]:ors_ab_ls[1]],axis=2)


orl_lb_hi = np.nansum(hi4pi.data[orl_vs[0]:orl_vs[1],:,:],axis=0)
orl_lv_hi = np.nansum(hi4pi.data[:,orl_bs[0]:orl_bs[1],:],axis=1)
orl_bv_hi = np.nansum(hi4pi.data[:,:,orl_ls[0]:orl_ls[1]],axis=2)


orl_lb_co = np.nansum(dht.data[orl_vs[0]:orl_vs[1],:,:],axis=0)
orl_lv_co = np.nansum(dht.data[:,orl_bs[0]:orl_bs[1],:],axis=1)
orl_bv_co = np.nansum(dht.data[:,:,orl_ls[0]:orl_ls[1]],axis=2)


fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

lbounds, bbounds, vbounds = [(-140,-165), (-25,0), (-5,20)]


axs[0,0].imshow(orab_lb_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])
axs[0,0].contour(lons,lats,orab_lb_hi,levels=[600],colors=['white'])
axs[0,0].contour(lons,lats,orab_lb_co,levels=[4],colors=['black'],linewidths=[0.8])


# axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='+',color='red')
axs[0,0].scatter(np.nanmean(orab['lon<deg>'],axis=1),np.nanmean(orab['lat<deg>'],axis=1),
            marker='x',color='orange')

axs[0,0].set(xlim=lbounds,ylim=bbounds)


axs[0,1].imshow(orab_bv_hi.transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
axs[0,1].contour(vels,lats,orab_bv_hi.transpose(),levels=[3800,6300,7500],colors=['white'])
axs[0,1].contour(vels,lats,orab_bv_co.transpose(),levels=[15],colors=['black'],linewidths=[0.8])

axs[0,1].scatter(np.nanmean(orab['vel<km/s>'],axis=1),np.nanmean(orab['lat<deg>'],axis=1),
            marker='x',color='orange')
# axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='+',color='red')
axs[0,1].set(xlim=vbounds, ylim=bbounds)


axs[1,0].imshow(orab_lv_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
axs[1,0].contour(lons,vels,orab_lv_hi,levels=[3500,4500],colors=['white'])
axs[1,0].contour(lons,vels,orab_lv_co,levels=[30],colors=['black'],linewidths=[0.8])

# axs[1,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['vel<km/s>'],axis=1),
#             marker='+',color='red')
axs[1,0].scatter(np.nanmean(orab['lon<deg>'],axis=1),np.nanmean(orab['vel<km/s>'],axis=1),
            marker='x',color='orange')

axs[1,0].set(xlim=lbounds, ylim=vbounds)


li,bi,vi = [-144,-22,0]
lr,br,vr = [-14,12,14]
axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

lbounds, bbounds, vbounds = [(-150,-180), (-25,0), (-10,20)]


axs[0,0].imshow(orl_lb_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])
axs[0,0].contour(lons,lats,orl_lb_hi,levels=[450,650,1000],colors=['white'])
axs[0,0].contour(lons,lats,orl_lb_co,levels=[1],colors=['black'],linewidths=[0.8])


# axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='+',color='red')
axs[0,0].scatter(np.nanmean(orl['lon<deg>']),np.nanmean(orl['lat<deg>']),
            marker='x',color='orange')

axs[0,0].set(xlim=lbounds,ylim=bbounds)


axs[0,1].imshow(orl_bv_hi.transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
axs[0,1].contour(vels,lats,orl_bv_hi.transpose(),levels=[3800,6300,7500],colors=['white'])
axs[0,1].contour(vels,lats,orl_bv_co.transpose(),levels=[15],colors=['black'],linewidths=[0.8])

axs[0,1].scatter(np.nanmean(orl['vel<km/s>']),np.nanmean(orl['lat<deg>']),
            marker='x',color='orange')
# axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#             marker='+',color='red')
axs[0,1].set(xlim=vbounds, ylim=bbounds)


axs[1,0].imshow(orl_lv_hi,origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
axs[1,0].contour(lons,vels,orl_lv_hi,levels=[2300,3500,4700],colors=['white'])
axs[1,0].contour(lons,vels,orl_lv_co,levels=[6],colors=['black'],linewidths=[0.8])

# axs[1,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['vel<km/s>'],axis=1),
#             marker='+',color='red')
axs[1,0].scatter(np.nanmean(orl['lon<deg>']),np.nanmean(orl['vel<km/s>']),
            marker='x',color='orange')

axs[1,0].set(xlim=lbounds, ylim=vbounds)


li,bi,vi = [-155,-20,-5]
lr,br,vr = [-19,14,20]
axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


ol0s,ol1s = np.array(others['l_pix']).transpose()
ob0s,ob1s = np.array(others['b_pix']).transpose()
ov0s,ov1s = np.array(others['v_pix']).transpose()


olbs_hi = [np.nansum(hi4pi.data[ov0s[i]:ov1s[i],:,:],axis=0) for i in range(9)]
olbs_co = [np.nansum(dht.data[ov0s[i]:ov1s[i],:,:],axis=0) for i in range(9)]

olvs_hi = [np.nansum(hi4pi.data[:,ob0s[i]:ob1s[i],:],axis=1) for i in range(9)]
olvs_co = [np.nansum(dht.data[:,ob0s[i]:ob1s[i],:],axis=1) for i in range(9)]

obvs_hi = [np.nansum(hi4pi.data[:,:,ol0s[i]:ol1s[i]],axis=2) for i in range(9)]
obvs_co = [np.nansum(dht.data[:,:,ol0s[i]:ol1s[i]],axis=2) for i in range(9)]


for i in [0]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

    lbounds, bbounds, vbounds = [(178,143), (-19,0), (-17,15)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[400,600,900],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[3.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[4300,5200,6000],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[25],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
#     axs[0,1].scatter(np.nanmean(per_tau['vel<km/s>'],axis=1),np.nanmean(per_tau['lat<deg>'],axis=1),
#                 marker='x',color='orange')
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=np.arange(1400,3000,500),colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[15],colors=['black'],linewidths=[0.8])

#     axs[1,0].scatter(np.nanmean(per_tau['lon<deg>'],axis=1),np.nanmean(per_tau['vel<km/s>'],axis=1),
#                 marker='+',color='red')
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [175,-16,-10]
    lr,br,vr = [-29,14,18]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


for i in [1]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(-40,-80), (-30,0), (-20,15)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[180,230],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[3.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[2000,3500,5000],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[25],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=np.arange(2000,3000,400),colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[15],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [-47,-23,-3]
    lr,br,vr = [-27,13,12]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


for i in [2]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(-155,-185), (-10,10), (-20,30)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[1400,1800],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[3.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[2000,3500,5000],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[25],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[2900,4300,5000],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[16],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
  

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [-159,-7,-9]
    lr,br,vr = [-23,13,29]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


for i in [3]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(65,35), (2,20), (-5,20)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[200,300,400],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[1.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[1800,2500,3000],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[1,5,10],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[1100,1200],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[1,5,15],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
    

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [61,5.5,-1]
    lr,br,vr = [-23,13,16]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


for i in [4]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(105,75), (-20,0), (-15,15)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[180,270,400],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[1.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[900,1400],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[1,5,10],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[2000,2300],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[1,5,15],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
    

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [98,-18,-7]
    lr,br,vr = [-21,10,15]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


for i in [5]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(-135,-152), (-10,5), (10,45)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[900,1150,1300],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[1.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[2500,3000],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[4,10],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[2000,2300],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[4,10],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
    

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [-138,-5,17]
    lr,br,vr = [-10,8,22]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


for i in [6]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(-150,-170), (-3,10), (-5,23)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=np.arange(400,1700,300),colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[1.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[1000,1500,2000,2500,3000,3500],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[4,10],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[800,1200,1600],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[4,10],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
    

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [-153,-1,-3]
    lr,br,vr = [-11,9,20]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


spec6 = np.array([np.nansum(hi4pi.data[:,ob0s[i]:ob1s[i],ol0s[i]:ol1s[i]],axis=(1,2)) for i in [6]])
spec6co = np.array([np.nansum(dht.data[:,ob0s[i]:ob1s[i],ol0s[i]:ol1s[i]],axis=(1,2)) for i in [6]])

plt.plot(vels,spec6[0])
plt.plot(vels,100*spec6co[0])

plt.vlines([-2,17],0,9e4)
plt.show()


for i in [7]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(130,110), (10,40), (-15,10)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[130,150,200],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[3.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[500,1000,1100,],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[4,10],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[1000,1600],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[4,10],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
    

    axs[1,0].set(xlim=lbounds, ylim=vbounds)
    
    axs[1,1].plot(vels,spec7hi[0])
    axs[1,1].plot(vels,30*spec7co[0])


    li,bi,vi = [128,22,-9]
    lr,br,vr = [-14,14,14]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


spec7hi = np.array([np.nansum(hi4pi.data[:,ob0s[i]:ob1s[i],ol0s[i]:ol1s[i]],axis=(1,2)) for i in [7]])
spec7co = np.array([np.nansum(dht.data[:,ob0s[i]:ob1s[i],ol0s[i]:ol1s[i]],axis=(1,2)) for i in [7]])

plt.plot(vels,spec7hi[0])
plt.plot(vels,30*spec7co[0])

plt.vlines([-9,5],0,9e4)
plt.show()


for i in [8]:
    fig,axs = plt.subplots(2,2,figsize=(10,7),layout='constrained')

#     lbounds, bbounds, vbounds = [(185,140), (-15,0), (-20,15)]
    lbounds, bbounds, vbounds = [(-143,-160), (-10,6), (-5,30)]


    axs[0,0].imshow(olbs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50])

    axs[0,0].contour(lons,lats,olbs_hi[i],levels=[700],colors=['white'])
    axs[0,0].contour(lons,lats,olbs_co[i],levels=[3.5],colors=['black'],linewidths=[0.8])


#     axs[0,0].scatter(np.nanmean(zcat['lon<deg>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')

    axs[0,0].set(xlim=lbounds,ylim=bbounds)


    axs[0,1].imshow(obvs_hi[i].transpose(),origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-100,100,-50,50])
    axs[0,1].contour(vels,lats,obvs_hi[i].transpose(),levels=[1000,1500,2000],colors=['white'])
    axs[0,1].contour(vels,lats,obvs_co[i].transpose(),levels=[4,10],colors=['black'],linewidths=[0.8])

    axs[0,1].scatter(np.nanmean(others['vel<km/s>'][i]),np.nanmean(others['lat<deg>'][i]),
                marker='x',color='orange')
    
#     axs[0,1].scatter(np.nanmean(zcat['vel<km/s>'],axis=1),np.nanmean(zcat['lat<deg>'],axis=1),
#                 marker='+',color='red')
    axs[0,1].set(xlim=vbounds, ylim=bbounds)


    axs[1,0].imshow(olvs_hi[i],origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100])
    axs[1,0].contour(lons,vels,olvs_hi[i],levels=[1000,1100,1300],colors=['white'])
    axs[1,0].contour(lons,vels,olvs_co[i],levels=[4,10],colors=['black'],linewidths=[0.8])

    
    axs[1,0].scatter(np.nanmean(others['lon<deg>'][i]),np.nanmean(others['vel<km/s>'][i]),
                marker='x',color='orange')
    

    axs[1,0].set(xlim=lbounds, ylim=vbounds)


    li,bi,vi = [-145,-8,3]
    lr,br,vr = [-12,9,23]
    axs[0,0].add_patch(Rectangle((li,bi),lr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[0,1].add_patch(Rectangle((vi,bi),vr,br,fill=None,edgecolor='white',linestyle='--'))
    axs[1,0].add_patch(Rectangle((li,vi),lr,vr,fill=None,edgecolor='white',linestyle='--'))


#when i come back: write the bounds of all the gmc-full HI clouds based on the traced out rectangles
#add in axis labels to all the plots
#maybe save the 2d zeroth moments? or the cubes defined by the bounds
#calculate masses etc for the gmc-full his


wco = np.nansum(dht.data,axis=0)


vchannel_ws = [23,14,20,22,15,16,29,12,18,14,20] #the width along v of the rectangles defining gmcfull his


np.nanmean(vchannel_ws),np.nanmedian(vchannel_ws) #characteristic number of velocity channels


firsts = np.arange(-100,100,18)
seconds = np.arange(-91,100,18)


firsts_int = ((firsts+100.)/1.3).astype('int')
seconds_int = ((seconds+100.)/1.3).astype('int')


first_slices = [(firsts_int[i],firsts_int[i+1]) for i in range(len(firsts)-1)]
second_slices = [(seconds_int[i],seconds_int[i+1]) for i in range(len(seconds)-1)]


import astropy.units as u
import astropy.constants as c
co_conv = ((u.K*1.3*u.km/u.s*1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
            * 2 * 1.36 * c.m_p).to(u.solMass/u.pc**2)).value
hi_conv = ((u.K*1.3*u.km/u.s*1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) * 1.36 * c.m_p).to(u.solMass/u.pc**2)).value


for i,ii in enumerate(first_slices):
    hilb = np.nansum(hi4pi.data[ii[0]:ii[1],:,:],axis=0)
    colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    plt.figure(figsize=(17,8))

    plt.imshow(hilb*hi_conv,origin='lower',norm='log',extent=[360,-360,-50,50])
    plt.contour(lons,lats,hilb*hi_conv,levels=[2.5,5,10,15,20],cmap='RdYlBu',zorder=0)
    plt.contourf(lons,lats,colb*co_conv,levels=[10,1e5],colors='black')
    
    
    plt.title(str(firsts[i])+' km/s to '+str(firsts[i]+18)+' km/s')
    plt.xlim(200,-200)
    plt.show()
    

for i,ii in enumerate(second_slices):
    hilb = np.nansum(hi4pi.data[ii[0]:ii[1],:,:],axis=0)
    colb = np.nansum(dht.data[ii[0]:ii[1],:,:],axis=0)
    
    plt.figure(figsize=(20,8))

    plt.imshow(hilb*hi_conv,origin='lower',norm='log',extent=[360,-360,-50,50])
    plt.contour(lons,lats,hilb*hi_conv,levels=[2.5,5,10,15,20],cmap='RdYlBu',zorder=0)
    plt.contourf(lons,lats,colb*co_conv,levels=[10,1e5],colors='black')
    
    
    plt.title(str(seconds[i])+' km/s to '+str(seconds[i]+18)+' km/s')
    plt.xlim(200,-200)
    plt.show()
    

fig,ax = plt.subplots(1,1,figsize=(20,))
for i in range(17):
