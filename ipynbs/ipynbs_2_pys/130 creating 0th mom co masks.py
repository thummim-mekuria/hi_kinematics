#!/usr/bin/env python
# coding: utf-8

# when i come back - formalize/finish co masks, then calculate masses, then go from there
# with acc radius, spectrum, also for co mass, keep both and report comparison? & explicitly state that the gmc
# mass used to calculate the gmc mass is pre-selected (ie not the one that is used when comparing with hi
# ... and its just for velocity gradient etc that the co also uses the hi vels? sticks to their own?
# yeah i think that it makes sense that each cloud (co/hi) have their own v bounds, since they alr have diff l b bounds


import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.io import fits

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import fun as fns
import datacube as dc


cat = QTable(fits.open('sighhplease.fits')[1].data)


coco = fits.open('h2_h2sized.fits')
hihi = fits.open('hi_hisized.fits')


vels = np.arange(155)*1.3-100


[(
plt.plot(50*np.nanmean(coco[i].data,axis=(1,2))),
plt.plot(np.nanmean(hihi[i].data,axis=(1,2))),
plt.show()
) for i in range(22)]


c0s = [np.log10(np.nansum(coco[i].data,axis=0)*1.3*1.97e20) for i in range(22)]


c2s = [np.sqrt(fns.get_momN(coco[i].data,vels,2)) for i in range(22)]
c2s = [np.where(c0s[i]<20,np.nan,c2s[i]) for i in range(22)]


c1s = [fns.get_mom1(coco[i].data,vels) for i in range(22)]
c1s = [np.where(c0s[i]<20,np.nan,c1s[i]) for i in range(22)]


[(
    plt.contourf(c0s[i],levels=np.arange(20,22,0.25)),
    plt.colorbar(),
    plt.title(cat['name'][i]),
#     plt.ylim(60,120),
    plt.show()
 ) for i in range(22)]


[(
    plt.contourf(c1s[i],levels=15),
    plt.colorbar(),
    plt.title(cat['name'][i]),
#     plt.ylim(60,120),
    plt.show()
 ) for i in range(22)]


[(
    plt.contourf(c2s[i],levels=15),
    plt.colorbar(),
    plt.title(cat['name'][i]),
#     plt.ylim(60,120),
    plt.show()
 ) for i in range(22)]


[(
plt.hist(np.ravel(c1s[i]),bins=np.arange(-30,60,0.65),histtype='step',density=True),
plt.hist(np.ravel(fns.get_mom1(hihi[i].data,vels)),bins=np.arange(-30,60,0.65),
         histtype='step',density=True),
plt.show()


) for i in range(22)]


from gaussdecomp import spectrum,fitter,cube,driver,utils


aqh0 = spectrum.Spectrum(np.nanmean(hihi[0].data,axis=(1,2)),vels)


aqhofit['par'].shape[0]/3


help(fitter.gaussfitter)


aqhofit


utils.gplot(vels,np.nanmean(hihi[0].data,axis=(1,2)),aqhofit['par'])


# num_pks = [fitter.gaussfitter(spectrum.Spectrum(spec,vels),silent=True)['par'].shape[0]/3 
#           for spec in hihi[0].data.reshape(155,-1).transpose()]

##^ this took too long, maybe put in a .py and run in the background


dust = fits.open('mean_and_std_lbd.fits')[1].data


bigc0 = np.log10(fits.open('../not_rn_8/allsky_0th_co.fits')[0].data)


old_mss = QTable(fits.open('lastlastLAST_cat.fits')[1].data)['mass_GMC<Msun>']


new_mss = [fns.mass_n_area(coco[i].data*u.K,
                 cat['distance<pc>'][i]*u.pc,
                 'CO')[0].to(u.solMass).value
 
 for i in range(22)]


clvs = [np.log10(np.nansum(coco[i].data,axis=1)*0.125*1.97e20) for i in range(22)]


cbvs = [np.log10(np.nansum(coco[i].data,axis=2)*0.125*0.5*1.97e20) for i in range(22)]


[(
    plt.contourf(cbvs[i],levels=np.arange(19.5,21.5,0.15)),
    plt.colorbar(),
    plt.title(cat['name'][i]),
    plt.ylim(55,105),
    plt.show()
 ) for i in range(22)]


[(
    plt.contourf(clvs[i],levels=np.arange(19.5,21.5,0.15)),
    plt.colorbar(),
    plt.title(cat['name'][i]),
    plt.ylim(55,105),
    plt.show()
 ) for i in range(22)]


from astrodendro import Dendrogram


[coco[i].data.size/1e5 for i in range(22)]


cmps = ['binary','binary','Purples','Purples','Oranges']


for tr in aqdnd:
    if tr.is_leaf:
        plt.contourf(np.log10(np.nansum(tr.get_mask(),axis=0)*1.3*1.97e20),
                     levels=np.arange(20,21.5,0.15),cmap=cmps[tr.idx])
#     plt.show()


for tr in aqdnd:
    if tr.is_leaf:
        plt.contourf(np.log10(np.nansum(tr.get_mask(),axis=1)*0.125*1.97e20),
                     levels=np.arange(20,21.5,0.15),cmap=cmps[tr.idx])
#     plt.show()
plt.ylim(60,100)


for tr in aqdnd:
    if tr.is_leaf:
        plt.contourf(np.log10(np.nansum(tr.get_mask(),axis=2)*0.5*0.125*1.97e20),
                     levels=np.arange(20,21,0.15),cmap=cmps[tr.idx])
#     plt.show()
plt.ylim(60,100)


np.array(['aquila_s',0, 
 'california',1,
 'camelopardalis',0,
 'chamaeleon',0,
 'cma_ob1',1,
 'corona_australis',0,
 'crossbones',1,
 'gemini_ob1',1,
 'hercules',0,
 'lacerta',0,
 'lupus',1,
 'maddalena',1, 
 'mon_ob1',1, 
 'mon_r2',1,
 'ophiuchus',1, 
 'orion_a',0,
 'orion_b',1,
 'orion_l',1, 
 'perseus',1, 
 'polaris',1,
 'rosette',1,
 'taurus',1,]).reshape(-1,2)


np.log10(4e20),np.log10(2e20),


lbls = ['aqu',
 'cal',
 'cam',
 'cha',
 'cma',
 'cor',
 'cro',
 'gem',
 'her',
 'lac',
 'lup',
 'mad',
 'mob',
 'mr2',
 'oph',
 'ora',
 'orb',
 'orl',
 'per',
 'pol',
 'ros',
 'tau']


dnds = [Dendrogram.load_from(fname) for fname in ['dend_'+lbls[i]+'.fits' for i in range(22)]]


help(dnds[1].plotter())


fig,ax = plt.subplots(1,1)
dnds[1].plotter().plot_tree(ax)


dnds[1].trunk[0].children[0].children


len(dnds[1])


dnds[1][20].children


plt.contourf(d0_cal,levels=np.arange(0.3,1,0.1))
plt.xlim(0,150)
plt.ylim(300,400)


dnds[1][24].descendants


for lf in [dnds[1][27]]:
    msk0 = np.nansum(lf.get_mask(),axis=0)
    plt.contourf(msk0,levels=np.arange(1,10))#,cmap=dc.cmap_min0('Pastel1'))
    plt.title(cat['name'][i]+' '+str(lf.idx))
    plt.show()
# plt.colorbar()
#         plt.show()


len(dnds[11])


for lf in [dnds[11][3]]:
    msk0 = np.nansum(lf.get_mask(),axis=0)
    plt.contourf(msk0,levels=np.arange(1,10))#,cmap=dc.cmap_min0('Pastel1'))
    plt.show()
# plt.colorbar()
#         plt.show()


len(dnds[12])


dnds[12][3].children


dnds[12].trunk


for lf in [dnds[12][7]]:
    msk0 = np.nansum(lf.get_mask(),axis=0)
    plt.contourf(msk0,levels=10)#,cmap=dc.cmap_min0('Pastel1'))
    plt.colorbar()
    plt.show()
# plt.colorbar()
#         plt.show()


msks = [c0s[i] > 20.3 for i in range(22)]


msks[1] = c0s[1]>20.5


from skimage.morphology import erosion


erosion(msks[6])


help(erosion)


for i in range(22):
    plt.imshow(c0s[i] > 20.3,origin='lower',cmap='magma')
#     plt.contourf(np.where(~erosion(c0s[i]>19.5),np.nan,c0s[i]),origin='lower',levels=15)
#     plt.contourf(np.where(~erosion(msks[i]),np.nan,c0s[i]),origin='lower',levels=15)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    
    plt.show()


for i in range(22):
#     plt.imshow(c0s[i] > 20.3,origin='lower')
    plt.contourf(c0s[i],origin='lower',levels=[20,20.5,21,22],cmap='binary')
    plt.contourf(np.where(~erosion(c0s[i]>20),np.nan,c0s[i]),origin='lower',levels=15)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    
    plt.show()


i = 1
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.3)
plt.imshow(((mskx*0.5 + msky)<63)+\
            ((mskx*0.6 + msky)>135)+\
           ((mskx*-0.6 + msky)<-70)+\
           ((mskx*-0.6 + msky)>70),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*0.5 + msky)<63)+\
            ((mskx*0.6 + msky)>135)+\
           ((mskx*-0.6 + msky)<-70)+\
           ((mskx*-0.6 + msky)>70))*msks[i]


i = 13
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.3)
plt.imshow(((mskx*-1.2+ msky)<-60)+\
            ((mskx*0.3 + msky)<10)+\
           ((mskx*-0.2 + msky)>30)+\
           ((mskx*0.6 + msky)>65),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*-1.2+ msky)<-60)+\
            ((mskx*0.3 + msky)<10)+\
           ((mskx*-0.2 + msky)>30)+\
           ((mskx*0.6 + msky)>65))*msks[i]


i = 12
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5)
plt.imshow(((mskx*-0.75+ msky)<-27)+\
            ((mskx*0.7 + msky)<18),#+\
#            ((mskx*-0.2 + msky)>30)+\
#            ((mskx*0.6 + msky)>65),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*-0.75+ msky)<-27)+\
            ((mskx*0.7 + msky)<18))*msks[i]


plt.contourf(bigc0,levels=np.linspace(20.7,22,10),extent=[360,-360,-50,50])
plt.xlim(180,150)
plt.ylim(-30,0)
plt.hlines([-12,-11],190,140)


c0s[i].shape


(cat['lic']-cat['loc'])[21],(cat['bic']-cat['boc'])[21]


coco[21].data.shape


(np.array((cat['bic'],cat['boc']))*0.125-50)


i = 21
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.7,cmap=dc.cmap_min0('viridis'))
plt.imshow(((mskx*0.5+ msky)<33),#+\
#             (mskx>25)*(mskx<65)*((mskx*0.4 + msky)>73)+\
#            (mskx>63)*((mskx*-0.15 + msky)>40),
#            ((mskx*-0.2 + msky)>30)+\
#            ((mskx*0.6 + msky)>65),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*0.5+ msky)<33)+\
            (mskx>25)*(mskx<65)*((mskx*0.35 + msky)>75)+\
           (mskx>63)*((mskx*-0.15 + msky)>40))*msks[i]


i = 6
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(((mskx*-1.9+ msky)<-36)+\
            ((mskx*-0.8 + msky)>42),#+\
#            ((mskx*-0.2 + msky)>30)+\
#            ((mskx*0.6 + msky)>65),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


# msks[i] = ~(((mskx*-1.9+ msky)<-33)+\
#             ((mskx*-0.8 + msky)>42))*msks[i]


i = 7
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(((mskx*-0.6+ msky)<-15)+\
            ((mskx*0.7 + msky)<30)+\
           ((mskx*1.8 + msky)>182),#+\
# #            ((mskx*0.6 + msky)>65),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*-0.6+ msky)<-15)+\
            ((mskx*0.7 + msky)<30)+\
           ((mskx*1.8 + msky)>182))*msks[i]


i = 18
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(
    ((mskx*-2.3+ msky)>25)+\
            ((mskx*1.5 + msky)<70)+\
           ((mskx*1.8 + msky)>182),#+\
# #            ((mskx*0.6 + msky)>65),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*-2.3+ msky)>25)+\
            ((mskx*1.5 + msky)<70)+\
           ((mskx*1.8 + msky)>182))*msks[i]


i = 16
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(
    ((mskx*-1.3+ msky)>6)+\
            ((mskx*-1.5 + msky)<-55)+\
           ((mskx*0.5 + msky)<19)+\
           ((mskx*0.6 + msky)>95),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~(((mskx*-1.3+ msky)>6)+\
            ((mskx*-1.5 + msky)<-55)+\
           ((mskx*0.5 + msky)<19)+\
           ((mskx*0.6 + msky)>95))*msks[i]


i = 10
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(
    ((mskx*-1.8+ msky)>20),#+\
#             ((mskx*-1.5 + msky)<-55)+\
#            ((mskx*0.5 + msky)<19)+\
#            ((mskx*0.6 + msky)>95),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~((mskx*-1.8+ msky)>20)*msks[i]


i = 14
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(
    ((mskx*-1.+ msky)<-135),#+\
# #             ((mskx*-1.5 + msky)<-55)+\
# #            ((mskx*0.5 + msky)<19)+\
# #            ((mskx*0.6 + msky)>95),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))


msks[i] = ~((mskx*-1.+ msky)<-135)*msks[i]


i = 19
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(
    ((mskx*-1.+ msky)<-65)+\
            ((mskx*1.2 + msky)<60),#+\
# #            ((mskx*0.5 + msky)<19)+\
# #            ((mskx*0.6 + msky)>95),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = ~(((mskx*-1.+ msky)<-65)+\
            ((mskx*1.2 + msky)<60))*msks[i]


i = 4
msky,mskx = np.indices(msks[i].shape)
plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
plt.imshow(
    ((mskx*-0.7+ msky)<-13),#+\
#             ((mskx*1.2 + msky)<60),#+\
# #            ((mskx*0.5 + msky)<19)+\
# #            ((mskx*0.6 + msky)>95),
           
            origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = ~((mskx*-0.7+ msky)<-13)*msks[i] #go through and edit
#the masks to look like this (ie with the added regions delineated)
#then make mom 0 maps, calculate mass, move on with my  life


for i in range(22):
    plt.contourf(msks[i],levels=[-1,0])
    plt.show()


for i in [1,4,7,11,16,18,19,21]:
    msks[i] = erosion(msks[i])


for i in range(22):
    plt.contourf(msks[i],levels=[-1,0])
    plt.show()


msks = [c0s[i]>20.6 for i in range(22)]


for i in [0,9,17]:
    msks[i] = c0s[i]>20.


for i in [2]: #so teeny tiny, that 25% of the pixels are 19-20 cm^-2
    msks[i] = c0s[i]>19


for i in range(22):
    plt.contourf(msks[i]*c0s[i],levels=np.arange(20,22,0.15),cmap='Paired')
    plt.colorbar()
    plt.show()


for i in [4,7,12,13,16,20,21]:
    msks[i] = c0s[i]>20.7


for i in range(22):
    plt.contourf(msks[i]*c0s[i],levels=np.arange(20,22,0.15),cmap='Paired')
    plt.colorbar()
    plt.show()


for i in [4]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(c0s[i]<20.7+\
        ((mskx*-0.7+ msky)<-13.5) +\
        ((mskx*-1.8+ msky)<-55) +\
        (msky > 40) +\
        ((mskx*1.2 + msky)<25)+\
        ((mskx*0.2 + msky)<10)+\
        (mskx<20)*(msky>15)*(msky<21)+\
        (msky>30)*(mskx>25)*(mskx<32),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
        ((mskx*-0.7+ msky)<-13.5) +\
        ((mskx*-1.8+ msky)<-55) +\
        (msky > 40) +\
        ((mskx*1.2 + msky)<25)+\
        ((mskx*0.2 + msky)<10)+\
        (mskx<20)*(msky>15)*(msky<21)+\
        (msky>30)*(mskx>25)*(mskx<32)))


for i in [1]:
#     plt.figure(dpi=200)
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        (mskx>30)*((mskx*0.56+ msky)<70)+\
        (mskx<32)*(msky<55) +\
        (mskx<25)*((mskx*0.2+ msky)<60)+\
        (mskx>86)*(mskx<104)*((mskx*1.8+ msky)>241)+\
        ((mskx*1.1+ msky)>210)+\
        ((mskx*-0.5 + msky)>74)+\
        ((mskx*-0.3 + msky)<-30)+\
        (msky<35)*(mskx>70)*(mskx<102)+\
        (msky>56)*(mskx>103)*(mskx<140)+\
        (msky>95),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~((mskx>30)*((mskx*0.56+ msky)<70)+\
        (mskx<32)*(msky<55) +\
        (mskx<25)*((mskx*0.2+ msky)<60)+\
        (mskx>86)*(mskx<104)*((mskx*1.8+ msky)>241)+\
        ((mskx*1.1+ msky)>210)+\
        ((mskx*-0.5 + msky)>74)+\
        ((mskx*-0.3 + msky)<-30)+\
        (msky<35)*(mskx>70)*(mskx<102)+\
        (msky>56)*(mskx>103)*(mskx<140)+\
        (msky>95)))


for i in [6]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*-2.3+ msky)<-48) +\
        (mskx>23)*(msky>10)*(msky<15)+\
        (mskx>21)*(msky>8)*(msky<13)+\
        (msky>44)*(mskx<10),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(((mskx*-2.3+ msky)<-48) +\
        (mskx>23)*(msky>10)*(msky<15)+\
        (mskx>21)*(msky>8)*(msky<13)+\
        (msky>44)*(mskx<10)))


for i in [11]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*0.5+ msky)<17) +\
        (mskx>50)+\
        (mskx>45)*(msky>11)*(msky<25)+\
        (mskx<22)*(mskx>10)*(msky>31)+\
        (mskx<42)*(mskx>33)*(msky>27)+\
        (mskx<5)*(msky>32),
#         (msky>44)*(mskx<10),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(((mskx*0.5+ msky)<17) +\
        (mskx>50)+\
        (mskx>45)*(msky>11)*(msky<25)+\
        (mskx<22)*(mskx>10)*(msky>31)+\
        (mskx<42)*(mskx>33)*(msky>27)+\
        (mskx<5)*(msky>32)))


for i in [12]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*0.7+ msky)<16) +\
        ((mskx*-0.5+ msky)<-15)+\
        ((mskx*-0.5+ msky)>25)+\
        (msky>30),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
  
        ((mskx*0.7+ msky)<16) +\
        ((mskx*-0.5+ msky)<-15)+\
        ((mskx*-0.5+ msky)>25)+\
        (msky>30)
))


for i in [13]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*0.25+ msky)<11) +\
        ((mskx*0.85+ msky)>78) +\
        ((mskx*-1.5+ msky)<-76)+\
        ((mskx*-0.3+ msky)>25)+\
        (msky>32)*(mskx>44)*(mskx<52)+\
        (msky>23)*(mskx>15)*(mskx<25),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
  
        ((mskx*0.25+ msky)<11) +\
        ((mskx*0.85+ msky)>78) +\
        ((mskx*-1.5+ msky)<-76)+\
        ((mskx*-0.3+ msky)>25)+\
        (msky>32)*(mskx>44)*(mskx<52)+\
        (msky>23)*(mskx>15)*(mskx<25)
))


for i in [7]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*-0.5+ msky)<-16) +\
        ((mskx*-1.5+ msky)<-85) +\
        ((mskx*-1.3+ msky)>25) +\
        ((mskx*0.6+ msky)<20) +\
        ((mskx*2.5+ msky)>217)*(msky>40) +\
        (mskx>75)+\
        (msky<7)*(mskx>24)*(mskx<31)+\
        (mskx<40)*(msky>50)*(msky<70)+\
        (mskx<5)*(msky>10)*(msky<50)+\
        (mskx>67)*(msky>36)*(msky<48)+\
        ((mskx*-0.6+msky)>27)*(msky>45)*(msky<54),
#         (mskx>21)*(msky>8)*(msky<13)+\
#         (msky>44)*(mskx<10),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
  ((mskx*-0.5+ msky)<-16) +\
        ((mskx*-1.5+ msky)<-85) +\
        ((mskx*-1.3+ msky)>25) +\
        ((mskx*0.6+ msky)<20) +\
        ((mskx*2.5+ msky)>217)*(msky>40) +\
        (mskx>75)+\
        (msky<7)*(mskx>24)*(mskx<31)+\
        (mskx<40)*(msky>50)*(msky<70)+\
        (mskx<5)*(msky>10)*(msky<50)+\
        (mskx>67)*(msky>36)*(msky<48)+\
        ((mskx*-0.6+msky)>27)*(msky>45)*(msky<54)
))


for i in [18]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*1.9+ msky)<50)+\
        ((mskx*2+ msky)<70)*(msky>30)+\
        ((mskx*1.5+ msky)<72)*(msky>50)+\
        ((mskx*1.5+ msky)>117)+\
        ((mskx*-0.75+ msky)<-17)+\
        (mskx<15)*(msky>40)*(msky<50)+\
        (mskx>43)*(msky>37)*(msky<48),#+\
# #         (mskx>21)*(msky>8)*(msky<13)+\
# #         (msky>44)*(mskx<10),

                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
    ((mskx*1.9+ msky)<50)+\
        ((mskx*2+ msky)<70)*(msky>30)+\
        ((mskx*1.5+ msky)<72)*(msky>50)+\
        ((mskx*1.5+ msky)>117)+\
        ((mskx*-0.75+ msky)<-17)+\
        (mskx<15)*(msky>40)*(msky<50)+\
        (mskx>43)*(msky>37)*(msky<48)
))


for i in [16]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*1+ msky)<20)+\
        (msky<8)*(mskx<40)+\
        (msky>65)+\
        ((mskx*-1.3+ msky)<-40)+\
        (mskx<32)*(msky>35)*(msky<52)+\
        (msky>27)*(mskx>7)*(mskx<28)+\
        (mskx<20)*(msky>20)*(msky<52),#+\


                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
    ((mskx*1+ msky)<20)+\
        (msky<8)*(mskx<40)+\
        (msky>65)+\
        ((mskx*-1.3+ msky)<-40)+\
        (mskx<32)*(msky>35)*(msky<52)+\
        (msky>27)*(mskx>7)*(mskx<28)+\
        (mskx<20)*(msky>20)*(msky<52)
))


for i in [21]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*-1.7+msky)<-103)*(mskx<90)*(msky<23)+\
        (mskx>80)*(msky<18)+\
        ((mskx*0.6+ msky)<35)+\
        ((mskx*-1.5+ msky)>65)+\
        ((mskx*-0.6+ msky)>40)*(mskx>45)+\
        ((mskx*0.4+ msky)>97)+\
        ((mskx*0.3+ msky)>75)*(mskx<28)+\
        (msky>68)*(mskx>20)*(mskx<47)+\
        ((mskx*0.7+ msky)>110)*(mskx>50)*(mskx<90)*(msky>60)+\
        (mskx>65)*(msky>55)*(msky<63)+\
        (msky>70)*(mskx>45)*(mskx<63),


                origin='lower',cmap=dc.cmap_min0('binary_r'))

msks[i] = msks[i]*(~(
    ((mskx*-1.7+msky)<-103)*(mskx<90)*(msky<23)+\
        (mskx>80)*(msky<18)+\
        ((mskx*0.6+ msky)<35)+\
        ((mskx*-1.5+ msky)>65)+\
        ((mskx*-0.6+ msky)>40)*(mskx>45)+\
        ((mskx*0.4+ msky)>97)+\
        ((mskx*0.3+ msky)>75)*(mskx<28)+\
        (msky>68)*(mskx>20)*(mskx<47)+\
        ((mskx*0.7+ msky)>110)*(mskx>50)*(mskx<90)*(msky>60)+\
        (mskx>65)*(msky>55)*(msky<63)+\
        (msky>70)*(mskx>45)*(mskx<63)
))

msks[i] = erosion(msks[i])


for i in [10]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
        ((mskx*-1.6+ msky)>35),
                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
    ((mskx*-1.6+ msky)>35)
))


for i in [14]:
    msky,mskx = np.indices(msks[i].shape)
    plt.imshow(c0s[i],origin='lower',vmin=20.5,cmap=dc.cmap_min0('viridis'))
    plt.imshow(
         ((mskx*-2.3+ msky)<-295)+\
        ((mskx*-1.4+ msky)<-175),
                origin='lower',cmap=dc.cmap_min0('viridis'))

msks[i] = msks[i]*(~(
    ((mskx*-2.3+ msky)<-295)+\
        ((mskx*-1.4+ msky)<-175)
))


for i in range(22):
    plt.contourf(msks[i]*c0s[i],levels=np.arange(20,22,0.15),cmap='Paired')
    plt.colorbar()
    plt.show()


fns.make_fits_from_list([msks[i]+1 for i in range(22)],'h2_0thmom_masks.fits')


#once again, finish masks, mass etc
#perseus, gemini, maybe a couple others just had their bounds tweaked so will need to reload c0s etc


mskd0s = [msks[i]*c0s[i] for i in range(22)]


fns.make_fits_from_list(mskd0s,'masked0thmom.fits')


#from here, look at spectra to decide velocity range
#calculate mass of gmc, then update radius_HI if needed 
#mass, zeroth moment mass
#then all the first moment stuff
#decide whether to use a mask (maybe look at number of emission peaks per spectrum) in veloctiy/first moment
#space in addition to / instead of using velocity range on spectrum to define 'co-moving hi'
