#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable

import fun as fns
import datacube as dc


hi4pi,dht = fits.open('both_same_grid_720d.fits')


h3sig35kms = (hi4pi.data*(hi4pi.data>0.043*3))[50:105]
d3sig35kms = (dht.data*(dht.data>0.2*3))[50:105]


h0 = np.nansum(h3sig35kms,axis=0)
c0 = np.nansum(d3sig35kms,axis=0)

h1 = fns.get_mom1(h3sig35kms,np.arange(55)*1.3-35)
c1 = fns.get_mom1(d3sig35kms,np.arange(55)*1.3-35)

h2 = fns.get_momN(h3sig35kms,np.arange(55)*1.3-35,2)
c2 = fns.get_momN(d3sig35kms,np.arange(55)*1.3-35,2)


zcat = QTable(fits.open('final_cat.fits')[1].data)


cham_scub_co = d3sig35kms[:,240:304,3312:3424]
cham_scub_hi = h3sig35kms[:,178:365,3274:3461]


def generate_patches(cube,params): 
    d,sidelen,lo,bo,ppc = params['d'],params['cloud_size_px'],\
                        params['lo'],params['bo'],params['patch_size_pc']

    pc10 = (ppc/d)*u.rad.to(u.deg)/0.125 
    #npix corresponding to 10pc at cloud distance

    npatch = sidelen//pc10 #how many patches you can make
    extra = sidelen%pc10 #extra pixels
    
    brdrds = [0]+ list(int(extra/2)+np.arange(1,npatch)*int(pc10)) +[sidelen]
    coords = np.array((brdrds[:-1],brdrds[1:])).transpose()

    slices_px = np.array([(i,j) for j in coords for i in coords])

    slices_lb = (slices_px.transpose().reshape(4,-1) + \
                np.broadcast_to(np.array([lo,bo,lo,bo]),
                (int(npatch**2),4)
                ).transpose()
                ).transpose().reshape(-1,2,2).swapaxes(1,2).astype('int')

    axcoords = np.flip(np.array(np.meshgrid(
                        np.arange(npatch),
                        np.arange(npatch))).transpose(),
                axis=0).reshape(-1,2).astype('int')
    
    
    specs = []
    for i in range(int(npatch**2)):
        lo,li = slices_lb[i][0]
        bo,bi = slices_lb[i][1]

        sp = np.nanmean(cube[:,bo:bi,lo:li],axis=(1,2))
        specs.append(sp)
    
    return specs,axcoords,ppc


ysos_cat = QTable(fits.open('J_A+A_674_A21_kyso.dat.gz.fits')[1].data)


#focus on finding sfr values in literature, it's okay if this information is not all coming from the same place
#start with the yso counts in the marton et al table, then calculate sfr from that number (look at lada lombardi
#alves)
#lso okay if the star formation tracer is something else other than ysos


plt.figure(figsize=(18,9))
plt.scatter(ysos_cat['GLON'],ysos_cat['GLAT'],s=1,color='red')
plt.scatter(ysos_cat['GLON']-360,ysos_cat['GLAT'],s=1,color='red')

plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],linewidths=0.6)

plt.xlim(180,0)
plt.ylim(-40,40)


plt.figure(figsize=(18,9))
plt.scatter(ysos_cat['GLON'],ysos_cat['GLAT'],s=1,color='red')
plt.scatter(ysos_cat['GLON']-360,ysos_cat['GLAT'],s=1,color='red')

plt.contour(c0,extent=[360,-360,-50,50],levels=[1.5],linewidths=0.6)

plt.xlim(0,-180)
plt.ylim(-40,40)


def generate_pars(cat,i):
    d = cat['distance<pc>'][i]
    lo = cat['l_pix_HI'][i][0]
    bo = cat['b_pix_HI'][i][0]
    
    csz = cat['l_pix_HI'][i][1] - cat['l_pix_HI'][i][0]
    
    patchsz = (zcat['radius_acc_HI<pc>'][i]//40 + 1)*10
    
    return {'d':d,'lo':lo,'bo':bo,'cloud_size_px':csz,\
            'patch_size_pc':patchsz}


zcat = QTable(fits.open('final_cat.fits')[1].data)


h0 = np.nansum(h3sig35kms,axis=0)
c0 = np.nansum(d3sig35kms,axis=0)


for i in range(17):
    hsps,axcs,ptchsz = generate_patches(h3sig35kms,generate_pars(zcat,i))
    csps,axcs,ptchsz = generate_patches(d3sig35kms,generate_pars(zcat,i))
    
    npatch = int(np.sqrt(len(hsps)))
    
    
    lo,li = zcat['l_pix_HI'][i]
    bo,bi = zcat['b_pix_HI'][i]
    
    plt.title(zcat['names'][i]+' '+str(ptchsz)+ ' pc patches')
    plt.imshow(h0[bo:bi,lo:li],origin='lower',
               extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50])
    plt.contourf(c0[bo:bi,lo:li],levels=np.arange(1.5,20,5),
           extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50])
    
    plt.contour(c0[bo:bi,lo:li],levels=[1.5],colors='white',
           extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50])
    
    
    fig,axs = plt.subplots(npatch,npatch,figsize=(2*npatch,2*npatch),
                           sharex=True,sharey=True)
    for j in range(len(hsps)):
        r,c = axcs[j]
        hmom2 = fns.get_momN(hsps[j],np.arange(55)*1.3-35,2)
        cmom2 = fns.get_momN(csps[j],np.arange(55)*1.3-35,2)

        axs[r,c].plot(np.arange(55)*1.3-35,hsps[j],
                 linewidth=0.7,alpha=0.8,color='gray',
                      label=round(np.sqrt(hmom2),2))

        axs[r,c].plot(np.arange(55)*1.3-35,20*csps[j],
                 linewidth=0.7,alpha=0.8,color='blue',
                     label = round(np.sqrt(cmom2),2))
                      
        axs[r,c].legend()

    for ax in axs.flat:
    #     ax.set_ylabel('T [K]')
    #     ax.set_xlabel('v [km/s]')
        ax.label_outer()
    # plt.plot(np.arange(55)*1.3-35,np.nanmean(cham_scub_hi,axis=(1,2)),
    #          color='black',linewidth=2.5)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()


plt.figure(figsize=(5,5))
plt.contourf(np.nansum(cham_scub_hi,axis=0),np.arange(100,1600,100),
           extent=[-49.25,-72.625,-27.75,-4.375])

plt.vlines([-53.125,-56.25,-59.375,-62.5,-65.625,-68.75],-27.75,-4.375,color='white')
plt.hlines([-23.875,-20.75,-17.625,-14.5,-11.375,-8.25],-49.25,-72.625,color='white')

plt.contourf(np.nansum(cham_scub_co,axis=0),extent=[-54,-68,-20,-12],levels=np.arange(1.5,20,5),
             cmap=dc.cmap_min0('binary_r'),zorder=10)

plt.xlim(-49.25,-72.625)


idxs = [0,31,56,81,106,131,156,187]


coords = np.array(([0,31],
[31,56],
[56,81],
[81,106],
[106,131],
[131,156],
[156,187]))


slices = np.array([(i,j) for j in coords for i in coords])


plt.figure(figsize=(5,5))
plt.contourf(np.nansum(cham_scub_hi,axis=0),np.arange(100,1600,100),
           extent=[-49.25,-72.625,-27.75,-4.375])

plt.vlines([-53.125,-56.25,-59.375,-62.5,-65.625,-68.75],-27.75,-4.375,color='white')
plt.hlines([-23.875,-20.75,-17.625,-14.5,-11.375,-8.25],-49.25,-72.625,color='white')

plt.contourf(np.nansum(cham_scub_co,axis=0),extent=[-54,-68,-20,-12],levels=np.arange(1.5,20,5),
             cmap=dc.cmap_min0('binary_r'),zorder=10)

plt.xlim(-49.25,-72.625)


slices_bothbig = (slices.transpose().reshape(4,-1) + 
 np.broadcast_to(np.array([3274,178,3274,178])\
                 ,(49,4)).transpose()).transpose().reshape(-1,2,2).swapaxes(1,2)


axcoords = np.array(np.meshgrid(np.arange(7),np.arange(7))).transpose()


axcoordsfxd = np.flip(axcoords.reshape(7,7,2),axis=0).reshape(-1,2)


#figure out how to match these two the l,b coords in the imshow


fig,axs = plt.subplots(7,7,figsize=(7,7),sharex=True,sharey=True)

for i in range(49):
    lo,li = slices_bothbig[i][0]
    bo,bi = slices_bothbig[i][1]
    
    
    r,c = axcoordsfxd[i]
    axs[r,c].plot(np.arange(55)*1.3-35,np.nanmean(h3sig35kms[:,bo:bi,lo:li],axis=(1,2)),
             linewidth=0.7,alpha=0.8,color='gray')
    
    axs[r,c].plot(np.arange(55)*1.3-35,20*np.nanmean(d3sig35kms[:,bo:bi,lo:li],axis=(1,2)),
             linewidth=0.7,alpha=0.8,color='blue')
    
for ax in axs.flat:
#     ax.set_ylabel('T [K]')
#     ax.set_xlabel('v [km/s]')
    ax.label_outer()
# plt.plot(np.arange(55)*1.3-35,np.nanmean(cham_scub_hi,axis=(1,2)),
#          color='black',linewidth=2.5)
plt.subplots_adjust(wspace=0, hspace=0)


dust = fits.open('mean_and_std_lbd.fits')[1].data


d0 = np.nansum(dust[150:300],axis=0)


plt.contourf(d0,norm='asinh',origin='lower',extent=[180,-180,-90,90],levels=np.arange(1e-2,2,5e-2))
plt.contour(c0,norm='asinh',origin='lower',extent=[360,-360,-50,50],levels=[0.5,1,1.5],linewidths=0.9)
# plt.contour(d0,norm='asinh',origin='lower',extent=[180,-180,-90,90],levels=[0.2],linewidths=0.9,colors='white')

plt.hlines([8,10],44,42)
plt.vlines([44,42],8,10)


# plt.colorbar()
plt.xlim(55,35)
plt.ylim(2,19)


dld = np.nansum(dust[:,420:485,:],axis=1)


dbd = np.nansum(dust[:,:,582:628],axis=2)


clv = np.nansum(d3sig35kms[:,432:544,:],axis=1)


cbv = np.nansum(d3sig35kms[:,:,2480:2560],axis=2)


plt.contourf(cbv,origin='lower',extent=[-50,50,-35,35],levels=[0.1,0.5,1,2])
plt.hlines(6,6,15)
plt.vlines([7,10],-10,30)


plt.colorbar()
plt.ylim(30,-10)
plt.xlim(5,15)


plt.contourf(clv,origin='lower',extent=[360,-360,-35,35],levels=[0.1,0.5,1,2])
plt.hlines(6,50,40)
plt.colorbar()
plt.ylim(35,-35)
plt.xlim(50,40)


dddl,ddddl = np.gradient(dld)


dddb,ddddb = np.gradient(dbd)


plt.contourf(dbd,extent=[-90,90,60,1200],levels=np.arange(7e-2,1,5e-2))
plt.contour(np.sqrt(dddb**2+ddddb**2),alpha=0.5,cmap='rainbow',linewidths=0.7,
             extent=[-90,90,60,1200],levels=[5e-3,0.01,0.02,0.05,0.07,0.1])
# plt.contourf(clv,extent=[360,-360,280,170],levels=[0.05,0.1,0.5,1])

# plt.hlines([195,230],6,12)
# plt.vlines([7,10],250,15)


plt.colorbar()
plt.ylim(190,260)
plt.xlim(5,20)


plt.contourf(dld,extent=[180,-180,60,1200],levels=np.arange(5e-2,5e-1,1e-2))
# plt.contour(np.sqrt(dddl**2+ddddl**2),alpha=0.5,cmap='rainbow',
#              extent=[180,-180,60,1200],levels=[5e-3,0.01,0.02,0.05,0.07,0.1])
plt.contourf(clv,extent=[360,-360,280,170],levels=[0.05,0.1,0.5,1])

plt.colorbar()
plt.ylim(190,260)
plt.xlim(49,38)


hlv = np.nansum(h3sig35kms[:,432:544,:],axis=1)
hbv = np.nansum(h3sig35kms[:,:,2480:2560],axis=2)


plt.contourf(hbv,extent=[-50,50,-35,35],levels=np.arange(500,8000,200))
plt.contour(cbv,extent=[-50,50,-35,35],levels=[1,2],linewidths=0.5)

plt.hlines(6,6,15)
plt.vlines([7,10],-10,30)


# plt.colorbar()
plt.ylim(30,-10)
plt.xlim(0,14)


plt.contourf(hlv,extent=[360,-360,-35,35],levels=np.arange(500,10000,500))
plt.colorbar()
plt.contour(clv,extent=[360,-360,-35,35],levels=[1,2],linewidths=0.5)

plt.xlim(50,40)
plt.ylim(30,-10)


#what is this saying about 1. the sudden widening then return to previous narrowness of the velocity spectra
#2. and a high drho/ds that corresponds to that
#would be interesting to look at this in bd/bv as well
#so ig whats going on is there is a filament of gas getting sucked in at 43deg?

#what can we learn by looking at velocity field x density gradient? -> ie in which direction is material
#building up and in which direction is it thinning out, ans how does that be a tracer for the motion / currents
#in the fluid? (think, seeing wind on top of water based on the wave it creates)

#do quivers for the drho/ds, then hw does that compare to lv/bv


(np.array([55,35])-180)/-0.223


(np.array([50,40])-360)/-0.125


(np.array([50,40])-180)/-0.223


(np.arange(786)*0.223-90)[393]


786/2


d0.shape


50/0.23


zcat[5]


top level narrative:
    revisiting the velocity gradients -> compare hi & co
    have a section for the clouds that are doing something interesting
    add in sfr & look at correlations with other properties
    look at correlation plots in nias angular momentum + m33 paper
    look at energetics (start by doing it for the cloud as a whole) -> look at how the line center is moving
    
    when thinking about associations -> look at the first moment maps of the entire complex
    (do the co separately, then for the hi that envelopes everybody) -> put the associations on the same
    correlation plots as the individual clouds and see if they follow the same pattern
    send hercules stuff to thavisha


for tomorrow/soon: 
    scatter plots looking at the specific angular momentum vs other stuff
    
    
priority 1 -> explain the velocity gradient stuff i writing
priority 2 -> get something to thavisha
priority 3 -> change patch sizes for the spectra
priority 4 -> angular momentum scatter plots
priority 5 -> channel maps -> approximate the ds9 scrolling through and seeing stuff move 
(choose the ones where interesting)
