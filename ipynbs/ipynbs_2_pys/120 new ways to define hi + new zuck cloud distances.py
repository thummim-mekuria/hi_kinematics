#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

import numpy as np
import pandas as pd
import datacube as dc
import fun as fns

import matplotlib.pyplot as plt


hi4pi,dht = fits.open('both_same_grid_720d.fits')


h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


h3sig.shape


hlv_posb = np.nansum(h3sig[:,400:,:],axis=1)*0.125*1.82e18
hlv_negb = np.nansum(h3sig[:,:400,:],axis=1)*0.125*1.82e18


plt.contourf(np.log10(hlv_posb),extent=[360,-360,-100,100],levels=np.arange(20,22,0.25))
plt.contour(np.log10(hlv_posb),extent=[360,-360,-100,100],levels=[20.4],colors='red')

plt.xlim(180,-180)


plt.contourf(np.log10(hlv_negb),extent=[360,-360,-100,100],levels=np.arange(20,22,0.25))
plt.contour(np.log10(hlv_negb),extent=[360,-360,-100,100],levels=[20.4],colors='red')

plt.xlim(180,-180)


np.arange(0,800,100)


for i in np.arange(0,800,100):
    this = np.nansum(h3sig[:,i:i+100,:],axis=1)*0.125*1.82e18
    
    plt.contourf(np.log10(this),extent=[360,-360,-100,100],levels=np.arange(20,22,0.25))
    plt.contour(np.log10(this),extent=[360,-360,-100,100],levels=[20.4],colors='red')

    plt.xlim(180,-180)
    plt.show()


hmom1 = fns.get_mom1(h3sig,np.arange(155)*1.3-100)


cat = QTable(fits.open('post_flash_cat.fits')[1].data)


cat


cath1s = fits.open('mom1_hi.fits')


for i in range(22):
    fig,axs = plt.subplots(1,3,figsize=(30,7))
    long = axs[0].contourf(hmom1,extent=[360,-360,-50,50],
                 levels=np.linspace(cat['vel<km/s>'][i][0],
                                    cat['vel<km/s>'][i][1],
                                    10),
                cmap='Paired')
    
    axs[0].set_xlim(cat['lon_GMC<deg>'][i])
    axs[0].set_ylim(cat['lat_GMC<deg>'][i])
    fig.colorbar(long,ax=axs[0])
    
    
    short = axs[1].contourf(cath1s[i].data,
                 levels=np.linspace(cat['vel<km/s>'][i][0],
                                    cat['vel<km/s>'][i][1],
                                    10),
                cmap='Paired')
    

    fig.colorbar(short,ax=axs[1])
    
#     axs[2].plot(np.nanmean())
    
    long2 = axs[2].contourf(hmom1,extent=[360,-360,-50,50],
                 levels=np.arange(-100,100,3),
                cmap='binary')
    
    fig.colorbar(long2,ax=axs[2])
    axs[2].set_xlim(cat['lon_GMC<deg>'][i])
    axs[2].set_ylim(cat['lat_GMC<deg>'][i])

    fig.suptitle(cat['name'][i])
    plt.show()


zuck20 = pd.read_csv('Handbook_Distances_Zucker2020.dat',delim_whitespace=True)


zuck20.columns


zuck20.iloc[176:180]


olds = ['Aquila_S', 'CMa_OB1','California','Cam','Chamaeleon',\
       'Corona_Australis','Draco','Gem_OB1','Hercules','Lacerta','Lupus','Maddalena',\
       'Mon_OB1_NGC2264','Mon_R2','Ophiuchus','Orion','Orion_Lam','Pegasus',\
       'Perseus','Polaris','Rosette','Spider','Taurus','Ursa_Major']


bglst = {}
for i in range(23):
    this_old = olds[i]
    nms = []
    for i in range(325):
        if zuck20['name'][i] == this_old:
            nms.append(i)
            
    bglst[this_old] = nms


for k in bglst.keys():
    if k != 'Orion':
        print(k,' ',round(np.nanmean(zuck20[['d16','d50','d84']].iloc[bglst[k]])))
#         print(k,bglst[k])


for k in bglst.keys():
    if k == 'Orion':
        for i in bglst[k]:
            print(zuck20[['d16','d50','d84']].iloc[i])
#             print(np.nanmean(zuck20[['d16','d50','d84']][i]))
#         print(k,bglst[k])


zuck20[np.nanmean(zuck20[['d16','d50','d84']],axis=1)<300][zuck20['l']<-130]


plt.figure(figsize=(20,10))
plt.contourf(dld_negb*2.8,levels=np.arange(5e-2,2,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_negb*2.8,levels=[0.2],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(326):
    if zuck20['b'][i] < 0:
        if zuck20['l'][i]>180:
            zuck20['l'][i] -=360
        plt.scatter(zuck20['l'][i],np.nanmean((zuck20['d16'][i],zuck20['d50'][i],zuck20['d84'][i])),
                                              color='red',marker='+')
#         plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
#                     cat['distance<pc>'][i]),color='orange')


plt.xlim(180,-180)
plt.ylim(60,1250)


plt.figure(figsize=(20,10))
plt.contourf(dld_posb*2.8,levels=np.arange(5e-2,2,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_posb*2.8,levels=[0.2],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(326):
    if zuck20['b'][i] > 0:
        if zuck20['l'][i]>180:
            zuck20['l'][i] -=360
        plt.scatter(zuck20['l'][i],np.nanmean((zuck20['d16'][i],zuck20['d50'][i],zuck20['d84'][i])),
                                              color='red',marker='+')
#         plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
#                     cat['distance<pc>'][i]),color='orange')


plt.xlim(180,-180)
plt.ylim(60,1250)


dlb = np.nansum(dust,axis=0)


plt.figure(figsize=(20,10))
plt.contourf(dlb*2.8,levels=np.arange(5e-2,5,1e-1),extent=[180,-180,-90,90],cmap='binary',norm='log')
plt.scatter(zuck20['l'],zuck20['b'],marker='+',
            c=np.nanmean(zuck20[['d16','d50','d84']],axis=1),
           cmap='Paired',vmax=2000)

plt.xlim(180,-180)
plt.colorbar()


clb = np.nansum(d3sig,axis=0)*1.3*1.97e20


plt.figure(figsize=(20,10))
plt.contourf(np.log10(clb),levels=np.arange(20,23,0.25),extent=[360,-360,-50,50],cmap='binary')
plt.colorbar()


plt.scatter(zuck20['l'],zuck20['b'],marker='+',s=400,
            c=np.nanmean(zuck20[['d16','d50','d84']],axis=1),
           cmap='Paired')

plt.xlim(180,-180)
plt.show()


hpk_vels = pd.read_csv('hpk_list.csv')[[str(i) for i in range(21)]]


cpk_vels = pd.read_csv('cpk_list.csv')[[str(i) for i in range(13)]]


nm_hpks = np.nansum(hpk_vels/hpk_vels,axis=1)


nm_cpks = np.nansum(cpk_vels/cpk_vels,axis=1)


nm_hpks_im = nm_hpks.reshape(800,5760)


nm_cpks_im = nm_cpks.reshape(800,5760)


in_dame = zuck20[zuck20['b']>-30][zuck20['b']<30]
one_kpc = in_dame[np.nanmean(in_dame[['d16','d50','d84']],axis=1)<1100]


plt.figure(figsize=(20,10))
plt.contourf(nm_hpks_im,extent=[360,-360,-50,50],levels=np.arange(0,7,1),cmap='Paired',extend='both')
plt.colorbar()

plt.scatter(one_kpc['l'],one_kpc['b'],marker='+',
           color='white',vmax=2000)

plt.xlim(180,-180)


plt.figure(figsize=(20,10))
plt.contourf(nm_cpks_im,extent=[360,-360,-50,50],levels=np.arange(0,7,1),cmap='Paired',extend='both')
plt.colorbar()

plt.scatter(one_kpc['l'],one_kpc['b'],marker='+',
           color='white',vmax=2000)

plt.xlim(180,-180)


plt.contourf(np.log10(clb),extent=[360,-360,-50,50])
plt.xlim(-130,-180)
plt.ylim(-40,-10)

plt.scatter(one_kpc['l'],one_kpc['b'],marker='+',
           color='red',vmax=2000)


plt.contourf(np.log10(clb),extent=[360,-360,-50,50],levels=np.arange(18,22,0.5),cmap='Paired')
plt.colorbar()
plt.xlim(-130,-180)
plt.ylim(-40,-10)

plt.scatter(one_kpc['l'],one_kpc['b'],marker='+',
           color='red',vmax=2000)


clbn = np.nansum(dht.data,axis=0)


plt.contourf(np.log10(clbn*1.97e20*1.3),
             extent=[360,-360,-50,50],
             levels=np.arange(18,22,0.5),cmap='Paired')
plt.colorbar()
plt.xlim(-130,-180)
plt.ylim(-40,-10)

plt.scatter(one_kpc['l'],one_kpc['b'],marker='+',
           color='red',vmax=2000)


und_or = zuck20.query('-130 > l > -180').query('b<-20.2')


in_dame_fr = one_kpc.query('~(-130 > l > -180 & b<-20.2)')


plt.figure(figsize=(20,10))
plt.contourf(nm_cpks_im,extent=[360,-360,-50,50],levels=np.arange(0,7,1),cmap='Paired',extend='both')
plt.colorbar()

plt.scatter(in_dame_fr['l'],in_dame_fr['b'],marker='+',
           color='white',vmax=2000)

plt.xlim(180,-180)


len(np.unique(in_dame_fr['name']))


for i in in_dame_fr.query('-5<b<5 & -90<l<90').index:
    plt.contourf(np.log10(clb),extent=[360,-360,-50,50],levels=np.arange(20,22.5,0.25))
    plt.colorbar()
    
    plt.xlim(in_dame_fr['l'][i]+5,in_dame_fr['l'][i]-5)
    plt.ylim(in_dame_fr['b'][i]-5,in_dame_fr['b'][i]+5)
    
    plt.scatter(in_dame_fr['l'][i],in_dame_fr['b'][i],color='red',marker='x')
    plt.title(in_dame_fr['name'][i])
    plt.show()
    
    
    lo,li = ((np.array((in_dame_fr['l'][i]+2,in_dame_fr['l'][i]-2)) - 360)/-0.125).astype('int')
    bo,bi = ((np.array((in_dame_fr['b'][i]-2,in_dame_fr['b'][i]+2)) + 50)/0.125).astype('int')
    
    sbc = d3sig[:,bo:bi,lo:li]
    sbh = h3sig[:,bo:bi,lo:li]
    
    sbc = np.nanmean(d3sig[:,bo:bi,lo:li],axis=(1,2))
    sbh = np.nanmean(h3sig[:,bo:bi,lo:li],axis=(1,2))
                    
#     plt.plot(np.nanmean(sbh,axis=(1,2)))
#     plt.plot(50*np.nanmean(sbc,axis=(1,2)))
    
    plt.plot(sbh)
    plt.plot(50*sbc)
    plt.title(in_dame_fr['name'][i])
    plt.show()
    
    
#     plt.contourf(np.log10(np.nansum(sbc,axis=1)*0.125*1.97e20),levels=np.arange(20,22.5,0.25))
#     plt.show()
    
#     plt.contourf(np.log10(np.nansum(sbc,axis=2)*0.125*1.97e20),levels=np.arange(20,22.5,0.25))
#     plt.show()


# for when i come back
# roughlt the idea i have is to determine the radial extent of the hi based on:
# let's say that the edges of atomic envelopes are around ~0.2 mag, 
# with 3d dust maps, we are able to separate out emission along the line of sight
# so that we can only integrate over a range of distances, then see what the footprint is
# of the edges defining the Av=0.2 contour in the 3d->2d projected map (only summing over the d, 
# when doing ld / bd by only looking at the ranges where you expect the cloud, give the Av=0.2 
# contour -> basically the iterative process i had been doing before, except the iterating happens
# on the dust map, which is patchier / more easily we can draw clsoed contours,
# then use that projected extent in lb for the hi cloud, 
# a useful check for this^ is: calculate total hydrogen mass from the dust, then calculate the hi mass
# within that area + the h2 within the borders of the cloud/envelope
# ... dust mass should => hi+h2 mass (the delta m would be co dark h2)
# and if dust mass is smaller, we are getting extra hi and this method doesnt work
# -> and for the hi, now theres the question of do we define the velocity just in the channels with co
# or wider? -> maybe the heuristic is, for clouds where the blending forces you to just choose the 
# channels with co (the peaks are too close together), don't include them
# -> another heuristic is, does Av=0.2 mag footprint roughly coincide with NH ~ 4e20? (look at KMT to double check
# that these are the same boundary)


np.unique(zuck20['name'])


len(np.unique(zuck20['name']))


dust = fits.open('mean_and_std_lbd.fits')[1].data


plt.scatter(cat['name'],cat['lat_GMC<deg>'].transpose()[0])
plt.scatter(cat['name'],cat['lat_GMC<deg>'].transpose()[1])

plt.xticks(rotation='vertical')

plt.hlines([-30,-15],'aquila_s','taurus',color='grey')
plt.hlines([-25,-10],'aquila_s','taurus',color='red')
plt.show()


dld_negb = np.nansum(dust[:,:393,:],axis=1)


plt.figure(figsize=(20,10))
plt.contourf(dld_posb*2.8,levels=np.arange(5e-2,2,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_posb*2.8,levels=[0.2],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(22):
    if np.nanmean(cat['lat_GMC<deg>'][i]) > 0:
        plt.scatter(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i],color='red',marker='+')
        plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i]),color='orange')

plt.scatter(90.5,247,color='red',marker='x')
plt.annotate('pegasus',(90.5,247),color='orange')

plt.xlim(180,-180)
plt.ylim(60,1250)


plt.figure(figsize=(20,10))
plt.contourf(dld_negb*2.8,levels=np.arange(5e-2,2,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_negb*2.8,levels=[0.2],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(22):
    if np.nanmean(cat['lat_GMC<deg>'][i]) < 0:
        plt.scatter(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i],color='red',marker='+')
        plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i]),color='orange')

plt.scatter(90.5,247,color='red',marker='x')
plt.annotate('pegasus',(90.5,247),color='orange')

plt.xlim(180,-180)
plt.ylim(60,1250)


stt = [70,200,300,400,600]
end = [300,400,600,600,1000]


for i in range(5):
    st = stt[i]-68
    en = end[i]-68
    lb = np.nansum(dust[st:en,:,:],axis=0)
    
    plt.figure(figsize=(20,10))
    plt.contourf(lb*2.8,extent=[180,-180,-90,90],levels=np.arange(5e-2,2,1e-1))
    plt.colorbar()
    plt.contour(lb*2.8,levels=[0.2],colors=['white'],extent=[180,-180,-90,90],linewidths=0.8)
    plt.contour(lb*2.8,levels=[1],colors=['orange'],extent=[180,-180,-90,90],linewidths=0.8)

    plt.xlim(180,-180)
    plt.show()


plt.figure(figsize=(20,10))
plt.contourf(dld_negb*2.8,levels=np.arange(5e-2,9e-1,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_negb*2.8,levels=[0.2,0.5],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(22):
    if np.nanmean(cat['lat_GMC<deg>'][i]) < 0:
        plt.scatter(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i],color='red',marker='+')
        plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i]),color='orange')

plt.scatter(90.5,247,color='red',marker='x')
plt.annotate('pegasus',(90.5,247),color='orange')

plt.xlim(180,140)
plt.ylim(350,600)

plt.vlines([169.5,152],350,600,linestyle='--',color='orange')
plt.hlines([411,547],180,140,linestyle='--',color='orange')


plt.figure(figsize=(20,10))
plt.contourf(dld_negb*2.8,levels=np.arange(5e-2,2,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_negb*2.8,levels=[0.2],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(22):
    if np.nanmean(cat['lat_GMC<deg>'][i]) < 0:
        plt.scatter(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i],color='red',marker='+')
        plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i]),color='orange')

plt.scatter(90.5,247,color='red',marker='x')
plt.annotate('pegasus',(90.5,247),color='orange')

plt.xlim(180,140)
plt.ylim(350,600)

plt.vlines([170,151],350,600,linestyle='--',color='orange')
plt.hlines([408,551],180,140,linestyle='--',color='orange')


(np.array((169.5,152))/-0.23 + 786).astype('int')


cal_bd = np.nansum(dust[:,:,49:125],axis=2)


plt.figure(figsize=(20,10))
plt.contourf(cal_bd*2.8,levels=np.arange(5e-2,9e-1,1e-1),extent=[-90,90,68,1250])
plt.colorbar()
plt.contour(cal_bd*2.8,levels=[0.2,0.5],colors=['white'],extent=[-90,90,68,1250],linewidths=0.8)

plt.ylim(350,650)
plt.xlim(-20,10)

plt.vlines([-13,-1.5],350,600,linestyle='--',color='orange')
plt.hlines([411,547],-20,10,linestyle='--',color='orange')


#more zoomed in look
#start w california
st = 411-68
en = 547-68
lb = np.nansum(dust[st:en,:,:],axis=0)

plt.figure(figsize=(20,10))
plt.contourf(lb*2.8,extent=[180,-180,-90,90],levels=np.arange(5e-2,2,1e-1),cmap='binary_r')
# plt.contour(np.log10(clb),extent=[360,-360,-50,50],levels=np.arange(20,22,0.25),cmap='magma',linewidths=0.6)

# plt.colorbar()
plt.contour(lb*2.8,levels=[0.25,0.5],colors=['white'],extent=[180,-180,-90,90],linewidths=0.8)
# plt.contour(lb*2.8,levels=[1],colors=['orange'],extent=[180,-180,-90,90],linewidths=0.8)

plt.xlim(180,140)
plt.ylim(-20,10)


plt.vlines([169.5,152],-20,0,linestyle='--',color='orange')
plt.hlines([-12,-1.5],180,140,linestyle='--',color='orange')


plt.show()


plt.figure(figsize=(20,10))
plt.contourf(dld_negb*2.8,levels=np.arange(5e-2,9e-1,1e-1),extent=[180,-180,68,1250])
plt.colorbar()
plt.contour(dld_negb*2.8,levels=[0.15,0.2,0.5],colors=['white'],extent=[180,-180,68,1250],linewidths=0.8)

for i in range(22):
    if np.nanmean(cat['lat_GMC<deg>'][i]) < 0:
        plt.scatter(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i],color='red',marker='+')
        plt.annotate(cat['name'][i],(np.nanmean(cat['lon_GMC<deg>'][i]),
                    cat['distance<pc>'][i]),color='orange')

plt.scatter(90.5,247,color='red',marker='x')
plt.annotate('pegasus',(90.5,247),color='orange')

plt.xlim(100,90)
plt.ylim(450,600)

plt.vlines([97,94.5],350,600,linestyle='--',color='orange')
plt.hlines([485,503,537],100,90,linestyle='--',color='orange')


#more zoomed in look
#start w california
st = 485-68
en = 537-68
lb = np.nansum(dust[st:en,:,:],axis=0)

# plt.figure(figsize=(20,10))
plt.contourf(lb*2.8,extent=[180,-180,-90,90],levels=np.arange(5e-2,1,1e-1),cmap='binary_r')
plt.contour(np.log10(clb),extent=[360,-360,-50,50],levels=np.arange(20,22,0.25),cmap='magma',linewidths=0.6)

# plt.colorbar()
plt.contour(lb*2.8,levels=[0.2,0.5],colors=['white'],extent=[180,-180,-90,90],linewidths=0.8)
# plt.contour(lb*2.8,levels=[1],colors=['orange'],extent=[180,-180,-90,90],linewidths=0.8)

plt.xlim(99,92)
plt.ylim(-17,-5)


plt.vlines([97,94.5],-20,0,linestyle='--',color='orange')
plt.hlines([-15,-9],100,90,linestyle='--',color='orange')


plt.show()


(np.array((97,94.5))/-0.23 + 786).astype('int')


lac_bd = np.nansum(dust[:,:,364:375],axis=2)


plt.figure(figsize=(20,10))
plt.contourf(lac_bd*2.8,levels=np.arange(5e-2,9e-1,1e-1),extent=[-90,90,68,1250])
plt.colorbar()
plt.contour(lac_bd*2.8,levels=[0.1,0.2,0.5],colors=['white'],extent=[-90,90,68,1250],linewidths=0.8)

plt.ylim(450,600)
plt.xlim(-20,-5)

# plt.vlines([-13,-1.5],350,600,linestyle='--',color='orange')
plt.hlines([485,503,537],-20,0,linestyle='--',color='orange')


20*u.arcsec.to(u.deg)


cmom1 = fns.get_mom1(d3sig,np.arange(155)*1.3-100)


plt.figure(figsize=(20,10))
plt.contourf(cmom1,extent=[360,-360,-50,50],levels=np.arange(-60,100,5),cmap='RdGy')
plt.xlim(180,-180)

plt.colorbar()


plt.figure(figsize=(20,10))
plt.contourf(hmom1,extent=[360,-360,-50,50],levels=np.arange(-60,100,5),cmap='RdGy')
plt.xlim(180,-180)

plt.colorbar()


plt.hist(np.ravel(hmom1),bins=np.arange(-100,100,1),alpha=0.4)
plt.hist(np.ravel(cmom1),bins=np.arange(-100,100,1),alpha=0.4)

plt.yscale('log')

plt.show()


lons,lats = np.meshgrid(np.arange(5760)*-0.125+360,np.arange(800)*0.125-50)


plt.figure(figsize=(20,20))
plt.scatter(np.ravel(cmom1),np.ravel(lats),s=1e-1)
plt.scatter(np.ravel(hmom1),np.ravel(lats),s=1e-2,alpha=0.2)


plt.figure(figsize=(20,20))
plt.scatter(np.ravel(cmom1),np.ravel(lons),s=1e-1)
plt.scatter(np.ravel(hmom1),np.ravel(lons),s=1e-2,alpha=0.2)
plt.ylim(180,-180)


plt.contourf(hmom1<-50,extent=[360,-360,-50,50])
plt.xlim(180,-180)


plt.contourf(hmom1>50,extent=[360,-360,-50,50])
plt.xlim(180,-180)


plt.contourf(cmom1>50,extent=[360,-360,-50,50])
plt.xlim(180,-180)


plt.contourf(cmom1<-50,extent=[360,-360,-50,50])
plt.xlim(180,-180)


plt.figure(figsize=(20,15))
plt.contourf(hmom1-cmom1,cmap='Paired',levels=np.linspace(5,15,10),extent=[360,-360,-50,50])
plt.colorbar()
plt.xlim(-130,-180)
plt.ylim(-30,10)


plt.figure(figsize=(20,15))
plt.contourf(cmom1,cmap='Paired',levels=np.linspace(-5,30,30),extent=[360,-360,-50,50])
plt.colorbar()
plt.xlim(-130,-180)
plt.ylim(-30,10)


cpk_vels.columns = ['c'+str(n) for n in np.arange(13)]


hpk_vels.columns = ['h'+str(n) for n in np.arange(21)]


hpk_vels['l'] = np.ravel(lons)
hpk_vels['b'] = np.ravel(lats)


cpk_vels['l'] = np.ravel(lons)
cpk_vels['b'] = np.ravel(lats)


bothpks = hpk_vels.merge(cpk_vels,how='outer')


for i in range(22):
    lo = cat[i]['loh'] - 200
    bo = cat[i]['boh'] - 180
    li = cat[i]['lih'] + 200
    bi = cat[i]['bih'] + 200
    
    vo = cat[i]['vo']
    vi = cat[i]['vi']
    
    bigh = h3sig[vo:vi,bo:bi,lo:li]
    
    h1 = fns.get_mom1(bigh,np.arange(vo,vi)*1.3-100)
    c1 = fits.open('mom1_co.fits')[i].data
    h1old = fits.open('mom1_hi.fits')[i].data
    
    h2 = np.sqrt(fns.get_momN(bigh,np.arange(vo,vi)*1.3-100,2))
    c2 = fits.open('mom2_co.fits')[i].data
    h2old = fits.open('mom2_hi.fits')[i].data
    
    plt.title(cat['name'][i])
    plt.hist(np.ravel(h1),alpha=0.5,density=True,bins=100,color='black')
    plt.hist(np.ravel(h1old),alpha=0.5,density=True,bins=70,color='blue')
    plt.hist(np.ravel(c1),alpha=0.5,density=True,bins=70,color='red')
    plt.yscale('log')
    plt.show()
    
    plt.title(cat['name'][i])
    plt.hist(np.ravel(h2),alpha=0.5,density=True,bins=100,color='black')
    plt.hist(np.ravel(h2old),alpha=0.5,density=True,bins=70,color='blue')
    plt.hist(c2[c2>0.5],alpha=0.5,density=True,bins=70,color='red')
    
    
    plt.yscale('log')
    
    
    plt.show()
#     print('\n\n')


for i in range(22):
    plt.figure(figsize=(6,6))

    lo = cat[i]['loh'] - 200
    bo = cat[i]['boh'] - 180
    li = cat[i]['lih'] + 200
    bi = cat[i]['bih'] + 200
    
    vo = cat[i]['vo']
    vi = cat[i]['vi']
    
    bigh = h3sig[vo:vi,bo:bi,lo:li]
    
    h1 = fns.get_mom1(bigh,np.arange(vo,vi)*1.3-100)
    c1 = fits.open('mom1_co.fits')[i].data
    h1old = fits.open('mom1_hi.fits')[i].data
    
#     h2 = np.sqrt(fns.get_momN(bigh,np.arange(vo,vi)*1.3-100,2))
#     c2 = fits.open('mom2_co.fits')[i].data
#     h2old = fits.open('mom2_hi.fits')[i].data
    
    
#     plt.contourf(h2/np.nanmean(c2),extent=[lo*-0.125+360,li*-0.125+360,
#                                           bo*0.125-50,bi*0.125-50],
#                  cmap='Paired',levels=np.arange(4,16))
#     plt.colorbar()

#     plt.contour(c2,extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                            cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                 cmap='magma',linewidths=0.7)

    
#     plt.show()

    
#     plt.contourf(h1,extent=[lo*-0.125+360,li*-0.125+360,
#                                           bo*0.125-50,bi*0.125-50],
#                  levels=10,cmap='RdGy')
#     plt.colorbar()

#     plt.contour(c1,extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                            cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                 cmap='magma',linewidths=0.7)

    
#     plt.show()
    
#     plt.contourf(h1,extent=[lo*-0.125+360,li*-0.125+360,
#                                           bo*0.125-50,bi*0.125-50],
#                  levels=np.linspace(np.nanmin(c1),np.nanmax(c1),10),cmap='RdGy')
# #     plt.colorbar()

#     plt.contour(c1,extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                            cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                 cmap='magma',linewidths=0.7)
    
#     plt.show()
#     plt.plot(np.nanmean(d3sig[vo:vi,bo:bi,lo:li],axis=(1,2)))

    
#     plt.show()
    
    
#     plt.contourf(h1old,extent=[lo*-0.125+360,li*-0.125+360,
#                                           bo*0.125-50,bi*0.125-50],
#                  levels=np.linspace(np.nanmin(c1),np.nanmax(c1),10),cmap='PiYG')
#     plt.colorbar()

#     plt.contour(c1,extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                            cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                 cmap='magma',linewidths=0.7)

    
#     plt.show()
    
    
#     plt.title(cat['name'][i])
#     plt.hist(np.ravel(h1),alpha=0.5,density=True,bins=100,color='black')
# #     plt.hist(np.ravel(h1old),alpha=0.5,density=True,bins=70,color='blue')
#     plt.hist(np.ravel(c1),alpha=0.5,density=True,bins=70,color='red')
#     plt.yscale('log')
#     plt.show()
    
    
#     plt.scatter(np.nanmean(c1)-np.nanmean(h1),cat['vel_grad_GMC<km/s/pc>'][i],color='black')
#     plt.scatter(np.nanmean(c1)-np.nanmean(h1old),cat['vel_grad_GMC<km/s/pc>'][i],color='red')
    
#     plt.annotate(cat['name'][i],(np.nanmean(c1)-np.nanmean(h1),cat['vel_grad_GMC<km/s/pc>'][i]))
    
# plt.yscale('log')
# plt.show()

#     plt.hist(np.ravel(c1),alpha=0.5,density=True,bins=70,color='red')
#     plt.yscale('log')
#     plt.show()

    plt.contourf(np.nanargmax(bigh,axis=0)*1.3-vo,levels=np.linspace(np.nanmin(c1),np.nanmax(c1),10),
                 cmap='Paired')
    plt.colorbar()
    plt.show()
    
    plt.contourf(h1,levels=np.linspace(np.nanmin(c1),np.nanmax(c1),10),cmap='Paired')
    plt.colorbar()
    plt.show()

    
grrr im still so stumped on how to define hi radii
new idea - some more looking at linewidths? go back to how / why ra was defined, and see if i could
take that reasining and squish it into shape for here
orrr see if the pixels co-inciding w co emission will have the same first moment values/distribution
then step outwards incremently until that stops being the case? would just do a circular aperture in that case?
alsoo -> realizing that there is def some sort of largescale velocity gradient/ distribution ... maybe see what 
happens if that is subtracted/normalized by? -> also the co looks consistently more blueshifted than the surrounding hi -> first make sure the linecenter is being caught by first moment map ... and whether this means anything
more dust?


plt.contourf(np.nanargmax(bigh,axis=0)*1.3-100,levels=np.arange(-40,20,5),cmap='Paired')
plt.colorbar()
plt.show()


plt.contourf(h1,levels=np.arange(-40,20,5),cmap='Paired')
plt.colorbar()
plt.show()


fttr = fitting.LinearLSQFitter()
init = models.Linear1D()


for i in range(22):
    lo = cat[i]['loh'] - 200
    bo = cat[i]['boh'] - 180
    li = cat[i]['lih'] + 200
    bi = cat[i]['bih'] + 200
    
    vo = cat[i]['vo']
    vi = cat[i]['vi']
    
    bigh = h3sig[vo:vi,bo:bi,lo:li]
    
#     h1 = fns.get_mom1(bigh,np.arange(vo,vi)*1.3-100)
#     c1 = fits.open('mom1_co.fits')[i].data
#     h1old = fits.open('mom1_hi.fits')[i].data
    
#     h2 = np.sqrt(fns.get_momN(bigh,np.arange(vo,vi)*1.3-100,2))
    c2 = fits.open('mom2_co.fits')[i].data
#     h2old = fits.open('mom2_hi.fits')[i].data
    
    nh = np.nansum(bigh,axis=2)*0.125*1.82e18
    
    plt.contourf(np.log10(nh),extent=[bo*0.125-50,bi*0.125-50,vo*1.3-100,vi*1.3-100],
                 cmap='Paired',levels=np.arange(20,22,0.15))
    plt.colorbar()
    
    plt.contour(np.log10(nh),extent=[bo*0.125-50,bi*0.125-50,vo*1.3-100,vi*1.3-100],
                 colors='white',levels=[np.log10(4e20)])

#     plt.contour(c2,extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                            cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                 cmap='magma',linewidths=0.7)

    
    plt.show()

    
plt.contourf(np.log10(hlv_posb),extent=[360,-360,-50,50],levels=np.arange(20,22,0.15))
plt.xlim(180,-180)
plt.contour(np.log10(hlv_posb),extent=[360,-360,-50,50],levels=[20.4,21],colors='white')


plt.contourf(np.log10(hlv_negb),extent=[360,-360,-50,50],levels=np.arange(20,22,0.15))
plt.xlim(180,-180)
plt.contour(np.log10(hlv_negb),extent=[360,-360,-50,50],levels=[20.4,21],colors='white')


hlb = np.nansum(h3sig,axis=0)*1.3*1.82e18


hlb2 = np.nansum(hi4pi.data,axis=0)*1.3*1.82e18


mn = np.broadcast_to(np.nanmean(hlb2,axis=0),(800,5760))


plt.figure(figsize=(17,7))
plt.contourf(np.log2(hlb2/mn),extent=[360,-360,-50,50],levels=10)#,levels=np.arange(20,23,0.25))
plt.colorbar()
plt.contour(np.log2(hlb2/mn),extent=[360,-360,-50,50],levels=[0],colors='white')

plt.xlim(180,-180)
plt.show()


mdn = np.broadcast_to(np.nanmedian(hlb2,axis=0),(800,5760))


plt.figure(figsize=(17,7))
plt.contourf(np.log2(hlb2/mdn),extent=[360,-360,-50,50],levels=10)#,levels=np.arange(20,23,0.25))
plt.colorbar()
plt.contour(np.log2(hlb2/mdn),extent=[360,-360,-50,50],levels=[0],colors='white')
plt.contourf(np.log10(clb),extent=[360,-360,-50,50],levels=[20.5,23],colors='orange')
plt.xlim(180,-180)
plt.show()


for i in range(22):
    plt.figure(figsize=(20,8))
    vo = cat['vo'][i]
    vi = cat['vi'][i]
    
    thislb = np.nansum(hi4pi.data[vo:vi,:,:],axis=0)*1.3*1.82e18
    
    plt.contourf(np.log10(thislb),extent=[360,-360,-50,50],levels=np.arange(20,22.5,0.25))
    plt.contour(np.log10(thislb),extent=[360,-360,-50,50],levels=[20.4],colors='white')
    plt.contourf(np.log10(clb),extent=[360,-360,-50,50],levels=[20.5,23],colors='orange',alpha=0.7)
    plt.title(cat['name'][i])
    plt.xlim(180,-180)
    plt.show()
    
    
plt.figure(figsize=(17,7))
plt.contourf(np.log10(hlb2),extent=[360,-360,-50,50],levels=np.arange(20,23,0.25))
plt.colorbar()
plt.contour(np.log10(hlb2),extent=[360,-360,-50,50],levels=[20.4],colors='white')

plt.xlim(180,-180)
plt.show()


plt.contourf(np.log10(hlb),extent=[360,-360,-50,50],levels=np.arange(20,23,0.25))
plt.colorbar()
plt.contour(np.log10(hlb),extent=[360,-360,-50,50],levels=[20.4],colors='white')

plt.xlim(180,-180)
plt.show()


plt.hist(cat['radius_GMC<pc>'],bins=20)
plt.show()


plt.scatter(cat['distance<pc>'],cat['radius_GMC<pc>'])
plt.loglog()
for i in range(22):
    plt.annotate(cat['name'][i],(cat['distance<pc>'][i],cat['radius_GMC<pc>'][i]))
plt.show()


c2s = [fits.open('mom2_co.fits')[i].data for i in range(22)]


plt.hist(np.concatenate([np.ravel(h2s[i])/np.nanmean(c2s[i])
 for i in range(22)]),bins=500)

plt.xlim(0.5,)
plt.ylim(0,1300)

# plt.yscale('log')

plt.show()


plt.figure(figsize=(10,10))
[plt.hist(np.ravel(h2s[i])/np.nanmean(c2s[i][c2s[i]>0.5]),
          bins=500,alpha=0.5,density=True) for i in range(22)]

plt.xlim(2,)
plt.ylim(0,6)

# plt.yscale('log')

plt.show()


[(cat['name'][i],np.nanmean(h2s[i])) for i in range(22)]


plt.hist(np.concatenate([np.ravel(h2s[i])
 for i in range(22)]),bins=300)

plt.xlim(0.5,)
# plt.ylim(0,1300)

# plt.yscale('log')

plt.show()


cat


cat[1]


(cat[1]['vel<km/s>']+100)/1.3


cat[1][['loh','lih','boh','bih']]


clhi = h3sig[50:93,50:600,1300:1900]


ch1 = fns.get_mom1(clhi,np.arange(43)*1.3-35)


ch2 = np.sqrt(fns.get_momN(clhi,np.arange(43)*1.3-35,2))


cc1 = fits.open('mom1_co.fits')[1].data


cc2 = fits.open('mom2_co.fits')[1].data


plt.imshow(cc2<7e-1,origin='lower')
plt.colorbar()


len(cc2[cc2>7e-1])


len(cc2[cc2<7e-1])


1498/(1498+831)


plt.hist(np.ravel(cc2),density=True,alpha=0.5,bins=100)
# plt.xlim(0,1)
plt.show()


plt.hist(np.ravel(ch2),alpha=0.5,bins=100,density=True)
plt.hist(cc2[cc2>7e-1],alpha=0.5,bins=100,density=True)

plt.yscale('log')
plt.show()


plt.hist(np.ravel(ch1),density=True,alpha=0.5)
plt.hist(np.ravel(cc1),density=True,alpha=0.5)

plt.show()


plt.figure(figsize=(15,15))
plt.contourf(ch2/np.nanmean(cc2),extent=[197.5,122.5,-43.75,25],cmap='Paired',levels=10)
plt.colorbar()

plt.contour(cc2,extent=[166.5,154,-13,-6],cmap='magma',linewidths=0.7)

plt.xlim(189,122)

# plt.xlim(175,140)
# plt.ylim(-30,5)


plt.figure(figsize=(15,15))

plt.contourf(ch1,extent=[197.5,122.5,-43.75,25],levels=10)
plt.colorbar()

plt.contourf(cc1,extent=[166.5,154,-13,-6],cmap='magma',levels=10)
plt.colorbar()

plt.xlim(189,122)


np.array([50,600])*0.125-50


np.array([1300,1900])*-0.125+360


#hmmm when im scrolling through this i keep seeing what look like ring structures
for i in range(0,1160,10):
    plt.figure(figsize=(25,10))
    plt.contourf(np.nansum(dust[i:i+10],axis=0),levels=np.arange(1e-2,1e-1,1e-2))
    plt.title(str(i+68))
    plt.colorbar()
    plt.show()


from astrodendro import Dendrogram


Dendrogram


plt.contourf(dld_posb,levels=np.arange(1e-2,5e-1,1e-2))
plt.vlines(1320,0,1100)


plt.contourf(dld_negb,levels=np.arange(1e-2,5e-1,1e-2))
plt.vlines(1320,0,1100)


dust[:,:,1320:].size/1e9


orish_dnd = Dendrogram.compute(dust[:,:,1320:],min_value=1e-2,min_delta=3e-2,min_npix=20)


dust.shape


1572/2


(1320-786)*-0.23


for tr in orish_dnd.trunk:
    sbcb = tr.get_mask()*dust[:,:,1320:]
    plt.contourf(np.nansum(sbcb,axis=0),levels=np.arange(1e-2,5e-1,3e-2),
                extent=[-122,-180,-90,90])
    
    plt.ylim(-30,30)
    plt.show()


orish_dnd.trunk


plt.imshow(np.nansum(orish_dnd[18].get_mask()*dust[:,:,1320:],axis=1),
           vmin=1e-2,origin='lower',cmap=dc.cmap_min0('viridis'))
plt.ylim(300,500)
plt.xlim(30,230)


plt.imshow(np.nansum(orish_dnd[18].get_mask()*dust[:,:,1320:],axis=0),
           vmin=1e-2,origin='lower',cmap=dc.cmap_min0('viridis'))
plt.ylim(300,350)
plt.xlim(30,230)


plt.imshow(np.nansum(orish_dnd[8].get_mask()*dust[:,:,1320:],axis=0),
           vmin=1e-2,origin='lower',cmap=dc.cmap_min0('viridis'))
plt.ylim(300,380)
plt.xlim(30,230)


plt.imshow(np.nansum(orish_dnd[8].get_mask()*dust[:,:,1320:],axis=1),
           vmin=1e-2,origin='lower',cmap=dc.cmap_min0('viridis'))
plt.ylim(300,500)
plt.xlim(30,230)
