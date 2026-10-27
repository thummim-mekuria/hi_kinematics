#!/usr/bin/env python
# coding: utf-8


from astrodendro import Dendrogram,periodic_neighbours
from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy_healpix import HEALPix as aphp
from astropy.modeling import models,fitting

import healpy as hp
from skimage.segmentation import find_boundaries

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import bootleg as btlg
import datacube as dc


# co_dend_new = Dendrogram.compute(dht.data[:,160:640,1360:4400],min_npix=500,
#                                  min_value=0.1,min_delta=2)
# btlg.dend_to_fits(co_dend_new,'decomposed_DHT_wcenter.fits')


hi4pi,dht = fits.open('both_same_grid_720d.fits')


hicat = pd.read_csv('hi4pi_dendleaves_5solMass.csv')
cocat = pd.read_csv('co_dendleaves.csv')


hi4pi = Dendrogram.load_from('decomposed_HI4PI.fits')


codend = Dendrogram.load_from('decomposed_DHT_wcenter.fits')


codend.data.shape


lfs = [lf.idx for lf in codend.leaves]


goodvs = [lf for lf in lfs
               if np.nanmean(codend[lf].indices()[0])*1.3-100 - 2*np.nanstd(codend[lf].indices()[0])*1.3 > -100.
               and np.nanmean(codend[lf].indices()[0])*1.3-100 + 2*np.nanstd(codend[lf].indices()[0])*1.3 < 100.]


goodbs = [codend[i].idx for i in goodvs
               if np.nanmean(codend[i].indices()[1])*0.125-50 - 2*np.nanstd(codend[i].indices()[1])*0.125 > -50.
               and np.nanmean(codend[i].indices()[1])*0.125-50 + 2*np.nanstd(codend[i].indices()[1])*0.125 < 50.]


goodlbvs = [codend[i].idx for i in goodbs
               if np.nanmean(codend[i].indices()[2])*-0.125+360 + 2*np.nanstd(codend[i].indices()[2])*-0.125 < 360.
               and np.nanmean(codend[i].indices()[2])*-0.125+360 - 2*np.nanstd(codend[i].indices()[2])*-0.125 > -360.]


cocat = pd.DataFrame(np.array(
    [(np.nanmean(codend[i].indices(),axis=1)+[0,160,1360])*(1.3,0.125,-0.125) + [-100,-50,360]
    for i in goodlbvs]),columns=['vel<km/s>','lat<deg>','lon<deg>'])


cocat['idx']=goodlbvs


cocat = cocat.sort_values('idx')
cocat = cocat[['idx','lon<deg>','lat<deg>','vel<km/s>']]


cocat['area<deg2>'] = [(0.125**2)*len(np.unique(codend[idx].indices()[1:],axis=1)[0]) 
 for idx in cocat['idx']]


cocat['radius<deg>'] = np.sqrt(cocat['area<deg2>']/np.pi)


cocat.to_csv('co_dendleaves.csv',index=False)


hi_dend_lbs = [np.nansum(hi4pi[i].get_mask()*hi4pi.data,axis=0) for i in hicat['idx']]


hi_dend_vspec = [np.nanmean(hi4pi[i].get_mask()*hi4pi.data,axis=(1,2)) for i in hicat['idx']]


hi_all_lb = np.nansum(hi_dend_lbs,axis=0)


co_dend_lbs = [np.nansum(codend[i].get_mask()*codend.data,axis=0) for i in cocat['idx']]


co_dend_vspec = [np.nanmean(codend[i].get_mask()*codend.data,axis=(1,2)) for i in cocat['idx']]


co_all_lb = np.nansum(co_dend_lbs,axis=0)


[(plt.plot(vels,spec)) for spec in hi_dend_vspec]
plt.plot(vels,0.1*np.nansum(hi_dend_vspec,axis=0),color='black')


fig,axs = plt.subplots(2,1,figsize=(15,18),dpi=200)

# plt.imshow(np.log10(hi_all_lb*1.828e18),extent=[190,-190,-50,50],origin='lower',
#              cmap=dc.cmap_min0('viridis'),vmin=np.log10(1e20))

axs[0].contour(np.arange(190,-190,-0.125),np.arange(-30,30,0.125),
            co_all_lb*1.97e20,levels=[1e20],colors='xkcd:crimson',linewidths=0.8)


axs[0].contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
            hi_all_lb*1.828e18,levels=[1e20,4e20],colors='black',
               linewidths=[0.5,0.8],linestyles=['--','solid'])

axs[0].set_xlim(185,-5)

axs[1].contour(np.arange(190,-190,-0.125),np.arange(-30,30,0.125),
            co_all_lb*1.97e20,levels=[1e20],colors='xkcd:crimson',linewidths=0.8)


axs[1].contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
            hi_all_lb*1.828e18,levels=[1e20,4e20],colors='black',
               linewidths=[0.5,0.8],linestyles=['--','solid'])

axs[1].set_xlim(5,-185)

for ax in axs.flat:
    ax.set(xlabel='galactic longitude [deg]',ylabel='galactic latitude [deg]')
    # plt.xlabel('galactic longitude [deg]')
    # plt.ylabel('galactic latitude [deg]')
# # plt.colorbar(orientation='vertical')


#hmmm from here it looks like a lot of the leaves dont have enough body/substance to integrate to larger
#column densities --> gotta choose branches / larger structures instead


#the next few series of cells is choosing structures from the dend, not based on
#being a leaf, but on having npix<5e5 and then removing redundancy by only choosing the 
#lowest level structure if multiple nested ones have npix<5e5
#but then the same constraints on centrral l,b,v and density will be enforced
#5e5 specificaly bc thats one we stop getting the structures that span the entire galaxy


d_lv_npx_lfs = np.array([(i.idx,i.level,i.get_npix()) for i in hi4pi.leaves]).transpose()


plt.scatter(id_lv_npx_al[0],id_lv_npx_al[1],marker='+')# if i.level==1]
plt.scatter(lf_levels[0],lf_levels[1],marker='x',color='orange')


plt.scatter(id_lv_npx_al[0],id_lv_npx_al[2],marker='+')
plt.scatter(d_lv_npx_lfs[0],d_lv_npx_lfs[2],marker='+')
# plt.yscale('log')
plt.loglog()


id_lv_npx_al = np.array([(i.idx,i.level,i.get_npix()) for i in hi4pi]).transpose()


smols = [i.idx for i in hi4pi if i.get_npix()<5e5]


big_parents = [smol for smol in smols if hi4pi[smol] \
               not in hi4pi.trunk and hi4pi[smol].parent.idx not in smols]


okay_qm = np.array([i for i in smols if i in big_parents or hi4pi[i] in hi4pi.trunk])


[(plt.imshow(np.nansum(hi4pi[i].get_mask(),axis=0),
             extent=[190,-190,-50,50],cmap=dc.cmap_min0('viridis')),
) for i in okay_qm]


[(plt.imshow(np.nansum(i.get_mask(),axis=0),
             extent=[190,-190,-50,50],cmap=dc.cmap_min0('viridis')),
) for i in hi4pi if i.get_npix()<5e5 and not i.is_leaf]


[(plt.imshow(np.nansum(i.get_mask(),axis=0),
             extent=[190,-190,-50,50],cmap=dc.cmap_min0('viridis')),
) for i in hi4pi if i.get_npix()>5e5]


hi_dend_lbs_bronchs = [np.nansum(hi4pi[i].get_mask()*hi4pi.data,axis=0) for i in okay_qm]


hi_all_lb_bronchs = np.nansum(hi_dend_lbs_bronchs,axis=0)


fig,axs = plt.subplots(2,1,figsize=(15,18),dpi=200)

# plt.imshow(np.log10(hi_all_lb*1.828e18),extent=[190,-190,-50,50],origin='lower',
#              cmap=dc.cmap_min0('viridis'),vmin=np.log10(1e20))

axs[0].contour(np.arange(190,-190,-0.125),np.arange(-30,30,0.125),
            co_all_lb*1.97e20,levels=[1e20],colors='xkcd:crimson',linewidths=0.8,alpha=0.5)


axs[0].contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
            hi_all_lb_bronchs*1.828e18,levels=[4e20],colors='black',
               linewidths=[0.8],linestyles=['solid'])

# axs[0].contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
#             hi_all_lb*1.828e18,levels=[2e20,4e20],colors='orange',
#                linewidths=[0.5,0.8],linestyles=['--','solid'])

axs[0].set_xlim(185,-5)


axs[1].contour(np.arange(190,-190,-0.125),np.arange(-30,30,0.125),
            co_all_lb*1.97e20,levels=[1e20],colors='xkcd:crimson',linewidths=0.8,alpha=0.5)


axs[1].contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
            hi_all_lb_bronchs*1.828e18,levels=[4e20],colors='black',
               linewidths=[0.8],linestyles=['solid'])

# axs[1].contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
#             hi_all_lb*1.828e18,levels=[2e20,4e20],colors='orange',
#                linewidths=[0.5,0.8],linestyles=['--','solid'])

axs[1].set_xlim(5,-185)

for ax in axs.flat:
    ax.set(xlabel='galactic longitude [deg]',ylabel='galactic latitude [deg]')
    # plt.xlabel('galactic longitude [deg]')
    # plt.ylabel('galactic latitude [deg]')
# # plt.colorbar(orientation='vertical')


his = fits.open('hi_clouds.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


zcat[10]


plt.figure(figsize=(10,10))
plt.imshow(np.nansum(his[10].data,axis=0)*1.82e18,origin='lower',vmin=4e20,cmap=dc.cmap_min0('viridis'),
           extent=[3980*-0.125+360,4163*-0.125+360,232*0.125-50,264*0.125-50])


plt.figure(figsize=(10,10))
plt.imshow(np.nansum(his[10].data[27:27+19,:,:],axis=0)*1.82e18,origin='lower',
           vmin=4e20,cmap=dc.cmap_min0('viridis'),
           extent=[3980*-0.125+360,4163*-0.125+360,232*0.125-50,264*0.125-50])


plt.imshow(np.sqrt(np.nansum(his[10].data,axis=1)),origin='lower',
           extent=[3980*-0.125+360,4163*-0.125+360,50*1.3-100,126*1.3-100])


h2s=fits.open('h2_clouds.fits')


plt.figure(figsize=(10,10))
plt.imshow(np.nansum(h2s[10].data,axis=0),origin='lower',
           extent=[-143,-155,-21,-17])


plt.imshow(np.sqrt(np.nansum(h2s[10].data,axis=1)),origin='lower',
           extent=[-143,-155,0,15])


#move this stuff (esp the functions!) to 89


vspec = np.nanmean(hi4pi.data,axis=(1,2))


bspec = np.nanmean(hi4pi.data,axis=(0,2))


ideal_bv = np.outer(vspec,bspec)


plt.imshow(np.outer(vspec,bspec),norm='asinh')


plt.imshow(np.nanmean(hi4pi.data,axis=0)/ideal_lb,norm='asinh')


ideal_lb = np.broadcast_to(bspec,(3040,800)).transpose() 


plt.imshow(ideal_lb,norm='asinh')


plt.imshow(np.nanmean(hi4pi.data,axis=2)/ideal_bv)


plt.imshow(np.nanmean(hi4pi.data,axis=1))
plt.xlim(500,2500)


plt.imshow(np.broadcast_to(vspec,(3040,155)).transpose(),norm='asinh')


np.con


models.Gaussian1D(20,0,20)


plt.plot(vels,vspec)
plt.plot(vels,models.Lorentz1D(22,0,20)(vels))


def lor(y):
    return models.Lorentz1D(20,0,20)(y)


lons = np.arange(360,-360,-0.125)*u.deg


vels = np.arange(-100,101,1.3)


lor(vels)


def sinu(x):
    return 100*np.sin(2*x)


def z(x,y):
    X,Y = np.meshgrid(x,y)
    return sinu(X)/(np.abs(X)*Y) + 1e-2*Y
    

plt.imshow((z(lons.to(u.rad).value,vels)),extent=[360,-360,-100,100],origin='lower',vmin=-2,vmax=2)
# plt.colorbar()


#hmmm i need to figure out what is the function (prob has to do with the rotation curve) that describes l,v


plt.imshow(hi4pi.data[:,400,:]*1.3-100,extent=[360,-360,-100,100],origin='lower')


plt.imshow(np.nansum(hi4pi.data,axis=1),extent=[190,-190,-100,100])


plt.plot(z(lons,))


plt.plot(lons,y(lons))


plt.imshow(np.array(np.meshgrid(lons,y(lons)))[1])


np.nanargmax(hi4pi.data[:,400,:],axis=0)


plt.plot(np.nanargmax(hi4pi.data[:,400,:],axis=0)*1.3-100)
plt.plot(50*np.sin(2*np.arange(190,-190,-0.125)*u.deg))
plt.xlim()


ideal_lqm = 60*np.sin(2*np.arange(190,-190,-0.125)*u.deg)+80


plt.imshow(np.vstack([ideal_lqm]*155))


plt.plot(ideal_lqm)


plt.imshow(np.outer(ideal_lqm,vspec).transpose())


plt.imshow(hi4pi.data[:,400,:])


from astropy.modeling import models


gl = models.Gaussian1D(1.,0.,1.)


def two_d_ize (f1,f2):
    

plt.plot(gl(np.arange(-10,10)))


gb = models.Gaussian1D(1.,0.,1.)


#(the ^ is the plot where i divided by the median
#and norm='asinh' shows two distinct populations)^maybe this is something? it is saying that most the places in the galaxy where the brightness is greater than
#the median at that lattitude is associated with molecular clouds -- ie if we were to say heres a smooth
#(turns out to be lorentzian (why not gaussian?)) distribution of gas as a function of latitude / distance above
#the midplane, but oh wait its not smoot there are bits that are clumpy / have higher concentration
#then those bits are also the ones where we see molecular gas

#also when comparing with planck map -- those are also the places where polarization fraction is lowest?
#and dust emission intensity is highest (so something about magnetic fields / polarization pushing
#gas to the midplane/ lowering the scale height? and no/less magnetic field->greater concentrations
#of hi -> moelcular gas?
#also, in this image,the places where there is high density hi without signficant h2 is at high b
#so could be explained gravitationally)
#offsets @ corona, chamealeon -> drifting from the dense hi? are the magnetic fields stronger there?
#alsoooo places with dense hi and significant dust but no h2 -> co dark h2??
#alsoo - make these "idealized" functions 2d -> eg for lv, its not just a gaussian lv across all l
#it's a sinusoid, and then for first moment also need to put together the two idealiziations
#... eventually create a toy rotataing milky way? where the model intensities as a function of b, 
#and as a function of l/v (coupled) & the observed deviations from that tell us something?


#next steps:
#still work on dendrogram decomposition?(see if pruning based on boundary pixel value works)
#add in planck 353 / polarization map into the mix
#start doing "science plots" -> in the same spirit as the co/halpha paper, correlations between bright spots /
#strong emission across the different things
#look more into velocity gradient method, turbulence, etc
#flows! eg california, other clouds with strong velocity gradients, maybe related to magnetic field stuff?

#whats going on with the hole near taurus/perseus? other enclosed/bay-like features?
#other ways of extracting 3d structure? 
#other tracers of molecular gas? 
#hii regions? stars? halpha?

#also put ^ together with ~rough~ reid distances
#and! correlating dust lbd with hi/h2 lbv
#also the turbulence / velocity stuff from convo with ruth
#also, still hisa to be explored (eg is that whats going on in the hole?)


def check_Sigma(struc):
    wstruc = np.nansum((struc.get_mask()*dend.data)[struc.indices()],axis=0)
    nh_mean_struc = np.nanmean((wstruc * hi_conv * u.K * u.km/u.s)/(1.36*c.m_p))
    return nh > 400
