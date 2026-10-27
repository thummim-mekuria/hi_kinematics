#!/usr/bin/env python
# coding: utf-8


from astrodendro import Dendrogram, ppv_catalog
import astropy.units as u
import astropy.constants as c
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datacube as dc
from skimage.segmentation import find_boundaries
from astropy.io import fits


hi4pi,dht = fits.open('both_same_grid_720d.fits')


hi4pi_dend = Dendrogram.load_from('decomposed_HI4PI.fits')


import multiprocessing


# from skimage.filters import butterworth 
# plt.figure(figsize=(20,10))
# plt.imshow(butterworth(whi),origin='lower',vmin=0,vmax=2e2,extent=[190,-190,-50,50])
# plt.contour(lons,lats,butterworth(whi),levels=[50],colors='orange',linewidths=0.9)

# plt.contour(lons,lats,wco,levels=[1],colors='white',linewidths=0.9)
# plt.xlim(190,50)
# plt.colorbar()
##may come in handy at some point -- filters for low or high frequency features
##basically what i was trying to do with the whole hi4pi_calib_flat etc 


whi = np.nansum(hi4pi.data,axis=0)
wco = np.nansum(dht.data,axis=0)


lons = np.arange(5760)*-0.125 + 360
lats = np.arange(800)*0.125 - 50


edenhoffer = fits.open('mean_and_std_lbd.fits')
dust_mean = edenhoffer[1]
dlons = (np.arange(1572)-787)*-0.23
dlats = (np.arange(786)-394)*0.23
dvels = (np.arange(1182)-1) + 69


wd = np.nansum(dust_mean.data,axis=0)


bspec = np.nanmean(whi,axis=1)
idealized_bspec = (np.concatenate([bspec]*5760,axis=0)).reshape(-1,800).transpose()


plt.figure(figsize=(20,14))
plt.imshow(whi/idealized_bspec,extent=[360,-360,-50,50],vmin=0,vmax=2,origin='lower',)
# plt.colorbar()

plt.contour(lons,lats,whi/idealized_bspec,levels=[1.5,2],colors='orange',linewidths=0.9)

plt.contour(lons,lats,wco,levels=[1],colors='black',linewidths=0.5)
# plt.contour(dlons,dlats,wd,levels=[0.2],colors='red',linewidths=0.7)

plt.xlim(180,-180)
plt.ylim(-50,50)


plt.figure(figsize=(20,14))
plt.imshow(whi-idealized_bspec,extent=[360,-360,-50,50],norm='asinh',
           cmap=dc.cmap_min0('viridis'),origin='lower')
# plt.colorbar()

plt.contour(lons,lats,wco,levels=[1],colors='orange',linewidths=0.9)
plt.contour(dlons,dlats,wd,levels=[0.2],colors='red',linewidths=0.7)

plt.xlim(180,-180)
plt.ylim(-50,50)


import fun as fns
mom1_hi = fns.get_mom1(hi4pi.data,[0,155])


idealized_vspec.shape


idealized_vspec2 = 


plt.figure(figsize=(20,14))
plt.imshow(mom1_hi/idealized_bspec,extent=[360,-360,-50,50],vmin=-0.05,vmax=0.05,origin='lower')
# plt.colorbar()

plt.contour(lons,lats,whi/idealized_bspec,levels=[1.5],colors='orange',linewidths=0.9)

plt.contour(lons,lats,wco,levels=[1],colors='black',linewidths=0.5)
# plt.contour(dlons,dlats,wd,levels=[0.2],colors='red',linewidths=0.7)

plt.xlim(180,-180)
plt.ylim(-50,50)


from astropy.table import QTable


hlpy.read_map('HFI_SkyMap_353_2048_R3.01_full.fits',field=0)[0]


hlpy.read_map('HFI_SkyMap_353_2048_R3.01_full.fits',field=1)[0]


np.nanmedian([0.067651086,0.0018415415,6.973286e-5,1512,8.9258965e-9,\
-1.8544616e-9,2.0536777e-9,5.6715805e-8,2.4524027e-10,9.193198e-8])


np.sqrt(50331648)


hlpy.mollzoom(hlpy.read_map('HFI_SkyMap_353_2048_R3.01_full.fits'))


fits.open('HFI_SkyMap_353_2048_R3.01_full.fits')


fits.open('HFI_SkyMap_353_2048_R3.01_full.fits')[1].header


QTable(fits.open('HFI_SkyMap_353_2048_R3.01_full.fits')[1].data)


import healpy as hlpy


from astropy_healpix import HEALPix as hlpx


fits.open('HFI_SkyMap_353_2048_R3.01_full.fits')[1].header


planck_hpx = hlpx(nside=2048,order='nested')


planck_dust = QTable(fits.open('HFI_SkyMap_353_2048_R3.01_full.fits')[1].data)


pol_frac = np.array(np.sqrt(planck_dust['Q_STOKES']**2 + planck_dust['U_STOKES']**2)/planck_dust['I_STOKES'])


planck_coords = np.array(planck_hpx.healpix_to_lonlat(np.arange(len(planck_dust))))


pd.DataFrame(np.array((planck_coords[0]*u.rad.to(u.deg),
                       planck_coords[1]*u.rad.to(u.deg),
                       pol_frac)).transpose())


plt.figure(figsize=(20,10))
plt.scatter((planck_coords[0]*u.rad.to(u.deg))[::100],
            (planck_coords[1]*u.rad.to(u.deg))[::100],
            c=100*pol_frac[::100],s=4,vmin=0,vmax=20)
plt.scatter((-360+planck_coords[0]*u.rad.to(u.deg))[::100],
            (planck_coords[1]*u.rad.to(u.deg))[::100],
            c=100*pol_frac[::100],s=4,vmin=0,vmax=20)

plt.colorbar()

plt.contour(lons,lats,wco,levels=[1],colors='white',linewidths=0.5)

plt.xlim(180,-180)
plt.show()


help(hlpy.mollview)


pol_frac_w = hlpy.cartview(pol_frac,min=0,max=0.5,nest=True,return_projected_map=True,xsize=100000)


np.min(planck_coords[0]*u.rad.to(u.deg)),np.max(planck_coords[0]*u.rad.to(u.deg))


pol_frac_w.shape


plt.imshow(pol_frac_w,norm='asinh',vmin=0,vmax=0.5)
plt.colorbar()


10**(0.1)


len(np.arange(0,360,0.028))


1.7177*u.arcmin.to(u.degree)


planck_hpx.pixel_resolution


320000/50331648


len(planck_coords[0])


plt.figure(figsize=(20,10))
# plt.scatter((planck_coords[0]*u.rad.to(u.deg))[::100],
#             (planck_coords[1]*u.rad.to(u.deg))[::100],
#             c=(planck_dust['I_STOKES']*pol_frac)[::100],s=4,vmin=1e-8,vmax=5e-4)
# plt.scatter((-360+planck_coords[0]*u.rad.to(u.deg))[::100],
#             (planck_coords[1]*u.rad.to(u.deg))[::100],
#             c=(planck_dust['I_STOKES']*pol_frac)[::100],s=4,vmin=1e-8,vmax=5e-4)

plt.imshow(pol_frac_w,extent=[180,-180,-90,90],vmax=0.5)

plt.colorbar()

plt.contour(lons,lats,wco,levels=[1],colors='white',linewidths=0.5)

plt.xlim(180,-180)
plt.show()


T_Snu = 2*c.k_B*(353*u.GHz)**2 / c.c**2


planck_hpx.npix


50331648/8191


hlpy.projector.CartesianProj.projmap(hlpy.projector.CartesianProj.vec2xy(),map=pol_frac)


import reproject as rpoj


from astropy.wcs import WCS


ssss = WCS()


len(np.arange(-90,90,0.125))


hi4pi_w_wcs = WCS({'CTYPE1' : 'GLON-CAR',
     'CTYPE2':'GLAT-CAR',
'CRVAL1' : 180., 
   'CRVAL2' :  -90.0  ,
'CRPIX1' : 0.0, 
     'CRPIX2' :  0.0  ,
'CDELT1' : -0.125,  
  'CDELT2': 0.125  ,
'NAXIS1' : 2880,
     'NAXIS2' : 1440})


hmm1,hmm2 = rpoj.reproject_from_healpix('HFI_SkyMap_353_2048_R3.01_full.fits',hi4pi_w_wcs)


plt.imshow(hmm1,norm='asinh',vmax=0.1)
plt.colorbar()
plt.xlim(0,1600)
plt.ylim(0,800)


np.sqrt(50331648/12)


fits.open('HFI_SkyMap_353_2048_R3.01_full.fits')[1].header


plt.figure(figsize=(20,10))
plt.scatter((planck_coords[0]*u.rad.to(u.deg))[::100],
            (planck_coords[1]*u.rad.to(u.deg))[::100],
            c=planck_dust['Q_STOKES'][::100],s=4,vmin=-5e-4,vmax=5e-4)
plt.scatter((-360+planck_coords[0]*u.rad.to(u.deg))[::100],
            (planck_coords[1]*u.rad.to(u.deg))[::100],
            c=planck_dust['Q_STOKES'][::100],s=4,vmin=-5e-4,vmax=5e-4)

plt.colorbar()

plt.contour(lons,lats,wco,levels=[1],colors='white',linewidths=0.5)

plt.xlim(180,-180)
plt.show()


plt.figure(figsize=(20,10))
plt.scatter((planck_coords[0]*u.rad.to(u.deg))[::100],
            (planck_coords[1]*u.rad.to(u.deg))[::100],
            c=planck_dust['U_STOKES'][::100],s=4,vmin=-5e-4,vmax=5e-4)
plt.scatter((-360+planck_coords[0]*u.rad.to(u.deg))[::100],
            (planck_coords[1]*u.rad.to(u.deg))[::100],
            c=planck_dust['U_STOKES'][::100],s=4,vmin=-5e-4,vmax=5e-4)

plt.colorbar()

plt.contour(lons,lats,wco,levels=[1],colors='white',linewidths=0.5)

plt.xlim(180,-180)
plt.show()


#^maybe this is something? it is saying that most the places in the galaxy where the brightness is greater than
#the median at that lattitude is associated with moelcular clouds -- ie if we were to say heres a smooth
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


lvhi = np.nansum(hi4pi.data,axis=1)


vspec = np.nanmean(lvhi,axis=1)
idealized_vspec = (np.concatenate([vspec]*5760,axis=0)).reshape(-1,155).transpose()


vels = np.arange(155)*1.3-100


lvco = np.nansum(dht.data,axis=1)


plt.figure(figsize=(20,10))
plt.imshow(idealized_vspec,origin='lower',extent=[360,-360,-100,100],norm='asinh')
# plt.contour(lons,vels,lvhi,levels=[2e3,5e3,1e4],colors='orange',linewidths=0.9)
#
plt.contour(lons,vels,lvco,levels=[1],colors='white',linewidths=0.9)
plt.xlim(180,-180)
# plt.colorbar()


plt.figure(figsize=(20,10))
plt.imshow(lvhi,origin='lower',extent=[360,-360,-100,100],norm='asinh')
# plt.contour(lons,vels,lvhi,levels=[2e3,5e3,1e4],colors='orange',linewidths=0.9)
#
plt.contour(lons,vels,lvco,levels=[1],colors='white',linewidths=0.9)
plt.xlim(180,-180)
# plt.colorbar()


plt.figure(figsize=(20,10))
plt.imshow(lvhi/idealized_vspec,origin='lower',extent=[360,-360,-100,100],norm='asinh')
plt.contour(lons,vels,lvhi/idealized_vspec,levels=[2.5],colors='orange',linewidths=0.9)
#
# plt.contour(lons,vels,lvco,levels=[1],colors='white',linewidths=0.9)
plt.xlim(180,-180)
# plt.colorbar()


plt.figure(figsize=(20,10))
plt.imshow(lvhi-idealized_vspec,origin='lower',extent=[360,-360,-100,100],norm='asinh')
# plt.contour(lons,vels,lvhi,levels=[2e3,5e3,1e4],colors='orange',linewidths=0.9)
#
plt.contour(lons,vels,lvco,levels=[1],colors='white',linewidths=0.9)
plt.xlim(180,-180)
# plt.colorbar()


co_dend = Dendrogram.load_from('../not_rn_3/COGAL_dendrogram.fits')


hicat = pd.read_csv('hi4pi_dendleaves_5solMass.csv')


cocat = pd.read_csv('co_dendleaves.csv')


cocoords = [np.nanmean(lf.indices(),axis=1)*[1.3,0.125,-0.125] + [-50,-30,180] for lf in coleavs]


cocat = pd.DataFrame(cocoords)


cocat['idx'] = [lf.idx for lf in coleavs]


cocat = cocat.sort_values('idx')


cocat.columns=['vel<km/s>','lat<deg>','lon<deg>','idx']


cocat = cocat[['idx','lon<deg>','lat<deg>','vel<km/s>']]


plt.figure(figsize=(20,10))
plt.scatter(hicat['lon<deg>'],hicat['lat<deg>'],alpha=0.3,s=300)
plt.scatter(cocat['lon<deg>'],cocat['lat<deg>'],alpha=0.3,s=300)


co_masks = np.array([co_dend[i].get_mask() for i in cocat['idx']])


cocat.to_csv('co_dendleaves.csv')


# hi_masks = np.array([hi4pi_dend[i].get_mask() for i in hicat['idx']])


hicat['area<deg2>'] = [(0.125**2)*len(np.unique(hi4pi_dend[i].indices()[1:],\
                                                axis=1)[0]) for i in hicat['idx']]


cocat['area<deg2>'] =[(0.125**2)*len(np.unique(co_dend[i].indices()[1:],axis=1)[0]) for i in cocat['idx']]


cocat['radius<deg>'] = np.sqrt(cocat['area<deg2>'].to_numpy()/np.pi)


hicat['radius<deg>'] = np.sqrt(hicat['area<deg2>'].to_numpy()/np.pi)


cocat.to_csv('co_dendleaves.csv',index=False)
hicat.to_csv('hi4pi_dendleaves_5solMass.csv',index=False)


plt.hist(cocat['radius<deg>'],alpha=0.5,bins=20,density=True)
plt.hist(hicat['radius<deg>'],alpha=0.5,bins=20,density=True)


#go through & look at / verify matches


matches=[(i,j)
 for i in range(len(hicat['idx']))
    for j in range(len(cocat['idx']))
if np.sqrt(
    (cocat['lon<deg>'][j] - hicat['lon<deg>'][i])**2
    + (cocat['lat<deg>'][j] - hicat['lat<deg>'][i])**2
    + (cocat['vel<km/s>'][j] - hicat['vel<km/s>'][i])**2


) < 5]


matches_dict={}


for i,pr in enumerate(matches):
    if pr[0] not in matches_dict.keys():
        matches_dict[pr[0]] = [pr[1]]
    else:
        matches_dict[pr[0]].append(pr[1])


for k in matches_dict:
    for v in matches_dict[k]:
        plt.imshow(np.nansum(co_dend[cocat['idx'][v]].get_mask()*co_dend.data,axis=0),
                               cmap=dc.cmap_min0('viridis'),origin='lower',
                     extent=[180,-180,-30,30])

  
    plt.contour(np.arange(190,-190,-0.125),np.arange(-50,50,0.125),
                 np.nansum(hi4pi_dend[hicat['idx'][k]].get_mask()*hi4pi_dend.data,axis=0),
                  levels=3,linewidths=0.7)

    plt.xlim(hicat['lon<deg>'][k]+20,hicat['lon<deg>'][k]-20),
    plt.ylim(hicat['lat<deg>'][k]-20,hicat['lat<deg>'][k]+20),
 
    plt.show()
                     
 
matches


p = hi4pi_dend.plotter()


import astrodendro.pruning as dend_prune


from astrodendro import periodic_neighbours


fig,ax = plt.subplots(1,1,figsize=(30,10))
p.plot_tree(ax)
# plt.yscale('log')
# plt.xlim(-10,100)


plt.hist(hicat['radius<deg>'],bins=20)
plt.hist(cocat['radius<deg>'],bins=20,alpha=0.7)


plt.show()


plt.figure(figsize=(20,10))
plt.scatter(cocat['lon<deg>'],cocat['lat<deg>'],s=50*cocat['radius<deg>'],alpha=0.5)
plt.scatter(hicat['lon<deg>'],hicat['lat<deg>'],s=50*hicat['radius<deg>'],alpha=0.5)

plt.scatter(cm['lon<deg>'],cm['lat<deg>'],s=300*cm['radius<deg>'],alpha=0.3)
plt.scatter(hm['lon<deg>'],hm['lat<deg>'],s=300*hm['radius<deg>'],alpha=0.3)


cm = cocat.iloc[np.unique(np.array(matches).transpose()[1])]


hm = hicat.iloc[list(matches_dict.keys())]


plt.imshow(np.nansum(hi4pi_dend[344].get_mask(),axis=0).astype('bool')
           ,cmap=dc.cmap_min0('viridis'),origin='lower')
plt.ylim(380,400)
plt.xlim(580,640)


co_dend.params


whi = np.nansum(hi4pi_dend.data,axis=0)


plt.imshow(find_boundaries(np.nansum(hi4pi_dend[344].get_mask(),axis=0).astype('bool')),
           cmap=dc.cmap_min0('viridis'),origin='lower')
plt.ylim(380,400)
plt.xlim(580,640)


plt.scatter(np.array(hi4pi_dend[344].indices())[2],np.array(hi4pi_dend[344].indices())[1])
plt.scatter(np.where(find_boundaries(np.nansum(hi4pi_dend[344].get_mask(),axis=0).astype('bool')) == True)[1],
            np.where(find_boundaries(np.nansum(hi4pi_dend[344].get_mask(),axis=0).astype('bool')) == True)[0])


for i in hicat['idx']:
    this=hi4pi_dend[i]
#     print(this)
    m2D = np.nansum(this.get_mask(),axis=0)
#     print(m2D)
    brdr_px = np.array(np.where(find_boundaries(m2D.astype('bool'))==True)).transpose()
#     print(brdr_px)
    brdr_vals = np.array([(m2D*whi)[cord[0],cord[1]] for cord in brdr_px])
    print(i,np.all(brdr_vals[brdr_vals>0]>400),np.nanmean(brdr_vals[brdr_vals>0]))


hmm=[i for i in range(len(hi4pi_dend)) if hi4pi_dend[i].level <=2]


hmmm= [i for i in range(len(hi4pi_dend))
                                 if hi4pi_dend[i].level >=2 if hi4pi_dend[i].level <=5]


hmmm1e22e2= [i for i in range(len(hi4pi_dend))
                                 if hi4pi_dend[i].level >=250]


plt.hist([hi4pi_dend[i].level for i in range(len(hi4pi_dend))])


len(hmmm1e22e2)


[(plt.imshow(np.nansum(hi4pi_dend.data*hi4pi_dend[idx].get_mask(),axis=0),
             origin='lower',cmap=dc.cmap_min0('viridis')),
  plt.show(),
  print(idx)


)for idx in hmmm1e22e2]


fig,ax = plt.subplots(1,1,figsize=(20,10))
p.plot_tree(ax)


np.array(hi4pi_dend[344].indices())[1:]


matches_dict


#next steps - check if all the co clouds have been assigned to an hi
#whats going on with oriona a+b being a leaf, but then b again later?
#plot the hi dend tree and see if mose of the leaves are at level 1 (ie they are going
#all thew ay down to the minimum as leaves, rather than joining up with other structures
#in whic case, maybe lower min value & re run?)
#also check that i am selecting the right leaf ids
#and what the max separation is allowed to be between the mean positions of 
#the co & hi leaf


plt.figure(figsize=(20,10))
plt.scatter(cocat['lon<deg>'],cocat['lat<deg>'],s=300*cocat['radius<deg>'],alpha=0.2)
plt.scatter(hicat['lon<deg>'],hicat['lat<deg>'],s=300*hicat['radius<deg>'],alpha=0.2)


hicat.to_csv('hi4pi_dendleaves_5solMass.csv')


plt.figure(figsize=(20,10))
plt.scatter(cocat['lon<deg>'],cocat['vel<km/s>'],s=300*cocat['radius<deg>'],alpha=0.3)
plt.scatter(hicat['lon<deg>'],hicat['vel<km/s>'],s=300*hicat['radius<deg>'],alpha=0.3)


plt.figure(figsize=(20,10))
plt.scatter(cocat['lat<deg>'],cocat['vel<km/s>'],s=300*cocat['radius<deg>'],alpha=0.3)
plt.scatter(hicat['lat<deg>'],hicat['vel<km/s>'],s=300*hicat['radius<deg>'],alpha=0.3)


hlons = np.arange(190,-190,-0.125)
hlats = np.arange(-50,50,0.125)
hvels = np.arange(-100,101,1.3)


clons = np.arange(2881)*-0.125 + 180
clats = np.arange(481)*0.125 - 30
cvels = np.arange(76)*1.3 - 49


comap = np.nansum(co_masks,axis=(0,1))


hi_masks_30 = np.array([hi4pi_dend[i].get_mask() 
                        for i in hicat['idx'][:30]])


hi_masks_60 = np.array([hi4pi_dend[i].get_mask() 
                        for i in hicat['idx'][30:60]])


hi_masks_90 = np.array([hi4pi_dend[i].get_mask() 
                        for i in hicat['idx'][60:90]])


hi_masks_120 = np.array([hi4pi_dend[i].get_mask() 
                        for i in hicat['idx'][90:120]])


hi_masks_150 = np.array([hi4pi_dend[i].get_mask()
                         for i in hicat['idx'][120:150]])


hi_masks_180 = np.array([hi4pi_dend[i].get_mask() 
                         for i in hicat['idx'][150:180]])


hi_masks_210 = np.array([hi4pi_dend[i].get_mask() 
                         for i in hicat['idx'][180:210]])


hi_masks_240 = np.array([hi4pi_dend[i].get_mask() 
                         for i in hicat['idx'][210:240]])


hi_masks_260 = np.array([hi4pi_dend[i].get_mask() 
                          for i in hicat['idx'][240:260]])


hi_masks_last = np.array([hi4pi_dend[i].get_mask() 
                          for i in hicat['idx'][260:272]])


plt.figure(figsize=(20,10))
plt.imshow(comap,extent=[180,-180,-30,30],cmap=dc.cmap_min0('viridis'),origin='lower',norm='asinh')
plt.contour(hlons,hlats,np.nansum(hi_masks_30,axis=(0,1)),
           levels=[10,50,100],colors='red')

# plt.xlim(1500,)


co_dend.params


plt.figure(figsize=(20,10))
plt.imshow(comap,cmap=dc.cmap_min0('viridis'),origin='lower',norm='asinh')
plt.xlim(0,1500)


#make ppv catalog for this^ and see how it compares to the lbv coords i extracted (in hi_dend.csv or sth smlr)


#then get the co dendrogram + catalogs & see if i can associate hi & co clouds (+also find a control group?)
