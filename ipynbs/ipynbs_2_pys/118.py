#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import matplotlib.pyplot as plt


import numpy as np


from astropy.table import QTable


from astropy.modeling import models,fitting
import datacube as dc
import pandas as pd
import astropy.units as u
import astropy.constants as c


cat = QTable(fits.open('post_flash_cat.fits')[1].data)


ocat = QTable(fits.open('re_expanded_cat.fits')[1].data)


plt.scatter(cat['Sigma_GMC<Msun/pc2>']*u.solMass/u.pc**2/c.m_p,cat['spec_ang_mom_GMC<cm2/s>'])

plt.xlabel(r'$\Sigma_{H_2}\ [M_\odot\ pc^{-2}]$')
plt.ylabel(r'$j_{GMC}\ [cm^2\ s^{-1}$]')

plt.loglog()

plt.show()


gmc_n = (cat['mass_GMC<Msun>']*u.solMass/(cat['radius_GMC<pc>']*u.pc)**3/c.m_p).cgs
hi_n = (cat['mass_HI<Msun>']*u.solMass/(cat['radius_HI<pc>']*u.pc)**3/c.m_p/1.36).cgs


plt.scatter(cat['radius_GMC<pc>'],cat['spec_ang_mom_GMC<cm2/s>'])
plt.xlabel(r'$R_{GMC}\ [pc]$')
plt.ylabel(r'$j_{GMC}\ [cm^2\ s^{-1}$]')


plt.fill_between(np.arange(1e-0,1e2,1e-1),
                 1e23*np.arange(1e-0,1e2,1e-1)**1.4,
                 1e23*np.arange(1e-0,1e2,1e-1)**1.7,alpha=0.3)

plt.loglog()

plt.show()


fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.)


fttr(init,np.log10(cat['mass_GMC<Msun>']),np.log10(cat['spec_ang_mom_GMC<cm2/s>']))


plt.scatter(cat['mass_GMC<Msun>'],cat['spec_ang_mom_GMC<cm2/s>'])
plt.xlabel(r'$M_{GMC}\ [M_\odot]$')
plt.ylabel(r'$j_{GMC}\ [cm^2\ s^{-1}$]')


plt.fill_between(np.arange(1e2,1e6,1e2),
                 10**(21.2)*np.arange(1e2,1e6,1e2)**0.5,
                 10**(21.2)*np.arange(1e2,1e6,1e2)**0.9,alpha=0.3)

plt.loglog()

plt.show()


fttr(init,np.log10(cat['vel_grad_HI<km/s/pc>']),np.log10(cat['radius_HI<pc>']))


fttr(init,np.log10(cat['vel_grad_GMC<km/s/pc>']),np.log10(cat['radius_GMC<pc>']))


fttr(init,np.log10(cat['radius_HI<pc>']),np.log10(cat['b_rot_HI']))


plt.scatter(cat['radius_HI<pc>'],cat['b_rot_HI'])


plt.xlabel(r'$\Omega_{HI}\ [pc\ km\ s^{-1}]$')
plt.ylabel(r'$R_{HI}\ [pc$]')


plt.loglog()

plt.show()


plt.scatter(cat['vel_grad_HI<km/s/pc>'],cat['radius_HI<pc>'])


plt.xlabel(r'$\Omega_{HI}\ [pc\ km\ s^{-1}]$')
plt.ylabel(r'$R_{HI}\ [pc$]')


plt.loglog()

plt.show()


plt.scatter(cat['vel_grad_GMC<km/s/pc>'],cat['radius_GMC<pc>'])

plt.xlabel(r'$\Omega_{GMC}\ [pc\ km\ s^{-1}]$')
plt.ylabel(r'$R_{GMC}\ [pc$]')


plt.loglog()

plt.show()


plt.scatter(gmc_n,cat['spec_ang_mom_GMC<cm2/s>'])
plt.xlabel(r'$n_{H_2}\ [cm^{-3}]$')
plt.ylabel(r'$j_{GMC}\ [cm^2\ s^{-1}$]')


plt.loglog()

plt.show()


geom_rad = np.sqrt(np.nanstd(cat['lon_GMC<deg>'],axis=1)*\
                   np.nanstd(cat['lat_GMC<deg>'],axis=1))*u.deg.to(u.rad)*\
                    cat['distance<pc>']


np.array((cat['name'],cat['radius_GMC<pc>'],geom_rad)).transpose()
#the geometric radius is almost always way bigger, because saying that the radius
#is only dependent on the area of CO emission is forgetting that there is probably
#H2 that isnt CO bright, that is still a part of the structure we care about


cat[['name','radius_GMC<pc>','vel_grad_GMC<km/s/pc>']][[18,15,13,20]]


plt.scatter(cat['name'],gmc_n)

plt.xticks(rotation='vertical')

plt.show()


(5511*u.solMass/(10.72*u.pc)**3/c.m_p).cgs


(7325*u.solMass/(9.75*u.pc)**3/c.m_p).cgs


plt.hist(gmc_n)
plt.hist(hi_n,alpha=0.5)

plt.show()


plt.hist(cat['radius_GMC<pc>'],bins=10,edgecolor='black')
plt.vlines([9.75,10.72],0,7,color='red')
plt.show()


plt.hist(np.log10(cat['mass_GMC<Msun>']),bins=10,edgecolor='black')
plt.vlines(np.log10([7325.43,5511]),0,4,color='red')
plt.show()


#the red lines are where pegasus and ursa major would fall in the distribution of cloud properties


#also one more thing to think about -> how does solid body rotation make sense next to clouds
#that are a bunch of small concentrations next to each other (ophiuchus, lupus, chamaeleon, aquila s etc)
#maybeee look at dust to see if there is a contour that encloses all of 
#them ie pointing to an underlying structure?
#and ni the case of lambda orionis, its def not a continuous structure
#makes looking at moments of inertia even more imperatve
#maybe use 3d structure for principal axis determination etc?
#or even if i don't do to ^ -> think about calculating the moment of inertia as an integral
#...will that give something too different from 0.2-0.5xmr^2 is a diff question
#but yea the sum should be p straightforward -> something like, once an axis of rotation
#is given by arctan(omegay/omegax), break up the cloud into x bins that go out perpendicularly
#from that axis, then do mr^2 for each 
#also figure out how to slice throug =h a datacube based on some arbitrarily defined line
#then do that and show those plots of position vs v to see if the line is a good fit
#in addition(? replacement?) to just the r vs v scatterplots
#bc ^ this might reveal something that the scatter plot might not? thinking of like inflows and stuff


fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.)


fttr(init,np.log10(cat['radius_GMC<pc>']),np.log10(cat['spec_ang_mom_GMC<cm2/s>']))


fttr(init,np.log10(cat['radius_HI<pc>']),np.log10(cat['spec_ang_mom_HI<cm2/s>']))


plt.scatter(cat['radius_GMC<pc>'],cat['spec_ang_mom_GMC<cm2/s>'])
plt.scatter(cat['radius_HI<pc>'],cat['spec_ang_mom_HI<cm2/s>'])

plt.fill_between(np.arange(1e-0,1e2,1e-1),
                 1e23*np.arange(1e-0,1e2,1e-1)**1.4,
                 1e23*np.arange(1e-0,1e2,1e-1)**1.7,alpha=0.3)
plt.loglog()
plt.xlim(1e-0,1e2)


ora=fits.open('h2_clouds.fits')[15].data


plt.imshow(np.nansum(ora,axis=0),cmap=dc.cmap_min0('viridis'),origin='lower',vmax=60)
plt.colorbar()

plt.scatter(np.arange(0,96),-0.25*np.arange(0,96)+25,color='red',s=10,alpha=0.5)
plt.scatter(np.arange(0,96),np.around(-0.25*np.arange(0,96)+25),color='orange',s=10,alpha=0.5)

# plt.xlim(35,50)
# plt.ylim(5,30)


lonscoords = np.meshgrid(np.arange(0,96),np.arange(0,32))[0]


latscoords = np.meshgrid(np.arange(0,96),np.arange(0,32))[1]


oradf = pd.DataFrame(np.array([np.ravel(lonscoords),
                       np.ravel(latscoords),
                       np.ravel(np.nansum(ora,axis=0))]).transpose(),columns=['lon','lat','im'])


plt.scatter(oradf['lon'],oradf['lat'],c=oradf['im'],cmap=dc.cmap_min0('viridis'))


plt.imshow(latscoords)
plt.colorbar()


plt.imshow(lonscoords)
plt.colorbar()


plt.imshow(np.around(25+lonscoords*-0.25))
plt.colorbar()


#hmm okay 1) decide if taking position/velocity slices through not the principal axes of the array
#is useful in the long term / worht the effort
#then b) figure out how to do it (it feels like there should be something otut in the world for this)


ora.shape


np.array((np.arange(20,70),np.around(-0.25*np.arange(20,70)+25))).transpose()


np.array((np.arange(20,70),-0.25*np.arange(20,70)+25)).transpose()


plt.contourf(np.nansum(ora,axis=1),vmin=10,cmap=dc.cmap_min0('viridis'))
plt.colorbar()


orl=fits.open('h2_clouds.fits')[17].data


plt.contourf(np.nansum(orl,axis=0),vmin=2,cmap=dc.cmap_min0('viridis'))
plt.colorbar()


plt.contourf(np.nansum(orl,axis=1),vmin=5,cmap=dc.cmap_min0('viridis'))
plt.colorbar()


plt.contourf(np.nansum(ora,axis=0),cmap=dc.cmap_min0('viridis'),levels=np.arange(5,100,2))
plt.colorbar()


plt.contourf(np.nansum(ora,axis=1),vmin=5,cmap=dc.cmap_min0('viridis'),levels=np.arange(5,100,2))
plt.colorbar()


plt.contourf(np.nansum(ora,axis=2).transpose(),vmin=5,cmap=dc.cmap_min0('viridis'),levels=np.arange(15,300,20))
plt.colorbar()


plt.contourf(np.nansum(orl,axis=2),vmin=5,cmap=dc.cmap_min0('viridis'))
plt.colorbar()


help(ax.scatter3D)


orl.shape


ora


orlcords = np.meshgrid(np.arange(88),np.arange(39),np.arange(108))


oracords = np.meshgrid(np.arange(32),np.arange(27),np.arange(96))


ora.shape


np.array(orlcords).shape


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')

ax.scatter3D(oracords[2],oracords[1],oracords[0],c=ora,norm='log',s=ora,cmap=dc.cmap_min0('viridis'))

ax.view_init(azim=200,elev=30)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')

ax.scatter3D(oracords[2],oracords[1],oracords[0],c=ora,norm='log',s=ora,cmap=dc.cmap_min0('viridis'))

ax.view_init(azim=200,elev=30)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')

ax.scatter3D(orlcords[2],orlcords[1],orlcords[0],c=orl,s=orl,cmap=dc.cmap_min0('viridis'))

ax.view_init(azim=90,elev=20)


dust = fits.open('mean_and_std_lbd.fits')[1].data


dlb = np.nansum(dust,axis=0)


plt.contourf(dlb[290:324,1370:1480])


dorld = np.nansum(dust[:,290:324,1370:1480],axis=1)


dorbd = np.nansum(dust[:,290:324,1370:1480],axis=2)


plt.contourf(dorbd,levels=np.arange(1e-2,7e-1,5e-2))
plt.colorbar()
plt.ylim(250,500)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dorbd,levels=np.arange(1e-2,7e-1,5e-2))
plt.colorbar()
plt.ylim(250,500)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dorld,levels=np.arange(1e-2,7e-1,5e-2))
plt.colorbar()
plt.ylim(250,500)
plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dorld,levels=np.arange(1e-2,7e-1,5e-2))
plt.colorbar()
plt.ylim(250,500)
plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dlb[335:375,10:120])


dcalld = np.nansum(dust[:,335:375,10:120],axis=1)


dcalbd = np.nansum(dust[:,335:375,10:120],axis=2)


plt.contourf(dcalbd.transpose(),levels=np.arange(1e-2,5e-1,3e-2))
plt.colorbar()
plt.xlim(250,600)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dcalbd,levels=np.arange(1e-2,5e-1,3e-2))
plt.colorbar()
plt.ylim(250,600)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dcalld,levels=np.arange(1e-2,3e-1,2e-2))
plt.colorbar()
plt.ylim(250,600)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dcalld,levels=np.arange(1e-2,3e-1,2e-2))
plt.colorbar()
plt.ylim(250,600)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dlb[260:340,70:120])


dperld = np.nansum(dust[:,260:340,70:120],axis=1)


dperbd = np.nansum(dust[:,260:340,70:120],axis=2)


plt.contourf(dperld,levels=np.arange(1e-2,9e-1,5e-2))
plt.colorbar()
plt.ylim(150,350)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dperbd.transpose(),levels=np.arange(1e-2,9e-1,5e-2))
plt.colorbar()
plt.xlim(150,350)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dlb[250:350,700:850],levels=np.arange(1e-2,9e-1,5e-2))


dcrald = np.nansum(dust[:,250:350,700:850],axis=1)


dcrabd = np.nansum(dust[:,250:350,700:850],axis=2)


plt.contourf(dcrabd.transpose(),levels=np.arange(3e-2,3e-1,2e-2))
plt.colorbar()
plt.xlim(100,200)
# plt.xlim(20,)
# plt.yscale('log')


plt.contourf(dcrald,levels=np.arange(3e-2,6e-1,3e-2))
plt.colorbar()
plt.ylim(100,200)
# plt.xlim(20,)
# plt.yscale('log')


cat[5]


plt.contourf(dcalld,levels=np.arange(1e-2,6e-1,3e-2))
plt.hlines(120,0,70,color='red')


plt.contourf(dorld,levels=np.arange(1e-2,6e-1,3e-2))
plt.hlines(120,0,70,color='red')


plt.contourf(dperld,levels=np.arange(3e-2,6e-1,3e-2))
plt.hlines(120,0,70,color='red')


plt.contourf(dcrald,levels=np.arange(3e-2,6e-1,3e-2))
plt.hlines(120,0,70,color='red')


nearby = np.nansum(dust[:120,:,:],axis=0)


plt.figure(figsize=(20,15))
plt.contourf(dlb,levels=np.arange(5e-2,1,5e-2),cmap='binary_r')
plt.contour(nearby,levels=np.arange(2e-2,5e-1,3e-2))


nearby2 = np.nansum(dust[120:300,:,:],axis=0)


plt.figure(figsize=(20,15))
plt.contourf(dlb,levels=np.arange(5e-2,1,5e-2),cmap='binary_r')
plt.contour(nearby2,levels=np.arange(2e-2,5e-1,3e-2))


plt.plot(np.nanmean(dust,axis=(1,2)))


far = np.nansum(dust[700:,:,:],axis=0)


plt.figure(figsize=(20,15))
plt.contourf(dlb,levels=np.arange(5e-2,1,5e-2),cmap='binary_r')
plt.contour(far,levels=np.arange(2e-2,5e-1,3e-2))


far2 = np.nansum(dust[900:,:,:],axis=0)


plt.figure(figsize=(20,15))
plt.contourf(dlb,levels=np.arange(5e-2,1,5e-2),cmap='binary_r')
plt.contour(far2,levels=np.arange(2e-2,5e-1,3e-2))


dbd = np.nansum(dust,axis=2)


plt.contourf(dbd.transpose(),levels=np.arange(3e-2,1,1e-2),norm='log',extend='both')


plt.contourf(dlb,levels=np.arange(5e-2,7e-1,5e-2))
plt.xlim(1010,1100)
plt.ylim(280,380)


dchmld = np.nansum(dust[:,280:380,1010:1100],axis=1)


dchmbd = np.nansum(dust[:,280:380,1010:1100],axis=2)


plt.contourf(dchmld,levels=np.arange(1e-1,1,8e-2))
# plt.xlim(35,120)
plt.ylim(60,165)


plt.contourf(dchmbd.transpose(),levels=np.arange(1e-1,1,8e-2))
# plt.xlim(35,120)
plt.xlim(60,165)


plt.contourf(dchmld,levels=np.arange(1e-1,1,8e-2))
plt.xlim(35,120)
plt.ylim(60,165)


plt.contourf(dchmbd.transpose(),levels=np.arange(1e-1,1,8e-2))
# plt.xlim(35,120)
plt.xlim(60,165)


cat['name','vel_grad_GMC<km/s/pc>']


plt.contourf(dlb,levels=np.arange(5e-2,9e-1,5e-2))
plt.xlim(1480,1550)
plt.ylim(310,370)


dorlld = np.nansum(dust[:,310:370,1480:1550],axis=1)
dorlbd = np.nansum(dust[:,310:370,1480:1550],axis=2)


plt.contourf(dorlld,levels=np.arange(5e-2,9e-1,5e-2))
plt.ylim(300,400)


plt.contourf(dorlbd.transpose(),levels=np.arange(5e-2,9e-1,5e-2))
plt.xlim(300,400)
