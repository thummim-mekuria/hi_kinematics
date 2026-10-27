#!/usr/bin/env python
# coding: utf-8

# in discussion, say: we left out 
# - serpens bc its too blended to the galactc plane
# - cepheus bc we know there stuff p next to each other that we know for sure is separated by 600pc, but they are hard to separate in velocity, so will cause issues with hi
# - pegasus & ursa bc they are outside the latitude range of the composite survey, theyre in the north but that has diff spatial (and velocity,but thats less impt?) resolution
# - we also used the big dame survey (instead of all the individual surveys with higher res) because the hi4pi resolution is 1.3km/s


import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import fitting,models

import numpy as np
import pandas as pd

import datacube as dc
import fun as fns
import interpolate as ip


hi4pi,dht = fits.open('both_same_grid_720d.fits')


cat = QTable(fits.open('post_flash_cat.fits')[1].data)


c0 = np.nansum(dht.data,axis=0)*1.3*1.97e20


north = fits.open('DT22_masked.fits')[0].data


north = north.transpose().swapaxes(1,2)


north_ip = ip.interpolate_vlb(north)


north0 = np.nansum(north_ip,axis=0)*1.3*1.97e20


nvls = np.arange(146)*0.65 - 47.


nlns = np.arange(1441)*0.25 - 60


nlts = (np.arange(677)-320)*0.25


plt.contourf(np.log10(north0),
             levels=np.arange(19,22,0.2),
             extent=[-60,300,-80,89])

plt.ylim(-50,-20)
plt.xlim(100,80)


plt.contourf(np.log10(north0),
             levels=np.arange(19,22,0.2),
             extent=[-60,300,-80,89])

# plt.ylim(20,50)
# plt.xlim(160,140)

plt.ylim(30,50)
plt.xlim(160,140)


plt.contourf(np.log10(c0),levels=np.arange(19,22,0.2),extent=[360,-360,-50,50])
plt.xlim(180,-180)

plt.xlim(130,70)
plt.ylim(-10,30)


(np.array((130,70))-360)/-0.125


(np.array((0,20))+50)/0.125


cph = dht.data[:,400:560,1840:2320]


plt.plot(np.nanmean(cph,axis=(1,2)))


for i in range(50,100):
    plt.contourf(cph[i],cmap=dc.cmap_min0('viridis'),
                 extent=[130,70,0,20],levels=np.arange(1e-2,6e-1,3e-2))
    plt.hlines([8,17],130,70)
    plt.vlines(100,5,30)
    
    plt.colorbar()
    plt.ylim(5,)
    plt.xlim(115,95)
    plt.show()


rcra = fits.open('DHT03_RCrA_mom.fits')


rcra[0].header


rcra = rcra[0].data.transpose().swapaxes(1,2)


rca0 = np.nansum(rcra.data,axis=0)*0.65*1.97e20


plt.contourf(np.log10(rca0),levels=np.arange(19,23,0.2))


cat[5]


(fns.mass_n_area(rcra*u.K,cat['distance<pc>'][5]*u.pc,'CO',dv=0.26,ds=0.25)[0]).to(u.solMass)


(fns.mass_n_area(rcra*u.K,cat['distance<pc>'][5]*u.pc,'CO',dv=0.26,ds=0.25)[1])


2012/1722


84.65/(np.pi*4.09**2)


aqs = fits.open('DHT06_Aquila-South_mom.fits')


aqs[0].header


aqs = aqs[0].data.transpose().swapaxes(1,2)


aqs0 = np.nansum(aqs,axis=0)*0.65*1.97e20


plt.contourf(np.log10(aqs0),levels=np.arange(19,22,0.25))


(fns.mass_n_area(aqs*u.K,cat['distance<pc>'][0]*u.pc,'CO',dv=0.65,ds=0.25)[0]).to(u.solMass)


cat[0]


1152/576


cam = fits.open('DHT19_Camelopardalis_bw_mom.fits')[0]


cam.header


cam = cam.data.transpose().swapaxes(1,2)


cam0 = np.nansum(cam,axis=0)*0.65*1.97e20


plt.contourf(np.log10(cam0),levels=np.arange(19,22,0.25))


np.nanmean((0.125,0.142))


fns.mass_n_area(cam*u.K,cat['distance<pc>'][2]*u.pc,'CO',dv=0.65,ds=0.25)[0].to(u.solMass)


3876/442


fns.mass_n_area(cam[:,0:35,30:]*u.K,cat['distance<pc>'][2]*u.pc,'CO',dv=0.65,ds=0.25)[0].to(u.solMass)


2698/442


cat[2]


for i in range(22):
    plt.scatter(cat['mass_GMC<Msun>'][i],cat['vel_dispersion_GMC<km/s>'][i])
    plt.annotate(cat['name'][i],
                 (cat['mass_GMC<Msun>'][i],cat['vel_dispersion_GMC<km/s>'][i]))
    
    
plt.xscale('log')


lac = fits.open('DHT12_Lacerta_mom.fits')[0]


lac.header


lac = lac.data.swapaxes(0,1).transpose()


lac0 = np.nansum(lac,0)*0.65*1.97e20


plt.contourf(np.log10(lac0),levels=np.arange(19,22,0.25))


fns.mass_n_area(lac*u.K,cat['distance<pc>'][9]*u.pc,'CO',dv=0.65)[0].to(u.solMass)


12613/3775


fns.mass_n_area(lac[:,30:,70:]*u.K,cat['distance<pc>'][9]*u.pc,'CO')[0].to(u.solMass)


10679/3775


len(lac0[lac0>4e20]),len(lac0[lac0>1e21]),lac0.size


len(cam0[cam0>4e20]),len(cam0[cam0>1e21]),cam0.size


#damn i was also mebedeling them by having a threshold of 1e21 to count towards mass
#esp cos these are the little guys, i was excluding proportionally way more of the material
#compared to doing this for a big cloud


#but still, look at what changing just the resolution causes 
#both spatially & velocity for mass, velocity for velocity dispersion
#and if the first moment/second moment values are the same between the 0.65 & 1.3 versions
#then stick to using the old dht ive been doing
#minimize redoing stuff! prioritise writing out what i have already done!

#also the masses are slightly diff bc of changes to velocity (most went down, hence greater 
#discrepancy with the more ones calculated rn)


np.array((['aqs','cam','cra','lac'],[1152.3884,2698.6607,2012.8439,10679.297])).transpose()


cat['mass_GMC<Msun>'][[0,2,5,9]]/[1152.3884,2698.6607,2012.8439,10679.297]


def calc_masses(cubes,distances,types,dvs,dss):
    num = len(cubes)
    masses = np.zeros(num)
    areas = np.zeros(num)
    ma = [fns.mass_n_area(cubes[i]*u.K,distances[i]*u.pc,
                                    types[i],dvs[i],dss[i]) 
                    for i in range(num)]
    
#     print([fns.mass_n_area(cubes[i]*u.K,distances[i]*u.pc,
#                                     types[i],dvs[i],dss[i]) 
#                     for i in range(num)])
    
#     [(print(cubes[i].shape,'\n',
#             distances[i],'\n',
#             types[i],'\n',
#             dvs[i],'\n',
#             dss[i],'\n\n\n'))
#                     for i in range(num)]
    return ma


vo,vi = (((cat['vel<km/s>']+np.array((0,1.3))).transpose() + 100)/1.3).astype('int')


cat['vo'] = vo
cat['vi'] = vi


cos = [
    dht.data[vo[i]:vi[i],
             cat['boc'][i]:cat['bic'][i],
             cat['loc'][i]:cat['lic'][i]]
      for i in range(22)]


inpt = calc_masses(cos,cat['distance<pc>'],
                   np.full(22,'CO'),np.full(22,1.3),
                   np.full(22,0.125))


[(plt.contourf(np.log10(np.nansum(co,axis=0)*1.3*1.97e20),
             levels = 10),
  plt.colorbar(),
plt.show())
             for co in cos]


[(plt.contourf(np.log10(np.nansum(co,axis=0)*1.3*1.97e20),
             levels = np.arange(20,22,0.25)),
plt.show())
             for co in cos]


newms, newas = np.array([(np.around(m.to(u.solMass).value,2),
  np.around(a.value,2))
 for i in range(22)
 for m,a in [inpt[i]]]).transpose()


newas/(np.pi*cat['radius_GMC<pc>']**2)


cos[7].shape


cat[7]


360-166.5


plt.contourf(np.log10(np.nansum(cos[7],axis=0)*1.3*1.97e20),
             levels=np.arange(20,22,0.15),cmap='Paired',extent=[-164.5,-173.5,-3.5,2.5])
plt.colorbar()
plt.hlines([-1.5,-1],-173,-165)
plt.vlines([188-360,194-360,-166.5,192-360],-3,2)
plt.xlim(-165,-173)


cat


cat['lon_GMC<deg>'][7] = [-166.5,-172]
cat['lat_GMC<deg>'][7] = [-1,2.5]


(cat['lon_GMC<deg>'][7]-360)/-0.125


(cat['lat_GMC<deg>'][7]+50)/0.125


cat[['loc','lic']][7] = (cat['lon_GMC<deg>'][7]-360)/-0.125


cat[['boc','bic']][7] = (cat['lat_GMC<deg>'][7]+50)/0.125


cat[['boc','bic']][7] = (cat['lat_GMC<deg>'][7]+50)/0.125


cat[7]


gemcube = (dht.data*(dht.data>0.2*3))[cat['vo'][7]:cat['vi'][7],
         392:420,
         4212:4256]


cat['boc'][7] = 392


cat['lic'][7] = 4256


cat['loc'][7] = 4212


cat.columns


newcat = cat[['name','lon_GMC<deg>','lat_GMC<deg>',
'distance<pc>','vel<km/s>','loc','lic','boc','bic','vo','vi']]


newcat.write('lastlastLAST_cat.fits')


gemom1 = fns.get_mom1(gemcube,np.arange(cat['vo'][7],cat['vi'][7])*1.3-100)


plt.contourf(np.nansum(gemcube,axis=1),levels=np.arange(1,55,5))
plt.colorbar()


plt.contourf(np.nansum(gemcube,axis=2),levels=np.arange(1,100,10))
plt.colorbar()


# plt.contourf(gemom1,np.arange(-5,5,3),cmap='Purples')
# plt.contourf(gemom1,np.arange(5,15,3),cmap='Oranges')

plt.contourf(gemom1,np.arange(-5,15,3),cmap='Spectral',levels=20)


# plt.colorbar()


plt.hist(np.ravel(gemom1),bins=60)
plt.show()


plt.plot(np.arange(cat['vo'][7],cat['vi'][7])*1.3-100,np.nanmean(gemcube,axis=(1,2)))


plt.contourf(np.log10(np.nansum(gemcube,axis=0)*1.3*1.97e20),levels=np.arange(20.5,23,0.25),
            extent=[-166.5,-172,-1,2.5])
plt.colorbar()

plt.contour(np.log10(np.nansum(gemcube,axis=0)*1.3*1.97e20),levels=[21],#,21.5],
            extent=[-166.5,-172,-1,2.5],colors='red')
plt.xlim(-167,-172)


plt.contourf(np.log10(np.nansum(gemcube,axis=0)*1.3*1.97e20),levels=np.arange(20.5,23,0.25))
plt.colorbar()


(dht.data*(dht.data>0.2*3))[cat['vo'][7]:cat['vi'][7],
         cat['boc'][7]:cat['bic'][7],
         cat['loc'][7]:cat['lic'][7]]


fns.mass_n_area(gemcube*u.K,cat['distance<pc>'][7]*u.pc,'CO')[0].to(u.solMass)/1e6


fns.mass_n_area(gemcube*u.K,cat['distance<pc>'][7]*u.pc,'CO')[1]


fns.mass_n_area(gemcube*u.K,cat['distance<pc>'][7]*u.pc,'CO')[0].to(u.solMass)/1e6


fns.mass_n_area(gemcube*u.K,cat['distance<pc>'][7]*u.pc,'CO')[1]


np.sqrt((9959.5911*u.pc**2)/np.pi)


1264224.2/808512.5584007167


cat['mass_GMC<Msun>'][7]


fns.mass_n_area(gemcube*u.K,cat['distance<pc>'][7]*u.pc,'CO')[0]


plt.contourf(np.log10(np.nansum(cos[7],axis=1)*0.125*1.97e20),
             levels=np.arange(20,22,0.15),cmap='Paired',extent=[-164.5,-173.5,-3.5,2.5])
plt.colorbar()
plt.vlines([188-360,194-360],-3,2)


plt.contourf(np.log10(np.nansum(cos[7],axis=2)*0.5*0.125*1.97e20),
             levels=np.arange(20,21.5,0.15),cmap='Paired',extent=[-3.5,2.5,-30,40])
plt.colorbar()
plt.contour(np.log10(np.nansum(cos[7],axis=2)*0.5*0.125*1.97e20),
            levels=[20.5,20.75],colors='black',extent=[-3.5,2.5,-30,40])

plt.vlines(-1.5,-30,40)


gemom1 = fns.get_mom1(cos[7],np.arange(53,108)*1.3-100)


np.gradient(np.histogram(np.ravel(gemom1),bins=np.arange(-10,30,0.65))[0])


cat[7]


gemom0 = np.log10(np.nansum(cos[7],axis=0)*1.3*1.97e20)


gemhi = hi4pi.data[cat['vo'][7]:cat['vi'][7],cat['boh'][7]:cat['bih'][7],cat['loh'][7]:cat['lih'][7]]


gemhimom1 = fns.get_mom1(gemhi,np.arange(53,108)*1.3-100)


plt.hist(np.ravel(gemhimom1),bins=np.arange(-10,30,0.65))

plt.yscale('log')
plt.show()


plt.hist(np.ravel(gemom1),bins=np.arange(-10,30,0.65/2),alpha=0.5,density=True)
plt.hist(np.ravel(gemhimom1),bins=np.arange(-10,30,0.65/2),alpha=0.5,density=True)

# plt.hist(np.ravel(gemom0*gemom1),bins=30)

# plt.yscale('log')
plt.show()


ext = [cat['lon_GMC<deg>'][7][0],cat['lon_GMC<deg>'][7][1],cat['lat_GMC<deg>'][7][0],cat['lat_GMC<deg>'][7][1]]

plt.contourf(gemom1,levels=np.arange(5,10,0.65),cmap=dc.cmap_min0('Oranges'),extent=ext)
plt.contourf(gemom1,levels=np.arange(10,30,0.65),cmap=dc.cmap_min0('viridis'),extent=ext)
plt.contourf(gemom1,levels=np.arange(0,5,0.65),cmap=dc.cmap_min0('binary'),extent=ext)
plt.contourf(gemom1,levels=np.arange(-10,0,0.65),cmap=dc.cmap_min0('Purples'),extent=ext)

plt.xlim(-165,-173)
plt.hlines([-1.5,-1],-173,-165)
plt.vlines([188-360,194-360,-166.5,192-360],-3.5,2)


plt.contourf(gemhimom1,levels=np.arange(0,7,0.65),cmap=dc.cmap_min0('binary'),extent=ext)
plt.contourf(gemhimom1,levels=np.arange(7,15,0.65),cmap=dc.cmap_min0('viridis'),extent=ext)


# plt.contour(gemom1,levels=np.arange(0,7,0.65),cmap=dc.cmap_min0('Oranges'),extent=ext,linewidths=0.7)
# plt.contour(gemom1,levels=np.arange(7,15,0.65),cmap=dc.cmap_min0('Purples'),extent=ext,linewidths=0.7)


newms/cat['mass_GMC<Msun>']


cat['mass_GMC<Msun>'][[0,2,5,9]]/[1152.3884,2698.6607,2012.8439,10679.297]


newms[[0,2,5,9]]


newms[[0,2,5,9]]/[1152.3884,2698.6607,2012.8439,10679.297]


#hmm camelopardalis is having a rough time with masses
#looking at velocity ranges now


fits.open('DHT03_RCrA_mom.fits')[0].header


cat[5]


np.arange(83)*0.26-0.95


#im actually integrating over more velocity, so makes sense that it is slightly higher (2%)


fits.open('DHT06_Aquila-South_mom.fits')[0].header


np.arange(46)*0.65-9.43


cat[0]


10--12, 20--9


#similar velocity range (not the same overlap, but adds up to being ~20+km/s)


aqsind = fits.open('DHT06_Aquila-South_mom.fits')[0].data.transpose().swapaxes(1,2)


aqsind.shape


cat['vel<km/s>'][0]


(np.arange(46)*0.65-9.43)[:35]


aq2 = fns.get_momN(aqsind[:35:,:],(np.arange(46)*0.65-9.43)[:35],2,dx=0.65)


np.nanmean(np.sqrt(aq2))


#hmmm this has a dispersion of 4.67 vs 0.86 for the one with dame data 
#maybe i should just use this ... idk
#double check by reducing the resolution and seeing if anything changes

#goes down when only looking up to +12 (to make up for the fact that i start at -12, but this at -9)


plt.plot(np.arange(46)*0.65-9.43,np.nanmean(aqsind,axis=(1,2)))


(np.arange(46)*0.65-9.43)[10:27]


aq2 = fns.get_momN(aqsind[10:27:,:],(np.arange(46)*0.65-9.43)[10:27],2,dx=0.65)


np.nanmean(np.sqrt(aq2))


#even more down
#so ig the fact that theres lots of close to 0 values at the edges still increases the standard
#deviation, bc its like, theres *some* data that landed over there, and for that to be possible the 
#true distribution must b wider than just where we're seeing the peak at the moment


randar = np.zeros(100)


randar[30:60] = models.Gaussian1D(3.,15.9,8.54)(np.arange(30))


plt.plot(randar,marker='+')


len(np.arange(-15,15.5,0.5))


fns.get_momN(models.Gaussian1D(3.,15.9,8.54)(np.arange(30)),np.arange(-15,15),2)


fns.get_mom1(models.Gaussian1D(3.,15.9,8.54)(np.arange(30)),np.arange(-15,15))


np.sqrt(48.868592284836126)


fns.get_momN(randar,np.arange(-50,50),2)


fits.open('DHT19_Camelopardalis_bw_mom.fits')[0].header


np.arange(256)*0.65-82.9


cat[2]


#Dv is 160 vs 40 km/s so makes sense that the mass from the ind survey is ~3-4 times


camind = fits.open('DHT19_Camelopardalis_bw_mom.fits')[0].data


camind.shape


(np.array((-20,20))+82.9)/0.65 #just integrating over the velocity range i chose


fns.mass_n_area((camind.transpose().swapaxes(1,2))[96:158,:,:]*u.K,
                cat['distance<pc>'][2]*u.pc,'CO',dv=0.65)[0].to(u.solMass)


newms[2] 


969/846.79 #yay this worked


plt.contourf(np.log10(np.nansum((camind.transpose().swapaxes(1,2))[96:158,:,:],axis=0)*0.65*1.97e20),
             levels=np.arange(20,22,0.1))


plt.contourf(np.log10(np.nansum((camind.transpose().swapaxes(1,2)),axis=0)*0.65*1.97e20),
             levels=np.arange(20,22,0.1))


#so weird that it makes such a big difference when the velocity channels outside of +/- 20
#dont even have that much added emission but whateves


cam2 = fns.get_momN((camind.transpose().swapaxes(1,2))[96:158,:,:],np.arange(96,158)*0.65-82.9,2,dx=0.65)


np.nanmean(np.sqrt(cam2))


cat[2]


#the velocity dispersion is the same! and we can trust the one with 0.65 bc ts larger than one channel width
#wooo


np.arange(96,158)*0.65-82.9


fits.open('DHT12_Lacerta_mom.fits')[0].header


np.arange(92)*0.65-29.58


cat[9]


30--30,


#the n of velocities for lacerta is only off by dv 8km/s for me, but its ~64 less, but prob
#this was bc those -30->-22 had significant emission
