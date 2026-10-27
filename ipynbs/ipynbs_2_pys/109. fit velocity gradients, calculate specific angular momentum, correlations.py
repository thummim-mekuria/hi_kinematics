#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting

import fun as fns
import datacube as dc

hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data


fcat = QTable(fits.open('re_expanded_cat.fits')[1].data)


lcenter = np.nanmean(fcat['lon_GMC<deg>'],axis=1)
bcenter = np.nanmean(fcat['lat_GMC<deg>'],axis=1)

sh = (fcat['radius_HI<pc>']/fcat['distance<pc>'])*u.rad.to(u.deg)

lpixhi = ((np.array([lcenter+sh,lcenter-sh]).transpose() -360)/-0.125).astype('int')
bpixhi = ((np.array([bcenter-sh,bcenter+sh]).transpose()+50)/0.125).astype('int')

lpixgmc = ((fcat['lon_GMC<deg>'] -360)/-0.125).astype('int')
bpixgmc = ((fcat['lat_GMC<deg>'] +50)/0.125).astype('int')


vo,vi = np.array([((fcat['vel<km/s>'][i]+100)/1.3).astype('int') for i in range(22)]).transpose()


loc,lic = np.array(lpixgmc).transpose()
boc,bic = np.array(bpixgmc).transpose()

loh,lih = np.array(lpixhi).transpose()
boh,bih = np.array(bpixhi).transpose()


vs = [np.arange((vi-vo)[i])*1.3 + fcat['vel<km/s>'][i][0] for i in range(22)]


h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


his = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]
cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


masses = np.zeros(22)
radii = np.zeros(22)
for i in range(22):
    cube = his[i]*u.K
        
    d = fcat['distance<pc>'][i]*u.pc
    
    m,a = fns.mass_n_area(cube,d,'21cm')
    
    masses[i] = m.to(u.kg).value
    radii[i] = np.sqrt(a.to(u.pc**2)/np.pi).value
    
    
fcat['mass_HI<kg>'] = masses
# fcat['radius_GMC<pc>'] = radii

# s, ra = fns.accumulation_radius(masses*u.kg,fcat['distance<pc>']*u.pc)
# sh = s.value
# fcat['radius_HI<pc>'] = ra
# fcat.write('re_expanded_cat.fits',overwrite=True)


fcat


np.array([np.arange(22),fcat['name']]).transpose()


fig,ax = plt.subplots(1,1,figsize=(20,7))
for i in range(22):
    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[loh,lih,boh,bih],cmap='Paired',alpha=0.5,levels=np.arange(20,22,0.2))
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))

plt.xlim(180,-180)
plt.ylim(-40,40)
plt.show()


fcat[['name','distance<pc>']]


his[0].shape


vs[0].shape


(((h0/u.cm**2)*1.36*c.m_p*h2*(u.km/u.s)**2/(fcat['radius_HI<pc>'][i]*u.pc))/c.k_B).to(u.K/u.cm**3)


plt.scatter(np.log10(np.ravel(c0)),np.log10(np.ravel(c_pmp.value)),s=1e-1)
plt.ylim(0,)


for i in range(22):
#     fig,axs = plt.subplots(1,2,figsize=(20,10))

    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
    
    h_turb = (((h0/u.cm**2)*1.36*c.m_p*h2*(u.km/u.s)**2/\
              (fcat['radius_HI<pc>'][i]*u.pc))/c.k_B).to(u.K/u.cm**3)
    c_turb = (((c0/u.cm**2)*2*1.36*c.m_p*c2*(u.km/u.s)**2/\
              (fcat['radius_GMC<pc>'][i]*u.pc))/c.k_B).to(u.K/u.cm**3)

    plt.scatter(np.ravel(np.sqrt(h2)),np.ravel(h_turb),s=1e-3)
    plt.scatter(np.ravel(np.sqrt(c2)),np.ravel(c_turb),s=5e-1)
    
    
#     print(fttr(init,np.ravel(h0),np.ravel(h_pmp)))
#     print(fttr(init,np.ravel(c0),np.ravel(c_pmp)))
#     print('\n\n')
    
    plt.loglog()

#     plt.ylim(1e-2)
#     imh0 = axs[0].contourf(h_pmp.value/1e4,vmin=0, #h2 is already sigma^2
#                     extent=[loh,lih,boh,bih],cmap='YlGnBu',levels=10)
    
#     imc0 = axs[1].contourf(c_pmp.value/1e4,vmin=0,
#                     extent=[loc,lic,boc,bic],cmap='Purples',levels=10)

#     axs[0].contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))
#     for ax in axs.flat:
#         ax.set_xlim(loh,lih)
#         ax.set_ylim(boh,bih)
        
#     for im,ax in zip([imh0,imc0],axs):
#     fig.colorbar(im,ax=ax)
    
    
    plt.ylim(1e2,)
    plt.xlim(1e-1,)
    plt.show()

#why is this giving nan????


for i in range(22):
#     fig,axs = plt.subplots(1,2,figsize=(20,10))

    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
    
    h_pmp = (((np.pi*c.G*((h0/u.cm**2)*1.36*c.m_p)**2)*1.5)/c.k_B).to(u.K/u.cm**3)
    c_pmp = (((np.pi*c.G*((c0/u.cm**2)*1.36*c.m_p)**2)*1.5)/c.k_B).to(u.K/u.cm**3)
    
    h_pturb = (((h0/u.cm**2)*1.36*c.m_p*h2*(u.km/u.s)**2/\
              (fcat['radius_HI<pc>'][i]*u.pc))/c.k_B).to(u.K/u.cm**3)
    c_pturb = (((c0/u.cm**2)*2*1.36*c.m_p*c2*(u.km/u.s)**2/\
              (fcat['radius_GMC<pc>'][i]*u.pc))/c.k_B).to(u.K/u.cm**3)

    plt.scatter(np.ravel(h_pturb),np.ravel(h_pmp),s=1e-3)
    plt.scatter(np.ravel(c_pturb),np.ravel(c_pmp),s=5e-1)
    
    
#     print(fttr(init,np.ravel(h0),np.ravel(h_pmp)))
#     print(fttr(init,np.ravel(c0),np.ravel(c_pmp)))
#     print('\n\n')
    
    plt.loglog()

#     plt.ylim(1e-2)
#     imh0 = axs[0].contourf(h_pmp.value/1e4,vmin=0, #h2 is already sigma^2
#                     extent=[loh,lih,boh,bih],cmap='YlGnBu',levels=10)
    
#     imc0 = axs[1].contourf(c_pmp.value/1e4,vmin=0,
#                     extent=[loc,lic,boc,bic],cmap='Purples',levels=10)

#     axs[0].contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))
#     for ax in axs.flat:
#         ax.set_xlim(loh,lih)
#         ax.set_ylim(boh,bih)
        
#     for im,ax in zip([imh0,imc0],axs):
#     fig.colorbar(im,ax=ax)
    
    
    plt.ylim(1e2,)
    plt.xlim(1e2,)
    plt.show()

#why is this giving nan????


SFR


fcat


fig,ax = plt.subplots(1,1,figsize=(20,10))
for i in [6,13]:
    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[loh,lih,boh,bih],cmap='Paired',alpha=0.5,levels=np.arange(20,22,0.2))
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))

# plt.xlim(-135,-175)
# plt.ylim(-31,-3)
plt.show()


fig,ax = plt.subplots(1,1,figsize=(20,10))
for i in [4,7,11,12]:
    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[loh,lih,boh,bih],cmap='Paired',alpha=0.5,levels=np.arange(20,22,0.2))
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))

# plt.xlim(-135,-175)
# plt.ylim(-31,-3)
plt.show()


fig,ax = plt.subplots(1,1,figsize=(15,10))
for i in [18,21,1]:
    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[loh,lih,boh,bih],cmap='Paired',alpha=0.5,levels=np.arange(20,22,0.2))
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

# plt.scatter(-156.3,-15)
# ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
# ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))

# plt.xlim(-135,-175)
# plt.ylim(-31,-3)
plt.show()


for the ones that overlap in the integrated maps but the mass from the total
cloud massses doesnt give a long enough accumulation cylinder to encompass them, 
see if taking the patch by patch spectra of the overlap regions (ughh for perseus x taurus)
can reveal which peaks / portions of the spectrum to assign to which cloud


oph&lupus
orions


lpixgmc[0]


tot_or_mass = fcat['mass_GMC<kg>'][15] + fcat['mass_GMC<kg>'][16] + fcat['mass_GMC<kg>'][17]


tot_tauper_mass = fcat['mass_GMC<kg>'][18] + fcat['mass_GMC<kg>'][21]


tot_percal_mass = fcat['mass_GMC<kg>'][18] + fcat['mass_GMC<kg>'][1]


tot_gemon_mass = fcat['mass_GMC<kg>'][7] + fcat['mass_GMC<kg>'][12]


tot_xbmr2_mass = fcat['mass_GMC<kg>'][6] + fcat['mass_GMC<kg>'][13]


tot_ophlup2_mass = fcat['mass_GMC<kg>'][10] + fcat['mass_GMC<kg>'][14]


fcat['distance<pc>'][[6,13]]


fns.accumulation_radius(tot_ophlup2_mass*u.kg,166.5*u.pc)


144+57


886-92


np.nanmean([144,189])


fcat['distance<pc>'][[14,10]]


np.nanmean([402,423,432])


fns.accumulation_radius(tot_or_mass*u.kg,470*u.pc)


np.nanmean(fcat['lon_GMC<deg>'][[15,16,17]])


np.nanmean(fcat['lat_GMC<deg>'][[15,16,17]])


from matplotlib.patches import Circle


419+(102/2)


fig,ax = plt.subplots(1,1,figsize=(15,10))
for i in [15,16,17]:
    loc,lic = lpixgmc[i]*-0.125 + 360
    boc,bic = bpixgmc[i]*0.125 -50
    
    loh,lih = lpixhi[i]*-0.125 + 360
    boh,bih = bpixhi[i]*0.125 -50
    
    
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    ax.contourf(np.log10(h0),extent=[loh,lih,boh,bih],cmap='Paired',alpha=0.5,levels=np.arange(20,22,0.2))
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    ax.contour(np.log10(c0),extent=[loc,lic,boc,bic],cmap='binary_r')

plt.scatter(-156.3,-15)
ax.add_patch(Circle((-156.3,-15.),15.9,fill=False,edgecolor='red'))
ax.add_patch(Circle((-156.3,-15.),12.45,fill=False,edgecolor='red'))

plt.xlim(-135,-175)
plt.ylim(-31,-3)
plt.show()


[loc*-0.125+360,lic*-0.125+360,boc*0.125-50,bic*0.125-50]


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()

plinit = models.Planar2D(1.,1.,0.)
plfttr = fitting.LevMarLSQFitter()


hvstds = np.zeros(22)
cvstds = np.zeros(22)

for i in range(22):
    this_hi = his[i]
    this_co = cos[i]
    
    hsp = np.nanmean(this_hi,axis=(1,2))
    csp = np.nanmean(this_co,axis=(1,2))

    hsig = np.sqrt(fns.get_momN(hsp,vs[i],2))
    csig = np.sqrt(fns.get_momN(csp,vs[i],2))
    
    hvstds[i] = hsig
    cvstds[i] = csig
#     fig,ax = plt.subplots(1,1,figsize=(7,7))
#     ax.plot(vs[i],hsp)
    
#     ax2 = ax.twinx()
#     ax2.plot(vs[i],csp)
#     plt.title(fcat['name'][i])
#     plt.show()

#from looking at spectra ^ maddalena, cma_ob1 have pretty significant contamination
#from double peaks -> better to dendrogram / something else to find the structures
#corresponding to the peaks? bc doing a velocity cut o its own is gonna lose a lot of information

#^ this change will also affect the calculated masses, velocity dispersion, velocity gradients, maybe everything


fcat['vel_dispersion_GMC<km/s>'] = cvstds
fcat['vel_dispersion_HI<km/s>'] = hvstds


np.array((fcat['name'],fcat['vel_grad_HI<km/s/pc>'])).transpose()


np.array((fcat['name'],fcat['vel_grad_GMC<km/s/pc>'])).transpose()


plt.scatter(OmCs,OmCs_wd)
plt.plot(OmCs_wd,OmCs_wd)
# plt.plot(OmHs,OmHs)


plt.scatter(OmHs,OmHs_wd)
plt.plot(OmHs_wd,OmHs_wd)
# plt.plot(OmHs,OmHs)


OmHs_wd = np.zeros(22)
OmCs_wd = np.zeros(22)

ThHs_wd = np.zeros(22)
ThCs_wd = np.zeros(22)

for i in range(22):
    this_hi = his[i]
    this_co = cos[i]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,vs[i])
    c1 = fns.get_mom1(this_co,vs[i])
    
    
    h2 = fns.get_momN(this_hi,vs[i],2)
    c2 = fns.get_momN(this_co,vs[i],2)
    
    
    xc,yc = np.meshgrid(np.arange(loc[i],lic[i])*-0.125+360,
                        np.arange(boc[i],bic[i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(loh[i],lih[i])*-0.125+360,
                        np.arange(boh[i],bih[i])*0.125-50)
    
    
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0),
            np.ravel(h1),
            np.ravel(np.sqrt(h2))
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    hdf = hdf_noisy[(hdf_noisy['w']>4e20) &\
                    (hdf_noisy['s']>0) & \
                    ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] #make sure that the weights are a valid number
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0),
            np.ravel(c1),
            np.ravel(np.sqrt(c2))
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>4e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]

    
    dvdxh = fttr(init,hxpc,hdf['v'],weights=1/hdf['s']) #converting to pc now makes life wayyy easier
    dvdyh = fttr(init,hypc,hdf['v'],weights=1/hdf['s'])
    
    dvdxc = fttr(init,cxpc,cdf['v'],weights=1/cdf['s'])
    dvdyc = fttr(init,cypc,cdf['v'],weights=1/cdf['s'])
      
    
    omh = np.sqrt(dvdxh.slope**2 + dvdyh.slope**2)
    omc = np.sqrt(dvdxc.slope**2 + dvdyc.slope**2)

    thh = np.arctan(dvdyh.slope/dvdxh.slope)*u.rad.to(u.deg)
    thc = np.arctan(dvdyc.slope/dvdxc.slope)*u.rad.to(u.deg)
    
    
    #correcting for arctan only returning |theta|<90
    if dvdxh.slope > 0 : #first & fourth quadrants
        thh = thh
    elif dvdxh.slope < 0: #second quadrant & third
        thh = thh + 180
        
    if dvdxc.slope > 0 : #first & fourth quadrants
        thc = thc  
    elif dvdxc.slope < 0: #second quadrant & third
        thc = thc - 180

    OmHs_wd[i] = omh
    OmCs_wd[i] = omc
    
    ThHs_wd[i] = thh
    ThCs_wd[i] = thc
    
        
    plt.scatter(cxpc,cdf['v'],s=1)
    plt.plot(cxpc,cxpc*dvdxc.slope + dvdxc.intercept,color='red')
    plt.show()
    
    
#     plt.imshow(h1,origin='lower',cmap='RdGy')
#     plt.show()
#     plt.imshow(c1,origin='lower',cmap='PiYG')
#     plt.show()
    
    
#     print(fcat['name'][i],dvdxh.slope.value,dvdyh.slope.value,
#           dvdxc.slope.value,dvdyc.slope.value,'\n\n')
    
#     dvdxh.slope.value * (0.125*u.deg.to(u.rad)*fcat['distance<pc>'] #km/s/deg * deg/pc
    
#     print(fcat['name'][i],hplane.slope_x,hplane.slope_y,
#           cplane.slope_x,cplane.slope_y,'\n\n\n\n\n')

    
#     plt.scatter(np.ravel(xh),np.ravel(h1),s=1)
#     plt.scatter(np.ravel(xc),np.ravel(c1),s=3)
#     plt.show()

    
OmHs = np.zeros(22)
OmCs = np.zeros(22)

ThHs = np.zeros(22)
ThCs = np.zeros(22)

for i in range(22):
    this_hi = his[i]
    this_co = cos[i]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,vs[i])
    c1 = fns.get_mom1(this_co,vs[i])
    
    
    h2 = fns.get_momN(this_hi,vs[i],2)
    c2 = fns.get_momN(this_co,vs[i],2)
    
    
    xc,yc = np.meshgrid(np.arange(loc[i],lic[i])*-0.125+360,
                        np.arange(boc[i],bic[i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(loh[i],lih[i])*-0.125+360,
                        np.arange(boh[i],bih[i])*0.125-50)
    
    
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0),
            np.ravel(h1),
            np.ravel(h2)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0),
            np.ravel(c1),
            np.ravel(c2)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]
    
    
    hxpc = hdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]

    
    dvdxh = fttr(init,hxpc,hdf['v']) #converting to pc now makes life wayyy easier
    dvdyh = fttr(init,hypc,hdf['v'])
    
    dvdxc = fttr(init,cxpc,cdf['v'])
    dvdyc = fttr(init,cypc,cdf['v'])
      
    
    omh = np.sqrt(dvdxh.slope**2 + dvdyh.slope**2)
    omc = np.sqrt(dvdxc.slope**2 + dvdyc.slope**2)

    thh = np.arctan(dvdyh.slope/dvdxh.slope)*u.rad.to(u.deg)
    thc = np.arctan(dvdyc.slope/dvdxc.slope)*u.rad.to(u.deg)
    
    
    #correcting for arctan only returning |theta|<90
    if dvdxh.slope > 0 : #first & fourth quadrants
        thh = thh
    elif dvdxh.slope < 0: #second quadrant & third
        thh = thh + 180
        
    if dvdxc.slope > 0 : #first & fourth quadrants
        thc = thc  
    elif dvdxc.slope < 0: #second quadrant & third
        thc = thc - 180

    OmHs[i] = omh
    OmCs[i] = omc
    
    ThHs[i] = thh
    ThCs[i] = thc
    
        
#     plt.imshow(h1,origin='lower',cmap='RdGy')
#     plt.show()
#     plt.imshow(c1,origin='lower',cmap='PiYG')
#     plt.show()
    
    
#     dvdxh = fttr(init,np.ravel(xh),np.ravel(h1))
#     dvdyh = fttr(init,np.ravel(yh),np.ravel(h1))
    
#     dvdxc = fttr(init,np.ravel(xc),np.ravel(c1))
#     dvdyc = fttr(init,np.ravel(yc),np.ravel(c1))
    
#     print(fcat['name'][i],dvdxh.slope.value,dvdyh.slope.value,
#           dvdxc.slope.value,dvdyc.slope.value,'\n\n')
    
#     dvdxh.slope.value * (0.125*u.deg.to(u.rad)*fcat['distance<pc>'] #km/s/deg * deg/pc
    
#     print(fcat['name'][i],hplane.slope_x,hplane.slope_y,
#           cplane.slope_x,cplane.slope_y,'\n\n\n\n\n')


    plt.scatter(cxpc,cdf['v'],s=1)
    plt.plot(cxpc,cxpc*dvdxc.slope + dvdxc.intercept,color='red')
    plt.show()
    
    
#     plt.scatter(np.ravel(xh),np.ravel(h1),s=1)
#     plt.scatter(np.ravel(xc),np.ravel(c1),s=3)
#     plt.show()

    
fcat['vel_grad_GMC<km/s/pc>'] = OmCs
fcat['vel_grad_HI<km/s/pc>'] = OmHs

fcat['rot_ax_GMC<deg>'] = ThCs-90
fcat['rot_ax_HI<deg>'] = ThHs-90

for i in range(22):
    if fcat['rot_ax_GMC<deg>'][i] < -180:
        fcat['rot_ax_GMC<deg>'][i] +=360
        
    if fcat['rot_ax_HI<deg>'][i] < -180:
        fcat['rot_ax_HI<deg>'][i] +=360


jgmc = 2/5*fcat['vel_grad_GMC<km/s/pc>']*fcat['radius_GMC<pc>']**2


fcat['spec_ang_mom_GMC<pc.km/s>'] = jgmc


jhi = 2/5*fcat['vel_grad_HI<km/s/pc>']*fcat['radius_HI<pc>']**2
fcat['spec_ang_mom_HI<pc.km/s>'] = jhi


plt.rcParams['font.size'] = 12. 


fttr(init,np.log10(jgmc),np.log10(jhi))


plt.scatter(jgmc,jhi,s=10)
plt.plot(jgmc,10**1.35677745*(jgmc)**(0.44091927),linewidth=0.8)

# plt.plot(jgmc,(41.02*jgmc)**(2.693))
plt.loglog()

plt.xlabel(r'$j_{GMC}$ [$pc\ km\ s^{-1}$]')
plt.ylabel(r'$j_{HI}$ [$pc\ km\ s^{-1}$]')


plt.show()


fttr(init,np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)),np.log10(fcat['radius_GMC<pc>']))


fttr(init,np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)),np.log10(fcat['radius_HI<pc>']))


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),fcat['radius_GMC<pc>'],s=10)
plt.scatter(fcat['mass_HI<kg>']*u.kg.to(u.solMass),fcat['radius_HI<pc>'],s=10)

plt.plot(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
         10**(-0.87036483)*(fcat['mass_GMC<kg>']*u.kg.to(u.solMass))**(0.46116072),linewidth=0.8)

plt.plot(fcat['mass_HI<kg>']*u.kg.to(u.solMass),
         10**(-0.17721721)*(fcat['mass_HI<kg>']*u.kg.to(u.solMass))**(0.35936346),linewidth=0.8)


plt.loglog()

plt.xlabel(r'M [$M_\odot$]')
plt.ylabel(r'R [pc]')


plt.show()


fttr(init,np.log10(fcat['radius_GMC<pc>']),np.log10(fcat['vel_dispersion_GMC<km/s>']))


fttr(init,np.log10(fcat['radius_HI<pc>']),np.log10(fcat['vel_dispersion_HI<km/s>']))


plt.scatter(fcat['radius_GMC<pc>'],fcat['vel_dispersion_GMC<km/s>'],s=10)
plt.scatter(fcat['radius_HI<pc>'],fcat['vel_dispersion_HI<km/s>'],s=10)

plt.plot(fcat['radius_GMC<pc>'],10**(-0.1383996)*(fcat['radius_GMC<pc>'])**(0.42418601),linewidth=0.8)
plt.plot(fcat['radius_HI<pc>'],10**(0.79641596)*(fcat['radius_HI<pc>'])**(0.09783737),linewidth=0.8)

plt.loglog()

plt.xlabel(r'R [pc]')
plt.ylabel(r'$\sigma_{v}$ [$km\ s^{-1}$]')


plt.show()


fttr(init,np.log10((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass)),
     np.log10(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)))


fttr(init,np.log10((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass)),
     np.log10(fcat['mass_HI<kg>']*u.kg.to(u.solMass)))


fcat


plt.scatter((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass),
            fcat['mass_GMC<kg>']*u.kg.to(u.solMass),s=10)

plt.scatter((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass),
            fcat['mass_HI<kg>']*u.kg.to(u.solMass),s=10)


plt.plot((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass),
         10**(-1.38454272)*((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>']
                            )*u.kg.to(u.solMass))**(1.07672233),linewidth=0.8)


plt.plot((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'])*u.kg.to(u.solMass),
         10**(-0.01431748)*((fcat['mass_GMC<kg>']+fcat['mass_HI<kg>'] 
         )*u.kg.to(u.solMass))**(0.9928546),linewidth=0.8)


plt.loglog()

plt.xlabel(r'$M_{GMC + HI}$ [$M_\odot$]')
plt.ylabel(r'$M_{GMC}$ [$M_\odot$]')

plt.show()


plt.hist(fcat['rot_ax_GMC<deg>'],alpha=0.5)
plt.hist(fcat['rot_ax_HI<deg>'],alpha=0.5)

plt.show()


plt.scatter(fcat['vel_dispersion_GMC<km/s>'],
            np.abs(fcat['rot_ax_GMC<deg>']-fcat['rot_ax_HI<deg>']),s=10)

plt.scatter(fcat['vel_dispersion_HI<km/s>'],
            np.abs(fcat['rot_ax_GMC<deg>']-fcat['rot_ax_HI<deg>']),s=10)

plt.loglog()

plt.show()


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),
            np.abs(fcat['rot_ax_GMC<deg>']-fcat['rot_ax_HI<deg>']),s=10)


plt.scatter(fcat['mass_HI<kg>']*u.kg.to(u.solMass),
            np.abs(fcat['rot_ax_GMC<deg>']-fcat['rot_ax_HI<deg>']),s=10)

plt.loglog()

plt.show()


plt.hist(fcat['vel_dispersion_GMC<km/s>'],bins=5,alpha=0.5)
plt.hist(fcat['vel_dispersion_HI<km/s>'],bins=5,alpha=0.5)

plt.show()


plt.scatter(OmCs,fcat['vel_dispersion_GMC<km/s>'],s=10)
plt.scatter(OmHs,fcat['vel_dispersion_HI<km/s>'],s=10)

plt.loglog()

plt.show()


plt.scatter(jgmc,fcat['vel_dispersion_GMC<km/s>'],s=10)
plt.scatter(jhi,fcat['vel_dispersion_HI<km/s>'],s=10)

plt.loglog()

plt.show()


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass),OmCs,s=10)
plt.scatter(fcat['mass_HI<kg>']*u.kg.to(u.solMass),OmHs,s=10)

plt.loglog()

plt.show()


plt.scatter(fcat['radius_GMC<pc>'],OmCs,s=10)
plt.scatter(fcat['radius_HI<pc>'],OmHs,s=10)

plt.loglog()

plt.show()


dispersion of the hi with itself -> split into the hi that is within the borders of the co
and the hi that is in the accumulation radius but 


consider all the hi together (masking out the parts that have co) -> so instead of assigning different pixels to 
one or another cloud, its just, whats going on at this pixel, whats going on at that pixel


to really focus on why better data is different? -> degrade hi4pi & compare to our results


significant part of the paper -> laying out the methodology


surface density of hi as a function of h2 mass


alpha_vir_GMC = ((5*fcat['radius_GMC<pc>']*u.pc*(fcat['vel_dispersion_GMC<km/s>']*u.km/u.s)**2)/(
    c.G*fcat['mass_GMC<kg>']*u.kg)).si

alpha_vir_HI = ((5*fcat['radius_HI<pc>']*u.pc*(fcat['vel_dispersion_HI<km/s>']*u.km/u.s)**2)/(
    c.G*fcat['mass_HI<kg>']*u.kg)).si


plt.hist(alpha_vir_GMC,bins=5,alpha=0.5)
plt.hist(alpha_vir_HI,bins=5,alpha=0.5)


plt.show()


plt.scatter(alpha_vir_GMC,fcat['mass_GMC<kg>']*u.kg.to(u.solMass),s=10)
plt.scatter(alpha_vir_HI,fcat['mass_HI<kg>']*u.kg.to(u.solMass),s=10)

plt.loglog()

plt.show()


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    plt.hist(np.ravel(c0),alpha=0.5)
    plt.hist(np.ravel(h0),alpha=0.5)
    plt.xlim(1e20,)
    plt.loglog()
    
    plt.show()


fcat.write('re_expanded_cat.fits',overwrite=True)


fcat


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0big = 
    
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
 
    plt.scatter(np.ravel(h2),np.ravel(h0),s=1e-2)
    
plt.ylim(1e19,)
plt.loglog()
plt.show()

#a cube with HI pix but from dht? ^ to make this happen


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)

    
    dvdxh,dvdyh = np.gradient(h1)
    dvdxc,dvdyc = np.gradient(c1)
    
    dvhmag = np.sqrt(dvdxh**2 + dvdyh**2)
    dvcmag = np.sqrt(dvdxc**2 + dvdyc**2)
    
    
    plt.scatter(np.ravel(c2),np.ravel(dvcmag),s=1e-2)
    
plt.xlim(1e-1,) 
plt.ylim(1e-3,) 

plt.loglog()
plt.show()


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)

    
    dvdxh,dvdyh = np.gradient(h1)
    dvdxc,dvdyc = np.gradient(c1)
    
    dvhmag = np.sqrt(dvdxh**2 + dvdyh**2)
    dvcmag = np.sqrt(dvdxc**2 + dvdyc**2)
    
    
    plt.scatter(np.ravel(h2),np.ravel(dvhmag),s=1e-2)
    
plt.xlim(1e1,) 
plt.ylim(1e-2,) 

plt.loglog()
plt.show()


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
 

    dvdxh,dvdyh = np.gradient(h1)
    dvdxc,dvdyc = np.gradient(c1)
    
    dndxh,dndyh = np.gradient(h0)
    dndxc,dndyc = np.gradient(c0)
    
    dvhmag = np.sqrt(dvdxh**2 + dvdyh**2)
    dvcmag = np.sqrt(dvdxc**2 + dvdyc**2)
    
    dnhmag = np.sqrt(dndxh**2 + dndyh**2)
    dncmag = np.sqrt(dndxc**2 + dndyc**2)
    
    plt.scatter(np.ravel(dvhmag),np.ravel(dnhmag),s=1e-2)
    
# plt.ylim(1e1,)
# plt.xscale('log')
plt.loglog()
plt.show()


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
 

    dvdxh,dvdyh = np.gradient(h1)
    dvdxc,dvdyc = np.gradient(c1)
    
    dndxh,dndyh = np.gradient(h0)
    dndxc,dndyc = np.gradient(c0)
    
    dvhmag = np.sqrt(dvdxh**2 + dvdyh**2)
    dvcmag = np.sqrt(dvdxc**2 + dvdyc**2)
    
    dnhmag = np.sqrt(dndxh**2 + dndyh**2)
    dncmag = np.sqrt(dndxc**2 + dndyc**2)
    
    plt.scatter(np.ravel(dvcmag),np.ravel(dncmag),s=1e-2)
    
plt.xlim(1e-3,)
plt.loglog()
plt.show()


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
 
    plt.scatter(np.ravel(h0),np.ravel(h2),s=1e-2)
    
plt.ylim(1e1,)
plt.loglog()
plt.show()


for i in range(22):
    h0 = np.nansum(his[i],axis=0)*1.3*1.82e18
    c0 = np.nansum(cos[i],axis=0)*1.3*1.97e20
    
    
    h1 = fns.get_mom1(his[i],vs[i])
    c1 = fns.get_mom1(cos[i],vs[i])
    
    
    h2 = fns.get_momN(his[i],vs[i],2)
    c2 = fns.get_momN(cos[i],vs[i],2)
 
    plt.scatter(np.ravel(c0),np.ravel(c2),s=1e-1)
    
plt.ylim(1e-1,)
plt.loglog()
plt.show()


plt.scatter(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_GMC<pc>']**2),
            fcat['mass_HI<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>']**2),
            s=10)

# plt.loglog()


plt.show()


plt.hist(fcat['mass_GMC<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_GMC<pc>']**2),bins=10,alpha=0.5)
plt.hist(fcat['mass_HI<kg>']*u.kg.to(u.solMass)/(np.pi*fcat['radius_HI<pc>']**2),bins=10,alpha=0.5)

plt.show()


fcat
