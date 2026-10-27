#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle,Ellipse

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


new = fcat[['name','lon_GMC<deg>','lat_GMC<deg>','distance<pc>']]


newvels = np.array(([-10,10],
[-35,20],
[-25,20],
[-20,20],
[6,30],
[-20,25],
[-10,30],
[-15,35],
[-10,40],
[-20,20],
[-10,20],
[20,35],
[-5,15],
[-10,25],
[-10,20],
[-15,20],
[-10,30],
[-20,30],
[-20,20],
[-30,20],
[-10,30],
[-15,20]))


new['vel<km/s>'] = newvels


loc,lic = ((new['lon_GMC<deg>']-360)/-0.125).astype('int').transpose()
boc,bic = ((new['lat_GMC<deg>']+50)/0.125).astype('int').transpose()


vo,vi = ((new['vel<km/s>']+100)/1.3).transpose().astype('int')


vo


vo


h3sig = hi4pi.data*(hi4pi.data>0.043*3)
d3sig = dht.data*(dht.data>0.2*3)


cos = [d3sig[vo[i]:vi[i],boc[i]:bic[i],loc[i]:lic[i]] for i in range(22)]


vels = [np.arange(vo[i],vi[i])*1.3-100 for i in range(22)]


masses = np.zeros(22)
radii = np.zeros(22)
for i in range(22):
    cube = cos[i]*u.K
        
    d = fcat['distance<pc>'][i]*u.pc
    
    m,a = fns.mass_n_area(cube,d,'CO')
    
    masses[i] = m.to(u.solMass).value
    radii[i] = np.sqrt(a.to(u.pc**2)/np.pi).value
    

new['mass_GMC<Msun>'] = masses
new['radius_GMC<pc>'] = radii


new['radius_HI<pc>'] = (fns.accumulation_radius(masses*u.solMass,\
                        fcat['distance<pc>']*u.pc)[1]).value


plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


lcenter = np.nanmean(new['lon_GMC<deg>'],axis=1)
bcenter = np.nanmean(new['lat_GMC<deg>'],axis=1)

sh = (new['radius_HI<pc>'].value/new['distance<pc>'])*u.rad.to(u.deg)


loh,lih = np.around((np.array([lcenter+sh,lcenter-sh]).transpose() \
            -360)/-0.125).astype('int').transpose()

boh,bih = np.around((np.array([bcenter-sh,bcenter+sh]).transpose()\
            +50)/0.125).astype('int').transpose()

his = [h3sig[vo[i]:vi[i],boh[i]:bih[i],loh[i]:lih[i]] for i in range(22)]


h0 = np.nansum(h3sig,axis=0)*1.82e18*1.3
c0 = np.nansum(d3sig,axis=0)*1.3*1.97e20


plt.rcParams['font.size']=40


hlv = np.nansum(h3sig,axis=1)*1.82e18*0.125
clv = np.nansum(d3sig,axis=1)*1.97e20*0.125


hbv = np.nansum(h3sig,axis=2)*1.82e18*0.125*0.5
cbv = np.nansum(d3sig,axis=2)*1.97e20*0.125*0.5


names = ['Aquila South', 'California', 'Camelopardalis', 'Chamaeleon',
       'Canis Major', 'Corona Australis', 'Crossbones', 'Gemini',
       'Hercules', 'Lacerta', 'Lupus', 'Maddalena', 'Monoceros OB1', 'Monoceros R2',
       'Ophiuchus', 'Orion A', 'Orion B', r'$\lambda$-Orionis', 'Perseus', 'Polaris',
       'Rosette', 'Taurus']


abvs = ['AqS','Cal','Cam','Chm','CMa','CrA','Xbn','Gem','Her','Lac','Lup','Mad','MOB','MR2',
        'Oph','OrA','OrB','$\lambda$Or','Per','Pol','Ros','Tau']


name_locs = np.array(([30., -15.5],#aqs
[153., -12],#cal
[150.0, 20.],#cam
[-66.0, -15.5],#chm
[-130,2],#cma
[-3, -21],#cra
[-133, -10.0],#xbn
[-169, 7],#gem
[42.0, 10],#her
[94, -12],#lac
[-26.0, 17],#lup
[-141.0,1],#mad
[-155, 6], #mob
[-134, -14.0], #mr2
[0, 27],#oph
[-146., -24],#ora
[-159, -19],#orb
[-163, -13.5], #lor
[150, -21.0], #per
[117.0, 27.5], #pol
[-155., -4.0], #ros
[175, -30])) #tau


fig,ax = plt.subplots(1,1,figsize=(43,12),dpi=200)

ax.contourf(np.log10(h0),
                 extent=[360,-360,-50,50],
                 levels = np.arange(20,22.5,0.2),
                cmap='YlGnBu_r')


for i in range(22):
    ax.add_patch(Circle((np.nanmean(new['lon_GMC<deg>'][i]),
                         np.nanmean(new['lat_GMC<deg>'][i])),
                       sh[i],fill=False,color='white',linewidth=1.5,linestyle='--'))
    
    ax.annotate(abvs[i],name_locs[i],color='white')


ax.contourf(np.log10(c0),
                 extent=[360,-360,-50,50],
                 levels = np.arange(20,22.5,0.2),
                cmap='magma')


plt.xlabel(r'Galactic Longitude')
plt.ylabel(r'Galactic Latitude')


plt.xticks([180,120,60,0,-60,-120,-180],[r'$180\degree$',r'$120\degree$',r'$60\degree$',r'$0\degree$',
                                        r'$-60\degree$',r'$-120\degree$',r'$-180\degree$'])

plt.yticks([-40,-20,0,20,40],[r'$-40\degree$',r'$-20\degree$',r'$0\degree$',r'$20\degree$',r'$40\degree$'])


plt.xlim(185,-185)
plt.ylim(-45,45)

plt.savefig('catalog_0th')
    
plt.show()


name_locs_lv = np.array(([37, -7],#aqs
[160., 7],#cal
[148.0, 0],#cam
[-63, 10],#chm
[-128,23],#cma
[-1, -8],#cra
[-133, 4],#xbn
[-174, 15],#gem
[39, 10],#her
[98, -12],#lac
[-15, 14],#lup
[-141.0,35],#mad
[-159, 16], #mob
[-142, 15.5], #mr2
[10, 12],#oph
[-146., -10],#ora
[-156, -5],#orb
[-165, -11], #lor
[163, -18], #per
[120.0, 0], #pol
[-152., 20], #ros
[170, 13])) #tau


fig,ax = plt.subplots(1,1,figsize=(40,12),dpi=200)

ax.contourf(np.log10(hlv),
                 extent=[360,-360,-100,100],
                 levels = np.arange(20.5,22,0.15),
                cmap='YlGnBu_r')


ax.contourf(np.log10(clv),
                 extent=[360,-360,-100,100],
                 levels = np.arange(20,22,0.15),
                cmap='magma')


for i in range(22):
    ax.add_patch(Ellipse((np.nanmean(new['lon_GMC<deg>'][i]),
                         np.nanmean(new['vel<km/s>'][i])),
                        2*sh[i],2*hsig[i],
                        fill=False,color='white',linewidth=1.5,linestyle='--'))
    ax.annotate(abvs[i],name_locs_lv[i],color='white')


plt.xlabel(r'Galactic Longitude')
plt.ylabel(r'LSR Velocity $[km\ s^{-1}]$')


plt.xticks([180,120,60,0,-60,-120,-180],[r'$180\degree$',r'$120\degree$',r'$60\degree$',r'$0\degree$',
                                        r'$-60\degree$',r'$-120\degree$',r'$-180\degree$'])

plt.yticks([-60,-30,-15,0,15,30,60])#,[r'$-60\degree$',r'$-40\degree$',r'$-20\degree$',r'$0\degree$',
#                                         r'$20\degree$',r'$40\degree$',r'$60\degree$'])


plt.xlim(185,-185)
plt.ylim(-75,75)

plt.savefig('catalog_lv')
    
plt.show()


# fiddle with:
# 1) lv annotation placement
# 2) crop the two images
# 3) put into overleaf
# 4) the yaxes for spectra


plt.rcParams['font.size']=20


fig,axs = plt.subplots(6,4,figsize=(30,30),dpi=200,sharex=True,sharey=True)

for ax,i in zip(axs.flat,range(22)):
    
    
    ax.plot(np.arange(155)*1.3-100,
            np.nanmean(h3sig[:,boh[i]:bih[i],loh[i]:lih[i]],
                       axis=(1,2)),color='black',linewidth=0.7,
                       linestyle=':',label='HI envelope')

    ax2 = ax.twinx()

    ax2.plot(np.arange(155)*1.3-100,
            np.nanmean(d3sig[:,boc[i]:bic[i],loc[i]:lic[i]],
                       axis=(1,2)),color='xkcd:algae',
                       linewidth=0.7,label='GMC')


    ax2.tick_params(axis='y',colors='xkcd:algae')
    ax2.label_outer()#set_yticks([])

    ax.set_title(names[i],loc='left',y=0.8,x=0.1)
    
    ax2.set_ylim(top=1.2)
        
    
    if i == 3:
        ax.legend(loc=(0.5,0.2),frameon=False)
        ax2.legend(loc=(0.5,0.1),frameon=False)
        
        
for ax in axs.flat:
#     ax.label_outer()
    ax.set_xlim(-40,80)
    ax.set_xticks([-30,0,30,60],[-30,0,30,60])
    
plt.subplots_adjust(wspace=0,hspace=0)

plt.show()


h2tops = [0.3,0.3,0.3,0.3,0.4,0.4,0.4,0.4]


#figure ^ this out later (ordering them so that y axes can be shared both for
#hi and h2)


idx = [0,3,9,19,1,2,5,6,8,14,10,12,13,15,16,21,4,7,11,20,17,18]


shorthi - 0, 5, 18, 19 
medhi - 1, 2, 3, 6    8,9,10,12  13,12,15,16  17,21
tallhi - 4, 7, 11, 20


shorth2 - 0, 2, 9, 10, 17
medh2
tallh2  - 13, 15, 16, 21


shorthi 0,1,2,3, 5,8,9,10, 19,
medhi 6,12,13,14, 15,16,17,18, 21
tallhi 4,7,11,20

shorth2 0,2,9,10, 19,17,
medh2 1,3,6,8, 13,14,18,
tall2h2 4,5,7,11, 12,15,16,20, 21,


fig,axs = plt.subplots(6,4,figsize=(20,30),dpi=200)

for ax,i in zip(axs.flat,range(22)):    
    l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
    b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
    l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
    b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
      
    ax.contourf(np.log10(np.nansum(his[i],axis=0)*1.3*1.82e18),
                 extent=[l0h,l1h,b0h,b1h],
                 levels = np.arange(20,22,0.1),
                cmap='YlGnBu_r',linewidths=0.5)
    
    
    ax.contour(np.log10(np.nansum(cos[i],axis=0)*1.3*1.97e20),
                 extent=[l0c,l1c,b0c,b1c],
                 levels = np.arange(19.5,22.5,0.1),
                cmap='magma',alpha=0.5)
    
#     plt.colorbar()

    xo = np.nanmax([l0h,l0c])
    xi = np.nanmin([l1h,l1c])
    
    ax.set_xlim(xo,xi)
plt.show()


masses_hi = np.zeros(22)
radii_hi = np.zeros(22)
for i in range(22):
    cube = his[i]*u.K
        
    d = fcat['distance<pc>'][i]*u.pc
    
    m,a = fns.mass_n_area(cube,d,'21cm')
    
    masses_hi[i] = m.to(u.solMass).value
    radii_hi[i] = np.sqrt(a.to(u.pc**2)/np.pi).value


plt.scatter(fcat['name'],masses)

plt.xticks(rotation='vertical')

plt.hlines([10,100,1e2,1e3,1e4,1e5,1e6],'aquila_s','taurus')
plt.yscale('log')

plt.show()


new['mass_HI<Msun>'] = masses_hi


fttr = fitting.LinearLSQFitter()
init = models.Linear1D()


OmHs = np.zeros(22)
OmCs = np.zeros(22)

ThHs = np.zeros(22)
ThCs = np.zeros(22)

for i in range(22):
    this_hi = his[i]
    this_co = cos[i]
    
    h0 = np.nansum(this_hi,axis=0)*1.82e18*1.3
    c0 = np.nansum(this_co,axis=0)*1.97e20*1.3
    
    
    h1 = fns.get_mom1(this_hi,np.arange(vo[i],vi[i])*1.3-100)
    c1 = fns.get_mom1(this_co,np.arange(vo[i],vi[i])*1.3-100)
    
    
    h2 = fns.get_momN(this_hi,np.arange(vo[i],vi[i])*1.3-100,2)
    c2 = fns.get_momN(this_co,np.arange(vo[i],vi[i])*1.3-100,2)
    
    
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
    hdf = hdf_noisy[(hdf_noisy['w']>3e20) &\
                    (hdf_noisy['s']>0) & \
                    ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] 
    #make sure that the linewidths are a valid number 
    
    
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
    
    #converting to pc now makes life wayyy easier
    hxpc = hdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*fcat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*fcat['distance<pc>'][i]

    
    dvdxh = fttr(init,hxpc,hdf['v']) 
    dvdyh = fttr(init,hypc,hdf['v'])
    
    dvdxc = fttr(init,cxpc,cdf['v'])
    dvdyc = fttr(init,cypc,cdf['v'])
      
    
    omh = np.sqrt(dvdxh.slope**2 + dvdyh.slope**2)
    omc = np.sqrt(dvdxc.slope**2 + dvdyc.slope**2)

    thh = np.arctan(dvdyh.slope/dvdxh.slope)*u.rad.to(u.deg)
    thc = np.arctan(dvdyc.slope/dvdxc.slope)*u.rad.to(u.deg)
    
#     print(new['name'][i],dvdyc.slope,dvdxc.slope)
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
    

#     plt.scatter(cxpc,cdf['v'],s=1)
#     plt.plot(cxpc,cxpc*dvdxc.slope + dvdxc.intercept,color='red')
#     plt.show()
    
    
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

    
new['vel_grad_GMC<km/s/pc>'] = OmCs
new['vel_grad_HI<km/s/pc>'] = OmHs

new['rot_ax_GMC<deg>'] = ThCs-90
new['rot_ax_HI<deg>'] = ThHs-90

for i in range(22):
    if new['rot_ax_GMC<deg>'][i] < -180:
        new['rot_ax_GMC<deg>'][i] +=360
        
    if new['rot_ax_HI<deg>'][i] < -180:
        new['rot_ax_HI<deg>'][i] +=360


new


jgmc = (new['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc*\
        ((new['radius_GMC<pc>']*u.pc)**2)).cgs
jhi = (new['vel_grad_HI<km/s/pc>']*u.km/u.s/u.pc*\
       ((new['radius_HI<pc>']*u.pc)**2)).cgs


Lgmc = (new['mass_GMC<Msun>']*u.solMass*jgmc).cgs
Lhi = (new['mass_HI<Msun>']*u.solMass*jhi).cgs


new['ang_mom_GMC<erg.s>'] = Lgmc
new['ang_mom_HI<erg.s>'] = Lhi

new['spec_ang_mom_GMC<cm2/s>'] = jgmc
new['spec_ang_mom_HI<cm2/s>'] = jhi


hsig = np.zeros(22)
csig = np.zeros(22)


for i in range(22):    
#     l0h,l1h = np.array((loh[i],lih[i]))*-0.125+360
#     b0h,b1h = np.array((boh[i],bih[i]))*0.125-50
    
#     l0c,l1c = np.array((loc[i],lic[i]))*-0.125+360
#     b0c,b1c = np.array((boc[i],bic[i]))*0.125-50
    
    hispec = np.nanmean(his[i],axis=(1,2))
    hspmom2 = np.sqrt(fns.get_momN(hispec,vels[i],2))
    hsig[i] = hspmom2
    
    cospec = np.nanmean(cos[i],axis=(1,2))
    cspmom2 = np.sqrt(fns.get_momN(cospec,vels[i],2))
    csig[i] = cspmom2
   

new['vel_dispersion_GMC<km/s>'] = csig
new['vel_dispersion_HI<km/s>'] = hsig


alphavir_gmc = (((new['vel_dispersion_GMC<km/s>']*u.km/u.s)**2)*\
(new['radius_GMC<pc>']*u.pc)/\
(c.G*new['mass_GMC<Msun>']*u.solMass)).cgs


alphavir_hi = (((new['vel_dispersion_HI<km/s>']*u.km/u.s)**2)*\
(new['radius_HI<pc>']*u.pc)/\
(c.G*new['mass_HI<Msun>']*u.solMass)).cgs


betarot_gmc = (((new['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)**2 *\
               (new['radius_GMC<pc>']*u.pc)**3)/(c.G*new['mass_GMC<Msun>']*u.solMass)).cgs

betarot_hi = (((new['vel_grad_HI<km/s/pc>']*u.km/u.s/u.pc)**2 *\
               (new['radius_HI<pc>']*u.pc)**3)/(c.G*new['mass_HI<Msun>']*u.solMass)).cgs


delthet = new['rot_ax_GMC<deg>'] - new['rot_ax_HI<deg>']

for i in range(22):
    if delthet[i] > 180:
        delthet[i] -= 360
            
    elif delthet[i] < -180:
        delthet[i] += 360


new['a_vir_GMC'] = alphavir_gmc
new['a_vir_HI'] = alphavir_hi

new['b_rot_GMC'] = betarot_gmc
new['b_rot_HI'] = betarot_hi


new['drot_ax<deg>'] = delthet


Sigmah = (new['mass_HI<Msun>']*u.solMass)/(np.pi*(new['radius_HI<pc>']*u.pc)**2)


Sigmah2 = (new['mass_GMC<Msun>']*u.solMass)/(np.pi*(new['radius_GMC<pc>']*u.pc)**2)


new['Sigma_HI<Msun/pc2>'] = Sigmah
new['Sigma_GMC<Msun/pc2>'] = Sigmah2


new['loh'] = loh
new['lih'] = lih
new['boh'] = boh
new['bih'] = bih

new['loc'] = loc
new['lic'] = lic
new['boc'] = boc
new['bic'] = bic


hi_hdus = fits.HDUList([fits.PrimaryHDU(his[0])]+[fits.ImageHDU(his[i]) for i in range(1,22)])


hi_hdus.writeto('hi_clouds.fits')


co_hdus = fits.HDUList([fits.PrimaryHDU(cos[0])]+[fits.ImageHDU(cos[i]) for i in range(1,22)])


co_hdus.writeto('h2_clouds.fits')


new['vo'] = vo
new['vi'] = vi


new.write('post_flash_cat.fits')
