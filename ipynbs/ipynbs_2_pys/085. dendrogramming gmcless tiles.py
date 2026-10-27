#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import datacube as dc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
import fun as fns
import astropy.units as u
import astropy.constants as c
from matplotlib.patches import Circle,Rectangle


hi4pi,dht = fits.open('both_same_grid_720d.fits')
zcat = QTable(fits.open('final_cat.fits')[1].data)


lons = np.arange(5760)*-0.125 + 360
lats = np.arange(800)*0.125 - 50
vels = np.arange(155)*1.3 - 100


hi_conv = (1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p)#.to(u.solMass/u.pc**2) 

co_conv = (1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p)#.to(u.solMass/u.pc**2) 


whi = (hi_conv*u.K*1.3*u.km/u.s*np.nansum(hi4pi.data,axis=0)).to(u.solMass/u.pc**2).value
wco = (co_conv*u.K*1.3*u.km/u.s*np.nansum(dht.data,axis=0)).to(u.solMass/u.pc**2).value


whole_gal_ds = pd.read_csv('whole_gal_dists.csv')
zcat = QTable(fits.open('final_cat.fits')[1].data)
gless_cat = QTable(fits.open('gmcless_cat.fits')[1].data)


h2s = fits.open('h2_clouds.fits')
his = fits.open('hi_clouds.fits')
gls = fits.open('gmcless_tiles.fits')


gl_mom0s = [(np.nansum(gls[i].data,axis=0)*hi_conv*u.K*1.3*u.km/u.s).to(u.solMass/u.pc**2).value
            for i in range(7)]
gl_mom1s = [fns.get_mom1(gls[i].data,[0,155]) for i in range(7)]
gl_mom2s = [np.sqrt(fns.get_momN(gls[i].data,[0,155],2)) for i in range(7)]
#just the second moment is the variance


gf_mom0s = [(np.nansum(his[i].data,axis=0)*hi_conv*u.K*1.3*u.km/u.s).to(u.solMass/u.pc**2).value
            for i in range(17)]
gf_mom1s = [fns.get_mom1(his[i].data,[0,155]) for i in range(17)]
gf_mom2s = [np.sqrt(fns.get_momN(his[i].data,[0,155],2)) for i in range(17)]
#just the second moment is the variance


[(plt.imshow(gf_mom0s[i],origin='lower'),
  plt.title(zcat['names'][i]),
  plt.colorbar(),
  plt.show())
 for i in range(17)]


[(plt.imshow(gf_mom1s[i],origin='lower'),
  plt.title(zcat['names'][i]),
  plt.colorbar(),
  plt.show())
 for i in range(17)]


[(plt.imshow(gf_mom2s[i],origin='lower',norm='asinh'),
  plt.title(zcat['names'][i]),
  plt.colorbar(),
  plt.show())
 for i in range(17)]


mom2_gal = fns.get_momN(hi4pi.data,[0,155],2)


mom1_gal = fns.get_mom1(hi4pi.data,[0,155])


dust_lbd = fits.open('mean_and_std_lbd.fits')
mean_dust, std_dust = dust_lbd[1],dust_lbd[2]
long_dust = np.concatenate((mean_dust.data[:,:,786:],mean_dust.data,mean_dust.data[:,:,:786]),axis=2)
# gas2dust = 1.87e21*(u.cm**-2)/u.mag #from imara & burkhart
gas2dust = 2e21*(u.cm**-2)/u.mag #from zhu et al: https://arxiv.org/pdf/1706.07109.pdf

dust_lons = np.arange(1572)*(-0.23) + 181
dust_lats = np.arange(786)*(0.23) -90
dust_ds = np.arange(1182) + 68

wd = np.nansum(mean_dust.data,axis=0)
wd_nh = wd*u.mag*gas2dust


fake_mom2 = np.nanstd(hi4pi.data,axis=0)


plt.figure(figsize=(20,10))
plt.contourf(lons,lats,whi,origin='lower',levels=np.arange(0,30,6),norm='log')#,6,8,10,15,20,25,30])
plt.colorbar(orientation='horizontal')
# plt.contour(lons,lats,mom1_gal,levels=[-10,10],
#             colors='black',linewidths=0.6)
plt.xlim(180,-180)


plt.figure(figsize=(20,10))
plt.contourf(lons,lats,fake_mom2,origin='lower',levels=[0,1,2,3,4,5,10,15,20,25,30],norm='log')
plt.colorbar(orientation='horizontal')


# plt.contourf(dust_lons,dust_lats,wd,levels=[0.2,0.5],colors='red',linewidths=0.7,alpha=0.5)

# plt.contour(lons,lats,mom1_gal,levels=[-10,10],
#             colors='black',linewidths=0.6)
plt.xlim(180,-180)
plt.ylim(-50,50)


plt.figure(figsize=(20,10))
plt.contourf(lons,lats,np.sqrt(mom2_gal),origin='lower',levels=np.arange(0,30,4),norm='log')
plt.colorbar(orientation='horizontal')


plt.contour(dust_lons,dust_lats,wd_nh.value,levels=[1e21],colors='red',linewidths=0.7)

# plt.contour(lons,lats,mom1_gal,levels=[-10,10],
#             colors='black',linewidths=0.6)
plt.xlim(180,-180)
plt.ylim(-50,50)


gless_cat['lat<deg>'][0][1],gless_cat['lat<deg>'][0][1]


[(plt.imshow(gl_mom0s[i],origin='lower',norm='asinh',
  extent=[gless_cat['lon<deg>'][i][0],gless_cat['lon<deg>'][i][1],
          gless_cat['lat<deg>'][i][0],gless_cat['lat<deg>'][i][1]]),
  plt.title(gless_cat['name'][i]),
#   plt.colorbar(orientation='horizontal'),
  plt.show())
 for i in range(7)]


plt.scatter(whole_gal_ds['lon<deg>'][whole_gal_ds['lat<deg>']<-20],
            whole_gal_ds['vel<km/s>'][whole_gal_ds['lat<deg>']<-20],
            c=whole_gal_ds['d1<kpc>'][whole_gal_ds['lat<deg>']<-20],norm='asinh')

plt.xlim(-130,-50)
plt.colorbar()


from copy import deepcopy


nu = deepcopy(whole_gal_ds)


nu['lon<deg>'] = whole_gal_ds['lon<deg>']-360


nu2 = deepcopy(whole_gal_ds)


nu2['lon<deg>'] = whole_gal_ds['lon<deg>']+360


gal_ds_long = pd.concat((nu2[nu2['lon<deg>']<360],whole_gal_ds,nu[nu['lon<deg>']>-360]))


[(plt.scatter(gal_ds_long['lon<deg>'],
            gal_ds_long['vel<km/s>'],
            c=gal_ds_long['d1<kpc>'],norm='asinh'),

plt.xlim(gless_cat['lon<deg>'][i][0],gless_cat['lon<deg>'][i][1]),
plt.colorbar(),
  plt.title(gless_cat['name'][i]),
#   plt.ylim(-30,30),
 plt.show()) for i in range(7)]


for i in range(7):
    plt.contourf(np.arange(gless_cat['lon<deg>'][i][0],gless_cat['lon<deg>'][i][1],-0.125),
                 np.arange(-100,101,1.3),
                 np.nansum(gls[i].data,axis=1),levels=np.arange(0,7e3,5e2))
    plt.title(gless_cat['name'][i])
    plt.ylim(-50,50)
#     plt.colorbar(orientation='horizontal')
    plt.show()


shorter_tiles = [gls[i].data[38:115,:,:] for i in range(7)]


from astrodendro import Dendrogram


dnds = [Dendrogram.compute(shorter_tiles[i],min_value=np.nanstd(shorter_tiles[i]),
                  min_delta=3*np.nanstd(shorter_tiles[i]),min_npix=500) for i in range(7)]


[len(dnds[i]) for i in range(7)]


def solm(thing):
    return (1.3*u.km/u.s*u.K*thing*hi_conv).to(u.solMass/u.pc**2)


for wo in st1dnd.trunk[0].descendants:
    if wo.is_branch:
        print(wo.idx,wo.descendants,'\n')


import bootleg as btlg


[btlg.dend_to_fits(dnds[i],'tile'+str(i+1)+'_dnd.fits') for i in range(7)]


[len(dnds[i].trunk) for i in range(7)]


[len(dnds[i].leaves) for i in range(7)]


dnds_chonk = [Dendrogram.compute(shorter_tiles[i],min_value=np.nanstd(shorter_tiles[i]),
                  min_delta=3*np.nanstd(shorter_tiles[i]),min_npix=3000) for i in range(7)]


[shorter_tiles[i].shape for i in range(7)]


[his[i].data.shape for i in range(17)]


np.nanmin(zcat['vdispersion_HI<km/s>'])


[shorter_tiles[i].size for i in range(7)]


5/1.3*3.6/0.125*3.6/0.125


dnds_taller_chonk = [Dendrogram.compute(shorter_tiles[i],min_value=np.nanstd(shorter_tiles[i]),
                  min_delta=5*np.nanstd(shorter_tiles[i]),min_npix=3000) for i in range(7)]


dnds_taller = [Dendrogram.compute(shorter_tiles[i],min_value=np.nanstd(shorter_tiles[i]),
                  min_delta=5*np.nanstd(shorter_tiles[i]),min_npix=500) for i in range(7)]


help(dnds_chonk[i])


dir(dnds_chonk[i].plotter())


(4e20*u.cm**-2 * (1.36*c.m_p)).to(u.solMass/u.pc**2)


((10*u.solMass/u.pc**2)/(1.36*c.m_p)).to(u.cm**-2)


gless_cat.add_row


QTable(np.array(([gless_cat[0]]*12 + [gless_cat[1]]*7 + [gless_cat[2]]*7 +\
[gless_cat[3]]*5 + [gless_cat[4]]*3 + [gless_cat[5]]*4 + [gless_cat[6]])))


gless_cat


(all in pix) l,b,r,v,


these = np.array(([250+360,50+450,50,35+38],[200+360,100+450,100,35+38],[550+360,100+450,75,38+38],\
[500+360,100+0,50,45+38],[1050+700,100,100,35+38],[1050+1500,150,80,34+38],\
[1050+400,100,100,45+38],[75+2200,115+150,100,40+38],[70+2200,115+50,100,43+38],\
[300+2200,100+150,100,50+38],[50+1700,50+150,50,38+38],[2150+20,500+20,30,40+38],\
[2150+150,500+20,50,40+38],[2150+350,500+30,100,40+38])).transpose()


gmcless_interm_cat = QTable(np.array((these[0]*-0.125 + 360, these[1]*0.125 - 50, these[2]*0.125, these[3]*1.3-100)).transpose(),
       names=['lon<deg>','lat<deg>','radius<deg>','vel<km/s>'])


gmcless_interm_cat.write('gless_interm_cat.fits',overwrite=True)


(147.5-360)/-0.125


(np.array([12.5 ,-18.75])+50)/0.125


(np.array([91.25, 35.0])-360)/-0.125


[(
plt.plot(np.nansum(l.get_mask()*shorter_tiles[i],axis=(1,2))),
  plt.show(),
  plt.imshow(np.nansum(l.get_mask()*shorter_tiles[i],axis=0),cmap=dc.cmap_min0('viridis')),
  plt.title(str(i)+str(l.idx)),
  plt.show())
#  print(i,l.idx))
 for i in range(7)
 for l in dnds_chonk[i].leaves
 if np.nanmax(solm(np.nansum(l.get_mask()*shorter_tiles[i],axis=0)).value) >=5

]


450:500,360:560


011,017,028,114,214,215,219,32,37,38,46,53,57,512,63


dnds_chonk[0][11]


gless_candidates = [dnds[i][l.idx]
 for i in range(7)
 for l in dnds_chonk[i].leaves
 if np.nanmax(solm(np.nansum(l.get_mask()*shorter_tiles[i],axis=0)).value) >=5

]


[(plt.plot(np.nansum(l.get_mask()*shorter_tiles[i],axis=(1,2))),
  plt.show(),
    plt.imshow(np.nansum(l.get_mask()*shorter_tiles[i],axis=0),cmap=dc.cmap_min0('viridis'),origin='lower'),
  plt.title(str(i)+str(l.idx)),
  plt.show())
 
 for i in range(7)
 for l in dnds_chonk[i].leaves
 
]


[(plt.imshow(np.nansum(l.get_mask()*shorter_tiles[i],axis=0),cmap=dc.cmap_min0('viridis'),origin='lower'),
  plt.title(str(i)+str(l.idx)),
  plt.show())
 
 for i in range(7)
 for l in dnds_taller[i].leaves
 
]


[(plt.imshow(np.nansum(l.get_mask()*shorter_tiles[i],axis=0),cmap=dc.cmap_min0('viridis'),origin='lower'),
  plt.title(str(i)+str(l.idx)),
  plt.show())
 
 for i in range(7)
 for l in dnds[i].leaves
 
]


#when the leaves go tho the edge of the tiile -> find the tile next to it and see if it can extent there


for i in range(7):
    plt.imshow(gl_mom0s[i],origin='lower')
    plt.title(gless_cat['name'][i])
    plt.colorbar(orientation='horizontal')
    plt.show()


for i in range(7):
    plt.imshow(gl_mom2s[i],origin='lower',
        extent=[gless_cat['lon<deg>'][i][0],gless_cat['lon<deg>'][i][1],
        gless_cat['lat<deg>'][i][0],gless_cat['lat<deg>'][i][1]])
    plt.title(gless_cat['name'][i])
    plt.colorbar(orientation='horizontal')
    plt.show()


for i in range(7):
    plt.imshow(gl_mom1s[i],origin='lower',norm='asinh',
        extent=[gless_cat['lon<deg>'][i][0],gless_cat['lon<deg>'][i][1],
        gless_cat['lat<deg>'][i][0],gless_cat['lat<deg>'][i][1]])
    plt.title(gless_cat['name'][i])
    plt.colorbar(orientation='horizontal')
    plt.show()


[(plt.imshow(gl_mom2s[i],origin='lower'),
  plt.title(gless_cat['name'][i]),
  plt.colorbar(orientation='horizontal'),
  plt.show())
 for i in range(7)]


gless_cat


gless_cat


[(plt.imshow(gl_mom0s[i],origin='lower'),
  plt.title(gless_cat['names'][i]),
  plt.colorbar(),
  plt.show())
 for i in range(17)]
