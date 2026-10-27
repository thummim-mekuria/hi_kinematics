#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
from astropy.table import QTable
import matplotlib.pyplot as plt

import numpy as np
import datacube as dc
import astropy.units as u
import astropy.constants as c


zcat = QTable(fits.open('final_cat.fits')[1].data)


zcat.write('interm_cat2.fits')


hi4pi,dht = fits.open('both_same_grid_720d.fits')


l0s,l1s=(((np.array(zcat['lon<deg>'])-360)/(-0.125)).astype('int')).transpose()
b0s,b1s=(((np.array(zcat['lat<deg>'])+50)/(0.125)).astype('int')).transpose()


zcat['l_pix'] = np.array((l0s,l1s)).transpose()
zcat['b_pix'] = np.array((b0s,b1s)).transpose()


co_cubes_long = [dht.data[:,b0s[i]:b1s[i],l0s[i]:l1s[i]] for i in range(19)]


specs = [np.nansum(co_cubes_long[i],axis=(1,2)) for i in range(19)]


v0s,v1s = np.array([(np.where(specs[i]>155*0.2)[0][0],\
                     np.where(specs[i]>155*0.2)[0][-1]) for i in range(19)]).transpose()


zcat['vel<km/s>'] = np.array((v0s*1.3-100.,v1s*1.3-100.)).transpose()
zcat['v_pix'] = np.array((v0s,v1s)).transpose()


co_cubes = [dht.data[v0s[i]:v1s[i],b0s[i]:b1s[i],l0s[i]:l1s[i]] for i in range(19)]


co_conv = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) * 2 * 1.36 * c.m_p
hi_conv = 1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) * 1.36 * c.m_p
px_size = 0.125*u.deg.to(u.rad)


wco   = [np.nansum(co_cubes[i],axis=0) for i in range(19)]
sigco = [(wco[i] *u.K *1.3*u.km/u.s *co_conv).to(u.solMass/u.pc**2) for i in range(19)]
mco   = [(np.nansum(sigco[i]* (px_size*zcat['distance<pc>'][i]*u.pc)**2)).to(u.kg) for i in range(19)]


zcat['mass_GMC<kg>'] = mco


rds = [np.sqrt(((px_size*zcat['distance<pc>'][i]*u.pc)**2 * len(np.where(sigco[i]>0)[0]))/np.pi) for i in range(19)]
rads = [np.sqrt((0.125*u.deg)**2 * len(np.where(sigco[i]>0)[0])/np.pi) for i in range(19)]


zcat['radius_GMC<pc>'] = rds
zcat['radius_GMC<deg>'] = rads


from astropy.modeling import models, fitting


def get_dispers (spec):
    fitter = fitting.LevMarLSQFitter()
#     print(np.nanmax(spec),np.nanargmax(spec),5)
    
#     ginit = models.Gaussian1D(1.,0.,1.)
    ginit = models.Gaussian1D(np.nanmax(spec),np.nanargmax(spec),5.)

    gfit = fitter(ginit,np.arange(len(spec)),spec)
#     print(gfit.parameters)
    return gfit


[(plt.plot(specs[i]),\
 plt.plot(get_dispers(specs[i])(np.arange(len(specs[i])))),\
  plt.title(zcat['names'][i]),
  plt.show(),
#   print(np.nanstd(specs[i])),
    
 ) for i in range(19)]


zcat['vdispersion_GMC<km/s>']=[get_dispers(specs[i]).stddev.value*1.3*u.km/u.s for i in range(19)]


zcat.remove_row(0)


zcat.remove_row(1)


zcat.write('final_cat.fits',overwrite=True)


for i in range(19):
    ii = np.array((v0s*1.3-100.,v1s*1.3-100.)).transpose()[i]
#     print(ii[0],ii[1])
    plt.hlines(zcat['names'][i],ii[0],ii[1],color='red',linestyle='--')
    
#     iii = zcat['vel<km/s>'][i]
#     plt.hlines(zcat['names'][i],iii[0],iii[1],color='orange',linestyle=':')

    
plt.vlines([-20,-10,0,10,20,30,40],'aquila_r','taurus',alpha=0.5)
plt.vlines([-15,-5,5,15,25,35],'aquila_r','taurus',alpha=0.5,linestyle='--')

plt.show()


lbs = [(np.nansum(hi4pi.data[v0s[i]:v1s[i],:,:],axis=0)*hi_conv*u.K*1.3*u.km/u.s,\
        np.nansum(dht.data[v0s[i]:v1s[i],:,:],axis=0)*co_conv*u.K*1.3*u.km/u.s) for i in range(19)]


lvs = [(np.nansum(hi4pi.data[:,b0s[i]:b1s[i],:],axis=1)*hi_conv*u.K*0.125*u.km/u.s,\
        np.nansum(dht.data[:,b0s[i]:b1s[i],:],axis=1)*co_conv*u.K*0.125*u.km/u.s) for i in range(19)]


bvs = [(np.nansum(hi4pi.data[:,:,l0s[i]:l1s[i]],axis=2)*hi_conv*u.K*0.125*u.km/u.s,\
        np.nansum(dht.data[:,:,l0s[i]:l1s[i]],axis=2)*co_conv*u.K*0.125*u.km/u.s) for i in range(19)]


lons = np.arange(5761)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50
vels = np.arange(155)*1.3 - 100


#going to drop aquila & camelopardalis


[(plt.title(zcat['names'][i]),\
  plt.imshow(lbs[i][0].to(u.solMass/u.pc**2).value,
             origin='lower',
             cmap=dc.cmap_min0('viridis'),
             extent=[360,-360,-50,50], 
             vmax=15),
  
#     plt.xlim(zcat['lon<deg>'][i]),
#     plt.ylim(zcat['lat<deg>'][i]),
    plt.xlim(360+(-0.125)*np.array((l0s[i]-50,l1s[i]+50))),
    plt.ylim(-50+(0.125)*np.array((b0s[i]-50,b1s[i]+50))),

    plt.contour(lons,lats,
                lbs[i][1].to(u.solMass/u.pc**2).value,
                levels=[10,30,50],cmap='binary_r',linewidths=[0.8]), 
    plt.show(),
  
  plt.imshow(lvs[i][0].to(u.solMass/u.pc**2).value,\
         origin='lower',
         cmap=dc.cmap_min0('viridis'),
         extent=[360,-360,-100,100], \
         vmax=2),\
  
#     plt.xlim(360+(-0.125)*np.array((l0s[i]-100,l1s[i]-100))),\
#     plt.ylim(-100+(1.3)*np.array((v0s[i]-15,v1s[i]-15))),\
#     plt.colorbar(),\
    plt.contour(lons,vels,
                lvs[i][1].to(u.solMass/u.pc**2).value,
                levels=[1e-2,0.1,0.5,1,2,5],cmap='binary_r',linewidths=0.8), 
  
    plt.show(),\
  
#   plt.imshow(bvs[i],
#              origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-50,50,-100,100], \
#              vmax=1500),\
  
#           plt.xlim(-50+(0.125)*np.array((b0s[i]-50,b1s[i]-50))),\
#           plt.ylim(-100+(1.3)*np.array((v0s[i]-8,v1s[i]-8))),\
# #   plt.colorbar(),
#   plt.contour(lats,vels,\
#              bvs[i],\
#              levels=[300,500,700,900]), 
  
 
 plt.show()) for i in np.concatenate(([1],np.arange(3,19)))]


[(plt.title(zcat['names'][i]),\
#   plt.imshow(lbs[i]*(u.K *1.3*u.km/u.s*hi_conv).to(u.solMass/u.pc**2).value,\
#              origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50], \
#              vmax=15),\
# #     plt.xlim(360+(-0.125)*np.array((l0s[i]-50,l1s[i]-50))),\
# #     plt.ylim(-50+(0.125)*np.array((b0s[i]-50,b1s[i]-50))),\
# #     plt.colorbar(),\
#     plt.contour(lons,lats,\
#              lbs[i]*(u.K *1.3*u.km/u.s*hi_conv).to(u.solMass/u.pc**2).value,\
#              levels=[1,3]), 
  
#     plt.show(),
#   print(zcat['lat<deg>'][i]),\

#   plt.imshow(lvs[i],\
#          origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100], \
#          vmax=2500),\
# #     plt.xlim(360+(-0.125)*np.array((l0s[i]-50,l1s[i]-50))),\
# #     plt.ylim(-100+(1.3)*np.array((v0s[i]-8,v1s[i]-8))),\
# #     plt.colorbar(),\
#     plt.contour(lons,vels,lvs[i],levels=[1000,1500,2000,3000]), 
  
#     plt.show(),\
  
 
  print(zcat['lon<deg>'][i]),\
  plt.imshow(bvs[i],
             origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-50,50,-100,100]), \
#              vmax=1500),\
  
#     plt.xlim(-50+(0.125)*np.array((b0s[i]-50,b1s[i]-50))),\
#     plt.ylim(-100+(1.3)*np.array((v0s[i]-8,v1s[i]-8))),\
#   plt.colorbar(),
    plt.contour(lats,vels, bvs[i],levels=[500],colors='white'),#,500,700,900]), 
  
 
 plt.show()) for i in range(19)]


bvranges = np.array(([-50,50],[-50,40],[-50,40],[-50,50],[-25,20],[-25,20],[-50,25],\
              [-50,50],[-25,15],[-40,15],[-50,30],[-30,10], [-40,50],[-40,25],\
             [-50,30],[-50,25],[-50,30],[-25,10],[-50,30])).transpose()


vranges = np.array(([-75,100],[-75,20],[-90,10],[-70,100],[-5,90],[-5,75],\
[-20,30],[-80,80],[-100,5],[-5,75],[-10,50],[-5,75],[-5,75],\
[-10,50],[-25,40],[-75,5],[-100,5],[0,60],[-60,20])).transpose()


plt.hlines(zcat['names'],bvranges[0],bvranges[1],color='orange',linestyle='--')
   
plt.hlines(zcat['names'],zcat['lat<deg>'].transpose()[0],zcat['lat<deg>'].transpose()[1],\
           color='black',linestyle=':')

plt.show()


plt.hlines(zcat['names'],vranges[0],vranges[1],color='orange',linestyle='--')
   
plt.hlines(zcat['names'],zcat['vel<km/s>'].transpose()[0],zcat['vel<km/s>'].transpose()[1],\
           color='black',linestyle=':')

plt.show()


[(plt.title(zcat['names'][i]),\
#   plt.imshow(lbs[i]*(u.K *1.3*u.km/u.s*hi_conv).to(u.solMass/u.pc**2).value,\
#              origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-50,50], \
#              vmax=15),\
# #     plt.xlim(360+(-0.125)*np.array((l0s[i]-50,l1s[i]-50))),\
# #     plt.ylim(-50+(0.125)*np.array((b0s[i]-50,b1s[i]-50))),\
# #     plt.colorbar(),\
#     plt.contour(lons,lats,\
#              lbs[i]*(u.K *1.3*u.km/u.s*hi_conv).to(u.solMass/u.pc**2).value,\
#              levels=[1,3]), 
  
#     plt.show(),
  print(zcat['lat<deg>'][i]),\

  plt.imshow(lvs[i],\
         origin='lower',cmap=dc.cmap_min0('viridis'),extent=[360,-360,-100,100], \
         vmax=2500),\
#     plt.xlim(360+(-0.125)*np.array((l0s[i]-50,l1s[i]-50))),\
#     plt.ylim(-100+(1.3)*np.array((v0s[i]-8,v1s[i]-8))),\
#   plt.ylim(vranges[0][i],vranges[1][i]),
  plt.vlines([230,100,230-360,100-360],-100,100,color='red'),
  plt.vlines([260,130,260-360,130-360],-100,100,color='red',linestyles='--'),

#     plt.colorbar(),\
    plt.contour(lons,vels,lvs[i],levels=[500,1000],colors=['white'],linewidths=[0.7,1]), 
    plt.xlim(360,0),
    plt.show(),\
  
 
#   print(zcat['lon<deg>'][i]),\
#   plt.imshow(bvs[i],
#              origin='lower',cmap=dc.cmap_min0('viridis'),extent=[-50,50,-100,100]), \
# #              vmax=1500),\
  
# #     plt.xlim(-50+(0.125)*np.array((b0s[i]-50,b1s[i]-50))),\
# #     plt.ylim(-100+(1.3)*np.array((v0s[i]-8,v1s[i]-8))),\
# #   plt.colorbar(),
#     plt.contour(lats,vels, bvs[i],levels=[500],colors='white'),#,500,700,900]), 
  
 
 plt.show()) for i in range(19)]


#when i come back -- figure out some sort of way to tell the program what 
#v range to integrate over, as a function of l & b (+ vice versa, ie what range of l for a given b/v etc)
#also, come back and convert wcointo  units of solmass/pc^2 to see at what threshold the CO closed contours are
#then possibly use that to determine the surface density threshold for the hi?
#for  the v range as a function of l ... do i want to pick the 'comomn/most frequented' v raanges, or the 
#'awayf from the data peaks' v ranges when looking for gmc-less? ... or, 
#first find them, then see how their v ranges compare to the global trend?
#also, what is the equivalent to 10 msol/pc^2 when the map is bv or lv? sth to do with *0.125 deg? 
#but then will xco work? or wilson/other people might have a way?
#alsoo, see where the gmcs are relative to the global lv trend around them


lvgal = np.nansum(hi4pi.data,axis=1)


bvgal = np.nansum(hi4pi.data,axis=2)


plt.figure(figsize=(20,7))
plt.imshow(bvgal,extent=[-50,50,-100,100],origin='lower',vmax=1e5)
plt.vlines(zcat['lon<deg>'].transpose()[0],bvranges[0],bvranges[1],colors='white')
# plt.vlines(zcat['lon<deg>'].transpose()[1],vranges[0],vranges[1],colors='white')
# for i in range(19):
#     plt.hlines(vranges[0][i],zcat['lon<deg>'][i][0],zcat['lon<deg>'][i][1],colors='white')
#     plt.hlines(vranges[1][i],zcat['lon<deg>'][i][0],zcat['lon<deg>'][i][1],colors='white')

# plt.vlines([-45,45],-100,100,colors='red')

# plt.xlim(180,-180)


plt.figure(figsize=(20,7))
plt.imshow(lvgal,extent=[360,-360,-100,100],origin='lower',vmax=10000)
plt.vlines(zcat['lon<deg>'].transpose()[0],vranges[0],vranges[1],colors='white')
plt.vlines(zcat['lon<deg>'].transpose()[1],vranges[0],vranges[1],colors='white')
for i in range(19):
    plt.hlines(vranges[0][i],zcat['lon<deg>'][i][0],zcat['lon<deg>'][i][1],colors='white')
    plt.hlines(vranges[1][i],zcat['lon<deg>'][i][0],zcat['lon<deg>'][i][1],colors='white')

plt.vlines([-45,45],-100,100,colors='red')

plt.xlim(180,-180)
