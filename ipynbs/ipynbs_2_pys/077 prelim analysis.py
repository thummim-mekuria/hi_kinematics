#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Ellipse

import datacube as dc

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

from photutils.aperture import CircularAnnulus

#------------------------------------------------------------
zcat = QTable(fits.open('final_cat.fits')[1].data)

hi4pi,dht = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


zcat['l_pix_GMC'][2] = ((zcat['lon_GMC<deg>'][2]-360)/(-0.125)).astype('int')


hi4pi,dht = fits.open('both_same_grid_720d.fits')


hi_conv = (1.3*u.km/u.s *1.828e18*((u.cm)**-2)/(u.K*(u.km/u.s)) \
           *1.36*c.m_p *u.K).to(u.solMass/u.pc**2) #last term is data

co_conv = (1.3*u.km/u.s *1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s))\
           *2*1.36*c.m_p *u.K).to(u.solMass/u.pc**2) #last term is data


px_size = 0.125*u.deg.to(u.rad)


hicubes = [hi4pi.data[zcat['v_pix_HI'][i][0]:zcat['v_pix_HI'][i][1],\
                      zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1],
                      zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]]
                      for i in range(17)]
whi   = [np.nansum(hicubes[i],axis=0) for i in range(17)]
sighi = [(whi[i] * hi_conv).to(u.solMass/u.pc**2) for i in range(17)]


cocubes = [dht.data[zcat['v_pix_GMC'][i][0]:zcat['v_pix_GMC'][i][1],\
                      zcat['b_pix_GMC'][i][0]:zcat['b_pix_GMC'][i][1],
                      zcat['l_pix_GMC'][i][0]:zcat['l_pix_GMC'][i][1]]
                      for i in range(17)]
wco   = [np.nansum(cocubes[i],axis=0) for i in range(17)]
siggmc = [(wco[i] * co_conv).to(u.solMass/u.pc**2) for i in range(17)]


lons = np.arange(5760)*(-0.125) + 360
lats = np.arange(800)*(0.125) - 50


for i in range(17):
    these_ls_hi = lons[zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]]
    these_bs_hi = lats[zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1]]
    
    these_ls_gmc = lons[zcat['l_pix_GMC'][i][0]:zcat['l_pix_GMC'][i][1]]
    these_bs_gmc = lats[zcat['b_pix_GMC'][i][0]:zcat['b_pix_GMC'][i][1]]
    
    
    plt.imshow(sighi[i].value,origin='lower',\
               extent=[these_ls_hi[0],these_ls_hi[-1],these_bs_hi[0],these_bs_hi[-1]])
    
    
    ctr = plt.contour(these_ls_hi,these_bs_hi,sighi[i].value,\
                 levels=[1,2,5,10,15,20,25,30,40,50,60,70,80,100])
    
    plt.contour(these_ls_gmc,these_bs_gmc,siggmc[i].value,\
                 levels=[1,5,10,15,20,25,30,50,70],colors=['orange'],\
                 linewidths=0.9) #might also want to explore logarithmic spacing in the contour levels
    
    fmt = {}
    for l in ctr.levels:
        fmt[l] = str(l)
    
    
    plt.clabel(ctr, ctr.levels[:-1], fmt=fmt,inline=True,fontsize=7)
    plt.scatter(np.nanmean(zcat['lon_GMC<deg>'][i]),np.nanmean(zcat['lat_GMC<deg>'][i]),
               color='orange',marker='x')
    
    plt.title(zcat['names'][i])
#     plt.xlim(these_ls[0],these_ls[-1])
    plt.show()


#for sigma as a function of r -- can be done pixel by pixel instead of averaging over annuli


c_coords_hi = [hicubes[i].shape[1]/2 for i in range(17)]
c_coords_co = [(cocubes[i].shape[1]/2,cocubes[i].shape[2]/2) for i in range(17)]


npix_totals_hi = np.array([np.nansum(sighi[i].shape[0]**2) for i in range(17)])
npix_totals_gmc = np.array([np.nansum(siggmc[i].shape[0]*siggmc[i].shape[1]) for i in range(17)])


all_px_ds_hi = np.array(([np.sqrt((i-c_coords_hi[k])**2 + (j-c_coords_hi[k])**2) 
                      for k in range(17) for i in range(sighi[k].shape[0]) 
                      for j in range(sighi[k].shape[1])]))


all_px_ds_gmc = np.array(([np.sqrt((i-c_coords_co[k][0])**2 + (j-c_coords_co[k][1])**2) 
                      for k in range(17) for i in range(siggmc[k].shape[0]) 
                      for j in range(siggmc[k].shape[1])]))


rarr_hi = [(all_px_ds_hi[
    np.nansum(npix_totals_hi[:i]):
    np.nansum(npix_totals_hi[:i])+npix_totals_hi[i]]).
    reshape(sighi[i].shape[0],-1) for i in range(17)]


rarr_gmc = [(all_px_ds_gmc[
    np.nansum(npix_totals_gmc[:i]):
    np.nansum(npix_totals_gmc[:i])+npix_totals_gmc[i]]).
    reshape(siggmc[i].shape[0],-1) for i in range(17)]


[print(rarr_gmc[i].shape,siggmc[i].shape) for i in range(17)]


from astropy.modeling import models, fitting


help(models.BrokenPowerLaw1D)


cal_init = models.Gaussian1D(1.,0.,1.)


fitter = fitting.LMLSQFitter()


def filter_zeros(x,y):
    finitex = x>0
    
    newx = x[finitex]
    newy = y[finitex]
    
    finitey = newy>0
    
    nnewx = newx[finitey]
    nnewy = newy[finitey]
    
    return np.array((nnewx,nnewy))
    

xn0,yn0 = filter_zeros(np.ravel(rarr_gmc[i])*0.125,np.ravel(siggmc[i].value))


fit


for i in range(17) :  
    x1,y1 = np.ravel(rarr_gmc[i])*0.125, np.ravel(siggmc[i].value)
#     x1,y1 = filter_zeros(np.ravel(rarr_gmc[i])*0.125, np.ravel(siggmc[i].value))


    x,y = np.log10(x1),np.log10(y1)
    plt.scatter(x1,y1,
              s=3e3/npix_totals_gmc[i])
    
    init = models.Gaussian1D(np.nanmax(y1),np.nanargmax(y1),np.nanstd(y1))

    fit = fitter(init,x1,y1)
    
    plt.plot(x1,fit(x1),color='red',linestyle='--')
    # plt.plot(np.log10(np.ravel(rarr_gmc[i])*0.125),(-0.75+np.log10(np.ravel(rarr_gmc[i])*0.125))*-45)

    plt.title(zcat['names'][i])
    plt.xlabel('log (radius [deg (change to pc later)])')
    plt.ylabel(r'log ($\Sigma\ [M_\odot\ pc^{-2}]$)')
#     plt.ylim(-2,4)
    plt.show()


#the fit is called log normal +
#power law tail (ie in not log space its a gaussian + exponential)


for i in range(17) :  
    x1,y1 = filter_zeros(np.ravel(rarr_gmc[i])*0.125, np.ravel(siggmc[i].value))

    x,y = np.log10(x1),np.log10(y1)
    plt.scatter(x,y,
              s=3e3/npix_totals_gmc[i])

    plt.plot(np.log10(x1),np.log10(fit(x1)),color='red',linestyle='--')
    # plt.plot(np.log10(np.ravel(rarr_gmc[i])*0.125),(-0.75+np.log10(np.ravel(rarr_gmc[i])*0.125))*-45)

    plt.title(zcat['names'][i])
    plt.xlabel('log (radius [deg (change to pc later)])')
    plt.ylabel(r'log ($\Sigma\ [M_\odot\ pc^{-2}]$)')
    plt.ylim(-2,4)
    plt.show()


for i in range(17):
#     plt.scatter(np.log10(np.ravel(rarr_hi[i])*0.125),
#               np.log10(np.ravel(sighi[i].value)),
#               s=3e3/npix_totals_hi[i])
    
    plt.scatter(np.log10(np.ravel(rarr_gmc[i])*0.125),
              np.log10(np.ravel(siggmc[i].value)),
              s=3e3/npix_totals_gmc[i])
    
    
    plt.title(zcat['names'][i])
    plt.xlabel('log (radius [deg (change to pc later)])')
    plt.ylabel(r'log ($\Sigma\ [M_\odot\ pc^{-2}]$)')
    plt.ylim(-3,)
    plt.show()


for i in range(17):
    rs_deg = np.ravel(rarr_hi[i])*0.125
    rs_pc = rs_deg*(u.deg.to(u.rad))*zcat['distance<pc>'][i]
    
    sigs = np.ravel(sighi[i].value)
    
    plt.scatter(np.log10(rs_pc), np.log10(sigs), s=3e3/npix_totals_hi[i])
    
#     plt.scatter(np.log10(np.ravel(rarr_gmc[i])*0.125),
#               np.log10(np.ravel(siggmc[i].value)),
#               s=3e3/npix_totals_gmc[i])
    plt.title(zcat['names'][i])
    plt.xlabel('log (radius [pc])')
    plt.ylabel(r'log ($\Sigma\ [M_\odot\ pc^{-2}]$)')
#   plt.loglog(),
    plt.show()


[(plt.plot(lats[zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1]],
           np.nanmean(sighi[i],axis=1)),
  
  plt.plot(lats[zcat['b_pix_HI'][i][0]:zcat['b_pix_HI'][i][1]],
           np.nanmedian(sighi[i],axis=1)),
  
  
  plt.xlabel('b [deg]'),
  plt.ylabel(r'$\Sigma_{HI}\ [M_\odot\ pc^{-2}]$'),
  plt.title(zcat['names'][i]),
            plt.show()) for i in range(17)]


[(plt.plot(lons[zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]],\
           np.nanmean(sighi[i],axis=0)),
  
  plt.plot(lons[zcat['l_pix_HI'][i][0]:zcat['l_pix_HI'][i][1]],\
           np.nanmedian(sighi[i],axis=0)),
  
  
  plt.xlabel('l [deg]'),
  plt.ylabel(r'$\Sigma_{HI}\ [M_\odot\ pc^{-2}]$'),
  plt.title(zcat['names'][i]),\
            plt.show()) for i in range(17)]


plt.scatter(np.log10(zcat['mass_HI<kg>']*u.kg.to(u.solMass)),
            zcat['mass_GMC<kg>']/zcat['mass_HI<kg>'],
            marker='+')

# plt.plot(np.arange(3.5,7,0.5),0.12*np.arange(3.5,7,0.5)-0.4,linestyle='--',alpha=0.5)
plt.xlabel(r'log ($M_{HI}\ [M_\odot]$)')
plt.ylabel(r'$M_{GMC}/M_{HI}$')

plt.show()


plt.scatter(zcat['radius_acc_HI<pc>'],
            zcat['radius_eff_GMC<pc>']/zcat['radius_acc_HI<pc>'],
            marker='+')

# plt.plot(np.arange(3.5,7,0.5),0.12*np.arange(3.5,7,0.5)-0.4,linestyle='--',alpha=0.5)
plt.xlabel(r'$R_{HI}\ [pc]$)')
plt.ylabel(r'$R_{GMC}/R_{HI}$')

plt.show()


plt.scatter(np.log10(zcat['mass_HI<kg>']*u.kg.to(u.solMass)),
            np.log10(zcat['radius_acc_HI<pc>']),
            marker='+')

plt.plot(np.arange(4.5,7.5,0.5),0.36*np.arange(4.5,7.5,0.5)-0.2,linestyle='--',alpha=0.5)
plt.xlabel(r'log ($M_{HI}\ [M_\odot]$)')
plt.ylabel(r'log ($R_{HI}\ [pc]$)')

plt.show()


plt.scatter(np.log10(zcat['mass_HI<kg>']*u.kg.to(u.solMass)),
            np.log10(zcat['mass_GMC<kg>']*u.kg.to(u.solMass)),
            marker='+')

plt.plot(np.arange(4.5,7.5,0.5),1.13*np.arange(4.5,7.5,0.5)-1.6,linestyle='--',alpha=0.5)
plt.xlabel(r'log ($M_{HI}\ [M_\odot]$)')
plt.ylabel(r'log ($M_{GMC}\ [M_\odot]$)')

plt.show()


plt.scatter(zcat['radius_acc_HI<pc>'],
            zcat['radius_eff_GMC<pc>'],
            marker='+')

plt.plot(np.arange(25,210),0.7*np.arange(25,210)-10,linestyle='--',alpha=0.5)
plt.xlabel(r'$R_{HI}\ [pc]$')
plt.ylabel(r'$R_{GMC}\ [pc]$')

plt.show()


plt.scatter(zcat['radius_acc_HI<pc>'],
            zcat['vdispersion_HI<km/s>'],
            marker='+')

# plt.plot(np.arange(4.5,7.5,0.5),1.13*np.arange(4.5,7.5,0.5)-1.6,linestyle='--',alpha=0.5)
plt.xlabel(r'$R_{HI}\ [pc]$')
plt.ylabel(r'$\Delta V_{HI}\ [km\ s^{-1}]$')
plt.show()


plt.scatter(zcat['vdispersion_HI<km/s>'],zcat['vdispersion_GMC<km/s>'],
            marker='+')

# plt.plot(np.arange(4.5,7.5,0.5),1.13*np.arange(4.5,7.5,0.5)-1.6,linestyle='--',alpha=0.5)
plt.xlabel(r'$\Delta V_{HI}\ [km\ s^{-1}]$')
plt.ylabel(r'$\Delta V_{GMC}\ [km\ s^{-1}]$')
plt.show()


plt.hist(zcat['radius_acc_HI<pc>'],bins=12,label='HI',color='grey')
plt.hist(zcat['radius_eff_GMC<pc>'],bins=12,fill=False,edgecolor='black',label='GMC')
plt.xlabel(r'$R_{cloud}\ [pc]$')

plt.legend()


plt.hist(np.log10(zcat['mass_HI<kg>']*u.kg.to(u.solMass)),
         bins=12,label='HI',color='grey')
plt.hist(np.log10(zcat['mass_GMC<kg>']*u.kg.to(u.solMass)),
         bins=12,fill=False,edgecolor='black',label='GMC')
plt.xlabel(r'$M_{cloud}\ [M_\odot]$')

plt.legend()


meanSigs = np.array([np.nanmean(sighi[i]).value for i in range(17)])
medSigs = np.array([np.nanmedian(sighi[i]).value for i in range(17)])
mxSigs = np.array([np.nanmax(sighi[i]).value for i in range(17)])


plt.scatter(meanSigs,np.log10(zcat['mass_GMC<kg>']*u.kg.to(u.solMass)), marker='+')

# plt.plot(np.arange(25,210),0.7*np.arange(25,210)-10,linestyle='--',alpha=0.5)
plt.xlabel(r'mean $\Sigma_{HI}\ [pc]$')
plt.ylabel(r'log ($M_{GMC}\ [M_\odot]$)')


plt.scatter(mxSigs,np.log10(zcat['mass_HI<kg>']*u.kg.to(u.solMass)), marker='+')

# plt.plot(np.arange(25,210),0.7*np.arange(25,210)-10,linestyle='--',alpha=0.5)
plt.xlabel(r'mean $\Sigma_{HI}\ [pc]$')
plt.ylabel(r'log ($M_{HI}\ [M_\odot]$)')


plt.scatter(np.nanmean(zcat['lat_GMC<deg>'],axis=1),meanSigs, marker='+')

plt.xlabel(r'b [deg]')
plt.ylabel(r'mean $\Sigma_{HI}\ [pc]$')


[(plt.hist(np.ravel(sighi[i].value),bins=np.arange(0,np.nanmax(sighi[i].value),5),
           log=True,label='HI',color='grey'),
  plt.hist(np.ravel(siggmc[i].value),bins=np.arange(0,np.nanmax(siggmc[i].value),5),
           log=True,fill=False,label='GMC'),
    plt.legend(),
    plt.xlabel(r'$\Sigma_{HI}\ [pc]$'),
    plt.ylabel('N'),
    plt.title(zcat['names'][i]),
    plt.show()) for i in range(17)]


[(plt.hist(np.ravel(sighi[i].value),bins=np.arange(0,np.nanmax(sighi[i].value),5),
           log=True,label=zcat['names'][i],edgecolor='xkcd:'+colors[i],fill=False,
          histtype='step',stacked=True),
  
#     plt.legend(),
    plt.xlabel(r'$\Sigma_{HI}\ [pc]$'),
    plt.ylabel('N'))
     for i in range(17)]


colors=['steel','carolina blue','flat blue','twilight','azure','greenblue','kiwi green','light seafoam',\
        'dark cream','sun yellow','orangeish','grapefruit','lipstick','bubblegum',\
        'deep magenta','light purple','liliac']


def vir(sig,r,m):
    return (5*(sig**2)*r)/(c.G*m)


gmc_virs = vir(zcat['vdispersion_GMC<km/s>']*u.km/u.s,zcat['radius_eff_GMC<pc>']*u.pc,zcat['mass_GMC<kg>']*u.kg).cgs


hi_virs = vir(zcat['vdispersion_HI<km/s>']*u.km/u.s,zcat['radius_acc_HI<pc>']*u.pc,zcat['mass_HI<kg>']*u.kg).cgs


plt.hist(hi_virs,color='grey')
plt.hist(gmc_virs,fill=False)
# plt.xscale('log')
plt.show()


#when i come back - -caluculate virial parameter , other properties ppl calculate in papers, 
#maybe create a utils.py or tools.py or aux.py that has all the functions i need that were between dc and cld
#eg automate mass calculation etc
#think about maybe doing the thing where i calculate some property at different sizes/bounds that define the cloud
#and see if there is a point at which things stabilizr -- ir something that points to a phyiscal boundary
#for this^ look at how lewis talks about the choice of wco for calculating alpha_co , 
#or nia calculating angular momentum vector as a function of size (highlighted pink)


[(plt.imshow(siggmc[i].value,origin='lower',cmap=dc.cmap_min0('viridis'),vmin=1),
  
plt.title(zcat['names'][i]),
plt.show()) for i in range(17)]


zcat.write('final_cat.fits',overwrite=True)
