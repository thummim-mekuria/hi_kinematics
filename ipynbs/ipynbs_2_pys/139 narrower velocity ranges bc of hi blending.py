#!/usr/bin/env python
# coding: utf-8


one = [fits.open('../rn/h2_0thmom_masks_DISPLAY.fits')[i].data-1 for i in range(22)]
two = [fits.open('../rn/h2_0thmom_masks.fits')[i].data-1 for i in range(22)]


[(np.nansum(two[i])/one[i].size) for i in range(22)]


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Ellipse
from matplotlib.colors import ListedColormap
from matplotlib import cm

import numpy as np
import pandas as pd

import datacube as dc
import fun as fns

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

from scipy.optimize import curve_fit


hi4pi,dht = fits.open('../rn/both_same_grid_720d.fits')


oldcat = QTable(fits.open('../rn/sighhplease.fits')[1].data)


coco = fits.open('../rn/h2_h2sized.fits')


h2_msks = [fits.open('../rn/h2_0thmom_masks_DISPLAY.fits')[i].data-1 for i in range(22)]


new_lbs = np.array([
       [  0,  64,   5,  95], #0
       [  5,  95,   0, 155],
       [  0,  24,   0,  48], #2
       [  5,  64,   2, 100],
       [  4,  36,   0,  44], #4
       [  0,  55,   3,  42],
       [  3,  47,   0,  32], #6
       [  0,  75,   0,  80],
       [  0,  31,   0,  64], #8
       [  10, 72,   4,  25],
       [  0, 150,   0, 136], #10
       [  0,  33,   3,  50],
       [  0,  32,   0,  56], #12
       [  5,  38,   2,  66],
       [  0, 178,   0, 170], #14
       [  20, 60,  10, 120],
       [  5,  68,   5,  72], #16
       [  10, 100,   0, 116],
       [  0, 100,   0,  56], #18
       [  3,  70,   8,  75],
       [  0,  24,   0,  32], #20
       [  0,  75,   0, 102]])


new_coco = [coco[i].data[:,new_lbs[i][0]:new_lbs[i][1],new_lbs[i][2]:new_lbs[i][3]]
for i in range(22)]

new_h2_msks = [h2_msks[i][new_lbs[i][0]:new_lbs[i][1],new_lbs[i][2]:new_lbs[i][3]].astype('bool')
for i in range(22)]


new_h2_mass = np.array([fns.mass_n_area(new_coco[i]*u.K,
                 oldcat['distance<pc>'][i]*u.pc,
                 'CO',new_h2_msks[i])[0].to(u.solMass).value for i in range(22)])


new_h2_rad = np.array([np.sqrt(fns.mass_n_area(new_coco[i]*u.K,
                 oldcat['distance<pc>'][i]*u.pc,
                 'CO',new_h2_msks[i])[1]/np.pi).value for i in range(22)])


new_hi_rad = (((3*new_h2_mass*u.solMass)/(4*np.pi*0.57*u.cm**-3*1.36*c.m_p))**(1/3)).to(u.pc).value


boc_etal = np.array([new_lbs[i]+np.array((oldcat['boc'][i],oldcat['boc'][i],
                     oldcat['loc'][i],oldcat['loc'][i])) for i in range(22)])


zc = [(
    np.sin(np.nanmean(boc_etal[i][:2]*0.125-50)*u.deg).value*oldcat['distance<pc>'][i])
 for i in range(22)]


np.nanmean(np.abs(zc)),np.nanmedian(np.abs(zc))


rad_hi_px = ((new_hi_rad/oldcat['distance<pc>'])*u.rad.to(u.degree)/0.125).astype('int')


bc,lc = np.nanmean(boc_etal.reshape(22,-1,2),axis=2).transpose().astype('int')


new_hihi = [hi4pi.data[:,bc[i]-rad_hi_px[i]:bc[i]+rad_hi_px[i],
             lc[i]-rad_hi_px[i]:lc[i]+rad_hi_px[i]] for i in range(22)]


[(plt.contourf(np.log10(np.nansum(new_hihi[i],axis=0)*1.3*1.82e18),levels=27,cmap='magma'),
  plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


[(plt.contourf(np.log10(np.nansum(new_hihi[i],axis=1)*0.125*1.82e18),
               levels=np.linspace(20,21.5,20),cmap='magma',
              ),
  plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


[(plt.contourf(np.log10(np.nansum(new_hihi[i],axis=2)*0.125*1.82e18),
               levels=np.linspace(20,21.5,20),cmap='magma'),
  plt.gca().set_aspect(1),
plt.show()
) for i in range(22)]


#after the break, transcribe my description of my method into the agenda
#add other agenda items (using spheres instead of cylinder, draw a little projection-y diagram to explain)
#hmm


plt.rcParams['font.size']=15


for i in range(22):
    
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,7))
    lv_hi = np.log10(np.nansum(new_hihi[i],axis=1)*0.125*1.82e18)
    bv_hi = np.log10(np.nansum(new_hihi[i],axis=2)*0.125*1.82e18)
    
    lv_h2 = np.log10(np.nansum(new_coco[i],axis=1)*0.125*1.97e20)
    bv_h2 = np.log10(np.nansum(new_coco[i],axis=2)*0.125*1.97e20)
    
    
    axs[0].contourf(lv_hi,
                   levels=np.linspace(20,np.nanmax(lv_hi),20),
                    cmap='magma',
                   extent=[(lc[i]-rad_hi_px[i])*-0.125+360,
                           (lc[i]+rad_hi_px[i])*-0.125+360,
                       -100,100]),
    
    axs[0].contourf(lv_h2,
                   levels=np.linspace(20,np.nanmax(lv_h2),10),
                    cmap='YlGnBu_r',
                   extent=[(boc_etal[i][2])*-0.125+360,
                           (boc_etal[i][3])*-0.125+360,
                       -100,100]),
    
    
#     axs[0].vlines(lc[i]*-0.125+360,-100,100)
    

    axs[1].contourf(bv_hi,
                   levels=np.linspace(19.5,np.nanmax(bv_hi),20),
                    cmap='magma',
                   extent=[(bc[i]-rad_hi_px[i])*0.125-50,
                    (bc[i]+rad_hi_px[i])*0.125-50,
                       -100,100]),
    
    
    axs[1].contourf(bv_h2,
                   levels=np.linspace(19.5,np.nanmax(bv_h2),10),
                    cmap='YlGnBu_r',
                   extent=[(boc_etal[i][0])*0.125-50,
                    (boc_etal[i][1])*0.125-50,
                       -100,100]),
    
#     axs[1].vlines(bc[i]*0.125-50,-100,100)
    

    axs[2].plot(np.nanmean(new_hihi[i].data,axis=(1,2)),np.arange(155)*1.3-100)
    axs[2].plot(100*np.nanmean(new_coco[i].data,axis=(1,2)),np.arange(155)*1.3-100)


    fig.supylabel('velocity [km/s]')
    axs[0].set_xlabel('longitude')
    axs[0].invert_xaxis()
    axs[1].set_xlabel('latitude')
    axs[2].set_xlabel('$T_B$ [K]')
    
    for ax in axs[1:]:
        ax.set_yticks([])

    fig.suptitle(oldcat['name'][i])
    plt.subplots_adjust(wspace=0.05)

    fig.tight_layout()

    plt.show()


status report for friday
pick a method for fitting a gaussian to the spectrum
spheres


from gaussdecomp import spectrum, fitter


hi_sps = [spectrum.Spectrum(
            np.nanmean(new_hihi[i],axis=(1,2)),
            np.arange(155)*1.3-100)
          
          for i in range(22)]


from astropy.modeling import models


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


hi_sps_gauspars


#fix the chosen ones, then apply those same velocity ranges to the three
#in notebook 141, then first moment maps & fit gradients and see if there is a big difference
#in those results
#also redo mass calculation with the narrower velocities and see if that makes a difference


chosen = [1,6,2,2,4,3,4,1,1,3,1,4,1,3,2,1,3,1,0,1,1,3]


for i in range(12,22):
    As,xos,sigs = hi_sps_gauspars[i]
#     print(np.abs(xos-fns.get_mom1(np.nanmean(coco[i].data,axis=(1,2)),np.arange(155)*1.3-100)))
    
    a,xo,sig = np.array([As,xos,sigs]).transpose()[chosen[i]]
    plt.plot(np.arange(155)*1.3-100,
             models.Gaussian1D(a,xo,sig)(np.arange(155)*1.3-100),
            color='red')
    
    for a,xo,sig in zip(As,xos,sigs):
        plt.plot(np.arange(155)*1.3-100,
             models.Gaussian1D(a,xo,sig)(np.arange(155)*1.3-100),
            linewidth=0.5,color='grey')
        
        
    plt.plot(np.arange(155)*1.3-100,hi_sps[i].flux,color='red',linewidth=0.3)
    
    plt.plot(np.arange(155)*1.3-100,50*np.nanmean(new_coco[i],axis=(1,2)),
             color='black')
    
    plt.xlabel('$v\mathrm{_{LSR}\ [km\ s^{-1}]}$')
    plt.title(oldcat['name'][i])
#     plt.xlim(-20,100)
    
    plt.show()


# for the angular momentum change -- look into torques

# iterate the decomposition -- for the locations closer to the galactic plane, decompose the spectra again, except restricted to the velocity range where the previous peak is
# 
# a sanity check would be to look at the clouds that are way off the galactic plane, and see how those gradients compare to what we get out of defining a window for the clouds closer to the galactic plane
# another sanity check - check against whats physically possible (what is reasonable timescales for the gas to accumulate given the range of velocities


for i in range(22):
    
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,7))
    lv_hi = np.log10(np.nansum(new_hihi[i],axis=1)*0.125*1.82e18)
    bv_hi = np.log10(np.nansum(new_hihi[i],axis=2)*0.125*1.82e18)
    
    lv_h2 = np.log10(np.nansum(new_coco[i],axis=1)*0.125*1.97e20)
    bv_h2 = np.log10(np.nansum(new_coco[i],axis=2)*0.125*1.97e20)
    

    axs[0].contourf(lv_hi,
                   levels=np.linspace(20,np.nanmax(lv_hi),20),
                    cmap='magma',
                   extent=[(lc[i]-rad_hi_px[i])*-0.125+360,
                           (lc[i]+rad_hi_px[i])*-0.125+360,
                       -100,100]),
    
    axs[0].contourf(lv_h2,
                   levels=np.linspace(20,np.nanmax(lv_h2),10),
                    cmap='YlGnBu_r',
                   extent=[(boc_etal[i][2])*-0.125+360,
                           (boc_etal[i][3])*-0.125+360,
                       -100,100]),
    
    axs[0].hlines(noo3vs[i],(lc[i]-rad_hi_px[i])*-0.125+360,
                           (lc[i]+rad_hi_px[i])*-0.125+360)
    
    
#     axs[0].vlines(lc[i]*-0.125+360,-100,100)
    

    axs[1].contourf(bv_hi,
                   levels=np.linspace(19.5,np.nanmax(bv_hi),20),
                    cmap='magma',
                   extent=[(bc[i]-rad_hi_px[i])*0.125-50,
                    (bc[i]+rad_hi_px[i])*0.125-50,
                       -100,100]),
    
    
    axs[1].contourf(bv_h2,
                   levels=np.linspace(19.5,np.nanmax(bv_h2),10),
                    cmap='YlGnBu_r',
                   extent=[(boc_etal[i][0])*0.125-50,
                    (boc_etal[i][1])*0.125-50,
                       -100,100]),
    
    axs[1].hlines(noo3vs[i],(bc[i]-rad_hi_px[i])*0.125-50,
                    (bc[i]+rad_hi_px[i])*0.125-50)
    
#     axs[1].vlines(bc[i]*0.125-50,-100,100)
    

    axs[2].plot(chosen_sps[i],np.arange(155)*1.3-100)
    axs[2].plot(100*np.nanmean(new_coco[i].data,axis=(1,2)),np.arange(155)*1.3-100)


    fig.supylabel('velocity [km/s]')
    axs[0].set_xlabel('longitude')
    axs[0].invert_xaxis()
    axs[1].set_xlabel('latitude')
    axs[2].set_xlabel('$T_B$ [K]')
    
    for ax in axs[1:]:
        ax.set_yticks([])
        
        
    for ax in axs:
        ax.set_ylim(noo3vs[i]+20*np.array([-1,1]))

    fig.suptitle(oldcat['name'][i])
    plt.subplots_adjust(wspace=0.05)

    fig.tight_layout()

    plt.show()


chosen_sps = [models.Gaussian1D(a,xo,sig)(np.arange(155)*1.3-100)
     
    for i in range(22)
    for a,xo,sig 
in [hi_sps_gauspars[i].transpose()[chosen[i]]]

]


[(plt.plot(chosen_sps[i]/np.nanmax(chosen_sps[i])),
  plt.vlines([np.where(chosen_sps[i]>0.13)[0][0],np.where(chosen_sps[i]>0.13)[0][-1]],
             0,1,color='black'),
  
  plt.plot(np.nanmean(new_coco[i],axis=(1,2))/np.nanmax(np.nanmean(new_coco[i],axis=(1,2))),
             color='orange'),
#   plt.title(ca)
  
  plt.xlim(np.where(chosen_sps[i]>0.1)[0][0]-20,np.where(chosen_sps[i]>0.1)[0][-1]+20),
 plt.show())
 for i in range(22)]


noonoovs = np.array([
    [np.where(chosen_sps[i]>0.1)[0][0],np.where(chosen_sps[i]>0.1)[0][-1]]
    for i in range(22)
])*1.3-100


noonoovs


[(
    plt.plot(np.arange(155)*1.3-100,np.nanmean(new_hihi[i],axis=(1,2))),
    plt.plot(np.arange(155)*1.3-100,50*np.nanmean(new_coco[i],axis=(1,2))),
    plt.vlines(noonoovs[i],0,100),
plt.show())
    for i in range(22)]


#just doing 4sigma
noo3vs = [(xo-2.*sig,xo+2.*sig)
     
    for i in range(22)
    for a,xo,sig 
in [hi_sps_gauspars[i].transpose()[chosen[i]]]

]


noo3vs


for i in range(22):
    
    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,7))
    lv_hi = np.log10(np.nansum(new_hihi[i],axis=1)*0.125*1.82e18)
    bv_hi = np.log10(np.nansum(new_hihi[i],axis=2)*0.125*1.82e18)
    
    lv_h2 = np.log10(np.nansum(new_coco[i],axis=1)*0.125*1.97e20)
    bv_h2 = np.log10(np.nansum(new_coco[i],axis=2)*0.125*1.97e20)
    
    
    axs[0].contourf(lv_hi,
                   levels=np.linspace(20,np.nanmax(lv_hi),20),
                    cmap='magma',
                   extent=[(lc[i]-rad_hi_px[i])*-0.125+360,
                           (lc[i]+rad_hi_px[i])*-0.125+360,
                       -100,100]),
    
    axs[0].contourf(lv_h2,
                   levels=np.linspace(20,np.nanmax(lv_h2),10),
                    cmap='YlGnBu_r',
                   extent=[(boc_etal[i][2])*-0.125+360,
                           (boc_etal[i][3])*-0.125+360,
                       -100,100]),
    axs[0].hlines(noonoovs[i],(lc[i]-rad_hi_px[i])*-0.125+360,
                          (lc[i]+rad_hi_px[i])*-0.125+360),
    
    axs[0].hlines(noo3vs[i],(lc[i]-rad_hi_px[i])*-0.125+360,
                          (lc[i]+rad_hi_px[i])*-0.125+360,linestyle='--'),
    
#     axs[0].vlines(lc[i]*-0.125+360,-100,100)
    

    axs[1].contourf(bv_hi,
                   levels=np.linspace(19.5,np.nanmax(bv_hi),20),
                    cmap='magma',
                   extent=[(bc[i]-rad_hi_px[i])*0.125-50,
                    (bc[i]+rad_hi_px[i])*0.125-50,
                       -100,100]),
    
    
    axs[1].contourf(bv_h2,
                   levels=np.linspace(19.5,np.nanmax(bv_h2),10),
                    cmap='YlGnBu_r',
                   extent=[(boc_etal[i][0])*0.125-50,
                    (boc_etal[i][1])*0.125-50,
                       -100,100]),
    axs[1].hlines(noonoovs[i],(bc[i]-rad_hi_px[i])*0.125-50,
                          (bc[i]+rad_hi_px[i])*0.125-50),
    
    axs[1].hlines(noo3vs[i],(bc[i]-rad_hi_px[i])*0.125-50,
                          (bc[i]+rad_hi_px[i])*0.125-50,linestyle='--'),
    

#     axs[1].vlines(bc[i]*0.125-50,-100,100)
    

    axs[2].plot(np.nanmean(new_hihi[i].data,axis=(1,2)),np.arange(155)*1.3-100)
    axs[2].plot(100*np.nanmean(new_coco[i].data,axis=(1,2)),np.arange(155)*1.3-100)
    axs[2].hlines(noonoovs[i],0,100),
    axs[2].hlines(noo3vs[i],0,100,linestyle='--'),


    fig.supylabel('velocity [km/s]')
    axs[0].set_xlabel('longitude')
    axs[0].invert_xaxis()
    axs[1].set_xlabel('latitude')
    axs[2].set_xlabel('$T_B$ [K]')
    
    for ax in axs[1:]:
        ax.set_yticks([])

    fig.suptitle(oldcat['name'][i])
    plt.subplots_adjust(wspace=0.05)

    fig.tight_layout()

    plt.show()


noo3vo,noo3vi = np.array((((np.array(noo3vs)+100)/1.3).transpose()[0].astype('int'),
np.around(((np.array(noo3vs)+100)/1.3).transpose()[1])))


for i in range(22):
    plt.figure(figsize=(4,4))
    plt.contourf(np.log10(np.nansum(
        new_hihi[i][int(noo3vo[i]):int(noo3vi[i]),:,:],
        axis=0)*1.3*1.82e18),levels=np.arange(20.6,22,0.1),cmap='magma')
    plt.colorbar()
    plt.gca().set_aspect(1)
    plt.show()


newcat = QTable.from_pandas(pd.DataFrame(boc_etal,columns=['boc','bic','loc','lic']))


newcat['loh'] = lc - rad_hi_px
newcat['lih'] = lc + rad_hi_px
newcat['boh'] = bc - rad_hi_px
newcat['bih'] = bc + rad_hi_px


noo3vo,noo3vi


newcat['vo'] = noo3vo.astype('int')
newcat['vi'] = noo3vi.astype('int')


rh_over_d = new_hi_rad/oldcat['distance<pc>']


100*(np.arctan(rh_over_d)*u.rad.to(u.deg)/(np.arctan(rh_over_d*np.sqrt(2))*u.rad.to(u.deg)))


np.sin(np.nanmean(oldcat['lat_GMC<deg>'],axis=1)*u.deg).value


# for a cylinder with height h and diameter D, seeing it from an angle o b=0 we would just see the lateral side (which has dimensions of h x D) and from an angle of b=90, we would just see just the circular face, which has a radius of r. from any other angle=


new_hi_rad


50/140


def angle_cyl (R,d,b):
    
    A = np.sqrt(d**2 + 2*R**2 - 2*np.sqrt(2)*R*d*np.cos((45+b)*u.deg))
    alpha = np.arcsin(np.sqrt(2)*R*np.sin((45+b)*u.deg)/A)
    
    
    B = np.sqrt(d**2 + 2*R**2 - 2*np.sqrt(2)*R*d*np.cos((135-b)*u.deg))
    beta = np.arcsin(np.sqrt(2)*R*np.sin((135-b)*u.deg)/B)
    
    return (alpha+beta).to(u.deg).value


theta_cyl = np.array([angle_cyl(new_hi_rad[i],oldcat['distance<pc>'][i],
                       np.nanmean(oldcat['lat_GMC<deg>'][i])) for i in range(22)])


theta_crc = 2*np.arctan(new_hi_rad/oldcat['distance<pc>'])*u.rad.to(u.deg)


scene2 = (np.abs(oldcat['distance<pc>']*np.sin(np.nanmean(oldcat['lat_GMC<deg>'],axis=1)*u.deg)) < new_hi_rad)


# plt.scatter(new_hi_rad/oldcat['distance<pc>'],
#          theta_cyl,label=r'$\theta_\mathrm{cylinder}$')

# plt.scatter(new_hi_rad/oldcat['distance<pc>'],
#          theta_crc,label=r'$\theta_\mathrm{circle}$')


plt.scatter(oldcat['name'],(np.abs(theta_crc-theta_cyl)/theta_crc),marker='+',color='black')

plt.scatter(oldcat['name'][scene2],(np.abs(theta_crc-theta_cyl)/theta_crc)[scene2],marker='+',color='red')
plt.xticks(rotation='vertical')
plt.hlines([0.1,0.2,0.3],'aquila_s','taurus')

plt.ylabel(r'$\Delta \theta\  /\ \theta_\mathrm{circle}$',)
plt.show()
# plt.legend()


np.array((oldcat['name'],
          100*np.abs(np.array((theta_crc/theta_cyl)).transpose()-1))).transpose()


bns = np.arange(0,100,4)


plt.hist(100*oldcat['err_M_GMC<Msun>']/oldcat['mass_GMC<Msun>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_Omg_GMC<km/s/pc>']/oldcat['vel_grad_GMC<km/s/pc>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_M_HI<Msun>']/oldcat['mass_HI<Msun>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_Omg_HI<km/s/pc>']/oldcat['vel_grad_HI<km/s/pc>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_j_GMC<cm2/s>']/oldcat['j_GMC<cm2/s>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_j_HI<cm2/s>']/oldcat['j_HI<cm2/s>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_R_GMC<pc>']/oldcat['radius_GMC<pc>'],alpha=0.2,bins=bns)
plt.hist(100*oldcat['err_R_HI<pc>']/oldcat['radius_HI<pc>'],alpha=0.2,bins=bns)


plt.hist(2*100*np.abs(theta_crc-theta_cyl)/theta_crc,bins=bns,histtype='step',color='black')

# plt.show()


plt.show()


beta_gmc = (((1/3)*(oldcat['vel_grad_GMC<km/s/pc>']*u.km/u.s/u.pc)**2*\
(oldcat['radius_GMC<pc>']*u.pc)**3)/(c.G*oldcat['mass_GMC<Msun>']*u.solMass)).cgs


beta_hi = (((1/3)*(oldcat['vel_grad_HI<km/s/pc>']*u.km/u.s/u.pc)**2*\
(oldcat['radius_HI<pc>']*u.pc)**3)/(c.G*oldcat['mass_HI<Msun>']*u.solMass)).cgs


plt.plot(beta_gmc,marker='.')
plt.plot(beta_hi,marker='.')

plt.grid()
plt.xticks(np.arange(0,22,2))

# plt.hlines(0.02,0,22)
plt.show()


alpha_gmc = (((5*oldcat['sigma_v_GMC<km/s>']*u.km/u.s)**2*
oldcat['radius_GMC<pc>']*u.pc)/(c.G*oldcat['mass_GMC<Msun>']*u.solMass)).cgs


(((5*oldcat['sigma_v_GMC<km/s>']*u.km/u.s)**2*
oldcat['radius_GMC<pc>']*u.pc)/(c.G*oldcat['mass_GMC<Msun>']*u.solMass)).cgs


(((5*oldcat['sigma_v_HI<km/s>']*u.km/u.s)**2*
oldcat['radius_HI<pc>']*u.pc)/(c.G*oldcat['mass_HI<Msun>']*u.solMass)).cgs


alpha_hi = (((5*oldcat['sigma_v_HI<km/s>']*u.km/u.s)**2*
oldcat['radius_HI<pc>']*u.pc)/(c.G*oldcat['mass_HI<Msun>']*u.solMass)).cgs


plt.plot(alpha_gmc,marker='.')
plt.plot(alpha_hi,marker='.')

plt.grid()
plt.xticks(np.arange(0,22,2))

# plt.hlines(0.02,0,22)
plt.show()


#hmmm this seems wrong, 
#when i come back, try to figure out whats going on here,
#then see if beta can tell us anything
#then read open tabs for ideas on how to calculate a timescale


def v_from_ld (l,d):
    Ro = 8.5 *u.kpc
    Vo = 220 *u.km/u.s
    
    Rg = np.sqrt(Ro**2 + d**2 -2*Ro*d*np.cos(l*u.deg))
    Reff = Ro*np.sin(l*u.deg)
    
    
    alpha = np.arcsin(Reff/Rg).to(u.deg).value
    
    beta = 90-alpha
    
    vlsr_l = np.cos(beta*u.deg)*Vo
    
    return vlsr_l-Vo*np.sin(l*u.deg)


[v_from_ld(oldcat['lon_GMC<deg>'][i],oldcat['distance<pc>'][i]*u.pc) for i in [3]]


ld_grid = np.array(np.meshgrid(np.arange(180,-180,-0.125),np.arange(0,2500,10))).transpose().reshape(-1,2)


vs_from_ld = [v_from_ld(l,d*u.pc).value for l,d in ld_grid]


ls,ds = ld_grid.transpose()


plt.figure(dpi=200)
plt.scatter(ds*np.sin(ls*u.deg),
            ds*np.cos(ls*u.deg),
            c=vs_from_ld,cmap=fns.segmented_spectral,s=1e-2)
plt.colorbar()
plt.gca().set_aspect(1)


plt.scatter(oldcat['distance<pc>']*np.sin(np.nanmean(oldcat['lon_GMC<deg>'],axis=1)*u.deg),
            oldcat['distance<pc>']*np.cos(np.nanmean(oldcat['lon_GMC<deg>'],axis=1)*u.deg),
                   marker='+',color='black')

plt.xlim(-2000,1000)
plt.ylim(-2000,500)


time_ish = ((((oldcat['j_HI<cm2/s>']-oldcat['j_GMC<cm2/s>'])*(u.cm**2/u.s))*\
(oldcat['radius_HI<pc>']*u.pc-oldcat['radius_GMC<pc>']*u.pc))/\
(c.G*(oldcat['mass_HI<Msun>']*u.solMass-oldcat['mass_GMC<Msun>']*u.solMass))).to(u.yr)


plt.plot(time_ish,marker='o')
plt.yscale('log')
plt.xticks(np.arange(0,22,2))
plt.grid()


sfrs_old = QTable(fits.open('../not_rn_7/re_expanded_cat.fits')[1].data)['SFR<Msun/yr>']


plt.scatter(time_ish,sfrs_old)
plt.loglog()


hihi = fits.open('../rn/hi_hisized.fits')


for i in range(22):
    fig,axs = plt.subplots(1,3,figsize=(15,7),dpi=200)
    axs[0].contourf(
        np.log10(np.nansum(hihi[i].data,axis=1)*0.125*1.82e18),
        levels=np.linspace(20,21.5,20),cmap='magma'),
    
    axs[1].contourf(
        np.log10(np.nansum(hihi[i].data,axis=2)*0.125*1.82e18),
        ),
    plt.show()


fns.make_fits_from_list(new_coco,'h2h2.fits')
fns.make_fits_from_list(new_hihi,'hihi.fits')


newcat


np.array((new_h2_mass/oldcat['mass_GMC<Msun>'])).transpose()


#something is going with these but i will just ignore it for now
#the 0.76923077 factor popping up over and over is weird tho
#so it must have been something systematic-ish that i did?
#maybe some difference between the display masks and the 
plt.plot(np.array((new_h2_mass/oldcat['mass_GMC<Msun>'])).transpose(),marker='.')
# plt.plot(np.array((new_h2_rad/oldcat['radius_GMC<pc>'])).transpose(),marker='.')
# plt.plot(np.array((new_hi_rad/oldcat['radius_HI<pc>'])).transpose(),marker='.')
