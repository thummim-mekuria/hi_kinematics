#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle

from astropy.io import fits
from astropy.table import QTable

import astropy.units as u
import astropy.constants as c


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


bigh1 = fits.open('rn/bigh1.fits')[0].data


bigh0 = fits.open('rn/bigh0.fits')[0].data


df = pd.read_csv('Untitled Folder/mromgs.csv')


cat = QTable(fits.open('rn/sighhplease.fits')[1].data)


hi_r_deg = np.arctan(np.array(df['hir<pc>'])/cat['distance<pc>'])*u.rad.to(u.deg)


ommaghi = np.sqrt(df['hio<m/s/pc>']/np.nanmax(df['hio<m/s/pc>']))#rankdata(df['hio<m/s/pc>'])/22.
ommagh2 = np.sqrt(df['h2o<m/s/pc>']/np.nanmax(df['h2o<m/s/pc>']))#rankdata(df['h2o<m/s/pc>'])/22.


lc = np.nanmean([cat['loc'],cat['lic']],axis=0)*-0.125+360
bc = np.nanmean([cat['boc'],cat['bic']],axis=0)*0.125-50


plt.rcParams['font.size'] = 40.


plt.figure(dpi=500,figsize=(50,25))

plt.contourf(np.where(bigh0<4e20,0.,bigh0),
             extent=[360,-360,-50,50],levels=30,cmap='binary_r',norm='log')
# plt.contourf(bigh1,cmap=redses,extent=[360,-360,-50,50],
#              levels=np.arange(-51,0,5),extend='min',alpha=0.5,antialiased=True)
# plt.contourf(bigh1,cmap=bluess,extent=[360,-360,-50,50],
#              levels=np.arange(1,60,5),extend='max',alpha=0.5,antialiased=True)


plt.xticks(np.arange(180,-182,-60),
          [str(tck)+'$\degree$' for tck in np.arange(180,-182,-60)])
plt.yticks(np.arange(-40,42,20),
          [str(tck)+'$\degree$' for tck in np.arange(-40,42,20)])

plt.xlim(185,-185)
plt.ylim(-45,45)
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.gca().set_aspect(1.)


for i,x,y,d,h2dx_pc,h2dy_pc,hidx_pc,hidy_pc in\
zip(range(22),lc,bc,cat['distance<pc>'],
                    h2_dvdx/df['h2o<m/s/pc>']*1e3,
                    h2_dvdy/df['h2o<m/s/pc>']*1e3, 
                    hi_dvdx/df['hio<m/s/pc>']*1e3,
                    hi_dvdy/df['hio<m/s/pc>']*1e3):
    
    h2dx_dg = h2dx_pc*d*u.deg.to(u.rad)
    h2dy_dg = h2dy_pc*d*u.deg.to(u.rad)
    
    
    hidx_dg = hidx_pc*d*u.deg.to(u.rad)
    hidy_dg = hidy_pc*d*u.deg.to(u.rad)
    
    plt.quiver(x,y,-hidx_dg,hidy_dg,
               pivot='mid',
               facecolor=[0.8,0.1,0.8,ommaghi[i]**2],
               edgecolor='white',linewidth=1.2,width=5e-3)
    
    
    plt.quiver(x,y,-h2dx_dg,h2dy_dg,pivot='mid',
               facecolor=[0.1,0.8,0.1,ommagh2[i]**2],#=cm.magma. ,→resampled(22)(int(22*ommagh2[i])),
               edgecolor='white',linewidth=1.2,width=2e-3
                
              )
    

plt.savefig('bigomg_galmom1',bbox_inches='tight')
plt.show()


h2_msks = [(fits.open('rn//h2_0thmom_masks.fits')[i].data-1).astype('bool')
           for i in range(22)]


h2_mom0 = fits.open('rn/h2_0s.fits')
# [np.where(~h2_msks[i],np.nan,
#                     np.nansum(h2h2[i].data[vos[i]:vis[i],:,:],axis=0)*1.3*1.97e20)
#             for i in range(22)]


plt.figure(dpi=500,figsize=(50,25))

plt.contourf(np.where(bigh0<4e20,0.,bigh0),
             extent=[360,-360,-50,50],levels=30,cmap='binary_r',zorder=-40,norm='log')


plt.xticks(np.arange(180,-182,-60),
          [str(tck)+'$\degree$' for tck in np.arange(180,-182,-60)])
plt.yticks(np.arange(-40,42,20),
          [str(tck)+'$\degree$' for tck in np.arange(-40,42,20)])

plt.xlim(185,-185)
plt.ylim(-45,45)
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.gca().set_aspect(1.)


[(
plt.contourf(np.where(~h2_msks[i],np.nan,np.log10(h2_mom0[i].data)),
                 levels=np.arange(20,22.5,0.1),
                 extent=[cat['loc'][i]*-0.125+360,
                         cat['lic'][i]*-0.125+360,
                         cat['boc'][i]*0.125-50,
                         cat['bic'][i]*0.125-50],
                extend='max',antialiased=True,
                cmap='magma',zorder=10),
    
    plt.contour(h2_msks[i],
                 extent=[cat['loc'][i]*-0.125+360,
                         cat['lic'][i]*-0.125+360,
                         cat['boc'][i]*0.125-50,
                         cat['bic'][i]*0.125-50],
                extend='max',linewidths=0.6,
                colors='white',zorder=8),
    
    plt.gca().add_patch(Circle((lc[i],bc[i]),hi_r_deg[i],
                               fill=False,color='white',zorder=-10)),
    
    plt.gca().add_patch(Circle((lc[i]+360,bc[i]),hi_r_deg[i],
                               fill=False,color='white',zorder=-10)),
    
    plt.gca().add_patch(Circle((lc[i]-360,bc[i]),hi_r_deg[i],
                               fill=False,color='white',zorder=-10))
    
    
)for i in range(22)]


plt.savefig('bigsig_galmom0',bbox_inches='tight')


plt.show()


hihi = fits.open('Untitled Folder/hihi.fits')
h2h2 = fits.open('Untitled Folder/h2h2.fits')


names = [
    'Aquila South',
    'California',
    'Camelopardalis',
    'Chamaeleon',
    'Canis Major OB1',
'Corona Australis',
'Crossbones',
'Gemini OB1',
'Hercules',
'Lacerta',
'Lupus',
'Maddalena',
'Monoceros OB1',
'Monoceros R2',
'Ophiuchus',
'Orion A',
'Orion B',
r'$\lambda$-Orionis',
'Perseus',
'Polaris',
'Rosette',
'Taurus']


lbvs = pd.read_csv('Untitled Folder/lbv_bounds_both.csv')


vo,vi = lbvs['vo'],lbvs['vi']


np.arange(-100,100,20)


fig,axs = plt.subplots(2,2,figsize=(10,10))

[(
#     plt.figure(dpi=500),
    
    axs.flat[j].plot(np.nanmean(hihi[i].data,axis=(1,2)),color=[0.8,0.1,0.8],linewidth=3.5,linestyle='--'),
    axs.flat[j].plot(np.arange(vo[i],vi[i]),
             np.nanmean(hihi[i].data,axis=(1,2))[vo[i]:vi[i]],color=[0.8,0.1,0.8],linewidth=3.5),
    
    
    axs.flat[j].plot(100*np.nanmean(h2h2[i].data,axis=(1,2)),color='black',linewidth=1.5,linestyle=':'), 
    axs.flat[j].plot(100*np.nanmean(h2h2[i].data,axis=(1,2)),
             color=[0.1,0.8,0.1],linewidth=0.8,linestyle=':',
            ), 
    
    axs.flat[j].plot(np.arange(vo[i],vi[i]),
             100*np.nanmean(h2h2[i].data,axis=(1,2))[vo[i]:vi[i]],color='black',linewidth=1.5), 
    axs.flat[j].plot(np.arange(vo[i],vi[i]),
             100*np.nanmean(h2h2[i].data,axis=(1,2))[vo[i]:vi[i]],color=[0.1,0.8,0.1],linewidth=0.8), 
    axs.flat[j].set_title(names[i]),
    
    axs.flat[j].set_xticks((np.arange(-40,40,10)+100)/1.3,labels=[-40,' ',-20,' ',0,' ',20,' ']),
    axs.flat[j].set_xlim((np.array([-35,35])+100)/1.3),
    axs.flat[j].set_ylim(-3,65),
    axs.flat[j].grid(alpha=0.2,zorder=-10),
    axs.flat[j].set_aspect(1/1.3),
    axs.flat[j].set_yticks(np.arange(0,61,10),labels=[0,' ',20,' ',40,' ',60]),
    axs.flat[j].label_outer()
    
     
)for i,j in zip([2,8,18,20],[0,1,2,3])]

fig.supylabel(r'$\overline{T_\mathrm{B}}\ [\mathrm{K}]$')
fig.supxlabel(r'$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')

fig.tight_layout()

plt.savefig('postage_spec',bbox_inches='tight')

plt.show()


hij = (np.array(df['hio<m/s/pc>'])*u.m/u.s/u.pc*\
      (np.array(df['hir<pc>'])*u.pc)**2).cgs


h2j = (np.array(df['h2o<m/s/pc>'])*u.m/u.s/u.pc*\
      (np.array(df['h2r<pc>'])*u.pc)**2).cgs


plt.hist(hij/h2j,edgecolor='white',linewidth=5.5)


np.arange(-40,40,10)


#make spectra with these


#break figure 6 into two, use colors: h2=[0.1,0.8,0.1], hi=[0.8,0.1,0.8]
#recolor figure 7, 8 (check if i can do it canva first)
#remake fig 11 to have them on top of each other

#look at word-y changes i need to make
#look at meqbetber


fake_Jh2Jmw = np.array([-0.9,-0.89,-0.85,
                        -0.7,-0.75,-0.78,-0.68,-0.67,
                       -0.13,0.13,
                       0.24,0.14,0.31,0.27,
                       0.37,0.45,
                       0.52,
                       0.68,0.71,0.77,
                       0.9,0.95])

fake_JhiJmw = np.array([-0.9,-0.89,-0.85,
                        -0.7,-0.56,
                        -0.45,-0.38,-0.43,-0.37,-0.48,
                       -0.24,-0.27,
                        -0.11,-0.05,-0.16,
                       0.13,0.14,
                       0.21,0.28,
                       0.52,0.58,
                       0.9])


plt.rcParams['font.size'] = 15


np.arange(-1,1.1,0.5/3)


plt.figure(dpi=200,figsize=(7,7))
bns = np.arange(-1,1.1,0.5/3)

plt.hist(fake_JhiJmw,bins=bns,alpha=0.7,
        color=[0.8,0.1,0.8],edgecolor=[(0.8,0.1,0.8,0.)],label='HI',zorder=12)
plt.xlabel('$\ \ \ \ \ \ \hat{J} \cdot \hat{J}_\mathrm{MW}$')

plt.hist(fake_Jh2Jmw,bins=bns,histtype='step',zorder=10,
        color=[0.1,0.8,0.1],linewidth=5,label=r'$\mathrm{H_2}$')
# plt.set_xlabel('$\hat{J_\mathrm{HI}} \cdot \hat{J_\mathrm{MW}}$')

plt.xlim(-1.05,1.05)
# plt.ylim(0,6.5)
#
# plt.vlines([-0.5,0.5,],0,6.5,color='red',linestyle='--')
plt.gca().add_patch(Rectangle((-1.1,0),0.6,11,alpha=0.1,color='black',linewidth=0.))
plt.gca().add_patch(Rectangle((-1.1,0),1.6,11,alpha=0.1,color='black',linewidth=0.))
# plt.gca().vlines([-0.71,0.71],0,7,alpha=0.7,linewidth=0.5,zorder=-20,color='black')


plt.legend()

plt.savefig('delphi_gals.png')

plt.show()


np.arange(-51,0,5)


np.arange(0.5,60,5)


plt.figure(dpi=200)
plt.contourf(bigh1,cmap='Reds_r',extent=[360,-360,-50,50],levels=np.arange(-50,0,5),extend='min')
plt.colorbar(ticks=np.arange(-50,0,10),label=r'$\overline{v_\mathrm{Gal}}\ \mathrm{[km\ s^{-1}]}$')

plt.contourf(bigh1,cmap='Blues',extent=[360,-360,-50,50],levels=np.arange(0,60,5),extend='max')


plt.xlim(180,-180)
plt.colorbar()


plt.figure(dpi=200)
plt.contourf(np.log10(bigh0),cmap='binary_r',
             extent=[360,-360,-50,50],levels=np.arange(20,22,0.15),extend='min')


plt.xlim(183,-183)
plt.gca().set_aspect(1.5)
