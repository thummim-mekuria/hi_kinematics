#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Ellipse
from matplotlib.colors import ListedColormap

import numpy as np
import pandas as pd

import datacube as dc
import fun as fns

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable

from scipy.optimize import curve_fit


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


delphi = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']
for i in range(22):
    if delphi[i] > 180:
        delphi[i] -=360
        
    if delphi[i] < -180:
        delphi[i] +=360


cat['del_phi<deg>'] = delphi


def line(x,m,b):
    return m*x+b


cat = QTable(fits.open('sighhplease.fits')[1].data)


goodman = pd.read_csv('goodman_1993_data.csv')


from astropy.modeling import models,fitting


linit = models.Linear1D(1.,0.)


fitter = fitting.LinearLSQFitter()


help(fitter)


cat


help(fitting.LinearLSQFitter)


def fit_line(x,y,yerr):
    raw = curve_fit(line,x,y,sigma=yerr)
    
    return raw[0],np.sqrt(np.diag(raw[1]))


fit_line(np.log10(cat['radius_GMC<pc>']),np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['err_j_GMC<cm2/s>']))


fit_line(np.log10(cat['radius_HI<pc>']),np.log10(cat['j_HI<cm2/s>']),np.log10(cat['err_j_HI<cm2/s>']))


fit_line(np.log10(goodman['R_eff<pc>']),np.log10(goodman['j<cm2/s>']),np.log10(goodman['err_j<cm2/s>']))


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))

axs[0].errorbar(cat['radius_GMC<pc>'],0.4*cat['j_GMC<cm2/s>'],
             xerr=cat['err_R_GMC<pc>'],yerr=0.4*cat['err_j_GMC<cm2/s>'],
             fmt='.',elinewidth=0.6,color='black',ms=3)

axs[0].loglog()

axs[0].plot(np.arange(1e0,5e2,5e0),np.arange(1e0,5e2,5e0)**1.69 * 10**(22.37),
         linewidth=0.75,color=fns.segmented_spectral(19),label = r'$j \propto R^{1.69 \pm 0.21}$')

axs[0].plot(np.arange(1e0,5e2,5e1),np.arange(1e0,5e2,5e1)**1.30 * 10**(22.91),
         linewidth=0.5,alpha=0.3,linestyle='--',color='black')


axs[0].fill_between(np.arange(1e0,5e2,5e0),np.arange(1e0,5e2,5e0)**(1.69+0.21) * 10**(22.37),
                np.arange(1e0,5e2,5e0)**(1.69-0.21) * 10**(22.37),alpha=0.1,
                 color=fns.segmented_spectral(19),linewidth=0.)


axs[0].set_xlabel(r'$R_\mathrm{H_2}\ [\mathrm{pc}]$')
axs[0].set_ylabel(r'$j\ [\mathrm{cm^2\ s^{-1}}]$')

axs[0].set_xlim(4e0,130)
# axs[0].set_ylim(3e23,2e26)
axs[0].set_ylim(1e23,6e25)


# axs[0].set_aspect(0.5)

axs[0].legend(fontsize='small')

# axs[0].set_title(r'$\mathrm{H_2}$')


axs[1].errorbar(cat['radius_HI<pc>'],0.4*cat['j_HI<cm2/s>'],
             xerr=cat['err_R_HI<pc>'],yerr=0.4*cat['err_j_HI<cm2/s>'],
             fmt='.',elinewidth=0.6,color='black',ms=3)

axs[1].loglog()

axs[1].plot(np.arange(1e1,5e2,5e1),np.arange(1e1,5e2,5e1)**1.30 * 10**(22.91),
         linewidth=0.75,color='black',label = r'$j \propto R^{1.30 \pm 0.38}$')

axs[1].plot(np.arange(1e0,5e2,5e0),np.arange(1e0,5e2,5e0)**1.69 * 10**(22.37),
         linewidth=0.5,alpha=0.5,linestyle='--',color=fns.segmented_spectral(19))


axs[1].fill_between(np.arange(1e1,5e2,5e1),np.arange(1e1,5e2,5e1)**(1.30+0.38) * 10**(22.91),
                np.arange(1e1,5e2,5e1)**(1.30-0.38) * 10**(22.91),alpha=0.1,
                 color='black',linewidth=0.)


axs[1].set_xlabel(r'$R_\mathrm{HI}\ [\mathrm{pc}]$')
# axs[0].set_ylabel(r'$j\ [\mathrm{cm^2\ s^{-1}}]$')

axs[1].set_xlim(1e1,3e2)
# axs[1].set_ylim(3e24,6e26)
axs[1].set_ylim(1e24,2e26)


# axs[1].set_aspect(0.6)
axs[1].legend(fontsize='small')

# axs[1].set_title(r'$\mathrm{HI}$')


fig.tight_layout()

plt.savefig('jR_H2hi_sep',bbox_inches='tight')


all_rs = np.concatenate((goodman['R_eff<pc>'],cat['radius_GMC<pc>'],cat['radius_HI<pc>']))


all_js = np.concatenate((0.4*goodman['j<cm2/s>'],0.4*cat['j_GMC<cm2/s>'],0.4*cat['j_HI<cm2/s>']))


all_ms = np.concatenate((goodman['M<Msun>'],cat['mass_GMC<Msun>'],cat['mass_HI<Msun>']))


all_err_js = np.concatenate((0.4*goodman['err_j<cm2/s>'],0.4*cat['err_j_GMC<cm2/s>'],0.4*cat['err_j_HI<cm2/s>']))


fit_line(np.log10(all_rs),np.log10(all_js),np.log10(all_err_js))


plt.figure(dpi=200)
plt.errorbar(goodman['R_eff<pc>'],0.4*goodman['j<cm2/s>'],
            yerr=0.4*goodman['err_j<cm2/s>'],fmt='.',elinewidth=0.6,
             color=fns.segmented_spectral(4),label='G+93')


plt.errorbar(cat['radius_GMC<pc>'],0.4*cat['j_GMC<cm2/s>'],
             xerr=cat['err_R_GMC<pc>'],yerr=0.4*cat['err_j_GMC<cm2/s>'],
             fmt='.',elinewidth=0.6,color=fns.segmented_spectral(19),label='$\mathrm{H_2}$')

plt.errorbar(cat['radius_HI<pc>'],0.4*cat['j_HI<cm2/s>'],
             xerr=cat['err_R_HI<pc>'],yerr=0.4*cat['err_j_HI<cm2/s>'],
             fmt='.',elinewidth=0.6,color='black',label='$\mathrm{HI}$')


plt.plot(np.arange(1e-2,1e3,1e-1),np.arange(1e-2,1e3,1e-1)**1.46 * 10**(22.64),linewidth=0.75,color='grey',label = r'$j\ \propto R^{1.46\pm0.04}$',linestyle='--')

plt.fill_between(np.arange(1e-2,1e3,1e-1),
                 np.arange(1e-2,1e3,1e-1)**(1.46-0.04) * 10**(22.64),
                 np.arange(1e-2,1e3,1e-1)**(1.46+0.04) * 10**(22.64),alpha=0.1,color='grey')

plt.loglog()
plt.legend()

plt.xlabel(r'$R\ [\mathrm{pc}]$')
plt.ylabel(r'$j\ [\mathrm{cm^2\ s^{-1}}]$')

plt.xlim(5e-2,3e2)
plt.ylim(5e20,5e26)


plt.savefig('jR_all',bbox_inches='tight')


plt.figure(dpi=200)
plt.plot(np.linspace(7e-2,1e0,10),np.full(10,1.63),color=fns.segmented_spectral(4),label='cores')
plt.fill_between(np.linspace(7e-2,1e0,10),np.full(10,1.63+0.18),np.full(10,1.63-0.18),alpha=0.2,color=fns.segmented_spectral(4),linewidth=0.)
         
plt.plot(np.linspace(4e0,1e2,10),np.full(10,1.69),color=fns.segmented_spectral(19),label='molecular clouds')
plt.fill_between(np.linspace(4e0,1e2,10),np.full(10,1.69+0.21),np.full(10,1.69-0.21),alpha=0.1,color=fns.segmented_spectral(19),linewidth=0.)

         
plt.plot(np.linspace(1e1,3e2,10),np.full(10,1.30),color='black',label='atomic envelopes')
plt.fill_between(np.linspace(1e1,3e2,10),np.full(10,1.30-0.38),np.full(10,1.30+0.38),alpha=0.1,color='black',linewidth=0.)


plt.plot(np.linspace(7e-2,3e2,10),np.full(10,1.46),color='grey',linestyle='--',label='all')
plt.fill_between(np.linspace(7e-2,3e2,10),np.full(10,1.46-0.04),np.full(10,1.46+0.04),alpha=0.1,color='grey',linewidth=0.)


plt.xscale('log')

plt.ylabel('index of $\ j-R\ $ relationship')
plt.xlabel('$R$ [pc]')

plt.legend(fontsize='x-small',loc='lower left')

plt.savefig('jR_index_coresH2hi',bbox_inches='tight')


np.nanmin(cat['radius_GMC<pc>']),np.nanmax(cat['radius_GMC<pc>'])


np.nanmin(cat['radius_HI<pc>']),np.nanmax(cat['radius_HI<pc>'])


fit_line(np.log10(cat['mass_GMC<Msun>']),np.log10(0.4*cat['j_GMC<cm2/s>']),np.log10(0.4*cat['err_j_GMC<cm2/s>']))


fit_line(np.log10(cat['mass_HI<Msun>']),np.log10(0.4*cat['j_HI<cm2/s>']),np.log10(0.4*cat['err_j_HI<cm2/s>']))


fit_line(np.log10(goodman['M<Msun>']),np.log10(goodman['j<cm2/s>']),np.log10(goodman['err_j<cm2/s>']))


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))

axs[0].errorbar(cat['mass_GMC<Msun>'],0.4*cat['j_GMC<cm2/s>'],
             xerr=cat['err_M_GMC<Msun>'],yerr=0.4*cat['err_j_GMC<cm2/s>'],
             fmt='.',elinewidth=0.6,color='black',ms=3)

axs[0].loglog()

axs[0].plot(np.arange(3e2,5e6,5e3),np.arange(3e2,5e6,5e3)**0.65 * 10**(21.59),
         linewidth=0.75,color=fns.segmented_spectral(19),label = r'$j \propto M^{0.65 \pm 0.10}$')


axs[0].plot(np.arange(5e2,2e7,5e4),np.arange(5e2,2e7,5e4)**0.56 * 10**(22.23),
         linewidth=0.5,alpha=0.3,linestyle='--',color='black')


axs[0].fill_between(np.arange(3e2,5e6,5e3),np.arange(3e2,5e6,5e3)**(0.65+0.10) * 10**(21.59),
                np.arange(3e2,5e6,5e3)**(0.65-0.10) * 10**(21.59),alpha=0.1,
                 color=fns.segmented_spectral(19),linewidth=0.)


axs[0].set_xlabel(r'$M_\mathrm{H_2}\ [\mathrm{M_\odot}]$')
axs[0].set_ylabel(r'$j\ [\mathrm{cm^2\ s^{-1}}]$')

axs[0].set_xlim(5e2,2e6)
axs[0].set_ylim(1e23,6e25)

# axs[0].set_aspect(1)

axs[0].legend(fontsize='small')

# axs[0].set_title(r'$\mathrm{H_2}$')


axs[1].errorbar(cat['mass_HI<Msun>'],0.4*cat['j_HI<cm2/s>'],
             xerr=cat['err_M_HI<Msun>'],yerr=0.4*cat['err_j_HI<cm2/s>'],
             fmt='.',elinewidth=0.6,color='black',ms=3)

axs[1].loglog()

axs[1].plot(np.arange(8e3,2e7,5e4),np.arange(8e3,2e7,5e4)**0.56 * 10**(22.23),
         linewidth=0.75,color='black',label = r'$j \propto M^{0.56 \pm 0.13}$')


axs[1].plot(np.arange(3e2,1e7,5e3),np.arange(3e2,1e7,5e3)**0.65 * 10**(21.59),
          linewidth=0.5,alpha=0.5,linestyle='--',color=fns.segmented_spectral(19))


axs[1].fill_between(np.arange(8e3,2e7,5e4),np.arange(8e3,2e7,5e4)**(0.56+0.13) * 10**(22.23),
                np.arange(8e3,2e7,5e4)**(0.56-0.13) * 10**(22.23),alpha=0.1,
                 color='black',linewidth=0.)


axs[1].set_xlabel(r'$M_\mathrm{HI}\ [\mathrm{M_\odot}]$')
# axs[0].set_ylabel(r'$j\ [\mathrm{cm^2\ s^{-1}}]$')

axs[1].set_xlim(1e4,1e7)
axs[1].set_ylim(1e24,2e26)

# axs[1].set_aspect(1)
axs[1].legend(fontsize='small',loc='upper left')

# axs[1].set_title(r'$\mathrm{HI}$')


fig.tight_layout()

plt.savefig('jM_H2hi_sep',bbox_inches='tight')


fit_line(np.log10(all_ms),np.log10(all_js),np.log10(all_err_js))


plt.figure(dpi=200)
plt.errorbar(goodman['M<Msun>'],0.4*goodman['j<cm2/s>'],
            yerr=0.4*goodman['err_j<cm2/s>'],fmt='.',elinewidth=0.6,
             color=fns.segmented_spectral(4),label='G+93')


plt.errorbar(cat['mass_GMC<Msun>'],0.4*cat['j_GMC<cm2/s>'],
             xerr=cat['err_M_GMC<Msun>'],yerr=0.4*cat['err_j_GMC<cm2/s>'],
             fmt='.',elinewidth=0.6,color=fns.segmented_spectral(19),label='$\mathrm{H_2}$')

plt.errorbar(cat['mass_HI<Msun>'],0.4*cat['j_HI<cm2/s>'],
             xerr=cat['err_M_HI<Msun>'],yerr=0.4*cat['err_j_HI<cm2/s>'],
             fmt='.',elinewidth=0.6,color='black',label='$\mathrm{HI}$')


plt.plot(np.arange(2e-1,2e7,1e1),np.arange(2e-1,2e7,1e1)**0.77 * 10**(20.98),linewidth=0.75,color='grey',label = r'$j\ \propto M^{0.77\pm0.03}$',linestyle='--')

plt.fill_between(np.arange(2e-1,2e7,1e1),
                 np.arange(2e-1,2e7,1e1)**(0.77-0.03) * 10**(20.98),
                 np.arange(2e-1,2e7,1e1)**(0.77+0.03) * 10**(20.98),alpha=0.1,color='grey')

plt.loglog()
plt.legend()

plt.xlabel(r'$M\ [\mathrm{M_\odot}]$')
plt.ylabel(r'$j\ [\mathrm{cm^2\ s^{-1}}]$')

plt.xlim(3e-1,2e7)
plt.ylim(5e20,5e26)

# plt.gca().set_aspect(0.6)
# plt.xticks([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7],minor=True)

plt.savefig('jM_all',bbox_inches='tight')


np.nanmin(goodman['M<Msun>']),np.nanmax(goodman['M<Msun>'])


np.nanmin(cat['mass_GMC<Msun>']),np.nanmax(cat['mass_GMC<Msun>'])


np.nanmin(cat['mass_HI<Msun>']),np.nanmax(cat['mass_HI<Msun>'])


plt.figure(dpi=200)
plt.plot(np.linspace(6e-1,8e2,10),np.full(10,0.40),color=fns.segmented_spectral(4),label='cores')
plt.fill_between(np.linspace(6e-1,8e2,10),np.full(10,0.40+0.12),np.full(10,0.40-0.12),alpha=0.2,color=fns.segmented_spectral(4),linewidth=0.,)
         
plt.plot(np.linspace(8e2,1e6,10),np.full(10,0.65),color=fns.segmented_spectral(19),label='molecular clouds')
plt.fill_between(np.linspace(8e2,1e6,10),np.full(10,0.65+0.10),np.full(10,0.65-0.10),alpha=0.1,color=fns.segmented_spectral(19),linewidth=0.)

         
plt.plot(np.linspace(1e4,7e6,10),np.full(10,0.56),color='black',label='atomic envelopes')
plt.fill_between(np.linspace(1e4,7e6,10),np.full(10,0.56-0.13),np.full(10,0.56+0.13),alpha=0.1,color='black',linewidth=0.,)


plt.plot(np.linspace(2e-1,2e7,10),np.full(10,0.77),linewidth=0.75,color='grey',linestyle='--',label='all')
plt.fill_between(np.linspace(2e-1,2e7,10),np.full(10,0.77-0.03),np.full(10,0.77+0.03),alpha=0.1,color='grey',linewidth=0.)


plt.xscale('log')

plt.ylabel('index of $\ j-M\ $ relationship')
plt.xlabel('$M\ [\mathrm{M_\odot}]$')

plt.xlim(3e-1,2e7)
plt.legend(fontsize='x-small',loc='lower right')


plt.savefig('jM_index_coresH2hi',bbox_inches='tight')


help(plt.xticks)


plt.figure(dpi=200)

plt.gca().add_patch(Rectangle((-60,0),120,5,alpha=0.2,color=fns.segmented_spectral(17),zorder=-20,linewidth=0.))
plt.gca().add_patch(Rectangle((-190,0),70,5,alpha=0.2,color=fns.segmented_spectral(6),zorder=-20,linewidth=0.))
plt.gca().add_patch(Rectangle((120,0),80,5,alpha=0.2,color=fns.segmented_spectral(6),zorder=-20,linewidth=0.))
plt.gca().add_patch(Rectangle((-120,0),60,5,alpha=0.2,color=fns.segmented_spectral(12),zorder=-20,linewidth=0.))
plt.gca().add_patch(Rectangle((60,0),60,5,alpha=0.2,color=fns.segmented_spectral(12),zorder=-20,linewidth=0.))


plt.hist(cat['del_phi<deg>'],bins=np.arange(-180,210,20),
         color=fns.segmented_spectral(19),edgecolor='black',alpha=0.8)

plt.hist(cat['del_phi<deg>'],bins=np.arange(-180,210,20),
         histtype='step',edgecolor='black',alpha=0.5)

# plt.vlines([-60,60,-120,120],0,5,color=fns.segmented_spectral(4),linestyle='--',linewidth=0.75)

plt.xticks(np.arange(-180,210,60),[str(tck)+'$\degree$' for tck in np.arange(-180,210,60)])
plt.xlim(-185,185)
plt.ylim(0,4.5)
plt.ylabel('N (clouds per bin)')
plt.xlabel('$\phi_\mathrm{H_2} - \phi_\mathrm{HI}$')

plt.savefig('delphi_hist',bbox_inches='tight')


plt.scatter(cat['name'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
plt.scatter(cat['name'],cat['jsum_HI<cm2/s>']/cat['jsum_GMC<cm2/s>'])

plt.xticks(rotation='vertical')
plt.yscale('log')
plt.hlines(1,'aquila_s','taurus')
plt.show()


plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
# plt.xticks(rotation='vertical')
# plt.scatter(cat['mass_HI<Msun>'],cat['jsum_HI<cm2/s>']/cat['jsum_GMC<cm2/s>'])

plt.loglog()
plt.show()


plt.scatter(cat['sigma_v_GMC<km/s>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
# plt.xticks(rotation='vertical')
plt.loglog()
plt.show()


curve_fit(line,np.log10(cat['sigma_v_HI<km/s>']/cat['sigma_v_GMC<km/s>']),np.log10(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>']))


plt.scatter(cat['sigma_v_HI<km/s>']/cat['sigma_v_GMC<km/s>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
# plt.xticks(rotation='vertical')
plt.loglog()
plt.show()


np.array((cat['name'],cat['jsum_HI<cm2/s>']/cat['j_HI<cm2/s>'])).transpose()


plt.hist(cat['jsum_HI<cm2/s>']/cat['j_HI<cm2/s>'],bins=10**np.arange(-1.,0.5,0.1),histtype='step',color='black',linewidth=1.3)
plt.hist(cat['jsum_GMC<cm2/s>']/cat['j_GMC<cm2/s>'],bins=10**np.arange(-1.,0.5,0.1),histtype='step',color=fns.segmented_spectral(19))

plt.xscale('log')

plt.show()


10**(-1.2),10**(-1.3),10**(-1.4)


plt.figure(dpi=200)
plt.hist(cat['J_HI<erg.s>']/cat['J_GMC<erg.s>'],bins=10**np.arange(-3.,3.,1/4),
        color=fns.segmented_spectral(19),edgecolor='black',alpha=0.8)

plt.vlines(1,0,6.5,linestyle='--',color=fns.segmented_spectral(4))

plt.xscale('log')
plt.xlim(1e-2,1e3)


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
axs[0].hist(cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'],bins=10**np.arange(-3.,3.,1/4),
        color=fns.segmented_spectral(19),edgecolor='black',alpha=0.8)

axs[0].vlines(1,0,6.5,linestyle='--',color=fns.segmented_spectral(4))

axs[0].set_xscale('log')
axs[0].set_xlim(1e-2,1e3)

axs[0].set_xlabel('$j_\mathrm{HI}\ /\ j_\mathrm{H_2}$')
axs[0].set_ylabel('N (clouds per bin)')
axs[0].set_ylim(0,6.5)


axs[1].hist(cat['J_HI<erg.s>']/cat['J_GMC<erg.s>'],bins=10**np.arange(-3.,3.,1/4),
        color=fns.segmented_spectral(19),edgecolor='black',alpha=0.8)

axs[1].vlines(1,0,6.5,linestyle='--',color=fns.segmented_spectral(4))

axs[1].set_xscale('log')
axs[1].set_xlim(1e-2,1e3)
axs[1].set_xlabel('$J_\mathrm{HI}\ /\ J_\mathrm{H_2}$')
axs[1].set_ylim(0,4.5)
axs[1].set_yticks([0,1,2,3,4])

plt.savefig('j_J_ratios',bbox_inches='tight')


plt.figure(dpi=200)
plt.hist(cat['mass_HI<Msun>']/cat['mass_GMC<Msun>'],bins=np.arange(0,20,1),edgecolor='black',color=fns.segmented_spectral(19),alpha=0.8)
# plt.xscale('log')

plt.xlabel('$M_\mathrm{HI}\ /\ M_\mathrm{H_2}$')
plt.ylabel('N (clouds per bin)')
plt.xlim(2,17)
plt.xticks(np.arange(3,16,3))
plt.show()


plt.scatter(cat['vel_grad_GMC<km/s/pc>']/cat['vel_grad_HI<km/s/pc>'], (cat['radius_GMC<pc>']/cat['radius_HI<pc>'])**2)
plt.loglog()


curve_fit(line,2*np.nanstd(cat['vel<km/s>'],axis=1),cat['sigma_v_HI<km/s>'])


curve_fit(line,2*np.nanstd(cat['vel<km/s>'],axis=1),cat['vel_grad_HI<km/s/pc>'])


0.2/np.sqrt(5.67e-5)


0.00237/np.sqrt(2.1216e-7)


curve_fit(line,2*np.nanstd(cat['vel<km/s>'],axis=1),cat['sigma_v_GMC<km/s>'])


0.014/np.sqrt(3.37e-4)


curve_fit(line,2*np.nanstd(cat['vel<km/s>'],axis=1),cat['vel_grad_GMC<km/s/pc>'])


0.001299/np.sqrt(8.36e-7)


sig_mom2s_hi = [np.nanmean(fits.open('hi_2s.fits')[i].data) for i in range(22)]
sig_mom2s_h2 = [np.nanmean(fits.open('h2_2s.fits')[i].data) for i in range(22)]


plt.hist(cat['sigma_v_GMC<km/s>']/sig_mom2s_h2,bins = np.arange(0,3,0.1),edgecolor='white')

plt.show()


plt.hist(cat['sigma_v_HI<km/s>']/sig_mom2s_hi,bins = np.arange(0,3,0.1),edgecolor='white')

plt.show()


fig,axs=plt.subplots(2,2,dpi=200,figsize=(10,10),sharex=True)
axs[0,0].scatter(2*np.nanstd(cat['vel<km/s>'],axis=1),cat['sigma_v_GMC<km/s>'],color=fns.segmented_spectral(19))
axs[0,1].scatter(2*np.nanstd(cat['vel<km/s>'],axis=1),cat['sigma_v_HI<km/s>'],color='black')


axs[1,0].scatter(2*np.nanstd(cat['vel<km/s>'],axis=1),cat['vel_grad_GMC<km/s/pc>'],color=fns.segmented_spectral(19))
axs[1,1].scatter(2*np.nanstd(cat['vel<km/s>'],axis=1),cat['vel_grad_HI<km/s/pc>'],color='black')


axs[0,0].set_ylabel('$\sigma_\mathrm{v}\ [\mathrm{km\ s^{-1}}]$')
# axs[0,1].set_ylabel('$\sigma_\mathrm{v,HI}\ [\mathrm{km\ s^{-1}}]$')


axs[1,0].set_ylabel('$\Omega_\mathrm{v}\ [\mathrm{km\ s^{-1}\ pc^{-1}}]$')
# axs[1,1].set_ylabel('$\Omega_\mathrm{v,HI}\ [\mathrm{km\ s^{-1}\ pc^{-1}}]$')

axs[1,0].set_yticks([0,0.1,0.2,0.3])
axs[1,1].set_yticks([0,0.07,0.14,0.21])
axs[0,1].set_yticks([6,11,16,21])


fig.supxlabel('$\Delta v\ \mathrm{[km\ s^{-1}]}$')


fig.tight_layout()

# for ax in axs.flat:
#     ax.loglog()
plt.savefig('delta_v_effects',bbox_inches='tight')


plt.scatter(cat['sigma_v_HI<km/s>'],cat['j_HI<cm2/s>']/cat['j_GMC<cm2/s>'])
plt.loglog()


curve_fit(line,np.log10(cat['j_GMC<cm2/s>']),np.log10(cat['j_HI<cm2/s>']))


0.31/np.sqrt(0.0368)


plt.scatter(cat['j_GMC<cm2/s>'],cat['j_HI<cm2/s>'])
plt.loglog()


plt.scatter(cat['J_GMC<erg.s>'],cat['J_HI<erg.s>'])
plt.loglog()


plt.hist(cat['J_GMC<erg.s>']+cat['J_HI<erg.s>'],bins=10**np.arange(61,69,0.5))
plt.xscale('log')


plt.scatter(np.nanstd(cat['vel<km/s>'],axis=1)*2, cat['vel_grad_GMC<km/s/pc>'])
# plt.scatter(np.nanstd(cat['vel<km/s>'],axis=1)*2, cat['vel_grad_HI<km/s/pc>'])

plt.loglog()


hi_1s = fits.open('hi_1s.fits')
h2_1s = fits.open('h2_1s.fits')


msks = [(fits.open('h2_0thmom_masks_DISPLAY.fits')[i].data-1).astype('bool') for i in range(22)]


plt.contourf(np.where(~msks[1],np.nan,h2_1s[1].data),levels=30,cmap=fns.segmented_spectral)


#13 : (msky<5)*(mskx>95)*(mskx<105)
#15 : (msky>100)*(msky<110)*(mskx>118)*(mskx<128)
#16 : (msky>35)*(msky<45)*(mskx>63)*(mskx<70)


plt.scatter(cat['name'],cat['vel_grad_HI<km/s/pc>'])
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.hlines([1e-2,2e-2,3e-2],'aquila_s','taurus')
plt.show()


#making masks for mr2, ora, orb bc there is a spot with weirdddd  first moment (and tbh also zeroth moment) values
for i in [13,15,16]:
    msky,mskx = np.indices(hi_1s[i].data.shape)
    
    plt.imshow(hi_1s[i].data<0,cmap='Paired',origin='lower')
    
    plt.imshow(((msky>37)*(msky<43)*(mskx>63)*(mskx<70) +\
               (msky>103)*(msky<109)*(mskx>118)*(mskx<128) +\
               (msky<3)*(mskx>97)*(mskx<103)),origin='lower',alpha=0.3)
    plt.colorbar()
    plt.title(i)
    plt.show()


hi_msks = [np.ones(hi_1s[i].data.shape).astype('bool') for i in range(22)]


msky,mskx = np.indices(hi_1s[13].data.shape)
hi_msks[13] = ~((msky<5)*(mskx>95)*(mskx<105))

msky,mskx = np.indices(hi_1s[15].data.shape)
hi_msks[15] = ~((msky>100)*(msky<110)*(mskx>118)*(mskx<128))

msky,mskx = np.indices(hi_1s[16].data.shape)
hi_msks[16] = ~((msky>35)*(msky<45)*(mskx>63)*(mskx<70))


fns.make_fits_from_list([hi_msks[i]+1 for i in range(22)],'hi_msks.fits')


h2_ss = np.array((0.6,0.2,1,0.5,0.8,0.8,0.8,0.6,0.6,0.8,0.25,0.7,0.6,0.5,0.05,0.25,0.3,0.25,0.25,0.25,0.7,0.25))


for i in range(22):
    ycoords2, xcoords2 = np.indices(h2_1s[i].data.shape)
    ycoordsi, xcoordsi = np.indices(hi_1s[i].data.shape)
    
    if i in [3,10,17,21]:
        vo2 = np.nanmean(h2_1s[i].data)
    else:
        vo2 = h2_1s[i].data[int(0.5*h2_1s[i].data.shape[0]),int(0.5*h2_1s[i].data.shape[1])]
        
    voi = np.nanmean(hi_1s[i].data)#hi_1s[i].data[int(0.5*hi_1s[i].data.shape[0]),int(0.5*hi_1s[i].data.shape[1])]
    
    
    ycoords2 = (ycoords2-np.nanmean(ycoords2))*0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    xcoords2 = (xcoords2-np.nanmean(xcoords2))*-0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    ycoordsi = (ycoordsi-np.nanmean(ycoordsi))*0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    xcoordsi = (xcoordsi-np.nanmean(xcoordsi))*-0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    vmod2 = np.where(~msks[i],np.nan,vo2 + xcoords2*cat['dvdx_GMC'][i] + ycoords2*cat['dvdy_GMC'][i])
    vraw2 = np.where(~msks[i],np.nan,h2_1s[i].data)
    
    vmodi = np.where(~hi_msks[i],np.nan,voi + xcoordsi*cat['dvdx_HI'][i] + ycoordsi*cat['dvdy_HI'][i])
    vrawi = np.where(~hi_msks[i],np.nan,hi_1s[i].data)
    
    
    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))

    rperp2 = (xcoords2*cat['dvdx_GMC'][i] + ycoords2*cat['dvdy_GMC'][i])/np.linalg.norm([cat['dvdx_GMC'][i],cat['dvdy_GMC'][i]],axis=0)
    rperpi = (xcoordsi*cat['dvdx_HI'][i] + ycoordsi*cat['dvdy_HI'][i])/np.linalg.norm([cat['dvdx_HI'][i],cat['dvdy_HI'][i]],axis=0)


    # plt.contourf(np.where(~msks[i],np.nan,vmod),cmap=fns.segmented_spectral,levels=30)

    axs[0].scatter(np.ravel(rperp2),np.ravel(vraw2),s=h2_ss[i],color=fns.segmented_spectral(19))
    axs[0].plot(np.ravel(rperp2),np.ravel(vmod2),color=fns.segmented_spectral(4))
#     axs[0].fill_between(np.ravel(rperp2),
#                         np.ravel(np.where(
#                             ~msks[i],np.nan,
#                             vo2 + xcoords2*(cat['dvdx_GMC'][i] + cat['sigma_dvdx_GMC'][i]) +\
#                                   ycoords2*(cat['dvdy_GMC'][i] + cat['sigma_dvdy_GMC'][i]))),
                                 
#                         np.ravel(np.where(
#                             ~msks[i],np.nan,
#                             vo2 + xcoords2*(cat['dvdx_GMC'][i] + cat['sigma_dvdx_GMC'][i]) +\
#                                   ycoords2*(cat['dvdy_GMC'][i] + cat['sigma_dvdy_GMC'][i]))),
#                         color=fns.segmented_spectral(27),alpha=0.1,zorder=-20)

    
    axs[1].scatter(np.ravel(rperpi),np.ravel(vrawi),s=np.cbrt(hi_1s[i].data.size)/1e3,color='black')
    axs[1].plot(np.ravel(rperpi),np.ravel(vmodi),color=fns.segmented_spectral(4))

    # plt.gca().set_aspect('equal')
    
    # fig.suptitle(i)
    fig.supxlabel('displacement from rotation axis [pc]')
    axs[0].set_ylabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
    fig.tight_layout()
    axs[0].set_title(i)
    
    plt.show()


hihi = fits.open('hi_hisized.fits')


coco = fits.open('h2_h2sized.fits')


(cat['boh'][14]+30)*0.125- 50


oph_hi.shape


342/10


tau_hi = hihi[21].data


np.nanstd(cat['vel<km/s>'][[1,21]],axis=1)


plt.figure(dpi=200,figsize=(7,10))
[(
    plt.plot(np.arange(155)*1.3-100,i*80 + np.nanmean(tau_hi[:,i*32:(i+1)*32,i*32:(i+1)*32],axis=(1,2)),
             color=fns.segmented_spectral(1-i/9),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][21],cat['vi'][21])*1.3-100,
             i*80 + np.nanmean(tau_hi[cat['vo'][21]:cat['vi'][21],i*32:(i+1)*32,i*32:(i+1)*32],axis=(1,2)),
             color=fns.segmented_spectral(1 - i/9),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(tau_hi[cat['vo'][21]:cat['vi'][21],i*32:(i+1)*32,i*32:(i+1)*32],axis=(1,2)),
                               np.arange(cat['vo'][21],cat['vi'][21])*1.3-100),i*80,
                color=fns.segmented_spectral(1 - i/9),s=10,linewidth=0.,zorder=30),
    
    plt.hlines(i*80,-100,105,linewidth=0.2,color='black',linestyle='--'),
    
    
    plt.plot(np.arange(155)*1.3-100,i*80 + 30*np.nanmean(coco[21].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
  ) for i in range(10)]

plt.xlim(-70,40)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[21].data),0,800,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])


plt.title('Taurus',fontsize='x-small',y = 0.96,x=0.1)
plt.savefig('tau_stacked_patch_spectra',bbox_inches='tight')


her_hi = hihi[8].data


her_hi.shape


134/2


#if i come back to redo these, would be great to add l,b coordinates to orient relative to the galactic plane
plt.figure(dpi=200,figsize=(7,10))
[(
    plt.plot(np.arange(155)*1.3-100,i*80 + np.nanmean(her_hi[:,i*13:(i+1)*13,60:75],axis=(1,2)),
             color=fns.segmented_spectral(1-i/9),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][8],cat['vi'][8])*1.3-100,
             i*80 + np.nanmean(her_hi[cat['vo'][8]:cat['vi'][8],i*13:(i+1)*13,60:75],axis=(1,2)),
             color=fns.segmented_spectral(1-i/9),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(her_hi[cat['vo'][8]:cat['vi'][8],i*13:(i+1)*13,60:75],axis=(1,2)),
                               np.arange(cat['vo'][8],cat['vi'][8])*1.3-100),i*80,
                color=fns.segmented_spectral(1-i/9),s=10,linewidth=0.,zorder=30),
    
    plt.hlines(i*80,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,i*80 + 40*np.nanmean(coco[8].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
  ) for i in range(10)]

plt.xlim(-20,45)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[8].data),0,760,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])


plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
plt.savefig('her_stacked_patch_spectra',bbox_inches='tight')


mad_hi.shape


plt.figure(dpi=200,figsize=(7,10))
[(
    plt.plot(np.arange(155)*1.3-100,i*80 + np.nanmean(mad_hi[:,2*i*6:(2*i+1)*6,18 + i*6:18+ (i+1)*6],axis=(1,2)),
             color=fns.segmented_spectral(1-i/300),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][11],cat['vi'][11])*1.3-100,
             i*80 + np.nanmean(mad_hi[cat['vo'][11]:cat['vi'][11],2*i*6:(2*i+1)*6,18 + i*6:18+ (i+1)*6],axis=(1,2)),
             color=fns.segmented_spectral(1 - i/9),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(mad_hi[cat['vo'][11]:cat['vi'][11],2*i*6:(2*i+1)*6,18 + i*6:18+ (i+1)*6],axis=(1,2)),
                               np.arange(cat['vo'][11],cat['vi'][11])*1.3-100),i*80,
                color=fns.segmented_spectral(1 - i/9),s=10,linewidth=0.,zorder=30),
    
    plt.plot(np.arange(155)*1.3-100,i*80 + 130*np.nanmean(coco[11].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5),
    
    
    plt.plot(np.arange(155)*1.3-100,i*80 + np.nanmean(mad_hi[:,6 + 2*i*6:6 + (2*i+1)*6,18 + i*6:18+ (i+1)*6],axis=(1,2)),
             color=fns.segmented_spectral(1-i/300),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][11],cat['vi'][11])*1.3-100,
             i*80 + np.nanmean(mad_hi[cat['vo'][11]:cat['vi'][11],6 + 2*i*6:6 + (2*i+1)*6,18 + i*6:18+ (i+1)*6],axis=(1,2)),
             color=fns.segmented_spectral(1 - i/9),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(mad_hi[cat['vo'][11]:cat['vi'][11],6 + 2*i*6:6 + (2*i+1)*6,18 + i*6:18+ (i+1)*6],axis=(1,2)),
                               np.arange(cat['vo'][11],cat['vi'][11])*1.3-100),i*80,
                color=fns.segmented_spectral(1 - i/9),s=10,linewidth=0.,zorder=30),
    
    
    plt.hlines(i*80,-100,105,linewidth=0.2,color='black'),
    
    plt.plot(np.arange(155)*1.3-100,i*80 + 130*np.nanmean(coco[11].data,axis=(1,2)),color='black',linewidth=0.6,alpha=0.5)
    
  ) for i in range(5)]

plt.xlim(-20,90)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[11].data),0,800,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])


plt.title('Maddalena',fontsize='x-small',y = 0.96,x=0.1)
# plt.savefig('tau_stacked_patch_spectra',bbox_inches='tight')


cat['rot_ax_HI<deg>'][8]


plt.plot(cat['loh']-cat['lih'],marker='+')


cat['radius_HI<pc>']


fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,10))
[(
    axs[0].plot(np.arange(155)*1.3-100,i*40 + np.nanmean(tau_hi[:,i*32:(i+1)*32,i*32:(i+1)*32],axis=(1,2)),
             color=fns.segmented_spectral(1-i/300),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    axs[0].plot(np.arange(cat['vo'][21],cat['vi'][21])*1.3-100,
             i*40 + np.nanmean(tau_hi[cat['vo'][21]:cat['vi'][21],i*32:(i+1)*32,i*32:(i+1)*32],axis=(1,2)),
             color=fns.segmented_spectral(1 - i/9),
            linewidth=1.5,),
    
    axs[0].scatter(fns.get_mom1(np.nanmean(tau_hi[cat['vo'][21]:cat['vi'][21],i*32:(i+1)*32,i*32:(i+1)*32],axis=(1,2)),
                               np.arange(cat['vo'][21],cat['vi'][21])*1.3-100),i*40,
                color=fns.segmented_spectral(1 - i/9),s=10,linewidth=0.,zorder=30)
    
  ) for i in range(10)]

axs[0].set_xlim(-70,40)
axs[0].set_xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

axs[0].set_aspect(0.4)

axs[0].set_xticks(np.arange(-60,60,30))

axs[0].vlines(np.nanmean(hi_1s[21].data),0,350,color='black',linestyle='--',zorder=-20,linewidth=0.7)
axs[0].set_yticks([])


cat['rot_ax_HI<deg>']


for i in range(0,328,33):
    # print((1-i/300))
    print(i/30)


plt.scatter(cat['radius_GMC<pc>'],cat['sigma_v_GMC<km/s>'])
# plt.scatter(cat['radius_HI<pc>'],cat['sigma_v_HI<km/s>'])

plt.loglog()


oph_hi = hihi[14].data
oph_hi_1row_specs = oph_hi[:,35:36,:].swapaxes(1,2).transpose().reshape(-1,155) #for the first row in pbp, 35:70, and chunk them in 34 in the next line
oph_hi_1row_specs = [np.nanmean(oph_hi_1row_specs[i*1:(i+1)*1],axis=0) for i in range(0,171,3)]


for i in range(10):
    plt.plot(np.arange(155)*1.3-100,oph_hi_1row_specs[i],color='black',linestyle=':',linewidth=3,alpha=0.8) #whole spectrum (-100,100)
    plt.plot(np.arange(155)*1.3-100,np.nanmean(hihi[14].data,axis=(1,2)),color='black',linewidth=0.8,alpha=0.2) #whole spectrum (-100,100)

    plt.plot(np.arange(cat['vo'][14],cat['vi'][14])*1.3-100,
             oph_hi_1row_specs[i][cat['vo'][14]:cat['vi'][14]],
           color='black',linewidth=4) #spectrum just within the selected range
    
    plt.plot(np.arange(155)*1.3-100,100*np.nanmean(coco[14].data,axis=(1,2)),color=fns.segmented_spectral(19),linewidth=1.5) #co spectrum
    
    plt.vlines(cat['vel<km/s>'][14],0,70,linestyle='--',linewidth=0.7,color='red')
    plt.xlim(-70,70)
    plt.show()


plt.scatter(cat['vel_grad_HI<km/s/pc>']/cat['vel_grad_GMC<km/s/pc>'], cat['radius_HI<pc>']/cat['radius_GMC<pc>'] )
plt.loglog()


plt.plot(cat['del_phi<deg>'],marker='+')
plt.hlines([0,-60,60],0,22)
plt.grid()
plt.xticks(np.arange(0,22,2))


cat['rot_ax_GMC<deg>'][13]+90,cat['rot_ax_HI<deg>'][13]+90


360-209.83


cat['vel_grad_HI<km/s/pc>'][13]


plt.plot(cat['vel_grad_GMC<km/s/pc>'],marker='+')


# for i in range(22):
#     fig,ax=plt.subplots(1,1,dpi=200)

#     ax.plot(np.arange(155)*1.3-100,
#             np.nanmean(hihi[i].data,axis=(1,2)),linestyle=':',linewidth=3,
#             alpha=0.8,color='black')

#     ax.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
#             np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),
#             linewidth=4,color='black')


#     if i in [20,12,11]:
#         ax.vlines(np.array([cat['vo'][i],cat['vi'][i]-1])*1.3-100,
#                   0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),
#                   linestyle='--',linewidth=0.7,color='red')

#     elif i in [8,4]:
#         ax.vlines(np.array([cat['vo'][i],cat['vi'][i]-1])*1.3-100,
#                   0,0.6*np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),
#                   linestyle='--',linewidth=0.7,color='red')
#     else:
#         ax.vlines(np.array([cat['vo'][i],cat['vi'][i]-1])*1.3-100,
#                   0,0.4*np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),
#                   linestyle='--',linewidth=0.7,color='red')

#     ax.set_xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

#     ax2 = ax.twinx()

#     ax2.plot(np.arange(155)*1.3-100,np.nanmean(coco[i].data,axis=(1,2)),
#              color=seg_spec(19),linewidth=1.5)


#     if i in [4,5,6,12,13,14,20]:
#         ax.set_title(names[i],loc='left',y=0.87,x=0.7,fontsize='small')
#     else:
#         ax.set_title(names[i],loc='left',y=0.87,x=0.05,fontsize='small')


#     ax.set_ylabel('$\mathrm{T_{B,HI}\ [K]}$')
#     ax.set_xlim(vbnds[i])
#     ax.set_yticks(np.around(np.linspace(0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),5),1),
#                    labels = np.around(np.linspace(0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2))),5),1))
    
#     ax2.set_ylabel('$\mathrm{T_{B,H_2}\ [K]}$',color=seg_spec(19))
#     ax2.set_yticks(np.around(np.linspace(0,np.nanmax(np.nanmean(coco[i].data,axis=(1,2))),5),3),
#                    labels = np.around(np.linspace(0,np.nanmax(np.nanmean(coco[i].data,axis=(1,2))),5),3),
#                    color=seg_spec(19))
#     # plt.savefig('spectrum'+lbls[i],bbox_inches='tight')

#     plt.show()


#fake for now
names = ['']*22
names[0] = 'Aquila South'
names[2] = 'Camelopardalis'
names[4] = 'Canis Major'
names[17] = r'$\lambda-$Orionis'


# fig,axs = plt.subplots
for i in [0,2,4,17,15,21,8]:
    plt.figure(dpi=200)
    plt.hist(np.ravel(np.where(~hi_msks[i],np.nan,hi_1s[i].data)),density=True,bins=70,histtype='step',color='black',linewidth=4,alpha=0.8)
    plt.hist(np.ravel(np.where(~msks[i],np.nan,h2_1s[i].data)),density=True,bins=50,histtype='step',color=fns.segmented_spectral(19),linewidth=1.5)
    
    
    plt.title(cat['name'][i])
    
    plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
    
    
    if i in [0,17]:
        xloc = 0.8
    else:
        xloc = 0.2
        
    plt.title(i)#names[i],y = 0.87,x=xloc)
    
    plt.ylabel(r'$P\ (\mathrm{v_{LSR}})$')
    
    # plt.savefig('mu1_pdf_'+cat['name'][i],bbox_inches='tight')
    plt.show()


cat['vel<km/s>'][5]


cat['distance<pc>'][6]


mad_hi = hihi[11].data


mad_hi_0row_specs = mad_hi[:,:6,:].swapaxes(1,2).transpose().reshape(-1,155) #for the 0th row in pbp, 0:6, and chunk them in 6 in the next line
mad_hi_0row_specs = [np.nanmean(mad_hi_0row_specs[i*6:(i+1)*6],axis=0) for i in range(9)]


for i in range(9):
    plt.plot(np.arange(155)*1.3-100,mad_hi_0row_specs[i],color='black',linestyle=':',linewidth=3,alpha=0.8) #whole spectrum (-100,100)
    plt.plot(np.arange(155)*1.3-100,np.nanmean(hihi[11].data,axis=(1,2)),color='black',linewidth=0.8,alpha=0.2) #whole spectrum (-100,100)

    plt.plot(np.arange(cat['vo'][11],cat['vi'][11])*1.3-100,
             mad_hi_0row_specs[i][cat['vo'][11]:cat['vi'][11]],
           color='black',linewidth=4) #spectrum just within the selected range
    
    plt.plot(np.arange(155)*1.3-100,100*np.nanmean(coco[11].data,axis=(1,2)),color=fns.segmented_spectral(19),linewidth=1.5) #co spectrum
    
    plt.vlines(cat['vel<km/s>'][11],0,70,linestyle='--',linewidth=0.7,color='red')
    plt.xlim(-30,100)
    plt.show()


plt.plot(cat['vel_grad_HI<km/s/pc>'],marker='+')
# plt.yscale('log')


fig,axs=plt.subplots(1,3,dpi=200,figsize=(12,5))
for i,ax in zip([0,3,6],axs.flat):
    ax.plot(np.arange(155)*1.3-100,mad_hi_0row_specs[i],color='black',linestyle=':',linewidth=3,alpha=0.8) #whole spectrum (-100,100)
    ax.plot(np.arange(155)*1.3-100,np.nanmean(hihi[11].data,axis=(1,2)),color='black',linewidth=0.8,alpha=0.2) #whole spectrum (-100,100)

    ax.plot(np.arange(cat['vo'][11],cat['vi'][11])*1.3-100,
             mad_hi_0row_specs[i][cat['vo'][11]:cat['vi'][11]],
           color='black',linewidth=4) #spectrum just within the selected range
    
    ax.plot(np.arange(155)*1.3-100,100*np.nanmean(coco[11].data,axis=(1,2)),color=fns.segmented_spectral(19),linewidth=1.5) #co spectrum
    
    fig.tight_layout()
    ax.set_xlim(-30,90)
    ax.vlines(cat['vel<km/s>'][11],0,70,linestyle='--',linewidth=0.7,color='red')
    
fig.supxlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$',y=-0.05)
fig.supylabel('$T\mathrm{_B\ [K]}$',x=-0.05)
plt.savefig('mad_0throw_patches036',bbox_inches='tight')
    # plt.show()


plt.figure(dpi=200)
plt.plot(np.log2(cat['vel_grad_HI<km/s/pc>']/0.05),marker='+')

plt.hlines([-1,1,0],0,22)
plt.grid()


cat['rot_ax_HI<deg>'][6]


np.tan(117*u.deg)


cat['lon_GMC<deg>'][11]


plt.contourf((np.nansum(hihi[11].data,axis=1)*0.125).transpose(),extent=[-100,100,-140.,-147,],cmap='magma',levels=20)
plt.xlim(-30,90)


plt.contourf((np.nansum(hihi[11].data,axis=2)*0.125).transpose(),extent=[-100,100,-5.,-0.5,],cmap='magma',levels=20)
plt.xlim(-30,90)


36/0.46


cat['lat_GMC<deg>'][4]


plt.contourf((np.nansum(hihi[4].data,axis=1)*0.125).transpose(),extent=[-100,100,-132.5,-138.5,],cmap='magma',levels=20)
plt.xlim(-30,90)


plt.contourf((np.nansum(hihi[4].data,axis=2)*0.125).transpose(),extent=[-100,100,-4,-1.5,],cmap='magma',levels=20)
plt.xlim(-30,90)


cat['lon_GMC<deg>'][20],cat['lat_GMC<deg>'][20]


plt.contourf((np.nansum(hihi[20].data,axis=1)*0.125).transpose(),extent=[-100,100,-151,-155,],cmap='magma',levels=20)
plt.xlim(-30,90)


plt.contourf((np.nansum(hihi[20].data,axis=2)*0.125).transpose(),extent=[-100,100,-3.5,-0.5,],cmap='magma',levels=20)
plt.xlim(-30,90)


plt.hist(cat['vel_grad_HI<km/s/pc>'],alpha=0.5)
# plt.hist(cat['vel_grad_HI<km/s/pc>']*cat['distance<pc>']*u.deg.to(u.rad),alpha=0.5)


plt.plot(cat['vel_grad_HI<km/s/pc>']*cat['distance<pc>']*u.deg.to(u.rad),marker='x')
# plt.plot(cat['vel_grad_HI<km/s/pc>'],marker='+')


for i in range(22):
    lb = np.nansum(hihi[i].data,axis=0)*1.3
    lv = np.nansum(hihi[i].data,axis=1)*0.125
    bv = np.nansum(hihi[i].data,axis=2)*0.125
    
    figs,axs = plt.subplots(1,3,dpi=200,figsize=(13,5))
    for mp,ax in zip([hi_1s[i].data,lv,bv],axs):
        if ax == axs[0]:
            cmp = fns.segmented_spectral
        else:
            cmp = 'magma'
        ax.contourf(mp,cmap=cmp)
        
    axs[1].set_ylim(cat['vo'][i]-5,cat['vi'][i]+5)
    axs[2].set_ylim(cat['vo'][i]-5,cat['vi'][i]+5)
        
    axs[0].set_title(cat['name'][i])
    axs[1].set_title(round(cat['vel_grad_HI<km/s/pc>'][i]*cat['distance<pc>'][i]*u.deg.to(u.rad),3))
    # axs[2].set_title(round(cat['vel_grad_HI<km/s/pc>'][i]*cat['distance<pc>'][i]*u.deg.to(u.rad)/(2*np.nanstd(cat['vel<km/s>'][i])),5))
    plt.show()


hi_1_gal = fits.open('bigh1.fits')[0].data


help(plt.quiver)


plt.figure(dpi=200)
plt.contourf(hi_1_gal,extent=[360,-360,-50,50],cmap = fns.segmented_spectral,levels=np.arange(-40,40))

plt.xlim(90,-90)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'],cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'],
          pivot='middle',width=4e-3,alpha=0.5)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'],cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'],
          pivot='middle',width=1e-3,headwidth=10,headlength=7)


plt.figure(dpi=200)
plt.contourf(hi_1_gal,extent=[360,-360,-50,50],cmap = fns.segmented_spectral,levels=np.arange(-40,60))

plt.xlim(270,90)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'],cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'],
          pivot='middle',width=4e-3,alpha=0.5)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'],cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'],
          pivot='middle',width=1e-3,headwidth=10,headlength=7)


plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)+360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'],cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'],
          pivot='middle',width=4e-3,alpha=0.5)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1)+360,np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'],cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'],
          pivot='middle',width=1e-3,headwidth=10,headlength=7)


halpha = fits.open('../not_rn_6/wham-ss-DR1-v161116-170912-grid.fits')[0]


halpha.header


(np.arange(1440)-720)*0.25


halpha.header


halpha_mom1 = fns.get_mom1(halpha.data,np.arange(162)*2.03 - 169.12,dx = 2.03)


plt.figure(dpi=200)
plt.contourf(halpha_mom1,cmap=fns.segmented_spectral,levels=np.linspace(-50,50,100),extent=[180,-180,-90,90])
# plt.colorbar()

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'],cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'],
          pivot='middle',width=4e-3,alpha=0.5)

# plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
#            -cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'],cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'],
#           pivot='middle',width=1e-3,headwidth=10,headlength=7)


plt.xlim(180,-180)
plt.ylim(-40,40)
plt.gca().set_aspect(1.7)


plt.figure(dpi=200)
plt.contourf(halpha_mom1,cmap=fns.segmented_spectral,levels=np.linspace(-50,50,100),extent=[180,-180,-90,90])
# plt.colorbar()

# plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
#            -cat['dvdx_HI']/cat['vel_grad_HI<km/s/pc>'],cat['dvdy_HI']/cat['vel_grad_HI<km/s/pc>'],
#           pivot='middle',width=4e-3,alpha=0.5)

plt.quiver(np.nanmean(cat['lon_GMC<deg>'],axis=1),np.nanmean(cat['lat_GMC<deg>'],axis=1),
           -cat['dvdx_GMC']/cat['vel_grad_GMC<km/s/pc>'],cat['dvdy_GMC']/cat['vel_grad_GMC<km/s/pc>'],
          pivot='middle',width=1e-3,headwidth=10,headlength=7)


plt.xlim(180,-180)
plt.ylim(-40,40)
plt.gca().set_aspect(1.7)


phillips99_data = pd.read_csv('J_A+AS_134_241_table1.dat.gz.txt',skiprows=9,delimiter='|')


phillips99_data.columns = ['Name', 'Note', 'Ah Am  RAs  Ed Em  DEs', 'l_Dist', 'Dist',
       'u_Dist', 'RV', 'DV   DV2', 'l_Size', 'Size', 'n_Size', 'Size2',
       'Trans', 'L_OMEGA', 'OMEGA__OMEGA2', 'l_GAMMA', 'GAMMA',
       'l_Mass', 'Mass__Mass2', 'MassEm', 'THETA_THETA2', 'OMEGAo', 'u_OMEGAo', 'Type',
       'r_']


#clean this up
#first just wrap up the fit lines through the three (gmn, h2, hi)


pOmegas = (np.array([(float(n.split()[0])) 
                     for i,n in enumerate(phillips99_data['OMEGA__OMEGA2'][1:-1])])\
           *1e-14/u.s).to(u.km/u.s/u.pc).value


pRs = phillips99_data['Size'][1:-1].astype('float')


pMs = np.zeros(len(phillips99_data['Mass__Mass2'][1:-1]))

for i,n in enumerate(phillips99_data['Mass__Mass2'][1:-1]):
    if len(n.split())==0.:
        pMs[i] = 0
    elif len(n.split())==1.:
        pMs[i] = float(n.split()[0])
        
    elif len(n.split())==3.:
        pMs[i] = np.nanmean([float(n.split()[0]),float(n.split()[2])])
              
    elif n.split()[0][-1]=='>':
        pMs[i] = np.nanmean([float(n.split()[0][:6]),float(n.split()[1])])
        
    elif n.split()[1][:2]=='->':
        pMs[i] = np.nanmean([float(n.split()[0]),float(n.split()[1][2:])])
        
    else:
        pMs[i] = np.nanmean([float(n.split()[0]),float(n.split()[1])])


phillips = pd.DataFrame(np.array([pMs,pRs,pOmegas]).transpose(),columns=['M<Msun>','R<pc>','Omega<km/s/pc>'])

phillips = phillips[~np.isinf(np.log10(pMs))]


phillips['j<cm2/s>'] = (np.array(phillips['Omega<km/s/pc>']*phillips['R<pc>']**2)*u.km/u.s*u.pc).cgs.value


curve_fit(line,np.log10(phillips['M<Msun>']),np.log10(phillips['j<cm2/s>']))


curve_fit(line,np.log10(phillips['R<pc>']),np.log10(phillips['j<cm2/s>']))


plt.scatter(phillips['R<pc>'],phillips['j<cm2/s>'],marker='+')
plt.scatter(cat['radius_GMC<pc>'],cat['j_GMC<cm2/s>'],marker='.')
plt.scatter(cat['radius_HI<pc>'],cat['j_HI<cm2/s>'],marker='x')

plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],marker='x')

plt.loglog()


1000*u.AU.to(u.pc)


plt.scatter(phillips['M<Msun>'],phillips['j<cm2/s>'],marker='+')
plt.scatter(cat['mass_GMC<Msun>'],cat['j_GMC<cm2/s>'],marker='.')
plt.scatter(cat['mass_HI<Msun>'],cat['j_HI<cm2/s>'],marker='x')

# plt.scatter(goodman['M<Msun>'],goodman['j<cm2/s>'],marker='x')

plt.loglog()
