#!/usr/bin/env python
# coding: utf-8


def beta_rot(omega,rad,mass):
    top = rad*(omega*rad)**2
    bottom = 3*c.G*mass
    
    return (top/bottom).cgs


# plt.plot(beta_rot(hi_omg,hi_rad,hi_mass),marker='.')
# plt.plot(beta_rot(h2_omg,h2_rad,h2_mass),marker='.')
plt.plot(h2_omg)
plt.grid()
plt.xticks(np.arange(0,22,3))


plt.show()


plt.plot(h2_rad/h2_rad_geom,marker='.')
plt.grid()
plt.xticks(np.arange(0,22,3))
plt.show()


[(plt.contourf(h2_msks[i]),
plt.show()
) for i in range(22)]


h2_galphi[17]


hi_galphi[17]


plt.figure(dpi=200)
plt.plot(beta_rot(hi_omg,hi_rad,hi_mass),marker='.')
plt.plot(beta_rot(h2_omg,h2_rad,h2_mass),marker='.')

# plt.grid()
plt.xticks(np.arange(0,22,3))
plt.hlines(0.1,0,22)
plt.yscale('log')
plt.ylim(1e-2,1e0)
plt.show()


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
from astropy.modeling import models,fitting

from scipy.optimize import curve_fit

from gaussdecomp import spectrum,fitter

plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


def line(x,m,b):
    return m*x+b


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
    'Taurus'  
]


hihi = fits.open('hihi.fits')
hi_msks = [(fits.open('hi_masks.fits')[i].data-1).astype('bool') for i in range(22)]

h2h2 = fits.open('h2h2.fits')
h2_msks = [(fits.open('h2_masks.fits')[i].data-1).astype('bool') for i in range(22)]


bounds = pd.read_csv('lbv_bounds_both.csv')


vos,vis = bounds['vo'],bounds['vi']


cat = QTable(fits.open('../rn/sighhplease.fits')[1].data)


h2_dvdx


h2_mass = np.array([fns.mass_n_area(h2h2[i].data*u.K,
                 cat['distance<pc>'][i]*u.pc,
                 'CO',
                 mask=h2_msks[i])[0].to(u.solMass).value
 
                 for i in range(22)])*u.solMass

hi_mass = np.array([fns.mass_n_area(hihi[i].data[vos[i]:vis[i],:,:]*u.K,
                 cat['distance<pc>'][i]*u.pc,
                 '21cm',
                 mask=hi_msks[i])[0].to(u.solMass).value
 
                 for i in range(22)])*u.solMass

#h2 is across the whole +/- 100km/s bc the idea is that we calculate that first,
#then we get r_hi from mh2, then we loko within that aperture to determine what the velocity bounds are


h2_rad = np.array([np.sqrt(fns.mass_n_area(h2h2[i].data*u.K,
                 cat['distance<pc>'][i]*u.pc,
                 'CO',
                 mask=h2_msks[i])[1]/np.pi).value
 
                 for i in range(22)])*u.pc

hi_rad = np.cbrt(0.75*h2_mass/(0.57/u.cm**3*1.36*c.m_p*np.pi)).to(u.pc)


lcs,bcs = np.nanmean(cat['lon_GMC<deg>'],axis=1), np.nanmean(cat['lat_GMC<deg>'],axis=1)

hi_r_deg = np.arctan(hi_rad.value/cat['distance<pc>'])*u.rad.to(u.deg)


h2_mom0 = [np.where(~h2_msks[i],np.nan,
                    np.nansum(h2h2[i].data[vos[i]:vis[i],:,:],axis=0)*1.3*1.97e20) 
           for i in range(22)]


hi_mom0 = [np.where(~hi_msks[i],np.nan,
                    np.nansum(hihi[i].data[vos[i]:vis[i],:,:],axis=0)*1.3*1.82e18) 
           for i in range(22)]


h2_mom1 = [np.where(~h2_msks[i],np.nan,
          fns.get_mom1(h2h2[i].data[vos[i]:vis[i],:,:],np.arange(vos[i],vis[i])*1.3-100))
 
 for i in range(22)]


h2_mom2 = [np.where(~h2_msks[i],np.nan,
          np.sqrt(fns.get_momN(h2h2[i].data[vos[i]:vis[i],:,:],np.arange(vos[i],vis[i])*1.3-100,2)))
 
 for i in range(22)]


hi_mom1 = [np.where(~hi_msks[i],np.nan,
          fns.get_mom1(hihi[i].data[vos[i]:vis[i],:,:],np.arange(vos[i],vis[i])*1.3-100))
 
 for i in range(22)]


hi_mom2 = [np.where(~hi_msks[i],np.nan,
          np.sqrt(fns.get_momN(hihi[i].data[vos[i]:vis[i],:,:],np.arange(vos[i],vis[i])*1.3-100,2)))
 
 for i in range(22)]


hi_mom1_fitpars = np.array([
    fns.fit_gradient(hi_mom1[i],
                     cat['distance<pc>'][i],
                     weights=hi_mom2[i]) for i in range(22)])


h2_mom1_fitpars = np.array([
    fns.fit_gradient(h2_mom1[i],
                     cat['distance<pc>'][i],
                     weights=h2_mom2[i]) for i in range(22)])


hi_dvdx,hi_dvdy = hi_mom1_fitpars.transpose()[:2]
h2_dvdx,h2_dvdy = h2_mom1_fitpars.transpose()[:2]

hi_omg = np.linalg.norm(hi_mom1_fitpars.transpose()[:2],axis=0)*u.km/u.s/u.pc
h2_omg = np.linalg.norm(h2_mom1_fitpars.transpose()[:2],axis=0)*u.km/u.s/u.pc


hi_sigmav = np.array([np.sqrt(fns.get_momN(np.nanmean(hihi[i].data[vos[i]:vis[i],:,:],axis=(1,2)),
            np.arange(vos[i],vis[i])*1.3-100, 2))

for i in range(22)])*u.km/u.s


h2_sigmav = np.array([np.sqrt(fns.get_momN(np.nanmean(h2h2[i].data[vos[i]:vis[i],:,:],axis=(1,2)),
            np.arange(vos[i],vis[i])*1.3-100, 2))

for i in range(22)])*u.km/u.s

hi_j = (hi_omg*hi_rad**2).cgs
h2_j = (h2_omg*h2_rad**2).cgs

hi_J = (hi_j*hi_mass).cgs
h2_J = (h2_j*h2_mass).cgs


sigd_1 = np.array([2,2, 8, 3,4, 3,4,4,1,5,9,6,3,3,2,2,2,1,2,1,5,2])
sigd_2 = np.array([7,24,11,9,60,8,44,89,11,25,9,104,37,39,7,22,21,20,15,8,65,7])

sigd = np.linalg.norm([sigd_1,sigd_2],axis=0)


sigdvx_hi,sigdvy_hi = np.array([fns.fit_gradient(hi_mom1[i],
            cat['distance<pc>'][i],
            weights=hi_mom2[i],
                 full=True)['1_sigma_errors'][:2] for i in range(22)]).transpose()*u.km/u.s/u.pc


sigdvx_h2,sigdvy_h2 = np.array([fns.fit_gradient(h2_mom1[i],
            cat['distance<pc>'][i],
            weights=h2_mom2[i],
                 full=True)['1_sigma_errors'][:2] for i in range(22)]).transpose()*u.km/u.s/u.pc


h2_sigr = (sigd/cat['distance<pc>'])*h2_rad
h2_sigm = (2*sigd/cat['distance<pc>'])*h2_mass
hi_sigr = (1/3*h2_sigm/h2_mass)*hi_rad
hi_sigm = (2*sigd/cat['distance<pc>'])*hi_mass

hi_sigomg = np.linalg.norm([sigdvx_hi*hi_dvdx,sigdvy_hi*hi_dvdy],axis=0)/hi_omg.value
h2_sigomg = np.linalg.norm([sigdvx_h2*h2_dvdx,sigdvy_h2*h2_dvdy],axis=0)/h2_omg.value

hi_sigj = np.linalg.norm([hi_sigomg/hi_omg, 2*hi_sigr/hi_rad],axis=0)*hi_j
h2_sigj = np.linalg.norm([h2_sigomg/h2_omg, 2*h2_sigr/h2_rad],axis=0)*h2_j


phillips99_data = pd.read_csv('../rn/J_A+AS_134_241_table1.dat.gz.txt',skiprows=9,delimiter='|')

phillips99_data.columns = ['Name', 'Note', 'Ah Am  RAs  Ed Em  DEs', 'l_Dist', 'Dist',
       'u_Dist', 'RV', 'DV   DV2', 'l_Size', 'Size', 'n_Size', 'Size2',
       'Trans', 'L_OMEGA', 'OMEGA__OMEGA2', 'l_GAMMA', 'GAMMA',
       'l_Mass', 'Mass__Mass2', 'MassEm', 'THETA_THETA2', 'OMEGAo', 'u_OMEGAo', 'Type',
       'r_']

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
        
pTypes = phillips99_data['Type'][1:-1]   
phillips = pd.DataFrame(np.array([pTypes,pMs,pRs,pOmegas]).transpose(),columns=['type','M<Msun>','R<pc>','Omega<km/s/pc>'])

phillips = phillips[~np.isinf(np.log10(pMs))] #only clouds with real masses

#selecting out geometries that arent rotationally symmetric
phillips = phillips[~(phillips['type']=='D/R ')]
phillips = phillips[~np.array([('F' in  nm) for nm in phillips['type']])]

phillips['j<cm2/s>'] = (np.array(phillips['Omega<km/s/pc>']*phillips['R<pc>']**2)*u.km/u.s*u.pc).cgs.value


all_rs = np.concatenate((phillips['R<pc>'],h2_rad.value,hi_rad.value)).astype('float')
all_js = np.concatenate((0.4*phillips['j<cm2/s>'],0.4*h2_j.value,0.4*hi_j.value)).astype('float')
all_ms = np.concatenate((phillips['M<Msun>'],h2_mass.value,hi_mass.value)).astype('float')


goodman = pd.read_csv('../rn/goodman_1993_data.csv')

chenx_rs = np.array([1703,1232,3530,1440,1272,1763,1329,2762])*u.AU.to(u.pc)
chenx_js = (np.array([0.12,0.12,0.45,0.13,0.24,0.12,0.42,0.45])*1e-3*u.km/u.s*u.pc).cgs

chenh_rs = np.array([0.038,0.033,0.025,0.032,0.026,0.026,
                     0.034,0.054,0.042,0.020,0.022,0.024,0.19,0.028,0.045,
                    0.044,0.050,0.048])

chenh_js = (np.array([1.3,1.3,0.36,1.3,0.8,0.52,0.59,2.5,1.,0.8,1.5,
                      1.6,6.,0.37,0.25,0.34,0.46,0.54])*1e-3 * u.km/u.s*u.pc).cgs

    
r_xs = np.linspace(-1.9,2.4,10)
m_xs = np.linspace(-1,6.5,10)


core_rxs = 10**np.linspace(-1.93,1.92,10)
core_mxs = 10**np.linspace(-1.,5.18,10)

h2_rxs = 10**np.linspace(0.67,1.99,10)
h2_mxs = 10**np.linspace(2.91,6.07,10)

hi_rxs = 10**np.linspace(1.33,2.39,10)
hi_mxs = 10**np.linspace(3.67,6.54,10)

all_rxs = 10**np.linspace(-1.93,2.39,10)
all_mxs = 10**np.linspace(-1.,6.54,10)


hi_theta = np.array([fns.dxdy_2_theta(hi_dvdx[i],hi_dvdy[i]) for i in range(22)])
h2_theta = np.array([fns.dxdy_2_theta(h2_dvdx[i],h2_dvdy[i]) for i in range(22)])


delphi = hi_theta-h2_theta

for i in range(22):
    if delphi[i] < -180:
        delphi[i] += 360
        
    if delphi[i] > 180:
        delphi[i] -= 360
        
        
hi_phi = hi_theta - 90
for i in range(22):
    if hi_phi[i] < -180:
        hi_phi[i] += 360
        
        
h2_phi = h2_theta - 90
for i in range(22):
    if h2_phi[i] < -180:
        h2_phi[i] += 360
        
        
hi_galphi = np.cos((hi_phi-90)*u.deg)*np.cos(bcs*u.deg)
h2_galphi = np.cos((h2_phi-90)*u.deg)*np.cos(bcs*u.deg)


[thing.unit for thing in [h2_mass, hi_mass,
h2_rad, hi_rad,
h2_omg, hi_omg]]


h2_J


quartnity = pd.DataFrame(np.array([
    h2_mass, hi_mass,
    h2_rad, hi_rad,
    h2_omg*1e3, hi_omg*1e3,
    h2_J,hi_J
]).transpose(),
    columns = ['h2m<msun>','him<msun>',
               'h2r<pc>','hir<pc>',
               'h2o<m/s/pc>','hio<m/s/pc>',
              'h2J<erg.s>','hiJ<erg.s>'])


quartnity.to_csv('ang_mom_stuff.csv',index=False)


h2_color = cm.YlGnBu_r.resampled(10)(6)
hi_color = cm.YlGnBu_r.resampled(10)(1)

err_color = fns.segmented_spectral(3)
dat_color = 'black'


err_kws = {'ms':3,'elinewidth':0.5,
           'ecolor':dat_color,'capsize':1}


ftline_kws = {'linewidth':0.7,'linestyle':'--'}
ftfill_kws = {'linewidth':0.,'alpha':0.3}


# plt.figure(dpi=200)

fig,axs=plt.subplots(1,2,dpi=200,figsize=(10,5))

x = np.linspace(3e0,4e2,20)


axs[0].plot(x,x**1.54 * 10**(22.56),color=dat_color,
            linewidth=0.75,linestyle=':',label=r'$j\propto R^{1.54\pm0.14}$')

# axs[0].plot(x,(x**1.65)*10**(22.56),color=hi_color,**ftline_kws)

axs[0].fill_between(x,x**(1.54-0.14) * 10**(22.56),
                      x**(1.54+0.14) * 10**(22.56),
                    color=h2_color,**ftfill_kws)

# axs[0].scatter(h2_rad,0.4*h2_j,color='black',zorder=20,s=10)
axs[0].errorbar(h2_rad,0.4*h2_j,
                xerr=h2_sigr,yerr=0.4*h2_sigj,
                fmt='o',**err_kws,color=dat_color)


axs[0].loglog()


axs[0].set_ylabel(r'$j\ \mathrm{[cm^2\ s^{-1}]}$')
axs[0].set_xlabel(r'$R\mathrm{_{H_2}\ [pc]}$')


axs[0].set_xlim(4e0,130)
axs[0].set_ylim(1.5e23,3e25)

axs[0].legend()


axs[1].plot(x,(x**1.65)*10**(21.98),color=dat_color,
            linewidth=0.75,linestyle='--',label=r'$j\propto R^{1.65\pm0.32}$')

# axs[1].plot(x,x**1.54 * 10**(21.98),color=h2_color,**ftline_kws)

axs[1].fill_between(x,x**(1.65-0.32) * 10**(21.98),
                      x**(1.65+0.32) * 10**(21.98),
                    color=cm.YlGnBu_r.resampled(10)(1),**ftfill_kws)


# axs[1].scatter(hi_rad,0.4*hi_j,color='black',zorder=20,s=10)
axs[1].errorbar(hi_rad,0.4*hi_j,
                xerr=hi_sigr,yerr=0.4*hi_sigj,
                fmt='s',**err_kws,color=dat_color)


axs[1].loglog()
axs[1].legend()

axs[1].set_xlim(17,3e2)
axs[1].set_ylim(1.5e23,1.5e26)

axs[1].set_xlabel(r'$R\mathrm{_{HI}\ [pc]}$')


fig.tight_layout()

plt.savefig('jR_hih2',bbox_inches='tight')


# plt.figure(dpi=200)

fig,axs=plt.subplots(1,2,dpi=200,figsize=(10,5))

x = np.linspace(5e2,1e7,20)


axs[0].plot(x,x**0.62 * 10**(21.72),color=h2_color,
            linewidth=0.75,linestyle='--',label=r'$j\propto M^{0.62\pm0.07}$')

axs[0].plot(x,(x**0.63)*10**(21.71),color=hi_color,**ftline_kws)

axs[0].fill_between(x,x**(0.62-0.07) * 10**(21.72),
                      x**(0.62+0.07) * 10**(21.72),
                    color=h2_color,**ftfill_kws)

# axs[0].scatter(h2_mass,0.4*h2_j,color='black',zorder=20,s=10,marker='o')
axs[0].errorbar(h2_mass,0.4*h2_j,xerr=h2_sigm,yerr=0.4*h2_sigj,fmt='o',**err_kws,color=dat_color)

axs[0].loglog()


axs[0].set_ylabel(r'$j\ \mathrm{[cm^2\ s^{-1}]}$')
axs[0].set_xlabel(r'$M\mathrm{_{H_2}\ [pc]}$')


axs[0].set_xlim(5e2,2e6)
axs[0].set_ylim(1.5e23,3e25)


axs[0].legend()


axs[1].plot(x,(x**0.63)*10**(21.72),color=hi_color,
            linewidth=0.75,linestyle='--',label=r'$j\propto M^{0.63\pm0.10}$')

axs[1].plot(x,x**0.62 * 10**(21.72),color=h2_color,**ftline_kws)

axs[1].fill_between(x,x**(0.63-0.10) * 10**(21.72),
                      x**(0.63+0.10) * 10**(21.72),
                    color=hi_color,**ftfill_kws)


# axs[1].scatter(hi_mass,0.4*hi_j,color='black',zorder=20,s=10,marker='s')
axs[1].errorbar(hi_mass,0.4*hi_j,xerr=hi_sigm,yerr=0.4*hi_sigj,fmt='s',**err_kws,color=dat_color)


axs[1].loglog()
axs[1].legend()

axs[1].set_xlim(3e3,8e6)
axs[1].set_ylim(1.5e23,1.5e26)

axs[1].set_xlabel(r'$M\mathrm{_{HI}\ [pc]}$')


fig.tight_layout()

# plt.savefig('jM_hih2',bbox_inches='tight')


curve_fit(line,np.log10(phillips['R<pc>'].astype('float')),
          np.log10(0.4*phillips['j<cm2/s>'].astype('float')))


curve_fit(line,np.log10(goodman['R_eff<pc>']),np.log10(goodman['j<cm2/s>']))[0]


np.sqrt(np.diag(curve_fit(line,np.log10(goodman['R_eff<pc>']),np.log10(goodman['j<cm2/s>']))[1]))


curve_fit(line,np.log10(chenx_rs),np.log10(chenx_js.value))


np.sqrt(np.diag(curve_fit(line,np.log10(chenx_rs),np.log10(chenx_js.value))[1]))


np.sqrt(np.diag(curve_fit(line,np.log10(chenh_rs),np.log10(chenh_js.value))[1]))


curve_fit(line,np.log10(chenh_rs),np.log10(chenh_js.value))


np.log10(1e3*u.AU.to(u.pc))


all_rs = np.concatenate((
    phillips['R<pc>'].astype('float'),
                         h2_rad.value,hi_rad.value,
                         chenh_rs,
                         chenx_rs,
                         goodman['R_eff<pc>']
))
    
all_js = np.concatenate((
    0.4*phillips['j<cm2/s>'].astype('float'),
                         0.4*h2_j.value,0.4*hi_j.value,
                         chenh_js.value,
                         chenx_js.value,
                         goodman['j<cm2/s>']
))


curve_fit(line,np.log10(goodman['R_eff<pc>']),np.log10(goodman['j<cm2/s>']))


np.sqrt(np.diag(curve_fit(line,np.log10(goodman['R_eff<pc>']),np.log10(goodman['j<cm2/s>']))[1]))


curve_fit(line,np.log10(all_rs),np.log10(all_js))


np.sqrt(np.diag(curve_fit(line,np.log10(all_rs),np.log10(all_js))[1]))


clumps = phillips.query('type == "CL  "')
condensations = phillips.query('type == "MIC " | type == "MSC "')
clouds = phillips.query('type == "MI  " | type == "MS  "')
complexes = phillips.query('type == "ML  " | type == "MLC "')


p_typs = [clumps,clouds,condensations,complexes],
plbls = ['Clumps [P99]','Clouds [P99]','Condensations [P99]','Complexes [P99]']
psmbls = ['v','D','2','x']


plt.figure(dpi=200,figsize=(7,4))

xs = np.linspace(1e-3,1e3,10)

plt.fill_between(xs,xs**(1.38+0.03)*10**22.87,
                 xs**(1.38-0.03)*10**22.87,
                 linestyle='--',linewidth=0.7,color='black')
plt.plot(xs,xs**1.5*10**22.87,linestyle='--',linewidth=0.7,color='red')


s = 5

plt.scatter(h2_rad,0.4*h2_j,marker='o',s=s,
             label='Molecular Clouds [This Work]',color=h2_color)

plt.scatter(hi_rad,0.4*hi_j,marker='s',s=s,label='Atomic Envelopes [This Work]',color=hi_color)

plt.scatter(chenh_rs,chenh_js,s=s,color='orange',label='Droplets [C19]')
plt.scatter(chenx_rs,chenx_js,s=s,color='blue',label='Protostellar Cores [C07]')

[plt.scatter(df['R<pc>'],df['j<cm2/s>'],
             marker=psmbl,s=s,color=dat_color,label=plbl) 
 for df,plbl,psmbl in zip([clumps,clouds,condensations,complexes],plbls,psmbls)]
   
plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],s=s,color='red',label='Cores [G93]')


# plt.scatter(h2_rad,h2_j,marker='.')
# plt.scatter(hi_rad,hi_j,marker='x')


plt.loglog()
plt.ylabel(r'$j\ \mathrm{[cm^2\ s^{-1}]}$')
plt.xlabel(r'$R\mathrm{\ [pc]}$')

plt.legend(markerscale=2,fontsize='x-small',loc='upper left',bbox_to_anchor=(-0.1, -0.2),ncols=3)
plt.xlim(4e-3,280)
plt.ylim(2e19,7e26)


help(plt.legend)


curve_fit(line,np.log10(h2_rad.value),np.log10(0.4*h2_j.value))


core_rxs = 10**np.linspace(-2.,1.,10)

h2_rxs = 10**np.linspace(0.67,1.99,10)

hi_rxs = 10**np.linspace(1.33,2.39,10)

cloud_rxs = 10**np.linspace(-1.41,0,10)

cndns_rxs = 10**np.linspace(-0.72,0.46,10)

all_rxs = 10**np.linspace(-3,3.5,10)


plt.figure(dpi=200)

xs = np.linspace(1e-2,1e3,10)


s = 10

# plt.plot(cloud_rxs,cloud_rxs**1.47*10**23.06,**ftline_kws,color=dat_color)
# plt.plot(cndns_rxs,cndns_rxs**1.40*10**23.06,**ftline_kws,color=dat_color)


# plt.plot(hi_rxs,hi_rxs**1.65*10**21.98,**ftline_kws,color=hi_color)
# plt.fill_between(hi_rxs,hi_rxs**(1.65+0.32)*10**21.98,
#                  hi_rxs**(1.65-0.32)*10**21.98,
#                    **ftfill_kws,color=hi_color)

# plt.plot(h2_rxs,h2_rxs**1.54*10**22.56,**ftline_kws,color=h2_color)
# plt.fill_between(h2_rxs,h2_rxs**(1.53+0.14)*10**22.56,
#                  h2_rxs**(1.53-0.14)*10**22.56,
#                    **ftfill_kws,color=hi_color)

# plt.fill_between(core_rxs,core_rxs**(1.52+0.06)*10**22.7,
#                          core_rxs**(1.52-0.06)*10**22.7,
#                 **ftfill_kws,color=dat_color)

plt.plot(all_rxs,all_rxs**1.5*10**23,linewidth=0.75,linestyle='-.',
         color=dat_color,label=r'$j \propto R^{1.5}$')

# plt.plot(all_rxs,all_rxs**1.52*10**23,**ftline_kws,color=dat_color)


plt.scatter(h2_rad,0.4*h2_j,marker='o',s=s,
             label='Molecular Clouds [This Work]',color=h2_color)

plt.scatter(hi_rad,0.4*hi_j,marker='s',s=s,
            label='Atomic Envelopes [This Work]',color=hi_color)


# [plt.scatter(df['R<pc>'],df['j<cm2/s>'],
#              marker=psmbl,s=s,color=dat_color,label=plbl) 
#  for df,plbl,psmbl in zip([clouds,condensations],
#                           ['Clouds [P99]','Condensations [P99]'],
#                           ['2','v'])]


# [plt.scatter(df['R<pc>'],df['j<cm2/s>'],
#              marker=psmbl,s=3,color=dat_color,label=plbl) 
#  for df,plbl,psmbl in zip([clouds],plbls[1],psmbls[1])]

    
# plt.scatter(h2_rad,h2_j,marker='.')
# plt.scatter(hi_rad,hi_j,marker='x')


plt.scatter(chenh_rs,chenh_js,s=s,color='orange',label='Droplets [C19]')
plt.scatter(chenx_rs,chenx_js,s=s,color='blue',label='Protostellar Cores [C07]')

[plt.scatter(df['R<pc>'],df['j<cm2/s>'],
             marker=psmbl,s=s,color=dat_color,label=plbl) 
 for df,plbl,psmbl in zip([clumps,clouds,condensations,complexes],plbls,psmbls)]
   
plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],s=s,color='red',label='Cores [G93]')


plt.loglog()
plt.ylabel(r'$j\ \mathrm{[cm^2\ s^{-1}]}$')
plt.xlabel(r'$R\mathrm{\ [pc]}$')

plt.legend(markerscale=1.5,fontsize='xx-small',loc='lower right')
# markerscale=2,fontsize='x-small',
#            bbox_to_anchor=(0.06, -0.15),
#            loc='upper left',ncols=2)
plt.xlim(4e-3,300)
plt.ylim(2e19,8e26)

plt.savefig('jR_all',bbox_inches='tight')


curve_fit(line,np.log10(clouds['R<pc>'].astype('float')),np.log10(clouds['j<cm2/s>'].astype('float')))


curve_fit(line,np.log10(condensations['R<pc>'].astype('float')),
          np.log10(condensations['j<cm2/s>'].astype('float')))


curve_fit(line,np.log10(clumps['R<pc>'].astype('float')),
          np.log10(clumps['j<cm2/s>'].astype('float')))


curve_fit(line,np.log10(complexes['R<pc>'].astype('float')),
          np.log10(complexes['j<cm2/s>'].astype('float')))


hi_rad


h0_lvls = np.array((
    [20.4,21.3],
    [20.8,21.8],
    [20.7,21.2],
    [20.7,21.3],
    [21.1,21.5],
    [20.5,21.],
    [21.1,21.5],
    [21.2,21.9],
    [20.7,21.5],
    [20.9,21.4],
    [20.7,22.],
    [20.7,21.6],
    [20.8,21.7],
    [21,21.6],
    [20.6,22.],
    [20,21.5],
    [20.7,21.8],
    [20.8,21.7],
    [20.6,21.4],
    [20.,21.2],
    [21.3,21.7],
    [20.7,21.8]))
    
    
h1_lvls = np.array((
    [-0.5,4.3], #0
    [-8,6.5],
    [-2.,-0.3],
    [-0.5,3.5], #3
    [14,18.3],
    [3.7,6.2],
    [8.7,13.], #6
    [2.,6.5],
    [3.7,11.5],
    [-1.3,1.7], #9
    [-1.,7],
    [19.7,28.1],
    [4.5,9], #12
    [7,12],
    [-0.5,6.],
    [5,9], #15
    [5.,10.5],
    [3.,8.5],
    [1.5,6.7],#18
    [-6.5,1],
    [10.5,14.5],
    [-0.5,9])) #21


c1_lvls = np.array((
    [1,4.5], #0
    [-8.5,6.5],
    [-2.5,0.5],
    [1,5.5],#3
    [11,19],
    [4.9,6.4],
    [10,16],#6
    [-3,9.5],
    [4.5,14],
    [-4.5,4],#9
    [2.5,9.7],
    [21.5,33],
    [0,10.5],#12
    [8,14],
    [-1,9],
    [3,12.5],#15
    [4,11.5],
    [-3,13.5],
    [-0.,9.3],#18
    [-5,-0.5],
    [7,19],
    [4,8.7])) #21
    

ytcks = [np.around(np.linspace(bounds['boh'][i]*0.125-50,
                               bounds['bih'][i]*0.125-50,3),1) 
                   for i in range(22)]

xtcks = [np.around(np.linspace(bounds['loh'][i]*-0.125+360,
                               bounds['lih'][i]*-0.125+360,3),1)
                   for i in range(22)]


for i in [7,12,14,15,17,]:
    h2_msks[i] = h2_msks[i]*(np.log10(h2_mom0[i])>20.)


[(
    plt.figure(dpi=200),
    plt.contourf(np.where(~hi_msks[i],np.nan,hi_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],27),
                extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                       bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50],
                 extend='both',
                ),
    
    plt.gca().set_aspect(1),
    plt.colorbar(shrink=0.6,ticks=np.around(np.linspace(h1_lvls[i][0],h1_lvls[i][1],5),1),
                label='$\overline{v\mathrm{_{LSR}}}\ \ \mathrm{[km\ s^{-1}]}$',
                pad=0.05),
    
    
#     plt.contour(np.log10(np.where(~h2_msks[i],np.nan,h2_mom0[i])),colors=dat_color,
#             levels=7,
#                 extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
#                        bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50],
#                  extend='both',linewidths=0.5
#                 ),
    
    
    plt.xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,),
    
#     plt.title(names[i],y=1.1),
    plt.xlabel('Galactic Longitude',y=-1.15),
    plt.ylabel('Galactic Latitude',x=-1.15),
    plt.gca().spines[:].set_visible(False),
    
    
    plt.yticks(ytcks[i],[str(tck)+'$\degree$' for tck in ytcks[i]]),
    plt.xticks(xtcks[i],[str(tck)+'$\degree$' for tck in xtcks[i]]),
    
    
    plt.gca().tick_params(size=1.5),
    
    
    plt.savefig(cat['name'][i]+'_hi_mom1_map',bbox_inches='tight'),
    plt.show()

) for i in range(22)]


lcs,bcs = np.nanmean(cat['lon_GMC<deg>'],axis=1), np.nanmean(cat['lat_GMC<deg>'],axis=1)


hi_r_deg = np.arctan(hi_rad.value/cat['distance<pc>'])*u.rad.to(u.deg)


[(
    plt.figure(dpi=200),
    plt.contourf(np.where(~hi_msks[i],np.nan,hi_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],27),
                extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                       bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50],
                 extend='both',alpha=0.3,antialiased=True,
                ),
    
    plt.contourf(np.where(~h2_msks[i],np.nan,h2_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(c1_lvls[i][0],c1_lvls[i][1],27),
                extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                       bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50],
                 extend='both',linewidths=0.5
                ),
    
    
    plt.gca().set_aspect(1),
    plt.colorbar(shrink=0.6,ticks=np.around(np.linspace(c1_lvls[i][0],c1_lvls[i][1],5),1),
                label='$\overline{v\mathrm{_{LSR}}}\ \ \mathrm{[km\ s^{-1}]}$',
                pad=0.05),
    
    
    plt.title(names[i],y=1.1),
    plt.xlabel('Galactic Longitude',y=-1.15),
    plt.ylabel('Galactic Latitude',x=-1.15),
    plt.gca().spines[:].set_visible(False),
    
    plt.gca().invert_xaxis(),
    
    
    plt.yticks(ytcks[i],[str(tck)+'$\degree$' for tck in ytcks[i]]),
    plt.xticks(xtcks[i],[str(tck)+'$\degree$' for tck in xtcks[i]]),
    
    
    plt.gca().tick_params(size=1.5),
    
    
#     plt.savefig(cat['name'][i]+'_h2_mom1_map',bbox_inches='tight'),
    plt.show()

    
) for i in range(22)]


#mentally prepare to change the scaling factor so it is custom for each cloud, but for now we are going with 100x

fig,axs=plt.subplots(4,3,dpi=200,figsize=(10,10),sharey=True,sharex=True)

[(
    axs.flat[i].plot(np.arange(155)*1.3-100,
                     np.nanmean(hihi[i].data,axis=(1,2)),
                     linestyle='--',color=hi_color,alpha=0.5,linewidth=2),
    
    axs.flat[i].plot(np.arange(vos[i],vis[i])*1.3-100,
                     np.nanmean(hihi[i].data[vos[i]:vis[i],:,:],axis=(1,2)),
                     color=hi_color,linewidth=2),
    
    axs.flat[i].plot(np.arange(155)*1.3-100,
                     40*np.nanmean(h2h2[i].data,axis=(1,2)),
                     linestyle=':',color=h2_color,alpha=0.5),
    
    axs.flat[i].plot(np.arange(vos[i],vis[i])*1.3-100,
                     40*np.nanmean(h2h2[i].data[vos[i]:vis[i],:,:],axis=(1,2)),
                     color=h2_color),
    
#     axs.flat[i-12].set_xlim((vos[i]-25)*1.3-100,(vis[i]+15)*1.3-100),
    axs.flat[i-12].set_title(names[i],y=0.8,x=0.25,fontsize='small')
    
    
) for i in range(12)]


for ax in axs.flat[-2:]:
#     ax.set_axis_off()
    ax.set_xlim(-30,45)
    
fig.supxlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
fig.supylabel('$T\mathrm{_B\ [K]}$')
    
fig.tight_layout()
# plt.savefig('postage_spectra_1',bbox_inches='tight')


#mentally prepare to change the scaling factor so it is custom for each cloud, but for now we are going with 100x

fig,axs=plt.subplots(4,3,dpi=200,figsize=(10,10),sharey=True,sharex=True)

[(
    axs.flat[i-12].plot(np.arange(155)*1.3-100,
                     np.nanmean(hihi[i].data,axis=(1,2)),
                     linestyle='--',color=hi_color,alpha=0.5,linewidth=2),
    
    axs.flat[i-12].plot(np.arange(vos[i],vis[i])*1.3-100,
                     np.nanmean(hihi[i].data[vos[i]:vis[i],:,:],axis=(1,2)),
                     color=hi_color,linewidth=2),
    
    axs.flat[i-12].plot(np.arange(155)*1.3-100,
                     40*np.nanmean(h2h2[i].data,axis=(1,2)),
                     linestyle=':',color=h2_color,alpha=0.5),
    
    axs.flat[i-12].plot(np.arange(vos[i],vis[i])*1.3-100,
                     40*np.nanmean(h2h2[i].data[vos[i]:vis[i],:,:],axis=(1,2)),
                     color=h2_color),
    
#     axs.flat[i-12].set_xlim((vos[i]-25)*1.3-100,(vis[i]+15)*1.3-100),
    axs.flat[i-12].set_title(names[i],y=0.8,x=0.25,fontsize='small')
    
    
) for i in range(12,22)]


for ax in axs.flat[-2:]:
    ax.set_axis_off()
    ax.set_xlim(-30,45)
    
fig.supxlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
fig.supylabel('$T\mathrm{_B\ [K]}$')
    
fig.tight_layout()
# plt.savefig('postage_spectra_2',bbox_inches='tight')


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))

bns = 10**np.arange(-1.,3.,1/4)


axs[0].hist(hi_j/h2_j,bins=bns,
        color=[0,0,0,0.85],edgecolor='white',linewidth=5,)

axs[0].vlines(1,0,6.5,linestyle='--',color=err_color,linewidth=2.5)

axs[0].set_xscale('log')
axs[0].set_xlim(3e-1,4e1)

axs[0].set_xlabel('$j_\mathrm{HI}\ /\ j_\mathrm{H_2}$')
axs[0].set_ylabel('N (clouds per bin)')
axs[0].set_ylim(0,6.5)


axs[1].hist(hi_J/h2_J,bins=bns,
        color=[0,0,0,0.85],edgecolor='white',linewidth=5)

axs[1].vlines(1,0,5.5,linestyle='--',color=err_color,linewidth=2.5)

axs[1].set_xscale('log')
axs[1].set_xlim(.9,3e2)
axs[1].set_xlabel('$J_\mathrm{HI}\ /\ J_\mathrm{H_2}$')
axs[1].set_ylim(0,5.5)
# axs[1].set_yticks([0,1,2,3,4,5,6])

# plt.savefig('jJ_ratios',bbox_inches='tight')


h2_ss = [5,1e1,3,5,3,3,3,5,5,3,1e1,5,5,5,5e1,5,5,5,5,5,3,1e1]
hi_ss = [1e2,1e3,1e1,1e2,1e1,1e2,1e1,1e1,1e2,1e1,1e3,1e1,1e2,1e2,1e3,1e2,1e2,1e2,1e2,1e2,1e1,1e3]


for i in range(22):
    ycoords2, xcoords2 = np.indices(h2_mom1[i].data.shape)
    ycoordsi, xcoordsi = np.indices(hi_mom1[i].data.shape)
    
    
    ycoords2 = (ycoords2-np.nanmean(ycoords2))*0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    xcoords2 = (xcoords2-np.nanmean(xcoords2))*-0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    ycoordsi = (ycoordsi-np.nanmean(ycoordsi))*0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    xcoordsi = (xcoordsi-np.nanmean(xcoordsi))*-0.125*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    vmod2 = np.where(~h2_msks[i],np.nan,h2_mom1_fitpars[i][2] + xcoords2*h2_dvdx[i] + ycoords2*h2_dvdy[i])
    vraw2 = np.where(~h2_msks[i],np.nan,h2_mom1[i])
    
    vmodi = np.where(~hi_msks[i],np.nan,hi_mom1_fitpars[i][2] + xcoordsi*hi_dvdx[i] + ycoordsi*hi_dvdy[i])
    vrawi = np.where(~hi_msks[i],np.nan,hi_mom1[i])
    
    
    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5),)

    rperp2 = (xcoords2*h2_dvdx[i] + ycoords2*h2_dvdy[i])/h2_omg[i]
    rperpi = (xcoordsi*hi_dvdx[i] + ycoordsi*hi_dvdy[i])/hi_omg[i]


    axs[0].scatter(np.ravel(rperp2),np.ravel(vraw2),s=np.log10(h2_mom1[i].size)/h2_ss[i],color=h2_color)
    axs[0].plot(np.ravel(rperp2),np.ravel(vmod2),color=dat_color,**ftline_kws)

    
    axs[1].scatter(np.ravel(rperpi),np.ravel(vrawi),s=np.log10(hi_mom1[i].size)/hi_ss[i],color=hi_color)
    axs[1].plot(np.ravel(rperpi),np.ravel(vmodi),color=dat_color,**ftline_kws)
    axs[1].plot(np.ravel(rperp2),np.ravel(vmod2),color=err_color,**ftline_kws,)


    fig.supxlabel('displacement from rotation axis [pc]',y=0.05)
    fig.supylabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
    fig.tight_layout()
    fig.suptitle(names[i],y=1.05)
    
#     plt.savefig(cat['name'][i]+'_rv',bbox_inches='tight')
    
    plt.show()


bigh0 = fits.open('../rn/bigh0.fits')[0].data
bigh1 = fits.open('../rn/bigh1.fits')[0].data


hi_dvdx, hi_dvdy = hi_mom1_fitpars.transpose()[:2]
h2_dvdx, h2_dvdy = h2_mom1_fitpars.transpose()[:2]


plt.plot([np.sqrt(np.abs(hi_omg[i].value/np.nanmax(hi_omg.value))) for i in range(22)],marker='.')
plt.plot([np.abs(hi_omg[i].value/np.nanmax(hi_omg.value)) for i in range(22)],marker='x')
plt.hlines([0,0.5,1],0,22)


plt.hist(hi_omg,bins=30,edgecolor='white')
plt.vlines([1e-2,5e-2,1e-1],0,3,color='red')
plt.show()


plt.plot(hi_omg,marker='.')
plt.hlines([0.05,0.01],0,22)
plt.yscale('log')


strength = [0.4,0.7,1]


st_color = [fns.segmented_spectral(0),
             fns.segmented_spectral(8),
             fns.segmented_spectral(16),
            fns.segmented_spectral(24),
           ]


plt.hist(hi_omg,edgecolor='white',bins=10,alpha=0.5)
# plt.vlines([1e-2,5e-2,1e-1],0,4,color='red')
plt.vlines(10**np.arange(-2.,0.,0.5),0,4,color='red')
plt.hist(h2_omg,edgecolor='white',bins=10,alpha=0.5)
# plt.vlines([6e-2,1e-1],0,4,color='red')

# plt.xscale('log')
plt.show()
# plt.show()


10**np.arange(-2.,0.,0.5)


hi_stal = np.zeros(22)
hi_stal[np.arange(22)[hi_omg.value >= 1e-1]] = 3
hi_stal[np.arange(22)[hi_omg.value < 1e-1]] = 2
hi_stal[np.arange(22)[hi_omg.value < 3e-2]] = 1
hi_stal[np.arange(22)[hi_omg.value < 1e-2]] = 0


hi_stal = hi_stal.astype('int')


h2_stal = np.zeros(22)
h2_stal[np.arange(22)[h2_omg.value >= 1e-1]] = 3
h2_stal[np.arange(22)[h2_omg.value < 1e-1]] = 2
h2_stal[np.arange(22)[h2_omg.value < 3e-2]] = 1
h2_stal[np.arange(22)[h2_omg.value < 1e-2]] = 0


h2_stal = h2_stal.astype('int')


plt.hist(h2_omg,edgecolor='white',bins=10)
plt.vlines([6e-2,1e-1],0,4,color='red')
plt.show()


plt.plot(h2_omg,marker='.')


st_color = [cm.YlGnBu_r.resampled(10)(1),
             cm.YlGnBu_r.resampled(10)(3),
             cm.YlGnBu_r.resampled(10)(6),
            cm.YlGnBu_r.resampled(10)(10),
           ]


cm.YlGnBu_r.resampled(10)


fns.segmented_spectral


help(plt.quiver)


plt.figure(dpi=300,figsize=(10,5))
# plt.contourf(bigh1,extent=[360,-360,-50,50],
#              levels=np.linspace(-40,40,27),
#              cmap='RdGy',extend='both',
#             antialiased=True)

plt.contourf(np.log10(bigh0),extent=[360,-360,-50,50],
             levels=np.linspace(20,22.5,20),
             cmap='binary_r',extend='both')


plt.xlim(90,-90)
plt.ylim(-40,50)

plt.gca().set_aspect(1)


[
    plt.gca().add_patch(Circle((l,b),r,fill=False,
                                color='black',linewidth=0.5,linestyle='--'))
    
    for l,b,r in zip(lcs,bcs,hi_r_deg)]


for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs,bcs,cat['distance<pc>'],hi_dvdx,hi_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
        
    plt.quiver(x,y,-dx_dg,dy_dg,pivot='mid',
               color=st_color[hi_stal[i]],
               zorder=10,)

    
for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs,bcs,cat['distance<pc>'],h2_dvdx,h2_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
        
    plt.quiver(x,y,-dx_dg,dy_dg,width=2.3e-3,pivot='mid',
               color='black',
               zorder=20,)
    
    
for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs,bcs,cat['distance<pc>'],h2_dvdx,h2_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
        
    plt.quiver(x,y,-dx_dg,dy_dg,width=2e-3,pivot='mid',
               color=st_color[h2_stal[i]],
               zorder=30,)
    

plt.scatter(lcs,bcs,s=300,marker='x',color='black',
            linewidth=0.2,zorder=15,alpha=0.7)  

    
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')

plt.xticks(np.arange(90,-120,-30),[str(tck)+'$\degree$' for tck in np.arange(90,-120,-30)])
plt.yticks(np.arange(-40,60,20),[str(tck)+'$\degree$' for tck in np.arange(-40,60,20)])

plt.savefig('catalog_omega_grad_galcent_TEST',bbox_inches='tight')

plt.show()


plt.contourf(np.log10(bigh0),extent=[360,-360,-50,50],
             levels=np.linspace(20,22.,15),
             cmap='magma',extend='both')
plt.xlim(90,-90)
plt.ylim(-40,40)


plt.figure(dpi=300,figsize=(10,5))
# plt.contourf(bigh1,extent=[360,-360,-50,50],
#              levels=np.linspace(-50,60,27),
#              cmap=fns.segmented_spectral,extend='both')

plt.contourf(np.log10(bigh0),extent=[360,-360,-50,50],
             levels=np.linspace(20,22.5,20),
             cmap='binary_r',extend='both')
# plt.colorbar()


plt.xlim(270,90)
plt.ylim(-40,40)

plt.gca().set_aspect(1)

for l,b,r in zip(lcs,bcs,hi_r_deg):
    plt.gca().add_patch(Circle((l,b),r,fill=False,color='black',
                               linewidth=0.3,linestyle='--',zorder=5) )
                        
    
for l,b,r in zip(lcs+360,bcs,hi_r_deg):
    plt.gca().add_patch(Circle((l,b),r,fill=False,
                               color='black',linewidth=0.3,linestyle='--',zorder=5) )   

    
for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs,bcs,cat['distance<pc>'],hi_dvdx,hi_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
    
        
    plt.quiver(x,y,-dx_dg,dy_dg,pivot='mid',
               color=st_color[hi_stal[i]],zorder=10)
    

for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs,bcs,cat['distance<pc>'],h2_dvdx,h2_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
    
    
    plt.quiver(x,y,-dx_dg,dy_dg,width=2.3e-3,pivot='mid',
               color='black',
              zorder=20)
    
for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs,bcs,cat['distance<pc>'],h2_dvdx,h2_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
    
    
    plt.quiver(x,y,-dx_dg,dy_dg,width=2e-3,pivot='mid',
               color=st_color[h2_stal[i]],
              zorder=20)
    

for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs+360,bcs,cat['distance<pc>'],hi_dvdx,hi_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
    
    
    plt.quiver(x,y,-dx_dg,dy_dg,pivot='mid',color=st_color[hi_stal[i]],zorder=10)
    

for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs+360,bcs,cat['distance<pc>'],h2_dvdx,h2_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
    
    
    plt.quiver(x,y,-dx_dg,dy_dg,width=2.3e-3,pivot='mid',
              zorder=20,color='black')
    
for i,x,y,d,dx_pc,dy_pc in zip(range(22),lcs+360,bcs,cat['distance<pc>'],h2_dvdx,h2_dvdy):
    dx_dg = dx_pc*d*u.deg.to(u.rad)
    dy_dg = dy_pc*d*u.deg.to(u.rad)
    
    
    plt.quiver(x,y,-dx_dg,dy_dg,width=2e-3,pivot='mid',
              zorder=20,color=st_color[h2_stal[i]])
    
# plt.scatter(lcs,bcs,s=1000,marker='+',color='white',linewidth=0.2,alpha=0.4)  
plt.scatter(lcs,bcs,s=300,marker='x',color='black',linewidth=0.2,
            alpha=0.7,zorder=15)  
plt.scatter(lcs+360,bcs,s=300,marker='x',color='black',linewidth=0.2,
            alpha=0.7,zorder=15)  

    
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')

plt.xticks(np.arange(270,60,-30),[str(tck)+'$\degree$' for tck in [-90,-120,-150,180,150,120,90]])
plt.yticks(np.arange(-40,60,20),[str(tck)+'$\degree$' for tck in np.arange(-40,60,20)])

plt.savefig('catalog_omega_grad_anticent_TEST',bbox_inches='tight')

plt.show()


xs = cat['distance<pc>']*np.sin(np.nanmean(cat['lon_GMC<deg>'],axis=1)*u.deg)
ys = cat['distance<pc>']*np.cos(np.nanmean(cat['lon_GMC<deg>'],axis=1)*u.deg)


dvdx_colors=['red','blue']


plt.hist(np.abs(hi_dvdy/hi_dvdx),bins=10**np.arange(-3,3,0.25),edgecolor='white',alpha=0.5)
plt.hist(np.abs(h2_dvdy/h2_dvdx),bins=10**np.arange(-3,3,0.25),edgecolor='white',alpha=0.5)

plt.xscale('log')
plt.show()


hi_dvdy[13]


hi_dvdx[13]


hi_mom1[8].data.shape


hi_mom1[13].data.shape


154/19


np.nanmean(hi_mom1[i])


hi_mom1_fitpars[8]


i = 8
vo,vi = bounds['vo'][i],bounds['vi'][i]
vels = np.arange(vo,vi)*1.3-100

spectra = np.zeros((19,vi-vo))
mu_1s = np.zeros(19)

for n in range(19):
    
    yo, yi = n*8,(n+1)*8
    xo, xi = 9*8, 10*8

    
    spectra[n] = np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    mu_1s[n] = fns.get_mom1(spectra[n],vels)
        

mu1_ys = [spectra[n][np.argmin(np.abs(mu_1s[n]-vels))] for n in range(19)]


mu_2s = [np.sqrt(fns.get_momN(spectra[n],vels,2)) for n in range(19)]


i=8


vels[1+np.argmax(spectra[9])]


6.9-5.3


normed_spectra = [spectrum/np.nanmax(spectrum) for spectrum in spectra]


spec_mom0_mass = np.nansum(spectra,axis=1)


opacities = spec_mom0_mass/np.nanmax(spec_mom0_mass)


plt.plot(opacities,marker='.')
plt.plot(opacities**1.7,marker='.')
# plt.plot(opacities**2,marker='.')


plt.yscale('log')
# plt.hlines([1,0.75,0.5,0.25],0,19)


hopla = 0.25


plt.figure(dpi=200,figsize=(5,9))
# plt.gca().set_aspect(15)
for n in range(19):
    yo, yi = n*8,(n+1)*8
    xo, xi = 9*8, 10*8
    
    plt.fill_between(vels,normed_spectra[n] + n*hopla,
                     n*hopla,
                     edgecolor='black',
                     color=cm.binary_r.resampled(19)(opacities[n]**1.2),
                     linewidth=0.,zorder=-n,
                    ),
    
    
    plt.plot(vels, n*hopla + normed_spectra[n],
             color=cm.Spectral_r.resampled(19)(n),
             linewidth=0.4,zorder=-n)
    
#     plt.plot(vels,
#              n*hopla+\
#              models.Gaussian1D(np.nanmax(spectra[n]),
#                                (n-9)*2*hi_omg[8].value+\
#                                vels[1+np.argmax(spectra[9])],
#                                3.5)(vels),
#             linewidth=1.7,color=cm.Spectral_r.resampled(19)(n),zorder=-n,
#             linestyle='--')
    
    
    plt.annotate('['+str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+', '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+']$\degree$',
                 (18,n*hopla+3e-2),fontsize='xx-small',
                    color=cm.Spectral_r.resampled(19)(n),)
    
    plt.scatter(mu_1s[n]-1.6,0.1,
               color=cm.Spectral_r.resampled(19)(n),
                s=15,alpha=0.5,linewidth=0.)
    
# for n in range(1,19):
#     plt.plot([mu_1s[n]-1.6,mu_1s[n-1]-1.6],[n*hopla,(n-1)*hopla],color='black',linewidth=0.5,zorder=-20) 
    
    
#     plt.plot([(n-9)*2*hi_omg[8].value+vels[1+np.argmax(spectra[9])],
#               (n-1-9)*2*hi_omg[8].value+vels[1+np.argmax(spectra[9])]],
#              [n*hopla,(n-1)*hopla],
#              color='black',linewidth=0.5)
    
# plt.scatter(vels[1+np.argmax(spectra[9])],
#                 np.nanmax(spectra[9]) + 9*hopla,
#                 marker='+',s=100,linewidth=3.,color='black')    

plt.xlim(-10.3,22.2)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')


plt.yticks([])
plt.ylim(-0.05,np.nanmax(normed_spectra[-1])+0.1+n*hopla)

plt.xticks(np.linspace(-10,20,5))
plt.ylabel(r'$T_\mathrm{B}$ [K]')

plt.savefig('hercules_stacked_TEST',bbox_inches='tight')


cm.binary


2.5e-1


mu1_ys = [74,
 62,
 50,
 44,
 44,
 44,
 42,
 30,
 32,
 42,
 45,
 40,
 46,
 40,
 32,
 22,
 20,
 23,
 23.]


mu_1s[9]


19*hi_omg[8].value


[(n-9)*2*hi_omg[8].value for n in range(19)]


mu_1s[9]


mu_2s[9]


models.Gaussian1D(mu1_ys[9],)


plt.plot(mu_1s,marker='.')
plt.hlines(6.88,0,19)
plt.plot(np.arange(19),np.arange(19)*-9.5*hi_omg[8].value+mu_1s[0]+mu_1s[9])


hi_sigmav[8]


np.nanmin(hi_mom2[8])


np.nanmin(mu_2s)


vels


ys = vels[np.argmin(np.abs(mu_1s[n]-vels)) for n in range(19)]


np.array((bounds['vo'][8],bounds['vi'][8]))*1.3-100


vels


plt.rcParams['font.size'] = 15


i = 13

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    yo, yi = n*13,(n+1)*13
    xo, xi = 5*13, 6*13

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    hopla = n*60
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
#     sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
#     h2sp = hopla + np.nanmean(h2h2[i].data[vo:vi,:,:],axis=(1,2))

        
    if thismu1 - np.nanmean(hi_mom1[i]) < 0:
        vo_clr = 0.5*(thismu1 - h1_lvls[i][0])/(np.nanmean(hi_mom1[i])-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
             (h1_lvls[i][1] - np.nanmean(hi_mom1[i]))
            
    
    plt.fill_between(vels,sp1,hopla,# color=fns.segmented_spectral(vo_clr),
            edgecolor='black',color=cm.binary.resampled(11)(n),
            linewidth=1.7,zorder=-n),
    
#     plt.scatter(thismu1,thisTB, color=fns.segmented_spectral(vo_clr),zorder=-n,),
    
    plt.plot(vels, sp1,
             color='black',linewidth=0.6,zorder=-n)
    
#     plt.plot(vels, h2sp,
#              color='black',linewidth=0.6)
    
    plt.vlines(fns.get_mom1(sp1,np.arange(vo,vi)*1.3-100),
               hopla,np.nanmax(sp1)+10,
               linewidth=0.5,zorder=-n,color='red')
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+20,
               linewidth=0.5,zorder=35,color='black')
    
       
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (15.5,hopla-9),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(3,17)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')


# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-15,np.nanmax(sp1)+15)


plt.savefig('mon_r2_stacked',bbox_inches='tight')


i = 8

plt.figure(dpi=200,figsize=(7,10))
for n in range(11):
    
    yo, yi = n*15,(n+1)*15
    xo, xi = 5*15, 6*15

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    hopla = n*40
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
#     sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
#     h2sp = hopla + np.nanmean(h2h2[i].data[vo:vi,:,:],axis=(1,2))

        
    if thismu1 - np.nanmean(hi_mom1[i]) < 0:
        vo_clr = 0.5*(thismu1 - h1_lvls[i][0])/(np.nanmean(hi_mom1[i])-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
             (h1_lvls[i][1] - np.nanmean(hi_mom1[i]))
            
    
    plt.fill_between(vels,sp1,hopla,# color=fns.segmented_spectral(vo_clr),
            edgecolor='black',color=cm.binary.resampled(11)(n),
            linewidth=1.7,zorder=-n),
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+20,
               linewidth=0.5,zorder=35,color='black')
    
    plt.vlines(fns.get_mom1(sp1,np.arange(vo,vi)*1.3-100),
               hopla,np.nanmax(sp1)+10,
               linewidth=0.5,zorder=-n,color='red')
    
    
#     plt.scatter(thismu1,thisTB, color=fns.segmented_spectral(vo_clr),zorder=-n,),
    
    plt.plot(vels, sp1,
             color='black',linewidth=0.6,zorder=-n)
    
#     plt.plot(vels, h2sp,
#              color='black',linewidth=0.6)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (-8,hopla+10),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-10,22)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')


# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,np.nanmax(sp1)+15)


plt.savefig('hercules_stacked',bbox_inches='tight')


np.nanmax(hi_omg.value*1e3)


plt.plot(h2_omg/hi_omg,marker='.')
plt.yscale('log')


hihi[16].data.shape


i = 16

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    yo, yi = n*20,(n+1)*20
    xo, xi = 5*20, 6*20

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    hopla = n*50
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
#     sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
#     h2sp = hopla + np.nanmean(h2h2[i].data[vo:vi,:,:],axis=(1,2))

        
    if thismu1 - np.nanmean(hi_mom1[i]) < 0:
        vo_clr = 0.5*(thismu1 - h1_lvls[i][0])/(np.nanmean(hi_mom1[i])-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
             (h1_lvls[i][1] - np.nanmean(hi_mom1[i]))
            
    
    plt.fill_between(vels,sp1,hopla,# color=fns.segmented_spectral(vo_clr),
            edgecolor='black',color=cm.binary.resampled(11)(n),
            linewidth=1.7,zorder=-n),
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+40,
               linewidth=0.5,zorder=35,color='black')
    
    plt.vlines(fns.get_mom1(sp1,np.arange(vo,vi)*1.3-100),
               hopla,np.nanmax(sp1)+10,
               linewidth=0.5,zorder=-n,color='red')
    
    
#     plt.scatter(thismu1,thisTB, color=fns.segmented_spectral(vo_clr),zorder=-n,),
    
    plt.plot(vels, sp1,
             color='black',linewidth=0.6,zorder=-n)
    
#     plt.plot(vels, h2sp,
#              color='black',linewidth=0.6)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (-0.5,hopla-10),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-1,15.5)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')


# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-15,np.nanmax(sp1)+15)


plt.savefig('orion_b_stacked',bbox_inches='tight')


hihi[20].data.shape


hihi[17].data.shape


i = 17

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    xo, xi = n*16,(n+1)*16
    yo, yi = 5*16, 6*16

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    hopla = n*50
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
#     sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
#     h2sp = hopla + np.nanmean(h2h2[i].data[vo:vi,:,:],axis=(1,2))

        
    if thismu1 - np.nanmean(hi_mom1[i]) < 0:
        vo_clr = 0.5*(thismu1 - h1_lvls[i][0])/(np.nanmean(hi_mom1[i])-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
             (h1_lvls[i][1] - np.nanmean(hi_mom1[i]))
            
    
    plt.fill_between(vels,sp1,hopla,# color=fns.segmented_spectral(vo_clr),
            edgecolor='black',color=cm.binary.resampled(11)(n),
            linewidth=1.7,zorder=-n),
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+40,
               linewidth=0.5,zorder=35,color='black')
    
    plt.vlines(thismu1,hopla,np.nanmax(sp1)+10,
               linewidth=0.5,zorder=-n,color='red')
    
    
#     plt.scatter(thismu1,thisTB, color=fns.segmented_spectral(vo_clr),zorder=-n,),
    
    plt.plot(vels, sp1,
             color='black',linewidth=0.6,zorder=-n)
    
#     plt.plot(vels, h2sp,
#              color='black',linewidth=0.6)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (-11,hopla+20),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-11.5,18)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')


# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,np.nanmax(sp1)+15)


plt.savefig('orion_l_stacked',bbox_inches='tight')


plt.plot(np.abs(hi_dvdy/hi_dvdx),marker='.')


plt.plot(np.abs(hi_dvdx/hi_dvdy),marker='.')


hi_theta = np.array([fns.dxdy_2_theta(hi_dvdx[i],hi_dvdy[i]) for i in range(22)])
h2_theta = np.array([fns.dxdy_2_theta(h2_dvdx[i],h2_dvdy[i]) for i in range(22)])

delphi = hi_theta-h2_theta

for i in range(22):
    if delphi[i] < -180:
        delphi[i] += 360
        
    if delphi[i] > 180:
        delphi[i] -= 360


np.array((cat['name'],np.around(hi_theta,2))).transpose()


np.array((cat['name'],np.around(h2_theta,2))).transpose()


np.cos(np.linspace(180,0,7)*u.deg).value


np.linspace(180,0,13)


fig,ax = plt.subplots(dpi=200)

bns = np.arange(-1,1.1,0.5/3)


ax.hist(np.cos(delphi*u.deg),bins=bns,
        color=[0,0,0,0.85],edgecolor='white',linewidth=5,)


ax.set_xlim(-1.05,1.05)
ax.set_ylim(0,6.5)


ax.vlines([-0.5,0.5,],0,6.5,color='red',linestyle='--')

ax.set_ylabel('N (clouds per bin)')
ax.set_xlabel(r'$\hat{J_\mathrm{H_2}} \cdot \hat{J_\mathrm{HI}}$')


ax.add_patch(Rectangle((-1.1,0),0.6,7,alpha=0.1,color='black'))
ax.add_patch(Rectangle((-1.1,0),1.6,7,alpha=0.1,color='black'))


ax2 = ax.twiny()
ax2.set_xticks([-1,-0.5,0,0.5,1],
               ['$180\degree$','$120\degree$','$90\degree$','$60\degree$','$0\degree$'])
ax2.set_xlim(-1.05,1.05)


plt.savefig('delphi',bbox_inches='tight')


plt.show()


hi_phi = hi_theta - 90
for i in range(22):
    if hi_phi[i] < -180:
        hi_phi[i] += 360
        
        
h2_phi = h2_theta - 90
for i in range(22):
    if h2_phi[i] < -180:
        h2_phi[i] += 360


hi_galphi = np.cos((hi_phi-90)*u.deg)*np.cos(bcs*u.deg)
h2_galphi = np.cos((h2_phi-90)*u.deg)*np.cos(bcs*u.deg)


fig,axs = plt.subplots(1,2,dpi=200,sharex=True,sharey=True,figsize=(10,5))

bns = np.arange(-1,1.1,0.5/3)


axs[0].hist(h2_galphi,bins=bns,
        color=h2_color,edgecolor='white',linewidth=5,)

axs[0].set_xlabel('$\hat{J_\mathrm{H_2}} \cdot \hat{J_\mathrm{MW}}$')


axs[1].hist(hi_galphi,bins=bns,
        color=hi_color,edgecolor='white',linewidth=5,)

axs[1].set_xlabel('$\hat{J_\mathrm{HI}} \cdot \hat{J_\mathrm{MW}}$')

for ax in axs:
    ax.set_xlim(-1.05,1.05)
    ax.set_ylim(0,6.5)
#     ax.vlines([-0.5,0.5,],0,6.5,color='red',linestyle='--')
    
    ax.add_patch(Rectangle((-1.1,0),0.6,7,alpha=0.1,color='black',linewidth=0.))
    ax.add_patch(Rectangle((-1.1,0),1.6,7,alpha=0.1,color='black',linewidth=0.))
    ax.vlines([-0.71,0.71],0,7,alpha=0.7,linewidth=0.5,zorder=-20,color=dat_color)


axs[0].set_ylabel('N (clouds per bin)',x=-0.1,)

# plt.xlabel('$\cos\ (\Delta \phi)$')


plt.subplots_adjust(wspace=0.15)


ax20 = axs[0].twiny()
ax20.set_xticks([-1,-0.5,0,0.5,1],
               ['$180\degree$','$120\degree$','$90\degree$','$60\degree$','$0\degree$'])
ax20.set_xlim(-1.05,1.05)
# ax20.set_xlabel('$\Delta \phi$',)


ax21 = axs[1].twiny()
ax21.set_xticks([-1,-0.5,0,0.5,1],
               ['$180\degree$','$120\degree$','$90\degree$','$60\degree$','$0\degree$'])
ax21.set_xlim(-1.05,1.05)
# ax21.set_xlabel('$\Delta \phi$',)


plt.savefig('delphi_gals',bbox_inches='tight')


plt.show()


np.cos(45*u.deg)


h2hi = fits.open('h2hi.fits')


from copy import deepcopy


h2_msks_hisz = deepcopy(hi_msks)
for i in range(22):
    h2_msks_hisz[i][bounds['boc'][i]-bounds['boh'][i]:bounds['bic'][i]-bounds['boh'][i],
                    bounds['loc'][i]-bounds['loh'][i]:bounds['lic'][i]-bounds['loh'][i]] = ~h2_msks[i] 

    h2_msks_hisz[i] = ~h2_msks_hisz[i]*hi_msks[i]


for i in range(22):
    
    plt.contourf(h2_msks_hisz[i])
#     plt.contourf(~h2_msks[i],extent=[bounds['loc'][i]-bounds['loh'][i],
#                                     bounds['lic'][i]-bounds['loh'][i],
#                                     bounds['boc'][i]-bounds['boh'][i],
#                                     bounds['bic'][i]-bounds['boh'][i]])
    plt.gca().set_aspect(1)
    plt.show()


h2_msks_hisz_cbs = [np.broadcast_to(h2_msks_hisz[i],
                    (155,h2_msks_hisz[i].shape[0],h2_msks_hisz[i].shape[1]))
                    for i in range(22)]

h2hi_mskd = [np.where(~h2_msks_hisz_cbs[i],np.nan,h2hi[i].data) for i in range(22)]


h2hi_mskd_0th = [np.nansum(h2hi_mskd[i],axis=0)*1.3*1.97e20 for i in range(22)]


plt.contourf(h2hi_mskd_0th[2])


#so say, the velocity gradients are p misaligned with the galaxy
#and this is the histogram showing that / what the break down is between eaech category
#the hi especially seems to favor the un/mis-aligneds, vs h2 is p uniform across the three
#so then we look at the patch by patch spectra
#and also the ratios of dvdy to dvdx
#and say okay so we should be aware that velocity gradients can be caused (at least for the hi)
#by blending rather than rotation
#the end


# plt.plot(np.log10(np.abs(h2_dvdy/h2_dvdx)),marker='.')
plt.plot(np.log10(np.abs(hi_dvdy/hi_dvdx)),marker='.')

plt.hlines([0,-1,-0.5,0.5,1],0,22)
plt.yticks([-1,0,1],[0.1,1,10])


plt.plot(np.log10(np.abs(h2_dvdy/h2_dvdx)),marker='.')
# plt.plot(np.log10(np.abs(hi_dvdy/hi_dvdx)),marker='.')

plt.hlines([0,-1,-0.5,0.5,1],0,22)
plt.yticks([-1,0,1],[0.1,1,10])


hihi[i].data.shape


h2h2[i].data.shape


loh


loc


lih,lic


for i in [8]:
    n = hihi[i].data.shape[1]
    
    num_ptch = 19
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    coords = np.array((inds[:-1],inds[1:])).transpose()
    slcs = np.array([(i,j) for j in coords for i in coords])
    
    
    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=300,figsize=(10,10))
    
    map_axs = np.flip(axs,axis=0).transpose().flat

        
    loh = bounds['loh'][i]*-0.125+360
    lih = bounds['lih'][i]*-0.125+360
    
    boh = bounds['boh'][i]*0.125-50
    bih = bounds['bih'][i]*0.125-50
    
    
    loc = bounds['loc'][i]*-0.125+360
    lic = bounds['lic'][i]*-0.125+360

    boc = bounds['boc'][i]*0.125-50
    bic = bounds['bic'][i]*0.125-50
    
    
    bkg_map = fig.add_subplot(111)#add_axes([0, 0, 1, 1])
    bkg_map.set_zorder(-10)
#     bkg_map.set_axis_off()
    bkg_map.contourf(hi_mom1[i],cmap=cm.Spectral.resampled(30),extend='both',#fns.segmented_spectral,extend='both',
                    extent=[loh,lih,boh,bih],
                    zorder=-20,alpha=0.6,antialiased=True,
                    levels=np.linspace(h1_lvls[i][0]-0.65,h1_lvls[i][1]))
    
    
#     bkg_map.contour(h2hi_mskd_0th[i],colors='black',#fns.segmented_spectral,extend='both',
#                     extent=[loh,lih,boh,bih],
#                     zorder=-20,alpha=0.6,
#                     levels=[20.5,21])
   
    
    bkg_map.set_xlim(loh,lih)
    bkg_map.set_ylim(boh,bih)
#     bkg_map.set_aspect('equal')
        
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(hihi[i].data[vos[i]:vis[i],bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(h2hi[i].data[vos[i]:vis[i],bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(vos[i],vis[i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.7,
                        linewidth=1,color='black')
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
#         map_axs[j].vlines(np.nanmedian(hi_mom1[i].data),
#                           0,1.5,color='black',
#                           linestyle='--',linewidth=1.5,alpha=0.2)
        
#         map_axs[j].scatter(np.nanmean(fns.get_mom1(hihi[i].data[vos[i]:vis[i],bo:bi,lo:li],
#                                    np.arange(vos[i],vis[i])*1.3-100)),0,
#                            s=10,color='black',alpha=0.5,linewidth=0)
        
        
#         if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
#         map_axs[j].plot(np.arange(vos[i],vis[i])*1.3-100,
#                         (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')
    
#         map_axs[j].plot(np.arange(vos[i],vis[i])*1.3-100,
#                         np.nanmean(h2h2[i].data[vos[i]:vis[i]],axis=(1,2))/np.nanmax(
#                             np.nanmean(h2h2[i].data[vos[i]:vis[i]],axis=(1,2))),
#                         linewidth=0.3,alpha=0.6,color='black')

#         map_axs[j].scatter(np.nanmean(fns.get_mom1(h2hi[i].data[vos[i]:vis[i],bo:bi,lo:li],
#                                np.arange(vos[i],vis[i])*1.3-100)),0,s=1,color='black')

#         map_axs[j].vlines(np.nanmedian(np.where(~h2_msks[i],np.nan,h2_mom1[i])),
#                               0,1.5,color='black',linestyle='--',linewidth=0.4,alpha=0.5)
        
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_ylim(-0.15,1.15)
        map_axs[j].set_axis_off()
    
    
    plt.subplots_adjust(wspace=0.07,hspace=0.07)
    plt.figaspect(1)
        
#     fig.suptitle(names[i],y=0.93)
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude',x=-0.01)
    
    
    bkg_map.set_xticks(np.around(np.linspace(loh,lih,5),1),
                       [str(tck)+'$\degree$' for tck in np.around(np.linspace(loh,lih,5),1)])
    bkg_map.set_yticks(np.around(np.linspace(boh,bih,5),1),
                       [str(tck)+'$\degree$' for tck in np.around(np.linspace(boh,bih,5),1)])
    
    
    plt.savefig(cat['name'][i]+'_pbp_TEST',bbox_inches='tight')
    
    
    plt.show()


hi_rad[8]


38/2


plt.rcParams['font.size'] = 25


for i in range(22):
    n = hihi[i].data.shape[1]
    
    num_ptch = 10
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    coords = np.array((inds[:-1],inds[1:])).transpose()
    slcs = np.array([(i,j) for j in coords for i in coords])
    
    
    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=300,figsize=(10,10))
    
    map_axs = np.flip(axs,axis=0).transpose().flat

        
    loh = bounds['loh'][i]*-0.125+360
    lih = bounds['lih'][i]*-0.125+360
    
    boh = bounds['boh'][i]*0.125-50
    bih = bounds['bih'][i]*0.125-50
    
    
    loc = bounds['loc'][i]*-0.125+360
    lic = bounds['lic'][i]*-0.125+360

    boc = bounds['boc'][i]*0.125-50
    bic = bounds['bic'][i]*0.125-50
    
    
    bkg_map = fig.add_subplot(111)#add_axes([0, 0, 1, 1])
    bkg_map.set_zorder(-10)
    bkg_map.set_axis_off()
    bkg_map.contourf(hi_mom1[i],cmap=fns.segmented_spectral,extend='both',
                    extent=[loh,lih,boh,bih],
                    zorder=-20,alpha=0.3,antialiased=True,
                    levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],26))
   
    
    bkg_map.set_xlim(loh,lih)
    bkg_map.set_ylim(boh,bih)
#     bkg_map.set_aspect('equal')
        
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(hihi[i].data[vos[i]:vis[i],bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(h2hi_mskd[i][vos[i]:vis[i],bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(vos[i],vis[i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.7,
                        linewidth=2,color='black')
        
        
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        
            
        map_axs[j].vlines(np.nanmedian(hi_mom1[i].data),
                          0,1.5,color='black',
                          linestyle='--',linewidth=1.5,alpha=0.2)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(hihi[i].data[vos[i]:vis[i],bo:bi,lo:li],
                                   np.arange(vos[i],vis[i])*1.3-100)),0,
                           s=10,color='black',alpha=0.5,linewidth=0)
        
        
#         if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
        map_axs[j].plot(np.arange(vos[i],vis[i])*1.3-100,
                        (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')

        map_axs[j].scatter(np.nanmean(fns.get_mom1(h2hi[i].data[vos[i]:vis[i],bo:bi,lo:li],
                               np.arange(vos[i],vis[i])*1.3-100)),0,s=1,color='black')

        map_axs[j].vlines(np.nanmedian(np.where(~h2_msks[i],np.nan,h2_mom1[i])),
                              0,1.5,color='black',linestyle='--',linewidth=0.4,alpha=0.5)
        
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_ylim(-0.15,1.15)
        map_axs[j].set_axis_off()
    
    
    plt.subplots_adjust(wspace=0.07,hspace=0.07)
    plt.figaspect(1)
        
    fig.suptitle(names[i],y=0.93)
#     plt.savefig(cat['name'][i]+'_pbp',bbox_inches='tight')

    plt.show()


#after these pick the clouds that i will do the stacked spectra for
#then do the zeroth moment maps
#then please just write oh my fucking god


np.cos(np.pi/4)


plt.plot(np.abs(hi_dvdx/hi_dvdy),marker='.')
plt.xticks(np.arange(0,22,3))
plt.yticks(np.arange(0,13))

plt.grid()
# plt.plot(hi_dvdy/hi_dvdx,marker='.')


np.cos(delphi*u.deg)[[15,16,17]]


parent_m = hi_mass + h2_mass
parent_J = hi_J + h2_J*np.cos(delphi*u.deg)

# parent_Etot = hi_Etot + h2_Etot #gonna ignore energy, bc 
#1) eskahun we didnt discuss it that much, and #2 there are some negative values (ie grav > the others)
#and i dont really want to deal with those sums rn


#okay so what we want to say instead is: here is the J of the parent, and we know the J of the molecular cloud
#if we take that delta j to be impulse, we can set it equal to force (or torque) times time
#ie if you integrate a torque over some time, you get change in angular momentum
#so i could do something thats more constrain-y than nail-down-y
#ie for a range of x torque values, we expect these ages
#or alternatively, for these ages guestimated from known/accepted gmc lifetimes, these torques can be ruled in/out


#the source of the torque could maybe be the change in density/viscosity between the two fluids?
#think of something that i can say is acting within the system (bc we still want to keep
#saying that J as a whole is conserved)
#as the gmc is starting to collapse, could some sort of accretion disk-y boundary layer be formed?
#bref this feels like the quickest set of logical leaps to get to a 1) satisfying-enough answer
#2) (more imp) in the next two days
#and if i cant think of a torque, (this is discussion, so can start writing before!)
#then thats okay ill just do the constraining thing


#and idt ill have time to redefine hi cubes & calculate velocity gradient mnamn, so move the "what if there 
#are more than two children" yemilewn to the discussion (or future work)


fig,axs = plt.subplots(1,2,dpi=200,sharex=True,sharey=True,figsize=(10,5))

bns = 10**np.arange(-2,2.,0.2)


axs[0].hist(np.abs(h2_dvdy/h2_dvdx),bins=bns,
        color=h2_color,edgecolor='white',linewidth=5,)

axs[0].set_xlabel(r'$[\Omega_y\ / \ \Omega_x]_\mathrm{H_2}$')


axs[1].hist(np.abs(hi_dvdy/hi_dvdx),bins=bns,
        color=hi_color,edgecolor='white',linewidth=5,)

axs[1].set_xlabel(r'$[\Omega_y\ / \ \Omega_x]_\mathrm{HI}$')

for ax in axs:
    ax.set_xscale('log')
    ax.set_xlim(0.05,30)
    ax.set_ylim(0,4.5)
    ax.vlines([1,],0,4.5,color='red',linestyle='--')
    ax.set_yticks(np.arange(5))
#     ax.add_patch(Rectangle((-1.1,0),0.6,7,alpha=0.1,color='black',linewidth=0.))
#     ax.add_patch(Rectangle((-1.1,0),1.6,7,alpha=0.1,color='black',linewidth=0.))


axs[0].set_ylabel('N (clouds per bin)',x=-0.1,)

# plt.xlabel('$\cos\ (\Delta \phi)$')


plt.subplots_adjust(wspace=0.15)


# plt.savefig('omg_rats',bbox_inches='tight')


plt.show()


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    
     axs.flat[i].contourf(np.where(~hi_msks[i],np.nan,hi_mom1[i]),cmap=fns.segmented_spectral,
                         levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],27),
                         extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                                 bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50,],
                         extend='both'),
    
    axs.flat[i].contour(np.where(~h2_msks[i],np.nan,np.log10(h2_mom0[i])),
                        linewidths=3e-1,colors=dat_color,levels = 7,
                                 extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                                         bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50,]),
    
    
    axs.flat[i].set_aspect(1),
    
    
    axs.flat[i].set_xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360),
    axs.flat[i].set_ylim(bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50),    

    
    axs.flat[i].tick_params(size=1.5),
    

) for i in range(12)]


[(
    axs.flat[i].set_title('['+str(h1_lvls[i][0])+', '+str(h1_lvls[i][1])+']',
                          
                          y=-0.02,x=0.13,color=dat_color,fontsize='small')
) for i in range(12)]


fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()

# fig.savefig('postage_himu1_1',bbox_inches='tight')
#


for i in [12,17,4,15,14,16]:
    h2_msks[i] = h2_msks[i]*np.log10(h2_mom0[i])>20.4


[(
    plt.hist(np.ravel(np.log10(h2_mom0[i])),bins=100),
    plt.title(i),
# plt.xlim(19,21),
plt.show())
    
    for i in [12,17,4,15,14,16,11]]


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    
     axs.flat[i-12].contourf(np.where(~hi_msks[i],np.nan,hi_mom1[i]),cmap=fns.segmented_spectral,
                         levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],27),
                         extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                                 bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50,],
                         extend='both'),
    
    axs.flat[i-12].contour(np.where(~h2_msks[i],np.nan,np.log10(h2_mom0[i])),
                        linewidths=3e-1,colors=dat_color,levels = 7,
                                 extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                                         bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50,]),
    
    
    axs.flat[i-12].set_aspect(1),
    
    
    axs.flat[i-12].set_xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360),
    axs.flat[i-12].set_ylim(bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50),    

    
    axs.flat[i-12].tick_params(size=1.5),
    

) for i in range(12,22)]


[(
    axs.flat[i-12].set_title('['+str(h1_lvls[i][0])+', '+str(h1_lvls[i][1])+']',
                          
                          y=-0.02,x=0.13,color=dat_color,fontsize='small')
) for i in range(12,22)]


for ax in axs.flat[-2:]:
    ax.set_axis_off()


fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()

# fig.savefig('postage_himu1_2',bbox_inches='tight')
#


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    axs.flat[i].contourf(np.where(~hi_msks[i],np.nan,np.log10(hi_mom0[i])),
                         levels=20,
                         extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                                 bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50,],
                         cmap='magma',extend='both',norm='asinh'),
    
    axs.flat[i].contour(np.where(~h2_msks[i],np.nan,np.log10(h2_mom0[i])),linewidths=7e-1,colors='black',
                        levels = 5,alpha=0.7,
                                 extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                                         bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50,]),
    
    
    axs.flat[i].invert_xaxis(),

    axs.flat[i].set_xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360),
    axs.flat[i].set_ylim(bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50)
)
for i in range(12)]


mnms = [np.around(np.log10(np.nanmin(hi_mom0[i])),1) for i in range(22)]
mxms = [np.around(np.log10(np.nanmax(hi_mom0[i])),1) for i in range(22)]

[(
    axs.flat[i].set_title('['+str(mnms[i])+', '+\
                             str(mxms[i])+']',
                          
                          y=-0.02,x=0.135,color=dat_color,fontsize='small')
) for i in range(12)]


# for ax in axs.flat[-2:]:
#     ax.set_axis_off()
    
fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()
# fig.savefig('postage_himu0_1',bbox_inches='tight')


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    axs.flat[i-12].contourf(np.where(~hi_msks[i],np.nan,np.log10(hi_mom0[i])),
                         levels=20,
                         extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                                 bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50,],
                         cmap='magma',extend='both',norm='asinh'),
    
    axs.flat[i-12].contour(np.where(~h2_msks[i],np.nan,np.log10(h2_mom0[i])),linewidths=7e-1,colors='black',
                        levels = 5,alpha=0.7,
                                 extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                                         bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50,]),
    
    
    axs.flat[i-12].invert_xaxis(),

    axs.flat[i-12].set_xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360),
    axs.flat[i-12].set_ylim(bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50)
)
for i in range(12,22)]


mnms = [np.around(np.log10(np.nanmin(hi_mom0[i])),1) for i in range(22)]
mxms = [np.around(np.log10(np.nanmax(hi_mom0[i])),1) for i in range(22)]

[(
    axs.flat[i-12].set_title('['+str(mnms[i])+', '+\
                             str(mxms[i])+']',
                          
                          y=-0.02,x=0.135,color=dat_color,fontsize='small')
) for i in range(12,22)]


for ax in axs.flat[-2:]:
    ax.set_axis_off()
    
fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()
# fig.savefig('postage_himu0_2',bbox_inches='tight')


plt.rcParams['font.size'] = 10


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    axs.flat[i].contourf(np.where(~hi_msks[i],np.nan,hi_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],27),
                extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                       bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50],
                 extend='both',alpha=0.3,antialiased=True,
                ),
    
    axs.flat[i].contourf(np.where(~h2_msks[i],np.nan,h2_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(c1_lvls[i][0],c1_lvls[i][1],27),
                extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                       bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50],
                 extend='both'
                ),
    
    
    axs.flat[i].invert_xaxis(),

    axs.flat[i].set_xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360),
    axs.flat[i].set_ylim(bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50)
)
for i in range(12)]


[(
    axs.flat[i].set_title('['+str(c1_lvls[i][0])+', '+\
                             str(c1_lvls[i][1])+']',
                          
                          y=-0.02,x=0.13,color=dat_color,fontsize='small')
) for i in range(12)]


# for ax in axs.flat[-2:]:
#     ax.set_axis_off()
    
fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()
# fig.savefig('postage_h2mu1_1',bbox_inches='tight')


fig,axs=plt.subplots(4,3,dpi=200,figsize=(8,10))#,sharey=True)

[(
    axs.flat[i-12].contourf(np.where(~hi_msks[i],np.nan,hi_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(h1_lvls[i][0],h1_lvls[i][1],27),
                extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                       bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50],
                 extend='both',alpha=0.3,antialiased=True,
                ),
    
    axs.flat[i-12].contourf(np.where(~h2_msks[i],np.nan,h2_mom1[i]),cmap=fns.segmented_spectral,
            levels=np.linspace(c1_lvls[i][0],c1_lvls[i][1],27),
                extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,
                       bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50],
                 extend='both'
                ),
    
    
    axs.flat[i-12].invert_xaxis(),

    axs.flat[i-12].set_xlim(bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360),
    axs.flat[i-12].set_ylim(bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50)
)
for i in range(12,22)]


[(
    axs.flat[i-12].set_title('['+str(c1_lvls[i][0])+', '+\
                             str(c1_lvls[i][1])+']',
                          
                          y=-0.02,x=0.13,color=dat_color,fontsize='small')
) for i in range(12,22)]


for ax in axs.flat[-2:]:
    ax.set_axis_off()
    
fig.supxlabel('Galactic Longitude $[\degree]$')
fig.supylabel('Galactic Latitude $[\degree]$')
    
fig.tight_layout()
# fig.savefig('postage_h2mu1_2',bbox_inches='tight')


hi_lvs = [np.nansum(hihi[i].data,axis=1)*0.125*1.82e18 for i in range(21)]
hi_bvs = [np.nansum(hihi[i].data,axis=2)*0.125*1.82e18 for i in range(21)]

h2_lvs = [np.nansum(h2h2[i].data,axis=1)*0.125*1.97e20 for i in range(21)]
h2_bvs = [np.nansum(h2h2[i].data,axis=2)*0.125*1.97e20 for i in range(21)]


for i in range(22):

    fig,axs = plt.subplots(1,3,dpi=200,figsize=(15,5),sharey=True)

    axs[0].contourf(np.log10(hi_lvs[i]),
                 levels=np.linspace(20,21.5,15),
                 cmap='magma',
                extent=[bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,-100,100])

    axs[0].contour(np.log10(h2_lvs[i]),levels=np.linspace(20,21.5,7),linewidths=0.5,
                 colors='black',
                  extent=[bounds['loc'][i]*-0.125+360,bounds['lic'][i]*-0.125+360,-100,100])

    axs[0].set_xlabel('Galactic Longitude')
    axs[0].set_ylabel('$v\mathrm{_{LSR}\ [km\ s^{-1}]}$')

    axs[0].hlines([bounds['vo'][i]*1.3-100,bounds['vi'][i]*1.3-100],
                 bounds['loh'][i]*-0.125+360,bounds['lih'][i]*-0.125+360,
                 linestyle='--',linewidth=0.7,color='black')


    axs[0].set_ylim(bounds['vo'][i]*1.3-100 - 20,bounds['vi'][i]*1.3-100 + 20)


    axs[1].contourf(np.log10(hi_bvs[i]),
                 levels=np.linspace(20,21.5,15),
                 cmap='magma',
                extent=[bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50,-100,100])

    axs[1].contour(np.log10(h2_bvs[i]),levels=np.linspace(20,21.5,7),linewidths=0.5,
                 colors='black',
                  extent=[bounds['boc'][i]*0.125-50,bounds['bic'][i]*0.125-50,-100,100])

    axs[1].set_xlabel('Galactic Latitude')
    axs[1].set_ylabel('$v\mathrm{_{LSR}\ [km\ s^{-1}]}$')


    axs[1].hlines([bounds['vo'][i]*1.3-100,bounds['vi'][i]*1.3-100],
                 bounds['boh'][i]*0.125-50,bounds['bih'][i]*0.125-50,
                 linestyle='--',linewidth=0.7,color='black')

    
    axs[2].plot(np.nanmean(hihi[i].data,axis=(1,2)),np.arange(155)*1.3-100,color=hi_color)
    axs[2].plot(50*np.nanmean(h2h2[i].data,axis=(1,2)),np.arange(155)*1.3-100,color=h2_color)
    
    axs[2].hlines([bounds['vo'][i]*1.3-100,bounds['vi'][i]*1.3-100],
                0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2)))+5,
                 linestyle='--',linewidth=0.7,color='black')
    
    axs[2].set_xlabel('$T\mathrm{_B\ [K]}$')
    axs[2].set_xlim(0,np.nanmax(np.nanmean(hihi[i].data,axis=(1,2)))+3)

    fig.tight_layout()
    
    plt.show()


hop = [25,25,25]
ptch_wdth = [int(hihi[i].data.shape[1]/11) for i in [8,13,16]]


np.log(1)


i = 17

plt.figure(dpi=200,figsize=(7,10))
for n in range(11):
    
    
    xo, xi = n*15,(n+1)*15
    yo, yi = 5*15, 6*15

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
#     sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
#     h2sp = np.nanmean(h2h2[i].data,axis=(1,2))[vo:vi]

        
    if thismu1 - np.nanmean(hi_mom1[i]) < 0:
        vo_clr = 0.5*(thismu1 - h1_lvls[i][0])/(np.nanmean(hi_mom1[i])-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
             (h1_lvls[i][1] - np.nanmean(hi_mom1[i]))
            
    
#     if thismu1 - np.nanmean(hi_mom1[i]) < 0:
#         vo_clr = 0.5*(thismu1 - np.nanmin(hi_mom1[i]))/\
#                  (np.nanmean(hi_mom1[i]) - np.nanmin(hi_mom1[i]))
    
#     else:
#         vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
#              (np.nanmax(hi_mom1[i]) - np.nanmean(hi_mom1[i]))
            

    plt.fill_between(vels,sp1,hopla,# color=fns.segmented_spectral(vo_clr),
            edgecolor='black',color=cm.binary.resampled(13)(n),
            linewidth=1.7,zorder=-n),
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+20,
               linewidth=0.5,zorder=-15,color='black')
    
    
    plt.scatter(thismu1,thisTB, color=fns.segmented_spectral(vo_clr),zorder=-n,),
    
#     plt.plot(vels, hopla + 10*h2sp/np.nanmax(h2sp),
#              color='black',linewidth=0.6,zorder=-n)
#     plt.plot(vels, hopla + sp2,
#              color='black',linewidth=0.6)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (-8,hopla+2),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-12,22)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,np.nanmax(sp1)+15)


# plt.savefig('lam_or2',bbox_inches='tight')


set up the simplest case: for each do a stacked spectra along omega hat, and along phi hat
    and then say, well the thing we *would* expect for solid(ish) body rotation is x, and the thing we see is y
    
start with just orion b and lambda orionis -> these two i will have talked about already, and also occupy
extremes-ish in terms of omx>omy, omy>omx
then show hercules and mr2 ? or tbh just mr2 bc its something that is within the field of the orions (will have
 already built familiarity with it, and it also demonstrates the line of sight distance thing)
the teaser can just be about all the things we can lear from dust, i dont have to give away hercules right now


#pick mr2,orionb,hercules as here are clouds that all have high dy/dx, 
#and also they represent high,middle,low omgmags (orb and ros are degenerate in being the bottom 1/3)
#the idea is we split the population into three : high, mid, low and put 7,8,7 clouds per bin


(hi_dvdy/hi_dvdx)[[8,13,16]]


[hi_mom1[i].shape[0] for i in [8,13,16]]
