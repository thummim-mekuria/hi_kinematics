#!/usr/bin/env python
# coding: utf-8


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


plt.rcParams['font.size'] = 5


phillips


phillips['type']


plt.figure(dpi=200)
plt.hist(phillips['type'],edgecolor='white')


# roughly split into the following categories [just for the aesthetiqueeee]
# - clumps (CL)
# - condensations (MIC,MSC
# - clouds (MI,MS)
# - complexes (MLC,ML)


np.unique(phillips['type'] )


clumps = phillips.query('type == "CL  "')
condensations = phillips.query('type == "MIC " | type == "MSC "')
clouds = phillips.query('type == "MI  " | type == "MS  "')
complexes = phillips.query('type == "ML  " | type == "MLC "')


condensations


curve_fit(line,np.log10(all_rs),np.log10(all_js))


np.sqrt(0.00172193)


curve_fit(line,np.log10(hi_rad.value),np.log10(hi_j.value))


curve_fit(line,np.log10(h2_rad.value),np.log10(h2_j.value))


curve_fit(line,np.log10(all_rs),np.log10(all_rs))


curve_fit(line,np.log10(phillips['R<pc>'].astype('float')),
          np.log10(phillips['j<cm2/s>'].astype('float')))


#okay the figure im making is:
#copy over the code from #143 for the jR all figr
#then change the exponent (+the label)
#make the four phillips datapoints same color diff shape
#say that when fit separately, phillips data has an index of 1.5 something
#hmm actually maybe just do a broken fit ie all of them with just yerasachewn indices


# plt.scatter(all_rs,all_js)
plt.figure(dpi=200)
plt.scatter(phillips['R<pc>'],phillips['j<cm2/s>'])
plt.scatter(hi_rad,hi_j)
plt.scatter(h2_rad,h2_j)

plt.plot(np.linspace(1e-2,1e3,10),10**22.97*np.linspace(1e-2,1e3,10)**1.5,color='red')
plt.plot(np.linspace(1e-2,1e3,10),10**22.97*np.linspace(1e-2,1e3,10)**1.24,
         color='red',linestyle='--')

plt.loglog()


plt.rcParams['font.size'] = 15.


plt.figure(dpi=200)
for df in [clumps,clouds,condensations,complexes]:
    plt.scatter(df['R<pc>'],df['j<cm2/s>'],label=df['type'].iloc[0])
    
plt.legend()
    
    
plt.loglog()


h2_color = fns.segmented_spectral(19)
hi_color = fns.segmented_spectral(24)

err_color = 'black'#fns.segmented_spectral(3)
err_kws = {'ms':3,'elinewidth':0.5,
           'ecolor':err_color,'capsize':1}

dat_color = 'black'

ftline_kws = {'linewidth':0.75,'linestyle':'--'}
ftfill_kws = {'linewidth':0.,'alpha':0.1}


phillips['R<pc>'].astype


np.log10(phillips['j<cm2/s>'])


curve_fit(line,np.log10(phillips['R<pc>'].astype('float')),
          np.log10(0.4*phillips['j<cm2/s>'].astype('float')))


curve_fit(line,np.log10(hi_rad.value),
          np.log10(0.4*hi_j.value))


np.sqrt(0.09710483)


curve_fit(line,np.log10(h2_rad.value),
          np.log10(0.4*h2_j.value))


np.sqrt(0.01953576)


p_typs = [clumps,clouds,condensations,complexes],
plbls = ['clumps','clouds','condensations','complexes']
psmbls = ['v','D','2','x']


plt.figure(dpi=200)

xs = np.linspace(1e-2,1e3,10)


plt.plot(core_rxs,core_rxs**1.53*10**23.06,**ftline_kws,color=dat_color)
plt.fill_between(core_rxs,core_rxs**(1.53+0.08)*10**23.06,
                 core_rxs**(1.53-0.08)*10**23.06,
                   **ftfill_kws,color=dat_color)


plt.plot(hi_rxs,hi_rxs**1.62*10**22.04,**ftline_kws,color=hi_color)
plt.fill_between(hi_rxs,hi_rxs**(1.62+0.31)*10**22.04,
                 hi_rxs**(1.62-0.31)*10**22.04,
                   **ftfill_kws,color=hi_color)

plt.plot(h2_rxs,h2_rxs**1.53*10**22.58,**ftline_kws,color=h2_color)
plt.fill_between(h2_rxs,h2_rxs**(1.53+0.14)*10**22.58,
                 h2_rxs**(1.53-0.14)*10**22.58,
                   **ftfill_kws,color=hi_color)


for df,plbl,psmbl in zip([clumps,clouds,condensations,complexes],plbls,psmbls):
    plt.scatter(df['R<pc>'],df['j<cm2/s>'],marker=psmbl,s=3,color=dat_color,
                label=plbl)
plt.errorbar(h2_rad,0.4*h2_j,xerr=h2_sigr,yerr=0.4*h2_sigj,fmt='o',**err_kws,
             label='molecular clouds',color=h2_color)

plt.errorbar(hi_rad,0.4*hi_j,xerr=hi_sigr,yerr=0.4*hi_sigj,fmt='s',
             **err_kws,label='atomic envelopes',color=hi_color)

# plt.scatter(h2_rad,h2_j,marker='.')
# plt.scatter(hi_rad,hi_j,marker='x')


plt.loglog()
plt.ylabel(r'$j\ \mathrm{[cm^2\ s^{-1}]}$')
plt.xlabel(r'$R\mathrm{\ [pc]}$')

plt.legend(markerscale=2,fontsize='small',bbox_to_anchor=(1, 0.77))
plt.xlim(3e-2,280)
plt.ylim(2e21,7e26)

# plt.savefig('jR_all',bbox_inches='tight')


# plt.figure(dpi=200)

fig,axs=plt.subplots(1,2,dpi=200,figsize=(10,5))

x = np.linspace(3e0,4e2,20)


axs[0].plot(x,x**1.53 * 10**(22.58),color=h2_color,
            linewidth=0.75,linestyle=':',label=r'$j\propto R^{1.53\pm0.14}$')

axs[0].plot(x,(x**1.62)*10**(22.58),color=hi_color,**ftline_kws)

axs[0].fill_between(x,x**(1.53-0.14) * 10**(22.58),
                      x**(1.53+0.14) * 10**(22.58),
                    color=h2_color,**ftfill_kws)

# axs[0].scatter(h2_rad,0.4*h2_j,color='black',zorder=20,s=10)
axs[0].errorbar(h2_rad,0.4*h2_j,xerr=h2_sigr,yerr=0.4*h2_sigj,fmt='o',**err_kws,color=dat_color)


axs[0].loglog()


axs[0].set_ylabel(r'$j\ \mathrm{[cm^2\ s^{-1}]}$')
axs[0].set_xlabel(r'$R\mathrm{_{H_2}\ [pc]}$')


axs[0].set_xlim(4e0,130)
axs[0].set_ylim(2e23,3e25)

# axs[0].legend()


axs[1].plot(x,(x**1.62)*10**(22.04),color=hi_color,
            linewidth=0.75,linestyle='--',label=r'$j\propto R^{1.62\pm0.31}$')

axs[1].plot(x,x**1.53 * 10**(22.04),color=h2_color,**ftline_kws)

axs[1].fill_between(x,x**(1.62-0.31) * 10**(22.04),
                      x**(1.62+0.31) * 10**(22.04),
                    color=hi_color,**ftfill_kws)


# axs[1].scatter(hi_rad,0.4*hi_j,color='black',zorder=20,s=10)
axs[1].errorbar(hi_rad,0.4*hi_j,xerr=hi_sigr,yerr=0.4*hi_sigj,fmt='s',**err_kws,color=dat_color)


axs[1].loglog()
# axs[1].legend()
axs[1].set_title(r'$j\propto R^{1.62\pm0.31}$',y=0.85,x=0.25)
axs[0].set_title(r'$j\propto R^{1.53\pm0.14}$',y=0.85,x=0.25)


axs[1].set_xlim(20,280)
axs[1].set_ylim(1e24,2e26)

axs[1].set_xlabel(r'$R\mathrm{_{HI}\ [pc]}$')


fig.tight_layout()

# plt.savefig('jR_hih2',bbox_inches='tight')


bigh0 = fits.open('../rn/bigh0.fits')[0].data


plt.figure(dpi=200)
plt.contourf(np.log10(bigh0),extent=[360,-360,-50,50],levels=15,cmap='magma')
plt.xlim(180,-180)
plt.gca().set_aspect(1)
plt.xlim(90,-90)


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
    [4.7,9.6],
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
    

plt.plot(hi_omg,marker='.')
plt.scatter([8,13,17],hi_omg[[8,13,17]],marker='x',color='red')


plt.xticks(np.arange(0,22,3),)
plt.grid()


plt.plot(h2_omg,marker='.')
plt.scatter([8,13,17],h2_omg[[8,13,17]],marker='x',color='red')


plt.xticks(np.arange(0,22,3),)
plt.grid()


plt.plot(hi_phi,marker='.')
plt.scatter([8,13,17],hi_phi[[8,13,17]],marker='x',color='red')

plt.xticks(np.arange(0,22,3),)
plt.grid()


# the [first chosen] one is hercules (wooo! it has the second highest hi gradient, highest h2 gradient)
#     then the next two are [monr2 and orion lam] > these three hi are at ~ 0,90,180 degrees from the galaxy
# 
# actually (maybe its a tale of two systems)?
# hercules & orion lam?
# 
# either way let me make the hercules spectra first, then gotta gotta GOTTA spend time writing eshi? please? 
# yesssss #teamwork
# 
# 


h2hi = fits.open('h2hi.fits')
from copy import deepcopy

h2_msks_hisz = deepcopy(hi_msks)
for i in range(22):
    h2_msks_hisz[i][bounds['boc'][i]-bounds['boh'][i]:bounds['bic'][i]-bounds['boh'][i],
                    bounds['loc'][i]-bounds['loh'][i]:bounds['lic'][i]-bounds['loh'][i]] = ~h2_msks[i] 

    h2_msks_hisz[i] = ~h2_msks_hisz[i]*hi_msks[i]
    
    
h2_msks_hisz_cbs = [np.broadcast_to(h2_msks_hisz[i],
                    (155,h2_msks_hisz[i].shape[0],h2_msks_hisz[i].shape[1]))
                    for i in range(22)]

h2hi_mskd = [np.where(~h2_msks_hisz_cbs[i],np.nan,h2hi[i].data) for i in range(22)]


np.arange(12)*14


plt.figure(dpi=200)
plt.contourf(np.nansum(h2hi_mskd[i],axis=0),levels=np.arange(0.5,11))
plt.gca().set_aspect(1)
plt.hlines(np.arange(12)*14,0,154)
plt.vlines(np.arange(12)*14,0,154)
plt.xlim(0,154)
plt.ylim(0,154)


#to play around with this somemore,
#it will make my life easier to precompute all the spectra, put them into an array
#then go from there, 
#but for now i will hardcode just a couple of the behavior that i want to talk about


i = 8

plt.figure(dpi=200,figsize=(7,10))
for n in range(11):
    
    
    hopla = n*25
    
    xo, xi =  5*14,6*14 
    #the rotation plane goes down the
    #middle of the cloud, which is the 6th "patch"
    #if you average out thee fov into 11x11 spectra
    yo, yi =  n*14,(n+1)*14

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
    sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
    h2sp = np.nanmean(h2h2[i].data,axis=(1,2))[vo:vi]

        
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
        
        
    plt.fill_between(vels,sp1,hopla, color=fns.segmented_spectral(vo_clr),
            edgecolor='black',
            linewidth=1.7,alpha=0.7,zorder=-n),
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+20,
               linewidth=0.5,zorder=-15,color='black')
    
    
    plt.scatter(thismu1,thisTB, color='black',zorder=-n,),
    
#     plt.plot(vels, hopla + 10*h2sp/np.nanmax(h2sp),
#              color='black',linewidth=0.6,zorder=-n)
#     plt.plot(vels, hopla + sp2,
#              color='black',linewidth=0.6)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (-8,hopla+2),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-10,22)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,np.nanmax(sp1)+15)


# plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
# plt.savefig('hercules_hi_psps',bbox_inches='tight')
# plt.savefig('herc_normed',bbox_inches='tight')


colors = [(0.5633986928104575, 0.0033986928104575, 0.0928104575163399, 1.0),
 (0.6633986928104575, 0.0633986928104575, 0.1928104575163399, 1.0),
 (0.7633986928104575, 0.1633986928104575, 0.2928104575163399, 1.0),
 (0.8758169934640523, 0.3045751633986928, 0.29411764705882354, 1.0),
 (0.9568627450980393, 0.42745098039215684, 0.2627450980392157, 1.0),
 (0.9803921568627452, 0.5973856209150326, 0.34117647058823525, 1.0),
 (0.9934640522875817, 0.7477124183006535, 0.4352941176470587, 1.0),
 (0.996078431372549, 0.8784313725490196, 0.5450980392156862, 1.0),
 (0.7477124183006538, 0.8980392156862746, 0.6274509803921569, 1.0),
 (0.5803921568627453, 0.8313725490196079, 0.6444444444444445, 1.0),
 (0.4, 0.7607843137254902, 0.6470588235294118, 1.0),
 (0.26405228758169935, 0.6091503267973857, 0.7098039215686275, 1.0),
 (0.25359477124183005, 0.45882352941176474, 0.7058823529411765, 1.0),
 (0.3686274509803922, 0.30980392156862746, 0.6352941176470588, 1.0),
 (0.2686274509803922, 0.20980392156862746, 0.5352941176470588, 1.0),
 (0.1686274509803922, 0.10980392156862746, 0.4352941176470588, 1.0)]


i = 8

plt.figure(dpi=200,figsize=(7,10))
for n in range(22):
    
    
    hopla = n*30
    
    xo, xi =  11*7,12*7 
    #the rotation plane goes down the
    #middle of the cloud, which is the 6th "patch"
    #if you average out thee fov into 11x11 spectra
    yo, yi =  n*7,(n+1)*7

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
    sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
    h2sp = np.nanmean(h2h2[i].data,axis=(1,2))[vo:vi]

        
#     if thismu1 - np.nanmean(hi_mom1[i]) < 0:
#         vo_clr = 0.5*(thismu1 - h1_lvls[i][0])/(np.nanmean(hi_mom1[i])-h1_lvls[i][0])
    
#     else:
#         vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
#              (h1_lvls[i][1] - np.nanmean(hi_mom1[i]))
        
        
    if thismu1 - np.nanmean(hi_mom1[i]) < 0:
        vo_clr = 0.5*(thismu1 - np.nanmin(hi_mom1[i]))/\
                 (np.nanmean(hi_mom1[i]) - np.nanmin(hi_mom1[i]))
    
    else:
        vo_clr = 0.5+ 0.5*(thismu1 - np.nanmean(hi_mom1[i]))/\
             (np.nanmax(hi_mom1[i]) - np.nanmean(hi_mom1[i]))
        
        
    plt.fill_between(vels,sp1,hopla, color=colors[int(10*vo_clr)],
            edgecolor='black',
            linewidth=1.7,zorder=-n),
    
    plt.vlines([np.nanmin(hi_mom1[i]),np.nanmedian(hi_mom1[i]),np.nanmax(hi_mom1[i])],
               hopla,np.nanmax(sp1)+20,
               linewidth=0.5,zorder=15,color='black')
    
    
    plt.scatter(thismu1,thisTB, color='black',zorder=-n,),
    
#     plt.plot(vels, hopla + 10*h2sp/np.nanmax(h2sp),
#              color='black',linewidth=0.6,zorder=-n)
#     plt.plot(vels, hopla + sp2,
#              color='black',linewidth=0.6)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+bounds['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+bounds['boh'][i])*0.125-50))+'$\degree$',
                 (-8,hopla+2),fontsize='xx-small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-10,22)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

# plt.vlines(np.nanmean(hi_mom1[i]),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,np.nanmax(sp1)+15)


# plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
# plt.savefig('hercules_hi_psps',bbox_inches='tight')
# plt.savefig('herc_normed',bbox_inches='tight')


fns.segmented_spectral


print(np.nanmin(hi_mom1[i]),
np.nanmean(hi_mom1[i]),
np.nanmedian(hi_mom1[i]),
np.nanmax(hi_mom1[i]))


plt.hist(np.ravel(hi_mom1[i]),bins=100)

plt.vlines((np.nanmin(hi_mom1[i]),
np.nanmean(hi_mom1[i]),
np.nanmedian(hi_mom1[i]),
np.nanmax(hi_mom1[i])),0,500,color='red')

plt.show()


plt.figure(dpi=200)
plt.plot(hi_phi,marker='.')
plt.scatter([8,9,13,15,16,17,20,21],hi_phi[[8,9,13,15,16,17,20,21]],color='red')

plt.xticks(np.arange(0,22,3))

plt.yticks([0,45,90,-135,135,-45,180,-180,-90])
plt.grid()


plt.figure(dpi=200)
plt.plot(hi_galphi,marker='.')
plt.scatter([8,9,13,15,16,17,20,21],hi_galphi[[8,9,13,15,16,17,20,21]],color='red')

plt.xticks(np.arange(0,22,3))

# plt.yticks([0,60,90,-60,-90,12a0,-120,180,-180])
plt.grid()


[(hihi[i].data.shape,hi_phi[i]) for i in [3,7,8,9,13,15,16,17,20,21]]


#qes bey, do this in pieces
#first, figure out where the main axis is / how wide in pixels, for 11 spectra, for each of the ^ above
#then makee the stacked spectra
#then choose two from each regime (+1:aligned, 0:unaligned, -1:misaligned)
#then describe the spectra


stacked_sps = [8,9,13,15,16,17,20,21] 


(hi_dvdy/hi_dvdx)[17]


plt.plot(hi_dvdy/hi_dvdx,marker='.')
plt.xticks(np.arange(0,22,3))
plt.ylim(-2,2)
# plt.ylim(-10,10)
plt.grid()


#from here ^, 9,15,21 are gonna be diagonals
#8,13,16,20 are gonna be vertical
#17 is the only one going horizontal


i = 17

plt.figure(dpi=200,figsize=(7,10))
for n in range(11):
    
    
    hopla = n*25
    
    
    xo, xi = n*15,(n+1)*15
    yo, yi = 5*15, 6*15

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
    sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
    h2sp = np.nanmean(h2h2[i].data,axis=(1,2))[vo:vi]

        
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
            

    plt.fill_between(vels,sp1,hopla, color=colors[int(10*vo_clr)],
            edgecolor='black',
            linewidth=1.7,alpha=0.7,zorder=-n),
    
    plt.vlines(np.nanmean(hi_mom1[i]),hopla,np.nanmax(sp1)+20,
               linewidth=0.5,zorder=-15,color='black')
    
    
    plt.scatter(thismu1,thisTB, color='black',zorder=-n,),
    
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


# plt.savefig('lam_or1',bbox_inches='tight')


cm.YlGnBu.resampled(11)


i = 17

# plt.figure(dpi=200,figsize=(7,10))
for n in range(11):
    
    
    hopla = n*25
    
    
    xo, xi = n*15,(n+1)*15
    yo, yi = 5*15, 6*15

    vo,vi = bounds['vo'][i],bounds['vi'][i]
    
    vels = np.arange(vo,vi)*1.3-100
    
    sp1 =   hopla + np.nanmean(hihi[i].data[vo:vi,yo:yi,xo:xi],axis=(1,2))
    
    thismu1 = fns.get_mom1(sp1,vels)
    thisclr = (thismu1 - h1_lvls[i][0])/(h1_lvls[i][-1]- h1_lvls[i][0])
    thisTB = sp1[np.argmin(np.abs(1 - vels/thismu1))]
    
    
    sp2 =   hopla + np.nanmean(h2hi_mskd[i][vo:vi,yo:yi,xo:xi],axis=(1,2))
    h2sp = np.nanmean(h2h2[i].data,axis=(1,2))[vo:vi]

        
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


plt.savefig('lam_or2',bbox_inches='tight')


plt.figure(dpi=200,figsize=(7,10))
for n in range(11):
    
    
    hopla = n*25
    
    xo, xi =  70,84 
    
    yo, yi =  n*14,(n+1)*14


plt.figure(dpi=200)
plt.plot(hi_omg,marker='.')
plt.plot(h2_omg,marker='.')

plt.xticks(np.arange(0,22,3))
plt.grid()


h2_mass
