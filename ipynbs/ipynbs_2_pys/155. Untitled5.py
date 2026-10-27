#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

from astropy.io import fits
import astropy.units as u
import astropy.constants as c

from scipy.optimize import curve_fit


import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


def line(x,m,b):
    return m*x + b


tab2 = pd.DataFrame(np.load('table2_2sighi.npy').transpose(),
             columns ='mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


# # comparisons wibnibb


def stats1(arr):
    return np.nanmin(arr),np.nanmax(arr)


def stats2(arr):
    return np.nanmean(arr),np.nanmedian(arr),np.nanstd(arr)


stats2(ibb['omgh2'])


ibb = pd.read_csv('ibb_tab1.csv',sep=',',nrows=33)


ibb_jhi = ((ibb['omghi'].to_numpy()*0.01*u.km/u.s/u.pc)*(70*u.pc)**2).cgs
ibb_jh2 = ((ibb['omgh2'].to_numpy()*0.01*u.km/u.s/u.pc)*(ibb['rh2'].to_numpy()*u.pc)**2).cgs


ibjh2 = (np.array((3.92,25.5,42.2,24.2,25.8))*u.pc*u.km/u.s).cgs
ibjhi = (np.array((25.2,120,107,60.3,73.4))*u.pc*u.km/u.s).cgs


ibjh2.value/tab2['jh2'].to_numpy()[[12,13,15,18,20]]


ibjhi.value/tab2['jhi'].to_numpy()[[12,13,15,18,20]]


stats2(ibb_jhi)


#one on one comparisons


#large-scale comparisons
#stil need to add in values from rosolowsky
val          min     max        av       med      std

mh2          808,    1.48e6,    1.40e5,  3.45e4,  3.07e5
mhi          2.28e4, 1.97e7,    2.17e6   6.38e5,  4.20e6
mhi/mh2      7.56,   33.98,     18.99,   17.35,   7.05

ibbmh2       2.9e4,  7.8e5,     1.58e5,  9.6e4,   1.49e5
ibbmhi       2e4,    5e5,       2.87e5,  2.8e5,   1.18e5
ibbmhi/ibbmh2

rh2          5.44,   141.3,     33.18,   26.65,   30.54
rhi          21.60,  264.3,     89.47,   75.44,   54.54
rhi/rh2      1.87,   5.06,      3.29,    3.24,    0.86

ibbrh2       #from rosolowsky table
ibbrhi       #all 70pc
ibbrhi/ibbrh2


omgh2        0.018,   0.25,     0.11,    0.087,   0.070
omghi        0.0073,  0.059,    0.030,   0.026,   0.015
omghi/omgh2  0.048    1.42,     0.42,    0.30,    0.32

ibbomgh2     0.069,    1.4,     0.74     0.79     0.35
ibbomhi      0.20,     1.31,    0.53     0.517    0.25
ibbohi/ibboh2
  
jh2          3.88e23, 3,21e26,  3.98e25, 1.32e25, 7.08e25
jhi          1.93e24, 4.21e26,  8.87e25, 5.56e25, 1.13e26
jhi/jh2      0.68,    10.9,     4.06,    3.51,    2.73     #maddalena and lambda-orionis

ibbjh2       #from rosolowsky table
ibbjhi       3.05e25, 1.98e26   8.07e25, 7.82e25, 3.81e25
ibbjhi/ibbjh2 #need jh2


ibrh2 = np.array([14.6,30.6,21.9,16.2,26.5,])
ibrhi = np.array([51.9,83.6,76.0,42.4,78.8])


iboh2 = np.array([4.6,6.68,22.,23.1,9.18])/100
ibohi = np.array([1.87,3.48,3.71,6.7,2.36])/100


ibjh2 = ((iboh2*u.km/u.s/u.pc)*(ibrh2*u.pc)**2).cgs
ibjhi = ((ibohi*u.km/u.s/u.pc)*(ibrhi*u.pc)**2).cgs


ibmh2 = np.array([2.5,3.0,2.4,1.4,4.0])*1e4*u.solMass


ib = pd.DataFrame(np.array((ibmh2,ibrh2,iboh2,ibjh2,ibrhi,ibohi,ibjhi)).transpose(),
                  columns='mh2,rh2,omgh2,jh2,rhi,omghi,jhi'.split(','))


in_ib = tab2.iloc[[12,13,15,18,20]]


np.full(80,len(ibb))


curve_fit(line,np.log10(ibb['rh2']),np.log10(1e-2*ibb['omgh2']))


curve_fit(line,np.log10(tab2['rh2']),np.log10(1e-3*tab2['omgh2']))


curve_fit(line,np.log10(tab2['rh2']),np.log10(tab2['jh2']))


curve_fit(line,np.log10(ibb['rh2']),np.log10(ibb['jh2']))


10**(-0.76)


plt.figure(dpi=200,figsize=(5,5))


plt.plot(np.linspace(2,200,5),0.29*np.linspace(2,200,5)**-0.41,color='black',
         linestyle='--',
         label='$\Omega_\mathrm{H_2} \propto R_\mathrm{H_2}^{-0.41}$')


plt.plot(np.linspace(2,200,5),0.17*np.linspace(2,200,5)**-0.34,color='green',
         linestyle='--',alpha=0.2,
         label='$\Omega_\mathrm{H_2} \propto R_\mathrm{H_2}^{-0.41}$')

plt.scatter(tab2['rh2'],1e-3*tab2['omgh2'],color='black',edgecolor='white',s=80,zorder=20)
# plt.scatter(tab2['rhi'],1e-3*tab2['omghi'],color=(0.8,0.,0.),edgecolor='white',s=80,zorder=20)
plt.scatter(ibb['rh2'],1e-2*ibb['omgh2'],color='green',alpha=0.2,edgecolor='white',s=80,zorder=20)

# plt.scatter(ibb['rh2'],np.full(len(ibb),70),color='green',edgecolor='white',s=80,zorder=-10,alpha=0.2)


plt.scatter(ib['rh2'],ib['omgh2'],color='grey',edgecolor='white',s=50,zorder=20)
[plt.plot([ib['rh2'].iloc[j],in_ib['rh2'].iloc[j]],
          [ib['omgh2'].iloc[j],1e-3*in_ib['omgh2'].iloc[j]],linewidth=0.75,zorder=21,label=j) 
 for j in range(5)]


plt.xlabel(r'$R_\mathrm{H_2}$ [pc]')
plt.ylabel(r'$\Omega_\mathrm{H_2}$ [pc]')
# plt.xlim(3,200)


plt.loglog()

plt.legend(fontsize='xx-small')
    
plt.savefig('Ro',bbox_inches='tight')


plt.figure(dpi=200,figsize=(5,5))


plt.plot(np.linspace(2,200,5),7.6*np.linspace(2,200,5)**0.72,color='black',
         linestyle='--',
         label='$R_\mathrm{HI} \propto R_\mathrm{H_2}^{0.72}$')

plt.scatter(tab2['rh2'],tab2['rhi'],color='black',edgecolor='white',s=80,zorder=20)

# plt.scatter(ibb['rh2'],np.full(len(ibb),70),color='green',edgecolor='white',s=80,zorder=-10,alpha=0.2)


plt.scatter(ib['rh2'],ib['rhi'],color='grey',edgecolor='white',s=50,zorder=20)
[plt.plot([ib['rh2'].iloc[j],in_ib['rh2'].iloc[j]],
          [ib['rhi'].iloc[j],in_ib['rhi'].iloc[j]],linewidth=0.75,zorder=21,label=j) for j in range(5)]


plt.xlabel(r'$R_\mathrm{H_2}$ [pc]')
plt.ylabel(r'$R_\mathrm{HI}$ [pc]')
# plt.xlim(3,200)


plt.loglog()

plt.legend(fontsize='small')
    
plt.savefig('RR',bbox_inches='tight')


ibb['jhi'] = ibb_jhi.value
ibb['jh2'] = ibb_jh2.value


ibb_omhi_lh2 = np.array([1.76,3.3,3.21,2.41,2.32,1.25,1.11,3.14,0.67,9.54,3.67,2.46,2.01,2.61,\
                           0.91,2.08,1.21,1.70,0.83,3.49,3.06,7.96,4.44,6.11,5.22,3.19,2.82,1.79,\
                           5.07,6.26,3.58,6.03,1.58,5.15,2.67,1.26,8.61,4.25,3.42,6.90,3.56,3.55])*0.01


ibb_mhi_lh2 = np.array([4.6,4.0,4.1,4.0,3.7,3.8,3.1,3.3,4.0,3.4,5.2,4.2,5.0,3.1,2.6,3.2,3.8,4.8,6.1,\
                         5.3,4.9,4.3,4.0,5.9,6.2,5.4,3.9,5.4,4.0,3.5,3.6,3.8,3.5,3.6,3.4,3.0,4.1,\
                         4.3,2.3,2.5,4.0,4.1])*0.01


ibb_jhi_lh2 = ((ibb_omhi_lh2*u.km/u.s/u.pc)*(70*u.pc)**2).cgs


plt.figure(dpi=200,figsize=(5,5))
bns=10**np.arange(23,27,1/7)
plt.hist(ibb['jh2'],bins=bns,color='green',histtype='step',zorder=20,linewidth=2)
plt.hist(tab2['jh2'],bins=bns,color='black',zorder=15)

plt.vlines(np.nanmean(ibb['jh2']),0,8.5,linestyle='--',color='green',alpha=0.5)
plt.vlines(np.nanmean(tab2['jh2']),0,8.5,linestyle='--',color='black')

plt.xscale('log')


# [plt.annotate(wch+' '+str(round(np.nanmean(jhi),3)),(2e24,4-0.25*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(2),['mw','m33',],[tab2['jh2'],ibb['jh2']],
#      ['black','green'])]


# [plt.annotate(wch+' '+str(round(np.nanstd(jhi),3)),(2e24,3-0.25*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(2),['mw','m33',],[tab2['jh2'],ibb['jh2']],
#      ['black','green'])]


# [plt.annotate(wch+' '+str(round(np.nanmedian(jhi),3)),(2e24,2-0.25*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(2),['mw','m33',],[tab2['jh2'],ibb['jh2']],
#      ['black','green'])]


plt.xlabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.ylabel('N (clouds)')
plt.xlim(2e23,5e26)
plt.ylim(0,5.5)
# plt.yticks(np.arange(5))

plt.savefig('jh2s',bbox_inches='tight')

plt.show()


plt.figure(dpi=200,figsize=(5,5))
bns=10**np.arange(24,27,1/7)
plt.hist(ibb_jhi_lh2.value,bins=bns,color='red',histtype='step',zorder=30,linewidth=2)
plt.hist(ibb['jhi'],bins=bns,color='green',histtype='step',zorder=30,linewidth=2,)
plt.hist(tab2['jhi'],bins=bns,color='black',zorder=25)

plt.vlines(np.nanmean(ibb_jhi_lh2.value),0,11,linestyle='--',color='red',alpha=0.5)
plt.vlines(np.nanmean(ibb['jhi']),0,11,linestyle='--',color='green',alpha=0.5)
plt.vlines(np.nanmean(tab2['jhi']),0,11,linestyle='--',color='black')


# [plt.annotate(wch+' '+str(round(np.nanmean(jhi),3)),(2e24,8-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],[tab2['jhi'],ibb['jhi'],ibb_jhi_lh2.value],
#      ['black','red','green'])]


# [plt.annotate(wch+' '+str(round(np.nanstd(jhi),3)),(2e24,5-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],[tab2['jhi'],ibb['jhi'],ibb_jhi_lh2.value],
#      ['black','red','green'])]


# [plt.annotate(wch+' '+str(round(np.nanmedian(jhi),3)),(2e24,2-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],[tab2['jhi'],ibb['jhi'],ibb_jhi_lh2.value],
#      ['black','red','green'])]

plt.xscale('log')

plt.xlabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.ylabel('N (clouds)')
plt.xlim(1e24,7e26)
plt.ylim(0,10.5)

plt.savefig('jhis',bbox_inches='tight')

plt.show()


curve_fit(line,np.log10(ibb_jh2.value),np.log10(ibb_jhi/ibb_jh2))


plt.figure(dpi=200,figsize=(5,5))

plt.plot(np.linspace(1e23,1e27,5),
            6e6*np.linspace(1e23,1e27,5)**-0.25,color='black',linestyle='--',
            label=r'$\frac{j_\mathrm{HI}}{j_\mathrm{H_2}} \propto j_\mathrm{H_2}^{-0.25}$')


plt.plot(np.linspace(1e23,1e27,5),
            6e27*np.linspace(1e23,1e27,5)**-1.07,color='green',linestyle='--',alpha=0.2,
            label=r'$\frac{j_\mathrm{HI}}{j_\mathrm{H_2}} \propto j_\mathrm{H_2}^{-1.07}$')


plt.scatter(tab2['jh2'],tab2['jhi']/tab2['jh2'],color='black',edgecolor='white',s=80)
plt.scatter(tab2['jh2'][badps],(tab2['jhi']/tab2['jh2'])[badps],edgecolor='black',color='white',s=80)


plt.scatter(ibb_jh2,ibb_jhi/ibb_jh2,color='green',edgecolor='white',s=80,zorder=-10,alpha=0.2)
plt.scatter(ibb_jh2,ibb_jhi/ibb_jh2,color='green',edgecolor='white',s=80,zorder=-10,alpha=0.2)


plt.scatter(ib['jh2'],ib['jhi']/ib['jh2'],color='grey',edgecolor='white',s=50,zorder=20)
[plt.plot([ib['jh2'].iloc[j],in_ib['jh2'].iloc[j]],
          [(ib['jhi']/ib['jh2']).iloc[j],(in_ib['jhi']/in_ib['jh2']).iloc[j]],
          linewidth=0.75,zorder=21,label=j) for j in range(5)]


plt.xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel(r'$j_\mathrm{HI}\ /\ j_\mathrm{H_2}$')

# plt.yticks([1,10],labels=['1','10'])
plt.xlim(3e23,5e26)
plt.ylim(5e-1,5e2)

plt.loglog()
plt.legend(fontsize='small')
plt.savefig('jjj',bbox_inches='tight')


# # tdiss


sigh2,sighi = np.load('meansigisg.npy').transpose()[1:]


def get_tdis(dat):
    egrav_hi = (c.G * dat['mhi'].to_numpy()*u.solMass / (dat['rhi'].to_numpy()*u.pc)).cgs
    erot_hi = (0.4*dat['jhi'].to_numpy()*u.cm**2/u.s*dat['omghi'].to_numpy()*u.km/u.s/u.pc).cgs
    ekin_hi = (0.5*(sighi*u.km/u.s)**2).cgs

    egrav_h2 = (c.G * dat['mh2'].to_numpy()*u.solMass / (dat['rh2'].to_numpy()*u.pc)).cgs
    erot_h2 = (0.4*dat['jh2'].to_numpy()*u.cm**2/u.s*dat['omgh2'].to_numpy()*u.km/u.s/u.pc).cgs
    ekin_h2 = (0.5*(sigh2*u.km/u.s)**2).cgs
    
    etot_h2 = erot_h2 + egrav_h2  #+ ekin_h2
    etot_hi = erot_hi + egrav_hi  #+ ekin_hi
    
    tdis = (dat['jhi']/etot_hi.value)*u.s.to(u.Myr)
    
    jinit = (dat['mh2']*dat['jh2'] + dat['mhi']*dat['jhi'])/(dat['mh2'] + dat['mhi'])

    
    return tdis,jinit


ibb['rhi']=np.full(len(ibb),70)


omgh2_old = np.array([87.6,51.0,65.5,40.7,94.1,111.2,65.8,19.3,323.1,176.8,67.0,\
                      50.6,107.77,50.6,79.2,70.6,44.5,160.6,159.7,63.8,87.6,120.3])
omghi_old = np.array([80.8,58.1,21.9,12.0,15.2,2.2,31.8,7.0,68.5,15.6,15.8,19.9,\
                      7.8,17.7,29.6,4.8,9.5,22.2,18.8,22.8,10.5,62.3])


rh2_old = np.array([5.45,38.66,5.70,32.43,8.92,4.74,19.72,98.36,7.25,13.34,16.32,\
                   57.88,24.66,29.72,15.04,24.40,17.91,25.99,15.50,17.15,30.71,10.16])
rhi_old = np.array([22.25, 123.67, 21.75,110.41,46.57,28.5,69.67,245.53,38.89,43.87,\
                   64.65,153.29,104.21,112.08,62.39,108.40,94.22,71.79,76.86,63.53,109.73,57.72])


jh2_old = ((omgh2_old*u.m/u.s/u.pc)*(rh2_old*u.pc)**2).cgs
jhi_old = ((omghi_old*u.m/u.s/u.pc)*(rhi_old*u.pc)**2).cgs


mhi_old = np.array([8.67,984.9,7.43,543.64,40.47,10.80,208.68,5105.2,67.14,48.53,\
                    139.80,1774.90,524.75,520.58,144.76,281.99,379.94,\
                    322.44,114.95,56.39,699.24,131.71])*1e3


mh2_old = np.array([0.88,151.78,0.83,108.02,8.11,1.86,27.14,11187.8,\
                    4.72,6.78,21.68,289.06,90.82,112.99,19.49,102.22,\
                   67.13,29.69,36.43,20.57,106.03,15.43])*1e3


get_tdis(tab2)


np.nanmin(tab2['omghi'])


plt.figure(dpi=200,figsize=(5,5))
bns = 10**np.arange(23,27,0.15)
plt.hist(oldtab2['jhi'],edgecolor='black',bins=bns,histtype='step',label='old')
# plt.hist(oldtab2['mh2'],color='green',bins=bns)

plt.hist(tab2['jhi'],histtype='step',edgecolor='red',
         bins=bns,label='new',linestyle='--')
# plt.hist(tab2['mh2'],histtype='step',edgecolor='green',bins=bns)

plt.xscale('log')
plt.legend()


plt.figure(dpi=200,figsize=(5,5))
bns = 10**np.arange(-1,3.5,0.15)
plt.hist(oldtab2['omghi'],edgecolor='black',bins=bns,histtype='step',label='old')
# plt.hist(oldtab2['mh2'],color='green',bins=bns)

plt.hist(tab2['omghi'],histtype='step',edgecolor='red',
         bins=bns,label='new',linestyle='--')
# plt.hist(tab2['mh2'],histtype='step',edgecolor='green',bins=bns)

plt.xscale('log')
plt.legend()


plt.figure(dpi=200,figsize=(5,5))
bns = 10**np.arange(3,7.5,0.25)
plt.hist(oldtab2['mhi'],edgecolor='black',bins=bns,histtype='step',label='old')
# plt.hist(oldtab2['mh2'],color='green',bins=bns)

plt.hist(tab2['mhi'],histtype='step',edgecolor='red',
         bins=bns,label='new',linestyle='--')
# plt.hist(tab2['mh2'],histtype='step',edgecolor='green',bins=bns)

plt.xscale('log')
plt.legend()


oldtab2 = pd.DataFrame(np.array([mh2_old,rh2_old,omgh2_old,jh2_old,
                                 mhi_old,rhi_old,omghi_old,jhi_old]).transpose(),
                       columns='mh2,rh2,omgh2,jh2,mhi,rhi,omghi,jhi'.split(','))


1/7


plt.figure(dpi=200,figsize=(5,5))
plt.hist(get_tdis(tab2)[1],color='black',zorder=20,
          bins=10**np.arange(24,27,1/7)) 

plt.hist(get_tdis(ibb)[1],color='green',linewidth=2,zorder=25,
          bins=10**np.arange(24,27,1/7),histtype='step') 


plt.hist(ibb_jhi_lh2.value,color='red',linewidth=2,zorder=25,
          bins=10**np.arange(24,27,1/7),histtype='step') 


plt.vlines(np.nanmean(get_tdis(ibb)[1]),0,12,linestyle='--',color='green',alpha=0.5)
plt.vlines(np.nanmean(ibb_jhi_lh2.value),0,12,linestyle='--',color='red',alpha=0.5)
plt.vlines(np.nanmean(get_tdis(tab2)[1]),0,12,linestyle='--',color='black')


[plt.annotate(wch+' '+str(round(np.nanmean(jhi),3)),(2e24,8-0.5*i),fontsize='small',color=clr) 
 for i,wch,jhi,clr in 
 zip(range(3),['mw','m33','m33woh2'],
     [get_tdis(tab2)[1],get_tdis(ibb)[1],ibb_jhi_lh2.value],
     ['black','green','red'])]


[plt.annotate(wch+' '+str(round(np.nanstd(jhi),3)),(2e24,5-0.5*i),fontsize='small',color=clr) 
 for i,wch,jhi,clr in 
 zip(range(3),['mw','m33','m33woh2'],
     [get_tdis(tab2)[1],get_tdis(ibb)[1],ibb_jhi_lh2.value],
     ['black','green','red'])]


[plt.annotate(wch+' '+str(round(np.nanmedian(jhi),3)),(2e24,2-0.5*i),fontsize='small',color=clr) 
 for i,wch,jhi,clr in 
 zip(range(3),['mw','m33','m33woh2'],
     [get_tdis(tab2)[1],get_tdis(ibb)[1],ibb_jhi_lh2.value],
     ['black','green','red'])]


plt.xlabel('pre-GMC $j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel('N(clouds)')
plt.xscale('log')
plt.ylim(0,11.5)
# plt.xlim()
# plt.loglog()
# plt.savefig('jinit',bbox_inches='tight')


# # scaling relation plots


p99_tab1 = pd.read_csv('../J_A+AS_134_241_table1.dat.gz.txt',skiprows=9,delimiter='|')
p99_tab1.columns = ['Name','Note','Ah Am RAs Ed Em DEs','l_Dist','Dist','u_Dist','RV','DV DV2',\
                   'l_Size','Size','n_Size','Size2','Line','L_OMEGA','OMEGA_OMEGA2','l_GAMMA',\
                   'GAMMA','l_Mass','Mass_Mass2','MassEm','THETA_THETA2','OMEGAo','u_OMEGAo','Type','r_']


#Q: which tracers are used for the data in phillips
p99_tab1['Line'].to_numpy()


pOmgs = np.array([ float(n.split()[0]) 
                  for i, n in enumerate(p99_tab1['OMEGA_OMEGA2'][1:-1])]) *\
                    1e-14*(1/u.s).to(u.km/u.s/u.pc).value

pRs = p99_tab1['Size'][1:-1].astype('float').to_numpy()

pMsraw = [pM.split() for pM in p99_tab1['Mass_Mass2'][1:-1].to_numpy()]

pTyps = p99_tab1['Type'][1:-1].to_numpy()


pMs = np.zeros(len(pMsraw))
for i,n in enumerate(pMsraw):
    if len(n)==0: # ""
        pMs[i] = 0
        
    elif len(n)==1: # "num"
        pMs[i] = float(n[0])
        
    elif len(n)==3: #"num","->","num"
        pMs[i] = np.nanmean(np.array(n)[[0,2]].astype('float'))
        
    elif n[0][-1]=='>':
        pMs[i] = np.nanmean(np.array((n[0][:6],n[1])).astype('float'))
        
    elif n[1][:2]=='->':
        pMs[i] = np.nanmean(np.array((n[0],n[1][2:])).astype('float'))
        
    else:
        pMs[i] = np.nanmean(np.array(pMsraw[126]).astype('float'))
        

p99_dat = pd.DataFrame(np.array([pTyps,pMs,pRs,pOmgs]).transpose(),
                      columns=['type','M<Msun>','R<pc>','Omega<km/s/pc>'])

p99_dat = p99_dat[~np.isinf(np.log10(pMs))] #M > 0
p99_dat = p99_dat[~(p99_dat['type']=='D/R ')] #no disks/rings
p99_dat = p99_dat[~np.array([('F' in nm) for nm in p99_dat['type']])] #no filaments

p99_dat['j<cm2/s>'] = ((p99_dat['R<pc>'].to_numpy()*u.pc)**2 *\
                        p99_dat['Omega<km/s/pc>'].to_numpy()*u.km/u.s/u.pc).cgs.value.astype('float')


chenx_rs = np.array([1703,1232,3530,1440,1272,1763,1329,2762])*u.AU.to(u.pc)
chenx_js = (np.array([0.12,0.12,0.45,0.13,0.24,0.12,0.42,0.45])*\
            1e-3*u.km/u.s*u.pc).cgs.value


chenh_rs = np.array([0.038,0.033,0.025,0.032,0.026,0.026,
                     0.034,0.054,0.042,0.020,0.022,0.024,0.19,0.028,0.045,
                    0.044,0.050,0.048])

chenh_js = (np.array([1.3,1.3,0.36,1.3,0.8,0.52,0.59,2.5,1.,0.8,1.5,\
 1.6,6.,0.37,0.25,0.34,0.46,0.54])*1e-3 * u.km/u.s*u.pc).cgs


curve_fit(line,np.log10(tab2['rh2']),np.log10(tab2['jh2']))


np.sqrt(curve_fit(line,np.log10(tab2['rh2']),np.log10(tab2['jh2']))[1][0,0])


np.log10(2e23)


curve_fit(line,np.log10(ibb['rh2'].to_numpy()),np.log10(ibb['jh2'].to_numpy()))


np.sqrt(curve_fit(line,np.log10(ibb['rh2'].to_numpy()),np.log10(ibb['jh2'].to_numpy()))[1][0,0])


plt.figure(dpi=200,figsize=(5,5))


plt.scatter(tab2['rh2'],tab2['jh2'],color=(0.,0.8,0.),
            s=100,zorder=30,edgecolor='white')


plt.scatter(ibb['rh2'],ibb['jh2'],color='black',alpha=0.5,
            s=50,edgecolor='white',zorder=-20)
plt.plot(xs, 5e22*xs**1.66,color='black',linestyle='--',linewidth=0.75,zorder=-30,alpha=0.5)


plt.xlim(4,200)
plt.ylim(3e23,4e26)
plt.loglog()

xs = np.linspace(1,200,10)
plt.plot(xs, 1e23*xs**1.59,color='black',linestyle='--')
plt.fill_between(xs,1e23*xs**(1.59+0.18),1e23*xs**(1.59-0.18),
                 alpha=0.2,color='black',linewidth=0.)


plt.xlabel(r'$R_\mathrm{H_2}$ [pc]')
plt.ylabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')

plt.annotate('$j_\mathrm{H_2} \propto R^{1.59 \pm 0.18}$',(4.5,2e26))
plt.annotate('$j_\mathrm{H_2} \propto R^{1.66 \pm 0.19}$',(4.5,8e25),alpha=0.2)


plt.savefig('jR_h2',bbox_inches='tight')


curve_fit(line,np.log10(tab2['rhi']),np.log10(tab2['jhi']))


np.sqrt(curve_fit(line,np.log10(tab2['rhi']),np.log10(tab2['jhi']))[1][0,0])


np.log10(6e22)


plt.figure(dpi=200,figsize=(5,5))


plt.scatter(tab2['rhi'],tab2['jhi'],color=(0.8,0.,0.),
            s=100,edgecolor='white',zorder=30)

# [plt.vlines(tab2['rhi'][i],jhi1[i]-jhi[i],jhi3[i]-jhi[i]) for i in range(22)]


plt.xlim(20,300)
plt.ylim(15e23,6e26)
plt.loglog()

xs = np.linspace(10,300,10)
plt.plot(xs, 1e22*xs**1.88,color='black',linestyle='--')
plt.fill_between(xs,1e22*xs**(1.88+0.19),1e22*xs**(1.88-0.19),
                 alpha=0.2,color='black',linewidth=0.,zorder=10)


plt.annotate('$j \propto R^{1.88 \pm 0.19}$',(22,4e26))

plt.vlines(70,1e23,7e26,zorder=-10,linestyle=':',color='black',alpha=0.5)

plt.xlabel(r'$R_\mathrm{HI}\ [\mathrm{pc}]$')
plt.ylabel(r'$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('jR_hi',bbox_inches='tight')


plt.figure(dpi=200,figsize=(5,5))
plt.hist(np.arccos(np.cos(np.mod((tab2['ththi']-tab2['thth2'])*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg),
        bins=np.arange(0,361,10),edgecolor='white',color='black')

plt.vlines(np.nanmean(np.arccos(np.cos(np.mod(ibb['deltht']*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg)),
          0,5.5,color=(0.8,0.,0.),linestyle='--',linewidth=0.75,zorder=-10)

plt.vlines(np.nanmean(np.arccos(np.cos(np.mod((tab2['ththi']-tab2['thth2'])*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg)),
          0,5.5,color='black',linestyle='--')
# plt.hist(np.arccos(np.cos(np.mod(ibb['deltht']*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg),
#         bins=np.arange(0,361,10),edgecolor='white',color=(0.,0.8,0.),alpha=0.3)

# plt.hist(-360+np.mod((old_thhi-old_thh2)*u.deg.to(u.rad),2*np.pi)*u.rad.to(u.deg),
#         bins=np.arange(-360,361,10),edgecolor='white',color='red')

plt.xticks(np.arange(-360,361,30),labels = [str(num)+r'$\degree$' for num in np.arange(-360,361,30)],
           )

plt.xlim(-10,190)
plt.ylim(0,5.5)
plt.ylabel('N (clouds)')
# plt.grid()
plt.xlabel(r'$| \theta_\mathrm{HI} - \theta_\mathrm{H_2} |$')


plt.savefig('theta_h2_hi',bbox_inches='tight')

plt.show()
# plt.xticks


cat_ed = np.array((
['AqS',  41,     28,    -20,   -12,   -5.1,  9.2, 133],
['Cal',  174,    154,   -14,   -1,    -25, 18.3, 470],
['Cam',  150,    144,    16,    19,   -8,  10.0,  213],
['CMj', -132.5, -138.5, -4,     1.5,   9.2,  25, 1209],
['Chm', -54,    -67,    -20,   -12,   -3.8,  9.2,  183],
['CrA',  4,       -1,    -23.5, -16,    0.1,  20, 151],
['Xbn',  -138,   -142,   -12.5, -6.5,   -10,  30, 886],
['Gem',  -164,   -174,   -4,     6,    -6.4,  30, 1786],
['Her',  49,     41,     7,     11.5, -10.3, 23.5, 227],
['Lac',  98,     94,    -17,   -7,    -10,  14.4, 503],
['Lup', -10,    -27,     5,     25.5, -7.7,  17.0, 189],
['Mad', -140,   -147,   -5,    -0.5,   20, 39.1, 2072],
['MOB', -155,   -162,   -0.5,   3.5,  -2.5,  15.7, 745],
['MR2',  -141,   -150,  -15,   -10,    0.,  18.3, 778],
['Oph',  12,     -12,    3,     28,   -6.4,  13.1, 144],
['OrA', -140,   -156,   -21,   -16.5,  -10,  25.4, 432],
['OrB', -149,   -158,   -18,   -9,    -1.2,  17,   423],
['LOr', -158.5, -173,   -19,   -6,    -11.6, 19.6, 402],
['Per',  162,    154,   -27.5, -16,   -10,  20.5, 294],
['Pol',  128,    117,    20,    30,   -12.9, 6.6,  352],
['Ros', -151,   -155,   -3.5,  -0.5,   1.4,  30, 1304],
['Tau',  177,    163,   -21,   -10,   -5.1,  17.0, 141],))


coordses = cat_ed[:,1:5].astype('float')


bc = np.nanmean(coordses[:,2:],axis=1)


galhi = np.arccos(np.cos(tab2['ththi']*u.deg)*np.cos(bc*u.deg))*u.rad.to(u.deg)


galh2 = np.arccos(np.cos(tab2['thth2']*u.deg)*np.cos(bc*u.deg))*u.rad.to(u.deg)


plt.figure(dpi=200,figsize=(6,6))
plt.hist(np.mod(np.array(galh2)*u.deg.to(u.rad),2*np.pi)*u.rad.to(u.deg),
        bins=np.arange(0,361,10),edgecolor='white',color=(0.0,0.8,0.))


plt.xticks(np.arange(-360,361,30),labels = [str(num)+r'$\degree$' for num in np.arange(-360,361,30)])

plt.xlim(-10,190)
plt.ylim(0,5.5)
plt.ylabel('N (clouds)')
# plt.grid()
plt.xlabel(r'$| \theta_\mathrm{MW} - \theta_\mathrm{H_2} |$')

plt.savefig('theta_mw_h2',bbox_inches='tight')

plt.show()
# plt.xticks


plt.figure(dpi=200,figsize=(6,6))
plt.hist(np.mod(np.array(galhi)*u.deg.to(u.rad),2*np.pi)*u.rad.to(u.deg),
        bins=np.arange(0,361,10),edgecolor='white',color=(0.8,0.,0.))

plt.xticks(np.arange(-360,361,30),labels = [str(num)+r'$\degree$' for num in np.arange(-360,361,30)])

plt.xlim(-10,190)
plt.ylim(0,5.5)
plt.ylabel('N (clouds)')
# plt.grid()
plt.xlabel(r'$| \theta_\mathrm{MW} - \theta_\mathrm{HI} |$')

# plt.savefig('theta_mw_hi',bbox_inches='tight')
plt.show()
# plt.xticks


goodman = pd.read_csv('../csvs/goodman_1993_data.csv')


tatmtsu_rs = np.array([
    0.049,0.088,0.083,0.084,0.059,0.099,0.091,0.049,0.098,0.095,0.052,\
    0.104,0.094,0.083,0.111,0.105,0.074,0.057,0.072,0.043,0.062,0.074,0.115,\
    0.083,0.084,0.062,0.074,0.081,0.062,0.052,0.083,0.157,0.057,0.072
])


tatmtsu_js = np.array([
    7.4e20,1.7e21,7.6e20,7.2e20,0,6e20,1.2e21,1.6e21,9.3e21,3.2e21,3e20,\
    2.3e21,1.8e21,2.2e21,2.1e21,1.9e21,2.2e21,0,0,8.4e20,1.5e21,0,1.9e21,4.4e21,\
    0,7.7e20,0,1.1e21,1.2e21,1.1e21,3e21,3.5e21,8.9e20,0
])


tatmtsu_rs = tatmtsu_rs[~np.isinf(np.log10(tatmtsu_js))]
tatmtsu_js = tatmtsu_js[~np.isinf(np.log10(tatmtsu_js))]


goodman


all_rs = np.concatenate((tab2['rh2'].to_numpy(),tab2['rhi'].to_numpy(),\
         p99_dat['R<pc>'].to_numpy().astype('float'),chenh_rs,chenx_rs,\
                         tatmtsu_rs,ibb['rh2'],np.full(len(ibb),70),np.full(len(ibb_jhi_lh2),70)))

all_js = np.concatenate((tab2['jh2'].to_numpy(),tab2['jhi'].to_numpy(),\
         p99_dat['j<cm2/s>'].to_numpy().astype('float'),chenh_js.value,chenx_js,\
                         tatmtsu_js,ibb['jh2'].to_numpy(),ibb['jhi'].to_numpy(),ibb_jhi_lh2.value))


curve_fit(line,np.log10(all_rs),np.log10(all_js))


np.sqrt(curve_fit(line,np.log10(all_rs),np.log10(all_js))[1][0,0])


plt.figure(dpi=200,figsize=(5,5))
plt.scatter(chenh_rs,chenh_js,s=70,
            color='black',marker='s',edgecolor='white',label='chenh')
plt.scatter(chenx_rs,chenx_js,s=70,
            color='black',marker='<',edgecolor='white',label='chenx')
plt.scatter(tatmtsu_rs,tatmtsu_js,s=70,
            color='black',marker='o',edgecolor='white',label='ttmts')
plt.scatter(p99_dat['R<pc>'],p99_dat['j<cm2/s>'],label='p99',
            color='black',marker='*',s=150,edgecolor='white')

plt.scatter(tab2['rh2'],tab2['jh2'],
            color=(0.,0.8,0.),s=75,edgecolor='white',label='h2')
plt.scatter(tab2['rhi'],tab2['jhi'],
            color=(0.8,0.,0.),s=75,edgecolor='white',label='hi')


plt.scatter(ibb['rh2'],ibb['jh2'],color='black',alpha=0.5,
            s=70,edgecolor='green',zorder=-20,label='ibb h2')

plt.scatter(np.full(len(ibb),70),ibb['jhi'],color='black',alpha=0.5,
            s=70,edgecolor='red',zorder=30,label='ibb hi')

plt.scatter(np.full(len(ibb_jhi_lh2),70),ibb_jhi_lh2,color='black',alpha=0.5,
            s=70,edgecolor='orange',zorder=20,label='ibb hi wo h2')
 
# plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],
#             color=(0.6,0.6,0.6),marker='^',s=30,)


plt.xlim(3e-3,4e2)
plt.ylim(2e19,8e26)
plt.loglog()


xs = np.linspace(1e-3,5e2,10)
plt.plot(xs, 1e23*xs**1.48,color='black',linestyle='--',zorder=-10)
plt.fill_between(xs,1e23*xs**(1.48+0.02),1e23*xs**(1.48-0.02),
                 alpha=0.2,color='black',linewidth=0.,zorder=-10)


plt.annotate('$j \propto R^{1.48 \pm 0.02}$',(5e-3,2e26))

plt.legend(fontsize='xx-small')


plt.xlabel('R [pc]')
plt.ylabel('j [$\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('jR_all',bbox_inches='tight')


ds = cat_ed[:,7].astype('float')


lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')


hi4pi, dht = fits.open('../fitses/both_same_grid_720d.fits')


rhi_px = (np.arctan(tab2['rhi']/ds)*u.rad.to(u.deg)*8).astype('int')


names = cat_ed[:,0]


names


h2maps, himaps = np.ravel([fits.open('../fitses/'+names[i]+'maps.fits') for i in range(22)]).reshape(22,2).transpose()


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den

def rvsv_pearson(i, show=False): 
    
    ccoords = coordses[i]
    d = ds[i]
    omgs = tab2[['omgh2','omghi']].to_numpy()[i]/1e3
    thts = tab2[['thth2','ththi']].to_numpy()[i]
    
    l0c,l1c,b0c,b1c = ccoords
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2s, his = h2maps[i], himaps[i]
    h2_0,h2_1,h2_2 = h2s.data
    hi_0,hi_1,hi_2 = his.data
    
    rhi_px = hi_0.shape[0]/2

    l0h,l1h = lc + rhi_px*0.125, lc - rhi_px*0.125
    b0h,b1h = bc - rhi_px*0.125, bc + rhi_px*0.125

    hlons,hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))
    clons,clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))
    
    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*d
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*d

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*d
    cys_pc = (clats-bc)*u.deg.to(u.rad)*d
    
    
    h2_dvdss, hi_dvdss = [(om*np.cos(th*u.deg).value, om*np.sin(th*u.deg).value) for om,th in zip(omgs,thts)]
    
    hbperp = (hi_dvdss[0]*hxs_pc + hi_dvdss[1]*hys_pc)/omgs[1]
    cbperp = (h2_dvdss[0]*cxs_pc + h2_dvdss[1]*cys_pc)/omgs[0]


    hv_model =  hi_vos[i] + hi_dvdss[0]*hxs_pc + hi_dvdss[1]*hys_pc
    cv_model =  h2_vos[i] + h2_dvdss[0]*cxs_pc + hi_dvdss[1]*cys_pc
    
#     print(hbperp.shape)
#     plt.contourf(cbperp,cmap='PiYG')
#     plt.contour(cv_model)
#     plt.contourf(cv_model)
#     plt.contourf(h2_1,cmap='RdGy')
#     plt.gca().set_aspect(1)
#     plt.show()
    
    
#     plt.contourf(hv_model,cmap="RdGy")
#     plt.gca().set_aspect(1)
#     plt.colorbar()
#     plt.show()
    
#     plt.contourf(hi_1,cmap="PiYG")
#     plt.gca().set_aspect(1)
#     plt.colorbar()
#     plt.show()

    
    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hi_2),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod',])

#     hdf_noisy = hdf
#     hdf = hdf_noisy[hdf_noisy['w']>hi_Nlims[i][0]]
    
    
    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0),
            np.ravel(h2_1),
            np.ravel(h2_2),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod'])
    

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()
    
#     print(cdf)


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    h2p, hip = pearson(meancdf['rper'], meancdf['v']),pearson(meanhdf['rper'], meanhdf['v'])

#     print(meancdf['v'],meancdf['vmod'])
    
    if show:
    
        fig,axs = plt.subplots(1,2,figsize=(12,5))
#         axs[0].scatter(cdf['rper'],cdf['v'])
#         axs[1].scatter(hdf['rper'],hdf['v'])
        axs[0].fill_between(meancdf['rper'],
                            meancdf['v']+2*sigcdf['v'],
                            meancdf['v']-2*sigcdf['v'],
                            color=(0.8,0.8,0.8,1),linewidth=0.5)

        axs[0].fill_between(meancdf['rper'],
                            meancdf['v']+sigcdf['v'],
                            meancdf['v']-sigcdf['v'],
                            color=(0.5,0.5,0.5,0.8),linewidth=0.5)


        axs[0].plot(meancdf['rper'],meancdf['vmod'],color='red',linewidth=2,linestyle='--')
        axs[0].plot(meancdf['rper'],meancdf['v'],color='black',marker='o',linewidth=0.5,ms=2)


        axs[1].fill_between(meanhdf['rper'],
                            meanhdf['v']+2*sighdf['v'],
                            meanhdf['v']-2*sighdf['v'],
                            color=(0.8,0.8,0.8,1),linewidth=0.5)

        axs[1].fill_between(meanhdf['rper'],
                            meanhdf['v']+sighdf['v'],
                            meanhdf['v']-sighdf['v'],
                            color=(0.5,0.5,0.5,0.8),linewidth=0.5)


        axs[1].plot(meanhdf['rper'],meanhdf['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[1].plot(meanhdf['rper'],meanhdf['vmod'],color='red',linewidth=2,linestyle='--')

        axs[0].annotate(str(np.around(h2p,2)),(meancdf['rper'].iloc[0],
                                               np.nanmax(meancdf['v'])+np.nanmean(sigcdf['v'])))

        axs[1].annotate(str(np.around(hip,2)),(meanhdf['rper'].iloc[0],
                                               np.nanmax(meanhdf['v'])+np.nanmean(sighdf['v'])))


        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')

        axs[0].set_ylabel(r'velocity [$\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names[i])

        fig.savefig('v_vs_r_'+names[i][:3],bbox_inches='tight')


        plt.show()
    
    return h2p,hip


fitpars = np.load('plane_fitpars.npy')/1e3


h2_vos, hi_vos = fitpars.swapaxes(0,1)[:,:,2]


[rvsv_pearson(i,show=True) for i in range(22)]


h2ps,hips = np.array([rvsv_pearson(i) for i in range(22)]).transpose()


badh2p = [5,10,11]
badhip = [2,5,]


plt.scatter(h2ps,hips)
[plt.annotate(names[i],(h2ps[i],hips[i])) for i in range(22)]
plt.hlines(0.5,-1,1)
plt.vlines(0.5,1,-1)
plt.show()


plt.hist(h2ps,edgecolor='white',color=(0.,0.8,0.),bins=np.arange(-1,1.1,0.1))

plt.show()


plt.hist(hips,edgecolor='white',color=(0.8,0.,0.),bins=np.arange(-1,1.1,0.1))
plt.show()


plt.figure(dpi=200,figsize=(7,7))


plt.scatter(tab2['rh2'],tab2['jh2'],color=(0.,0.8,0.),
            s=100,zorder=30,edgecolor='white')


plt.scatter(tab2['rh2'][badh2p],tab2['jh2'][badh2p],edgecolor=(0.,0.8,0.),
            s=100,zorder=30,color='white')

plt.xlim(4,200)
plt.ylim(3e23,4e26)
plt.loglog()

xs = np.linspace(1,200,10)
plt.plot(xs, 1e23*xs**1.59,color='black',linestyle='--')
plt.fill_between(xs,1e23*xs**(1.59+0.18),1e23*xs**(1.59-0.18),
                 alpha=0.2,color='black',linewidth=0.)


plt.xlabel(r'$R_\mathrm{H_2}$ [pc]')
plt.ylabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')

plt.annotate('$j_\mathrm{H_2} \propto R^{1.59 \pm 0.18}$',(4.5,2e26))

plt.savefig('jR_h2',bbox_inches='tight')


plt.figure(dpi=200,figsize=(7,7))


plt.scatter(tab2['rhi'],tab2['jhi'],color=(0.8,0.,0.),
            s=100,edgecolor='white',zorder=30)

plt.scatter(tab2['rhi'][badhip],tab2['jhi'][badhip],edgecolor=(0.8,0.,0.),
            s=100,color='white',zorder=30)

# [plt.vlines(tab2['rhi'][i],jhi1[i]-jhi[i],jhi3[i]-jhi[i]) for i in range(22)]


plt.xlim(20,300)
plt.ylim(15e23,6e26)
plt.loglog()

xs = np.linspace(10,300,10)
plt.plot(xs, 1e22*xs**1.88,color='black',linestyle='--')
plt.fill_between(xs,1e22*xs**(1.88+0.19),1e22*xs**(1.88-0.19),
                 alpha=0.2,color='black',linewidth=0.)


plt.annotate('$j_\mathrm{HI} \propto R^{1.88 \pm 0.19}$',(22,4e26))


plt.xlabel(r'$R_\mathrm{HI}\ [\mathrm{pc}]$')
plt.ylabel(r'$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('jR_hi',bbox_inches='tight')


curve_fit(line,np.log10(tab2['jh2']),np.log10(tab2['jhi']/tab2['jh2']))


curve_fit(line,np.log10(tab2['rh2']),np.log10(tab2['rhi']))


plt.scatter(tab2['rh2'],tab2['rhi'])
plt.loglog()


plt.scatter(tab2['jh2'],tab2['jhi']*(tab2['rh2']/tab2['rhi']))
plt.plot(tab2['jh2'],tab2['jh2'])
plt.loglog()

plt.xlabel('jh2')
plt.ylabel('jhi rescaled')


plt.scatter(tab2['jh2'],tab2['jhi'])
plt.loglog()

plt.xlabel('jh2')
plt.ylabel('jhi')
plt.title('jhi propto jh2 ^ 0.75')


10**(0.88)


10**(6.8)


fig,axs = plt.subplots(1,2,dpi=200,figsize=(11,5))


axs[0].plot(np.linspace(2,200,5),
            7.6*np.linspace(2,200,5)**0.72,color='red',label='$R_\mathrm{HI} \propto R_\mathrm{H_2}^{0.72}$')
axs[0].scatter(tab2['rh2'],tab2['rhi'],color='black',edgecolor='white',s=80,zorder=20)


axs[0].set_xlabel(r'$R_\mathrm{H_2}$ [pc]')
axs[0].set_ylabel(r'$R_\mathrm{HI}$ [pc]')
axs[0].set_xlim(3,200)


axs[1].plot(np.linspace(1e23,1e27,5),
            6e6*np.linspace(1e23,1e27,5)**-0.25,color='red',
            label=r'$\frac{j_\mathrm{HI}}{j_\mathrm{H_2}} \propto j_\mathrm{H_2}^{-0.25}$')


axs[1].scatter(tab2['jh2'],tab2['jhi']/tab2['jh2'],color='black',edgecolor='white',s=80)
axs[1].scatter(tab2['jh2'][badps],(tab2['jhi']/tab2['jh2'])[badps],edgecolor='black',color='white',s=80)


axs[1].set_xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
axs[1].set_ylabel(r'$j_\mathrm{HI}\ /\ j_\mathrm{H_2}$')

axs[1].set_yticks([1,10],labels=['1','10'])
axs[1].set_xlim(3e23,5e26)

for ax in axs:
    ax.loglog()
    ax.legend()

    
plt.savefig('RRjj',bbox_inches='tight')


oldrhis = [22.25,123.67,21.74,110.41,46.57,28.5,69.67,245.53,38.89,43.87,64.65,\
          153.29,104.21,112.08,62.39,108.40,94.22,71.79,76.86,63.53,109.73,57.72]
