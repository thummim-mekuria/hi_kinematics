#!/usr/bin/env python
# coding: utf-8


# %%time

import datacube as dc
import clouds as cld

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram
from astropy.visualization import ZScaleInterval as zscale


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf

DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 


#----------------------
cat = pd.read_csv('Mekuria_catalog.csv',sep=',',skiprows=5)
names = cat['name']


l_unit = u.deg
b_unit = u.deg
d_unit = u.pc


#----------------------


def read_DHT(ns='all'):
    output={}
    if ns == 'all':
        ns = ['03','06','09','12','16','19','21','26','27','37']
        #for now, figure out DHT01 (bad data) 
        #& DHT18 (two ceph ds, beeg) 
        #& DHT08 (beeg) 
        #& DHT24, DHT34 (lam or & chamaeleon are disjointed down to the really low sig levels)
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output


#---------------------
surveys = read_DHT()
for s in surveys:
    k = np.where(cat['survey']==s)[0]

    these = cat['name'][k]
    
    nums = np.zeros(len(these))
    for i, this in enumerate(these):
        nums[i] = np.where(cat['name']==this)[0]
        
    nums = nums.astype(np.int64)
    cube = surveys[s]
    surveys[s] = [cube, nums]


for s in surveys:
    clob = surveys[s][0]
    sig = np.nanstd(clob.cube)

    clouds = surveys[s][1]
    
    if (s not in ['DHT03','DHT06', 'DHT34','DHT37']) and (s not in ['DHT21']) and (s != 'DHT08'):
        d = Dendrogram.compute(clob.cube, min_value=3*sig,min_delta=sig, min_npix=400)
        
    elif s in ['DHT03','DHT06', 'DHT34','DHT37']:
        d = Dendrogram.compute(clob.cube, min_value=sig,min_delta=sig, min_npix=500)
        
    elif s == 'DHT21':
        d = Dendrogram.compute(clob.cube, min_value=6*sig,min_delta=sig, min_npix=1000)

        
    for cl in clouds:
        l_c = cld.ind_ang('get index', clob.lons, [cat['l_c__deg'][cl]])
        b_c = cld.ind_ang('get index', clob.lats, [cat['b_c__deg'][cl]])


        for tr in d.trunk:
            inds = np.array(tr.indices())
            coords = np.transpose(inds[1:]).tolist()            
            if [b_c, l_c] in coords:
#                 print(s, cl,'success')
                m = tr.get_mask()

                new_cube = clob.cube * m
                new_clob = clob.clone()
                new_clob.cube = new_cube

                cat['clob'][cl] = new_clob


# play = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
# ps = np.nanstd(play.cube)
# pd = Dendrogram.compute(play.cube, min_value=ps,min_delta=ps,min_npix=400)


tb1_long = pd.read_csv('Tb1_ParketAl_Carina.txt', delim_whitespace=True, skipfooter=4,header=[0,1])
tb2_long = pd.read_csv('Tb2_ParketAl_Carina.txt', delim_whitespace=True, skipfooter=5,header=[0,1])


temp = tb1_long.drop(3)
temp = temp.drop(8)
tb1 = temp.drop(21)


#CO bright vs #CO dark clouds
HIs1 = tb1['HI']['#']
HIs2 = tb2['HI']['#']

 
MsHI = (tb1['M_HI']['(M_sun)']).to_numpy() * u.M_sun
MsH2 = (tb2['M_H2']['(M_sun)']).to_numpy() * u.M_sun


Msum = np.zeros_like(MsHI)


for h, hi in enumerate(HIs1):
    npH2 = np.array(HIs2)
    here = np.where(npH2==hi)[0]
    add = np.nansum(MsH2[here])
    Msum[h] = MsHI[h] + add


less = ((Msum.value-MsHI.value)/MsHI.value==0)


plt.scatter(HIs1,(Msum.value-MsHI.value)/MsHI.value)
plt.xticks(rotation='vertical')
plt.hlines(0,'H1','H29')
plt.xlabel('HI Clouds')
plt.ylabel(r'Added Cloud Mass [$M_\odot$]')
plt.show()


len(HIs1)


len(MsHI)


plt.scatter(HIs1,MsHI.value,color='red')
plt.scatter(HIs1*less,MsHI.value*less,color ='black')

plt.xticks(rotation='vertical')
plt.hlines(1e6,'H1','H29')
plt.hlines(1e7,'H1','H29')

plt.xlabel('HI Clouds')
plt.ylabel(r'Added Cloud Mass [$M_\odot$]')
plt.yscale('log')
plt.show()


plt.scatter((Msum.value-MsHI.value)/MsHI.value *~less,MsHI.value *~less,color='red')
plt.scatter((Msum.value-MsHI.value)/MsHI.value * less,MsHI.value * less,color = 'black')

plt.xticks(rotation='vertical')


plt.ylabel('$M_{HI}$')
plt.xlabel(r'$M_{H_2}/M_{HI}$')
plt.yscale('log')
plt.show()


#another diagnostic that will be useful to look at in order to understand the dynamics of what is going
#on with the HI is the Mach number (look @ imara burkhart paper)
#if there are multiplt molecular clouds under a single name, it might not matter
#if they all under the same HI  cloud
#look up erik rosolowsky & alyssa goodman's work to see how they go about 
#applying dendrograms to solar neighborhood clouds / how they approach thresholding


RsHI = (tb1['r_geo_phys']['(pc)']).to_numpy() * u.pc


v_VT = np.sqrt(((c.G * Msum)/(4*RsHI))).to(u.km/u.s)


plt.scatter(Msum*~less, v_VT*~less,color='red')
plt.scatter(Msum*less, v_VT*less,color='black')

plt.xscale('log')
plt.xlabel(r'$M_{HI + H_2}$ [$M_\odot$]')
plt.ylabel('Virial Dispersion [$km\ s^{-1}$]')
plt.show()


sigs = (tb1['sigma_v']['(km/s)']).to_numpy() * u.km/u.s


plt.scatter(Msum*~less, sigs*~less,color='red',label='GMC-ful')
plt.scatter(Msum*less, sigs*less,color='black',label='GMC-less')

plt.xscale('log')
plt.xlabel(r'$M_{HI + H_2}$ [$M_\odot$]')
plt.ylabel(r'$\sigma_v\ [km\ s^{-1}]$')
# plt.hlines(1/np.sqrt(0.5),np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 0.5$')
# plt.hlines(1,np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 1$',linestyle='--')
# plt.hlines(1/(np.sqrt(2)),np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 2$',linestyle='-.')
# plt.hlines(1/2,np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 4$',linestyle=':')

plt.ylim(4,12)
plt.legend()
plt.show()


sigs = (tb1['sigma_v']['(km/s)']).to_numpy() * u.km/u.s


plt.scatter(Msum*~less, (v_VT/sigs)*~less,color='red')
plt.scatter(Msum*less, (v_VT/sigs)*less,color='black')

plt.xscale('log')
plt.xlabel(r'$M_{HI + H_2}$ [$M_\odot$]')
plt.ylabel(r'$v_{VT}/\sigma_v$')
plt.hlines(1/np.sqrt(0.5),np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 0.5$')
plt.hlines(1,np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 1$',linestyle='--')
plt.hlines(1/(np.sqrt(2)),np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 2$',linestyle='-.')
plt.hlines(1/2,np.nanmin(Msum.value), np.nanmax(Msum.value),label=r'$\alpha = 4$',linestyle=':')

plt.ylim(0,2)
plt.legend()
plt.show()


dsHI = (tb1['d']['(kpc)']).to_numpy() * 1e3 * u.pc
lsHI = (tb1['l_c']['(deg)']).to_numpy() * u.deg


sep = np.zeros_like(dsHI)
for i,inner in enumerate(dsHI[1:-1]):
    
    db = dsHI[i]
    da = dsHI[i+2]
    dell = (lsHI[i+2] - lsHI[i]).to(u.rad)
    
    s = 0.5*np.sqrt( da**2 + db**2 \
                    - 2*da*db*np.cos(dell))
    sep[i+1] = s
    
sepkpc = sep/1e3
bs = np.arange(0, np.nanmax(sepkpc.value),0.1)


mass = MsHI > 5e6*c.M_sun

plt.hist(sepkpc.value,bins=bs, edgecolor = 'goldenrod',facecolor='none')
plt.hist(sepkpc.value*mass,bins=bs, edgecolor = 'red',facecolor='none')

plt.xlabel('Separation [kpc]')
plt.ylabel('Number')
plt.xlim(0.1,2)
plt.ylim(0,10)
plt.xticks((0.5,1,1.5,2))
plt.show()
