#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))


# moving forward:
# start by reading in more detail the Carina Arm paper to look at their methodoloy
# look into clumpfind & decomposing the HI4PI dataset before settling on a specific approach for identifying the controld HI group
# 
# & in parallel start examining the GMC-ful HI clouds and see if i can reproduce plots / results that are already knonwn
# take a closer look at the methodology of papers that use the Dame CO data & use that as a way to determine a better way to fdefine the contours of the GMCs
# compare accumulation radius bounds with nia's / other published values as a sanity check


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
        ns = ['03','04','06','09','12','16','18','19','21','24','26','27','34','37']#figure out DHT01
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output


cat_wLupus = cat
cat = cat_wLupus.drop(10)
#lupus is from DHT01 which is not downloading properly


surveys = read_DHT()


cat.sort_values('survey')
#03, 04, 06, 09, 12, 16, 19, 24, 26, 34, 37 only have one cloud
#18 has two (but both cepheus? maybe this one should also be excluded bc of ambiguity)
#21 has 3
#27 has 4


DHT21 = surveys['DHT21']
DHT27 = surveys['DHT27']


np.nanstd(DHT27.cube)


d27 = Dendrogram.compute(DHT27.cube, min_value=3*0.489, min_npix=400) #this is the one i will go with
for i in d27.trunk:
    m = i.get_mask()
    plt.imshow(np.nansum(m*d27.data, axis = 0),**kw_imshow)
    plt.show()


print(len(d27.trunk))


d21 = Dendrogram.compute(DHT21.cube, min_value=6*0.382, min_npix=1000) #this is the one i will go with
for i in d21.trunk:
    m = i.get_mask()
    plt.imshow(np.nansum(m*d21.data, axis = 0),**kw_imshow)
    plt.show()


d21 = Dendrogram.compute(DHT21.cube, min_value=7*0.382, min_delta=5*0.382, min_npix=1000)
for i in d21.trunk:
    m = i.get_mask()
    plt.imshow(np.nansum(m*d21.data, axis = 0),**kw_imshow)
    plt.show()


#plots showing selection process -> narrowing down from Zucker


Z = pd.read_csv('zucker_ds',sep=',',skiprows=5)
Z_names = Z['Cloud']
Z_ds = Z['d']
Z_bs = (Z['b0']+Z['b1'])*0.5

d_mask = Z_ds.to_numpy() < 1000
b_mask = np.abs(Z_bs.to_numpy()) < 30
solar_neighbors = np.arange(1,len(Z)+1) * d_mask * b_mask
short = solar_neighbors[solar_neighbors!=0]
keep_bd = short-1


plt.scatter(Z_names, np.abs(Z_bs), marker = '+', color = 'black')
plt.xticks(rotation='vertical')
plt.hlines(30,Z_names[0], Z_names[26], color = 'grey', linestyle='--', label = r'$b$, DHT Survey')
plt.ylabel(r'Galactic Latitude [$\degree$]')
plt.legend()
plt.show()


plt.scatter(Z_names, Z_ds, marker = '+', color = 'black')
plt.xticks(rotation='vertical')
plt.hlines(1000,Z_names[0], Z_names[26], color = 'grey', linestyle='--', label = 'Solar Neighborhood boundary')
plt.ylabel('Heliocentric distance [pc]')
plt.legend()
plt.show()


surveys = read_DHT()


surveys


for s in surveys:
    name = s
    cl = surveys[s]
     
    if s!= 'DHT18':
        sig = np.nanstd(cl.cube)
        d = Dendrogram.compute(cl.cube,min_value=sig, min_delta=sig,min_npix=20)
        if len(d.trunk) < 10 and len(d.trunk) >=5:
            trs = d.trunk
            fig, axs = plt.subplots(1, len(trs))
            if len(trs) == 1:
                m = trs[0].get_mask()
                im = np.nansum(cl.cube * m, axis = 0)
                plt.imshow(im,**kw_imshow)
                plt.title(s)
            else:
                for t, tr in enumerate(trs):
                    m = tr.get_mask()
                    im = np.nansum(cl.cube * m, axis = 0)
                    axs[t].imshow(im,**kw_imshow)
                    axs[t].set_title(s+'  '+str(t))
            plt.show()
        else:
            print(s, 'has', len(d.trunk),'trunks')


for s in surveys:
    name = s
    cl = surveys[s]
     
    if s!= 'DHT18':
        sig = np.nanstd(cl.cube)
        d = Dendrogram.compute(cl.cube,min_value=sig, min_delta=sig,min_npix=20)
        if len(d.trunk) < 5:
            trs = d.trunk
            fig, axs = plt.subplots(1, len(trs))
            if len(trs) == 1:
                m = trs[0].get_mask()
                im = np.nansum(cl.cube * m, axis = 0)
                plt.imshow(im,**kw_imshow)
                plt.title(s)
            else:
                for t, tr in enumerate(trs):
                    m = tr.get_mask()
                    im = np.nansum(cl.cube * m, axis = 0)
                    axs[t].imshow(im,**kw_imshow)
                    axs[t].set_title(s+'  '+str(t))
            plt.show()
        else:
            print(s, 'has', len(d.trunk),'trunks')


for s in surveys:
    name = s
    cl = surveys[s]
     
    if s!= 'DHT18':
        sig = np.nanstd(cl.cube)
        d = Dendrogram.compute(cl.cube,min_value=3*sig, min_delta=3*sig,min_npix=20)
        print(len(d.trunk))
#     im = 
#     plt.imshow(im,**kw_imshow)
#     plt.title(s)
#     plt.show()


for s in surveys:
    name = s
    cl = surveys[s]
     
    if s!= 'DHT18':
        sig = np.nanstd(cl.cube)
        d = Dendrogram.compute(cl.cube,min_value=sig, min_delta=sig,min_npix=20)
        print(len(d.trunk))
#     im = 
#     plt.imshow(im,**kw_imshow)
#     plt.title(s)
#     plt.show()


for s in surveys:
    name = s
    cl = surveys[s]
     
    if s!= 'DHT18':
        d = Dendrogram.compute(cl.cube)
        print(len(d.trunk))
#     im = 
#     plt.imshow(im,**kw_imshow)
#     plt.title(s)
#     plt.show()


ns = ['03','04','06','09','12','16','18','19','21','24','26','27','34','37']#figure out DHT01

sigs = [0.30,]


for s in surveys:
    name = s
    cl = surveys[s]
     
    if s!= 'DHT18':
        d = Dendrogram.compute(cl.cube, min_npix=20)
        print(s, len(d.trunk))
#         if len(d.trunk) < 5:
#             trs = d.trunk
#             fig, axs = plt.subplots(1, len(trs))
#             if len(trs) == 1:
#                 m = trs[0].get_mask()
#                 im = np.nansum(cl.cube * m, axis = 0)
#                 plt.imshow(im,**kw_imshow)
#                 plt.title(s)
#             else:
#                 for t, tr in enumerate(trs):
#                     m = tr.get_mask()
#                     im = np.nansum(cl.cube * m, axis = 0)
#                     axs[t].imshow(im,**kw_imshow)
#                     axs[t].set_title(s+'  '+str(t))
#             plt.show()
#         else:
#             print(s, 'has', len(d.trunk),'trunks')
#     im = 
#     plt.imshow(im,**kw_imshow)
#     plt.title(s)
#     plt.show()


s = cat['survey']


d = Dendrogram.compute(cl.cube)


cl = surveys[s[0]]


m = d.trunk[0].get_mask()


fig, axs = plt.subplots(2,2,figsize=(10,5),dpi=100)
lb = [cl.lons[0], cl.lons[-1],cl.lats[0],cl.lats[-1]]
lv = [cl.lons[0], cl.lons[-1],cl.vels[0],cl.vels[-1]]

axs[0,0].imshow(np.nansum(m * cl.cube, axis = 0),**kw_imshow, extent = lb)
axs[0,1].imshow(np.nansum(m * cl.cube, axis = 1),**kw_imshow, extent = lv)

axs[1,0].imshow(np.nansum(cl.cube, axis = 0),**kw_imshow, extent = lb)
axs[1,1].imshow(np.nansum(cl.cube, axis = 1),**kw_imshow, extent = lv)

# plt.xlabel(r'l [$\degree$]')
# plt.ylabel(r'b [$\degree$]')

plt.show()


head = fits.open('DHT19_Camelopardalis_qdeg_mom.fits')[0].header


data = fits.open('DHT19.fits')[0].data


chead['CTYPE1'] = 'VELO-LSR'


fits.writeto('DHT19.fits', cdata,  chead, overwrite=True)


c1 = dc.unpack_fits_cube('DHT19.fits','DHT')


c2 = dc.unpack_fits_cube('DHT19_Camelopardalis_qdeg_mom.fits','DHT')


plt.imshow(np.nansum(c2.cube*c2.dv.value,axis=0),**kw_imshow)


plt.imshow(np.nansum(c1.cube*c1.dv.value,axis=0),**kw_imshow)


pl1.db


pl2.db


plt.imshow(np.nansum(pl1.cube*pl1.dv.value,axis=0),**kw_imshow)


or_mom = dc.unpack_fits_cube('DHT27_Orion_mom.fits','DHT')
or_raw = dc.unpack_fits_cube("DHT27_Orion_raw.fits", 'DHT')
or_int = dc.unpack_fits_cube('DHT27_Orion_interp.fits', 'DHT')


fig, axs = plt.subplots(1,3,figsize=(10,3), dpi=100)

axs[0].imshow(np.nansum(or_mom.cube*or_mom.dv.value, axis =0), **kw_imshow)
axs[1].imshow(np.nansum(or_raw.cube*or_raw.dv.value, axis =0), **kw_imshow)
axs[2].imshow(np.nansum(or_int.cube*or_int.dv.value, axis =0), **kw_imshow)

# plt.colorbar()
plt.show()


len(or_mom.cube), len(or_raw.cube), len(or_int.cube)


np.nanstd(or_mom.cube), np.nanstd(or_raw.cube), np.nanstd(or_int.cube)


cat


HI4PI.dv.to(u.km/u.s)


type(DHT_CO.data_header['T_RMS'])


CO_Ms = np.zeros(len(cat))*u.kg
CO_As = np.zeros(len(cat))*u.pc**2

HI_Ms = np.zeros(len(cat))*u.kg
HI_As = np.zeros(len(cat))*u.pc**2
HI_RAs = np.zeros(len(cat))*u.pc

HI_dispers = np.zeros(len(cat))*u.km/u.s


# %%time
for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    
#     plot_dend_mask(CO_clob, 3*sigma, sigma, 10)
    
    CO_spectrum = dc.spectrum1D(CO_clob.cube)
    
    dSig_CO, dA_CO = cld.surface_density(CO_clob, 'CO', show_Sigma=False)
    M_CO, A_CO = cld.mass_and_area(n, cat, dSig_CO, dA_CO)
    CO_Ms[n] = M_CO
    CO_As[n] = A_CO


    d = cat['d__pc'][n] * d_unit
    RA, s = cld.calculate_accumulation_radius(M_CO, d)
    HI_RAs[n] = RA
    
    peak_v_CO = CO_clob.vels[np.nanargmax(CO_spectrum)]
    
    HI_vels_to_CO_vels = HI4PI.units['vels'].to(DHT_CO.units['vels'])
    HI_vels_COunit = HI4PI.vels * HI_vels_to_CO_vels
    peak_v_HI = np.where(HI_vels_COunit==peak_v_CO)[0]

    b_cent = 0.5*(b0+b1)
    l_cent = 0.5*(l0+l1)
    

    HI_clob = HI4PI.create_subcube([b_cent-s, b_cent+s, l_cent-s, l_cent+s], v_guess=peak_v_HI) 
    HI_spectrum = dc.spectrum1D(HI_clob.cube)
    HI_sig = dc.dispersion_1sigma(HI_spectrum) * HI_clob.dv.to(u.km/u.s)
#     plt.plot(HI_clob.vels, HI_spectrum)
#     plt.show()
    HI_dispers[n] = HI_sig
    
    dSig_HI, dA_HI = cld.surface_density(HI_clob, '21cm', show_Sigma=False)
    M_HI, A_HI = cld.mass_and_area(n, cat, dSig_HI, dA_HI)
    
    HI_Ms[n] = M_HI
    HI_As[n] = A_HI

    
#     plt.plot(CO_clob.vels,20*CO_spectrum, marker = '+', label='CO')
#     plt.plot(HI_clob.vels * HI_vels_to_CO_vels,HI_spectrum, marker = '+',label = 'HI')
#     plt.legend()
    
#     plt.title(names[n])
#     plt.xlabel(r'Velocity [$km\ s^{-1}$]')
#     plt.ylabel(r'T [$K$]')
#     plt.show()

    
    fig, axs = plt.subplots(1,2,figsize=(15,3), dpi=200)
 
    im = dSig_HI.value
    imx = HI_clob.lons 
    imy = HI_clob.lats


    vmin,vmax = zlims(im) 

    axs[0].imshow(im, vmin=vmin, vmax=vmax, \
                extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)
    

    im2 = dSig_CO.value
    im2x = CO_clob.lons
    im2y = CO_clob.lats
    
    vmin2,vmax2 = zlims(im2)

    axs[1].imshow(im2, vmin=vmin2, vmax=vmax2, \
                extent = [im2x[0],im2x[-1], im2y[0],im2y[-1]],**kw_imshow)
#     plt.colorbar(orientation='vertical')#, label = r'$W_{HI}$ $(K\times km/s)$')

    axs[0].contour(im2x, im2y, im2, colors = 'orange', levels = [1.5*np.nanmin(im2), 0.1* np.nanmax(im2),0.5*np.nanmax(im2)], \
                linewidths = [0.3, 0.3, 1], linestyles=['--','solid','solid'])


    plt.ylabel('Galactic Latitude')
    plt.xlabel('Galactic Longitude')


    plt.show()

    
# %%time
for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    
#     plot_dend_mask(CO_clob, 3*sigma, sigma, 10)
    
    CO_spectrum = dc.spectrum1D(CO_clob.cube)
    
    dSig_CO, dA_CO = cld.surface_density(CO_clob, 'CO', show_Sigma=False)
    M_CO, A_CO = cld.mass_and_area(n, cat, dSig_CO, dA_CO)
    CO_Ms[n] = M_CO
    CO_As[n] = A_CO


    d = cat['d__pc'][n] * d_unit
    RA, s = cld.calculate_accumulation_radius(M_CO, d)
    HI_RAs[n] = RA
    
    peak_v_CO = CO_clob.vels[np.nanargmax(CO_spectrum)]
    
    HI_vels_to_CO_vels = HI4PI.units['vels'].to(DHT_CO.units['vels'])
    HI_vels_COunit = HI4PI.vels * HI_vels_to_CO_vels
    peak_v_HI = np.where(HI_vels_COunit==peak_v_CO)[0]

    b_cent = 0.5*(b0+b1)
    l_cent = 0.5*(l0+l1)
    

    HI_clob = HI4PI.create_subcube([b_cent-s, b_cent+s, l_cent-s, l_cent+s], v_guess=peak_v_HI) 
    HI_spectrum = dc.spectrum1D(HI_clob.cube)
    HI_sig = dc.dispersion_1sigma(HI_spectrum) * HI_clob.dv.to(u.km/u.s)
#     plt.plot(HI_clob.vels, HI_spectrum)
#     plt.show()
    HI_dispers[n] = HI_sig
    
    dSig_HI, dA_HI = cld.surface_density(HI_clob, '21cm', show_Sigma=False)
    M_HI, A_HI = cld.mass_and_area(n, cat, dSig_HI, dA_HI)
    
    HI_Ms[n] = M_HI
    HI_As[n] = A_HI

    
#     plt.plot(CO_clob.vels,20*CO_spectrum, marker = '+', label='CO')
#     plt.plot(HI_clob.vels * HI_vels_to_CO_vels,HI_spectrum, marker = '+',label = 'HI')
#     plt.legend()
    
#     plt.title(names[n])
#     plt.xlabel(r'Velocity [$km\ s^{-1}$]')
#     plt.ylabel(r'T [$K$]')
#     plt.show()

    
    fig, axs = plt.subplots(1,2,figsize=(15,3), dpi=200)
 
    im = np.nansum(HI_clob.cube*HI_clob.dv.value, axis = 0)#dSig_HI.value
    imx = HI_clob.lons 
    imy = HI_clob.lats


    vmin,vmax = zlims(im) 

    axs[0].imshow(im, vmin=vmin, vmax=vmax, \
                extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)
    

    im2 = np.nansum(CO_clob.cube*CO_clob.dv.value, axis = 0)#dSig_CO.value
    im2x = CO_clob.lons
    im2y = CO_clob.lats
    
    vmin2,vmax2 = zlims(im2)

    axs[1].imshow(im2, vmin=vmin2, vmax=vmax2, \
                extent = [im2x[0],im2x[-1], im2y[0],im2y[-1]],**kw_imshow)
#     plt.colorbar(orientation='vertical')#, label = r'$W_{HI}$ $(K\times km/s)$')

    axs[0].contour(im2x, im2y, im2, colors = 'orange', levels = [1.5*np.nanmin(im2), 0.1* np.nanmax(im2),0.5*np.nanmax(im2)], \
                linewidths = [0.3, 0.3, 1], linestyles=['--','solid','solid'])


    plt.ylabel('Galactic Latitude')
    plt.xlabel('Galactic Longitude')


    plt.show()

    
dc.dispersion_1sigma(HI_spectrum)*HI_clob.dv


plt.scatter(HI_clob.vels,HI_spectrum)


HI_dispers


xs = np.arange(3,5.5,0.1)


ys = 3*xs


plt.figure(figsize=(8,8))
plt.scatter(np.log(np.sqrt(HI_As/np.pi).value), np.log(HI_Ms.to(c.M_sun).value), marker = '+')
plt.ylabel(r'$log$ [Cloud Mass / $M_\odot$]')
plt.xlabel(r'$log$ [Cloud Size / pc]')
# plt.xscale('log')
# plt.yscale('log')
plt.plot(xs, ys, linestyle = '--', color = 'red', label = 'Fit line: $R = M^3$')
plt.legend()
plt.show()


plt.figure(figsize=(8,8))
plt.scatter(np.log(np.sqrt(HI_As/np.pi).value), HI_alphas, marker = '+')
# plt.ylabel(r'$log$ [Cloud Mass / $M_\odot$]')
# plt.xlabel(r'$log$ [Cloud Size / pc]')
# plt.xscale('log')
# plt.yscale('log')
# plt.plot(xs, ys, linestyle = '--', color = 'red', label = 'Fit line: $R = M^3$')
plt.legend()
plt.show()


def virial_parameter(sig, M, R):
    num = 5*R*sig**2
    den = c.G * M
    return num/den


HI_alphas = (virial_parameter(HI_dispers, HI_Ms, np.sqrt(HI_As/np.pi))).si


plt.hist(HI_alphas)


np.sqrt(HI_As/np.pi)


plt.figure(figsize=(8,8))
plt.scatter(HI_RAs, np.sqrt(HI_As/np.pi), marker = '+')
plt.xlabel(r'$H_2$ mass [$M_\odot$]')
plt.ylabel(r'HI mass [$M_\odot$]')
# plt.xscale('log')
# plt.yscale('log')
plt.show()


plt.scatter(names, CO_Ms.to(c.M_sun)/1e4, marker = '+')
plt.xticks(rotation='vertical')
plt.show()


names[17], CO_Ms[17].to(c.M_sun)/1e4, np.sqrt(CO_As[17]/np.pi)


names[18], CO_Ms[18].to(c.M_sun)/1e4, np.sqrt(CO_As[18]/np.pi)


names[24], CO_Ms[24].to(c.M_sun)/1e4, np.sqrt(CO_As[24]/np.pi)


names[20], CO_Ms[20].to(c.M_sun)/1e4, np.sqrt(CO_As[20]/np.pi)


names[10], CO_Ms[10].to(c.M_sun)/1e4, np.sqrt(CO_As[10]/np.pi)


names[14], CO_Ms[14].to(c.M_sun)/1e4, np.sqrt(CO_As[14]/np.pi)


names[2], CO_Ms[2].to(c.M_sun)/1e4, np.sqrt(CO_As[2]/np.pi)


names[21], CO_Ms[21].to(c.M_sun)/1e4, np.sqrt(CO_As[21]/np.pi)


names[16], CO_Ms[16].to(c.M_sun)/1e4, np.sqrt(CO_As[16]/np.pi)


# R_As = np.zeros(len(names))
# Ms = np.zeros(len(names))


np.nanstd(HI_clob.cube)


orion_dend_L = Dendrogram.compute(HI_clob.cube,\
                               min_value = 3*np.nanstd(HI_clob.cube),\
                               min_delta = np.nanstd(HI_clob.cube),\
                               min_npix = 10)


orion_dend = Dendrogram.compute(HI_clob.cube,\
                               min_value = 3*np.nanstd(HI_clob.cube),\
                               min_delta = np.nanstd(HI_clob.cube),\
                               min_npix = 10)


masked_L = orion_dend_L.trunk[0].get_mask() * HI_clob.cube


orion_dend_bb_L = Dendrogram.compute(CO_clob.cube,\
                               min_value = 3*np.nanstd(CO_clob.cube),\
                               min_delta = np.nanstd(CO_clob.cube),\
                               min_npix = 10)


orion_dend_bb = Dendrogram.compute(CO_clob.cube,\
                               min_value = 3*np.nanstd(CO_clob.cube),\
                               min_delta = np.nanstd(CO_clob.cube),\
                               min_npix = 10)


for k in orion_dend_bb_L.trunk:
    print(k.get_npix())


masked2 = orion_dend_bb.trunk[1].get_mask() * CO_clob.cube


masked2_L = orion_dend_bb_L.trunk[2].get_mask() * CO_clob.cube


HI_clob.cube


fig, axs = plt.subplots(1, 2, figsize=(10,3), dpi=200)


im = np.nansum(masked_L*HI_clob.dv.value, axis = 0)#dSig_HI.value
imp = np.nansum(masked*HI_clob.dv.value, axis = 0)#dSig_HI.value

imx = HI_clob.lons 
imy = HI_clob.lats


vmin,vmax = zlims(im) 

# axs[0].imshow(im, vmin=vmin, vmax=vmax, \
#             extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


# # axs[0].contour(imx, imy, imp, colors = 'orange', levels = [np.nanmin(im2)], \
# #             linewidths = [1])


# axs[1].imshow(imp, vmin=vmin, vmax=vmax, \
#             extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)

# axs[1].contour(imx, imy, im, colors = 'orange', levels = [np.nanmin(im2)], \
#             linewidths = [1])

# plt.show()
im2 = np.nansum(masked2*CO_clob.dv.value, axis = 0)#dSig_CO.value
im2L = np.nansum(masked2_L*CO_clob.dv.value, axis = 0)#dSig_CO.value

im2x = CO_clob.lons
im2y = CO_clob.lats

vmin,vmax = zlims(im2)

axs[0].imshow(im2, vmin=vmin, vmax=vmax, \
            extent = [im2x[0],im2x[-1], im2y[0],im2y[-1]],**kw_imshow)

axs[1].imshow(im2L, vmin=vmin, vmax=vmax, \
            extent = [im2x[0],im2x[-1], im2y[0],im2y[-1]],**kw_imshow)
# plt.colorbar(orientation='vertical')#, label = r'$W_{HI}$ $(K\times km/s)$')

# plt.contour(im2x, im2y, im2, colors = 'orange', levels = [np.nanmin(im2)], \
#             linewidths = [1])


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')


plt.show()


df = comparison.dropna()


del_l_1 = df['l0_Z'] - df['l1_L']


del_l_2 = df['l0_L'] - df['l1_Z']


del_b_1 = df['b0_Z'] - df['b0_L']


del_b_1


del_b_2


del_b_2 = df['b1_Z'] - df['b1_L']


del_l_1[7] +=360


df


plt.figure(figsize=(8,8))
plt.scatter(df['d_Z'], del_b_1, marker = '+')
plt.scatter(df['d_L'], del_b_1, marker = 'x')
plt.scatter(df['d_L'], del_b_2, marker = '+')
plt.scatter(df['d_Z'], del_b_2, marker = 'x')
plt.show()


plt.figure(figsize=(8,8))
# plt.scatter(df['d_Z'], del_l_1%360, marker = '+')
plt.scatter(df['d_L'], del_l_1%360, marker = 'x')
plt.scatter(df['d_L'], del_l_2%360, marker = '+')
# plt.scatter(df['d_Z'], del_l_2%360, marker = 'x')
plt.show()


-1/np.abs(-1)


x = -371902347820


(np.abs(x)+6)* (x/np.abs(x))


-371902347820/np.abs(371902347820)


del_l_2


np.median(del_l_2%360)


np.nanmax(masked)


np.nanmax(masked2)


# %%time
for n in range(len(Zucker_df)):
    b0 = Zucker_df['b0'][n]
    b1 = Zucker_df['b1'][n]
    
    l0 = Zucker_df['l0'][n]
    l1 = Zucker_df['l1'][n]
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    CO_spectrum = dc.spectrum1D(CO_clob.cube)
    
    dSig_CO, dA_CO = cld.surface_density(CO_clob, 'CO', show_Sigma=False)
    M, A = cld.mass_and_area(n, Zucker_df, dSig_CO, dA_CO)
    Ms[n] = M.to(c.M_sun).value

    d = Zucker_df['d'][n] * u.Unit(Zucker_df['d_unit'][n])
    RA, s = cld.calculate_accumulation_radius(M, d)
    R_As[n] = RA.value
    
    peak_v_CO = CO_clob.vels[np.nanargmax(CO_spectrum)]
    
    HI_vels_to_CO_vels = HI4PI.units['vels'].to(DHT_CO.units['vels'])
    HI_vels_COunit = HI4PI.vels * HI_vels_to_CO_vels
    peak_v_HI = np.where(HI_vels_COunit==peak_v_CO)[0]

    b_cent = 0.5*(b0+b1)
    l_cent = 0.5*(l0+l1)
    

    HI_clob = HI4PI.create_subcube([b_cent-s, b_cent+s, l_cent-s, l_cent+s], v_guess=peak_v_HI) 
    HI_spectrum = dc.spectrum1D(HI_clob.cube)
    
    
#     plt.plot(CO_clob.vels,20*CO_spectrum, marker = '+', label='CO')
#     plt.plot(HI_clob.vels * HI_vels_to_CO_vels,HI_spectrum, marker = '+',label = 'HI')
#     plt.legend()
    
#     plt.title(names[n])
#     plt.xlabel(r'Velocity [$km\ s^{-1}$]')
#     plt.ylabel(r'T [$K$]')
#     plt.show()

    dSig, dA = cld.surface_density(HI_clob, '21cm', show_Sigma=False)
    
    
#     plt.figure(figsize=(15,3), dpi=200)
    
 
#     im = dSig.value
#     imx = HI_clob.lons 
#     imy = HI_clob.lats


#     vmin,vmax = zlims(im) 

#     plt.imshow(im, vmin=vmin, vmax=vmax, \
#                 extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)
    

#     im2 = dSig_CO.value
#     im2x = CO_clob.lons
#     im2y = CO_clob.lats
    
#     vmin,vmax = zlims(im2)

#     plt.imshow(im2, vmin=vmin, vmax=vmax, \
#                 extent = [im2x[0],im2x[-1], im2y[0],im2y[-1]],**kw_imshow)
#     plt.colorbar(orientation='vertical')#, label = r'$W_{HI}$ $(K\times km/s)$')

#     plt.contour(im2x, im2y, im2, colors = 'orange', levels = [1.5*np.nanmin(im2), 0.1* np.nanmax(im2),0.5*np.nanmax(im2)], \
#                 linewidths = [0.3, 0.3, 1], linestyles=['--','solid','solid'])


#     plt.ylabel('Galactic Latitude')
#     plt.xlabel('Galactic Longitude')


#     plt.show()

    
# for n in range(len(names)):
#     print(names[n], R_As[n], Ms[n]/1e4)


comparison = pd.read_csv('comparing_lb_coords',sep=',',skiprows=5)


df = comparison


for k in df.index:
    plt.figure(figsize=(5,5))
#     print(df['b0_Z'][k], df['l0_Z'][k], df['l1_Z'][k])
#     print(df['b0_L'][k], df['l0_L'][k], df['l1_L'][k])
#     print('\n')
    plt.hlines(df['b0_Z'][k], df['l0_Z'][k], df['l1_Z'][k], label = 'Z', color = 'black')
    plt.hlines(df['b1_Z'][k], df['l0_Z'][k], df['l1_Z'][k], color = 'black')
    
    plt.vlines(df['l0_Z'][k], df['b0_Z'][k], df['b1_Z'][k], color = 'black')
    plt.vlines(df['l1_Z'][k], df['b0_Z'][k], df['b1_Z'][k], color = 'black')
    
    
    plt.hlines(df['b0_L'][k], df['l0_L'][k], df['l1_L'][k], label = 'L', color = 'red')
    plt.hlines(df['b1_L'][k], df['l0_L'][k], df['l1_L'][k], color = 'red')
    
    plt.vlines(df['l0_L'][k], df['b0_L'][k], df['b1_L'][k], color = 'red')
    plt.vlines(df['l1_L'][k], df['b0_L'][k], df['b1_L'][k], color = 'red')
    
    plt.title(df['Cloud'][k])
    plt.legend()
    
#     plt.xlim(-180, 360)
#     plt.ylim(-40,40)
#     plt.xticks(np.arange(-180, 360, 10))
    plt.show()


R_As[17]


l0


Zucker_df['l0'][24], Zucker_df['l1'][24]


names[24]


l_cent


s


171+22-360


HI_clob.lons[0], HI_clob.lons[-1]
