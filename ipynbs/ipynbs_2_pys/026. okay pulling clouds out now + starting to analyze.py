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

DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')


#----------------------
cat = pd.read_csv('Mekuria_catalog.csv',sep=',',skiprows=5)
names = cat['name']


l_unit = u.deg
b_unit = u.deg
d_unit = u.pc


def extend_bounds(x, add):
    sign = x/np.abs(x)
    extended = np.abs(x) + add
    return extended*sign


def plot_dend_mask(clob, min_v, min_del, min_npx):
    d = Dendrogram.compute(clob.cube, min_value = min_v, min_delta = min_del, min_npix = min_npx)
    px_count = np.zeros(len(d.trunk))
#     print('there are', len(d.trunk),'trunks')
    for i in range(len(d.trunk)):
        px_count[i] = d.trunk[i].get_npix()
    plt.plot(px_count, marker = '+')
    plt.show()
#     print('and they have',px_count,'\n')
    
#     this = np.nanargmax(px_count)
    
#     mask = d.trunk[this].get_mask()
#     masked_cube = clob.cube*mask
    
#     imx = clob.lons
#     imy = clob.lats
#     plt.imshow(np.nansum(masked_cube*clob.dv.value, axis = 0),extent = [imx[0],imx[-1], imy[0],imy[-1]], **kw_imshow)
#     plt.show()


d


def spot_test(n):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6

    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6

    clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    d = Dendrogram.compute(clob.cube, 3*0.28, 5*0.28, 20)

    for i in range(len(d.trunk)):
        print(i)
        m = d.trunk[i].get_mask()
        md = clob.cube*m
        plt.imshow(np.nansum(md*clob.dv.value, axis = 0), **kw_imshow)
        plt.show()
        
    return d


spot_test(17)


dir(od)


for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    print(n, names[n])
    plot_dend_mask(CO_clob, sigma, sigma, 10)
   

for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    print(n, names[n])
    plot_dend_mask(CO_clob, 3*sigma, 5*sigma, 20)
   

for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    print(n, names[n])
    plot_dend_mask(CO_clob, sigma, sigma, 10)
   

for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
        
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    x,y,z = CO_clob.cube.shape
    print(n,x,y,z, x*y*z)


for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    print(n, names[n])
    plot_dend_mask(CO_clob, 3*sigma, sigma, 100)
   

for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    print(n, names[n])
    plot_dend_mask(CO_clob, 3*sigma, sigma, 10)
   

# %%time
for n in range(len(cat)):
    b0 = cat['b_c__deg'][n] + 6
    b1 = cat['b_c__deg'][n] - 6
    
    l0 = cat['l_c__deg'][n] + 6
    l1 = cat['l_c__deg'][n] - 6
    
    
    CO_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    sigma = 0.28
    
    plot_dend_mask(CO_clob, 3*sigma, sigma, 10)
    
    CO_spectrum = dc.spectrum1D(CO_clob.cube)
    
    dSig_CO, dA_CO = cld.surface_density(CO_clob, 'CO', show_Sigma=False)
    M, A = cld.mass_and_area(n, cat, dSig_CO, dA_CO)
#     Ms[n] = M.to(c.M_sun).value

    d = cat['d__pc'][n] * d_unit
    RA, s = cld.calculate_accumulation_radius(M, d)
#     R_As[n] = RA.value
    
    peak_v_CO = CO_clob.vels[np.nanargmax(CO_spectrum)]
    
    HI_vels_to_CO_vels = HI4PI.units['vels'].to(DHT_CO.units['vels'])
    HI_vels_COunit = HI4PI.vels * HI_vels_to_CO_vels
    peak_v_HI = np.where(HI_vels_COunit==peak_v_CO)[0]

    b_cent = 0.5*(b0+b1)
    l_cent = 0.5*(l0+l1)
    

    HI_clob = HI4PI.create_subcube([b_cent-s, b_cent+s, l_cent-s, l_cent+s], v_guess=peak_v_HI) 
    HI_spectrum = dc.spectrum1D(HI_clob.cube)
    dSig_HI, dA_HI = cld.surface_density(HI_clob, '21cm', show_Sigma=False)

    
    plt.plot(CO_clob.vels,20*CO_spectrum, marker = '+', label='CO')
    plt.plot(HI_clob.vels * HI_vels_to_CO_vels,HI_spectrum, marker = '+',label = 'HI')
    plt.legend()
    
    plt.title(names[n])
    plt.xlabel(r'Velocity [$km\ s^{-1}$]')
    plt.ylabel(r'T [$K$]')
    plt.show()

    
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

    
comparison.dropna()


Zucker_df


Zucker_df['l0']


Zucker_df['l0']%360


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
