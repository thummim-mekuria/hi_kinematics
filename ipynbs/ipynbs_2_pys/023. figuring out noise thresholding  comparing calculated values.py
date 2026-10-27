#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


# questions for nia:
# 1. how do you know when a spectrum is sampled well (across how many velocity channels?)
# 2. issues i ran into when calculating R_A
#     - related: how do i handle conversion between physical lengths (pc) and angular (esp when considering that there are maybe 2 layers of projection? from the flat plane of infinity -> spherical sky -> cylindrical coordinates?)
#     #factor of 2
#     #for mass calculations, look at the reported values in other papers & compare to the value i calculated
#     #for the mass values, because Zucker's distances are more updated, trust the mass values that are calculated using those
# 3. help with talking through kinematic distance calculation


# %%time

import datacube as dc
import clouds as cld

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c

from astropy.visualization import ZScaleInterval as zscale


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')

DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')


#----------------------

Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)


# def subcube_sigma(clob, coords):
#     b0, b1, l0, l1 = coords
  
#     dlon_buffer = 2*np.abs(clob.dl.value) 
#     lons_dg = clob.lons * clob.units['lons'].to(u.degree) 

#     if np.any(lons_dg > 180 + dlon_buffer): #glon range is (0, 360)
#         glon_con = '0-360' 
#     else: #glon range is +/- 180 
#         glon_con = '+/- 180'
        

#     #2) check database convention
#     if glon_con == '+/- 180':
#         if l0 > 180 + dlon_buffer:
#             l0 -= 360

#         if l1 > 180 + dlon_buffer:
#             l1 -= 360

#     if glon_con == '0-360':
#         if l0 < 0  - dlon_buffer:
#             l0 += 360

#         if l1 < 0 - dlon_buffer:
#             l1 += 360

#     #convert physical coords -> index locations within the larger cube
#     lat_inds = cld.ind_ang('get index', clob.lats, (b0,b1))
#     lon_inds = cld.ind_ang('get index', clob.lons, (l0,l1))

#     #fixes index orders in case data is in decreasing order
#     lat_inds.sort()
#     lon_inds.sort()

#     b0, b1 = lat_inds #get each element separately
#     l0, l1 = lon_inds

#     scube = clob.cube[:, b0:b1, l0:l1]
#     sublats = clob.lats[b0:b1]
#     sublons = clob.lons[l0:l1]
    
#     sigma = np.nanstd(scube)
    
#     return sigma


def compare_areas(clob, n, database):
#     #------------------ 3D (the .pys handle this) -----------------------
    b0, b1 = database['b0'][n], database['b1'][n]
    l0, l1 = database['l0'][n], database['l1'][n]
    coords = [b0, b1, l0, l1]
    
    
    with3D = clob.create_subcube(coords)
    
    dSig3D, dA3D = cld.surface_density(with3D, 'CO', show_Sigma=False)
    M3D, A3D = cld.mass_and_area(n, database, dSig3D, dA3D)
    
    
    #------------------------------------ 2D ------------------------------------------
    #------------------ setting up ------------------------

    b0, b1, l0, l1 = coords
  
    dlon_buffer = 2*np.abs(clob.dl.value) 
    lons_dg = clob.lons * clob.units['lons'].to(u.degree) 

    if np.any(lons_dg > 180 + dlon_buffer): #glon range is (0, 360)
        glon_con = '0-360' 
    else: #glon range is +/- 180 
        glon_con = '+/- 180'
        

    #2) check database convention
    if glon_con == '+/- 180':
        if l0 > 180 + dlon_buffer:
            l0 -= 360

        if l1 > 180 + dlon_buffer:
            l1 -= 360

    if glon_con == '0-360':
        if l0 < 0  - dlon_buffer:
            l0 += 360

        if l1 < 0 - dlon_buffer:
            l1 += 360

    #convert physical coords -> index locations within the larger cube
    lat_inds = cld.ind_ang('get index', clob.lats, (b0,b1))
    lon_inds = cld.ind_ang('get index', clob.lons, (l0,l1))

    #fixes index orders in case data is in decreasing order
    lat_inds.sort()
    lon_inds.sort()

    b0, b1 = lat_inds #get each element separately
    l0, l1 = lon_inds

    scube = clob.cube[:, b0:b1, l0:l1]
    sublats = clob.lats[b0:b1]
    sublons = clob.lons[l0:l1]
    
    #------------------ now do thresholding ------------------------

    sigma_whole = np.nanstd(clob.cube)
    
    zeroth = np.nansum(scube*clob.dv.value, axis = 0)
    
    mask2D = zeroth > sigma_whole * 3 * clob.dv.value #3sigma cutoff, same as other one
    
    pxs_with_signal = np.ones_like(zeroth) * mask2D
    d = database['d'][n]
    dA2D = d**2 * pxs_with_signal * clob.dl.to(u.rad) * clob.db.to(u.rad)
    
    A2D = np.sum(dA2D)
    

    return A3D,  A2D


names = Zucker_df['Cloud']


As3D = np.zeros(len(Zucker_df))
As2D = np.zeros(len(Zucker_df))

for n in range(len(Zucker_df)):
    A3D, A2D = compare_areas(DHT_CO, n,  Zucker_df)
    As3D[n] = A3D.value
    As2D[n] = A2D.value


fig, axs = plt.subplots(1,2, figsize = (15, 5))

axs[0].vlines(names, As2D, As3D, color = 'orange', alpha = 0.5)
axs[0].vlines(names, 10, As3D, color = 'orange', linestyle = '--', alpha = 0.3)


axs[0].scatter(names, As3D , marker = '+', color = 'black', label = r'$\sigma$ cut on ppv')
axs[0].scatter(names, As2D , marker = '+', color = 'green', label = r'$\sigma$ cut on 0th moment')


axs[0].set_xticks(names, labels = names, rotation = 'vertical')

axs[0].set_ylabel(r'cloud area [$pc^{2}$]')

axs[0].set_yscale('log')
axs[0].legend()

# ------------------
res = 100* np.abs(As3D - As2D)/As3D

axs[1].vlines(names, 0, res, color = 'orange', linestyle = '--', alpha = 0.3)
axs[1].hlines(10, names[0], names[24], color = 'purple', linestyle = '--', alpha = 0.3, label = '10%')
axs[1].hlines(50, names[0], names[24], color = 'purple', linestyle = '-.', alpha = 0.3, label = '50%')


axs[1].scatter(names, res , marker = '+', color = 'black')

axs[1].set_xticks(names, labels = names, rotation = 'vertical')

axs[1].set_ylabel(r'percentage difference')

# axs[1].set_yscale('log')

# axs[1].legend()

plt.show()


#look through nia's papers to find a section that details the steps of how the
#sigma cutoff is calculated
#look at observational papers by john lewis, lada, lombardi, adam leroy, erik rosolowsky
#& look at their data analysis / data reduction sections
#papers that explicitly talk about identifying molecular clouds are a good place to look


#to help troubleshoot what's going on with the sigma cut-off, write out what i 
#am doing in words


# sigmas = np.zeros(len(Zucker_df))
# for n in range(len(Zucker_df)):
#     b0 = np.array(Zucker_df['b0'][n])
#     b1 = np.array(Zucker_df['b1'][n])


#     l0 = np.array(Zucker_df['l0'][n])
#     l1 = np.array(Zucker_df['l1'][n])
    
#     this_sigma = subcube_sigma(DHT_CO, [b0, b1, l0, l1])
    
#     sigmas[n] = this_sigma

# # --------------------------------

# plt.scatter(names, sigmas, marker = '+', color = 'black', \
#             label = 'calculated from subcube ')
# # plt.hlines(0.25, 0, len(names), linestyle = '--', color = 'orange', \
# #            label = 'median value from Dame 2001')
# plt.hlines(np.nanstd(DHT_CO.cube), 0, len(names), linestyle = '--',\
#            color = 'green', label = 'calculated from full DHT cube')
# plt.xticks(rotation = 'vertical')

# plt.legend(loc='upper left')
# plt.ylabel( r'$\sigma_{rms}$ [K]')

# plt.title(r'$\sigma_{rms}$ values')

# plt.show()


#for this ^ how to go about choosing range of velocities where there is data if the thresholding method
#we're using involves first calculating 0th moment of the subcube
#ie how do you create a mask for the 3D cube ??
#also the 3D way of doing it retains less data -> is this good ('cleaner'?) or bad(losing important info)


#calculating distances kinematically to compare to literature values
def R_poly (R, l, b, v_lsr, theta_sun, r_sun, a1, a2, a3):

    term1 = R/r_sun
    term2 = (R/r_sun)**a2
    
    coeff1 = (v_lsr/np.cos(b))+ (theta_sun*np.sin(l))
    coeff2 = a1*theta_sun*np.sin(l)
    coeff3 = a3*theta_sun*np.sin(l)
             
    final = coeff1*term1 - coeff2*term2 - coeff3
    
    return final


# def dist(R, r_sun, l, b):
#     print(R, r_sun)
#     term1 = r_sun*np.cos(l)/np.cos(b)
    
#     term2_1 = R**2
#     term2_2 = (r_sun*np.sin(l))**2
#     term2 = np.sqrt(term2_1 - term2_2)/np.cos(b)
    
    
#     final = term1 + term2
#     if final

# #     if R > r_sun:
# #         final = term1 + term2
    
# #     else:
# #         print('in else')
# #         final = (term1 + term2, term1 - term2)
    
#     return final


def dist_yarn(R, r_sun, l, b): #this and dist give pretty similar results, <0.2% differences
    term1 = r_sun*np.cos(b)*np.cos(l)
    
    term2_1 = term1**2
    term2_2 = (r_sun**2)
    term2_3 = R**2
    
    term2 = np.sqrt(term2_1 - term2_2 + term2_3)
    
    
    final = term1 - term2
    
    return final


#------------------------------------------------

# #values from Russeil et al 2017
# Theta_sun = 240*u.km/u.s
# R_sun = 8.34*1000*u.pc
# a1 = 1.022
# a2 = 0.0803
# a3 = 0.


#values from Brand & Blitz 1993
Theta_sun = 220*u.km/u.s
R_sun = 8.5*1000*u.pc
a1 = 1.00767
a2 = 0.0394
a3 = 0.00712

#------------------------------------------------

possible_Rs = np.arange(0, 17, 0.05)*1000*u.pc


my_Zucker_dists_near = np.zeros(len(Zucker_df))
for n in range(len(Zucker_df)):
    b0 = np.array(Zucker_df['b0'][n])
    b1 = np.array(Zucker_df['b1'][n])


    l0 = np.array(Zucker_df['l0'][n])
    l1 = np.array(Zucker_df['l1'][n])
    
    this_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    vels = this_clob.vels * this_clob.units['vels']
    L = np.nanmedian(this_clob.lons) * this_clob.units['lons']
    B = np.nanmedian(this_clob.lats) * this_clob.units['lats']
    
    r_out = np.zeros((len(possible_Rs), len(vels)))
    
    for r, R in enumerate (possible_Rs):
        this_R = R
        row = np.zeros(len(vels))
        
        for v, V in enumerate(vels):
            r_value = R_poly(R, L, B, V, Theta_sun, R_sun, a1, a2, a3)
            row[v] = r_value.value
            
        r_out[r] = row
        
    
    r_out_med = np.nanmedian(r_out, axis = 1)

    closest_R = possible_Rs[np.argmin(np.abs(0-r_out_med))]

    d = dist_yarn(closest_R, R_sun, L, B)
    
    my_Zucker_dists_near[n] = d.value
            

my_Zucker_dists = np.zeros(len(Zucker_df))
for n in range(len(Zucker_df)):
    b0 = np.array(Zucker_df['b0'][n])
    b1 = np.array(Zucker_df['b1'][n])


    l0 = np.array(Zucker_df['l0'][n])
    l1 = np.array(Zucker_df['l1'][n])
    
    this_clob = DHT_CO.create_subcube([b0, b1, l0, l1])
    
    vels = this_clob.vels * this_clob.units['vels']
    L = np.nanmedian(this_clob.lons) * this_clob.units['lons']
    B = np.nanmedian(this_clob.lats) * this_clob.units['lats']
    
    r_out = np.zeros((len(possible_Rs), len(vels)))
    
    for r, R in enumerate (possible_Rs):
        this_R = R
        row = np.zeros(len(vels))
        
        for v, V in enumerate(vels):
            r_value = R_poly(R, L, B, V, Theta_sun, R_sun, a1, a2, a3)
            row[v] = r_value.value
            
        r_out[r] = row
        
    
    r_out_med = np.nanmedian(r_out, axis = 1)

    closest_R = possible_Rs[np.argmin(np.abs(0-r_out_med))]

    d = dist_yarn(closest_R, R_sun, L, B)
    
    my_Zucker_dists[n] = d.value
            

Zucker_ds = Zucker_df['d']


fig, axs = plt.subplots(1,2, figsize = (15, 5))

axs[0].vlines(names, my_Zucker_dists/1000, Zucker_ds/1000, color = 'orange', alpha = 0.5)
axs[0].vlines(names, -2, Zucker_ds/1000, color = 'orange', linestyle = '--', alpha = 0.3)

axs[0].scatter(names, my_Zucker_dists/1000, marker = '+', color = 'black', label = r'kinematic distance calculated from data')
axs[0].scatter(names, Zucker_ds/1000 , marker = '+', color = 'green', label = r'distance from Zucker et al')
axs[0].scatter(names, my_Zucker_dists_near/1000 , marker = '+', color = 'red', label = r'near kinematic distance')


axs[0].set_xticks(names, labels = names, rotation = 'vertical')

axs[0].set_ylabel(r'heliocentric distance [kpc]')

# axs[0].set_yscale('log')
axs[0].legend()

# # ------------------
# res = 100* np.abs(As3D - As2D)/As3D

# axs[1].vlines(names, 0, res, color = 'orange', linestyle = '--', alpha = 0.3)
# axs[1].hlines(10, names[0], names[24], color = 'purple', linestyle = '--', alpha = 0.3, label = '10%')
# axs[1].hlines(50, names[0], names[24], color = 'purple', linestyle = '-.', alpha = 0.3, label = '50%')


# axs[1].scatter(names, res , marker = '+', color = 'black')

# axs[1].set_xticks(names, labels = names, rotation = 'vertical')

# axs[1].set_ylabel(r'percentage difference')

# # axs[1].set_yscale('log')

# # axs[1].legend()

# plt.show()


# # plt.scatter(names, dists/Zucker_ds, marker = '+', color = 'green', \
# #             label = r'$\frac{kinematic}{published}$')

# plt.scatter(names, 100*(dists - Zucker_ds)/dists, marker = '+', color = 'green', \
#             label = r'percent variation')

# plt.xticks(rotation = 'vertical')
# # plt.ylabel('heliocentric distance [pc]')
# plt.ylabel('percent variation from kinematic distance')

# plt.legend()
# # plt.yscale('log')
axs[0].set_ylim(-2, 15)

plt.show()


#one way to figure out the distances of GMC-less HI clouds would be to find GMC-ful HI clouds
#that are peaking at the same velocity range & say they are located at the same kinematic distance
#and note whenever kinematic distance ambiguities come up

#but also, keep looking into the discrepancies between kinematic & parallax distance & see if there
#is an explanation
#read this in full: https://arxiv.org/pdf/2006.13654.pdf
#go back and read Zucker more closely to see if they comment on this discrepancy as well

#for deciding what the acceptable error bars on the distance are, look at 
#https://arxiv.org/pdf/2309.04532.pdf 


#^^ what this plot says is that, whenever you do the subcube first then calculate the std, 
#there is a wide variation based on the l, b coordinates you are looking at
#however, if you calculate the std of the entire cube, it is pretty similar to the published values
#so can use that 

#also for sigma thresholding -> will need to do two separate ones? ie one where you take 0th moment cube, then
#pick the values that are above the sigma_cube * dv in the 0th moment, vs the other where you make a 
#mask over the cube that keeps everything above 3sigma
#for mass calculation, either method will work (and doing the thresholding on 2D retains more pixels than
#thresholding on the 3D) .... but, if we want to be selecting a velocity range that corresponds to 
#the cloud, then we want to be selecting in the 3D (once you collapse then make a mask, its unclear how you
# can then apply that mask back onto the cube in order to select velocity channels)
#and you cant really use the subcube's sigma at that point bc that only works to reduce noise if you already
#know that most of the pixels in the cube contain noise (true for the big cube but not for the spatially
#isolated ones) otherwise a 3sigma threshold is well above where it should be 

#talk to nia ab this conundrum


#also, the kinematic distance i calculated using rotation curve is 
#way different than the published results to the GMCs (2x for Orion, up to 100x 
#for aquila, chamaeleon, maddalena, etc)
#(double check kinematic distance calculation using 
#trey venger's code https://github.com/tvwenger/kd )


len(clouds)


len(gaia_ds)


len(kin_ds)


#look at table 1 in yan et al (https://arxiv.org/pdf/2006.13654.pdf)
#- there seems to be a big discrepancy between 
#the gaia (parallax) & kinematic distances

clouds = ['G026.9−03.5', 'G027.8−02.1','G028.8−01.9', 'G029.6+03.7','G030.3+01.2', \
         'G030.6−03.6', 'G035.0+00.6','G036.6−00.9','G037.0−01.4', 'G037.5+04.0', \
         'G038.3−00.1', 'G039.4−02.7', 'G040.7−04.1','G041.3+03.2', 'G041.5+02.3',\
         'G042.0−00.9', 'G043.2+02.1', 'G043.3+03.1', 'G044.5+02.6', 'G044.8+04.0',\
         'G045.5−04.3', 'G046.2+03.1', 'G046.2−01.6', 'G046.5−03.3', 'G046.9−02.9', \
         'G048.3−01.7',' G048.9+02.3', 'G049.0+03.4']

gaia_ds = np.array((425, 487, 261, 547, 428, 317, 606, 588, 528, 515, 645, \
          492, 377, 1071, 905, 719, 872, 731, 772, 859, 938, 937, \
          394, 1348, 1300, 1134, 483, 480)) * u.pc

kin_ds = np.array((1.16, 1.26, 0.21, 0.49, 0.32, 0.62, 0.84, 0.83, 1.2, 1.01,\
         1.05, 0.96, 0.49, 1.73, 1.15, 0.69, 1.51, 0.67, 0.97, 1.33,\
         1.24, 1.80, 0.47, 1.84, 1.78, 1.47, 0.47, 0.44)) * 1000*u.pc


# plt.scatter(clouds, gaia_ds, marker = '+', label = 'gaia', color = 'pink')
# plt.scatter(clouds, kin_ds, marker = '+', label = 'kinematic', color = 'black')

plt.scatter(clouds, 100*(kin_ds - gaia_ds)/kin_ds, marker = '+', label = 'percent variation', color = 'green')


# plt.ylabel('distance [pc]')
plt.ylabel('percent variation from kinematic distance')
plt.xlabel('clouds')
plt.xticks(rotation='vertical')
plt.show()


MD_catalog = pd.read_csv('MD2017_catalog.dat', delim_whitespace=True)


MD_catalog['Rgal']


MD_Ls = MD_catalog['GLON']
MD_Bs = MD_catalog['GLAT']
MD_vlsrs = MD_catalog['Vcent']


len(MD_Bs), len(MD_Ls), len(MD_vlsrs)


#comparing the kinematic distances i calculate from l, b, v with the ds reported by MD
Rs = np.zeros((len(MD_Bs), len(possible_Rs)))
              
for i in range(len(MD_Bs)):
    l = MD_Ls[i]*u.deg
    b = MD_Bs[i]*u.deg
    v = MD_vlsrs[i]*u.km/u.s
              
    row = np.zeros(len(possible_Rs))
    for r, R in enumerate(possible_Rs):
        y = R_poly(R, l, b, v, Theta_sun, R_sun, a1, a2, a3)
        row[r] = y.value
              
    Rs[i] = row
    
    
R_solutions = np.zeros(len(Rs))
for i, arr in enumerate(Rs):
    zeroth = np.argmin(np.abs(0-arr))
    R_solutions[i] = possible_Rs[zeroth].value
    
my_MD_ds = np.zeros_like(R_solutions)  

for i, ith in enumerate(R_solutions):
    l = MD_Ls[i]*u.deg
    b = MD_Bs[i]*u.deg
    
    R_gal = ith*u.pc
    
    d = dist_yarn(R_gal, R_sun, l, b)
    
    my_MD_ds[i] = d.value
    

Rs.shape


dist = dist_yarn(closest_R, R_sun, L, B)


Rs.shape


my_MD_ds_kpc = my_MD_ds/1000


MD_nums = np.arange(0, len(MD_catalog))
plt.scatter(MD_nums, 100*(np.abs(my_MD_ds_kpc - MD_catalog['Dfar']))/my_MD_ds_kpc, s = 0.05)
# plt.scatter(MD_nums, 100*(np.abs(my_MD_ds/1000 - MD_catalog['Dnear']))/my_MD_ds_kpc, s = 0.05, alpha = 0.5)


plt.ylabel('percent difference in distance')


# plt.yscale('log')
plt.ylim(-1,20)

plt.show()


MD_nums


MD_catalog['Dfar']


MD_ds


len(np.isnan(MD_ds)==True)


for i in range(len(MD_ds)):
    print('far', MD_catalog['Dfar'][i], 'near', MD_catalog['Dnear'][i], 'mine', MD_ds[i] )


len(R_solutions)


MD_ds = dists


# since it seems like there's a wide range of R values, maybe the thing to do is some sort of adaptive search
# algorithm that does a pretty coarse array of possible Rs, then iterates n number of times, each time finding the 
# R that minimizes the polynomial and looking prev_interval/x around it to find the direction
# of decrease, and repeat (think leapfrog + binary search (the one where you keep dividing your
#                                                         total set in half))


for i in range(0,10):
    plt.plot(Rs[i])
    plt.show()


plt.plot(r_out_med)


r_value = R_poly(R, L, B, V, Theta_sun, R_sun, a1, a2, a3)


np.zeros((3, 10))


for n in range(len(Zucker_df)):
    print(Zucker_df['Cloud'][n])
    b0 = np.array(Zucker_df['b0'][n])
    b1 = np.array(Zucker_df['b1'][n])


    l0 = np.array(Zucker_df['l0'][n])
    l1 = np.array(Zucker_df['l1'][n])


    DHT_CO.create_subcube([b0, b1, l0, l1], 'CO')


#from this it doesn't look like theres a big (enough) difference in the resultig sub-cube, based on whether the threshold is applied on the cube
#or on the 0th moment (integrated intensity) map ==> ie using the noise

for n in range(len(Zucker_df)):
    print(Zucker_df['Cloud'][n])
    b0 = np.array(Zucker_df['b0'][n])
    b1 = np.array(Zucker_df['b1'][n])


    l0 = np.array(Zucker_df['l0'][n])
    l1 = np.array(Zucker_df['l1'][n])


    DHT_CO.create_subcube([b0, b1, l0, l1], 'CO')


0.22*np.sqrt(DHT_CO.dv)


0.22*DHT_CO.dv


for n in (15, 17, 20, 22):
    print(Zucker_df['Cloud'][n])
    b0 = np.array(Zucker_df['b0'][n])
    b1 = np.array(Zucker_df['b1'][n])


    l0 = np.array(Zucker_df['l0'][n])
    l1 = np.array(Zucker_df['l1'][n])
    
    cld.get_parent_HI(n, Zucker_df, DHT_CO, HI4PI, show_0th=False, show_spec=False)
    print('\n')
#     co_clob = DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'a', 'CO')
#     sig, a = cld.surface_density(co_clob, 'CO', show_Sigma=False).values()
# #     print(cld.calculate_mass(n, Zucker_df, sig, a))
#     HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'HI')
    
    
#there s a lot of math that needs to be checked because comparing these values with nias ones in
#tables 1 & 2 of this paper:https://iopscience.iop.org/article/10.1088/0004-637X/732/2/78/pdf ... 
#the inconsistencies between what she has and what i get are not systematic ....
#something is wrong but idk what


# my values 
#                   Zucker num       d            M          R_MC       R_A
# Perseus             20           294pc      2.3e4 Msun    10.9pc
# Orion A             17           432pc      6.6e4 Msun    16.3pc
# Rosette             22           1304pc     7.5e4 Msun    27.4pc
# MonR2               15           778pc      7.3e4 Msun    23.4pc
# 


Zucker_df['d'][15], Zucker_df['d'][17], Zucker_df['d'][20], Zucker_df['d'][22]


figuring_out_values['my d'] = [294, 432, 1304, 778]
figuring_out_values["Nia's d"] = [320, 414, 1600, 830]


figuring_out_values["my M"] = [1.7, 3.2, 7.2, 4.3]#[2.3, 6.6, 7.5, 7.3]
figuring_out_values["Nia's M"] = [1.2, 6.9, 7.7, 9.0]


figuring_out_values["my R_MC"] = [8.5, 10.4, 26.7, 17.2]#[10.9, 16.3, 27.4, 23.4]
figuring_out_values["Nia's R_MC"] = [16.2, 21.9, 26.5, 30.6]


figuring_out_values["my R_A"] = [67.1, 94.6, 98.5, 97.9]
figuring_out_values["Nia's R_A"] = [42.4, 76.0, 78.8, 83.6]


new_Zucker = pd.DataFrame(columns=Zucker_df.columns)


new_Zucker['d'] = figuring_out_values["Nia's d"]
new_Zucker['Cloud'] = figuring_out_values["Cloud"]


new_Zucker['d_unit'] = ['pc', 'pc', 'pc', 'pc']
new_Zucker['lb_unit'] = ['deg', 'deg', 'deg', 'deg']


new_Zucker['l0'] = [Zucker_df['l0'][20], Zucker_df['l0'][17], Zucker_df['l0'][22], Zucker_df['l0'][15]]
new_Zucker['l1'] = [Zucker_df['l1'][20], Zucker_df['l1'][17], Zucker_df['l1'][22], Zucker_df['l1'][15]]
new_Zucker['b0'] = [Zucker_df['b0'][20], Zucker_df['b0'][17], Zucker_df['b0'][22], Zucker_df['b0'][15]]
new_Zucker['b1'] = [Zucker_df['b1'][20], Zucker_df['b1'][17], Zucker_df['b1'][22], Zucker_df['b1'][15]]


for n in range(len(new_Zucker)):
    print(new_Zucker['Cloud'][n])
    b0_phy = np.array(new_Zucker['b0'][n])
    b1_phy = np.array(new_Zucker['b1'][n])


    l0_phy = np.array(new_Zucker['l0'][n])
    l1_phy = np.array(new_Zucker['l1'][n])
    
    cld.get_parent_HI(n, new_Zucker, DHT_CO, HI4PI, show_0th=False, show_spec=False)
    print('\n')
#     co_clob = DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'a', 'CO')
#     sig, a = cld.surface_density(co_clob, 'CO', show_Sigma=False).values()
# #     print(cld.calculate_mass(n, Zucker_df, sig, a))
#     HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'HI')
    
    
#there s a lot of math that needs to be checked because comparing these values with nias ones in
#tables 1 & 2 of this paper:https://iopscience.iop.org/article/10.1088/0004-637X/732/2/78/pdf ... 
#the inconsistencies between what she has and what i get are not systematic ....
#something is wrong but idk what


#from new Zucker, ie all my same data except using the d values from Nia's paper
recalculated_Ms = [2.8, 6.1, 11.3, 8.5]
recalculated_RMCs = [11.9, 15.6, 33.7, 25.0]
recalculated_RAs = [71.0, 92.0, 112.9, 102.2]


rho = M/V
sigma = M/A

rho = sigma/R


def calc_RA(M, rho):
    RA = (M/(2*np.pi*rho))**(1/3)
    return RA


((M_pers)/(2*np.pi*(42.4*u.pc)**(3))).to(c.M_sun/u.pc**3)


rho_HI = 0.025*c.M_sun/(u.pc**3)


#redoing the R_A calculation, this time using the mass values from Nia's paper
M_pers = 1.2e4*c.M_sun
M_ora = 6.9e4*c.M_sun
M_ros = 7.7e4*c.M_sun
M_mon = 9.0e4*c.M_sun


calc_RA([M_pers.value, M_ora.value, M_ros.value, M_mon.value], rho_HI.value)


figuring_out_values


rerecalculated_RAs = [ 53.5, 95.8, 99.3, 104.6]


#re-redoing the R_A calculation, this time using the 13CO (before was 12CO) values from Nia's paper
Ms_from13CO = np.array((1.4, 2.4, 3.0, 4.0))*1e4*c.M_sun.value


calculated_from_13COM = calc_RA(Ms_from13CO, rho_HI.value)


#calculating distances kinematically to compare to literature values
def R_poly (R, l, b, v_lsr, theta_sun, r_sun, a1, a2, a3):

    term1 = R/r_sun
    term2 = (R/r_sun)**a2
    
    coeff1 = (v_lsr/np.cos(b))+ (theta_sun*np.sin(l))
    coeff2 = a1*theta_sun*np.sin(l)
    coeff3 = a3*theta_sun*np.sin(l)
             
    final = coeff1*term1 - coeff2*term2 - coeff3
    
    return final


# def dist(R, r_sun, l, b):
#     print(R, r_sun)
#     term1 = r_sun*np.cos(l)/np.cos(b)
    
#     term2_1 = R**2
#     term2_2 = (r_sun*np.sin(l))**2
#     term2 = np.sqrt(term2_1 - term2_2)/np.cos(b)
    
    
#     final = term1 + term2
#     if final

# #     if R > r_sun:
# #         final = term1 + term2
    
# #     else:
# #         print('in else')
# #         final = (term1 + term2, term1 - term2)
    
#     return final


def dist_yarn(R, r_sun, l, b): #this and dist give pretty similar results, <0.2% differences
    term1 = r_sun*np.cos(b)*np.cos(l)
    
    term2_1 = term1**2
    term2_2 = (r_sun**2)
    term2_3 = R**2
    
    term2 = np.sqrt(term2_1 - term2_2 + term2_3)
    
    final = term1 + term2
    
    return final


# #values from Russeil et al 2017
# Theta_sun = 240*u.km/u.s
# R_sun = 8.34*1000*u.pc
# a1 = 1.022
# a2 = 0.0803
# a3 = 0.


#values from Brand & Blitz 1993
Theta_sun = 220*u.km/u.s
R_sun = 8.5*1000*u.pc
a1 = 1.00767
a2 = 0.0394
a3 = 0.00712


possible_Rs = np.arange(1.7, 17, 0.01)*1000*u.pc


possible_Rs


b0_pers = np.array(Zucker_df['b0'][0])
b1_pers = np.array(Zucker_df['b1'][0])


l0_pers = np.array(Zucker_df['l0'][0])
l1_pers = np.array(Zucker_df['l1'][0])


pers_clob = DHT_CO.create_subcube([b0_pers, b1_pers, l0_pers, l1_pers], 'CO')


#Dame = 0.25


b0_pers = np.array(new_Zucker['b0'][0])
b1_pers = np.array(new_Zucker['b1'][0])


l0_pers = np.array(new_Zucker['l0'][0])
l1_pers = np.array(new_Zucker['l1'][0])


pers_clob = DHT_CO.create_subcube([b0_pers, b1_pers, l0_pers, l1_pers], 'CO')


#------------------------------------------------------------
pers_vels = pers_clob.vels * pers_clob.units['vels']
persL = np.nanmedian(pers_clob.lons)
persB = np.nanmedian(pers_clob.lats)

r_out = np.zeros((len(possible_Rs),len(pers_vels)))

for r, R in enumerate(possible_Rs):
    this_R = R
    row = np.zeros(len(pers_vels))
    
    
    for v, V in enumerate(pers_vels):
        r_value = R_poly(R, persL, persB, V, Theta_sun, R_sun, a1, a2, a3)
        row[v] = r_value.value
        
#         print(row[v])
#         row[v] = (R_poly(R, persL, persB, V, Theta_sun, R_sun, a1, a2, a3)).value()
        
    r_out[r] = row


#------------------------------------------------------------
r_out_med = np.nanmedian(r_out, axis =1)

closest_R = possible_Rs[np.argmin(np.abs(0-r_out_med))]

dist = dist_yarn(closest_R, R_sun, persL, persB)

print(dist)


plt.imshow(pers_clob.calculate_0th_moment(), **kw_imshow)


for i in range(len(new_Zucker)):
    print(i)


kinematic_dists = np.zeros(len(new_Zucker))

for i in range(len(new_Zucker)):

    b0 = np.array(new_Zucker['b0'][i])
    b1 = np.array(new_Zucker['b1'][i])

    l0 = np.array(new_Zucker['l0'][i])
    l1 = np.array(new_Zucker['l1'][i])


    clob = DHT_CO.create_subcube([b0, b1, l0, l1], 'CO')


    #------------------------------------------------------------
    vels = clob.vels * clob.units['vels']
    L = np.nanmedian(clob.lons)
    B = np.nanmedian(clob.lats)

    r_out = np.zeros((len(possible_Rs),len(vels)))

    for r, R in enumerate(possible_Rs):
        this_R = R
        row = np.zeros(len(vels))


        for v, V in enumerate(vels):
            r_value = R_poly(R, L, B, V, Theta_sun, R_sun, a1, a2, a3)
            row[v] = r_value.value

    #         print(row[v])
    #         row[v] = (R_poly(R, persL, persB, V, Theta_sun, R_sun, a1, a2, a3)).value()

        r_out[r] = row


    #------------------------------------------------------------
    r_out_med = np.nanmedian(r_out, axis =1)

    closest_R = possible_Rs[np.argmin(np.abs(0-r_out_med))]

#     dist = 
    d = dist_yarn(closest_R, R_sun, L, B)
#     print(dist(closest_R, R_sun, L, B))
    kinematic_dists[i] = d.value

#     print(dist)


kinematic_dists


for n in range(len(new_Zucker)):
    print(new_Zucker['Cloud'][n])
    b0_phy = np.array(new_Zucker['b0'][n])
    b1_phy = np.array(new_Zucker['b1'][n])


    l0_phy = np.array(new_Zucker['l0'][n])
    l1_phy = np.array(new_Zucker['l1'][n])
    
    cld.get_parent_HI(n, new_Zucker, DHT_CO, HI4PI, show_0th=False, show_spec=False)
    print('\n')
#     co_clob = DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'a', 'CO')
#     sig, a = cld.surface_density(co_clob, 'CO', show_Sigma=False).values()
# #     print(cld.calculate_mass(n, Zucker_df, sig, a))
#     HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'HI')
    
    
#there s a lot of math that needs to be checked because comparing these values with nias ones in
#tables 1 & 2 of this paper:https://iopscience.iop.org/article/10.1088/0004-637X/732/2/78/pdf ... 
#the inconsistencies between what she has and what i get are not systematic ....
#something is wrong but idk what


fig, axs = plt.subplots(2, 2, figsize = (15, 15))

axs[0,0].scatter(figuring_out_values['Cloud'], figuring_out_values['my d'], \
                 marker = '+', label = 'my values (from Zucker)')
axs[0,0].scatter(figuring_out_values['Cloud'], figuring_out_values["Nia's d"], \
                 marker = 'x', label = "Nia's values")
# axs[0,0].scatter(new_Zucker['Cloud'], kinematic_dists, \
#                  marker = '+', label = "kinematic distances calculated from data", color = 'xkcd:fuchsia')

axs[0,0].set_ylabel('distance to GMC [pc]')
# axs[0,0].set_yscale('log')
axs[0,0].legend()


#--------------------------------------------------------------------------------------------------------------
axs[0,1].scatter(figuring_out_values['Cloud'], figuring_out_values['my M'], \
                 marker = '+', label = 'my values', s = 40)
axs[0,1].scatter(figuring_out_values['Cloud'], figuring_out_values["Nia's M"], \
                 marker = 'x', label = "Nia's values", s= 40)
axs[0,1].scatter(new_Zucker['Cloud'], recalculated_Ms, \
                 marker = '+', label = "recalculated using Nia's d", s = 40)

axs[0,1].set_ylabel(r'mass of GMC [$10^4$$M_\odot$]')
axs[0,1].legend()

#--------------------------------------------------------------------------------------------------------------

axs[1,0].scatter(figuring_out_values['Cloud'], figuring_out_values['my R_MC'], \
                 marker = '+', label = 'my values')
axs[1,0].scatter(figuring_out_values['Cloud'], figuring_out_values["Nia's R_MC"], \
                 marker = 'x', label = "Nia's values")
axs[1,0].scatter(new_Zucker['Cloud'], recalculated_RMCs, \
                 marker = '+', label = "recalculated using Nia's d")

axs[1,0].set_ylabel('effective radius of GMC [pc]')
axs[1,0].legend()


#--------------------------------------------------------------------------------------------------------------
# axs[1,1].scatter(figuring_out_values['Cloud'], figuring_out_values['my R_A'], \
#                  marker = '+', label = 'my values')
axs[1,1].scatter(figuring_out_values['Cloud'], figuring_out_values["Nia's R_A"], \
                 marker = 'x', label = "Nia's values")
# axs[1,1].scatter(new_Zucker['Cloud'], recalculated_RAs, \
#                  marker = '+', label = "recalculated using Nia's d")
axs[1,1].scatter(new_Zucker['Cloud'], rerecalculated_RAs, \
                 marker = '+', label = r"recalculated using Nia's $M_{GMC}$")
# axs[1,1].scatter(new_Zucker['Cloud'], calculated_from_13COM, \
#                  marker = '+', label = r"calculated from $^{13}CO$ mass")


#this plot is confusing => i used eq 12 in nia's paper 1, but even when using the same M values as listed
#in the paper, i am getting a different value of R_A
#can't tell if it's because the value of rho i should be using is diff from 0.125Msun/pc**3
#or whether i need to be using a different formula


axs[1,1].set_ylabel('accumulation radius [pc]')
axs[1,1].legend()


plt.show()


#look in the literature for different boundaries people use to define the clouds
#catherine zucker, charlie lada


#the nia way does rms -> calculate the sigma of the cube, 
#then collapse the cube into a surface density map, then apply that rms
#and the relationship between rms of the cube and rms of the 2D (0th moment?) 
#is something to do with square root of velocity dispersion


#a careful way to calculate the rms is to go to a 
#region of the sky where you expect there to not be signal
#then calculate the rms value there & use it as an estimation of sigma_rms


Zucker_df


new_Zucker


# nia's values
#              d         M            R_MC        R_A
# Perseus    320pc     1.2e4 Msun    16.2pc      42.4pc
# Orion A    414pc     6.9e4 Msun    21.9pc      76.0pc
# Rosette    1600pc    7.7e4 Msun    26.5pc      78.8pc
# MonR2      830pc     9.0e4 Msun    30.6pc      83.6pc


b_cs = np.nanmean((new_Zucker['b0'], new_Zucker['b1']), axis = 0)
nias_ms = np.array(figuring_out_values["Nia's M"])*1e4*c.M_sun
nias_ds = np.array(figuring_out_values["Nia's d"])*u.pc

cld.calculate_accumulation_bounds(nias_ms,b_cs, nias_ds)


for n in range(len(Zucker_df)):
    print(Zucker_df['Cloud'][n])
    b0_phy = np.array(Zucker_df['b0'][n])
    b1_phy = np.array(Zucker_df['b1'][n])


    l0_phy = np.array(Zucker_df['l0'][n])
    l1_phy = np.array(Zucker_df['l1'][n])
    
    cld.get_parent_HI(n, Zucker_df, DHT_CO, HI4PI, show_0th=False, show_spec=False)
    print('\n')
#     co_clob = DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'a', 'CO')
#     sig, a = cld.surface_density(co_clob, 'CO', show_Sigma=False).values()
# #     print(cld.calculate_mass(n, Zucker_df, sig, a))
#     HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'HI')
    

for n in range(len(Zucker_df)):
    print(Zucker_df['Cloud'][n])
    b0_phy = np.array(Zucker_df['b0'][n])
    b1_phy = np.array(Zucker_df['b1'][n])


    l0_phy = np.array(Zucker_df['l0'][n])
    l1_phy = np.array(Zucker_df['l1'][n])
    
    cld.get_parent_HI(n, Zucker_df, DHT_CO, HI4PI, show_0th=False, show_spec=False)
    print('\n')
#     co_clob = DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'a', 'CO')
#     sig, a = cld.surface_density(co_clob, 'CO', show_Sigma=False).values()
# #     print(cld.calculate_mass(n, Zucker_df, sig, a))
#     HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'HI')
    

np.mean((3, 4))


for n in range(len(Zucker_df)):
    print(n)
    b0_phy = np.array(Zucker_df['b0'][n])
    b1_phy = np.array(Zucker_df['b1'][n])


    l0_phy = np.array(Zucker_df['l0'][n])
    l1_phy = np.array(Zucker_df['l1'][n])
    
    DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'CO')
    HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p', 'HI')
    

Zucker_df


orA_spec=dc.spectrum1D(oriona_CO.cube)


orA_spec


rms(dc.spectrum1D(oriona_CO.cube))


np.nanstd(dc.spectrum1D(oriona_CO.cube))


whole_orion_cube = dc.unpack_fits_cube('DHT27_Orion_raw.fits', 'DHT')


whole_orion_cube.cube.shape


48*57 + 32*(73-57) + 16*(121-73)


153*209 - (48*57 + 32*(73-57) + 16*(121-73))


2*27961


(77*153 + 77*209 + 153*209)-(48*57 + 32*(73-57) + 16*(121-73))


print(whole_orion_cube.lats[0], whole_orion_cube.lats[-1])
print(whole_orion_cube.lons[0], whole_orion_cube.lons[-1])


-3.125-(-22.125)


226.125-200.125


19*26


plt.imshow(np.nanmean(whole_orion_cube.cube, axis=0), **kw_imshow)
plt.hlines(121, 0, 210)
plt.hlines(73, 0, 210)
plt.hlines(57, 0, 210)

plt.vlines(16, 0, 145)
plt.vlines(32, 0, 145)
plt.vlines(48, 0, 145)
plt.vlines(208, 0, 145)


77*153*209


2,462,229


from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import snr


snr(orA_spec)


oriona_CO = DHT_CO.create_subcube([-20.5, -18.5, 207.0, 215.0], 'p')


Zucker_df


oriona_CO = DHT_CO.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p')


oriona_HI = HI4PI.create_subcube([b0_phy, b1_phy, l0_phy, l1_phy], 'p')


from scipy.signal import peak_widths, find_peaks

def surface_density (clob, data_type, W_thresh = 1.5, v_channel_width = 25, show_Sigma = True):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from extract_cloud/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    and whether to assume that the spectrum is noisy / will need to isolate velocity peak or not (this is mostly
    for HI clouds)
    3) W_thresh (using value from lewis, lada, dame, but think about whether i wanna change it)

    '''
    spec = dc.spectrum1D(clob.cube)
    plt.plot(spec)
    plt.show()
    
    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/2203.07480.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        thresh = np.nanmin(spec)
        

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        thresh = np.nanmean(spec) + 3 * np.nanstd(spec)
    
        
    p = find_peaks (spec, height = thresh) 
    #something counts as a peak only if its above the mean-3sigma value 
    #this is because the assumption here is that the majority of the velocity channels
    #are noisy (ie a smaller fraction of the channels have data vs noise. this means that 
    #the mean of the values across the spectral axis is telling you about noise, and 
    #any value that is within the 3sigma of the noise mean is also noise. ie what counts as
    #a velocity channel containing signal is one whose value is > mean + 3sigma)
    #but then we run into a problem when every pixel contains data / the noise has already been clipped away
    print(p)
    peak = p[0] #array of peaks locations
    peak_info = p[1] #array of peak heights
    
    if len(peak) > 1:
        i = np.argmax(peak_info['peak_heights']) #if there are multiple peaks, pick the tallest one #this is problematic though if there is blending
        peak = [peak[i]]
    vel_center = peak[0]

    #this assumes that vel_center is at a place along the spectral axis that there is
    # still data +/- v_channel_width/2 around it ... seems like an ok assumption for now but might 
    # run into problems at some point
    v_channel_start =  vel_center - int(v_channel_width/2)
    v_channel_end = vel_center + int(v_channel_width/2)

    subcube = (clob.cube[v_channel_start:v_channel_end+1, :, :])*clob.units['data']

    #calculate integrated intensity
    W = np.nansum(subcube*clob.dv, axis=0)
    # W = np.nansum(cube*dv, axis=0)

    
    #use W_thresh to create a mask:
    W_mask = W.value > W_thresh
    W_old = W
    W = W_old * W_mask
    
    #calculate column density from the integrated intensity & conversion factor
    N = W * conversion_factor
    
    #define surface density here
    dSig = mmw * (c.m_p.to(u.g)) * N
    
    
    #calculate on sky area & cloud area at each pixel
    dAsky = np.abs(clob.dl * clob.db)
    

    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        kw_imshow = {'cmap': cmap0, 'origin':'lower'}

    
        dSig_Msunpc2 = dSig.to(c.M_sun/u.pc**2)

        fig, ax = plt.subplots (figsize = (15, 10), dpi = 200)
        Sig_plot = ax.imshow(dSig_Msunpc2.value, extent = [clob.lons[0], clob.lons[-1], \
                                                clob.lats[0], clob.lats[-1]], **kw_imshow)


        cbar = fig.colorbar(Sig_plot, ax=ax, label = Sig_label +r' [$M_\odot \ pc^{-2}$]', orientation = 'vertical')

        if data_type == '21cm':
            KMT_thresh = plt.contour(clob.lons, clob.lats, dSig_Msunpc2.value,  linewidths = 0.7, \
            levels=[10, 100], colors = ['xkcd:red', 'xkcd:orange']) #atomic-to-molecular transition for HI maps 
            cbar.add_lines(KMT_thresh)


        ax.set_xlabel(r'Galactic Longitude [$\degree$]')
        ax.set_ylabel(r'Galactic Latitude [$\degree$]')
        ax.set_title(Sig_label)


        plt.show()

    
    return {'surface density per pixel': dSig, 'pixel area': dAsky}


def mass_and_area (cloud_num, database, surface_density_map, dA_perpx):
    """
    the first two help with retrieving the distance to the cloud (this is)
    """
    x, y = surface_density_map.shape
    num_pix = x*y

    d = database['d'][cloud_num] * u.Unit(database['d_unit'][cloud_num])
    A_rad = num_pix*dA_perpx.value*((area.unit).to(u.rad**2))
    
    A_total = d**2 * A_rad #total area, this assumes that all the pixels in surface_density_map contain mass we care about
    M_total = np.nansum(surface_density_map) * A_total
    M_total_kg = M_total.to(u.kg)


    return {'M':M_total_kg, 'A':A_total}


orion


mass_and_area(17, Zucker_df, Sig_oriona, area)


x, y = Sig_oriona.shape


x


surface_density(oriona_CO, 'CO', show_Sigma=True).values()


surface_density(oriona_HI, '21cm', show_Sigma=True).values()


surface_density(oriona_HI, '21cm', show_Sigma=True).values()


surface_density(oriona_CO, 'CO', show_Sigma=True).values()


### realizing what a 3sigma cutoff means !!! ie that you are finding the mean, 
#then setting the lower limt as 3*sigma below the mean
#should change how i use signal to extract the spectrum / make masks in the first place
#go look at create_subcube


Sig_oriona.shape


(Sig_oriona.flat())


cld.surface_density(oriona_CO, 1, 'CO', show_Sigma=True)


sig, a = (cld.surface_density(oriona_CO, 1, 'CO', show_Sigma=False)).values()


sig.shape


plt.imshow(sig.value, **kw_imshow)
# plt.colorbar()


plt.imshow(oriona_CO.calculate_0th_moment(), **kw_imshow)
# plt.colorbar()


plt.imshow(oriona_HI.calculate_0th_moment())


cloud_num, database, \
larger_cloud, p_padding = 0, v_padding = 0, v_guess = None, sigma = 3


cloud, parent = cld.get_parent_HI(1, Zucker_df, DHT_CO, HI4PI, \
                      CO_padding=20, HI_padding=40,  v_padding=10, \
                      show_0th = False, show_spec = False, sigma = 3).values()


#surface density maps for solar neighborhood clouds + their HI envelopes
for i in range(len(Zucker_df)):
    cloud, parent = cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, \
                      CO_padding=20, HI_padding=40,  v_padding=10, \
                      show_0th = False, show_spec = False, sigma = 3).values()


#0th moment maps & velcoity spectra for solar neighborhood clouds + their HI envelopes

for i in range(len(Zucker_df)):
    cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, \
                      CO_padding=20, HI_padding=40,  v_padding=10, \
                      show_0th = True, show_spec = True, sigma = 3)


# (HI_cube.size * HI_cube.itemsize)/1e9 #the size of HI_cube in gigs = 35.3


HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 
# ^^ this takes too long to run, so instead of making a 0th moment map of the whole
# thing, i will just cut out the regions around the solar neighborhood clouds


# ### GALAXY 0TH MOMENT (in HI) ####
# plt.figure(figsize=(15,3), dpi=200)

# im = HI_0th
# imx = HI_lons
# imy = HI_lats


# vmin,vmax = zlims(im) 


# plt.imshow(im, vmin=vmin, vmax=vmax, \
#             extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


# plt.ylabel('Galactic Latitude')
# plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')


# plt.title('Galaxy in HI: 0th moment')
# plt.show()


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
