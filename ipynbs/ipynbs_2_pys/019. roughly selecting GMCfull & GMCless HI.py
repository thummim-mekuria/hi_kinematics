#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


from astropy.utils import data
from spectral_cube import SpectralCube

import datacube as dc
import clouds as cld

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c

from astropy.visualization import ZScaleInterval as zscale


from transform import transform_1D
from transform import transform_2D

from scipy.signal import peak_widths, find_peaks


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}
kw_colorb = {'label':r'$W_{CO}$ (K*km/s)', 'orientation':'horizontal'}


#-----------------------------------------------------------------------------------
def rms(thing):
    return np.sqrt((np.nanmean(thing**2)))
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
#all-sky map from HI4PI
LAB = dc.unpack_fits_cube('labh.fits', 'LAB')
LAB_cube, LAB_vels, LAB_lats, LAB_lons, LAB_units, LAB_head = LAB.values()

LAB_dv = ((LAB_vels[1] - LAB_vels[0])*LAB_units['vels']).to(u.km/u.s)
LAB_db = (LAB_lats[1] - LAB_lats[0])*LAB_units['lats']
LAB_dl = (LAB_lons[1] - LAB_lons[0])*LAB_units['lons']
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
#DHT CO survey
DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
CO_cube, CO_vels, CO_lats, CO_lons, CO_units, CO_head = DHT_CO.values()

#assume arrays have linear spacing
CO_dv = ((CO_vels[1] - CO_vels[0])*CO_units['vels']).to(u.km/u.s)
CO_db = (CO_lats[1] - CO_lats[0])*CO_units['lats']
CO_dl = (CO_lons[1] - CO_lons[0])*CO_units['lons']

CO_0th = np.nansum(CO_cube*CO_dv.value, axis=0)
#-----------------------------------------------------------------------------------


def calculate_subcube_linewidth (b0, l0, s):
    subcube = LAB_cube[:, b0:b0+s, l0:l0+s]
    
    spectrum = subcube.mean(axis=(1,2))
    
    linewidth = dc.dispersion_1sigma(spectrum) * LAB_dv
    
    return linewidth


starting_bs = np.arange(0,360,20)
starting_ls = np.arange(0,720,20)


linewidths = []
coords = []
for b in starting_bs:
    for l in starting_ls:
        sigma = calculate_subcube_linewidth(b, l, 20)
        linewidths.append(sigma)
        coords.append([b,l])
        
        
linewidths = np.array(linewidths)


plt.hist(linewidths, 100)
plt.vlines(15, 0, 40, linestyle='--', color = 'red', label = '92% below')
plt.vlines(20, 0, 40, linestyle='--', color = 'orange', label = '96% are below')
plt.vlines(25, 0, 40, linestyle='--', color = 'green', label = '98% are below')
plt.xlabel('Linewidths [km/s]')
plt.ylabel('Apertures')

plt.legend()
plt.show()


coords[301], coords[345], coords[346], coords[538], \
coords[539], coords[540], coords[575], coords[577], coords[580]


linewidths[301], linewidths[345], linewidths[346], linewidths[538], \
linewidths[539], linewidths[540], linewidths[575], linewidths[577], linewidths[580]


count_sub15 = 0
for i in linewidths:
    if i < 15:
        count_sub15 +=1


count_sub20 = 0
for i in linewidths:
    if i < 20:
        count_sub20 +=1


count_sub25 = 0
for i in linewidths:
    if i < 25:
        count_sub25 +=1


(count_sub15/648)*100, (count_sub20/648)*100, (count_sub25/648)*100


linewidths_2D = linewidths.reshape(18,36)

plt.figure(figsize=(15,3), dpi=200)

im = linewidths_2D
imx = LAB_lons
imy = LAB_lats


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'Linewidths [$km \ s^{-1}$]')

plt.title('Linewidths across the galaxy')
plt.show()


LAB_vels_kms = LAB_vels*(LAB_units['vels']).to(u.km/u.s)


def velocity_peaks(b0,l0, s):
    subcube = LAB_cube[:, b0:b0+s, l0:l0+s]
    
    spectrum = dc.spectrum1D(subcube)
    
    
    sigm = np.nanstd(spectrum)
    p = find_peaks (spectrum, height = 3* sigm) 
    #from visually inspecting the spectra, it looks like the 3*sigm line 
    #is still above a good number of bumps that maybe should count?? think more ab this
    
    peak = p[0] #array of peaks locations
    peak_info = p[1] #array of peak heights
        
    num_peaks = len(peak)
    
#     if num_peaks > 1:
#         print('peaks are at', peak)
#         dc.dispersion_1sigma(spectrum)
#         plt.plot(spectrum)
#         plt.hlines(3*sigm, 0,800, linestyle = '--')
#         plt.show()
    if len(peak) == 1:
        peak_locs = peak[0]
        
    else:
        peak_locs = int(np.nanmean(peak)) #not sure how to represent multiple, so averaging for now
#         print(int(np.nanmean(peak)))
#         print(peak_locs)
        
    peak_locs_kms = LAB_vels_kms[peak_locs]
        
    return num_peaks, peak_locs_kms


vel_peaks = []
for b in starting_bs:
    for l in starting_ls:
        num_peaks = velocity_peaks(b, l, 20)[0]
        vel_peaks.append(num_peaks)
        
        
vel_peaks = np.array(vel_peaks)

#okay so it looks like 'find peaks' thinks a peak is any local maximum,
#but i want to be thinking of 'peaks' in the sense of separate emission lines
# in the spectra that correspond to separate clouds


values, base = np.histogram(vel_peaks, bins=5)
cumulative = np.cumsum(values)
cumulative100 = cumulative/(cumulative[-1])

plt.plot(base[:-1], cumulative100)
plt.vlines(1, 0.7, 1, linestyle='--', color = 'red', label = '70% have 1')
plt.vlines(2, 0.7, 1, linestyle='--', color = 'orange', label = '92% have 2,1')
plt.vlines(3, 0.7, 1, linestyle='--', color = 'green', label = '97% have 3, 2, 1')
plt.xlabel('Number of velocity peaks')
plt.ylabel('Apertures')

plt.legend(loc='lower right')
plt.show()


peak_vels = []
for b in starting_bs:
    for l in starting_ls:
        peak_ind = velocity_peaks(b, l, 20)[1]
        peak_vels.append(peak_ind)
        
# print(peak_vels)
peak_vels = np.array(peak_vels)


count_1 = 0
for i in vel_peaks:
    if i == 1:
        count_1 +=1


count_2 = 0
for i in vel_peaks:
    if i ==2:
        count_2 +=1


count_3 = 0
for i in vel_peaks:
    if i == 3:
        count_3 +=1


count_4 = 0
for i in vel_peaks:
    if i == 4:
        count_4 +=1


count_5 = 0
for i in vel_peaks:
    if i == 5:
        count_5 +=1


(count_1/648)*100, (count_2/648)*100, (count_3/648)*100, (count_4/648)*100, (count_5/648)*100


num_peaks_2D = vel_peaks.reshape(18,36)

plt.figure(figsize=(15,3), dpi=200)

im = num_peaks_2D
imx = LAB_lons
imy = LAB_lats


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower', cmap='Paired')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical')

plt.title('Number of velocity peaks across the galaxy')
plt.show()


#for the velocity spectra that are wider than 25km/s -> 
#go back and plot them and see what they look like
#for all the numbers / assumptions i am making for now: keep in 
#mind that they are reasonable, but not the final value
#how can we make choices as agnostically as possible?


# peak_vels_flat = np.nanmean(peak_vels, axis = 1)
peak_vels_2D = peak_vels.reshape(18,36)

plt.figure(figsize=(15,3), dpi=200)

im = peak_vels_2D
imx = LAB_lons
imy = LAB_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin = vmin, vmax = vmax,\
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower')#, **kw_imshow)#, cmap='Paired')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.title('Location of velocity peaks across the galaxy')
plt.show()


#look into the specifics of how the data are taken to figure 
#out why there are really low negative velocities at higher latitudes
#could just be noise, but maybe they could be high-velocity-clouds (HVCs)


#repeat across-galaxy plots with better resolution (smaller apertures)
b0s_10 = np.arange(0,360,10)
l0s_10 = np.arange(0,720,10)


linewidths_10 = []
for b in b0s_10:
    for l in l0s_10:
        sigma = calculate_subcube_linewidth(b, l, 10)
        linewidths_10.append(sigma)
        
        
linewidths_10 = np.array(linewidths_10)


plt.hist(linewidths_10, 100)
plt.vlines(15, 0, 160, linestyle='--', color = 'red', label = '91% below')
plt.vlines(20, 0, 160, linestyle='--', color = 'orange', label = '96% are below')
plt.vlines(25, 0, 160, linestyle='--', color = 'green', label = '98% are below')
plt.xlabel('Linewidths [km/s]')
plt.ylabel('Apertures')

plt.legend()
plt.show()


len(linewidths_10)


count_sub25 = 0
for i in linewidths_10:
    if i < 25:
        count_sub25 +=1


(count_sub25/2592)*100 #, (count_sub20/648)*100, (count_sub25/648)*100


linewidths_10_2D = linewidths_10.reshape(36,72)

plt.figure(figsize=(15,3), dpi=200)

im = linewidths_10_2D
imx = LAB_lons
imy = LAB_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin = vmin, vmax = vmax,\
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'Linewidths [$km \ s^{-1}$]')

plt.title('Linewidths across the galaxy')
plt.show()


vel_peaks_10 = []
for b in b0s_10:
    for l in l0s_10:
        num_peaks = velocity_peaks(b, l, 10)[0]
        vel_peaks_10.append(num_peaks)
        
        
vel_peaks_10 = np.array(vel_peaks_10)


#----------


peak_vels_10 = []
for b in b0s_10:
    for l in l0s_10:
        peak_ind = velocity_peaks(b, l, 10)[1]
        peak_vels_10.append(peak_ind)
        
        
peak_vels_10 = np.array(peak_vels_10)


#make sure to go back and figure out a way to capture the smaller peaks (especially if there is multiple peaks)
#for find peaks, better to define what counts as a peak as height > noise, rather than height > 75%


values, base = np.histogram(vel_peaks_10, bins=5)
cumulative = np.cumsum(values)
cumulative100 = cumulative/(cumulative[-1])

plt.plot(base[:-1], cumulative100)
plt.vlines(1, 0.7, 1, linestyle='--', color = 'red', label = '83% have 1')
plt.vlines(2, 0.7, 1, linestyle='--', color = 'orange', label = '96% have 2,1')
plt.vlines(3, 0.7, 1, linestyle='--', color = 'green', label = '98% have 3, 2, 1')
plt.xlabel('Number of velocity peaks')
plt.ylabel('Apertures')

# plt.legend(loc='lower right') #change labels to reflect percentages first
plt.show()


count_3 = 0
for i in vel_peaks_10:
    if i == 3:
        count_3 +=1


(count_3/2592)*100 


num_peaks_10_2D = vel_peaks_10.reshape(36,72)

plt.figure(figsize=(15,3), dpi=200)

im = num_peaks_10_2D
imx = LAB_lons
imy = LAB_lats


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower',cmap='Paired')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical')

plt.title('Number of velocity peaks across the galaxy')
plt.show()


peak_vels_10_2D = peak_vels_10.reshape(36,72)

plt.figure(figsize=(15,3), dpi=200)

im = peak_vels_10_2D
imx = LAB_lons
imy = LAB_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin = vmin, vmax = vmax,\
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.title('Location of velocity peaks across the galaxy')
plt.show()


plt.figure(figsize=(15,3), dpi=200)

im = peak_vels_10_2D
imx = LAB_lons[l0s_10]
imy = LAB_lats[b0s_10]


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower', cmap = 'twilight')

plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')


# plt.contour(imx, imy, im, levels = np.arange(-45,45,5), cmap = 'twilight', linewidths = 0.7)
plt.contour(imx, imy, im, levels = np.arange(-45,45,5), cmap = 'twilight', linewidths = 1.5)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')

plt.title('Location of velocity peaks across the galaxy')
plt.show()


subcube = LAB_cube[:,130:230, 250:300]


from copy import deepcopy
comp = deepcopy(LAB)


#repeat velocity peak plots with even better resolution (smaller apertures)
#because getting the velocity peak accurate for each line of sight is
#important for selecting the right velocity channels to integrate -> rigt now
#using the average across the whole cube so whole chunks of the sky are getting 
#integrated where there isn't info

b0s_5 = np.arange(0,360,5)
l0s_5 = np.arange(0,720,5)


peak_vels_5 = []
for b in b0s_5:
    for l in l0s_5:
        peak_ind = velocity_peaks(b, l, 5)[1]
        peak_vels_5.append(peak_ind)
        
        
peak_vels_5 = np.array(peak_vels_5)


peak_vels_5_2D = peak_vels_5.reshape(72,144)

plt.figure(figsize=(15,3), dpi=200)

im = peak_vels_5_2D
imx = LAB_lons[l0s_5]
imy = LAB_lats[b0s_5]


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower', cmap = 'twilight')

plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.contour(imx, imy, im, levels = np.arange(-45,45,5), cmap = 'twilight', linewidths = 0.7)
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

# plt.contour(imx, imy, im, levels = np.arange(-45,45,20), cmap = 'twilight')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.title('Location of velocity peaks across the galaxy')
plt.show()


def get_bulk_surface_density(b0, l0):
    subcube = LAB_cube[:, b0:b0+5, l0:l0+5]
    
    comp = deepcopy(LAB)
    comp['cube'] = subcube
    comp['lats'] = LAB_lats[b0:b0+5]
    comp['lons'] = LAB_lons[l0:l0+5]

    sig_per_pix, a_per_pix = (cld.surface_density(comp, 1, '21cm', show_Sigma=False)).values()
    sig_per_pix_Msunpc2 = sig_per_pix.to(c.M_sun/u.pc**2)
    
    if np.all(np.isnan(sig_per_pix_Msunpc2)):
        print('found all nan slice at', b0, l0)
    bulk_Sigma = np.nanmean(sig_per_pix_Msunpc2)
    
    return bulk_Sigma.value
    

Sigmas_5 = []
for b in b0s_5:
    for l in l0s_5:
        Sig_at_px = get_bulk_surface_density(b, l)
        Sigmas_5.append(Sig_at_px)
        
        
Sigmas_5 = np.array(Sigmas_5)


Sigmas_5_2D = Sigmas_5.reshape(72,144)

mask = Sigmas_5_2D > 5

plt.figure(figsize=(15,3), dpi=200)

im = Sigmas_5_2D
imx = LAB_lons[l0s_5]
imy = LAB_lats[b0s_5]


vmin,vmax = zlims(im) 


plt.imshow(im*mask, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower', cmap = cmap0)

plt.colorbar(orientation='vertical', label = r'$\Sigma_{HI}$[$M_{\odot} \ pc^{-2}$]')

# plt.contour(imx, imy, im, levels = [5], linewidths = 0.7)
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

# plt.contour(imx, imy, im, levels = np.arange(-45,45,20), cmap = 'twilight')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.title(r'Surface Density Map ($2\degree.5$ $\times 2\degree.5$)')
plt.show()


def rms_noise(thing):
    mean = np.nanmean(thing)
    noise_square = (thing - mean)**2
    
    return np.sqrt(np.nanmean(noise_square)) 


CO_cubic_mask = np.zeros_like(CO_cube)

for c, ch in enumerate(CO_cube):
    noise_per_ch =  rms_noise(ch)
    ch_mask = ch > 3*noise_per_ch
    
    CO_cubic_mask[c] = ch_mask


CO_cubic_mask_inv = np.zeros_like(CO_cube)

for c, ch in enumerate(CO_cube):
    noise_per_ch =  rms_noise(ch)
    ch_mask = ch < 3*noise_per_ch
    
    CO_cubic_mask_inv[c] = ch_mask


CO_cubic_mask_norm_0th = (np.nansum(CO_cubic_mask, axis = 0))/(np.nansum(CO_cubic_mask, axis = 0))
CO_cubic_mask_norm_0th_inv = (np.nansum(CO_cubic_mask_inv, axis = 0))/(np.nansum(CO_cubic_mask_inv, axis = 0))


mask_from_CO_to_HI = transform_2D(CO_cubic_mask_norm_0th, CO_lats, CO_lons, LAB_lats[b0s_5], LAB_lons[l0s_5])
mask_from_CO_to_HI_inv = transform_2D(CO_cubic_mask_norm_0th_inv, CO_lats, CO_lons, LAB_lats[b0s_5], LAB_lons[l0s_5])


plt.figure(figsize=(15,3), dpi=200)

im = Sigmas_5_2D*mask*mask_from_CO_to_HI
# im2 = mask_from_CO_to_HI
imx = LAB_lons[l0s_5]
imy = LAB_lats[b0s_5]


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower', cmap = cmap0)

plt.colorbar(orientation='vertical', label = r'$\Sigma_{HI}$[$M_{\odot} \ pc^{-2}$]')

# plt.contour(imx, imy, im, levels = [5], linewidths = 0.7)
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.title(r'Surface Density Map ($2\degree.5$ $\times 2\degree.5$)')
plt.show()


np.all(mask_from_CO_to_HI==mask_from_CO_to_HI_inv)


plt.figure(figsize=(15,3), dpi=200)

im = Sigmas_5_2D*mask*mask_from_CO_to_HI_inv
# im2 = mask_from_CO_to_HI
imx = LAB_lons[l0s_5]
imy = LAB_lats[b0s_5]


vmin,vmax = zlims(im) 


plt.imshow(im, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]], origin='lower', cmap = cmap0)

plt.colorbar(orientation='vertical', label = r'$\Sigma_{HI}$[$M_{\odot} \ pc^{-2}$]')

# plt.contour(imx, imy, im, levels = [5], linewidths = 0.7)
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
# plt.colorbar(orientation='vertical', label = r'Velocity peak [$km \ s^{-1}$]')

plt.title(r'Surface Density Map ($2\degree.5$ $\times 2\degree.5$)')
plt.show()


def transform_2D(initial_data, initial_ax0_coord, initial_ax1_coord, target_ax0_coord, target_ax1_coord):
    #initial_data is a 2-d array, and its physical coordinates along axis = 0 & axis = 1 are given
    #as well as the physical coordinates of the target arlong axis = 0 and axis = 1


    #start out by saving important values
    init_ax0_len = len(initial_ax0_coord)
    init_ax1_len = len(initial_ax1_coord)
    tar_ax0_len = len(target_ax0_coord)
    tar_ax1_len = len(target_ax1_coord)
    

    #create an array of indices that point the coordinate arrays -> this is the data that's 
    #going to be transformed first, because then it will say which elements along each axis
    #of initial_data are going to be kept / moved over to initial_data_trans
    initial_ax0_coord_ins = np.arange(0, init_ax0_len)
    initial_ax1_coord_ins = np.arange(0, init_ax1_len)

    #then create the intermediate/useful arrays (the ones that contain the indices to be kept)
    init_ax0_ins_trans = transform_1D(initial_ax0_coord_ins, initial_ax0_coord, target_ax0_coord)
    init_ax1_ins_trans = transform_1D(initial_ax1_coord_ins, initial_ax1_coord, target_ax1_coord)


    #create empty array that the data will be transferred to
    initial_data_trans = np.zeros((tar_ax0_len, tar_ax1_len))

    #then go through and populate the new array with elements selected from initial_data 
    #the data are selected according to the transformed index arrays

    for i in range(tar_ax0_len-1):
        for j in range(tar_ax1_len-1):
            init_0th_coord = init_ax0_ins_trans[i]
            init_1th_coord = init_ax1_ins_trans[j]

            if np.isnan(init_0th_coord) or np.isnan(init_1th_coord): #ie padding had to be added to expand the initial to the target
                initial_data_trans[i][j] = False 
                #FALSE FOR NOW BECAUSE I WANT TO TRANSFORM A MASK, BUT WILL BE BETTER TO HAVE THIS (THE
                #VALUE THAT GETS SET WHERE THERE'S NO DATA) TO BE AN INPUT

            else: #real coordinate values
                init_0th_coord = int(init_0th_coord) #should already be an integer, but fixing the type
                init_1th_coord = int(init_1th_coord) #should already be an integer, but fixing the type

                initial_data_trans[i][j] = initial_data[init_0th_coord][init_1th_coord]

    return initial_data_trans


#repeat velocity peak plots with even better resolution (smaller apertures)
#because getting the velocity peak accurate for each line of sight is
#important for selecting the right velocity channels to integrate -> rigt now
#using the average across the whole cube so whole chunks of the sky are getting 
#integrated where there isn't info

b0s_2 = np.arange(0,360,2)
l0s_2 = np.arange(0,720,2)


Sigmas_2 = []
for b in b0s_2:
    for l in l0s_2:
        Sig_at_px = get_bulk_surface_density(b, l)
        Sigmas_2.append(Sig_at_px)
        
        
Sigmas_2 = np.array(Sigmas_2)


Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)


orion = cld.extract_cloud(17, Zucker_df, DHT_CO, p_padding=60)


#surface density maps for solar neighborhood clouds + their HI envelopes FROM LAB
for i in range(len(Zucker_df)):
    cloud, parent = cld.get_parent_HI(i, Zucker_df, DHT_CO, LAB, \
                      CO_padding=20, HI_padding=5,  v_padding=10, \
                      show_0th = False, show_spec = False, sigma = 3).values()
    
    cld.surface_density(cloud['comp'],1.5, 'CO', show_Sigma=True)
    cld.surface_density(parent['comp'],1.5, '21cm', show_Sigma=True)


#0th moment maps & velcoity spectra for solar neighborhood clouds + their HI envelopes FROM LAB

for i in range(len(Zucker_df)):
    cld.get_parent_HI(i, Zucker_df, DHT_CO, LAB, \
                      CO_padding=20, HI_padding=5,  v_padding=10, \
                      show_0th = True, show_spec = True, sigma = 3)


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
