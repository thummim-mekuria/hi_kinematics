#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


import astropy.units as u
from astropy.utils import data
from spectral_cube import SpectralCube

import matplotlib.pyplot as plt

import datacube as dc
import clouds as cld
import numpy as np


HI4PI_dc['header']


HI4PI_dc = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')


HI4PI_dv = (HI4PI_dc['vels'][1] - HI4PI_dc['vels'][0])*HI4PI_dc['units']['vels']


HI4PI_cube = SpectralCube.read('HI4PI_CAR.fits')


HI4PI_subcube = HI4PI_cube[:, 1000:1050, 1000:1050]


HI4PI_subcube


sample_spec = HI4PI_cube[:, 1000:1050, 1000:1050].mean(axis=(1,2))


plt.plot(sample_spec)


dc.dispersion_1sigma(e.value) * HI4PI_dv #examle velocity dispersion / line width along a line of sight


#okay so one thing i can do is: make an evenly spaced sample 
#of apertures across the entire galaxy, then take spectra along each of them, 
#then take the average?median? width of the spectra (if there is one) as a characteristic
#width to use for the velocity channels (or make a histogram of the spectra fwhms ... hmmm)
#what would be a good aperture size though?


#this is the 0th moment map of the subcube that i just took the spectrum of ... 
#the aperture feels too small ... 
plt.imshow(HI4PI_subcube.moment(order=0).value, origin = 'lower')
plt.colorbar()

plt.show()


# went back to look at the HI clouds associated with the GMCs, and it seems like their (=the GMCs) sizes are: 
# 
# (==> will have to do some sort of math equivalent to the acumulation radius to get from here to HI sizes, and also need to know something about distance in order to turn aperture sizes into physical extents in pc)
# 
# (these are by eye, form looking at the 0th moment maps plotted by get_parent_HI)
# 
# aquila:     b = 10 x l = 7
# cma:        b = 7  x l = 8
# california: b = 7  x l = 10
# cam:        b = 7  x l = 8
# cepheus:    b = 14 x l = 21
# chamaeleon: b = 10 x l = 14
# corona aus: b = 10 x l = 9
# crossbones: b = 7  x l = 7
# gemini:     b = 7  x l = 10
# hercules:   b = 7  x l = 9
# lacerta:    b = 8  x l = 7
# lupus:      b = 18 x l = 15
# maddalena:  b = 8  x l = 7
# mon_ob1:    b = 7  x l = 9
# mon_r2:     b = 7  x l = 8
# ophiuchus:  b = 9  x l = 9
# orion A:    b = 7  x l = 12
# orion B:    b = 8  x l = 10
# orion lam:  b = 13 x l = 12
# perseus:    b = 10 x l = 9
# polaris:    b = 12 x l = 12
# rosette:    b = 6  x l = 8
# serpens:    b = 8  x l = 9
# taurus:     b = 9  x l = 14
# 


bs = np.array((10,7,7,7,14,10,10,7,7,7,8,18,8,7,7,9,7,8,13,10,12,6,8,9))


ls = np.array((7,8,10,8,21,14,9,7,10,9,7,15,7,9,8,9,12,10,12,9,12,8,9,14))


np.median(bs), np.mean(bs)


np.median(ls), np.mean(ls)


# the median (& mean) latitude extents are 8 (& 9)
# the median (& mean) longitude extents are 9 (& 10)
# 
# db = 0.083 degrees & dl = -0.0833 degrees (!!! 6x better resolution than LAB)
# so the median aperture would be: 96 x 108 indices
# & the mean aperture would be: 108 x 120 indices
# 
# 
# lets say we go with apertures for sampling that are ~100 x 100 indices (which again, technically is the size of the GMCS not the HI, but it is a starting point for the sample velocity spectra)
# 
# 
# the HI4PI cube is 2163 (b) x 4323 (l), so if we were to sample every single 100x100 aperture across it, it would be ~21 x ~43 = ~900 spectra .... which feels doable? 
# 
# the thing that i then want to do is ... (wow this is going to be a big for loop probably)
# 1) go through and create an array of b starting positions, and l starting positions,
# 2) slice out subcubes at those coordinates
# 3) take velocity spectra of each subcube 
# 4) crunch the spectra through dispersion1sigma to calculate the linewidths
# 5) add those to an array
# 6) plot the histogram of the array to see if there is a 'characteristic' linewidth


HI4PI_db = (HI4PI_dc['lats'][1] - HI4PI_dc['lats'][0])*HI4PI_dc['units']['lats']
HI4PI_dl = (HI4PI_dc['lons'][1] - HI4PI_dc['lons'][0])*HI4PI_dc['units']['lons']


starting_bs = np.arange(0,2100,100)
starting_ls = np.arange(0,4300,100)


def calculate_HI4PIsubcube_linewidth (b0, l0):
    subcube = HI4PI_cube[:, b0:b0+100, l0:l0+100]
    
    spectrum = subcube.mean(axis=(1,2))
    
    linewidth = dc.dispersion_1sigma(spectrum.value) * HI4PI_dv
    
    return linewidth


linewidths = []


for b in starting_bs:
    for l in starting_ls:
        sigma = calculate_HI4PIsubcube_linewidth(b, l)
        linewidths.append(sigma)


len(linewidths) #noice that this is 903 long as expected 


linewidths = np.array((linewidths))


linewidths


count_sub10 = 0
for i in linewidths:
    if i < 10000:
        count_sub10 +=1


count_sub10


753/903


plt.hist(linewidths, 100)
plt.vlines(10000, 0, 65, linestyle='--', color = 'red')
plt.show()


count_sub15 = 0
for i in linewidths:
    if i < 15000:
        count_sub15 +=1


count_sub15


847/903


# almost 94 % of the spectra have linewidths of 15km/s or less, so that seems like an okay (for now) velocity channel width to use (hmmmmm in nia's m33 paper, she used 50km/s -- the wider velocity dispersion maybe has to do with larger angular extent?))
# 
# 
# also just realized that this is for HI4PI ... but not really (right now) able to plot a surface density map for that, so re-do for LAB
# 
# for LAB, because resolution is 0.5 degrees, going to use an aperture of 20 x 20 indices (10 degrees square)


LAB = dc.unpack_fits_cube('labh.fits', 'LAB')


LAB_cube = LAB['cube']
LAB_dv = (LAB['vels'][1] - LAB['vels'][0]) * LAB['units']['vels']
LAB_db = (LAB['lats'][1] - LAB['lats'][0]) * LAB['units']['lats']
LAB_dl = (LAB['lons'][1] - LAB['lons'][0]) * LAB['units']['lons']


def calculate_LABsubcube_linewidth (b0, l0):
    subcube = LAB_cube[:, b0:b0+20, l0:l0+20]
    
    spectrum = subcube.mean(axis=(1,2))
    
#     plt.plot(spectrum)
#     plt.show()
    
    linewidth = dc.dispersion_1sigma(spectrum) * LAB_dv
    
    return linewidth


LAB_cube.shape


starting_bs_LAB = np.arange(0,360,20)
starting_ls_LAB = np.arange(0,720,20)


linewidths_LAB = []


for b in starting_bs_LAB:
    for l in starting_ls_LAB:
        sigma = calculate_LABsubcube_linewidth(b, l)
        linewidths_LAB.append(sigma)


linewidths_LAB = np.array(linewidths_LAB)


len(linewidths_LAB)


plt.hist(linewidths_LAB, 100)
plt.vlines(10000, 0, 40, linestyle='--', color = 'red')
plt.show()


count_sub15 = 0
for i in linewidths_LAB:
    if i < 15000:
        count_sub15 +=1


count_sub20 = 0
for i in linewidths_LAB:
    if i < 20000:
        count_sub20 +=1


count_sub25 = 0
for i in linewidths_LAB:
    if i < 25000:
        count_sub25 +=1


count_sub25


# about 98% of the spectra in LAB have a velocity dispersion <= 25km/s, so sounds like a pretty okay velocity channel width to use for making the surface density maps
# 
# will need to edit the function so it finds the velocity peak, then only integrates 25 velocity channels centered on the peak (the index with the peak, + 12 to the -v and 12 to the +v 


(25/LAB_dv)*1000


subcube = LAB_cube[:, 200:220, 500:520]


spec = dc.spectrum1D(subcube)


np.where(spec == np.max(spec))


plt.plot(spec)
plt.vlines(439, 0, 25, color = 'red', linestyle = '--')


from scipy.signal import peak_widths, find_peaks

maxim = np.nanmax(spec)
    
p = find_peaks(spec, height = 0.75 * maxim)
#first unpack what find_peaks returns
peak = p[0]
peak_info = p[1] 


# so the thing i want do is find the peak of the velocity spectrum, then find the location (in velocity indices) of where that peak is, then define the velocity window as +/- 12 channels from the location of the channel where the max value is, then do the integration to calculate surface density just across those 25 channels


peak


peak_info


peak_widths(spec, peak) 
#returns 0 = fwhm, 1 = 0.5 * peak height (where fwhm was evaluated), 2, 3 = left, right inds where fwhm was evaluated


if len(peak) > 1:#ie if there are multiple peaks found that are above the threshold
    #look at the peak heights & find the location of the tallest one
    i = np.argmax(peak_info['peak_heights']) 
    #then use that tallest peak for fwhm purposes
    peak = [peak[i]] #use a list (even tho 1 element) because peak_widths wants an iterable

fwhm = peak_widths(spectrum, peak)[0] #the first from peak_width output is fwhm
sigma = fwhm/sqrt2ln2x2 


#sum across narrow velocity channel
cld.surface_density(LAB, 1.5,'21cm' )


#sum across narrow velocity channel
cld.surface_density(LAB, 1.5,'21cm' )


#sum across entire cube
cld.surface_density(LAB, 1.5,'21cm' )


638/648


#also this will be a side project -- figure out whether it is worth the effort to edit the output from datacube
#to be an object, similar to spectral cube ... or to just straight-up use spectral cube


0.5/0.083


HI4PI_0th = HI4PI_cube.moment(order=0)


plt.imshow(HI4PI_0th.value, origin = 'lower')
plt.colorbar()

plt.show()


from astropy.convolution import Gaussian1DKernel


kernel = Gaussian1DKernel(2.5)


new_cube = HI4PI_cube.spectral_smooth(kernel)


HI4PI_sample_spec = HI4PI_cube[]


HI4PI_cube


HI4PI_wheel = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')


HI4PI_wheel['vels'][0], HI4PI_wheel['vels'][-1]


cube = SpectralCube.read(fn)
print(cube)
SpectralCube with shape=(7, 4, 3) and unit=Jy / beam:
n_x:      3  type_x: RA---ARC  unit_x: deg    range:    52.231466 deg:   52.231544 deg
n_y:      4  type_y: DEC--ARC  unit_y: deg    range:    31.243639 deg:   31.243739 deg
n_s:      7  type_s: VRAD      unit_s: m / s  range:    14322.821 m / s:   14944.909 m / s

# extract the subcube between 98 and 100 GHz
slab = cube.spectral_slab(98 * u.GHz, 100 * u.GHz)  

# Ignore elements fainter than 1 Jy/beam
masked_slab = slab.with_mask(slab > 1 Jy/beam)  

# Compute the first moment and write to file
m1 = masked_slab.moment(order=1)  
m1.write('moment_1.fits') 


def transform_1D(initial_data, initial, target):

    #initial & target here contain coordinate information -> initial_data is the actual information
    #that will be transformed to fit / be in the new coordinates defined by target

    #save important values
    init_first = initial[0]
    tar_first = target[0]
    init_last = initial[-1]
    tar_last = target[-1]

    #all of this is also assuming that d1 is the spacing  everywhere in initial, 
    #and d2 is the spacing everywhere in target
    #for now, not worrying about the case in which this isn't true --
    # but might need to revisit this
    d1 = initial[1] - initial[0]
    d2 = target[1] - target[0]
    d2d1 = d2/d1 #turns out to be a number that comes up a lot


    ###################################################

    #declare all the important booleans as false
    #they turn true when conditions are fulfilled

    #flags for check 1
    target_starts_before_initial = False 
    target_ends_after_initial = False
    they_overlap = False

    # flags for check 2
    nice_spacing = False
    need_interp = False

    #flags for check 3
    init_reversed = False

    #flags for check 5
    in_phase = False

    ###################################################


     #------- step 1: check whether the arrays ascending in the same direction ------- 
            
    if d2d1 < 0: #they are going in opposite directions
        initial_orig = initial 
        initial = np.flip(initial_orig) #reverse initial array to go in the same direction as target

        d2d1_orig_2 = d2d1 #this is a different instance of saving a value for d2d1
        d2d1 = -d2d1_orig_2
        init_reversed = True
                    
    #----------------------------------- step 1: -------------------------------------- 


    #-------------- step 2: the most important check: do the arrays overlap & to what extent ---------------------      
    
    #because of the previous check, we know that the two arrays are now increasing in the same direction
    #as each other -- but we still want to know if both are in ascending or descending order
    if d2 > 0: #whatever the sign is of d2, that's the direction in which both arrays are now going
       
        # case 1: initial starts first or they start at the same number
        if (init_first <= tar_first): 
             if (tar_first < init_last):
                they_overlap = True

                if (tar_last > init_last):
                    target_ends_after_initial = True
                    
            #else they don't overlap


        #case 2: target starts first (tar_first < init_first)
        else: 
            target_starts_before_initial = True 
    
            if (init_first < tar_last):
                they_overlap = True
                  
                if (tar_last > init_last):  
                    target_ends_after_initial = True
                    #might need to modify this a little based on Vv (case of the lons), but works for now
                    #this checks that the end of initial (or at least some portion of last_subset_init) are inside the end of target
                    #the way this is defined is kind of circular, but the main idea is that you want to find the closest overlap and 
                   #that might be a few in-to initial, and the last element might still overshoot (look at the case of CO_lons & HI_lons)
                #else initial ends last (encompasses target) 

            #else they don't overlap         


    else: #d2 < 0:

        # case 1: initial starts first or they start at the same number
        if (init_first >= tar_first):
            if (tar_first > init_last):
                they_overlap = True

                if (tar_last < init_last):
                    target_ends_after_initial = True

                #else initial ends last (encompasses target)

            #else they don't overlap

       
        #case 2: target starts first (tar_first > init_first)
        else:
            target_starts_before_initial = True

            if (init_first > tar_last):
                they_overlap = True

                if (tar_last < init_last):
                    target_ends_after_initial = True

                #else initial ends after target (initial is not encompassed)

            #else they don't overlap

    #------------------------------------------------- step 2 ------------------------------------------------------

    
    #----------------------------------- step 3: second most important check: how do the spacings relate to each other --------------
    if (d2%d1 == 0): 
        nice_spacing = True
        d2d1 = int(d2d1) 
        #if the above two conditions are met, it should already be an integer, this is just for type purposes

        #the spacing of the initial is denser, so no need to calculate interpolated values
        #this is the easiest scenario
    
    elif (d1%d2 == 0): 
        d2d1_orig = d2d1
        d2d1 = int(1/d2d1_orig) #use the ratio that's of spacings that's a whole number 

        nice_spacing = True
        need_interp = True # ie the spacing of the initial is sparser, will need to interpolate 

    #else their spacing doesn't fit neatly 
        #------------------------------------------------------------ step 3 -----------------------------------------------------------


    if they_overlap: #the transformation will have physical meaning only if this is true
    #now begins the work of setting up the transformation
        #---------------------  step 4: set transformation starting and end points --------------------------------------- 
        if target_starts_before_initial:
            start = 0 #where to start itierating inside initial
        else: 
        #initial starts before target, so will need to pick a subset from the middle of initial
            start = np.argmin(np.abs(initial - target[0]))


        if target_ends_after_initial:
            end = len(initial)
        else:
        #initial's end is outside of target so we need to find the last element within intial to iterate to
            end = 1 + np.argmin(np.abs(initial - target[-1]))
        #--------------------------------------- step 4: ----------------------------------------------------------------- 


        #------------------------------------- step 5: find point(s) of closest approach/intersection ---------------------------------------- 
        if nice_spacing: 
            #if this is true, there are only two options: 1) they're in phase, ie they have elements that are equal to each 
            #other periodically by index, or 2) they are out of phase ie none of their elements are equal to one another

            #to test in/out-of-phase-ness, take a subset of the denser array & compare to the other 
            #if they're in phase (and we already know they overlap), there will be at least one exact element from the subset of the 
            #denser array that is equal to an element form the other, and if they are out of phase, none of the elements from the
            #subset of the denser array will match an element from the sparser

            if need_interp: #target is the denser array & interpolation needed (d1 > d2)
                denser = target #array from which to take subset
                sparser = initial
                subset_size = d2d1

                #start and end points hadn't been set for target like they were for initial (in step 4, so doing that now)
                #for the subset searching for the last matching pair
                if target_starts_before_initial:
                    first_subset_start = np.argmin(np.abs(initial[0] - target))
                else:
                    first_subset_start = 0 

                #for the subset searching for the last matching pair   
                if target_ends_after_initial:
                    last_subset_end = 1 + np.argmin(np.abs(initial[-1] - target))
                else:
                    last_subset_end = len(target)

                
            elif not need_interp: #initial is the denser array, interpolation not needed (d2 > d1)
                denser = initial
                sparser = target
                
                #for the subset searching for the last matching pair
                first_subset_start = start 
                #for the subset searching for the last matching pair   
                last_subset_end = end

                if d2d1 == 1: #special case, in which we need to select d2 number of elements in subset
                    subset_size =  int(d2)+1 #in case d2 is a decimal, +1 accounts for: int rounds down but need to be rounding up 
                else:
                    subset_size = d2d1


            first_subset = denser[first_subset_start : first_subset_start + subset_size]


            for i, ith in enumerate(first_subset):
                compare = np.abs(ith - sparser)
                j = np.argmin(compare)

                if ith == sparser[j]: #the subset element is exactly equal to an element in the sparser array
                    i += first_subset_start #update i since it was only indexing within subset
                    in_phase = True
                    break


            if in_phase: #at some point in the for loop, we found an exact match, which is just the the first matching pair 
                #of elements between the two arrays, & since we know now that they're in phase, 
                #we also want to pin down the last matching pair

                last_subset = denser [-subset_size + last_subset_end: last_subset_end]
                for k, kth in enumerate(last_subset):
                    compare = np.abs(kth - sparser)
                    m = np.argmin(compare)

                    if (kth == sparser[m]):
                        k += last_subset_end - subset_size #update k, since before it was only indexing within last_subset
                        break


                #finally we can construct the new array
                if need_interp: 
                    #so now we have the following information: i & k (= positions of first/last matching pair in denser),
                    #j & m (= positions of first/last matching pair in sparser), d2d1 -> factor by which desnser has more 
                    #elements (ie how many spaces in denser between the elements that match every single one in sparser)
                    
                    #create an array thats the same length as how many indices are between the first & 
                    #last matching pairs, and that is as wide as how many interpolations will need to happen
                    #SOMETHING IS WRONG WITH HOW THIS LONG IS GETTING INITIALIZED TO BE / INDEXING INTO DENSER & SPARER
                    
                    init_data_trans = np.zeros(((m-j +1),d2d1)) 
                    for s in range(j, m+1):
                        interp_padding = np.full((d2d1-1), np.nan)
                        
                        #this is to fill up the space bc initial is being expanded to be the same size as target
                        #might be a good idea to look into actual interpolation math between the data points instead of just nans

                        interp_unit = np.append(initial_data[s], interp_padding)
                        
                        init_data_trans[s] = interp_unit
                    
                    init_data_trans = init_data_trans.flatten()
    
                    if target_starts_before_initial:
                        padding_start = np.full(len(target[:i]), np.nan) #create an array of nans the same length as how many items are in padding before initial 
                        init_data_trans = np.append(padding_start, init_data_trans) 


                    if target_ends_after_initial:
                        padding_end = np.full(len(target[k+d2d1-1:]), np.nan) #don't want to include kth through kth+d2d1-1 elements cos that was accounted for by interp_padding
                        init_data_trans = np.append(init_data_trans, padding_end) #add the stuff thats at the end of target past initial's end

        
                else: #interpolation isn't needed -> ie initial = denser & target = sparser
  
                    init_data_trans = initial_data[i: k+1 :d2d1] #subset of initial that matches with target, selecting only the nth elements

                    if target_starts_before_initial:
                        padding_start = np.full(len(target[:j]), np.nan) #create an array of nans the same length as how many items are in padding before initial 
                        init_data_trans = np.append(padding_start, init_data_trans) 


                    if target_ends_after_initial:
                        padding_end = np.full(len(target[m+1:]), np.nan) #don't want to include mth element cos that still exists in initial
                        init_data_trans = np.append(init_data_trans, padding_end) #add the stuff thats at the end of target past initial's end

                        
                return init_data_trans
            
            else: #for loop completed without exact match => they're completely out of phase 
                print("they're out of phase, physically not sure what this would entail, but that's a problem for later")


        else: #not nice spacing
            #DECIDE IF THIS NEEDS TO ALSO BE BROKEN DOWN INTO INTERP NEEDED / NOT NEEDED, OR IF ONE FORMULA WILL HANDLE BOTH CASES
            print('not nice spacing')
            print('d1 is',d1,'d2 is',d2)
            print('d2d1 is', d2d1)
            print('start is',start)
            print('end is', end)
            
            print('in target, start corresponds to',np.argmin(np.abs(initial[start] - target)))
            print('in target, end corresponds to',np.argmin(np.abs(initial[end-1] - target)))
            for i, init in enumerate(initial):
                j = np.argmin(np.abs(init - target))
        

        #--------------------------------------------------------------- step 5 --------------------------------------------------------- 

    
    else:
        print('the two arrays do not coincide physically, the transformation has no meaning')
 

HI_vels_km = HI_vels*HI_units['vels'].to(CO_units['vels'])


from scipy.ndimage import binary_dilation


mask_CO = CO_0th > 3*rms(CO_0th)


plt.imshow(mask_CO)


mask_CO_dil = binary_dilation(mask_CO, iterations=10)


im = np.array(((2,0,0),(4,0,0)))


imdil = binary_dilation(im)


plt.imshow(im)


imdil


plt.imshow(imdil)


# question: how much resolution can i afford to lose when transforming from CO_vels to HI_vels
# ^^ the main thing that will determine this is the width of the smallest velocity profile that i want to resolve
# because if the smallest velocity peak is ~10km.s => decide what percentage of that (10%, 5%) is okay to lose, ie a loss of data at the level of 0.5-1km.s would be okay
# 
# so maybe the thing that i need to do first is go through and look at the velocity profiles along different lines of sight to see if there is a 'characteristic' velocity channel width (<= this will also help with figuring out surface density mapping)


plt.imshow(mask_CO_dil)


HI_vels_km[755]


CO_vels[492]


HI_vels_km[134]


CO_vels[0]


493-0


755-134


493/621


transform_1D(CO_vels, CO_vels, HI_vels_km)


transform_1D(HI_lons, HI_lons, CO_lons)


r = np.full((3, 5), 8)


np.append(np.array((np.nan, np.nan)), r)


r


CO_lons[-4:]


f.shape


CO_lons.shape


m = [None]


m.append([0])


m


def transform_2D(initial_data, initial_ax0_coord, initial_ax1_coord, target_ax0_coord, target_ax1_coord):
    #initial_data is a 2-d array, and its physical coordinates along axis = 0 & axis = 1 are given
    #as well as the physical coordinates of the target arlong axis = 0 and axis = 1


    #start out by saving important values
    init_ax0_len = len(initial_ax0_coord)
    init_ax1_len = len(initial_ax1_coord)

    tar_ax0_len = len(target_ax0_coord)
    tar_ax1_len = len(target_ax1_coord)


    initial_ax0_coord_ins = np.arange(0, init_ax0_len)
    initial_ax1_coord_ins = np.arange(0, init_ax1_len)

    #then create the intermediate/useful arrays
    init_ax0_ins_trans = transform_1D(initial_ax0_coord_ins, initial_ax0_coord, target_ax0_coord)
    init_ax1_ins_trans = transform_1D(initial_ax1_coord_ins, initial_ax1_coord, target_ax1_coord)


    #create empty array that the data will be transferred to
    initial_data_trans = np.zeros((tar_ax0_len, tar_ax1_len))

    #then go through and populate the new array with elements selected from initial_data 
    #the data are selected according to the transformed coordinate arrays

    #AS IT IS, THIS ONLY WORKS FOR DATA THAT DOESN'T NEED INTERPOLATION (INITIAL IS DENSER THAN TARGET)
    #BECAUSE IF INTERPOLATION IS NEEDED, WILL HAVE TO DO MORE MATH -- NOT TACKLING THAT FOR NOW THO
    print(init_ax0_ins_trans)#[0], init_ax0_ins_trans[-1])

    for i in range(tar_ax0_len):
        for j in range(tar_ax1_len):
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


CO_l_ins = np.arange(0, len(CO_lons))
CO_b_ins = np.arange(0, len(CO_lats))


def transform_jptr(mask_CO, CO_lats, CO_lons, HI_lats, HI_lons):

    CO_l_ins_trans = transform_1D(CO_l_ins, CO_lons, HI_lons)
    CO_b_ins_trans = transform_1D(CO_b_ins, CO_lats, HI_lats)

    new_mask = np.zeros_like(HI_0th)


    for b in range(len(new_mask)):
        l_row = new_mask[b]

        for l in range(len(l_row)):
            new_lat = CO_b_ins_trans[b]
            new_lon = CO_l_ins_trans[l]

            if np.isnan(new_lat) or np.isnan(new_lon):
                new_mask[b][l] = False #np.nan
            else:
                new_lat = int(new_lat)
                new_lon = int(new_lon)
                new_mask[b][l] = mask_CO[new_lat][new_lon]
            
            
    return new_mask


mask_CO_inv = ~mask_CO


mask_CO = CO_0th > 0.5*rms(CO_0th)


plt.imshow(new_mask_f_inv)


new_mask_f_inv = transform_2D(mask_CO_inv, CO_lats, CO_lons, HI_lats, HI_lons)


new_mask_f = transform_2D(mask_CO, CO_lats, CO_lons, HI_lats, HI_lons)


mask_CO.shape


CO_0th_masked = mask_CO * CO_0th


new_mask_inv = np.zeros_like(new_mask)


for b in range(len(new_mask_inv)):
    l_row = new_mask_inv[b]
    
    for l in range(len(l_row)):
        new_lat = CO_b_ins_trans[b]
        new_lon = CO_l_ins_trans[l]
        
        if np.isnan(new_lat) or np.isnan(new_lon):
            new_mask_inv[b][l] = True #np.nan
        else:
            new_lat = int(new_lat)
            new_lon = int(new_lon)
            new_mask_inv[b][l] = ~mask_CO[new_lat][new_lon]


for b in range(len(new_mask)):
    l_row = new_mask[b]
    
    for l in range(len(l_row)):
        new_lat = CO_b_ins_trans[b]
        new_lon = CO_l_ins_trans[l]
        
        if np.isnan(new_lat) or np.isnan(new_lon):
            new_mask[b][l] = False #np.nan
        else:
            new_lat = int(new_lat)
            new_lon = int(new_lon)
            new_mask[b][l] = mask_CO[new_lat][new_lon]


np.zeros((2,3))


np.max(new_mask_inv)


np.all(new_mask == new_mask_f)


HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = new_mask_f_inv*HI_0th
imx = HI_lons
imy = HI_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')

# im2 = CO_0th_masked
# im2x = CO_lons
# im2y = CO_lats 

# plt.contour(im2x,im2y,im2,levels=[0.1], linewidths = 0.7, colors='xkcd:pale orange')


plt.title('Galaxy in HI: 0th moment')
plt.show()


HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = new_mask
imx = HI_lons
imy = HI_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')

im2 = CO_0th_masked
im2x = CO_lons
im2y = CO_lats 

plt.contour(im2x,im2y,im2,levels=[0.1], linewidths = 0.7, colors='xkcd:pale orange')


plt.title('Galaxy in HI: 0th moment')
plt.show()


plt.figure(figsize=(15,3), dpi=200)

im = new_mask
imx = HI_lons
imy = HI_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin = vmin, vmax =vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.colorbar(orientation='vertical', label = r'$W_{CO}$ $(K\times km/s)$')

# plt.contour(imx,imy,im,levels=[0.1], linewidths = 0.7, colors='xkcd:pale orange')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')


plt.title('Galaxy in CO: 0th moment')
plt.show()


### GALAXY 0TH MOMENT (in CO) ####
plt.figure(figsize=(15,3), dpi=200)

im = CO_0th_masked
imx = CO_lons
imy = CO_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin = vmin, vmax =vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)

# plt.contour(imx, imy, im, levels = [0.1])

plt.colorbar(orientation='vertical', label = r'$W_{CO}$ $(K\times km/s)$')

plt.contour(imx,imy,im,levels=[0.1], linewidths = 0.7, colors='xkcd:pale orange')


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')


plt.title('Galaxy in CO: 0th moment')
plt.show()


CO_0th = np.nansum(CO_cube*CO_dv.value, axis=0) 


### GALAXY 0TH MOMENT (in CO) ####
plt.figure(figsize=(15,3), dpi=200)

im = CO_0th
imx = CO_lons
imy = CO_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{CO}$ $(K\times km/s)$')


plt.title('Galaxy in CO: 0th moment')
plt.show()


HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = HI_0th
imx = HI_lons
imy = HI_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')

im2 = CO_0th_masked
im2x = CO_lons
im2y = CO_lats 

plt.contour(im2x,im2y,im2,levels=[0.1], linewidths = 0.7, colors='xkcd:pale orange')


plt.title('Galaxy in HI: 0th moment')
plt.show()


from astropy import units as u
import astropy.constants as c


def transform_1D(initial_data, initial_coord, target):

    #naming this shorter for brevity's sake 
    # will need to refer to it lots 
    #target is only the coordinate information, the actual data isn't needed
    initial = initial_coord


    #save important values
    init_first = initial[0]
    tar_first = target[0]
    init_last = initial[-1]
    tar_last = target[-1]

    #all of this is also assuming that d1 is the spacing  everywhere in initial, 
    #and d2 is the spacing everywhere in target
    #for now, not worrying about the case in which this isn't true --
    # but might need to revisit this
    d1 = initial[1] - initial[0]
    d2 = target[1] - target[0]
    d2d1 = d2/d1 #turns out to be a number that comes up a lot


    ###################################################

    #declare all the important booleans as false
    #they turn true when conditions are fulfilled

    #flags for check 1
    target_starts_before_initial = False 
    target_ends_after_initial = False
    they_overlap = False

    # flags for check 2
    nice_spacing = False
    need_interp = False

    #flags for check 3
    init_reversed = False

    #flags for check 5
    in_phase = False

    ###################################################


     #------- step 1: check whether the arrays ascending in the same direction ------- 
            
    if d2d1 < 0: #they are going in opposite directions
        initial_orig = initial 
        initial = np.flip(initial_orig) #reverse initial array to go in the same direction as target

        d2d1_orig_2 = d2d1 #this is a different instance of saving a value for d2d1
        d2d1 = -d2d1_orig_2
        init_reversed = True
                    
    #----------------------------------- step 1: -------------------------------------- 


    #-------------- step 2: the most important check: do the arrays overlap & to what extent ---------------------      
    
    #because of the previous check, we know that the two arrays are now increasing in the same direction
    #as each other -- but we still want to know if both are in ascending or descending order
    if d2 > 0: #whatever the sign is of d2, that's the direction in which both arrays are now going
       
        # case 1: initial starts first or they start at the same number
        if (init_first <= tar_first): 
             if (tar_first < init_last):
                they_overlap = True

                if (tar_last > init_last):
                    target_ends_after_initial = True
                    
            #else they don't overlap


        #case 2: target starts first (tar_first < init_first)
        else: 
            target_starts_before_initial = True 
    
            if (init_first < tar_last):
                they_overlap = True
                  
                if (tar_last > init_last):  
                    target_ends_after_initial = True
                    #might need to modify this a little based on Vv (case of the lons), but works for now
                    #this checks that the end of initial (or at least some portion of last_subset_init) are inside the end of target
                    #the way this is defined is kind of circular, but the main idea is that you want to find the closest overlap and 
                   #that might be a few in-to initial, and the last element might still overshoot (look at the case of CO_lons & HI_lons)
                #else initial ends last (encompasses target) 

            #else they don't overlap         


    else: #d2 < 0:

        # case 1: initial starts first or they start at the same number
        if (init_first >= tar_first):
            if (tar_first > init_last):
                they_overlap = True

                if (tar_last < init_last):
                    target_ends_after_initial = True

                #else initial ends last (encompasses target)

            #else they don't overlap

       
        #case 2: target starts first (tar_first > init_first)
        else:
            target_starts_before_initial = True

            if (init_first > tar_last):
                they_overlap = True

                if (tar_last < init_last):
                    target_ends_after_initial = True

                #else initial ends after target (initial is not encompassed)

            #else they don't overlap

    #------------------------------------------------- step 2 ------------------------------------------------------

    
    if they_overlap: #the transformation will have physical meaning only if this is true
        #----------------------------------- step 3: second most important check: how do the spacings relate to each other --------------
        if (d2%d1 == 0): 
            nice_spacing = True
            d2d1 = int(d2d1) 
            #if the above two conditions are met, it should already be an integer, this is just for type purposes

            #the spacing of the initial is denser, so no need to calculate interpolated values
            #this is the easiest scenario
        
        elif (d1%d2 == 0): 
            d2d1_orig = d2d1
            d2d1 = int(1/d2d1_orig) #use the ratio that's of spacings that's a whole number 

            nice_spacing = True
            need_interp = True # ie the spacing of the initial is sparser, will need to interpolate 

        #else their spacing doesn't fit neatly 
        #------------------------------------------------------------ step 3 -------------------------------------------------------------


        #---------------------  step 4: set transformation starting and end points --------------------------------------- 
        if target_starts_before_initial:
            start = 0 #where to start itierating inside initial
        else: 
        #initial starts before target, so will need to pick a subset from the middle of initial
            start = np.argmin(np.abs(initial - target[0]))


        if target_ends_after_initial:
            end = len(initial)
        else:
        #initial's end is outside of target so we need to find the last element within intial to iterate to
            end = 1 + np.argmin(np.abs(target[-1]-initial))
        #--------------------------------------- step 4: ----------------------------------------------------------------- 


        #------------------------------------- step 5: find point(s) of closest approach/intersection ---------------------------------------- 
        if nice_spacing: 
            #if this is true, there are only two options: 1) they're in phase, ie they have elements that are equal to each 
            #other periodically by index, or 2) they are out of phase ie none of their elements are equal to one another

            #to test in/out-of-phase-ness, take a subset of the denser array & compare to the other 
            #if they're in phase(and we already know they overlap), there will be at least one exact element from the subset of the 
            #denser array that is equal to an element form the other 
            #if they are out of phase, none of the elements from the subset of the denser array will match an element from the other

            if need_interp: #target is the denser array & interpolation needed (d1 > d2)
                if target_starts_before_initial:
                    subset_start = np.argmin(np.abs(initial[0] - target))
                else:
                    subset_start = 0 

                denser = target #array from which to take subset
                sparser = initial
                subset_end = subset_start + d2d1 
            

            elif not need_interp: #initial is the denser array, interpolation not needed (d2 > d1)
                denser = initial
                sparser = target
                subset_start = start 


                if d2d1 == 1: #special case, in which we need to select d2 number of elements in subset
                    subset_end = start + int(d2)+1 #in case d2 is a decimal, +1 accounts for: int rounds down but need to be rounding up 
                else:
                    subset_end = start + d2d1


            subset = denser[subset_start : subset_end]


            for i, ith in enumerate(subset):
                compare = np.abs(ith - sparser)
                j = np.argmin(compare)

                if ith == sparser[j]: #the subset element is exactly equal to an element in the sparser array
                    i += subset_start #update i since it was only indexing within subset
                    in_phase = True
                    break

            if in_phase: #at some point in the for loop, we found an exact match
                #the exact match should be the first matching pair 
                #we also want to pin down the last matching pair

                subset_size = subset_end - subset_start

                #finally we can construct the new array

                if need_interp: 
                    print('still working on interpolatin bruv')

                else: #interpolation isn't needed -> ie initial = denser & target = sparser
                    subset_last = denser [-subset_size + end: end]
                    for k, kth in enumerate(subset_last):
                        compare = np.abs(kth - sparser)
                        m = np.argmin(compare)

                        if (kth == sparser[m]):
                            k += end - subset_size #update k, since before it was only indexing within subset_last
                            break


                    init_data_trans = initial_data[i: k+1 :d2d1] #subset of initial that matches with target, selecting only the nth elements

                    if target_starts_before_initial:
                        padding_start = np.array([np.nan] * len(target[:j])) #create an array of nans the same length as how many items are in padding before initial starts
                        init_data_trans = np.append(padding_start, init_data_trans) 


                    if target_ends_after_initial:
                        padding_end = np.array([np.nan]*len(target[m+1:])) #don't want to include mth element  
                        init_data_trans = np.append(init_data_trans, padding_end) #add the stuff thats at the end of target past initial's end

                        
                return init_data_trans
            
            else: #for loop completed without exact match => they're completely out of phase 
                print("they're out of phase, physically not sure what this would entail, but that's a problem for later")


        else: #not nice spacing
            #DECIDE IF THIS NEEDS TO ALSO BE BROKEN DOWN INTO INTERP NEEDED / NOT NEEDED, OR IF ONE FORMULA WILL HANDLE BOTH CASES
            print('not nice spacing')
            for i, init in enumerate(initial):
                j = np.argmin(np.abs(init - target))
        

        #if this is true, then we will need to know whether they are in phase or not (ie whether there are
        #elements in one array that are equal to elements in the other)       

        #find the first elements that match between the two arrays


        # #find last pair of matching elements
        # last_subset_init = initial[end-d2d1 : end]

        # for k, last_init in enumerate(last_subset_init):
        #     k = k - d2d1 + end

        #     compare = np.abs(last_init - target)
        #     m = np.argmin(compare)
        #     last_tar = target[m]
        #     if (last_init == last_tar):
        #         print('k is ', k, 'and end is ', end)
        #         #now we also have our last_init & last_tar, as well as their positions k & m
        #         break

            
        #--------------------------------------------------------------- step 5 --------------------------------------------------------- 

    
    else:
        print('the two arrays do not coincide physically, the transformation has no meaning')
 

#then go through the CO_cube, and transform the data points to new positions in the empty HI
for i in range(CO_cube):
    CO_cube[i] = vel_ax
    for j in range(vel_ax):
        vel_ax[j] = lat_ax
        for k in range(lat_ax):
            vbl_data = lat_ax[k]
            

HI_vels_kms = HI_vels*(u.m.to(u.km))


CO_vels[0], CO_vels[-1]


rms_CO = rms(CO_cube)


inter_mask = CO_cube < rms_CO


inter_mask.shape


plt.pcolormesh(inter_mask[200])


#plot surface density
cld.surface_density(LAB, 1.5, '21cm', show_Sigma=True)


Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)


orion = cld.extract_cloud(17, Zucker_df, DHT_CO, p_padding=60)


orion_cube = orion['cube']
orion_head = orion['header']
orion_units = orion['units']


plt.imshow(np.nansum(orion_cube, axis = 0), origin = 'lower')


orion_dend = Dendrogram.compute(orion_cube,\
                               min_value = dend_sigma,\
                               min_delta = 2*dend_sigma,\
                               min_npix = 3*dend_ppb)


orion_mtdt = {}
orion_mtdt['data_unit'] = u.Jy #needs to be u.Jy in order for this to work -- why not K ??
orion_catalog = ppv_catalog(orion_dend, orion_mtdt)


orion_dclust = SpectralCloudstering(orion_dend, orion_catalog, orion_head)


trunks_flat = np.nanmedian(orion_dclust.trunks_asgn.data, axis = 0)


trunks_flat.shape


orion_dclust.trunks_asgn.data.shape


plt.imshow(trunks_flat, origin = 'lower')
plt.colorbar()


plt.pcolormesh(orion_dclust.trunks_asgn.data[6])
plt.colorbar()

plt.show()


HI_dend = Dendrogram.compute(HI_cube, \
                                min_value=dend_sigma,\
                                min_delta=2*dend_sigma,\
                                min_npix=3*dend_ppb)


HI_mtdt = {} #metadata
HI_mtdt['data_unit'] = u.Jy
HI_catalog = ppv_catalog(HI_dend, HI_mtdt)


HI_dclust = SpectralCloudstering(HI_dend, HI_catalog, HI_head)


HI_dclust.clusters_asgn.data.shape


plt.pcolormesh(HI_dclust.clusters_asgn.data[80])


#following along tutorial (https://scimes.readthedocs.io/en/latest/tutorial.html)
#for calculating dendograms
orion_hdu = fits.open('DHT27_Orion_raw.fits')[0]
orion_data = orion_hdu.data
orion_head = orion_hdu.header


orion_dend = Dendrogram.compute(orion_data, \
                                min_value=dend_sigma,\
                                min_delta=2*dend_sigma,\
                                min_npix=3*dend_ppb)


metadata = {}
metadata['data_unit'] = u.Jy #hmm isn't it u.K ??
orion_cat = ppv_catalog(orion_dend, metadata)


orion_dclust = SpectralCloudstering(orion_dend, orion_cat, orion_head)

orion_dclust.showdendro()

plt.pcolormesh(orion_dclust.clusters_asgn.data[80])


# would be interesting/helpful to isolate the different velocity peaks & see how they map onto spatial features


HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = HI_0th
imx = HI_lons
imy = HI_lats


vmin,vmax = zlims(im) 


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow)


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')


plt.title('Galaxy in HI: 0th moment')
plt.show()


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
