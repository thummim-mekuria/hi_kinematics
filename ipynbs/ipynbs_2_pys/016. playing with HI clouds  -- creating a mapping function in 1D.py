#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


#-----------------------------------------------------------------------------------
two90 = np.arange(2, 90, 2)
reversetwo90 = np.flip(two90)

eighteen72 = np.arange(18,72,4)
reverse1872 = np.flip(eighteen72)
#-----------------------------------------------------------------------------------


transform_array(tens, twos)


-2%-10


-10%2


10%-2


tens = np.array((10, 20, 30, 40, 50))


two5s = np.array((7, 9, 11, 13, 15, 17, 19))


twos = np.array((6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36))


transform_array(two5s, twos)


d2 = twos[1] - twos[0]
d1 = tens[1] - tens[0]


d2d1 = d2/d1


CO_lats[-4]


len(CO_lats)


481-4


CO_lats[100:-4+1]


CO_lons[-1], HI_lons[-1]


two90[8:35]


eighteen72[13]


len(CO_cube), len(CO_cube[0]), len(CO_cube[0][0])


CO_lon_collapsed = np.nanmean(np.nanmean(CO_cube, axis = 0), axis=0)


HI_lon_collapsed = np.nanmean(np.nanmean(HI_cube, axis = 0), axis=0)


rms_CO = rms(CO_lon_collapsed)


mask_CO = CO_lon_collapsed > rms_CO


len(mask_CO)


mask_CO_transformed = transform_array(CO_lons, mask_CO, HI_lons)


inverted_mask_transformed = ~mask_CO_transformed


plt.subplots(figsize = (15,8))

# plt.plot(CO_lons, CO_lon_collapsed, label = 'CO orig')
# plt.plot(HI_lons, HI_lon_collapsed, label = 'HI')
# # plt.plot(HI_lons, CO_lons_transformed, label = 'CO transformed')
# plt.plot(CO_lons, mask_CO*CO_lon_collapsed, alpha = 0.5, label = 'mask')
# plt.hlines(rms_CO, -180, 180)


# plt.plot(CO_lon_collapsed, label = 'CO orig')
plt.plot(HI_lon_collapsed, label = 'HI')
plt.plot(inverted_mask_transformed*HI_lon_collapsed, label = 'mask transformed')


plt.plot(CO_lons_transformed, label = 'CO transformed')
# plt.plot(mask_CO, alpha = 0.5, label = 'mask orig')
plt.plot(mask_CO_transformed*CO_lons_transformed, label = 'mask transformed')
plt.hlines(rms_CO, 0, 700)


plt.legend()
plt.show


#YAY IT WORKS


CO_lons_transformed = transform_array(CO_lons, CO_lon_collapsed, HI_lons, HI_lon_collapsed)


transform_array(CO_lats, HI_lats)


thrts = np.array((1.3, 2.6, 3.9))


thrts[1] - thrts[0]


a = 2.6


b = 1.3


a%b


import datacube as dc
import clouds as cld

import numpy as np


from astropy import units as u
import astropy.constants as c


def transform_array(initial_coord, initial_data, target):

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


                    init_data_trans = initial_data[i:k+1 :d2d1] #subset of initial that matches with target, selecting only the nth elements

                    if target_starts_before_initial:
                        padding_start = np.zeros_like(target[:j])
                        init_data_trans = np.insert(init_data_trans, 0, padding_start) #add the stuff thats at the beginning of target before initial starts

                    if target_ends_after_initial:
                        padding_end = np.zeros_like(target[m:])
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


       #ALSO EVENTUALLY WILL NEED TO FIGURE OUT HOW TO HANDLE PAIRS OF ARRAYS - WHERE ONE TELLS YOU THE COORDINATE, AND THE OTHER TELLS YOU THE DATA
       #BECAUSE RIGHT NOW, WE'RE ASSUMING THE NUMBERS INSIDE THE ARAY ARE THE THING WE'RE STRETCHING AND SQUISHING, BUT ACTUALLY THE TRANSFORMATION MATTERS ON
       #THE COORDINATE, BUT THE ACTUAL THING WE NEED TO BE TRANSOFRMING IS THE DATAR -> THIS MIGHT ALSO PLAY INTO HOW INTERPOLATION WORKS (DO YOU LEAVE 0S WHERE THERE'S
       #NO DATA OR DO YOU FIND THE AVERAGE OR WOT)

    
    else:
        print('the two arrays do not coincide physically, the transformation has no meaning')
 

#works -----    all are the same sign, direction doesn't matter
transform_array(two90, eighteen72)            #+ascending, +ascending            same sign, same dir #int
transform_array(two90, reverse1872)           #+ascending, +descending           same sign, diff dir #int
transform_array(-two90, -eighteen72)          #-descending, -descending          same sign, same dir #int
transform_array(-two90, -reverse1872)         #-descending, -ascending           same sign, diff dir #int
transform_array(reversetwo90, eighteen72)     #+descending, +ascending           same sign, diff dir #int
transform_array(reversetwo90, reverse1872)    #+descending, +descending          same sign, same dir #int
transform_array(-reversetwo90, -eighteen72)   #-ascending, -descending           same sign, diff dir #int
transform_array(-reversetwo90, -reverse1872)  #-ascending, -ascending            same sign, same dir #int


# and all the times where it did work was when the two arrays did intersect
#bundle up all the times it didn't work as being the times where the two arrays weren't intersecing


#one extra element ----- all are diff sign, direction doesn't matter, +/- eighteen72
# all of these were flagged as target starts after initial 


transform_array(two90, -eighteen72)           #+ascending, -descending           diff sign, diff dir FALSE
transform_array(-two90, eighteen72)           #-descending, +ascending           diff sign, diff dir TRUE
transform_array(reversetwo90, -eighteen72)    #+descending, -descending          diff sign, same dir FALSE
transform_array(-reversetwo90, eighteen72)    #-ascending, +ascending            diff sign, same dir TRUE


#error didn't assign k ----- all are diff sign, direction doesn't matter, +/- reversed1872
# all of these were flagged as target starts before initial 

transform_array(two90, -reverse1872)          #+ascending, -ascending            diff sign, same dir TRUE
transform_array(-two90, reverse1872)          #-descending, +descending          diff sign, same dir FALSE
transform_array(reversetwo90, -reverse1872)   #+descending, -ascending           diff sign, diff dir TRUE
transform_array(-reversetwo90, reverse1872)   #-ascending, +descending           diff sign, diff dir FALSE


#conclusion -- whenever they're the same sign, it works
#whenever they're different signs, something goes wrong
#when its eighteen72 & -eighteen72 the thing that goes wrong is one more element
#when its reversed1872 & -reversed1872 the thing that goes wrong is k doesn't get assignd


#as expected -- says needs interpolating
transform_array(eighteen72, two90)
transform_array(-eighteen72, two90)
transform_array(reverse1872, two90)
transform_array(-reverse1872, two90)

#as expected -- says needs interpolating

transform_array(-eighteen72, two90)
transform_array(-eighteen72, -two90)
transform_array(-eighteen72, reversetwo90)
transform_array(-eighteen72, -reversetwo90)


#as expected -- says needs interpolating
transform_array(reverse1872, two90)
transform_array(reverse1872, -two90)
transform_array(reverse1872, reversetwo90)
transform_array(reverse1872, -reversetwo90)


#as expected -- says needs interpolating
transform_array(-reverse1872, two90)
transform_array(-reverse1872, -two90)
transform_array(-reverse1872, reversetwo90)
transform_array(-reverse1872, -reversetwo90)


#mix and match that all the permutations of two90, eighteen72, -two90, -eighteen72, 
#and their respective reverses give the correct transformation


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
