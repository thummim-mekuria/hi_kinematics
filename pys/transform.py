
import datacube as dc
import clouds as cld

import numpy as np


from astropy import units as u
import astropy.constants as c

import numpy as np



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

        #--------------------------------------------------------------- step 5 --------------------------------------------------------- 

        




        #----------------------------------------------------------- step 6: finally we can construct the new array ----------------------------------------------------------- 

                
                if need_interp: 
                    #so now we have the following information: i & k (= positions of first/last matching pair in denser),
                    #j & m (= positions of first/last matching pair in sparser), d2d1 -> factor by which desnser has more 
                    #elements (ie how many spaces in denser between the elements that match every single one in sparser)
                    
                    #create an array thats the same length as how many indices are between the first & 
                    #last matching pairs, and that is as wide as how many interpolations will need to happen
                    init_data_trans = np.zeros(((m-j +1),d2d1)) 
                    for s in range(0, m-j +1):
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


        #--------------------------------------------------------------------------------- step 6 -------------------------------------------------------------------------------------- 

            
            else: #for loop completed without exact match => they're completely out of phase 
                print("they're out of phase, physically not sure what this would entail, but that's a problem for later")




        else: #not nice spacing
            #DECIDE IF THIS NEEDS TO ALSO BE BROKEN DOWN INTO INTERP NEEDED / NOT NEEDED, OR IF ONE FORMULA WILL HANDLE BOTH CASES
            print('not nice spacing')
            print('d1 is',d1,'d2 is',d2)
            for i, init in enumerate(initial):
                j = np.argmin(np.abs(init - target))
        

       

    
    else:
        print('the two arrays do not coincide physically, the transformation has no meaning')
 





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

    #AS IT IS, THIS ONLY WORKS FOR DATA THAT DOESN'T NEED INTERPOLATION (INITIAL IS DENSER THAN TARGET)
    #BECAUSE IF INTERPOLATION IS NEEDED, WILL HAVE TO DO MORE MATH -- NOT TACKLING THAT FOR NOW THO





