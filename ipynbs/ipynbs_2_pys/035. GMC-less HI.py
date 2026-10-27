#!/usr/bin/env python
# coding: utf-8


# %%time

import datacube as dc
import clouds as cld

import numpy as np
# %matplotlib

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


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf

DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')
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
        ns = ['03','06','09','12','16','19','21','24','26','27','34','37']
        #for now, figure out DHT01 (bad data) 
        #& DHT18 (two ceph ds, beeg) 
        #& DHT08 (beeg) 
        #& DHT24, DHT34 (lam or & chamaeleon are disjointed down to the really low sig levels)
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output

#----------------------

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


cat


DHT_CO.create_subcube([207,218,])


import transform as tr


tr.transform_2D


HI4PI.lons[0], HI4PI.lons[-1]


HI4PI.lats[0], HI4PI.lats[-1]


DHT_CO.lats[0], DHT_CO.lats[-1]


np.argmin(np.abs(HI4PI.vels-(37000)))


HI4PI.lats[718:723]


HI4PI.lats[1440:1444]


b0DHT, b1DHT = 600, 1600
v0DHT, v1DHT = 400, 550


HI4PI.vels[400]/1e3,HI4PI.vels[550]/1e3


1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]


DHT_CO.vels[0], DHT_CO.vels[-1]


DHT_lats_vels = HI4PI.clone()


DHT_lats_vels.cube = HI4PI.cube[v0DHT:v1DHT,b0DHT:b1DHT,:]


DHT_lats_vels.lats = HI4PI.lats[b0DHT:b1DHT]
DHT_lats_vels.vels = HI4PI.vels[v0DHT:v1DHT]


y = DHT_lats_vels.lats
x = DHT_lats_vels.lons
z = DHT_lats_vels.vels


imHI = np.nansum(DHT_lats_vels.cube, axis=0)


plt.imshow(imHI, extent=[x[0],x[-1],y[0],y[-1]])
plt.show()


DHTdSig, DHTdA = cld.surface_density(DHT_CO,'CO',show_Sigma=False)


DHTdN = (DHTdSig/(c.m_p * 2 * 1.36)).cgs


plt.hist(DHTdN.flat)
plt.xscale('log')
plt.yscale('log')
plt.vlines(1e21,0,1e6,color='red')


plt.figure(figsize=(10,10))
plt.imshow(DHTdN.value * (DHTdN.value > 1e21) ,**kw_imshow)
# plt.imshow(DHTdN.value,**kw_imshow)

plt.show()


CO_mask = DHTdN > 1e21/u.cm**2


initial_data, initial_ax0_coord, initial_ax1_coord, target_ax0_coord, target_ax1_coord


LABHI = dc.unpack_fits_cube('labh.fits','LAB')


LABHI.cube.shape


DHT_CO.cube.shape


HImask = transform_2D(~binary_dilation(CO_mask,iterations=4),DHT_CO.lats,DHT_CO.lons, LABHI.lats, LABHI.lons)


from scipy.ndimage import binary_dilation


newm = binary_dilation(HImask)


np.nanmax(HImask)


HImask.shape


DHT_lats_vels.cube.shape


HIdSig, HIdA = cld.surface_density(DHT_lats_vels,'21cm',show_Sigma=False)


HIdSig = HIdSig.to(c.M_sun/u.pc**2)


from matplotlib.patches import Rectangle


HImask


np.argmin(np.abs(LABHI.vels-(-10e3)))


np.argmin(np.abs(LABHI.vels-(22e3)))


newim = np.nansum(LABHI.cube[436:467,:,:],axis=0)


fig,ax = plt.subplots(2,1,figsize=(20,10))
ax[0].imshow(HIdSig.value,**kw_imshow,extent=[x[0],x[-1],y[0],y[-1]])
# ax.colorbar()
ax[0].contour(x,y,HIdSig.value,levels=[10],colors=['white'],linewidths=0.8)
# ax[0].contour(x,y,HIdSig.value,levels=[100],colors=['white'],linewidths=0.8)

ax[0].contour(LABHI.lons,LABHI.lats,HImask,linewidths=0.7,colors='red',levels=[0.0])


ax[0].add_patch(Rectangle((55,5),20,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-120,-12),20,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-40,-10),20,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((78,-15),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-105,5),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-65,5),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((167,10),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-90,-15),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((128,-13),10,5,alpha=0.5,color='orange'))


ax[0].hlines(5,100,-100,color='orange',linestyle='--')
ax[0].hlines(-5,100,-100,color='orange',linestyle='--')
ax[0].vlines(100,-5,5,color='orange',linestyle='--')

("")

# plt.xlim(180,100)
ax[0].set_ylim(-30,30)
ax[0].set_xlim(180,-10)
ax[0].set_ylabel(r'galactic latitude [$\degree$]')
ax[0].set_xlabel(r'galactic longitude [$\degree$]')


#----------------------------------------------------------------------------
ax[1].imshow(HIdSig.value,**kw_imshow,extent=[x[0],x[-1],y[0],y[-1]])
# ax.colorbar()
ax[1].contour(x,y,HIdSig.value,levels=[10],colors=['white'],linewidths=0.8)
# ax[1].contour(x,y,HIdSig.value,levels=[100],colors=['white'],linewidths=0.8)

ax[1].contour(LABHI.lons,LABHI.lats,HImask,linewidths=0.7,colors='red',levels=[0.0])

ax[1].add_patch(Rectangle((55,5),20,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-130,-12),30,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-40,-10),20,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((78,-15),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-105,6),18,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-65,5),30,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((167,10),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-90,-15),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((128,-13),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-145,2),35,5,alpha=0.5,color='orange'))


ax[1].hlines(5,100,-100,color='orange',linestyle='--')
ax[1].hlines(-5,100,-100,color='orange',linestyle='--')
ax[1].vlines(-100,-5,5,color='orange',linestyle='--')
#solar neighbors
ax[1].add_patch(Rectangle((115,20),15,10,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-70,-20),20,10,color='cyan',fill=False,linestyle='--'))
ax[1].add_patch(Rectangle((150,-25),30,20,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-160,-25),25,20,color='cyan',fill=False))
ax[1].add_patch(Rectangle((40,5),10,10,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-10,10),10,10,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-5,-25),10,10,color='cyan',fill=False,linestyle='--'))
ax[1].add_patch(Rectangle((-175,-20),15,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((90,-20),10,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-165,-5),15,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((25,-25),15,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((145,15),5,5,color='cyan',fill=False))


# plt.xlim(180,100)
ax[1].set_ylim(-30,30)
ax[1].set_xlim(-180,10)
ax[1].set_ylabel(r'galactic latitude [$\degree$]')
ax[1].set_xlabel(r'galactic longitude [$\degree$]')


plt.show()


# moving forward after this ^:
# - do the 0th moment on a narrower velocity range, look at how that affects where the CO & HI thresholds are
# - compare dendrogram vs using the m33 approach for the orange boxes (using median accumulation radii)
# - start of the averaged-over-area spectra, but might be useful to go back & do pixel-by-pixel
# - send nia papers on the HISA
# - compare the HI masses, surface densities etc by calculating using 1) kinematic distance & 2) zucker distances


newim


fig,ax = plt.subplots(2,1,figsize=(20,10))
ax[0].imshow(HIdSig.value,**kw_imshow,extent=[x[0],x[-1],y[0],y[-1]])
# ax.colorbar()
ax[0].contour(x,y,newim)#,levels=[10],colors=['white'],linewidths=0.8)
# ax[0].contour(x,y,HIdSig.value,levels=[100],colors=['white'],linewidths=0.8)

ax[0].contour(LABHI.lons,LABHI.lats,HImask,linewidths=0.7,colors='red',levels=[0.0])


ax[0].add_patch(Rectangle((55,5),20,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-120,-12),20,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-40,-10),20,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((78,-15),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-105,5),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-65,5),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((167,10),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((-90,-15),10,5,alpha=0.5,color='orange'))
ax[0].add_patch(Rectangle((128,-13),10,5,alpha=0.5,color='orange'))


ax[0].hlines(5,100,-100,color='orange',linestyle='--')
ax[0].hlines(-5,100,-100,color='orange',linestyle='--')
ax[0].vlines(100,-5,5,color='orange',linestyle='--')

#solar neighbors
ax[0].add_patch(Rectangle((115,20),15,10,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((-70,-20),20,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((150,-25),30,20,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-160,-25),20,20,color='cyan',fill=False))
ax[0].add_patch(Rectangle((40,5),10,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-10,10),10,10,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-5,-25),10,10,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((-175,-20),15,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((90,-20),10,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((-165,-5),15,15,color='cyan',fill=False))
ax[0].add_patch(Rectangle((25,-25),15,15,color='cyan',fill=False,linestyle='--'))
ax[0].add_patch(Rectangle((145,15),5,5,color='cyan',fill=False,linestyle='--'))


# plt.xlim(180,100)
ax[0].set_ylim(-30,30)
ax[0].set_xlim(180,-10)
ax[0].set_ylabel(r'galactic latitude [$\degree$]')
ax[0].set_xlabel(r'galactic longitude [$\degree$]')


#----------------------------------------------------------------------------
ax[1].imshow(HIdSig.value,**kw_imshow,extent=[x[0],x[-1],y[0],y[-1]])
# ax.colorbar()
ax[1].contour(x,y,HIdSig.value,levels=[10],colors=['white'],linewidths=0.8)
# ax[1].contour(x,y,HIdSig.value,levels=[100],colors=['white'],linewidths=0.8)

ax[1].contour(LABHI.lons,LABHI.lats,HImask,linewidths=0.7,colors='red',levels=[0.0])

ax[1].add_patch(Rectangle((55,5),20,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-130,-12),30,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-40,-10),20,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((78,-15),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-105,6),18,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-65,5),30,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((167,10),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-90,-15),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((128,-13),10,5,alpha=0.5,color='orange'))
ax[1].add_patch(Rectangle((-145,2),35,5,alpha=0.5,color='orange'))


ax[1].hlines(5,100,-100,color='orange',linestyle='--')
ax[1].hlines(-5,100,-100,color='orange',linestyle='--')
ax[1].vlines(-100,-5,5,color='orange',linestyle='--')
#solar neighbors
ax[1].add_patch(Rectangle((115,20),15,10,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-70,-20),20,10,color='cyan',fill=False,linestyle='--'))
ax[1].add_patch(Rectangle((150,-25),30,20,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-160,-25),25,20,color='cyan',fill=False))
ax[1].add_patch(Rectangle((40,5),10,10,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-10,10),10,10,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-5,-25),10,10,color='cyan',fill=False,linestyle='--'))
ax[1].add_patch(Rectangle((-175,-20),15,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((90,-20),10,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((-165,-5),15,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((25,-25),15,15,color='cyan',fill=False))
ax[1].add_patch(Rectangle((145,15),5,5,color='cyan',fill=False))


# plt.xlim(180,100)
ax[1].set_ylim(-30,30)
ax[1].set_xlim(-180,10)
ax[1].set_ylabel(r'galactic latitude [$\degree$]')
ax[1].set_xlabel(r'galactic longitude [$\degree$]')


plt.show()


#a test to check


DHT21


-360+194.5


cat


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


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
    
#     print('created empt',initial_data_trans.shape)

    #then go through and populate the new array with elements selected from initial_data 
    #the data are selected according to the transformed coordinate arrays

    #AS IT IS, THIS ONLY WORKS FOR DATA THAT DOESN'T NEED INTERPOLATION (INITIAL IS DENSER THAN TARGET)
    #BECAUSE IF INTERPOLATION IS NEEDED, WILL HAVE TO DO MORE MATH -- NOT TACKLING THAT FOR NOW THO
#     print(init_ax0_ins_trans)#[0], init_ax0_ins_trans[-1])

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
 

DHT27 = surveys['DHT27']
#DHT27[0].calculate_rms_noise() is 0.4888


d27 = Dendrogram.compute(DHT27[0].cube,5*0.49,0.49,100)


# fig,ax = plt.subplots(1,1,figsize=(8,8))
# d27.plotter().plot_tree(ax)


bigm = d27.trunk[0].get_mask()
for t in d27.trunk:
    bigm += t.get_mask()


fig,axs = plt.subplots(1,2,figsize=(10,10))
lb0 = np.nansum(bigm*DHT27[0].cube,axis=0)
lv0 = np.nansum(bigm*DHT27[0].cube,axis=1)

lb0lims = DHT27[0].lons[0], DHT27[0].lons[-1],DHT27[0].lats[0],DHT27[0].lats[-1]
lv0lims = DHT27[0].lons[0], DHT27[0].lons[-1],DHT27[0].vels[0],DHT27[0].vels[-1]

axs[0].imshow(lb0,extent=lb0lims,**kw_imshow)
axs[0].set_xlim(222,200)
axs[0].set_xlabel(r'Galactic Longitude [$\degree$]')
axs[0].set_ylabel(r'Galactic Latitude [$\degree$]')

axs[1].imshow(lv0,extent=lv0lims,**kw_imshow)
axs[1].set_xlabel(r'Galactic Longitude [$\degree$]')
axs[1].set_ylabel(r'$V_{LSR}$ [$km\ s^{-1}$]')
axs[1].set_ylim(-5,20)
axs[1].set_xlim(222,200)

plt.show()


orb = d27.trunk[0].get_mask()*DHT27[0].cube
orbclob = DHT27[0].clone()
orbclob.cube = orb


ora = d27.trunk[1].get_mask()*DHT27[0].cube
oraclob = DHT27[0].clone()
oraclob.cube = ora


mr2 = d27.trunk[6].get_mask()*DHT27[0].cube
mr2clob = DHT27[0].clone()
mr2clob.cube = mr2


crs = d27.trunk[7].get_mask()*DHT27[0].cube
crsclob = DHT27[0].clone()
crsclob.cube = crs


d27s = [crsclob,mr2clob,oraclob,orbclob]
ns = DHT27[1]


M27s = np.zeros_like(d27s)
R27s = np.zeros_like(d27s)
R_ac_27s = np.zeros_like(d27s)
s_ac_27s = np.zeros_like(d27s)


for cl, cl27 in enumerate(d27s):
    dSig, dA = cld.surface_density(cl27, 'CO',show_Sigma=False)
    M,A = cld.mass_and_area(ns[cl],cat,dSig,dA)
    RA,s = cld.calculate_accumulation_radius(M,cat['d__pc'][ns[cl]]*u.pc, cat['b_c__deg'][ns[cl]])
    M27s[cl] = M.to(c.M_sun)
    R27s[cl] = np.sqrt(A/np.pi)
    R_ac_27s[cl] = RA
    s_ac_27s[cl] = s


l0i27,l1i27 = cld.ind_ang('get index',HI4PI.lons,[-360+DHT27[0].lons[0],-360+DHT27[0].lons[-1]])
b0i27,b1i27 = cld.ind_ang('get index',HI4PI.lats,[DHT27[0].lats[0],DHT27[0].lats[-1]])
v0i27,v1i27 = cld.ind_ang('get index',1e-3*HI4PI.vels,[DHT27[0].vels[0],DHT27[0].vels[-1]])


HI27c = HI4PI.cube[v0i27:v1i27, b0i27:b1i27, l0i27:l1i27]


from matplotlib.patches import Circle,Rectangle


colors27 = ['#caa0ff', '#bbcc00', '#ffb7c8','#ffcc00']


fig,ax = plt.subplots(1,1,figsize=(10,10))

hilims = [HI4PI.lons[l0i27]+360, 360+HI4PI.lons[l1i27], HI4PI.lats[b0i27], HI4PI.lats[b1i27]]
ax.imshow(np.nansum(HI27c*HI4PI.dv.value/1e3, axis=0),**kw_imshow,extent=hilims)

for n in range(4):
    clob = d27s[n]
    cube = clob.cube
    ncat = ns[n]
    ax.contour(d27s[n].lons,clob.lats,np.nansum(cube*clob.dv.value,axis=0),levels=[1],colors=colors27[n])
    ax.add_patch(Circle((cat['l_c__deg'][ncat], cat['b_c__deg'][ncat]),\
                    s_ac_27s[n],fill=False,edgecolor=colors27[n]))

    
ax.add_patch(Rectangle((-149+360,-16),2,2,color="#ff0000",fill=False))


overlap = HI27c[:,74:98,157:181] # the red square
overlap_vs = HI4PI.vels[v0i27:v1i27]


overlap_T = overlap.transpose().swapaxes(0,1)


plt.imshow(np.nansum(overlap*HI4PI.dv.value/1e3,axis=0),**kw_imshow,extent = [-147,-149,-16,-14])
plt.colorbar(label=r'$W_{HI}$ [$K\ km\ s^{-1}$]')
plt.xticks(rotation='vertical')
plt.ylabel(r'b [$\degree$]')
plt.xlabel(r'l [$\degree$]')

plt.show()


pk_vels = np.zeros(4)


plt.figure(figsize=(10,10))
plt.imshow(np.nansum(overlap*HI4PI.db.value,axis=1),**kw_imshow,extent=[-147,-149,-1.02,3.87])
for p, pk in enumerate(pk_vels):
    plt.hlines(pk/10, -147, -149,color=colors27[p],linestyle='--',label=cat['name'][ns[p]])
plt.colorbar(label=r'$W_{HI}$ [$K\ \degree$]')
plt.xlabel(r'l [$\degree$]')
plt.ylabel(r'$V_{LSR}$ [$10\times km\ s^{-1}$]')
plt.legend()
plt.show()


plt.figure(figsize=(10,10))
plt.plot(overlap_vs/1e3,dc.spectrum1D(overlap),label='averaging the whole square')

for n in range(4):
    sp = dc.spectrum1D(d27s[n].cube)
    plt.plot(d27s[n].vels,5e2*sp,label=cat['name'][ns[n]],color=colors27[n])
    pk = d27s[n].vels[np.nanargmax(sp)]
    pk_vels[n] = pk
    plt.vlines(pk,\
               0,5e2*np.nanmax(sp),color=colors27[n],linestyle='--')
    
    
plt.xlabel(r'$V_{LSR}$ [$km\ s^{-1}$]')
plt.ylabel(r'$T_B$ [$K$]')
plt.legend()


plt.show()


for i in range(0,24,12):
    for k in range(0,24,12):
        plt.plot(overlap_vs/1e3,overlap_T[i][k],alpha=0.4,color='grey')

for n in range(len(ns)):
    plt.plot(d27s[n].vels,5e2*dc.spectrum1D(d27s[n].cube),label=cat['name'][ns[n]],color=colors27[n])
    
plt.xlabel(r'$V_{LSR}$ [$km\ s^{-1}$]')
plt.ylabel(r'$T_B$ [$K$]')
plt.legend()


plt.show()


fig,axs=plt.subplots(4,4,figsize=(10,10))

for i in range(4):
    for k in range(4):
        axs[i,k].plot(overlap_vs/1e3,overlap_T[2*k+i*3][k*6],color='grey')

for ax in axs.flat:
        for n in range(len(ns)):
            ax.plot(d27s[n].vels,5e2*dc.spectrum1D(d27s[n].cube),\
                    label=cat['name'][ns[n]],alpha=0.5,color=colors27[n])
    
# plt.xlabel(r'$V_{LSR}$ [$km\ s^{-1}$]')
# plt.ylabel(r'$T_B$ [$K$]')
# plt.legend()

# for ax in fig.get_axes():
#     ax.label_outer()
plt.show()


dover = Dendrogram.compute(overlap,)


fig,ax = plt.subplots(1,1,figsize=(8,8))
dover.plotter().plot_tree(ax)


M27s


R27s


R_ac_27s


plt.imshow(np.nansum(d27.trunk[8].get_mask()*DHT27[0].cube,axis=0),**kw_imshow)


fig,ax = plt.subplots(1,1,figsize=(10,10))
hi0 = np.nansum(hicube*HI4PI.dv.value,axis=0)
limshi = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
vnhi,vxhi = zlims(hi0)
plt.imshow(hi0,**kw_imshow,vmin=vnhi,vmax=vxhi,extent=limshi)
plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],co0,levels=[4],colors=['white'])
ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2,label='crossbones'))
ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2,label='orion a'))
ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2,label='orion b'))
ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2,label='mon r2'))
plt.legend()
plt.show()


for t in hod.trunk[2].descendants:
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    limshi = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
#     if t.is_leaf:
    ax.imshow(np.nansum(t.get_mask()*hicube,axis=0),**kw_imshow,extent=limshi)
    ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
    ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
    ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
    ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))

    plt.show()


# calculated accumulation radii & got:
# crossbones = 4 degrees
# mon_r2 = 7 degrees
# orion_a = 11 degrees
# orion_b = 11 degrees


from matplotlib.patches import Circle
fig,ax = plt.subplots(1,1,figsize=(10,10))
limsco = [DHT_CO.lons[2500], DHT_CO.lons[2750],DHT_CO.lats[50],DHT_CO.lats[220]]
co0 = np.nansum(cocube*DHT_CO.dv.value,axis=0)
vnco,vxco=zlims(co0)
ax.imshow(co0,**kw_imshow,extent=limsco,vmin=vnco,vmax=vxco)
# plt.colorbar()
ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
ax.add_patch(Circle((-149,-19),11,edgecolor='white',fill=False,linewidth=2))
ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))

plt.show()


206-360


cod = Dendrogram.compute(cocube,cos,cos,100)


fig,ax = plt.subplots(1,1,figsize=(8,8))
cod.plotter().plot_tree(ax)


hod = Dendrogram.compute(hicube,5*hos,hos,1000)


fig,ax = plt.subplots(1,1,figsize=(8,8))
hod.plotter().plot_tree(ax)


cat


oraclob = DHT_CO.create_subcube([-21,-17,-143,-155])


dSigora, dAora = cld.surface_density(oraclob,'CO')


Mora, Aora = cld.mass_and_area(11,cat,dSigora,dAora)


Mora, Aora


Morb, Aorb


cld.calculate_accumulation_radius(Mora, cat['d__pc'][11]*u.pc)


plt.imshow(oraclob.calculate_0th_moment().value,**kw_imshow)


cat['l_c__deg'][11]-360-6


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[7],colors='red')
plt.scatter(cat['l_c__deg'][11]-360,cat['b_c__deg'][11],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][11]-360+11,cat['l_c__deg'][11]-360-11)
plt.ylim(cat['b_c__deg'][11]-11,cat['b_c__deg'][11]+11)


orbclob = DHT_CO.create_subcube([-18,-10,-150,-158])


dSigorb, dAorb = cld.surface_density(orbclob,'CO')


Morb, Aorb = cld.mass_and_area(12,cat,dSigorb,dAorb)


Morb


cld.calculate_accumulation_radius(Morb,cat['d__pc'][12]*u.pc)


cat


mr2clob = DHT_CO.create_subcube([-15, -10,-143,-149])


dSigmr2, dAmr2 = cld.surface_density(mr2clob,'CO')


Mmr2, Amr2 = cld.mass_and_area(9, cat, dSigmr2,dAmr2)


cld.calculate_accumulation_radius(Mmr2,cat['d__pc'][9]*u.pc)


cat['b_c__deg'][9]-3


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[4.5],colors='red')
plt.scatter(cat['l_c__deg'][9]-360,cat['b_c__deg'][9],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][9]-360+7,cat['l_c__deg'][9]-360-7)
plt.ylim(cat['b_c__deg'][9]-7,cat['b_c__deg'][9]+7)


219-360


crsclob = DHT_CO.create_subcube([-13,-6,-138,-142])


dSigcrs, dAcrs = cld.surface_density(crsclob,'CO')


Mcrs, Acrs = cld.mass_and_area(5, cat, dSigcrs, dAcrs)


cld.calculate_accumulation_radius(Mcrs,cat['d__pc'][5]*u.pc)


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[5],colors='red')
plt.xlim(-135,-160)


hicube.shape


cocube.shape


DHT_CO.lats[50],DHT_CO.lats[220]


for i in range(0,200,5):
    plt.pcolormesh(hicube[i],cmap=cmap0)
#     plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
#                 np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[7],colors='red')
    plt.show()


#crossbones - 886pc
#mpn_r2 - 778pc
#orion_a - 432pc
#orion_b - 423pc


hicube_T = hicube.transpose().swapaxes(0,1)


from scipy.signal import find_peaks

find_peaks(dc.spectrum1D(hicube[:,200:210,300:310]),height=25)


plt.plot(dc.spectrum1D(hicube[:,200:210,300:310]))


specs = np.zeros((490,700))
for b,latslc in enumerate(hicube_T):
    for l,lonslc in enumerate(latslc):
        H = 10
        p = find_peaks(lonslc,height=H)
        nps = len(p[0])
        specs[b][l] = nps


HI4PI.lats[600],HI4PI.lats[1050]


1050-600


plt.figure(figsize=(10,10))
lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
plt.imshow(np.nansum(hicube_T,axis=2),origin='lower',cmap=cmap0,extent=lims)
plt.colorbar()
# plt.ylim(0,450)
plt.show()


hicube = HI4PI.cube[400:600,600:1090,3600:4300]


HI4PI.lons[4300]


cut2 = hicube[:,370:430,180:300].transpose().swapaxes(0,1)


cut2.shape


plt.plot(HI4PI.vels[400:600]/1e3,dc.spectrum1D(hicube[:,370:400,180:190]))
plt.plot(oraclob.vels,10*dc.spectrum1D(oraclob.cube),label='ora')
plt.plot(orbclob.vels,10*dc.spectrum1D(orbclob.cube),label='orb')
plt.plot(mr2clob.vels,10*dc.spectrum1D(mr2clob.cube),label='mr2')
plt.plot(crsclob.vels,10*dc.spectrum1D(crsclob.cube),label='crs')
plt.xlim(-20,20)


for k in range(0,60,5):
    for m in range(0,120,10):
        plt.plot(HI4PI.vels[400:600]/1e3,cut2[k][m],label='this px')
        plt.plot(oraclob.vels,10*dc.spectrum1D(oraclob.cube),label='ora')
        plt.plot(orbclob.vels,10*dc.spectrum1D(orbclob.cube),label='orb')
        plt.plot(mr2clob.vels,10*dc.spectrum1D(mr2clob.cube),label='mr2')
        plt.plot(crsclob.vels,10*dc.spectrum1D(crsclob.cube),label='crs')
        plt.legend()
        plt.xlim(-10,25)
        plt.show()


fig,ax = plt.subplots(figsize=(10,10))
# lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
plt.imshow(specs[370:430,180:300],origin='lower',cmap='Paired')#,extent=lims)
plt.colorbar()
# ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
# ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
# ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
# ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))
plt.show()


fig,ax = plt.subplots(figsize=(10,10))
lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
ax.imshow(specs,origin='lower',cmap='Paired',extent=lims)
# plt.colorbar()
ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))
# plt.xlim(-135,-145)
# plt.ylim(-10,-6)
plt.show()


cutout = hicube[:,100:301,300:502]


cutout_T = cutout.transpose().swapaxes(0,1)


plt.plot(dc.spectrum1D(oraclob.cube))


# when coming back -> look more at promise of HISA? esp where multiple peaks are reported
# also, look at taurus/california/perseus dht survey to see if similar pattern (CO peaks at the valleys of HI)
# is replicated
# look at HISA papers some more, look into gaussian decomposition
# also -- look at physical sizes of accumulation radii and if it can be justified to smoosh orions / mon / xbns tgt?
# ^ or if some are actually separate, do axis=1 np.nansum on the HI dendrogram to see how the leafs spread out in v
# ^ or also, gaussian decomposition??


for k in range(0,200,20):
    for m in range(0,200,20):
        plt.plot(HI4PI.vels[400:600]/1e3,cutout_T[k][m],label='this px')
        plt.plot(oraclob.vels,10*dc.spectrum1D(oraclob.cube),label='ora')
        plt.plot(orbclob.vels,10*dc.spectrum1D(orbclob.cube),label='orb')
        plt.plot(mr2clob.vels,10*dc.spectrum1D(mr2clob.cube),label='mr2')
        plt.plot(crsclob.vels,10*dc.spectrum1D(crsclob.cube),label='crs')
        plt.legend()
        plt.xlim(-10,25)
        plt.show()


cutout_T.shape


hicube.shape


fig,ax = plt.subplots(figsize=(10,10))
# lims = [HI4PI.lons[3600],HI4PI.lons[4300],HI4PI.lats[600],HI4PI.lats[1090]]
plt.imshow(specs[100:301,300:502],origin='lower',cmap='Paired')#,extent=lims)
plt.colorbar()
# ax.add_patch(Circle((-141,-9),4,edgecolor='red',fill=False,linewidth=2))
# ax.add_patch(Circle((-149,-19),11,edgecolor='blue',fill=False,linewidth=2))
# ax.add_patch(Circle((-154,-14),11,edgecolor='orange',fill=False,linewidth=2))
# ax.add_patch(Circle((-146,-12.5),7,edgecolor='black',fill=False,linewidth=2))
plt.show()


plt.imshow(np.nansum(hicube,axis=0),origin='lower')


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[2.5],colors='red')
plt.scatter(cat['l_c__deg'][5]-360,cat['b_c__deg'][5],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][5]-360+4,cat['l_c__deg'][5]-360-4)
plt.ylim(cat['b_c__deg'][5]-4,cat['b_c__deg'][5]+4)


plt.contour(DHT_CO.lons[2500:2750],DHT_CO.lats[50:220],\
                np.nansum(cocube*DHT_CO.dv.value,axis=0),levels=[7],colors='red')
plt.scatter(cat['l_c__deg'][12]-360,cat['b_c__deg'][12],marker='+',color='blue')
plt.xlim(cat['l_c__deg'][12]-360+11,cat['l_c__deg'][12]-360-11)
plt.ylim(cat['b_c__deg'][12]-11,cat['b_c__deg'][12]+11)
