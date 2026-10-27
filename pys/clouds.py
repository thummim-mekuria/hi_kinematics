#===================================================================================
import numpy as np
import pandas as pd

import astropy.units as u 
import astropy.constants as c

import datacube as dc
from scipy.signal import peak_widths, find_peaks

from copy import deepcopy

import cloud_object as clob
#===================================================================================


#===================================================================================
def ind_ang (arrow, array, inp):
    """
    converts between indices & coordinate values along that axis(e.g., vel,lat,lon) 
    
    arrow:'get index' or 'get physical'
    array: eg velocity, latitude or longitude axes that come out of unpack_fits_cube
    
    num: an array of integers; paired with 'get physical'
         an array of floats/integers; paired with 'get index'
    """
    
    if arrow == 'get angular':
        output = array[inp]
  
   
    if arrow == 'get index':
        """
        1) define a buffer, i.e., when to decide that a physical value 
        is not in the bounds of array
        
        2) assuming all values in array are equally spaced, the buffer equals two 
        times this separation
            
        3) use the sign of delta to determine whether the array is increasing
         or decreasing then set bounds accordingly
        """
        output = np.zeros_like(inp) 
        
        delta = 2*(array[1] - array[0])
        if delta > 0:
            lowest  = array[0]  - delta
            highest = array[-1] + delta
            
        if delta < 0:
            lowest  = array[-1] - delta
            highest = array[0]  + delta
        
        
       
        for i, coord in enumerate(inp):
            #if inp is within these bounds, return its indices
            if (coord >= lowest) and (coord <= highest):
                diff = abs(array - coord)
                closest = np.argmin(diff)

                output[i] = closest

       
             #else, return the edge of array closest to num
            else: 
                if (coord <= lowest):
                    output[i] = np.argmin(array)
                elif coord >= highest:
                    output[i] = np.argmax(array)
            
            
        output = output.astype(int)



    return output
#===================================================================================


#===================================================================================
def surface_density (cube, data_type):
    '''
    parameters:
    1) a clob (=cloud object)  (the output from create_subcube/unpack_fits_cube)
            assume that clob.cube is moment-masked / is mostly signal
    2) string indicating whether the data is CO or 21cm 
            (to decide different which conversion factor to use)
    '''

    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 * c.m_p #2 bc molecule, 1.36 bc accounti for he 
        dv = 1.3 * u.km/u.s
        dA_px = (0.125)**2 * (u.deg**2).to(u.rad**2)

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) 
        #^ from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36 * c.m_p
        dv = 1.288 * u.km/u.s
        # dA_px = (0.125)**2 * (u.deg**2).to(u.rad**2) #based on dl/db for HI4PI

    W = np.nansum(cube * u.K * dv, axis=0) 
    N = W * conversion_factor 
    dSig = mmw * N 
    
    return dSig,dA_px
#===================================================================================


#===================================================================================
def mass(dSig, d, ang_area):
    """
    """
    Sig = np.nansum(dSig)

    #UGHH so when ang_area plugged in here is just the area of one pixel (db*dl),then
    #the calculated mass is correct, but if the ang_area plugged in is the whole area
    #of the cloud then it doesnt => whyyyyyyy (in my brain -- shouldnt Sigma = M/A, so
    #then sum(Sigma) = sum(M)/sum(A??))
    
    a_phys =  (d**2) * ang_area #to get rid of 'rad' as a unit
    
    m_total = (Sig * a_phys).to(u.kg)

    return m_total
#===================================================================================


#===================================================================================
def accumulation_radius(m, d, b, rho_HI = 0.025*c.M_sun/(u.pc**3)):
    """
    calculated using equation 12 from imara & blitz angular momentum paper
    """
    #in physical units
    r_a = ((m/(2*np.pi*rho_HI))**(1/3)).to(u.pc)

    #in sky units
    s = ((r_a/d) * u.rad).to(u.deg) * np.cos(b.to(u.rad))
    return r_a, s
#===================================================================================


#===================================================================================
def accumulation_prisms (vs,bs,ls,r_as,v_rmss,clob):
    big_mask = np.zeros_like(clob.cube, dtype=bool) 
    for v,b,l,r,vrms in zip(vs,bs,ls,r_as,v_rmss):
        #4 sigma wide, to recover > 99.7% of the velocity channels with the spectrum
        v0,v1 = v-(2*vrms), v+(2*vrms)

        #square with l=2xr_a to enclose area ~ circle defined by accumulation radius 
        b0,b1 = b-r, b+r 
        l0,l1 = l-r, l+r
        # print(l0,l1,b0,b1,v0,v1)

        #turn physical locations into index locations
        ls = ind_ang('get index', clob.lons,[l0,l1])
        bs = ind_ang('get index', clob.lats,[b0,b1])
        vs = ind_ang('get index', clob.vels/1e3,[v0, v1]) #bc alr know about m/s->km/s
        
        ls.sort()
        bs.sort()
        vs.sort()
        
        l0,l1 = ls
        b0,b1 = bs
        v0,v1 = vs
        # print(l0,l1,b0,b1,v0,v1)
        # print(big_mask.shape)
        small_mask = np.ones_like(big_mask[v0:v1,b0:b1,l0:l1],dtype=bool)
        # print(small_mask)

        # print(np.all(small_mask==True))
        big_mask[v0:v1,b0:b1,l0:l1] = small_mask
        # print(big_mask[v0:v1,b0:b1,l0:l1])
        # print('\n\n\n')
        # print(np.all(big_mask==False))

        
    return big_mask
#===================================================================================

