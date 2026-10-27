
import numpy as np
import pandas as pd
import astropy.units as u 
import astropy.constants as c
import matplotlib.pyplot as plt
import datacube as dc
from scipy.signal import peak_widths, find_peaks

from copy import deepcopy

import cloud_object as clob



#values from Russeil et al 2017
Theta_sun = 240*u.km/u.s
R_sun = 8.34*1000*u.pc
A1 = 1.022
A2 = 0.0803
A3 = 0.


# #values from Brand & Blitz 1993
# Theta_sun = 220*u.km/u.s
# R_sun = 8.5*1000*u.pc
# A1 = 1.00767
# A2 = 0.0394
# A3 = 0.00712



def R_poly(R, l, b, v_lsr, theta_sun=Theta_sun, r_sun=R_sun, a1=A1, a2=A2, a3=A3):

    term1 = R/r_sun
    term2 = (R/r_sun)**a2
    
    coeff1 = (v_lsr/np.cos(b))+ (theta_sun*np.sin(l))
    coeff2 = a1*theta_sun*np.sin(l)
    coeff3 = a3*theta_sun*np.sin(l)
             
    final = coeff1*term1 - coeff2*term2 - coeff3
    
    return final



def dist(R, l, b, r_sun=R_sun):
    term1 = r_sun*np.cos(l)/np.cos(b)
    
    term2_1 = R**2
    term2_2 = (r_sun*np.sin(l))**2
    term2 = np.sqrt(term2_1 - term2_2)/np.cos(b)
    
    if R<=r_sun:
        near = term1 - term2
        far = term1 + term2
        
        return [near,far]
    else:
        d = term1 + term2
        return [d]


#===========================================
possible_Rs = np.arange(1.7, 17, 0.1)*1000*u.pc

poly_out = np.zeros(len(possible_Rs))
def calculate_kd(l,b,v_lsr):
    l = l*u.deg
    b = b*u.deg
    v_lsr = v_lsr * u.km/u.s
    for i, posR in enumerate(possible_Rs):
        p = R_poly(posR,l,b,v_lsr)

        poly_out[i] = p.value

    this = np.argmin(np.abs(0-poly_out))

    R_gal = possible_Rs[this]

    d = dist(R_gal,l,b)

    return d


