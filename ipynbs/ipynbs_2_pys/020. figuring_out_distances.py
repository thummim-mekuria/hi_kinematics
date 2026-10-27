#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as math
import scipy.optimize as sci_op 

import astropy.units as u
import astropy.constants as c
from astropy.io import fits

import datacube as dc
import matplotlib

from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits


Orion_clob = dc.unpack_fits_cube('DHT27_Orion_raw.fits', 'DHT')


zeroth = Orion_clob.calculate_0th_moment()


im = zeroth
vmin, vmax = zlims(im)
plt.imshow(im, vmin=vmin, vmax=vmax, origin = 'lower')
plt.colorbar()
plt.show()


im = zeroth[0:50, 40:180]
vmin, vmax = zlims(im)
imx = Orion_clob.lons[40:180]
imy = Orion_clob.lats[0:50]

plt.imshow(im, vmin=vmin, vmax=vmax, origin = 'lower', \
           extent = [imx[0], imx[-1], imy[0], imy[-1]])
plt.colorbar()
plt.show() 


oriona = Orion_clob.cube[10:45,0:50,40:180]
oriona_lons = Orion_clob.lons[0:50]
oriona_lats = Orion_clob.vels[10:45]


plt.subplots(figsize = (40,10))
im = np.nansum(oriona*Orion_clob.dv, axis=1)
imx = oriona_lons
imy = oriona_lats
vmin, vmax = zlims(im)

plt.imshow(im.value,  origin = 'lower',\
          extent = [imx[0], imx[-1], imy[0], imy[-1]])
plt.colorbar()
plt.show() 


orionb = Orion_clob.cube[:,35:95, 125:200]


im = np.nansum(orionb*Orion_clob.dv, axis=0)
vmin, vmax = zlims(im)

plt.imshow(im.value,  origin = 'lower',\
          )
plt.colorbar()
plt.show() 


im = np.transpose(np.nansum(orionb*Orion_clob.dv, axis=2))
vmin, vmax = zlims(im)
imx = Orion_clob.vels
imy = Orion_clob.lats[35:95]

plt.imshow(im.value,  origin = 'lower',\
          extent = [imx[0], imx[-1], imy[0], imy[-1]])
plt.colorbar()
plt.show() 


#values from Brand & Blitz 1993
Theta_sun = 220*u.km/u.s
R_sun = 8.5*1000*u.pc
a1 = 1.00767
a2 = 0.0394
a3 = 0.00712


#values from Russeil et al 2017
Theta_sun = 240*u.km/u.s
R_sun = 8.34*1000*u.pc
a1 = 1.022
a2 = 0.0803
a3 = 0.


def R_poly (R, l, b, v_lsr, theta_sun, r_sun, a1, a2, a3):

    term1 = R/r_sun
    term2 = (R/r_sun)**a2
    
    coeff1 = (v_lsr/np.cos(b))+ (theta_sun*np.sin(l))
    coeff2 = a1*theta_sun*np.sin(l)
    coeff3 = a3*theta_sun*np.sin(l)
             
    final = coeff1*term1 - coeff2*term2 - coeff3
    
    return final


def dist(R, r_sun, l, b):
    term1 = r_sun*np.cos(l)/np.cos(b)
    
    term2_1 = R**2
    term2_2 = (r_sun*np.sin(l))**2
    term2 = np.sqrt(term2_1 - term2_2)/np.cos(b)
    
    final = term1 + term2
    
    return final


def dist_yarn(R, r_sun, l, b):
    term1 = r_sun*np.cos(b)*np.cos(l)
    
    term2_1 = term1**2
    term2_2 = (r_sun**2)
    term2_3 = R**2
    
    term2 = np.sqrt(term2_1 - term2_2 + term2_3)
    
    final = term1 + term2
    
    return final


oriona_vlsrs = (Orion_clob.vels[19:35]*Orion_clob.units['vels']).to(u.km/u.s)
oriona_ls = (oriona_lons)*u.deg.to(u.rad)
oriona_bs = (oriona_lats)*u.deg.to(u.rad)


orionaB = np.nanmedian(oriona_bs)
orionaL = np.nanmedian(oriona_ls)


possible_Rs = np.arange(1.7, 17, 0.1)*1000*u.pc


r_out = np.zeros((len(possible_Rs),len(oriona_vlsrs)))


for r, R in enumerate(possible_Rs):
    this_R = R
    row = np.zeros_like(oriona_vlsrs)
    for v, V in enumerate(oriona_vlsrs):
        row[v] = R_poly(R, orionaL, orionaB, \
                          V, Theta_sun, R_sun, a1, a2, a3)
    r_out[r] = row.value
        

r_out_med = np.nanmedian(r_out, axis =1)


plt.plot(possible_Rs,r_out_med)
plt.ylabel("polynomial, [km/s]")
plt.xlabel("possible Rs [pc]")

plt.hlines(0,possible_Rs[0].value, possible_Rs[-1].value, linestyle='--')

plt.show()


p_possible_Rs = possible_Rs[70:80]


dists = dist(p_possible_Rs, R_sun, orionaL, orionaB)
dists_y = dist_yarn(p_possible_Rs, R_sun, orionaL, orionaB)


((dists-dists_y)/dists_y)*100


plt.plot(dists, linestyle = '--')
plt.plot(dists_y, linestyle = '--')

plt.show()


864.43*np.sin(orionaB)*0.5


def select_bs_from_ds (d, H):
    #the equation for scale height can be given by H = d * np.sin(b)
    b = np.arcsin(H/d)
    return b


dists = dist(possible_Rs, R_sun, orionaL, orionaB)


bs_scale_height = (select_bs_from_ds(dists, 100*u.pc)).to(u.deg)


np.argmin(np.abs(bs_scale_height - 53*u.deg))


np.nanmin(bs_scale_height)


plt.hist(bs_scale_height, bins = 100)
plt.xlim(-10,10)

plt.show()


bs_scale_height[64:80]


dists[64:80]


possible_Rs[64:80]


orionaL*u.rad.to(u.deg)


r_out.shape


possible_Rs.shape


np.argmin(np.abs(0-r_out_med))


np.arctan(13/8500)


13/(np.sin(np.arctan(13/8500)))


r_out_med[70:80]


r_out_med[70:80]


possible_Rs[70:80]


p_possible_Rs = possible_Rs[70:80]


r_out_med[70:80]


p_possible_Rs


dist_yarn(p_possible_Rs, r_sun, orionaL, orionaB)


dists = dist(p_possible_Rs, r_sun, orionaL, orionaB)


dists


dist(possible_Rs[75], r_sun, orionaL, orionaB)


grid = np.meshgrid(oriona_bs, oriona_ls, oriona_vlsrs)


orionaB = np.nanmedian(oriona_bs)
orionaL = np.nanmedian(oriona_ls)


orionB


-0.33*360/(2*3.14)


len(grid[0][0][0])


len(oriona_vlsrs), len(oriona_bs), len(oriona_ls)


grid = np.array(grid)


212-360


np.sin(-148)


np.sin(212)


-148*u.deg.to(u.rad)


212*u.deg.to(u.rad)


np.sin(3.7000980142279785), np.sin(-2.5830872929516078)


#from table 1 of https://arxiv.org/pdf/0905.0723.pdf
rdv_ls = np.array((53.59, 29.89, 49.49, 18.89, 30.49, \
                          35.14, 34.24, 19.94, 23.44, 38.94, 23.44,\
                         30.79, 30.29, 53.14, 22.44, 24.49, 49.49, 19.39, \
                          34.74, 23.04, 18.69, 18.19, 25.64, 24.79,\
                         30.44, 23.24, 19.89, 22.04, 24.39, 30.54))*(u.deg.to(u.rad))

rdv_bs = np.array((0.04, -0.06, -0.41, -0.51, -0.36, -0.76, 0.14,\
                  -0.81, -0.21, -0.46, -0.21, -0.06, -0.21, 0.04,\
                   0.34,0.49, -0.26, -0.01, -0.66, -0.41, -0.06,\
                  -0.31, -0.11, 0.09, -0.26, -0.36, -0.56, 0.19,\
                   0.14, -0.06))*(u.deg.to(u.rad))

rdv_vs = np.array((23.74, 100.68, 56.90, 65.82, 12.26, 35.22, 57.75,\
                  42.87, 101.10, 41.59, 103.65, 94.73, 104.50, 22.04,\
                  84.52, 102.38, 50.94, 26.72, 46.69, 74.32, 45.42,\
                  50.09, 93.88, 110.45, 103.65, 77.30, 44.14, 50.94, \
                  113.43, 93.45))*(u.km/u.s)


r_out = np.zeros((len(rdv_bs), len(possible_Rs)))


r_out = np.zeros_like(rdv_bs)


row[v] = R_poly(R, orionaL, orionaB, \
                          V, theta, r_sun, a1, a2, a3)


rdv_ls


for c in range(len(rdv_bs)):
    its_b = rdv_bs[c]
    its_l = rdv_ls[c]
    its_v = rdv_vs[c]
    
    p = np.zeros(len(possible_Rs))
    for r, R in enumerate(possible_Rs):
        p[r] = R_poly(R, its_l, its_b, its_v, theta, r_sun, a1, a2, a3).value
    
    most_prob = np.argmin(np.abs(0.0 - p))
    r_out[c] = possible_Rs[most_prob].value
        

r_out = r_out*u.pc


dist_yarn(r_out, r_sun, rdv_ls, rdv_bs)


yarn_kpc = (dist_yarn(r_out, r_sun, rdv_ls, rdv_bs))/1000

BB_kpc = (dist(r_out, r_sun, rdv_ls, rdv_bs))/1000


yarn_kpc


yarn_kpc


BB_kpc


BB_kpc


(dist_MD(r_out, r_sun, rdv_ls, rdv_bs))/1000


(dist_MD(r_out, r_sun, rdv_ls, rdv_bs))/1000


def dist_MD(R, r_sun, l, b):
    term1 = r_sun*np.cos(l)
    
    term2_1 = R**2
    term2_2 = (r_sun*np.sin(l))**2
    term2 = np.sqrt(term2_1 - term2_2)
    
    final = term1 - term2
    
    return final


def dist(R, r_sun, l, b):
    term1 = r_sun*np.cos(l)/np.cos(b)
    
    term2_1 = R**2
    term2_2 = (r_sun*np.sin(l))**2
    term2 = np.sqrt(term2_1 - term2_2)/np.cos(b)
    
    final = term1 + term2
    
    return final


def dist_yarn(R, r_sun, l, b):
    term1 = r_sun*np.cos(b)*np.cos(l)
    
    term2_1 = term1**2
    term2_2 = (r_sun**2)
    term2_3 = R**2
    
    term2 = np.sqrt(term2_1 - term2_2 + term2_3)
    
    final = term1 - term2
    
    return final


def draw_bounds (x1, y1, length, lcolor = 'xkcd:grey', lwidth = 10):
    kwgs = {'color':lcolor, 'linewidth' : lwidth}
    
    #define 3 x & 5 y positions, all the lines go between coords 
    #that are made up by mixing & matching these x & ys
    
    spacing = length/20
    
    x1 = x1 
    x1r = x1 + spacing
    
    x2 = x1 + length
    x2l, x2r = x2 - spacing, x2 + spacing
    
    x3 = x2 + length
    x3l = x3 - spacing
    x3r = x3 + spacing
    
    
    y1 = y1
    y1a = y1 + spacing
    
    y2 = y1 + length
    y2a = y2 + spacing
    y2b = y2 - spacing
    
    y3 = y1 + 2*length
    y3a = y3 + spacing
    y3b = y3 - spacing
    
    y4 = y2 + 2*length
    y4a = y4 + spacing
    y4b = y4 - spacing
    
    y5 = y3 + 2*length
    y5a = y5 + spacing
    y5b = y5 - spacing
    

    plt.vlines(x1, y1a, y3b, **kwgs)
    plt.vlines(x2, y1a, y2b, **kwgs)
    plt.vlines(x2, y2a, y3b, **kwgs)
    plt.vlines(x3, y1a, y3b, **kwgs)
    
    plt.vlines(x1, y3a, y5b, **kwgs)
    plt.vlines(x2, y3a, y4b, **kwgs)
    plt.vlines(x2, y4a, y5b, **kwgs)
    plt.vlines(x3, y3a, y5b, **kwgs)
    
    plt.hlines(y1, x1r, x2l, **kwgs)
    plt.hlines(y1, x2r, x3l, **kwgs)
    plt.hlines(y3, x1r, x2l, **kwgs)
    plt.hlines(y3, x2r, x3l, **kwgs)
    plt.hlines(y5, x1r, x2l, **kwgs)
    plt.hlines(y5, x2r, x3l, **kwgs)
    
    plt.plot((x1r, x2l),(y3b, y2), **kwgs)
    plt.plot((x1r, x2l),(y3a, y4), **kwgs)
    plt.plot((x3l, x2r),(y3a, y4), **kwgs)
    plt.plot((x3l, x2r),(y3b, y2), **kwgs)

    
#     plt.vlines(l0, b0, b1, label = name, **kwgs)
#     plt.vlines(l1, b0, b1, **kwgs)
    
#     plt.hlines(b0, l0, l1, **kwgs)
#     plt.hlines(b1, l0, l1, **kwgs)


plt.subplots(figsize = (10,20))
draw_bounds(0,0,100)

plt.show()


d = 500*u.pc
R = 1*u.pc
T = 8000*u.K
n = 2e3/(u.cm**3)

omega = ((R/d)**2)
S_coeff = 2*c.k_B*T/c.c**2

tau_coeff_old = (4*(c.e.emu**6)*(n**2)*R*((2*np.pi)**0.5))/(3*c.m_p*c.c*c.k_B*T*(3*c.k_B*c.m_p*T)**0.5)


def f (nu, tau, S, om):
    exp = (tau*nu**-2)
    term1 = 1-(np.e**exp)
    term2 = S*om*(nu**2)
    
    return term1*term2
    

nus = (c.c/(np.linspace(1*u.um, 100*u.cm, 10000))).cgs


fs = (f(nus,tau_coeff, S_coeff, omega))
fs_Jy = fs.to(u.Jy)


fs


plt.plot(nus, fs_Jy)
# plt.ylim(-0.5,0.5)
plt.xlim(1e8,1e12)
plt.xscale('log')
# plt.yscale('log')

plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$F_\nu$ [Jy]')

plt.show()
