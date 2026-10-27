#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:85% !important; }</style>"))


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


# %%time

#all-sky map from HI4PI
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI_cube, HI_vels, HI_lats, HI_lons, HI_units = HI4PI

HI_dv = ((HI_vels[1] - HI_vels[0])*HI_units['vels']).to(u.km/u.s)
HI_db = (HI_lats[1] - HI_lats[0])*HI_units['lats']
HI_dl = (HI_lons[1] - HI_lons[0])*HI_units['lons']


# %%time

#DHT CO survey
DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
CO_cube, CO_vels, CO_lats, CO_lons, CO_units = DHT_CO

#assume arrays have linear spacing
CO_dv = ((CO_vels[1] - CO_vels[0])*CO_units['vels']).to(u.km/u.s)
CO_db = (CO_lats[1] - CO_lats[0])*CO_units['lats']
CO_dl = (CO_lons[1] - CO_lons[0])*CO_units['lons']

CO_0th = np.nansum(CO_cube*CO_dv.value, axis=0)


Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=5)


#these had problems with multiple peaks but got fixed with extract_cloud v3
for i in [1, 3, 4, 5, 6, 11, 13, 14, 19, 22]:
    cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, CO_padding=15,\
                      HI_padding=40,  v_padding = 10, \
                      show_0th = False, show_spec = True, sigma = 3)


# these had problems with multiple peaks
for i in [1, 22, 3, 4, 5, 19, 9, 11, 13, 14]:
    cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, CO_padding=15,\
                      HI_padding=40,  v_padding = 10, \
                      show_0th = False, show_spec = True, sigma = 3)


# cld.get_parent_HI(1, Zucker_df, DHT_CO, HI4PI, CO_padding=15, HI_padding=30, show_0th = False, show_spec = True)


for i in range(len(Zucker_df)):
    cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, CO_padding=15, HI_padding=40,  v_padding=10, show_0th = False, show_spec = True)


# something is going on with these clouds (shows up in both the 0th moment maps & velocity spectra)
#     - Rosette
#     - Orion Lam
#     - Mon OB1
#     - Maddalena
#     - Lacerta
#     - Gem OB1
#     - Cepheus (both near & far)
#     - Cam
#     - CMa OB1
#     
# * the trend here seems to be that all the faulty clouds are from the same surveys as each other (12, 18, 19, 24, 31)
# * and all the ones that work well are from the same surveys as each other (1, 3, 4, 6, 8, 16, 21, 27, 30, 37)
# 
# maybe the solution is to use a different sigma cut-off to define the signal, based on which category of surveys the data came from
#     --- okay weird, v_padding = 0 throws an error, but v_padding = 10 seems to work ?


#doesn't apply anymore -- was used to figue out if the systematic errors 
#(all now resolved by extract_cloud v3) had to do with which surveys things came from
data_qual = ['ok', 'ok', 'ok', 'ok', 'not', 'not', 'ok', 'ok',\
             'ok', 'ok', 'ok','not', 'ok', 'ok','ok', 'ok', 'ok', \
             'ok', 'ok', 'not', 'ok', 'ok','ok', 'ok', 'ok']

#----
Zucker_df['data_qual']= data_qual
#----
by_data_qual_Zucker_df = Zucker_df.sort_values('data_qual')
#----
l_av = 0.5*(by_data_qual_Zucker_df['l0']+by_data_qual_Zucker_df['l1'])
b_av = 0.5*(by_data_qual_Zucker_df['b0']+by_data_qual_Zucker_df['b1'])
#----
l_av_not_ok = l_av[:4]
b_av_not_ok = b_av[:4]
#----
l_av_ok = l_av[4:]
b_av_ok = b_av[4:]
#----
def draw_bounds (l0, l1, b0, b1, lcolor, lwidth, name):
    kwgs = {'color':lcolor, 'linewidth' : lwidth}
    plt.vlines(l0, b0, b1, label = name, **kwgs)
    plt.vlines(l1, b0, b1, **kwgs)
    
    plt.hlines(b0, l0, l1, **kwgs)
    plt.hlines(b1, l0, l1, **kwgs)
    
#----    
    
    
# fig, axs = plt.subplots(1,1, figsize = (15,3))

# plt.scatter(l_av_ok,b_av_ok, marker = '+', color = 'green', label = 'ok')
# # plt.scatter(l_av_not_ok,b_av_not_ok, marker = 'x', color = 'red', label = 'not ok')

# plt.xlabel(r'Galactic Longitude ($\degree$)')
# plt.ylabel(r'Galactic Latitude ($\degree$)')


# plt.vlines(0,-45,45, linestyle = ':', color = 'grey', label = 'l-b bounds')
# plt.vlines(360,-45,45, linestyle = ':', color = 'grey')

# plt.hlines(-30,-5,365, linestyle = ':', color = 'grey')
# plt.hlines(30,-5,365, linestyle = ':', color = 'grey')

# # plt.hlines(10,-5,365, linestyle = '--', color = 'grey', label = 'Galactic Plane')
# # plt.hlines(-10,-5,365, linestyle = '--', color = 'grey')


# draw_bounds(15,0,-5,20,'xkcd:grey', 0.5,'#01: Superbeam')
# draw_bounds(360,210,-5,20,'xkcd:grey', 0.5,'#01: Superbeam, pt 2')
# draw_bounds(4,0,-16,-24,'xkcd:orange', 0.5,'#03: R CrA')
# draw_bounds(360,358,-16,-24,'xkcd:orange', 0.5,'#03: R CrA, pt 2')
# draw_bounds(36,4,18,4,'xkcd:pink', 0.5,'#04: Aquila Rift')
# draw_bounds(40,29,-20,-12,'xkcd:lavender', 0.5,r'#6: Aquila South')
# draw_bounds(75,15,15,-10,'xkcd:spring green', 0.5,r'#8: First Quadrant')
# # draw_bounds(107,95,-20,-9,'xkcd:lavender', 0.5,r'#12: Lacerta')
# draw_bounds(130,118,18,34,'xkcd:orange', 0.5,'#16: Polaris Flare')
# # draw_bounds(210,70,-10,20,'xkcd:cyan', 0.5,'#18: Second Quadrant, 1/4 $\degree$ sampled')
# # draw_bounds(156,140,-10,20,'xkcd:red', 0.5,'#19: Camelopardalis')
# draw_bounds(180,155,-5,-30,'xkcd:pink', 0.5,'#21: Taurus')
# # draw_bounds(201,191,-18,-7,'xkcd:red', 0.5,r'#24: $\lambda$-Orionis')
# draw_bounds(200,226,-2,-22,'xkcd:goldenrod', 0.5,'#27: Orion')
# draw_bounds(294,306,-12,-20,'xkcd:red', 0.5,'#30: Chamaeleon')
# # draw_bounds(190,300,2,-5,'xkcd:red', 0.5,'#31: Third Quadrant')
# draw_bounds(12,0,6,26,'xkcd:cyan', 0.5,'#37: Ophiuchus')
# draw_bounds(360,344,6,28,'xkcd:cyan', 0.5,'#37: Ophiuchus, pt 2')


# plt.xlim(365,-5)
# plt.ylim(-35,35)

# # plt.legend()

# plt.show()

#----    

fig, axs = plt.subplots(1,1, figsize = (15,3))
plt.scatter(l_av_ok,b_av_ok, marker = '+', color = 'green', label = 'ok')
plt.scatter(l_av_not_ok,b_av_not_ok, marker = 'x', color = 'red', label = 'not ok')

plt.xlabel(r'Galactic Longitude ($\degree$)')
plt.ylabel(r'Galactic Latitude ($\degree$)')


plt.vlines(0,-45,45, linestyle = ':', color = 'grey', label = 'l-b bounds')
plt.vlines(360,-45,45, linestyle = ':', color = 'grey')

plt.hlines(-30,-5,365, linestyle = ':', color = 'grey')
plt.hlines(30,-5,365, linestyle = ':', color = 'grey')

# plt.hlines(10,-5,365, linestyle = '--', color = 'grey', label = 'Galactic Plane')
# plt.hlines(-10,-5,365, linestyle = '--', color = 'grey')


draw_bounds(15,0,-5,20,'xkcd:green', 0.5,'#01: Superbeam')
draw_bounds(360,210,-5,20,'xkcd:green', 0.5,'#01: Superbeam, pt 2')
draw_bounds(4,0,-16,-24,'xkcd:green', 0.5,'#03: R CrA')
draw_bounds(360,358,-16,-24,'xkcd:green', 0.5,'#03: R CrA, pt 2')
draw_bounds(36,4,18,4,'xkcd:green', 0.5,'#04: Aquila Rift')
draw_bounds(40,29,-20,-12,'xkcd:green', 0.5,r'#6: Aquila South')
draw_bounds(75,15,15,-10,'xkcd:green', 0.5,r'#8: First Quadrant')
draw_bounds(130,118,18,34,'xkcd:green', 0.5,'#16: Polaris Flare')
draw_bounds(180,155,-5,-30,'xkcd:green', 0.5,'#21: Taurus')
draw_bounds(200,226,-2,-22,'xkcd:green', 0.5,'#27: Orion')
draw_bounds(294,306,-12,-20,'xkcd:green', 0.5,'#30: Chamaeleon')
draw_bounds(12,0,6,26,'xkcd:green', 0.5,'#37: Ophiuchus')
draw_bounds(360,344,6,28,'xkcd:green', 0.5,'#37: Ophiuchus, pt 2')
draw_bounds(190,300,2,-5,'xkcd:green', 0.5,'#31: Third Quadrant')
draw_bounds(156,140,-10,20,'xkcd:green', 0.5,'#19: Camelopardalis')


draw_bounds(107,95,-20,-9,'xkcd:red', 0.5,r'#12: Lacerta')
draw_bounds(210,70,-10,20,'xkcd:red', 0.5,'#18: Second Quadrant, 1/4 $\degree$ sampled')
draw_bounds(201,191,-18,-7,'xkcd:red', 0.5,r'#24: $\lambda$-Orionis')


plt.xlim(365,-5)
plt.ylim(-35,35)

plt.xticks(np.arange(-10, 370, 20))
plt.yticks(np.arange(-40, 40, 10))


# plt.legend()

plt.show()


# ## alright so:
# color-coded based on whether the 3-sigma sigal threshold worked or not --- next step => look at Rice / Miville-Deschenes papers to see if they talk about needing to use different signal thresholds based on the inidiviual survey that an on-sky area came from (Rice talked about quadrants a lot?? was this a part of it? or was the focus on more real-world physical differences between the different parts imaged?)
# 
# -- not sure if there is any physical meaning to this -- first check the data side of things, but if that doesn't do anything, also possible that the thing to look at is a wider velocity range for HI etc...
#     - from looking at the data (this table: https://lweb.cfa.harvard.edu/rtdc/CO/Table1.jpg & table 1 in Rice) ... not seeing anything systematic that demarcates the green surveys from the red
#     
#     
#     
# - another option: check if the HI4PI data has anything systematic going on ... ? => look at the pairs of 0th moments to see if all the CO data is awry in similar ways as each other for the clouds
# 
# 
# -- looked again, and from the velocity spectra, something that is showing up is that all the ones that look correct are centered on ~0, and all the ones that look wonky are at high velocities
#     - check padding to see what's going on
#     - gonna plot the entire spectra of these guys to see if they only have signal in spurts, rather than a plateau


for i in [1, 22, 3, 4, 5, 19, 9, 11, 13, 14]:
    cld.get_parent_HI(i, Zucker_df, DHT_CO, HI4PI, CO_padding=15,\
                      HI_padding=40,  v_padding = 10, \
                      show_0th = True, show_spec = False)


# summary of what i've done on (07.07) --- 
# * wrote a function that takes in the cloud number (in the Zucker database), the CO cube, the HI cube, the CO padding, the HI padding, then makes a pair plots the CO emission from a molecular cloud + the HI region around the cloud (defined by the cloud location + HI padding) with the molecular cloud contour overplotted
# * wrote a function that takes in the cloud cube & the HI cube in the direction of the cloud (output by the first function) and plots the velocity spectra of the CO and HI emissions 
# 
# * tested both functions out on the Zucker catalog + confirmed that i get two plots + two spectra for each (except Ursa Major, something is going on with the CO data there)
# 
# * pattern i've noticed so far, the molecular clouds are not at the concentrated regions of HI, but just next to them ?? increase HI padding to see if this is true
#         - for Perseus, California, Taurus far the moelcular clouds seem to be @ around W = 400
#         
# 

# note:
# * Pegasus is out of the bounds of the DHT survey (|b| > 30) so i took out of my version of the Zucker catalog
# 

# next steps:
# * most of the clouds only have one HI gaussian along the line of sight, but some have multiple so figure out a robust way to do velocity cuts that only considers the HI gaussian overlayed on the CO gaussian
# * right now, the vels arrays are the entire length of the CO / HI vels arrays -> figure out ways to save the v0_noise & v1_noise values to use for plotting spectra (look at extract_cloud)
# * look at population level characteristics of the parent HI clouds
# * start thinking about how to define ^ these clouds (ie if i didnt know the CO emission was there, what could i have done to extract them) then use that definition in parts of the sky that don't have CO emission
# * read Kalberla paper


CO_0th_pos = np.abs(CO_0th)

### GALAXY 0TH MOMENT (in CO) ####
plt.figure(figsize=(15,3), dpi=200)

im = CO_0th_pos
imx = CO_lons
imy = CO_lats


vmin,vmax = zlims(im)   

plt.imshow(im, vmin=vmin, vmax=vmax, \
          extent = [imx[0],imx[-1], imy[0],imy[-1]],  **kw_imshow)

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{CO}$ $(K\times km/s)$')

plt.contour(imx,imy,im,levels=[5], linewidths=0.7 ,colors='xkcd:pale orange')
# plt.xlim(800,900)
# plt.ylim(-30,-5)

plt.title('Galaxy in CO: 0th moment')
plt.show()


# (HI_cube.size * HI_cube.itemsize)/1e9 #the size of HI_cube in gigs = 35.3


# HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0) 
#^^ this takes too long to run, so instead of making a 0th moment map of the whole
#thing, i will just cut out the regions around the solar neighborhood clouds


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


#3D maps based on HI4PI
MerPhan_3D_1 = fits.open('HI_dens_mean_BEG03.fits')[0].data
MerPhan_3D_2 = fits.open('HI_dens_mean_SBM15.fits')[0].data

vmin,vmax = zlims(MerPhan_3D_1)


MerPhan_xs = np.linspace(-24,24,len(MerPhan_3D_1))
MerPhan_ys = np.linspace(-24,24,len(MerPhan_3D_1[0]))
MerPhan_zs = np.linspace(-2,2,len(MerPhan_3D_1[0][0]))


MerPhan_3D_1_xy = np.nanmedian(MerPhan_3D_1, axis=2)
MerPhan_3D_1_yz = np.nanmedian(MerPhan_3D_1, axis=1)
MerPhan_3D_1_xz = np.nanmedian(MerPhan_3D_1, axis=0)


MerPhan_3D_2_xy = np.nanmedian(MerPhan_3D_2, axis=2)
MerPhan_3D_2_yz = np.nanmedian(MerPhan_3D_2, axis=1)
MerPhan_3D_2_xz = np.nanmedian(MerPhan_3D_2, axis=0)


fig, axs = plt.subplots(2,3, figsize=(10,15), dpi=100)

kws = {'origin':'lower', 'vmin':vmin, 'vmax':vmax}

axs[0,0].imshow(MerPhan_3D_1_xy, **kws,\
                extent=[MerPhan_xs[0],MerPhan_xs[-1], \
                        MerPhan_ys[0],MerPhan_ys[-1]])
axs[0,0].set_title('Face-On, BEG03')
axs[0,0].set_xlabel('X [kpc]')
axs[0,0].set_ylabel('Y [kpc]')


axs[0,1].imshow(MerPhan_3D_1_yz, **kws,\
                extent=[MerPhan_zs[0],MerPhan_zs[-1], \
                        MerPhan_ys[0],MerPhan_ys[-1]])
axs[0,1].set_title('Edge-On, BEG03')
axs[0,1].set_xlabel('Z [kpc]')
axs[0,1].set_ylabel('Y [kpc]')

axs[0,2].imshow(MerPhan_3D_1_xz, **kws,\
                extent=[MerPhan_zs[0],MerPhan_zs[-1], \
                        MerPhan_xs[0],MerPhan_xs[-1]])
axs[0,2].set_title('Edge-On, BEG03')
axs[0,2].set_xlabel('Z [kpc]')
axs[0,2].set_ylabel('X [kpc]')


axs[1,0].imshow(MerPhan_3D_2_xy, **kws,\
                extent=[MerPhan_xs[0],MerPhan_xs[-1], \
                        MerPhan_ys[0],MerPhan_ys[-1]])
axs[1,0].set_title('Face-On, SBM15')
axs[1,0].set_xlabel('X [kpc]')
axs[1,0].set_ylabel('Y [kpc]')

axs[1,1].imshow(MerPhan_3D_2_yz, **kws,\
                extent=[MerPhan_zs[0],MerPhan_zs[-1], \
                        MerPhan_ys[0],MerPhan_ys[-1]])

axs[1,1].set_title('Edge-On, SBM15')
axs[1,1].set_xlabel('Z [kpc]')
axs[1,1].set_ylabel('Y [kpc]')

axs[1,2].imshow(MerPhan_3D_2_xz, **kws,\
                extent=[MerPhan_zs[0],MerPhan_zs[-1], \
                        MerPhan_xs[0],MerPhan_xs[-1]])
axs[1,2].set_title('Edge-On, SBM15')
axs[1,2].set_xlabel('Z [kpc]')
axs[1,2].set_ylabel('X [kpc]')


plt.show()


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
