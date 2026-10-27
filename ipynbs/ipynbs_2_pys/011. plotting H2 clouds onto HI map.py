#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# # notes for check in 05.10
# * plotted Rice, MD, Solar Neighborhood clouds & CO emission on top of HI
# * notice that CO emission is concentrated to a smaller rage of b (-20,20) than the 21cm emission (-50,50)
# * both Rice & MD clouds are clustered around b = 0
#     * more restricted in b than the solar neighborhood clouds
#      * plotted the clouds over just the CO emission & it looks like the Rice/MD catalogs don't pick out CO emission > ~b=50 
#      --> this could be a feature of my plotting though, need to plot the rectangles (l,b + size) in order to really know
#      
# * plotted some empirical relations (just with Rice & MD, then adding in Solar Neighborhood)
#     - noticed that there were 'near' & 'far' values fiven for R, M, d, etc but on the whole the patterns are consistent whether what's being plotted is near or far quantities
# 
# * (from plots / histograms) can see that the Solar Neighborhood clouds are both more massive & have larger R than the catalog clouds 
#         -- the M-R fit that works for Rice/MD has the same slope, but needs to be moved vertically to fit the Solar Neighborhood clouds
#     
#  
# * haven't calculated virial parameters for MD clouds yet -- but the pattern with Rice & Zucker clouds is kind of similar
#     
#     
# ## moving forward notes from 05.10
# * read the methods of Rice & MD to see what they dd to identify the clouds
#     - [Rice+16 paper](https://iopscience.iop.org/article/10.3847/0004-637X/822/1/52/pdf)
#     - [MD+17 paper](https://arxiv.org/pdf/1610.05918.pdf)
# * look into the l/b positions in Rice & MD catalogs and see if they identify the solar neighborhood clouds // read the paper to see if they mention anything about off-plane clouds
# * look at the details of the LAB vs other HI full sky surveys
# 
#     
#     
# * are we probing all the way to the center of the galaxy
# * look at 13C0 surveys
#     
# #### moving forward notes from 05.05
# 
# * take the spatial locations of the H2 clouds -- and overplot onto HI survey
#     - use both Rice & MD
#     - overplot the solar neighborhood clouds i already know about
#     
# * see if MD comments on why they have more clouds than Rice
#     - look at the Lada & Dame paper to see what they say about both ^
# 
# 
# 
# * pin for later: looks like virial parameters are > 1 for a lot of the clouds
#     - check the literature to see if anyone has noticed this
#     - what are the implications of this, physically?
#     
#     
#     
#   
# #### keep in back of mind: what are bits of information we want to know about the HI clouds -- what do papers have to say?
# 


import datacube as dc
import clouds as cld


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c


from astropy.visualization import ZScaleInterval as zscale
zlims = zscale().get_limits

cmap0 = dc.cmap_grey0('viridis')


kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}


# # 1.sanity check: read in the CO & HI maps + plot 0th moments

# ### 1a. Read in LAB HI cube


#LAB fits -- has been superseded but this is the one that works so will use it for now
LAB_HI = dc.unpack_fits_cube('labh.fits','DHT')
HI_cube = LAB_HI[0]
HI_vels = LAB_HI[1]
HI_lats = LAB_HI[2]
HI_lons = LAB_HI[3]
HI_units = LAB_HI[4] 
HI_units['vels'] = u.m/u.s
#had to go in and change units of vels to m/s, so need to edit dc.unpack_fits_cube
#to make sure it reads the header for units instea of hard-coding 

HI_dv = ((HI_vels[1] - HI_vels[0])*HI_units['vels']).to(u.km/u.s)
HI_db = (HI_lats[1] - HI_lats[0])*HI_units['lats']
HI_dl = (HI_lons[1] - HI_lons[0])*HI_units['lons']

HI_0th = np.nansum(HI_cube*HI_dv.value, axis=0)


# ### 1b. Read in DHT CO/H2 cube


#DHT CO survey
DHT_CO = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
CO_cube = DHT_CO[0]
CO_vels = DHT_CO[1]
CO_lats = DHT_CO[2]
CO_lons = DHT_CO[3]
CO_units = DHT_CO[4]

#assume arrays have linear spacing
CO_dv = (CO_vels[1] - CO_vels[0])*CO_units['vels']
CO_db = (CO_lats[1] - CO_lats[0])*CO_units['lats']
CO_dl = (CO_lons[1] - CO_lons[0])*CO_units['lons']

CO_0th = np.nansum(CO_cube*CO_dv.value, axis=0)


# ### 1c. Plot CO 0th moment map


CO_0th_pos = np.abs(CO_0th)


(np.abs(CO_db*CO_dl)).to(u.arcmin**2)


### GALAXY 0TH MOMENT (in CO) ####
plt.figure(figsize=(15,3), dpi=200)

im = CO_0th_pos
imx = CO_lons
imy = CO_lats


vmin,vmax = zlims(im)   

plt.imshow(im,  vmin=vmin, vmax=vmax,\
          extent = [imx[0],imx[-1], imy[0],imy[-1]],  **kw_imshow)

plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{CO}$ $(K\times km/s)$')

plt.contour(imx,imy,im,levels=[5], linewidths=0.7 ,colors='xkcd:pale orange')
# plt.xlim(800,900)
# plt.ylim(-30,-5)

plt.title('Galaxy in CO: 0th moment')
plt.show()


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

plt.contour(imx,imy,im,levels=[5], linewidths = 0.7, colors='xkcd:pale orange')
# plt.xlim(800,900)
# plt.ylim(350,450)

plt.title('Galaxy in CO: 0th moment')
plt.show()


# ### 1d. Plot HI 0th moment map


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(15,3), dpi=200)

im = HI_0th
imx = HI_lons
imy = HI_lats


vmin,vmax = zlims(im) 


plt.imshow(im,  \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow,\
           norm = mpl.colors.LogNorm())


plt.ylabel('Galactic Latitude')
plt.xlabel('Galactic Longitude')
plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$')


plt.title('Galaxy in HI: 0th moment')
plt.show()


# # 2. Overplot Rice & MD clouds on HI

# ### 2a. Read in Rice+16 catalog & put in pandas dataframe


#the catalog is in the form of a table, 
#so have to select the [1] item instead of the [0] for images
Rice_header = fits.open('Rice+16_cloud_catalog.fits')[1].header
Rice_data_raw = fits.open('Rice+16_cloud_catalog.fits')[1].data

Rice_data_raw_list = Rice_data_raw.tolist()


#should probably package this into a function in datacube
num_cols = Rice_header['TFIELDS']
columns = [None]*num_cols
units = [None]*num_cols
for i in range(num_cols): 
    col_name = Rice_header['TTYPE'+str(i+1)] #do it this way because the header fields start at 1
    columns[i] = col_name
    
    #check if a unit is associated with a column
    if 'TUNIT'+str(i+1) in Rice_header.keys(): 
        unit = Rice_header['TUNIT'+str(i+1)]
        units[i] = unit
    else:
        units[i] = ''
        

#initialize the dataframe
Rice_df_raw = pd.DataFrame(columns=columns)
for r in Rice_data_raw_list:
    #assign position in the database based on 'cloud index'
    ind = r[0]-1  #'cloud index' starts at 1, but python starts counting at 0
    
    #then add each cloud as a row
    Rice_df_raw.loc[ind] = r


#drop Nans along the mass column
Rice_df = Rice_df_raw.dropna(subset=['mass'])


Rice_df


# ### 2b. Read in MD+17 catalog in put in pandas dataframe


MD_df = pd.read_csv('MD2017_catalog.dat', delim_whitespace=True)
#MD2017_catalog.dat is a file i downloaded from Vizier 
#+ edited by adding column headings according to the ReadMe from the tar unzip


MD_df


# ### 2c. Read in Zucker catalog and put in pandas frame


Zucker_df = pd.read_csv('local_molecular_clouds_zucker',sep=',',skiprows=4)


Zucker_df


# ### 2d. pull out lat/lon values for the Rice clouds


Rice_lons = Rice_df['l (deg)'].to_numpy()
Rice_lats = Rice_df['b (deg)'].to_numpy()

#^ these are in 0-360 form, so convert to -/+ 180
Rice_lons_180 = np.zeros_like(Rice_lons)
for i in range(len(Rice_lons)):
    l = Rice_lons[i]
    if l > 180:
        l -= 360
    Rice_lons_180[i] = l


# ### 2e. pull out lat/lon values for the MD clouds


MD_lons = MD_df['GLON'].to_numpy() #already in -/+180 blesssss
MD_lats = MD_df['GLAT'].to_numpy()


# ### 2f. pull out lat/lon values for the Zucker clouds


#because the l&b bounds are what is saved, will need to pull out the two values & average them
Zucker_lons = 0.5*(Zucker_df['l0'] + Zucker_df['l1'])
Zucker_lats = 0.5*(Zucker_df['b0'] + Zucker_df['b1'])


#then convert to -/+ 180
Zucker_lons_180 = np.zeros_like(Zucker_lons)
for i in range(len(Zucker_lons)):
    l = Zucker_lons[i]
    if l > 180:
        l -= 360
    Zucker_lons_180[i] = l


# ### 2f. then plot the clouds as points on top of HI map


Zucker_df


### GALAXY 0TH MOMENT (in HI) ####
plt.figure(figsize=(20,10), dpi=200) #h
# plt.figure(figsize=(10,15), dpi=100) #v


########################### 21cm emission ###########################


im = HI_0th #h
imx = HI_lons
imy = HI_lats

im_trans = np.transpose(HI_0th) #v


vmin,vmax = zlims(im) 
#


plt.imshow(im, vmin=vmin, vmax=vmax, \
            extent = [imx[0],imx[-1], imy[0],imy[-1]],**kw_imshow) #h

# # plt.imshow(im_trans, \
#             extent = [imy[0],imy[-1], imx[0],imx[-1]],**kw_imshow) #v


########################### 21cm emission ###########################


########################### CO emission ###########################
im2 = CO_0th
im2x = CO_lons
im2y = CO_lats


vmin2,vmax2 = zlims(im2)    


# plt.contour(im2x,im2y,im2,levels=[5], colors='xkcd:neon blue',linewidths=1)


# plt.imshow(im2, vmin=vmin2, vmax=vmax2,\
#             extent = [im2x[0],im2x[-1], im2y[0],im2y[-1]],**kw_imshow, alpha = 1)


########################### CO emission ###########################


########################### clouds ###########################


#------- MD clouds -------
plt.scatter(MD_lons, MD_lats, s=3, color='xkcd:powder blue') #h
# # plt.scatter(MD_lats,MD_lons, s=0.5, color='white', alpha = 0.5) #v

# #------- Rice clouds -------
plt.scatter(Rice_lons_180, Rice_lats, s=3, color='xkcd:pinkish red') #h
# # # plt.scatter(Rice_lats,Rice_lons_180, s=0.5, color='red', alpha = 0.5) #v


#------- Zucker clouds -------
plt.scatter(Zucker_lons_180, Zucker_lats, s=400, color='xkcd:bright lavender',\
            marker='+', linewidths=4) #h
plt.scatter(Zucker_lons_180[17], Zucker_lats[17], s=600, color='xkcd:bright pink',\
            marker='+', linewidths=4) #orion
plt.scatter(17, 0, s=600, color='xkcd:bright pink', marker='+', linewidths=4) #m16/eagle

# plt.scatter(Zucker_lats,Zucker_lons_180, s=0.5, color='red', alpha = 0.5) #v

########################### clouds ###########################

plt.ylabel('Galactic Latitude') #h
plt.xlabel('Galactic Longitude') #h
# plt.colorbar(orientation='vertical', label = r'$W_{HI}$ $(K\times km/s)$') #h

# plt.ylabel('Galactic Longitude') #v
# plt.xlabel('Galactic Latitude') #v
# plt.colorbar(orientation='horizontal', label = r'$W_{HI}$ $(K\times km/s)$') #v

# plt.ylim(-20,20) #h
# plt.xlim(-20,20) #v


# plt.title('Galaxy in HI: 0th moment')
plt.show()


# # 3. empirical relations

# ### 3a. Pull M & R from MD


#MD only has Mnear & Mfar, no plain M
#so use both
MD_M_n = (MD_df['Mnear'].to_numpy())*u.M_sun
MD_M_f = (MD_df['Mfar'].to_numpy())*u.M_sun

MD_R_n = (MD_df['Rnear'].to_numpy())*u.pc
MD_R_f = (MD_df['Rfar'].to_numpy())*u.pc


# ### 3b. Pull out M & R from Rice


#Rice has near_mass, far_mass, mass -- for the sake of consistency
#just use near_mass & far_mass
R_M_n = (Rice_df['near_mass'].to_numpy())*u.M_sun
R_M_f = (Rice_df['far_mass'].to_numpy())*u.M_sun 

R_R_n = (Rice_df['near_size'].to_numpy())*u.pc
R_R_f = (Rice_df['far_size'].to_numpy())*u.pc


# ### 3c. Calculate M & R for Zucker clouds


#create arrays to hold quantities to be calculated
Z_M = np.zeros_like(Zucker_df['l0'])
Z_A = np.zeros_like(Z_M)
Z_R = np.zeros_like(Z_M)
Z_Sigma = np.zeros_like(Z_M)
for i in range(26):    
    
    #add in a check that scans that lat & lon values of clouds & first confirms that 
    #they're in the galaxy cube, before calling cld.extract
    #something along the lines of 
    #if l0 in range(np.min(glat),np.max(glat))
    
    
    #extract the cloud
    extracted = cld.extract_cloud(i, Zucker_df, DHT_CO)
    
    cube = extracted[0]
    
    #0th moment map
    W = np.nansum(cube*CO_dv, axis=0)
    
    #velocity spectrum
    spec = dc.spectrum1D(cube)
    
    #velocity dispersion
    #convert units from pixel/index to velocity ones(based on gvels<->index) + add units

    sigma = CO_dv * dc.dispersion_1sigma(spec)
    sigma = sigma.to(u.m/u.s)
    Z_Sigma[i] = sigma.value
    sigma_unit = sigma.unit
    
    #c.G is in m^3/kg * s^2 so convert units to work
    #mass & area
    mass, area = cld.calc_mass_and_area(i, Zucker_df,DHT_CO,1.5)
    mass = mass.to(u.kg)
    Z_M[i] = mass.value
    mass_unit = mass.unit

    
    area = area.to(u.m**2)
    Z_A[i] = area.value
    area_unit = area.unit

    radius = np.sqrt(abs(area/np.pi))
    Z_R[i] = radius.value

# #once all the quantities have been calculated, add to database
# #eventually = write a function so units are saved in column names 
# CLOUDS.insert(len(CLOUDS.columns),'mass',masses)
# m_unit = [str(mass_unit)]*len(masses)
# CLOUDS.insert(len(CLOUDS.columns),'m_unit',m_unit)

# CLOUDS.insert(len(CLOUDS.columns),'area',areas)
# a_unit = [str(area_unit)]*len(areas)
# CLOUDS.insert(len(CLOUDS.columns),'a_unit',a_unit)

# CLOUDS.insert(len(CLOUDS.columns),'velocity_dispersion',sigmas)
# sig_unit = [str(sigma_unit)]*len(sigmas)
# CLOUDS.insert(len(CLOUDS.columns),'sig_unit',sig_unit)


# ### 3c. log-ify all of them


logMD_M_n = np.log10(MD_M_n.value)
logMD_M_f = np.log10(MD_M_f.value)
logMD_R_n = np.log10(MD_R_n.value)
logMD_R_f = np.log10(MD_R_f.value)

logR_M_n = np.log10(R_M_n.value)
logR_M_f = np.log10(R_M_f.value)
logR_R_n = np.log10(R_R_n.value)
logR_R_f = np.log10(R_R_f.value)

logZ_M = np.log10(Z_M)
logZ_R = np.log10(Z_R)


# ### 3d. calculate fit


def linear_fit (x, m, b):
    y = m*x + b
    
    return y


#set up the radius array to be nice input for the fit
fit_x = np.arange(-2,3,1e-1)
# fit_x = np.arange(-2,20,1e-1) #fit with Solar Neighborhood clouds


MRfit_y = linear_fit(fit_x, 2.3, 1.5)
# MRfit_y = linear_fit(fit_x, 2.3, -5.5) #fit with Solar Neighborhood clouds


# ### 3e. Plot M vs R with fit


## M vs R

#----MD----
plt.scatter(logMD_R_n, logMD_M_n,s=0.5,label='MD near', alpha=0.7, color = 'xkcd:bright blue')
plt.scatter(logMD_R_f, logMD_M_f,s=0.5,label='MD far', alpha=0.7, color = 'xkcd:powder blue')


#----Rice----
plt.scatter(logR_R_n, logR_M_n,s=0.5,label='Rice near', alpha=0.7, color = 'xkcd:crimson')
plt.scatter(logR_R_f, logR_M_f,s=0.5,label='Rice far', alpha=0.7, color = 'xkcd:pinkish red')


# # #----Zucker----
# plt.scatter(logZ_R, logZ_M,s=0.5,label='Solar Neighborhood', \
#             alpha=0.7, color = 'xkcd:pastel purple')


#----fit----
plt.plot(fit_x, MRfit_y, linestyle='--', color = 'black', alpha = 0.4, label=r'log M=2.3(log R)+1.5')

plt.ylabel(r'log Masses [$M_\odot$]')
plt.xlabel('log R [pc]')

plt.title(r'Mass-Radius relationship of $H_2$ clouds from Rice+16 & MD+17')
plt.legend()
plt.show()


# # 4. Histograms

# ### 4a. Masses


bin_num = 100

#----MD----
plt.hist(MD_M_n, bins=bin_num, histtype='step', label='MD near', color='xkcd:bright blue')
plt.hist(MD_M_f, bins=bin_num, histtype='step', label='MD far', color='xkcd:powder blue')

# ----Rice----
plt.hist(R_M_n, bins=bin_num, histtype='step', label='Rice near', color='xkcd:crimson')
plt.hist(R_M_f, bins=bin_num, histtype='step', label='Rice far', color='xkcd:pinkish red')

#----Zucker----
plt.hist(Z_M, bins=bin_num, histtype='step', label='Solar Neighborhood', color='xkcd:pastel purple')


plt.xscale('log')
plt.yscale('log')


plt.xlabel(r'Masses [$M_\odot$]')
plt.ylabel('N[Clouds]')

plt.legend()

plt.show()


# ### 4b. Sizes


bin_num = 100

#----MD----
plt.hist(MD_R_n, bins=bin_num, histtype='step', label='MD near', color='xkcd:bright blue')
plt.hist(MD_R_f, bins=bin_num, histtype='step', label='MD far', color='xkcd:powder blue')

#----Rice----
plt.hist(R_R_n, bins=bin_num, histtype='step', label='Rice near', color='xkcd:crimson')
plt.hist(R_R_f, bins=bin_num, histtype='step', label='Rice far', color='xkcd:pinkish red')

#----Zucker----
plt.hist(Z_R, bins=bin_num, histtype='step', label='Solar Neighborhood', color='xkcd:pastel purple')


plt.xscale('log')
plt.yscale('log')


plt.xlabel(r'Sizes [pc]')
plt.ylabel('N[Clouds]')

plt.legend()

plt.show()


# ### 4c. Virial Parameters


#pull out from Rice Catalog
R_vir = Rice_df['virial_alpha']


#calculate for MD & Zucker clouds
Z_vir = np.zeros_like(Z_M)
for i in range(len(Z_vir)):
     Z_vir[i] = (5*Z_R[i]*((Z_Sigma[i])**2))/(c.G.value*Z_M[i])
        
        
# CLOUDS.insert(len(CLOUDS.columns),'virial_parameter',vir_pars)


bin_num = 100


#----Rice----
plt.hist(R_vir, bins=bin_num, histtype='step', label='Rice', color='xkcd:crimson')

#----Zucker----
plt.hist(Z_vir, bins=bin_num-80, histtype='step', label='Solar Neighborhood', color='xkcd:pastel purple')


plt.xlim(0,20)
plt.xlabel(r'Virial $\alpha$s')

plt.yscale('log')

plt.legend()
plt.show()


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
