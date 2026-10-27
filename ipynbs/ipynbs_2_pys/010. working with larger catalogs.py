#!/usr/bin/env python
# coding: utf-8


from IPython.display import display, HTML
display(HTML("<style>.container { width:95% !important; }</style>"))


# # notes for check-in 5.05
# * started with [Rice+16 catalog](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/EWY90X)
#     - [Rice+16 paper](https://iopscience.iop.org/article/10.3847/0004-637X/822/1/52/pdf)
#     - plotted some empirical relations as a sanity check
#     - plotted some histograms
#     
# 
# * still working on figuring out data format for [MD+17 catalog](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/834/57#/article)
#     - [MD+17 paper](https://arxiv.org/pdf/1610.05918.pdf)
#     
# * also started to look at HI all sky survey for comparison
#     - [HI 4_Pi survey](https://lambda.gsfc.nasa.gov/product/foreground/fg_hi4pi_get.html)
#     
# * noticed that the analyses are population level (eg here is where the clouds are in the galaxy as a function of distance from galactic center), but we want to do more seeing stuff *within* each cloud? 
# * catalog from MD is ~8000 & Rice is ~1000, which one would be more useful for our purposes? (the total mass identified by MD is only ab ~10 Rice's, so their cloud finding algorithm seems more granular, rather than finding clouds Rice ddin't)
# * question: next step from here could be: looking at surface density v r of the cloud? (just bulk quantities are recorded though, so what's the best way to get more resolution/detail)
#     
# ## moving forward notes from 05.05
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
#     
#     
#     
# #### moving forward notes from 4.19:
# - look at catalogs from: https://arxiv.org/pdf/2006.08632.pdf Rice & MD's catalogs and whats going on there 
# - possible: see if there are any surveys that are more complete, but the majority of clouds in the galaxy are near the plane so Dame's does a pretty good job
#     - 13CO survey: Heyer, FCRAO
# - what can 13CO tell us? (in comparison to 12CO)
# - take a look @ https://www.nature.com/articles/s41550-023-01901-5 -- how do molecular clouds form? uses a CII
#     
# 
# #### keep in back of mind: what are bits of information we want to know about the HI clouds -- what do papers have to say?
# 


import datacube as dc
import clouds as cld


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits

from astropy import units as u
import astropy.constants as c


#open the fits file
HI_survey = fits.open('NHI_HPX.fits')


# ### HI survey
# 
# #### CHECKPOINT 1: documentation timestamp: near the end of coworking (~5pm) on 05/08
# so it looks like the way the data was stored is in 2D
# where all the pixels are along one axis
# and the second axis is information ab each pixel
#     - six pieces of information:
#     i) HPXINDEX => which pixel
#    ii) RA2000 => RA coordinates
#   iii) DEC2000 => Dec coordinates
#    iv) GLON => galactic lon coordinate
#     v) GLAT => galactic lat coordinate
#    vi) NHI => column density in 1/cm^2
#    
#    
#    so the pieces of information i need in order to make an all-sky map* of this are:
#    glon, glat, nhi
#    [3]   [4]   [5]
#      
# * the main goal being to overplot the H2 clouds from Rice/MD's catalogs onto an HI map of the galaxy in order to identify the parent atomic clouds 

# #### CHECKPOINT 2: documentation timestamp: 9:14pm
# 
# okay so what i have done so far is open up the fits bin table that HI_survey is stored as
# realized that thed data is stored in the ^^ above format
# 
# ##### ~~ understanding what's going on ~~
# * also realized that the way the pixels are ordered is by galactic latitude
# * also realized that i can access certain 'columns' in the table by doing .[column name]
#     - so did that to pull out glons glats (instead of by index like i had maybe started to do by writing [3][4][5] earlier)
# * i iterated through the glats array, and wrote a cell that populates two arrays: nums & inds
#     - nums saves the number of latitude values that are the same as each other (i.e. how many of the pixels correspond to a single row in latitude) 
#     - and inds saves the start and end indices of these rows in glats
#         - the start is the index of the first pixel that has that latitude value, the end is the index of the first pixel that *doesn't* have that latitude value
#             - did it this way instead of end equaling the index of the last pixel that has that value, bc when you index in python (ie list[start:end]) the subset that you get is the elements of list up to end-1
#             
#             
# ##### ~~ creating a 2D map/ image/ visualization  ~~
#    * this is when i will need to output the .NHI column and use those values to populate the image            
#    * so the length of nums tells me how many distinct values of latitude there are => 
#        - should be equal to the how big the data is along the glat axis
#    * the np.max of nums (ie the most number of pixels that have the same latitude value) tells me how many distinct longitude values there are
#        - ie this is the longest latitude row
#    * using inds, i pulled out the range in glons that corresponds to this longest_lat_row
#        - the spacing between the longitude values in this row is the smallest spacing, so i will use them to set dlon
#            - know this because the shortest latitude rows (at the beginning & end) have big jumps in longitude
#            - and as the latitude arrays get longer, they sample the longitude axis more frequently/better resolution
#    * dlon is consistent to ~9 decimal points (and the dlons that are not the same as the others round to that same value)
#        - so this will be what i use to construct the longitude axis 
#            - ig, its not necessary in the constructing step, but will be useful when figuring out where to assign the pixel values from the spars-er latitude rows
#        - **next step:figure out if the jump in latitudes is consistent
#    * **the goal from here is to then construct a lat/lon grid & populate it with the value of the HI column density at each spot
#             
#             

# #### CHECKPOINT 3: documentation timestamp: 9:58pm
# 
# * i plotted nums and can see that: there is a good number of rows that are all the longest length (haha)
#     used:
#     xscatter = np.arange(0,len(nums),1)
#     plt.scatter(xscatter, nums, s=1e-2)
#     to plot how nums (ie lenght of latitude rows) changes over latitude
#     
#     - but!! their longitude values arent neatly lined up ... so whatever lon/lat grid i make will be kinda patchy
#         - can see this by looking at the first longitude value + dlon of each row thats the same length as the max length
#         
# 
#     
# * also: i had assumed that W x H (where W = length of longest latitude row (okay alliteration) & H = number of latitudes) would == npix (ie the length of the fits bintable)
#     - but it does nottttt, npix > W x H, so gotta figure out what's going on there
#     
#     
# * hmm this feels like its way too convoluted to just read out some data from a fits table -- other people have used this dataset so gonna pivot & try to find their papers to see what they did / whether they describe how they extracted the data
#     - another option is: topcat -- figure out how to actually install the downloaded .dmg (or whatever format)
# 


# pull out the fits bintable
rawHI = HI_survey[1].data


rawHI


# pull out the columns i care about
pxs = rawHI.HPXINDEX
glons = rawHI.GLON
glats = rawHI.GLAT
nhis = rawHI.NHI


# go through glats & figure out how long each 'row' is, + save the start/end indices of each row
nums = []
inds = []
it = glats[0]
num = 0
start = 0
for i in range(1,len(glats)):
    if glats[i] == it:
        num += 1
    else:
        #save the most recent iteration
        end = i
        inds.append((start,end))
        nums.append(num)
        
        #reset things
        num = 0
        it = glats[i]
        start = i

#save the final things
inds.append((start,len(glats)))
nums.append(num)


# doing the stuff covered in checkpoint 2 --- before realizing that theres multiple longest latitude rows + that the longitudes arent lined up

#find the longest row, pull out its start and end
longest_glat = np.argmax(nums)
s,e = inds[longest_glat]

#then save the longitude values along this row
longest_glons = glons[s:e]

#calculate the dlon between each pixel along glon 
diffs = np.zeros_like(longest_glons)
for i in range(len(longest_glons)-1):
    diffs[i] = longest_glons[i+1] - longest_glons[i]
diffs = diffs[:-1]


"""
after having already looked at it, i know that dlon ~= diffs[0] 
(or tbh any element in diffs because they're all the same to within 9 decimal points)

did the checking using:
    for i in range(len(diffs)):
        diffs[i] = round(diffs[i],9)
    np.all(diffs == diffs[0]) 
        ^ if this is True, means that after rounding to 9 decimal points, all the elements in diffs are ==

"""

dlon = diffs[0]


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
        

plt.imshow(h1)


#initialize the dataframe
Rice_df_raw = pd.DataFrame(columns=columns)
for r in Rice_data_raw_list:
    #assign position in the database based on 'cloud index'
    ind = r[0]-1  #'cloud index' starts at 1, but python starts counting at 0
    
    #then add each cloud as a row
    Rice_df_raw.loc[ind] = r


Rice_df_raw #okay now the quantities are in but not their units ...


#drop Nans along the mass column
Rice_df = Rice_df_raw.dropna(subset=['mass'])


columns


Rice_df


# ## empirical relations


columns


masses = (Rice_df['mass'].to_numpy())*u.M_sun #in Msun
sizes = (Rice_df['size'].to_numpy())*u.pc #in pc

logM = np.log10(masses.value)
logR = np.log10(sizes.value)

#there is a cloud with size=0pc


# #idea -- but need to fiddle with numpy, for figuring out empirical relationship 
# model = np.polyfit(logR, logM, 2)


def linear_fit (x, m, b):
    y = m*x + b
    
    return y


#set up the radius array to be nice input for the fit
fit_x = logR.copy()
#so linear fit has an ordered input 
fit_x.sort()


#curious to see if this matches the values reported by Rice

total_mass = np.nansum(masses)
total_mass


MRfit_y = linear_fit(fit_x, 2.2, 1.5)
#by eye got to 2.5, pretty close! MD find 2.2


## M vs R

#----data----
plt.scatter(logR,logM,s=5,label='data', alpha=0.7)

plt.plot(fit_x,MRfit_y, label = 'fit', color='red', linestyle='--')

plt.ylabel(r'log Masses [$M_\odot$]')
plt.xlabel('log R [pc]')
plt.title('Mass-Radius')
plt.legend()
plt.show()


plt.hist(masses, bins=5000)

plt.xscale('log')

plt.xlabel(r'Masses [$M_\odot$]')
plt.show()


plt.hist(sizes, bins=100)

plt.xlabel(r'Sizes [pc]')
plt.show()


clouds = Rice_df['cloud index']
vir_pars = Rice_df['virial_alpha']


plt.hist(vir_pars, bins=200,edgecolor='white')

plt.xlim(0,20)
plt.xlabel(r'Virial $\alpha$s')
plt.show()


# ### post-conversation with enrico
# something to explain this^^ either, the masses that we are inferring through the beam is too small
# or the clouds are unbound (because the sigmas that we are seeing are > the escape velocities)
# 
# another thing  to do is: start with virial parameter == 1 and infer the mass
# 


## 2nd Larson's Law
plt.scatter(clouds,vir_pars,s=5)
plt.hlines(1,0,1100, color='red',linestyle='--', label=r'$\alpha$=1')
plt.ylim(0,20) #outliers upto 70
plt.xlabel('Cloud #s')
plt.ylabel(r'Virial $\alpha$')


plt.legend()
plt.show()


#interesting sinusoidal shape although that 
#could just be from how things are sorted by index


line_width = (Rice_df['v_rms'].to_numpy())*u.km/u.s #not sure if this is the right thing
logSigma = np.log10(line_width.value)


RSig_fit_y = linear_fit(fit_x,0.5, -0.35)


## 1st Larson's Law

#----data----

plt.scatter(logR,logSigma,s=10,label='data')


plt.plot(fit_x, RSig_fit_y, color='red',linestyle='--',label='fit')

plt.ylabel('Velocity Dispersion [km/s]')
plt.xlabel('R [pc]')
plt.title('Linewidth-Size relation')

plt.legend()
plt.show()


# https://github.com/brkent/datacube-display/blob/main/animateNGC2841.gif
# may come in handy at some pt in the future
