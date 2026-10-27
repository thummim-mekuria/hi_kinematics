#!/usr/bin/env python
# coding: utf-8


import pandas as pd
from astropy.io import fits
from astropy.table import QTable
import numpy as np
import matplotlib.pyplot as plt
import datacube as dc
import interpolate as ip


cat = pd.read_csv('solar_neighborhood_clouds.csv',delim_whitespace=True)


dht = fits.open('DHT_CO.fits')[0]


# clouds near the galactic plane / in danger of breaking reid (start with |l|<45)
# - ophiuc x
# - corona x
# - aquilr 
# - lupus  
# - 
# 


np.concatenate([np.arange(9),np.arange(8)])


5//3


np.concatenate(([7],np.arange(9)))


aqr y 167
aqs y 248 249 217
cma 
cal y
cam y
cph 
cfr 
gem y
her y
lac y
lps   271 281 238 198 293
mdl  
mob y
mr2 y
oph y
ora y
orb y
orl y
per y
pol y
ros 
tau y
xbn y


columns=['name','lon','lat','dist']
data = [['AqR',np.array((35,11)),np.array((4,14)),484],
['AqS',np.array((42,28)),     np.array((-21,-11)), 133],
['CMa',np.array((-133,-140)), np.array((-6,2) )  , 1209],
['Cal',np.array((170,155)) ,  np.array((-13,-6) ), 470],
['Cam',np.array((156,140))   ,np.array((15,26)) ,  213],
['Cph',np.array((120,95) ) ,  np.array((7,25))   , 352],
['CFr',np.array((120,95))   , np.array((7,25))   , 923],
['Chm',np.array((-70,-50))   ,np.array((-20,-12)), 183],
['CrA',np.array((4,-3.5) )  ,np.array( (-25,-16)), 151],
['Gem',np.array((-163,-172)),np.array((-5,5))  ,  1786],
['Her',np.array((49,40) )  , np.array( (6,13) )  , 227],
['Lac',np.array((98,93)  ) ,np.array( (-17,-8)) , 503],
['Lps',np.array((-10,-28) )  ,np.array((2,24))   , 189],
['Mdl',np.array((-138,-147)) ,np.array((-5,3) )  ,2072],
['MOB',np.array((-153,-165) ),np.array((-2,5) )  , 745],
['MR2',np.array((-142,-150)) ,np.array((-15,-8)) , 778],
['Oph',np.array((8,-12) )   , np.array((12,23) )  ,144],
['OrA',np.array((-141,-157)) ,np.array((-22,-17)) ,432],
['OrB',np.array((-149,-160) ),np.array((-17,-10)), 423],
['OrL',np.array((-155,-174)) ,np.array((-20,-5))  ,402],
['Per',np.array((163,155))  , np.array((-25,-15)) ,294],
['Pol',np.array((130,118))   ,np.array((20,32))   ,352],
['Ros',np.array((-147,-160)),np.array((-6,3)),1304],
['Tau',np.array((177,163)),np.array((-21,-11)),141],
['Xbn',np.array((-138,-142)),np.array((-13,-6)),886]]

zcat = QTable.from_pandas(pd.DataFrame(data,columns=columns))


zcat.to_pandas().to_csv('back2zucker_cat.csv',index=False)


help(zcat.write)


zcat.write('back2zucker_cat.csv',format='csv',overwrite=True)


dht = fits.open('DHT_CO.fits')[0]


dht2 = fits.open('COGAL_all_mom.fits')[0]


vbl = dht2.data.astype('float32').transpose().swapaxes(1,2)


interp_dht = ip.interpolate_vlb(vbl)


w2 = np.nansum(interp_dht,axis=0)*1.3


w = np.nansum(dht.data,axis=0)*1.3


plt.figure(figsize=(20,6))
plt.imshow(w,origin='lower',cmap=dc.cmap_min0('viridis'),vmax=45,vmin=2,extent=[180,-180,-30,30])
    

plt.colorbar()
plt.xlim(40,-40)


for i in range(len(zcat)):
    plt.figure(figsize=(8,8))
    plt.imshow(w,extent=[180,-180,-30,30],\
               origin='lower',cmap=dc.cmap_min0('viridis'),vmin=1,vmax=45)
    plt.xlim(zcat['lon'][i])
    plt.ylim(zcat['lat'][i])
    plt.title(zcat['name'][i])
    plt.colorbar()
    plt.show()


#clean up the dht data , plot all the clouds -- 0th moment, velocity spectra, etc
#then start calculating masses / surface density profiles / virial paramteres / maybe larson's laws
