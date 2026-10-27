#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
import datacube as dc
import numpy as np


hi4pi = fits.open('HI4PI.fits')[0]


lb = np.nansum(hi4pi.data,axis=0)
lv = np.nansum(hi4pi.data,axis=1)
bv = np.nansum(hi4pi.data,axis=2)


hilons = np.arange(4323)*(-0.083) + 180.166
hilats = np.arange(1200)*(0.083) -49.99
hivels = np.arange(156)*1.288 - 100.


lb_neg = np.nansum(hi4pi.data[:155,:,:],axis=0)
lb_pos = np.nansum(hi4pi.data[155:,:,:],axis=0)


plt.figure(figsize=(20,8))
plt.imshow(lb_neg,origin='lower',norm='log',extent=[180,-180,-50,50])
plt.contour(hilons,hilats,lb_neg,colors='white',levels=[200,400,800])


plt.figure(figsize=(20,8))

plt.imshow(lb_pos,origin='lower',norm='log',extent=[180,-180,-50,50])
plt.contour(hilons,hilats,lb_pos,colors='white',levels=[200,400,800])


plt.figure(figsize=(20,8))
plt.imshow(lv,origin='lower',norm='log',extent=[180,-180,-200,200])
plt.contour(hilons,hivels,lv,colors='white')


firsts = np.array((-200,-160,-120,-80,-40,0,40,80,120,160,200))
seconds = np.array((-180,-140,-100,-60,-20,20,60,100,140,180))


firsts_int = ((firsts+200.)/1.288).astype('int')
seconds_int = ((seconds+200.)/1.288).astype('int')


first_slices = [(firsts_int[i],firsts_int[i+1]) for i in range(len(firsts)-1)]
second_slices = [(seconds_int[i],seconds_int[i+1]) for i in range(len(seconds)-1)]


for i,ii in enumerate(second_slices):
    thislb = np.nansum(hi4pi.data[ii[0]:ii[1],:,:],axis=0)
    plt.figure(figsize=(20,8))

    plt.imshow(thislb,origin='lower',norm='log',extent=[180,-180,-50,50])
    plt.contour(hilons,hilats,thislb,colors=['red','purple'],
                levels=[100,400,800],linewidths=[0.8])
    plt.title(str(hivels[ii[0]])+' to '+str(hivels[ii[1]-1]))
    plt.show()
    

for i,ii in enumerate(first_slices):
    thislb = np.nansum(hi4pi.data[ii[0]:ii[1],:,:],axis=0)
    plt.figure(figsize=(20,8))

    plt.imshow(thislb,origin='lower',norm='log',extent=[180,-180,-50,50])
    plt.contour(hilons,hilats,thislb,colors=['red','purple'],
                levels=[100,400,800],linewidths=[0.8])
    plt.title(str(hivels[ii[0]])+' to '+str(hivels[ii[1]-1]))
    plt.show()
    

hi4pi = hi4pi[:-1,:,:]


shorter_ = np.array(([-100,-70],[-90,-60],[-70,-40],[-60,-30],\
          [-40,-10],[-30,0],[-10,10],[0,30],[10,40],\
          [30,60],[40,70],[60,90],[70,100]))


shorter_slices = ((shorter_+100)/1.288).astype('int')


hi4pi.shape


fits.PrimaryHDU(hi4pi).writetoteto('HI4PI_shorter.fits',overwrite=True)


for i,ii in enumerate(shorter_slices):
    thislb = np.nansum(hi4pi[ii[0]:ii[1],:,:],axis=0)
    plt.figure(figsize=(20,8))

    plt.imshow(thislb,origin='lower',norm='log',extent=[180,-180,-50,50])
    plt.contour(hilons,hilats,thislb,colors=['red','purple'],
                levels=[100,400,800],linewidths=[0.8])
    plt.title(str(hivels[ii[0]])+' to '+str(hivels[ii[1]-1]))
    plt.show()
    

from scipy import ndimage


0.125/0.083


coer_arser_hi4pi = ndimage.zoom(hi4pi,(1,0.083/0.125,0.083/0.125))


plt.imshow(np.nansum(coer_arser_hi4pi,axis=0),norm='log')


def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and 
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndnew = ndarray.reshape(new_shape)
#     ndarray = ndarray.reshape(flattened)
#     for i in range(len(new_shape)):
#         op = getattr(ndarray, operation)
#         ndarray = op(-1*(i+1))
    return ndnew


coarser_hi = ndimage.gaussian_filter(hi4pi,sigma=(1,0.125/0.083,0.125/0.083),\
                                    mode=('nearest','nearest','wrap'))


coarserer = np.resize(coarser_hi,(155,800,2881))


plswork = np.repeat(np.repeat(coarser_hi,5,axis=1),3,axis=2)[:,::2,::2]


plswork.shape


plt.imshow(np.nansum(plswork,axis=0),origin='lower',vmax=600,vmin=0,cmap='Paired')
# plt.xlim(700,1250)
# plt.ylim(100,310)
plt.colorbar()


plt.imshow(np.nansum(coarser_hi,axis=0),origin='lower',vmax=600,vmin=0,cmap='Paired')
plt.xlim(700,1250)
plt.ylim(100,310)
plt.colorbar()


plt.imshow(np.nansum(hi4pi,axis=0),origin='lower',vmax=600,vmin=0,cmap='Paired')
plt.xlim(700,1250)
plt.ylim(100,310)
plt.colorbar()
