#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import datacube as dc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import QTable
import fun as fns


lons = np.arange(180,-180,-0.125)
lats = np.arange(-50,50,0.125)
vels = np.arange(-100,101,1.3)
coords = np.array(np.meshgrid(lons,lats,vels)).transpose().reshape(-1,3).transpose()
df = pd.DataFrame(coords.transpose(),columns=['lons','lats','vels'])
df['vunc'] = np.full(len(df),0.6)
df['pfar'] = np.full(len(df),0.5)
df['!names']= np.arange(len(df))
df = df[['!names','lons', 'lats', 'vels', 'vunc', 'pfar']]
df.loc[0] = ['!','','','','','']
df.to_csv('sources_info.inp',index=False)


pd.read_csv('reid/sources_info.inp',nrows=100)


hi4pi,dht = fits.open('both_same_grid_720d.fits')


zcat = QTable(fits.open('final_cat.fits')[1].data)


cal_hi = hi4pi.data[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1],
                 zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1],
                 zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]]


cal_co = dht.data[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1],
                 zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1],
                 zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]]


lons = np.arange(360,-360,-0.125)
lats = np.arange(-50,50,0.125)
vels = np.arange(-100,101,1.3)


# plt.imshow(np.nansum(cal_hi,axis=0),origin='lower',cmap='viridis')
plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=0),levels=np.arange(200,3000,200))
plt.colorbar()


plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            np.nansum(cal_co,axis=0),levels=np.arange(5,3000,5),
            colors='white',linewidths=0.5,alpha=0.8)

plt.xlim(172,152)


# plt.imshow(np.nansum(cal_hi,axis=0),origin='lower',cmap='viridis')
plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=1),levels=np.arange(200,10000,1000))
plt.colorbar()


plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_co,axis=1),levels=np.arange(5,3000,20),
            colors='white',linewidths=0.5,alpha=0.8)

plt.xlim(172,152)


# plt.imshow(np.nansum(cal_hi,axis=0),origin='lower',cmap='viridis')
plt.contourf(lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=2),levels=np.arange(200,15000,1000))
plt.colorbar()


plt.contour(lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_co,axis=2),levels=np.arange(5,3000,40),
            colors='white',linewidths=0.5,alpha=0.8)


# plt.imshow(np.nansum(cal_hi,axis=0),origin='lower',cmap='viridis')


plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            np.nansum(cal_co,axis=0),levels=np.arange(0,30,2),
            )

plt.colorbar()

plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=0),linewidths=0.5,alpha=0.8,
            levels=np.arange(200,3000,150),colors='white')

plt.xlim(172,152)


# plt.imshow(np.nansum(cal_hi,axis=0),origin='lower',cmap='viridis')


plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_co,axis=1),levels=np.arange(0,300,20),
            )

plt.colorbar()

plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=1),linewidths=0.5,alpha=0.8,
            levels=np.arange(200,9000,500),colors='white')

plt.xlim(172,152)


# plt.imshow(np.nansum(cal_hi,axis=2),origin='lower',cmap='viridis')

plt.contourf(lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_co,axis=2),levels=np.arange(0,200,30),
            )

plt.colorbar()

plt.contour(lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            vels[zcat['v_pix_HI'][0][0]:zcat['v_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=2),linewidths=0.5,alpha=0.8,
            levels=np.arange(200,10000,700),colors='white')

# plt.xlim(172,152)


# plt.imshow(fns.get_mom1(cal_co,zcat['v_pix_GMC'][0]),origin='lower',cmap='viridis',norm='asinh')


plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            fns.get_mom1(cal_co,zcat['v_pix_HI'][0]),levels=np.arange(-10,10,2),
            cmap='YlGnBu')

plt.colorbar()


plt.xlim(172,152)


# plt.imshow(fns.get_mom1(cal_co,zcat['v_pix_GMC'][0]),origin='lower',cmap='viridis',norm='asinh')

plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            fns.get_mom1(cal_hi,zcat['v_pix_HI'][0]),levels=np.arange(-20,10,3),
            cmap='YlGnBu')

plt.colorbar()


plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            fns.get_mom1(cal_co,zcat['v_pix_HI'][0]),levels=np.arange(-10,10,2),
            cmap='YlGnBu',alpha=0.8)

plt.xlim(172,152)


# plt.imshow(fns.get_mom1(cal_co,zcat['v_pix_GMC'][0]),origin='lower',cmap='viridis',norm='asinh')


plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            fns.get_momN(cal_co,zcat['v_pix_HI'][0],2),levels=np.arange(0,15,1),
            cmap='viridis')

plt.colorbar()


plt.xlim(172,152)


# plt.imshow(fns.get_mom1(cal_co,zcat['v_pix_GMC'][0]),origin='lower',cmap='viridis',norm='asinh')


plt.contourf(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            fns.get_momN(cal_hi,zcat['v_pix_HI'][0],2),levels=np.arange(5,15,1),
            cmap='viridis')

plt.colorbar()

plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            fns.get_momN(cal_co,zcat['v_pix_HI'][0],2),levels=np.arange(0,9,1),
            colors='white',alpha=0.8,linewidths=0.5)


plt.xlim(172,152)


co_spcs = np.array([np.ravel(cal_co[:,b:b+1,l:l+1]) for b in range(169) for l in range(169)])
hi_spcs = np.array([np.ravel(cal_hi[:,b:b+1,l:l+1]) for b in range(169) for l in range(169)])


spc_coords = np.array([((zcat['b_pix_HI'][0][0]+b)*0.125 - 50,
                           (zcat['l_pix_HI'][0][0]+l)*-0.125 + 360) 
                          for b in range(169) for l in range(169)])


len(hi_spcs[4])


from scipy.signal import find_peaks


def calc_dvtv(array):
    n = len(array)
    dvtv = np.array([array[0]]+[(array[i+1] - array[i-1])/2 
                                for i in range(1,n-1)]+[array[-1]])
    return dvtv
    

def find_hisa(hi_arr,co_arr):
    hi_1dvtv = calc_dvtv(hi_arr)
    hi_2dvtv = calc_dvtv(hi_1dvtv)
    
    co_1dvtv = calc_dvtv(co_arr)
    co_2dvtv = calc_dvtv(co_1dvtv)
    
    
    hi_valleys = find_peaks(-1*hi_arr)[0]
    co_peaks = find_peaks(co_arr)[0]
    
    hi2_peaks = find_peaks(hi_2dvtv)[0]
    co2_valleys = find_peaks(-1*co_2dvtv)[0]
    
    
    hi_inds = np.intersect1d(hi_valleys,hi2_peaks)
    co_inds = np.intersect1d(co_peaks,co2_valleys)
    
    hisa = np.intersect1d(hi_inds,co_inds)
    
    return hisa


len(find_hisa(hi_spcs[4],co_spcs[4]))


cal_hisa = [i for i in range(169**2) if len(find_hisa(hi_spcs[i],co_spcs[i])) >= 1 ]


pxHISA = np.zeros(len(spc_coords))


pxHISA[cal_hisa] = 1


plt.imshow(pxHISA.reshape(169,169),extent = [lons[zcat['l_pix_HI'][0][0]],
                                             lons[zcat['l_pix_HI'][0][1]],
                                             lats[zcat['b_pix_HI'][0][0]],
                                             lats[zcat['b_pix_HI'][0][1]]])
plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            np.nansum(cal_co,axis=0),levels=np.arange(5,3000,5),
            colors='white',linewidths=0.5,alpha=0.8)


plt.contour(lons[zcat['l_pix_HI'][0][0]:zcat['l_pix_HI'][0][1]],
            lats[zcat['b_pix_HI'][0][0]:zcat['b_pix_HI'][0][1]],
            np.nansum(cal_hi,axis=0),levels=np.arange(200,3000,200))

plt.imshow(pxHISA.reshape(169,169),extent = [lons[zcat['l_pix_HI'][0][0]],
                                             lons[zcat['l_pix_HI'][0][1]],
                                             lats[zcat['b_pix_HI'][0][0]],
                                             lats[zcat['b_pix_HI'][0][1]]])


spc_coords[cal_hisa]


[(plt.plot(co_spcs[i]/np.nanmax(co_spcs[i])),
 plt.plot(hi_spcs[i]/np.nanmax(hi_spcs[i])),
 plt.title(spc_coords[i]),
 plt.show())
 
for i in range(169**2) if len(find_hisa(hi_spcs[i],co_spcs[i])) >= 1 ]


[(plt.plot(co_spcs[i]/np.nanmax(co_spcs[i])),
 plt.plot(hi_spcs[i]/np.nanmax(hi_spcs[i])),
 plt.title(spc_coords[i]),
 plt.show())
 
for i in range(169**2)]


#different scenarios
#their peaks are aligned
