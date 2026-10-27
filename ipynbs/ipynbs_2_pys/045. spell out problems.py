#!/usr/bin/env python
# coding: utf-8


import datacube as dc
import clouds as cld
import interpolate as ip
import bootleg as btlg

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import math

# %matplotlib ipympl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram


import decomp_DHT_2 as dd2


pd.read_csv('sources_info.inp',delim_whitespace=True)


pd.read_csv('sources_from_bfr.csv',delim_whitespace=True)


pd.read_csv('summary_from_bfr.csv')


dd2.lf_catalog


(dd2.lf_catalog['area_exact'].data)/(DHT_CO.db*DHT_CO.dl)


for i,lf in enumerate(dd2.lf_catalog['structure'].data):
    m = lf.get_mask()
    m1 = np.nansum(m,axis=0)
    m4 = (m1!=0)
    m2 = np.ones_like(m1) * m1
    m3 = m1/m1
    m5 = np.ones_like(m1) * m4
    print(np.nansum(m2),np.nansum(m3),np.nansum(m4),\
          np.nansum(m5),((dd2.lf_catalog['area_exact'].data)/(DHT_CO.db*DHT_CO.dl))[i])
#     print(np.nansum(m2)/np.nansum(m3))
#     print(m2.shape,m3.shape)


for i,lf in enumerate(dd2.lf_catalog['structure'].data):
    m = lf.get_mask()
    m1 = np.nansum(m,axis=0)
    m4 = (m1!=0)
    m2 = np.ones_like(m1) * m1
    m3 = m1/m1
    m5 = np.ones_like(m1) * m4
    print(np.nansum(m2),np.nansum(m3),np.nansum(m4),\
          np.nansum(m5),((dd2.lf_catalog['area_exact'].data)/(DHT_CO.db*DHT_CO.dl))[i])
#     print(np.nansum(m2)/np.nansum(m3))
#     print(m2.shape,m3.shape)


pd.re


mys = [mys[i].value for i in range(len(mys))]


plt.scatter(mys,ppvs,marker='+')


lfcat = pd.read_csv('COGAL_decomposed_catalog_lvs_only.csv')


dend = Dendrogram.load_from('COGAL_dendrogram.fits')


lvs = [dend[i] for i in range(len(dend)) if dend[i].is_leaf]


def surface_density (clob, data_type, n_sigma = 3, show_Sigma = False):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from create_subcube/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    '''

    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/astro-ph/0411089.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        rms_noise = 0.25 *u.K #mean calculated from values in Table 1 of the paper announcing the full DHT dataset (https://iopscience.iop.org/article/10.1086/318388/pdf)

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        rms_noise = 0.043 *u.K #from Table 1 of the paper announcing the full HI4PI dataset (https://arxiv.org/pdf/1610.06175.pdf)
    
   #assume that the clob has already been processed such that 
   #the noise (values < 3 sigma_rms) in each channel has been suppressed & there is only *one* peak along the velocity axis
   #this works fine for now because working with CO-bright clouds & we can specify what velocity range
   #might need to change for later if CO-dark clouds need this function
   
    #lets assume (at least for now) that the clob 
    #has already been processed through dendrogram.trunk[i].get_mask()

    W = np.nansum(clob.cube * clob.units['data'] * clob.dv, axis=0) #calculate integrated intensity
    N = W * conversion_factor #calculate column density from the integrated intensity & conversion factor
    dSig = mmw * (c.m_p.to(u.g)) * N #define surface density here
    
    mask = (W!=0)
    dAsky = mask * np.ones_like(W.value) 
    A_tot_ang = np.nansum(dAsky) * (clob.dl * clob.db).to(u.rad**2)

    
    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        kw_imshow = {'cmap': cmap0, 'origin':'lower'}

        dSig_Msunpc2 = dSig.to(c.M_sun/u.pc**2)

        fig, ax = plt.subplots (figsize = (15, 10), dpi = 200)
        Sig_plot = ax.imshow(dSig_Msunpc2.value, extent = [clob.lons[0], clob.lons[-1], \
                                                clob.lats[0], clob.lats[-1]], **kw_imshow)


        cbar = fig.colorbar(Sig_plot, ax=ax, label = Sig_label +r' [$M_\odot \ pc^{-2}$]', orientation = 'vertical')

        if data_type == '21cm':
            KMT_thresh = plt.contour(clob.lons, clob.lats, dSig_Msunpc2.value,  linewidths = 0.7, \
            levels=[10, 100], colors = ['xkcd:red', 'xkcd:orange']) #atomic-to-molecular transition for HI maps 
            cbar.add_lines(KMT_thresh)


        ax.set_xlabel(r'Galactic Longitude [$\degree$]')
        ax.set_ylabel(r'Galactic Latitude [$\degree$]')
        ax.set_title(Sig_label)


        plt.show()

    # return {'surface density per pixel': dSig, 'angular area per pixel': dAsky_signal}
    return dSig, A_tot_ang


src = pd.read_csv('sources_info.inp',delim_whitespace=True,skiprows=1)


sm = pd.read_csv('summary.prt',delim_whitespace=True,skiprows=2)


sm = sm.loc[1:]


src.columns = ['id','lon','lat','vel','pfar','vunc']


sm = sm[['!','Long.','Lat.','+/-']]


sm.columns = ['lon','lat','vel','dist']


src = src[['id','lon','lat','vel']]


sm


src['id'] = src['id'].astype('int')


src['lon'] = round(src['lon'],2)


src['lat'] = round(src['lat'],2)
src['vel'] = round(src['vel'],1)


src


help(sm.merge)


src.columns


sm.columns


sm['lon'] = sm['lon'].astype('float64')
sm['lat'] = sm['lat'].astype('float64')
sm['vel'] = sm['vel'].astype('float64')


mg = src.merge(sm,on=['lon','lat','vel'])


have_ds= [i for i in range(179) if src['id'].loc[i] in mg['id'].to_list()]


lvs_w_ds = np.array(lvs)[have_ds]


mg['dist'] = mg['dist'].astype('float64')


DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')

DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 

smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole)
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641]
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]

DHT_CO.cube = wo_gal_cent
#===========================


plt.scatter(np.arange(172),np.array(oldnums)-np.array(newnums),marker='+')


badold = oldnums[28:]


badnew = newnums[28:]


badold.sort()


badnew.sort()


cocat = pd.read_csv('co_catalog.csv',skiprows=26,delim_whitespace=True)


dend = Dendrogram.load_from('COGAL_dendrogram.fits')


import astropy.units as u


for n,num in enumerate(cocat['_idx']):
    lf = dend[num]
    m = lf.get_mask()
    m1 = np.nansum(m,axis=0)
    m2 = m1!=0
    m3 = np.nansum(m2)
    print(m3)
    mya = cocat['A'].loc[n]
    d = 1e3*cocat['Dist.'].loc[n]
    mya_rad = mya/(d**2)
    mya_deg = (mya_rad*u.rad**2).to(u.deg**2)
    mya_npx = mya_deg/((0.125*u.deg)**2)
    print(mya_npx)
    print('\n')


for i in range(144):
    if badold[i]!=badnew[i]:
        print(i,badold[i],badnew[i])


for i in range(172):
    print(i,oldnums[i],newnums[i])


newnums=[840,
661,
270,
145,
149,
318,
387,
302,
192,
816,
1663,
1912,
623,
256,
672,
123,
325,
1097,
289,
543,
462,
180,
212,
353,
429,
190,
767,
188,
601,
354,
226,
1703,
804,
201,
666,
498,
139,
165,
724,
359,
494,
178,
430,
323,
718,
195,
546,
1430,
554,
1037,
1250,
841,
511,
1202,
1102,
147,
313,
672,
400,
5340,
306,
685,
498,
711,
766,
933,
260,
137,
815,
267,
268,
1134,
249,
175,
317,
225,
257,
299,
322,
1857,
233,
227,
215,
197,
268,
258,
474,
1693,
236,
221,
619,
1161,
304,
314,
143,
952,
418,
1135,
676,
294,
703,
173,
269,
244,
1209,
268,
1007,
773,
1061,
1469,
807,
276,
1501,
851,
550,
1403,
924,
116,
545,
650,
137,
242,
433,
1486,
1090,
1176,
2977,
293,
773,
416,
247,
686,
2815,
919,
520,
197,
411,
436,
326,
670,
168,
291,
243,
123,
163,
308,
892,
565,
211,
480,
1193,
319,
628,
370,
180,
412,
527,
820,
209,
276,
358,
300,
167,
966,
122,
489,
560,
786,
163,
158,
277,
264]


oldnums = [840,
661,
270,
145,
149,
318,
387,
302,
192,
816,
1663,
1912,
623,
256,
672,
123,
325,
1097,
289,
543,
462,
180,
212,
353,
429,
190,
767,
188,
302,
601,
354,
226,
1703,
804,
201,
666,
498,
139,
165,
1092,
724,
359,
494,
178,
263,
430,
323,
718,
195,
546,
270,
1430,
554,
1037,
1250,
841,
511,
1202,
1102,
467,
147,
313,
672,
400,
5340,
306,
685,
498,
711,
766,
933,
260,
137,
815,
267,
268,
1134,
249,
175,
317,
246,
225,
257,
299,
322,
1857,
233,
227,
215,
197,
268,
258,
474,
1693,
236,
221,
619,
1161,
304,
314,
143,
952,
418,
1135,
676,
294,
703,
173,
269,
244,
1209,
268,
1007,
773,
1061,
1469,
807,
276,
1501,
851,
550,
1403,
1867,
924,
116,
545,
650,
137,
242,
433,
1486,
1090,
1176,
2977,
293,
773,
416,
247,
686,
2815,
919,
520,
197,
411,
436,
326,
670,
168,
291,
243,
123,
163,
308,
892,
565,
211,
480,
1193,
319,
628,
370,
180,
412,
527,
820,
209,
276,
358,
300,
167,
966,
122]


help(pd.Series())


from astropy.table import Table


help(Table.Column)


# %%time
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]

x,y = np.meshgrid(HI4PI.lons, HI4PI.lats)

DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


cocat = pd.read_csv('co_catalog.csv',delim_whitespace=True,skiprows=32)


pd.read_csv('co_catalog.csv',skiprows=4,nrows=28)


ang_R_A = (cocat['R_A']/(cocat['Dist.']*1e3) * np.cos(cocat['lat'] * u.deg.to(u.rad))) *\
            u.rad.to(u.deg) 


vrms_3x = 3*cocat['v_rms']


whi = np.nansum(HI4PI.cube,axis=0)


# fig,ax = plt.subplots(1,1,figsize=(15,10))

# lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]
# ax.imshow(whi,extent=lims)


# for l,v,r,vrms in zip(cocat['lon'],cocat['vel'],ang_R_A,4*vrms_3x):
#     ax.add_patch(Ellipse((l,v),2*r,vrms,fill=False,edgecolor='white',alpha=0.5))
    
    
# ax.set_ylim(-50,50)


# fig,ax = plt.subplots(1,1,figsize=(15,10))

# lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
# ax.imshow(whi,extent=lims)

# for l,b,r in zip(cocat['lon'],cocat['lat'],ang_R_A):
#     ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5))
# #     print(l,b,r)
    
    
# # ax.set_xlim(180,130)


parents = np.zeros_like(HI4PI.cube,dtype=bool)

for i in range(len(cocat)):#[130,69,121,118,64]:
    l = cocat['lon'][i]*u.deg #add unit
    b = cocat['lat'][i]*u.deg #add unit
    v = cocat['vel'][i]*u.km/u.s #add unit
    r = ang_R_A[i]*u.deg #add unit
    vrms = vrms_3x[i]*DHT_CO.dv #DHT px -> km/s
    parents = cld.add_accumulation_prism(v,b,l,r,vrms,parents)

# parents = HI4PI.cube * parents


# fig,ax = plt.subplots(1,1,figsize=(15,10))
# lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]

# ax.imshow(np.nansum(parents,axis=0),extent=lims,origin='lower')
# for l,b,r in zip(cocat['lon'],cocat['lat'],ang_R_A):
#     ax.add_patch(Circle((l,b),r,fill=False,edgecolor='white',alpha=0.5,linewidth=2))
# #     print(l,b,r)


# fig,ax = plt.subplots(1,1,figsize=(15,10))
# lims = [HI4PI.lons[0],HI4PI.lons[-1],1e-3*HI4PI.vels[0],1e-3*HI4PI.vels[-1]]

# ax.imshow(np.nansum(parents,axis=1),extent=lims,origin='lower')

# ax.set_ylim(-50,50)
# for l,v,r,vrms in zip(cocat['lon'],cocat['vel'],ang_R_A,vrms_3x):
#     #everything else is in physical units so convert vrms from px -> km/s
#     vrms = vrms*1e-3*HI4PI.dv.value
#     ax.add_patch(Ellipse((l,v),2*r,4*vrms,\
#                          fill=False,edgecolor='white',alpha=0.5,linewidth=2))
    
# # ax.set_xlim(80,50)
    

import cc3d
parents_masked_overlap,N = cc3d.connected_components(parents,connectivity=26,return_N = True)
#connectivity=6 and 18 also give the same N


from numba import njit

@njit
def separate_clouds(labels_map,data,N):
    x,y,z = data.shape
    out = np.zeros((N,x,y,z))
    for n in range(1,N+1):
        m = (labels_map==n)
        mskd = m*data
        out[n-1] = mskd
        
    return out


# #come back and figure out whats wrong with this / the function ^ 
# distinct_hi = separate_clouds(parents_masked_overlap,HI4PI.cube,N)


# lims = [HI4PI.lons[0], HI4PI.lons[-1], HI4PI.lats[0], HI4PI.lats[-1]]
# for i in range(N):
#     plt.imshow(np.nansum(HI4PI.cube*distinct_hi[i],axis=0),\
#                cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)
# #     plt.colorbar()
#     plt.title(str(i+1))
#     plt.show()


# lims = [HI4PI.lons[0], HI4PI.lons[-1], 1e-3*HI4PI.vels[0], 1e-3*HI4PI.vels[-1]]
# for i in range(N):
#     plt.imshow(np.nansum(HI4PI.cube*distinct_hi[i],axis=1),\
#                cmap=dc.cmap_grey0('viridis'),origin='lower',extent=lims)
# #     plt.colorbar()
#     plt.title(str(i+1))
#     plt.show()


# #once i figure out getting distinct_hi, save it as a fits
# HI_parents = fits.PrimaryHDU(distinct_hi)
# HI_parents.writeto('decomp_HI_parents_frm_RA.fits',overwrite=True)


fig,ax = plt.subplots(1,1,figsize=(15,10))

lims = [HI4PI.lons[0],HI4PI.lons[-1],HI4PI.lats[0],HI4PI.lats[-1]]
ax.imshow(whi,extent=lims,origin='lower')


for l,b,a_bb,r_A,ind in zip(cocat['lon'],cocat['lat'],cocat['area_exact'],ang_R_A,cocat.index):
    r_bb = np.sqrt(a_bb/np.pi)
#     r_A = r #1 degree -> 8pixels bc 0.125deg/pix
    if r_bb > r_A:
        ax.add_patch(Circle((l,b),r_A,fill=False,edgecolor='white'))
        ax.add_patch(Circle((l,b),r_bb,fill=False,edgecolor='red'))

#         print(ind)
# ax.set_xlim(180,130)


myAs = cocat['A'].to_numpy() * u.pc**2


ppvAs = ((u.pc*1e3*cocat['Dist.'].to_numpy())**2)*\
            (cocat['area_exact'].to_numpy()*(u.deg**2).to(u.rad**2))


dend = Dendrogram.load_from('COGAL_dendrogram.fits')


ind_lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


lvs_w_ds = [ind_lvs[i] for i in range(172) if i in cocat['!id']]
