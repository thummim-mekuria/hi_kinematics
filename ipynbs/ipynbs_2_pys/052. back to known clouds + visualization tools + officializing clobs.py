#!/usr/bin/env python
# coding: utf-8


#01.19 meeting notes
#re-evaluate HI velocity dispersion
#the main point is not to get all the known gmc clouds ever, but just to choose a sample, 
#then base the selection of the non-gmc hi clouds based on ^
#it will be helpful to build in ways to: inspect spectra, show 0th moment maps (in a way where you can 
#overlap the CO & HI data, select different velocity ranges in real time, etc)
#deliverable for next meeting: ^ 
#also for next time: revisit the timelines for the paper, 
#within the next two meetings, have the techniques nailed down for
#selecting the gmc-full & gmc-less hi clouds
#in order to move towards quantitatively analyzing the different populations
#add a column for surface density (in cat)


#look at the literature for why the calculated distances for clouds near l~0 might be 
#larger than the others


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
from astrodendro import Dendrogram,periodic_neighbours


def display(thing_3d,order=0,cmp = dc.cmap_grey0('viridis'),**kwargs):
    if type(thing_3d) == u.Quantity:
        thing_3d = thing_3d.value
    plt.imshow(np.nansum(thing_3d,axis=order),origin='lower',cmap=cmp,**kwargs)
    plt.show()


DHT_CO.data_header


DHT_CO.data_header.update(CTYPE2='GLAT-CAR',CTYPE3='GLON_CAR')


DHT_CO.data_header.update(NAXIS1=len(DHT_CO.vels),NAXIS2=len(DHT_CO.lats),NAXIS3=len(DHT_CO.lons))


DHT_CO.data_header['UNIT1'] = 'km/s'


DHT_CO.data_header['UNIT3'] = 'deg'


np.nanmin(DHT_CO.cube)


np.nanmax(DHT_CO.cube)


hdu = fits.PrimaryHDU(data=DHT_CO.cube,header=DHT_CO.data_header)


hdu.writeto('DHT_CO.fits')


DHT_CO.data_header.update(CDELT2=0.125,CDELT3=-0.125)


DHT_CO.data_header.update(CRVAL1=-49.4)


DHT_CO.data_header['CRVAL1'] = DHT_CO.vels[0]
DHT_CO.data_header['CRVAL2'] = DHT_CO.lats[0]
DHT_CO.data_header['CRVAL3'] = DHT_CO.lons[0]


DHT_CO.vels[0],DHT_CO.vels[1],DHT_CO.vels[2],


-49.4+(1.3*75)


DHT_CO.lats[-1]


-30.125+(0.125*480)


DHT_CO.data_header


DHT_CO.data_header.comments['CRVAL3'] = 'value at pixel 0'


DHT_CO.data_header.comments['UNIT3'] = 'axis 3 unit fourth'


DHT_CO.data_header.insert('CTYPE2',('CRVAL1',DHT_CO.data_header['CRVAL1']),after=True)


hdu=fits.PrimaryHDU(data=DHT_CO.cube,header=DHT_CO.data_header)


hdu.writeto('DHT_CO.fits',overwrite=True)


HI4PI.cube.shape


np.nanstd(DHT_CO.cube)


DHT_CO.data_header


help(DHT_CO.data_header.rename_keyword)


DHT_CO.data_header.update(BITPIX=32)


DHT_CO.data_header


HI4PI.data_header.update(BITPIX=32)


HI4PI.data_header.remove('DATASUM')


hduhi = fits.PrimaryHDU(data=HI4PI.cube,header=HI4PI.data_header)


hduhi.writeto('HI4PI.fits')


HI4PI.cube.dtype


DHT_CO.cube.dtype


fits.open('COGAL_all_mom.fits')[0].data.dtype


fits.open('HI4PI_CAR.fits')[0].data.dtype


HI4PI.lons[0]


HI4PI.data_header


HI4PI.data_header


spectral_cube.SpectralCube.read('DHT_CO.fits')


import spectral_cube


'Hi'.isupper()


help(list.insert)


comments = []
for i,ii in enumerate(hdr.split()):
    if ii =='/':
        for j,jj in enumerate(hdr.split()[i:]):
            if jj.isupper():
                comment = hdr.split()[i:j+i+1]
                print(hdr.split()[i:j+i+1])
                comments.append(' '.join(comment))
                break


h = fits.Header.fromstring('NAXIS = 1 \n NAXIS1 = 3',sep='\n')


h['NAXIS1']


hs = hdr.split()


hs.insert(1,2)


hs = hdr.split()


for i in range(len(hs)-1,-1,-1):
    ii = hs[i]
    
    if ii =='=':
        hs.insert(i,'      ')
        hs.insert(i-1,'\n')
        
#         for j,jj in enumerate(hdr.split()[i:]):
#             if jj.isupper():# and jj!='END' and hs[j+i+1] == '=':
#                 print(jj)
# #                 hs.insert(j+i+1,'\n')
# #                 comment = hdr.split()[i:j+i+1]
# #                 print(hdr.split()[i:j+i+1])
# #                 comments.append(' '.join(comment))
#                 break


hs


hdr_str_hopefullythisoneworks=' '.join(hs[1:-1])


hdr_str_hopefullythisoneworks


co_head = fits.Header.fromstring(hdr_str_hopefullythisoneworks,sep='\n')


co_head


help(list.insert)


' '.join(hs)


help(fits.Header.fromstring)


fits.Header.fromstring(hdr)


dir(fits.open('DHT_CO.fits')[0].header)


spectral_cube.SpectralCube.read('DHT_CO.fits')


help(spectral_cube.SpectralCube.read)


DHT_CO.data_header['CRVAL1']


DHT_CO.data_header.update(UNIT1='km/s')


DHT_CO.data_header.tofile('head.txt')


hdstr = pd.read_csv('head.txt',header=None).to_string(index=False)


fits.Header.fromfile('head.txt')


np.array(hdr.split())


hdr


hdr = '    '.join(hdr)


hdr = hdstr.split()[1:]


help(DHT_CO.data_header.tofile)


DHT_CO.data_header


hdstr = DHT_CO.data_header.tostring()


hdstr.split(' ')


hdr_list = hdstr.split()


hdr_ar = np.array(hdr_list)


hdr_ar.reshape(43,4)


len(hdr_ar)


172/8


dir(DHT_CO.data_header)


DHT_CO = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')

DHT_CO.cube = DHT_CO.cube[209:285,:,:] #trim to narrower velocity range
DHT_CO.vels = DHT_CO.vels[209:285]

smooth_whole = ip.interpolate_vlb(DHT_CO.cube)

wo_gal_cent = np.zeros_like(smooth_whole) #interpolate bc 'mom' is staccatoed
wo_gal_cent[:,:,:641]        = smooth_whole[:,:,:641] #take out |b|<5 for |l|<100
wo_gal_cent[:,:,2241:]       = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]

DHT_CO.cube = wo_gal_cent


HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
trimmed_HI4PI_cube = HI4PI.cube[317:627,482:1682,:]
HI4PI.cube = trimmed_HI4PI_cube.astype('float32')
HI4PI.vels = HI4PI.vels[317:627]
HI4PI.lats = HI4PI.lats[482:1682]


dnd = Dendrogram.load_from('COGAL_dendrogram.fits')
cat = pd.read_csv('physcl_prprts_catalog.csv')


mekcat = pd.read_csv('Mekuria_catalog.csv',skiprows=2)


lvs_map = fits.open('lvs_ind_map.fits')[0].data


cat_strut = [dnd[i] for i in range(len(dnd)) if i in cat['_idx'].to_numpy()] #dend structures found in cat


cat_masks = [c.get_mask() for c in cat_strut]


cat_masks = np.array(cat_masks)


cat_masks_collapsed = np.nansum(cat_masks,axis=0)


just_lvs = cat_masks_collapsed*DHT_CO.cube


just_lvs.shape


hdu = fits.PrimaryHDU(just_lvs)


hdu.writeto('structures_w_ds.fits',overwrite=True)


plt.imshow(cat_masks_collapsed[40],cmap=dc.cmap_grey0('viridis'),origin='lower')


display(cat_strut[0].get_mask())


m = cat_strut[0].get_mask()
v,b,l = np.where(m!=0)
v0,v1 = np.min(v),np.max(v)
b0,b1 = np.min(b),np.max(b)
l0,l1 = np.min(l),np.max(l)


for st in cat_strut:
    m = st.get_mask()
    v,b,l = np.where(m!=0)
    v0,v1 = np.min(v),np.max(v)
    b0,b1 = np.min(b),np.max(b)
    l0,l1 = np.min(l),np.max(l)
    display(m[v0:v1,b0-5:b1+5,l0-5:l1+5])


for cl in mekcat.values:
    print(cl,len(cl))


cld.ind_ang('get angular',DHT_CO.vels,[v0,v1])


v0,v1


mekcat['ehm'] = [None]*len(mekcat)


def find_cloud_lvs(cloud_cat,lf_cat,dnd,clob):
    #this only finds the structures that are at the central coordinate of the named cloud
    #how to improve so it can find structures along the whole cloud? (eg california, chamaleon)
    
    struts = [dnd[i] for i in range(len(dnd)) if i in lf_cat['_idx'].to_numpy()] 
    data_locs = [np.where(strt.get_mask()!=0) for strt in struts]
    
    
    phys_coords = [(cld.ind_ang('get angular',clob.vels,[np.min(d[0]),np.max(d[0])]),\
      cld.ind_ang('get angular',clob.lats,[np.min(d[1]),np.max(d[1])]),\
      cld.ind_ang('get angular',clob.lons,[np.min(d[2]),np.max(d[2])]))
     for d in data_locs] 
    
    matched = [(k,cat['_idx'][i],i) for k in range(len(cloud_cat)) \
               for i in range(len(phys_coords)) \
               if phys_coords[i][1][0]<=cloud_cat['b_c<deg>'][k]<=phys_coords[i][1][1]\
               and phys_coords[i][2][1]<=cloud_cat['l_c<deg>'][k]<=phys_coords[i][2][0]]
    return matched
#     cloud_cat['leaves']=[None]*len(cloud_cat)
#     for i,cl in enumerate(cloud_cat.values):
#         lc,bc = cl[1],cl[2]
#         associates = []
    
#         for st in struts:
#             m = st.get_mask()
#             v,b,l = np.where(m!=0)
#             v0,v1 = np.min(v),np.max(v)
#             b0,b1 = np.min(b),np.max(b)
#             l0,l1 = np.min(l),np.max(l)

#             vmn,vmx = cld.ind_ang('get angular',clob.vels,[v0,v1])
#             bmn,bmx = cld.ind_ang('get angular',clob.lats,[b0,b1])
#             lmx,lmn = cld.ind_ang('get angular',clob.lons,[l0,l1])
            
#             if (lmn<=lc<=lmx) and (bmn<=bc<=bmx):
#                 associates.append(st.idx)
                
#         cloud_cat['leaves'][i] = [associates]
        
#         print(cl[0],associates)   
#     return cloud_cat
                

#                 print(name,str(st.idx))
#                 display(m[v0:v1,b0-5:b1+5,l0-5:l1+5])
            

matches = find_cloud_lvs(mekcat,cat,dnd,DHT_CO) #list comprehension version


data_locs = [np.where(strt.get_mask()!=0) for strt in cat_strut]
    
phys_coords = [(cld.ind_ang('get angular',DHT_CO.vels,[np.min(d[0]),np.max(d[0])]),\
  cld.ind_ang('get angular',DHT_CO.lats,[np.min(d[1]),np.max(d[1])]),\
  cld.ind_ang('get angular',DHT_CO.lons,[np.min(d[2]),np.max(d[2])]))
 for d in data_locs] 


matched_vbl_locs = [(m[0],phys_coords[m[2]]) for m in matches]


np.where(mekcat['name']=='Orion_A')


mekcat['name']


matches


[i for i in range(len(matches)) if matches[i][0] == 0]


len(phys_coords)


np.where(cat['_idx']==316)


def plot_this(cloud,how='0',**kwargs):
    if type(cloud) == str:
        this = np.where(mekcat['name']==cloud)[0][0] 
        #assumes theres only one element ie only one row in the database with this name 
        ns = [i for i in range(len(matches)) if matches[i][0] == this]
        for n in ns:
    #         print(n)
            ar = matched_vbl_locs[n]
            vels = ar[1][0]
            lats = ar[1][1]
            lons = ar[1][2]

            idx = matches[n][1]

            v0,v1 = cld.ind_ang('get index',DHT_CO.vels,vels)
            b0,b1 = cld.ind_ang('get index',DHT_CO.lats,lats)
            l0,l1 = cld.ind_ang('get index',DHT_CO.lons,lons)

            m = dnd[idx].get_mask() * DHT_CO.cube
            print(name)
            lims = {'0':[lons[0],lons[1],lats[0],lats[1]],\
                    '1':[lons[0],lons[1],vels[0],vels[1]],\
                    '2':[lats[0],lats[1],vels[0],vels[1]]}
            if len(how)==1: 
                #eventually add an else that lets multiple subplots when e.g. how = (0,1)
                #^ will need to edit display or make a display_multiple
                
                display(m[-2+v0:v1+2,-2+b0:b1+2,-2+l0:l1+2],order=int(how),extent=lims[how],**kwargs)

    if type(cloud) == int: #assume what is being given is an idx
        this = np.where(cat['_idx']==cloud)[0][0]
        ar = phys_coords[this]
        print(ar)
        
        vels = ar[0]
        lats = ar[1]
        lons = ar[2]
        
        v0,v1 = cld.ind_ang('get index',DHT_CO.vels,vels)
        b0,b1 = cld.ind_ang('get index',DHT_CO.lats,lats)
        l0,l1 = cld.ind_ang('get index',DHT_CO.lons,lons)

        m = dnd[cloud].get_mask() * DHT_CO.cube
        
        lims = {'0':[lons[0],lons[1],lats[0],lats[1]],\
                    '1':[lons[0],lons[1],vels[0],vels[1]],\
                    '2':[lats[0],lats[1],vels[0],vels[1]]}
        if len(how)==1: 
            display(m[-2+v0:v1+2,-2+b0:b1+2,-2+l0:l1+2],order=int(how),extent=lims[how],**kwargs)


for a,ar in enumerate(matched_vbl_locs):
    nm    = mekcat['name'][ar[0]]
    vels = ar[1][0]
    lats = ar[1][1]
    lons = ar[1][2]
    
    idx = matches[a][1]
    
    v0,v1 = cld.ind_ang('get index',DHT_CO.vels,vels)
    b0,b1 = cld.ind_ang('get index',DHT_CO.lats,lats)
    l0,l1 = cld.ind_ang('get index',DHT_CO.lons,lons)
    
    m = dnd[idx].get_mask() * DHT_CO.cube
    print(nm)
    display(m[-2+v0:v1+2,-2+b0:b1+2,-2+l0:l1+2],extent=[lons[0],lons[1],lats[0],lats[1]])
    display(m[-2+v0:v1+2,-2+b0:b1+2,-2+l0:l1+2],order=1,extent=[lons[0],lons[1],vels[0],vels[1]])
    

find_cloud_lvs(mekcat,cat,dnd,DHT_CO) #nested for loops version


for st in cat_strut:
    m = st.get_mask()
    v,b,l = np.where(m!=0)
    v0,v1 = np.min(v),np.max(v)
    b0,b1 = np.min(b),np.max(b)
    l0,l1 = np.min(l),np.max(l)
    
    vmn,vmx = cld.ind_ang('get angular',DHT_CO.vels,[v0,v1])
    bmn,bmx = cld.ind_ang('get angular',DHT_CO.lats,[b0,b1])
    lmx,lmn = cld.ind_ang('get angular',DHT_CO.lons,[l0,l1])
    
    for cl in mekcat.values:
        name,lc,bc,d = cl
        if (lmn<=lc<=lmx) and (bmn<=bc<=bmx):
            print(name,str(st.idx))
            display(m[v0:v1,b0-5:b1+5,l0-5:l1+5])


aquila: 217(102),248(121),249(122) -- d: x2
california: 173(75),133(52), 136(53) -- d:x0.5
chamaeleon: 218(103),216(101) -- d:xBIG
corona: 275(145) -- d:xBIG
crossbones:304(171) -- d:1-1.5x
hercules: 257(129) -- d:1x
lacerta: 152(64) -- d:1.5x
monr2: 301(168) -- d:1.5x
ophiuchus: 260(131) -- d:xBIG
ora: 246(120) -- d:x2
orb: 250(123) -- d:x2    
perseus: 214(100) -- d:x1
polaris: 157(66) -- d:1-1.5x
scissors: 202(93)
taurus: 264(134) -- d:x1.5

      
cam
cma_ob1
gem_ob1

lupus
maddalena
mon_ob1 - 189?
orion_lam
rosette


cat.loc[64]


mekcat


lon_ids = cld.ind_ang('get index',DHT_CO.lons,mekcat['l_c<deg>'].to_numpy())
lat_ids = cld.ind_ang('get index',DHT_CO.lats,mekcat['b_c<deg>'].to_numpy())


for i,b, l in zip(mekcat.index,lat_ids,lon_ids):
    cb = dnd.data[:,b-80:b+80,l-80:l+80]
    plt.imshow(np.nansum(cb,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'))
    plt.title(mekcat['name'][i])
#     plt.scatter(l,b,color='orange',marker='x')
    plt.show()
#     sp = dc.spectrum1D(cb)
#     plt.plot(sp)
#     plt.show()
#     print(b,l,DHT_CO.lats[b],DHT_CO.lons[l])


lvs.shape


lims_lb_CO = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]


lims_lv_CO = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.vels[0],DHT_CO.vels[-1]]


plt.figure(figsize=(20,5))
plt.imshow(np.nansum(lvs,axis=1),origin='lower',cmap=dc.cmap_grey0('viridis'))#,extent=lims_lv_CO)
# plt.scatter(mekcat['l_c<deg>'],mekcat['b_c<deg>'],color='orange')
plt.xlim(0,500)


# plt.figure(figsize=(20,5))
plt.figure(figsize=(10,10))
plt.imshow(lvs[42],origin='lower',\
           cmap=dc.cmap_grey0('Paired'),extent=lims_lb_CO,vmin=190,vmax=290)
plt.colorbar()
plt.scatter(mekcat['l_c<deg>'],mekcat['b_c<deg>'],color='red',marker='+')
plt.xlim(180,150)


# plt.figure(figsize=(20,5))
plt.figure(figsize=(10,10))

plt.imshow(np.nansum(lvs,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims_lb_CO)
plt.scatter(mekcat['l_c<deg>'],mekcat['b_c<deg>'],color='orange')
plt.xlim(180,150)


plt.scatter(cat['lon<deg>'],cat['lat<deg>'],marker='+',c=cat['_idx'],cmap='Paired')
plt.colorbar()
plt.scatter(mekcat['l_c<deg>'],mekcat['b_c<deg>'],marker='x',color='black')
plt.xlim(180,150)
# plt.colorbar()


def gaussian(xs, peaks, FWHMs, means):
    """Return a Gaussian function."""
    gs = []
    for p,fwhm,mn in zip(peaks,FWHMs,means):
        sigma = fwhm / 2.354820045  # (2 * sqrt( 2 * ln(2)))
        g = p * np.exp(-(xs - mn)**2 / 2. / sigma**2)
        gs.append(g)
        
    return gs


def gauss_fit(spec,thresh=None):
    from scipy.signal import find_peaks,peak_widths
    if thresh==None:
        thresh=np.nanstd(spec)
        
    N = len(spec)
    xs = np.arange(N)
    
    fnd_pk_out = find_peaks(spec,height=thresh)
    pk_xs = fnd_pk_out[0] #peak positions along x
    pk_ys = fnd_pk_out[1]['peak_heights'] #peak positions along y
    
    pk_wdth_out = peak_widths(spec,pk_xs)
    wdths = pk_wdth_out[0] #fwhms
    
    decomposed = gaussian(xs,pk_ys,wdths,pk_xs)
    decomposed = np.array(decomposed)
    
    nump_comps = len(decomposed)
    
    return decomposed
