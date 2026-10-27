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


dnd = Dendrogram.load_from('COGAL_dendrogram.fits')
# cat = pd.read_csv('physcl_prprts_catalog.csv')


cat = pd.read_csv('physcl_prprts_catalog.csv')


hi_dv = hi4pi.header['CDELT1'] * u.Unit(hi4pi.header['CUNIT1'])
hi_nv = hi4pi.header['NAXIS3']
hi_0v = hi4pi.header['CRVAL1'] * u.Unit(hi4pi.header['CUNIT1'])

hi_db = hi4pi.header['CDELT2'] * u.Unit(hi4pi.header['CUNIT2'])
hi_nb = hi4pi.header['NAXIS2']
hi_0b = hi4pi.header['CRVAL2'] * u.Unit(hi4pi.header['CUNIT2'])

hi_dl = hi4pi.header['CDELT3'] * u.Unit(hi4pi.header['CUNIT3'])
hi_nl = hi4pi.header['NAXIS1']
hi_0l = hi4pi.header['CRVAL3'] * u.Unit(hi4pi.header['CUNIT3'])


hi4pi = fits.open('HI4PI.fits')[0]


Aquila_S [248, 249]
California [173]
Cam [182]
Chamaeleon [218]
CMa_OB1 [306]
Corona_Australis [275]
Crossbones [304]
Gem_OB1 [208]
Hercules [257]
Lacerta [152]
Lupus [271]
Maddalena [311, 313]
Mon_OB1 [189]
Mon_R2 [301]
Ophiuchus [242, 260]
Orion_A [246]
Orion_B [250]
Orion_Lam []
Perseus [214]
Polaris [157]
Rosette [226]
Serpens_AqR []
Taurus [264]


mek_idx = [249,173,182,218,306,275,304,208,257,152,271,311,189,301,260,246,250,214,157,226,264]


mek_idx.sort()


cat_lons = u.deg*cat['lon<deg>'].to_numpy()


cat_lats = u.deg*cat['lat<deg>'].to_numpy()


ls_c = np.array([int((lon - hi_0l)/hi_dl) for lon in cat_lons])


bs_c = np.array([int((lat - hi_0b)/hi_db) for lat in cat_lats])


ras = round(cat['accumulation_radius<deg>']/hi_db.value).astype('int').to_numpy()


cat_ras = u.deg*cat['accumulation_radius<deg>'].to_numpy()


mekcat = pd.read_csv('Mekuria_catalog.csv',skiprows=2)


hi4pi.data.dtype


hi4pi_cube = hi4pi.data.astype('float32')


rng = np.random.default_rng(0)
x = np.linspace(-5., 5., 200)
y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
y += rng.normal(0., 0.2, x.shape)


plt.plot(x,y)
plt.plot(x,g(x))


from astropy.modeling import models, fitting


g_init = models.Gaussian1D(amplitude=20., mean=0, stddev=10.)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, x, y)


g


g


g


g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)


g(spc)


from specutils.fitting import fit_lines


from specutils import Spectrum1D


one52 = Spectrum1D(spectral_axis=hi_vs*u.km/u.s,flux=spc*u.K)


g_init_q = models.Gaussian1D(amplitude=1.*u.K, mean=0*u.km/u.s, stddev=1.*u.km/u.s)


sg = fit_lines(one52,g_init_q)


sg(hi_vs*u.km/u.s)


for idx,l,b,r in zip(mek_idx,ls_c[mek_ns],bs_c[mek_ns],ras[mek_ns]):
    sub = hi4pi_cube[:,b-r:b+r,l-r:l+r]
    spc = dc.spectrum1D(sub)
    g = fit_g(g_init, hi_vs, spc)
    plt.plot(hi_vs,spc)
    plt.plot(hi_vs,sg(hi_vs*u.km/u.s))
    plt.title(str(idx))
    plt.show()
    break


mek_ns[-7]


ls_c[134], ras[134]


mek_ns=[i for i in range(186) if cat['_idx'].to_numpy()[i] in mek_idx]


cat['_idx'].to_numpy()[mek_ns]


hi4pi.header


hi_vs = 1.288*np.arange(310) - 199.572


dht = fits.open('DHT_CO.fits')[0]


lvs = [dnd[i] for i in range(len(dnd)) if dnd[i].is_leaf]


ms = [lf.get_mask() for lf in lvs]


bigm = np.nansum(np.array(ms),axis=0)


wm = np.nansum(bigm,axis=0)


plt.imshow(wm,origin='lower',cmap=dc.cmap_grey0('viridis'))


dht.header['CDELT2']


dht.header['NAXIS3']


2881*(-0.125) + 180.125


lims = [180.125,-180. ,-30.125,29.875]


wco = np.nansum(dht.data,axis=0)


for thing in mekcat.values:
    name, lc, bc, d = thing
    plt.imshow(wco,origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims)
    plt.xlim(lc+15, lc-15)
    plt.ylim(bc-15, bc+15)
    plt.title(name)
    plt.show()


for thing in mekcat.values:
    name, lc, bc, d = thing
    plt.imshow(wm,origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims)
    plt.xlim(lc+5, lc-5)
    plt.ylim(bc-5, bc+5)
    plt.title(name)
    plt.show()


newhdu = fits.PrimaryHDU(data=dht.data,header=newhead)


dht.header['NAXIS1']=76


dht.header['NAXIS3']=2881


newhead = dht.header


newhdu.writeto('DHT_CO.fits',overwrite=True)


px01 = dht.header['CRVAL1']*u.Unit(dht.header['CUNIT1'])
# px02
# px03
d1 = dht.header['CDELT1']*u.Unit(dht.header['CUNIT1'])
# d2
# d3


dir(dht)


dir(dnd)


mekcat


hd = fits.open('DHT_CO.fits')[0].header


lons = mekcat['l_c<deg>'].to_numpy()


lons


[int((lon - hd['CRVAL3'])/hd['CDELT3']) for lon in lons]


[(hd['CRVAL3'] + (l-1)*hd['CDELT3']) * u.Unit(hd['CUNIT3']) for l in lons]


hd['CRVAL1']


from astrodendro import Dendrogram


dnd = Dendrogram.load_from('COGAL_dendrogram.fits')


lvs_map = fits.open('lvs_ind_map.fits')[0].data


cat_strut = [dnd[i] for i in range(len(dnd)) if i in cat['_idx'].to_numpy()] #dend structures found in cat


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


import astropy.units as u


fits.open('DHT_CO.fits')[0].header['CDELT1'] * u.Unit(fits.open('DHT_CO.fits')[0].header['CUNIT1'])


lon_ids = cld.ind_ang('get index',DHT_CO.lons,mekcat['l_c<deg>'].to_numpy())
lat_ids = cld.ind_ang('get index',DHT_CO.lats,mekcat['b_c<deg>'].to_numpy())


from astrodendro import ppv_catalog


help(ppv_catalog)


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
