#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import clouds as cld
import astropy.units as u
import astropy.constants as c


def fix_colnames(filename,num):
    df = pd.read_csv(filename,skiprows=num,delim_whitespace=True,header=None)
#     print(df.head())
    head = pd.read_csv(filename,nrows=num-1).loc[1].to_list()[0].split()[1:] 
    #asssume the format of the file is summary.prt from BDC
    
    df.columns = head
    
    return df


import decomp_DHT_2 as dd2


dd2


src = dd2.sources
DHT_CO = dd2.DHT_CO


src.columns = [src.columns.to_list()[i][0] for i in range(len(src.columns))]


src[['lon','lat']] = round(src[['lon','lat']],2)
src['vel'] = round(src['vel'],1)


sm = fix_colnames('summary.prt',4)


src.columns.to_list()[1:4] = sm.columns[0:3]


new_cols = src.columns[1:4].to_list()+sm.columns[3:8].to_list()+\
           (' '.join(sm.columns[8:].to_list()).upper()).split() #want to differentiate column names between d1&d2


sm.columns = new_cols


mg = sm.merge(src,on=['lon','lat','vel'])
neighbors = mg[mg['Dist.'] <=3]


ind_lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


lvs_w_ds = [ind_lvs[i] for i in range(179) if i in mg['!id']]


len(lvs_w_ds)


len(mg)


dd2.DHT_CO


noo = DHT_CO.clone()


noo


noo.cube = lvs_w_ds[-1]


np.all(np.nan_to_num(noo.cube) == np.nan_to_num(lvs_w_ds[-1]))


mg['Dist.'][70]


mg['Dist.'].to_numpy()[70]


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
    dAsky = mask * np.ones_like(mask) 
    A_tot_ang = np.nansum(dAsky) * (clob.dl * clob.db)#.to(u.rad**2)
    print(np.nansum(dAsky))
    
    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        
    # return {'surface density per pixel': dSig, 'angular area per pixel': dAsky_signal}
    return dSig, A_tot_ang

# ######### END: SURFACE DENSITY #########


mg


lf_nums = [dd2.d[i].idx for i in range(len(dd2.d)) if dd2.d[i].is_leaf]


lf_nums


pd.read_csv('COGAL_decomposed_catalog.csv')['area_exact']/(DHT_CO.db*DHT_CO.dl)


pd.read_csv('COGAL_decomposed_catalog_wclust.csv')['area_exact']


As_1*u.pc**2


As - As_1


aaaa = (CO_clouds_pckd['area_exact'].to(u.rad**2).value)*(1e3*u.pc*(CO_clouds_pckd['Dist.'].value))**2


ughhh = mg.iloc[25:]


okayy = mg.iloc[:25]


o_npix = [dd2.d[i].get_npix() for i in range(len(dd2.d)) if dd2.d[i].idx in okayy['_idx']]
u_npix = [dd2.d[i].get_npix() for i in range(len(dd2.d)) if dd2.d[i].idx in ughhh['_idx']]


o_mks = [dd2.d[i].get_mask() for i in range(len(dd2.d)) if dd2.d[i].idx in okayy['_idx']]
u_mks = [dd2.d[i].get_mask() for i in range(len(dd2.d)) if dd2.d[i].idx in ughhh['_idx']]


oarea = np.nansum(np.array(o_mks),axis=1)
uarea = np.nansum(np.array(u_mks),axis=1)


onuma = oarea


(oarea/oarea).shape


o_mks[0].shape


import datacube as dc


for i in ughhh['_idx']:
    plt.imshow(np.nansum(dd2.d[i].get_mask(),axis=0),cmap=dc.cmap_grey0('viridis'))
    plt.show()


plt.figure(figsize=(20,10))
for i in ughhh['_idx']:
    plt.imshow(np.nansum(dd2.d[i].get_mask(),axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower')
plt.show()


plt.figure(figsize=(20,10))
for i in okayy['_idx']:
    plt.imshow(np.nansum(dd2.d[i].get_mask(),axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower')
plt.show()


for i in okayy['_idx']:
    plt.imshow(np.nansum(dd2.d[i].get_mask(),axis=0),cmap=dc.cmap_grey0('viridis'))
    plt.show()


plt.pcolormesh((o_mks[14][30]),cmap=dc.cmap_grey0('viridis'))


# california? = 165
# perseus? = 212
# orionb? 244


mg.loc[115]


plt.hist(o_npix,density=True,bins=10,alpha=0.7)
plt.hist(u_npix,density=True,bins=20,alpha=0.7)
plt.yscale('log')
plt.show()


plt.hist(okayy['lat'],density=True,bins=10,alpha=0.7)
plt.hist(ughhh['lat'],density=True,bins=20,alpha=0.7)
plt.show()


plt.hist(okayy['lon'],density=True,bins=10,alpha=0.7)
plt.hist(ughhh['lon'],density=True,bins=20,alpha=0.7)
plt.show()


plt.hist(okayy['Dist.'],density=True,bins=10,alpha=0.7)
plt.hist(ughhh['Dist.'],density=True,bins=20,alpha=0.7)
plt.show()


for i in range(172):
    if np.abs((aaaa[i]/1e6 - As[i]*u.pc**2/1e6)/(aaaa[i]/1e6)) > 0.5:
        plt.scatter(i,(aaaa[i]/1e6 - As[i]*u.pc**2/1e6)/(aaaa[i]/1e6),marker='+')
#         plt.xlim(-1,30)


dd2.catalog


Ms_1*u.kg


Ms*u.kg


plt.hist(np.log10(np.sqrt(As/np.pi)),bins=30)#,color='blue',alpha=0.5)
plt.hist(np.log10(np.sqrt(As[neighbors.index]/np.pi)),bins=30)#,color='orange',alpha=0.5)

plt.xlabel(r'log$_{10}$ (effective radius [pc])')
plt.show()


plt.hist(np.log10(Ms*units['M'].to(c.M_sun)),bins=30)
plt.hist(np.log10(Ms[neighbors.index]*units['M'].to(c.M_sun)),bins=30)

plt.xlabel(r'log$_{10}$ (mass [$M_\odot$])')

plt.show()


plt.hist(np.log10(R_As),bins=30)
plt.hist(np.log10(R_As[neighbors.index]),bins=30)

plt.xlabel(r'log$_{10}$ (accumultion radii [pc])')

plt.show()


plt.scatter(np.log10(np.sqrt(As/np.pi)),np.log10(R_As),marker = '+')
plt.plot(np.arange(0,3.5,0.1),0.8 + 0.67*np.arange(0,3.5,0.1),color='red',linestyle='--',label='fit: 0.8 + 0.67R')
plt.scatter(np.log10(np.sqrt(As[neighbors.index]/np.pi)),\
            np.log10(R_As[neighbors.index]),marker = '+')


plt.ylabel('log(accumulation radius [pc])')
plt.xlabel('log(effective radius [pc])')
plt.legend()
plt.show()


cat = dd2.cat

cat_lon = cat['l_c__deg'].to_numpy()
for l,ll in enumerate(cat_lon):
    if ll > 180:
        cat_lon[l] -= 360
        
cat_lat = cat['b_c__deg'].to_numpy()


neigh_lvs = [lvs_w_ds[i] for i in range(len(lvs_w_ds)) if i in neighbors['!id']]


lf_nums = np.array([dd2.d[i].idx for i in range(len(dd2.d)) if dd2.d[i].is_leaf])


lf_nums


lf_idx_w_ds = lf_nums[mg['!id'].astype('int').to_list()]


lf_idx_w_ds


tbl_lfs_w_ds


#when i come back - match lon / lat with x_cen/y_cen to see if it makes sense to merge mg & tbl the way i did


mg


tbl = pd.read_csv('COGAL_decomposed_catalog.csv')
tbl_lfs_w_ds = tbl.loc[lf_idx_w_ds]


mg['_idx'] = lf_idx_w_ds
mg_n_tbl = mg.merge(tbl_lfs_w_ds, on='_idx')


mg_n_tbl['M'] = Ms
mg_n_tbl['A'] = As
mg_n_tbl['R_A'] = R_As


mg_n_tbl.columns


from astropy.table import Table


CO_clouds_pckd[115]


np.sqrt((CO_clouds_pckd[115]['A'] * u.pc**2)/np.pi)


(ppvAs[115].value*((u.deg**2).to(u.rad**2)))*(1e3*u.pc*CO_clouds_pckd[115]['Dist.'])**2


np.sqrt((5529.8738*u.pc**2)/np.pi)


myAs


ppvAs = ppvAs.value*((u.deg**2).to(u.rad**2))*(1e3*u.pc*CO_clouds_pckd['Dist.'].value)**2


def calcr(a):
    r = np.sqrt(a/np.pi)
    return r


diff_r = (100*np.abs(calcr(myAs) - calcr(ppvAs))/calcr(myAs))


dfdiff = pd.DataFrame(diff_r)


biggest_a_diff = [80,88,124,150,59]


dfdiff.sort_values(0)


ughhh


for i in range(3):
    print(dd2.d[i])


for st in dd2.d


for bad in bad_idx:
    struct = dd2.d[bad]
    
    msk = struct.get_mask()
    w = np.nansum(msk,axis=0)
    m = (w!=0)
    mskdd = np.ones_like(m) * m 
    datot = np.nansum(mskdd) * (DHT_CO.db*DHT_CO.dl).to(u.rad**2)
#     print(np.nansum(m))
    print(datot)
    
    idcs = struct.indices()
#     print(len(idcs[0]))
    aridcs = np.array(idcs)
    t_idcs = aridcs.transpose()
    sp = aridcs[1:]
    tsp = sp.transpose()
#     print(sp.transpose()[0])
    no_dbs = np.unique(tsp,axis=0)
    n = no_dbs.shape[0]
    adeg = n*DHT_CO.db*DHT_CO.dl
#     print(n)
    print(adeg.to(u.rad**2))
    print('\n')


lvs =[dd2.d[i] for i in range(len(dd2.d)) if dd2.d[i].is_leaf]


ind_lvs = []
for lf in lvs:
    m = lf.get_mask()
    mskd = m*DHT_CO.cube
    ind_lvs.append(mskd)
    

CO_clouds_pckd['area_exact'][biggest_a_diff]


CO_clouds_pckd['A'][biggest_a_diff]


lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


np.where(arxt==83.4375)


189,202,255,283,160


from astrodendro import Dendrogram


dend = Dendrogram.load_from("COGAL_dendrogram.fits")


len(dend)


catalog = pd.read_csv('COGAL_decomposed_catalog.csv')


catalog['structs'] = dend


for i in range(309):
    if dend[i].is_leaf:
        print(i)


ctwcl = pd.read_csv('COGAL_decomposed_catalog_wclust.csv')


pd.read_csv('COGAL_decomposed_catalog.csv') 


lvs


ctwcl


lvs = [dend[i] for i in range(len(dend)) if dend[i].is_leaf]


len(np.unique(np.array(lvs[-1].indices()[1:]).transpose(),axis=0))


np.nansum(np.ones_like(np.nansum(lvs[-1].get_mask(),axis=0)) * (np.nansum(lvs[-1].get_mask(),axis=0) != 0 ))


cb = fits.open('COGAL_decomposed_leaves.fits')[0].data[-1]


w1 = np.nansum(cb*1.3,axis=0)


m1 = (w1 !=0)


m = w!=0


np.nansum(np.ones_like(m1) * m1)


cb


np.nansum(np.ones_like(m) * m)


catalog


for i in [189,202,255,283,160]:
    print(dend[i])


lvs = [dend[i].idx for i in range(len(dend)) if dend[i].is_leaf]


np.array(lvs)[[80,88,124,150,59]]


sig = DHT_CO.calculate_rms_noise()
minn = np.nanmin(DHT_CO.cube)

kw_bool = {'save_all_leaves':True}
d = compute_dend(minn,20*sig,500)


arxt = pd.read_csv('COGAL_decomposed_catalog.csv')['area_exact']


len(lvs)


1857, 2
 (1693, 2)
(1090, 2)
(1193, 2)
 (5340, 2)


help(np.unique)


ughhh.loc[biggest_a_diff]


biggest_a_diff


ughhh.loc[biggest_a_diff]


bad_idx = ughhh.loc[[80,88,124,150,59]]['_idx'].to_numpy()


bad_idx


ppvAs[biggest_a_diff]


myAs[biggest_a_diff]


plt.hist(diff_r,bins=50)
# plt.xlim(-1,20)
plt.show()


CO_clouds_pckd = mg_n_tbl[['!id', '_idx','lon', 'lat', 'vel', 'Dist.', 'Integrated', 'Arm',
       'area_ellipse', 'area_exact', 'flux', 'major_sigma', 'minor_sigma',
       'position_angle', 'radius', 'v_cen', 'v_rms', 'x_cen', 'y_cen',
       'M', 'A', 'R_A']]

CO_clouds_pckd = Table.from_pandas(CO_clouds_pckd)


CO_clouds_pckd['Dist.'].unit = 1e3*u.pc


CO_clouds_pckd['M'].unit = u.kg
CO_clouds_pckd['A'].unit = u.pc**2
CO_clouds_pckd['R_A'].unit = u.pc


CO_clouds_pckd['vel'].unit = u.km/u.s


CO_clouds_pckd['lon'].unit = u.deg
CO_clouds_pckd['lat'].unit = u.deg


for name in CO_clouds_pckd.columns:
    if name in dd2.catalog.columns:
        CO_clouds_pckd[name].unit = dd2.catalog[name].unit 


CO_clouds_pckd.write('co_catalog.csv',format='ascii.ecsv',overwrite=True) 
#want the 'e' otherwise units arent saved


myAs = np.array(CO_clouds_pckd['A'].tolist())*u.pc**2


ppvAs = np.array(CO_clouds_pckd['area_exact'].tolist())*u.deg**2


pd.read_csv('co_catalog.csv',nrows=28,skiprows=4)#this is where the unit info goes


pd.read_csv('co_catalog.csv',delim_whitespace=True,skiprows=32)


#from here, 1: use accumulation radii to find HI parents 
#&& try to do background /continuum subtraction to make closed contours emerge


Table.write.help()
