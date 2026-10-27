#!/usr/bin/env python
# coding: utf-8


#this is mostly a testbed for developing .pys


# import decomp_DHT_2 as dd2
import numpy as np
import matplotlib.pyplot as plt
import bootleg as bdend
import clouds as cld
import pandas as pd


#when i come back -- find that one functino that reads in + cleans up summary .prt (maybe  in the drive trash)
#then use the distances to calculate M A etc, you know what to do after


def fix_colnames(filename,num):
    df = pd.read_csv(filename,skiprows=num,delim_whitespace=True,header=None)
#     print(df.head())
    head = pd.read_csv(filename,nrows=num-1).loc[1].to_list()[0].split()[1:] 
    #asssume the format of the file is summary.prt from BDC
    
    df.columns = head
    
    return df


sm = fix_colnames('summary.prt',4)


sm


src = pd.read_csv('sources_info.inp',delim_whitespace=True,header=[0,1])


src.columns = [src.columns.to_list()[i][0] for i in range(len(src.columns))]


src[['lon','lat']] = round(src[['lon','lat']],2)
src['vel'] = round(src['vel'],1)


src.columns.to_list()[1:4] = sm.columns[0:3]


new_cols = src.columns[1:4].to_list()+sm.columns[3:8].to_list()+\
           (' '.join(sm.columns[8:].to_list()).upper()).split() #want to differentiate column names between d1&d2


sm.columns = new_cols


sm


mg = sm.merge(src,on=['lon','lat','vel'])


import datacube as dc


plt.figure(figsize=(20,10))
plt.scatter(mg['lon'],mg['lat'],marker='+',s=1e2,c = mg['Dist.'],cmap=dc.cmap_grey0('RdYlGn_r'))
plt.colorbar()
plt.show()


ids_w_ds = mg['!id'].to_numpy().astype('int')


decomclus = fits.open('COGAL_decomposed_clustered.fits')[0].data


clusters_with_ds = [decomclus[i] for i in range(94) if i in mg['!id']]


from numba import njit


DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


mg


import astropy.units as u


units


neighbors = mg[mg['Dist.'] <=3]


mg[mg['Dist.'] <=3]


mg[mg['lon'] <-140]


plt.imshow(np.nansum(clusters_with_ds[16:17],axis=(1,0)),cmap=dc.cmap_grey0('rainbow'),origin='lower')
plt.xlim(2500,2880)


Table.write.help('pandas.csv')


from astropy.table import Table
dat = Table([[1, 2], [3, 4]], names=('a', 'b'))
# dat.write('table.dat', format='ascii')


lvs = fits.open('COGAL_decomposed_leaves.fits')[0].data


l55 = [dend[i].idx for i in range(309) if dend[i].is_leaf]


l55==lvs_nums


l22 = [dend[i] for i in range(309) if dend[i].is_leaf]


lvs_nums = [l22[i].idx for i in range(len(l22))]


ind_lvs=np.array(ind_lvs)


np.all(np.nansum(ind_lvs,axis=0) == lvs)


ind_lvs = []
for lnum in lvs_nums:
    mskd = lvs*(lvs==lnum)
    ind_lvs.append(mskd)
    plt.imshow(np.nansum(mskd,axis=0),cmap=dc.cmap_grey0('rainbow'))
    plt.show()


lvs


max([l22[i].get_npix() for i in range(179)])


nums =  pd.read_csv('cluster_nums.csv')['cluster nums'].to_list()


nums.sort()


len(l33),len(l44)


len(l44), len(l33)


plt.pcolormesh(lvs[20])
cb=plt.colorbar()


dir(cb)


plt.hist(lvs.flatten(),bins=500)
plt.yscale('log')
plt.xlim(30,50)


l44 = nums


l33 = dend_cat[dend_cat['parent']==-1]['_idx'].to_list()


for i in range(0,55,5):
    fig, axs = plt.subplots(1,2,figsize=(15,5))
    axs[0].pcolormesh(dend.index_map[i],cmap=dc.cmap_grey0('Paired'))
    axs[0].set_xlim(0,500)
    axs[0].set_ylim(100,300)
    axs[1].pcolormesh(lvs[i],cmap=dc.cmap_grey0('Paired'))
    axs[1].set_xlim(0,500)
    axs[1].set_ylim(100,300)
    plt.show()


dend_cat = pd.read_csv('COGAL_decomposed_catalog.csv')


dend_cat[dend_cat['structure_type']=='L']


len(decomclus),len(clusters_with_ds)


dend_cat[dend_cat['structure_type']=='T']['_idx']


dend_cat


dend.trunk


from astrodendro import Dendrogram


dend = Dendrogram.load_from('COGAL_dendrogram.fits')


c


# plt.figure(figsize=(15,3))
plt.pcolormesh(dend.index_map[35],cmap=dc.cmap_grey0('Paired'))
plt.xlim(2500,2880)


clus_labels = [i**2 for i in range(9)]


clus_labels = np.array(clus_labels)


clusters_pd = pd.DataFrame({'cluster nums':clus_labels})


clusters_pd


clusters_pd.to_csv('cluster_nums.csv')


clusters_pd = pd.DataFrame({'cluster nums':clus_labels})
clusters_pd.to_csv('cluster_nums.csv')


dir(dend)


dend_cat[dend_cat['structure_type']=='T']['_idx']


dat.write('table.fits', format='fits')


fits.open('table.fits')[1].data


pd.read_csv('table.csv')


Table.write.help()


fits.TableHDU(dat)


help(fits.TableHDU.writeto)


from scimes import SpectralCloudstering


help(SpectralCloudstering)


cat.sort_values('l_c__deg')


plt.scatter(mg[mg['Dist.'] >10]['lon'], mg[mg['Dist.'] >10]['lat'],\
            c=mg[mg['Dist.'] >10]['Dist.'],marker='+',cmap='Paired')
plt.colorbar()


Ms[neighbors.index]


plt.hist(np.log10(np.sqrt(As/np.pi)),bins=30)
plt.hist(np.log10(np.sqrt(As[neighbors.index]/np.pi)),bins=30)

plt.xlabel(r'log$_{10}$ (effective radius [pc])')
plt.show()


plt.hist(np.log10(Ms*units['M'].to(c.M_sun)),bins=30)
plt.hist(np.log10(Ms[neighbors.index]*units['M'].to(c.M_sun)),bins=30)

plt.xlabel(r'log$_{10}$ (mass [$M_\odot$])')

plt.show()


plt.hist(mg['Dist.'],bins=30)
plt.hist(neighbors['Dist.'],bins=30)

plt.xlabel('distance [kpc]')
# plt.xscale('log')
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


plt.scatter(np.log10(np.sqrt(As/np.pi)),np.log10(R_As),marker = '+')
plt.plot(np.arange(0,3.5,0.1),0.8 + 0.67*np.arange(0,3.5,0.1),color='red',linestyle='--',label='fit: 0.8 + 0.67R')
plt.scatter(np.log10(np.sqrt(As[neighbors.index]/np.pi)),\
            np.log10(R_As[neighbors.index]),marker = '+')


plt.ylabel('log(accumulation radius [pc])')
plt.xlabel('log(effective radius [pc])')
plt.legend()
plt.show()


plt.figure(figsize=(10,6))
# plt.scatter(mg['lon'],mg['lat'],marker='+')
plt.scatter(neighbors['lon'],neighbors['lat'],marker='x',c=R_As[neighbors.index])
plt.show()


len(clusters_with_ds)


neigh_clus = [clusters_with_ds[i] for i in range(91) if i in neighbors['!id'].astype('int')]


len(neigh_clus)


cat = pd.read_csv('Mekuria_catalog.csv',sep=',',skiprows=5)


for cl in neigh_clus:
    plt.imshow(np.nansum(cl,axis=0),origin='lower',cmap=dc.cmap_grey0('rainbow'))
    plt.show()


neighbors_cutout = cutout(neighbors['lon'].to_numpy(),neighbors['lat'].to_numpy(),np.nansum(neigh_clus,axis=0))


cat_lon = cat['l_c__deg'].to_numpy()
for l,ll in enumerate(cat_lon):
    if ll > 180:
        cat_lon[l] -= 360
        
cat_lat = cat['b_c__deg'].to_numpy()
  

cats_in_the_neighborhood = cutout(cat_lon,cat_lat,np.nansum(neigh_clus,axis=0))


len(neighbor_cat)


zero_neighbors = np.nansum(neigh_clus,axis=1)


zero_neighbors.shape


neighbors['lon'].to_numpy()[i]-5


zero_mg = np.nansum(clusters_with_ds,axis=1)


# for i in range(79):
#     lims = [neighbors['lon'].to_numpy()[i]-1, neighbors['lon'].to_numpy()[i]+1,\
#             neighbors['lat'].to_numpy()[i]-1,neighbors['lat'].to_numpy()[i]+1]
lms = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]
plt.imshow(np.nansum(zero_mg,axis=0),origin='lower',cmap=dc.cmap_grey0('rainbow'),extent=lms)
plt.imshow(np.nansum(zero_neighbors,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lms)

plt.title(str(i))
# plt.xlim(500,1000)
plt.show()


cat


for nc, cc in enumerate(cats_in_the_neighborhood):
    lims = [cat_lon[nc]-5, cat_lon[nc]+5, cat_lat[nc]-5,cat_lat[nc]+5]
    plt.imshow(np.nansum(cc,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims)
    plt.title(cat['name'][nc])
#     plt.scatter(neighbors['lon'].to_numpy()[nc],neighbors['lat'].to_numpy()[nc],color='red')
    plt.show()
#     plt.imshow(np.nansum(cat_n_cloud[c],axis=0),**kw)
#     plt.show()


for nc, cc in enumerate(neighbors_cutout):
    lims = [neighbors['lon'].to_numpy()[nc]-5, neighbors['lon'].to_numpy()[nc]+5,\
            neighbors['lat'].to_numpy()[nc]-5,neighbors['lat'].to_numpy()[nc]+5]
    plt.imshow(np.nansum(cc,axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'),extent=lims)
    
    plt.scatter(neighbors['lon'].to_numpy()[nc],neighbors['lat'].to_numpy()[nc],color='red')
    plt.show()
#     plt.imshow(np.nansum(cat_n_cloud[c],axis=0),**kw)
#     plt.show()


import pyreadstat


dir(pyreadstat)


neighbors


import astropy.constants as c


compute = dd2.compute_dend


clob = dd2.DHT_CO
cat = dd2.cat


cat.head()


sig = clob.calculate_rms_noise()
minn = np.nanmin(clob.cube)


# sig_sig_50 = compute(sig, sig,50)
# min_min_50 = compute(minn, minn,50)
# min_5sig_50 = compute(minn,5*sig,50)
# min_5sig_100 = compute(minn,5*sig,100)


def count_pxs(lst):
    return [lst[i].get_npix() for i in range(len(lst))]


lons = [cat['l_c__deg'][i] for i in range(18)]
for l, ll in enumerate(lons):
    if ll > 180:
        lons[l] -=360


lats = [cat['b_c__deg'][i] for i in range(18)]


# @njit
def cutout(ls,bs,data):
    output = []
    for i, l in enumerate(ls):
        b = bs[i]
        
        lsign = l/np.abs(l)
        bsign = b/np.abs(b)
        
        l0 = l-10
        l1 = l+10

        b0 = b-10
        b1 = b+10
                
        l0,l1 = cld.ind_ang('get index', clob.lons,[l0,l1])
        b0,b1 = cld.ind_ang('get index', clob.lats,[b0,b1])
        
        l0,l1 = np.sort((l0,l1))
        b0,b1 = np.sort((b0,b1))
            
        sub = data[:,b0:b1,l0:l1]
    
        output.append(sub)
    return output


cat_n_cloud = cutout(lons,lats,clob.cube)


kw = dd2.kw_imshow


min_20sig_500 = compute(minn,20*sig,500)


compute_c = dd2.compute_cluster


kw_bool = {'save_all_leaves':True}


m20s500 = compute_c(min_20sig_500)


m20s500 = compute_c(min_20sig_500,**kw_bool)


m20s500 = compute_c(min_20sig_500,**kw_bool)


bigm = min_20sig_500.trunk[0].get_mask()
for tr in min_20sig_500.trunk:
    msk = tr.get_mask()
    bigm += msk


bigm_l = min_20sig_500.leaves[0].get_mask()
for lf in min_20sig_500.leaves:
    msk = lf.get_mask()
    bigm_l += msk


dir(m20s500)


clus_labels = m20s500.clusters_asgn.data
clus_cube = m20s500.dendrogram.data


ind_clust = []
for kl in m20s500.clusters:
    mskd = clus_cube*(clus_labels==kl)
    ind_clust.append(mskd)
    
ind_clust = np.array(ind_clust)


np.array(ind_clust).shape


for ind in ind_clust:
    plt.imshow(np.nansum(ind,axis=0),**kw)
    plt.show()


clust_n_clouds_2 = cutout(lons,lats,m20s500.leaves_asgn.data)


clust_n_clouds = cutout(lons,lats,m20s500.clusters_asgn.data)


for cnum in m20s500.clusters:
    dat = m20s500.clusters_asgn.data
    mskd = dat*(dat==cnum)
    plt.imshow(np.nansum(mskd,axis=0),**kw)
    plt.show()


m20s500.clusters_asgn.data


for c,cc in enumerate(clust_n_clouds):
    plt.imshow(np.nansum(cc,axis=0),**kw)
    plt.show()
#     plt.imshow(np.nansum(cat_n_cloud[c],axis=0),**kw)
#     plt.show()


for c,cc in enumerate(clust_n_clouds):
    plt.imshow(np.nansum(cc,axis=0),**kw)
    plt.show()
#     plt.imshow(np.nansum(cat_n_cloud[c],axis=0),**kw)
#     plt.show()


for c,cc in enumerate(clust_n_clouds):
    plt.imshow(np.nansum(cc,axis=0),**kw)
    plt.show()
    plt.imshow(np.nansum(cat_n_cloud[c],axis=0),**kw)
    plt.show()
    plt.imshow(np.nansum(clust_n_clouds_2[c],axis=0),**kw)
    plt.show()


m20s500.sa


plt.imshow(np.nansum(m20s500.clusters_asgn.data,axis=0),**kw)
plt.colorbar()


plt.imshow(np.nansum(m20s500.leaves_asgn.data,axis=0),**kw)
plt.colorbar()


plt.imshow(np.nansum(m20s500.leaves_asgn.data,axis=0),**kw)
plt.colorbar()


plt.imshow(np.nansum(bigm_l*clob.cube,axis=0),**kw)
plt.colorbar()


plt.imshow(np.nansum(bigm*clob.cube,axis=0),**kw)
plt.colorbar()


for l in min_20sig_500.leaves:
    m = l.get_mask()
    plt.imshow(np.nansum(m*clob.cube,axis=0),**kw)
    plt.show()
#     plt.colorbar()


ind_clust.shape


@njit
def weighted_coords(data, vbl_axes, ax_ress): #prob very imprecise, since the clusters are spread out in l, b, v space
    #might need to do multiple distance calculations???

    N = len(data)
    IDs, ls, bs, vs = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
#     print(IDs,ls,bs,vs)
    for i in range(N):
        cb = data[i]
#         print(i)

        v,b,l = dc.weighted_vbl(cb,vbl_axes,ax_ress)
#         print(i, ' ',v,b,l)
        IDs[i] = i
        ls[i] = l
        bs[i] = b
        vs[i] = v


    output = np.array(IDs, vs, bs, ls)

    return output


cat


smol = ind_clust[23:26]


smol.shape


out = weighted_coords(ind_clust, (clob.vels, clob.lats, clob.lons),\
                      (clob.dv.value,clob.db.value,clob.dl.value))
