#!/usr/bin/env python
# coding: utf-8


# %%time

import datacube as dc
import clouds as cld
import transform as tr
import kd as kd

import numpy as np
# %matplotlib

import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy import units as u
import astropy.constants as c
from astrodendro import Dendrogram
from astropy.visualization import ZScaleInterval as zscale

from matplotlib.patches import Rectangle


zlims = zscale().get_limits
cmap0 = dc.cmap_grey0('viridis')
kw_imshow = {'cmap':cmap0, 'origin':'lower'}#,  'norm': mpl.colors.LogNorm()}


#----------------------
HI4PI = dc.unpack_fits_cube('HI4PI_CAR.fits', 'HI4PI')
HI4PI.data_header['T_rms'] = 43*1e-3 #u.K, from Table 1 of HI4PI paper: https://arxiv.org/pdf/1610.06175.pdf
LAB = dc.unpack_fits_cube('labh.fits','LAB')

DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')
DHT_CO2 = dc.unpack_fits_cube('COGAL_all_interp.fits','DHT')
DHT_CO3 = dc.unpack_fits_cube('COGAL_all_mom.fits','DHT')


DHT_CO.data_header['T_rms'] = 0.28 #u.K, median value from Table 1 of Dame paper: 


#----------------------
cat = pd.read_csv('Mekuria_catalog.csv',sep=',',skiprows=5)
names = cat['name']


l_unit = u.deg
b_unit = u.deg
d_unit = u.pc

#----------------------
def read_DHT(ns='all'):
    output={}
    if ns == 'all':
        ns = ['03','06','09','12','16','19','21','24','26','27','34','37','18']
        #for now, figure out DHT01 (bad data) 
        #& DHT18 (two ceph ds, beeg) 
        #& DHT08 (beeg) 
        #& DHT24, DHT34 (lam or & chamaeleon are disjointed down to the really low sig levels)
     
    for s in ns:
#         print(s)
        s_clob = dc.unpack_fits_cube('DHT'+s+'.fits','DHT')
        output['DHT'+s] = s_clob
        
    return output

#----------------------

surveys = read_DHT()
for s in surveys:
    k = np.where(cat['survey']==s)[0]

    these = cat['name'][k]
    
    nums = np.zeros(len(these))
    for i, this in enumerate(these):
        nums[i] = np.where(cat['name']==this)[0]
        
    nums = nums.astype(np.int64)
    cube = surveys[s]
    surveys[s] = [cube, nums]


ss = pd.read_csv('sources_info.inp',delim_whitespace=True)


ss['l [deg]'] = round(ss['l [deg]'],2)


ss['b [deg]'] = round(ss['b [deg]'],2)


ss['v_lsr [km/s]'] = round(ss['v_lsr [km/s]'],1)


ss


pd.read_csv('sources_info.inp',delim_whitespace=True)


sm = pd.read_csv('summary.prt',delim_whitespace=True,skiprows=2)


ss


mg = ss.merge(sm,on=['l [deg]', 'b [deg]'])#, 'v_lsr [km/s]'])


len(mg),len(sm)


pd.read_csv('summary.prt',delim_whitespace=True,skiprows=2)


len(ss),len(sm)


not_f = [i for i in range(len(ss)) if ss['l [deg]'][i] not in sm['l [deg]'].to_list()]


sm['l [deg]'][20],sm['b [deg]'][20], sm['v_lsr [km/s]'][20]


116.81 in ss['l [deg]'].to_list()


177*2


newss = ss[['!ID','l [deg]','b [deg]','v_lsr [km/s]','!v_unc']]


newss.to_csv('sources_info.inp',index=False)


ss['v_unc'] = np.full(len(ss),2.)


ss[['!ID','l [deg]','b [deg]','v_lsr [km/s]','v_unc','pfar']].to_csv('sources_info.inp',index=False,sep=' ')


np.full(len(ss),2.)


ss.loc[not_f]


ss.loc[not_f].to_csv('sources_info.inp')


ss['l [deg]'].to_list().index(123.5)


ss['l [deg]'][40],ss['b [deg]'][40], ss['v_lsr [km/s]'][40]


plt.figure(figsize=(20,5))

s = surveys['DHT18']
clob = s[0]
lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
plt.imshow(np.nansum(np.abs(s[0].cube),axis=0),**kw_imshow,extent=lims)
    
plt.xlim(230,130)
plt.ylim(-30,10)
plt.show()


plt.figure(figsize=(20,5))
for s in surveys.values():
    
    clob = s[0]
    lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
    plt.imshow(np.nansum(np.abs(s[0].cube),axis=0),**kw_imshow,extent=lims)
    
# plt.contour(DHT_CO.lons,DHT_CO.lats,np.nansum(np.abs(DHT_CO.cube),axis=0),colors='orange',linewidths=0.7)
# plt.contour(DHT_CO.lons+360,DHT_CO.lats,np.nansum(np.abs(DHT_CO.cube),axis=0),colors='orange',linewidths=0.7)

plt.xlim(230,130)
plt.ylim(-30,10)
plt.show()


plt.figure(figsize=(20,5))
for s in surveys.values():
    clob = s[0]
    lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
    plt.imshow(np.nansum(np.abs(s[0].cube),axis=0),**kw_imshow,extent=lims)
    
# plt.contour(DHT_CO.lons,DHT_CO.lats,np.nansum(np.abs(DHT_CO.cube),axis=0),colors='orange',linewidths=0.7)
# plt.contour(DHT_CO.lons+360,DHT_CO.lats,np.nansum(np.abs(DHT_CO.cube),axis=0),colors='orange',linewidths=0.7)

plt.xlim(230,130)
plt.ylim(-30,10)
plt.show()


plt.figure(figsize=(20,5))
lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]
plt.imshow(np.nansum(np.abs(DHT_CO.cube),axis=0),**kw_imshow,extent=lims)
plt.xlim(-130,-190)
plt.ylim(-30,10)
plt.show()


vmn,vmx=zlims(np.nansum(np.abs(DHT_CO.cube),axis=0))


DHT_CO.cube.shape


import interpolate as ip


ip.interpolate(dc.spectrum1D(DHT_CO.cube[:,240:241,960:961]))


plt.plot(dc.spectrum1D(DHT_CO.cube[:,240:241,960:961]))
plt.plot(ip.interpolate(dc.spectrum1D(DHT_CO.cube[:,240:241,960:961])))


dc.spectrum1D(DHT_CO.cube[:,320:321,1280:1281])


32*40


for x in range (32,35):
    plt.plot(dc.spectrum1D(DHT_CO.cube[:,10*x:10*x+1,40*x:40*x+1]),marker='x')
    plt.hlines(0,0,50)
    plt.show()


plt.figure(figsize=(20,5))
lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]
plt.imshow(np.nansum(np.abs(DHT_CO.cube),axis=0),origin='lower',\
           cmap=dc.cmap_grey0('Paired'))#,vmin=vmn,vmax=vmx)
# plt.xlim(130,125)
# plt.ylim(5,10)
# plt.xlim(140,120)
# plt.ylim(5,15)
plt.xlim(200,300)
plt.ylim(280,281)
plt.colorbar()
plt.show()


x1 = np.nansum(np.abs(DHT_CO.cube),axis=0)[280:281,200:300]


len(x1[0])


plt.figure(figsize=(10,10))
plt.scatter(np.arange(0,100),x1[0],marker='x',color='red',alpha=0.5)
plt.plot(np.arange(0,100),ip.interpolate(x1[0]),marker='+')


plt.figure(figsize=(10,10))
plt.scatter(np.arange(0,100),x1[0],marker='x',color='red',alpha=0.5)
plt.scatter(np.arange(0,100),ip.interpolate(x1[0]),marker='+')


x = np.abs(DHT_CO.cube[23:24,206:207,280:288])[0][0]


0.5*(0.15159798+0.00243568)


plt.scatter(np.arange(0,8),x)
plt.scatter(np.arange(0,8),ip.interpolate(x),alpha=0.3)


DHT_CO = dc.unpack_fits_cube('COGAL_local_mom.fits','DHT')


plt.scatter(np.arange(0,55),np.nanmean(DHT_CO.cube,axis=(1,2)),marker='x')
plt.scatter(np.arange(0,55),ip.interpolate1D(np.nanmean(DHT_CO.cube,axis=(1,2))),marker='+')


dht27 = surveys['DHT27'][0]


for k in range(len(DHT_CO.cube)):
    plt.imshow(np.abs(DHT_CO.cube[:,200:400,300:500][k]),**kw_imshow)
    # plt.ylim(100,200)
    # plt.xlim(50,150)
    plt.colorbar()
    plt.show()


plt.plot(np.nanmean(DHT_CO.cube[:,270:350,300:301],axis=1),marker='+')
# plt.plot(ip.interpolate1D(np.nanmean(DHT_CO.cube[:,300:301,313:314],axis=(1,2))),marker='+')
plt.show()


#pick up here -- currently working on a way to create a function to interpolate missing values
#first thing to work on when i come back: use max_gap to nail down how many indexes at each end of the line are just themselves
#the end goal is to use the local moment masked co dht cube, interpolate it myself, then when dendrogramming, can just do sigma? or not even noise threshold?


0 == 0.


ip.int


len(np.transpose(DHT_CO.cube[:,200:400,300:500]))


np.swapaxes(specsmooth,1,2).shape


vbl lbv 


len(smooth)


specsmooth = ip.interpolate_vlb(DHT_CO.cube[:,200:400,300:500])


smooth = ip.interpolate_vlb(DHT_CO.cube[:,200:400,300:500])


for k in range(len(smooth)):
    fig,axs=plt.subplots(1,2,figsize=(10,7))
    
    axs[0].imshow(smooth[k])#,cmap=dc.cmap_grey0('viridis'))
    axs[1].imshow(np.abs(DHT_CO.cube[:,200:400,300:500])[k])#,cmap=dc.cmap_grey0('viridis'))
    plt.show()
    # plt.xlim(170,200)
# axs[1].pcolormesh([k])


plt.imshow(np.nansum(smooth,axis=0),cmap=dc.cmap_grey0('viridis'))


plt.imshow(np.nansum(np.abs(DHT_CO.cube),axis=0),cmap=dc.cmap_grey0('viridis'))
plt.xlim(200,500)
plt.ylim(100,400)


smooth_whole = ip.interpolate_vlb(DHT_CO.cube)


388/2


plt.plot(np.nanmean(smooth_whole,axis=(1,2)))
plt.hlines(0.2,0,55)


plt.imshow(smooth_whole[30]*(smooth_whole[30]>0.2),cmap=dc.cmap_grey0('viridis'))


plt.imshow(smooth_whole[30],cmap=dc.cmap_grey0('viridis'))


DHT_CO.db


DHT_CO.dl


5/DHT_CO.dv


d_sm_wh_nolims = Dendrogram.compute(wo_gal_cent, 0.1, 0.1,50)


len(d_sm_wh_nolims)


d_sm_wh = Dendrogram.compute(wo_gal_cent, 0.1, 0.1, 125)


p = d_sm_wh.plotter()


np.nanstd(wo_gal_cent)


fig,ax = plt.subplots(1,1)
p.plot_tree(ax)
plt.hlines(0.1,0,700,color='red')
plt.ylim(0,3)


len(d_sm_wh.trunk)


npx = np.zeros(208)
for t,tr in enumerate(d_sm_wh.trunk):
    xx = tr.get_npix()
    npx[t] = xx


plt.imshow(np.nansum(smooth_whole*d_sm_wh[10].get_mask(),axis=2),cmap=dc.cmap_grey0('viridis'))
plt.imshow(np.nansum(smooth_whole*d_sm_wh[27].get_mask(),axis=2),cmap=dc.cmap_grey0('binary'),alpha=0.3)

# plt.xlim(0,500)


plt.imshow(np.nansum(smooth_whole*d_sm_wh[114].get_mask(),axis=0),cmap=dc.cmap_grey0('viridis'))
# plt.xlim(0,500)


d_sm_wh.trunk[10].children[1].children


d_sm_wh.trunk[27]


dor = Dendrogram.compute(smooth_whole[:,70:300,2500:2750],0.18,0.18,60)


np.nanstd(smooth_whole)


len(dor),len(dor.trunk)


len(dor),len(dor.trunk)


len(dor),len(dor.trunk)


for k,kk in enumerate(dor.trunk):
#     if not kk.is_leaf:
    plt.imshow(np.nansum(smooth_whole[:,70:300,2500:2750]*kk.get_mask(),axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'))
    plt.show()
#         print(k,len(kk.sorted_leaves()))


plt.imshow(np.nansum(smooth_whole[:,:,2500:2750],axis=2),cmap=dc.cmap_grey0('viridis'),origin='lower')
# plt.xlim(2500,2750)
plt.xlim(70,300)


plt.imshow(np.nansum(smooth_whole[:,70:300,:],axis=1),cmap=dc.cmap_grey0('viridis'),origin='lower')
plt.xlim(2500,2750)
# plt.ylim(70,300)


plt.imshow(np.nansum(smooth_whole,axis=0),cmap=dc.cmap_grey0('viridis'),origin='lower')
plt.xlim(2500,2750)
plt.ylim(70,300)


for k,kk in enumerate(d_sm_wh_nolims.trunk):
#     if not kk.is_leaf:
    plt.imshow(np.nansum(smooth_whole*kk.get_mask(),axis=0),origin='lower',cmap=dc.cmap_grey0('viridis'))
    plt.show()
#         print(k,len(kk.sorted_leaves()))


len(d_sm_wh.so)


npx = np.zeros(208)
for t,tr in enumerate(d_sm_wh.trunk):
    xx = tr.get_npix()
    npx[t] = xx


plt.hist(npx,bins=500)
plt.yscale('log')
plt.xscale('log')

plt.show()


len(d_sm_wh.trunk[0].children[1].descendants)


import astropy.units as u


(DHT_CO.db*DHT_CO.dl).to(u.steradian).value*(1e3*u.pc)**2


(DHT_CO.dl).to(u.rad).value*(1e3*u.pc)


np.shape(smooth_whole[30])


np.sum(smooth_whole[30]>0.2)


481*2881


86996/1385761


plt.plot(np.nansum(smooth_whole,axis=(1,2)))
plt.hlines(0.2,0,55)


np.nanstd(smooth_whole)


plt.imshow(np.nansum(np.abs(smooth_whole),axis=0),cmap=dc.cmap_grey0('viridis'))
plt.xlim(200,500)
plt.ylim(100,400)


plt.imshow(np.abs(DHT_CO.cube[:,200:400,300:500])[30],cmap=dc.cmap_grey0('viridis'))
plt.colorbar()
# axs[1].pcolormesh(np.transpose(specsmooth)[k])


plt.imshow(np.transpose(specsmooth)[30],cmap=dc.cmap_grey0('viridis'))
plt.colorbar()
# plt.xlim(170,200)
# axs[1].pcolormesh([k])


for k in range(30,40):
    fig,axs=plt.subplots(1,2,figsize=(10,6))
    axs[0].pcolormesh(np.transpose(np.transpose(DHT_CO.cube[:,200:400,300:500]))[k])
    axs[1].pcolormesh(np.transpose(specsmooth)[k])
    plt.show()


np.transpose(specsmooth).shape


plt.pcolormesh(np.transpose(DHT_CO.cube[:,200:400,300:500])[3])


plt.pcolormesh(specsmooth[2])


plt.figure(figsize=(20,5))
lims = [DHT_CO.lons[0],DHT_CO.lons[-1],DHT_CO.lats[0],DHT_CO.lats[-1]]
plt.imshow(np.nansum(np.abs(DHT_CO.cube),axis=0),**kw_imshow,extent=lims)
plt.xlim(150,100)
plt.ylim(0,30)
plt.show()


wo_gal_cent = np.zeros_like(smooth_whole)
wo_gal_cent[:,:,:641] = smooth_whole[:,:,:641]
wo_gal_cent[:,:,2241:] = smooth_whole[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = smooth_whole[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = smooth_whole[:,281:,641:2241]


for x in range(len(wo_gal_cent)):
    plt.figure(figsize=(10,4))
    plt.pcolormesh(wo_gal_cent[x],cmap=dc.cmap_grey0('viridis'))
    plt.show()


# plt.figure(figsize=(10,4))
plt.pcolormesh(wo_gal_cent[35],cmap=dc.cmap_grey0('viridis'))
plt.xlim(2600,2700)
plt.ylim(50,200)
plt.show()


# DHT_CO.lats[201],DHT_CO.lats[281]
# DHT_CO.lons[641],DHT_CO.lons[2241]
wo_gal_cent = np.zeros_like(DHT_CO.cube)
wo_gal_cent[:,:,:641] = DHT_CO.cube[:,:,:641]
wo_gal_cent[:,:,2241:] = DHT_CO.cube[:,:,2241:]
wo_gal_cent[:,:201,641:2241] = DHT_CO.cube[:,:201,641:2241]
wo_gal_cent[:,281:,641:2241] = DHT_CO.cube[:,281:,641:2241]


np.nanstd(wo_gal_cent)


#eq 6 from MD
def weighted_l(cube,dl,lons,dv):
    wco = np.nansum(cube*dv,axis=0) #include data units at somepoint
    wco_tot = np.nansum(wco)
    lfront = np.transpose(cube)
    num = 0
    for l,ll in enumerate(lons):
        this_slice = lfront[l] 
        wcoi = np.nansum(this_slice*dv*ll) #include data units at somepoint
        num += wcoi
        
    lc = num/wco_tot
    return lc


#eq 7 from MD
def weighted_b(cube,db,lats,dv):
    wco = np.nansum(cube*dv,axis=0) #include data units at somepoint
    wco_tot = np.nansum(wco)
    bfront = np.swapaxes(cube,0,1)
    num = 0
    for b,bb in enumerate(lats):
        this_slice = bfront[b] 
        wcoi = np.nansum(this_slice*dv*bb) #include data units at somepoint
        num += wcoi
        
    bc = num/wco_tot
    return bc


#eq 11 from MD
def weighted_v(cube,dv,vels):
    #except instead of the form in 11, it looks more like the forms in 6 and 7 (T B is not from gaussian decomp)
    wco = np.nansum(cube*dv,axis=0) #include data units at somepoint
    wco_tot = np.nansum(wco)
    
    num = 0
    for v,vv in enumerate(vels):
        this_slice = cube[v] 
        wcoi = np.nansum(this_slice*dv*vv) #include data units at somepoint
        num += wcoi
        
    vc = num/wco_tot
    return vc


def weighted_vbl(clob):
    l = weighted_l(clob.cube*clob.units['data'], clob.dl, clob.lons, clob.dv)
    b = weighted_b(clob.cube*clob.units['data'], clob.db, clob.lats, clob.dv)
    v = weighted_v(clob.cube*clob.units['data'], clob.dv, clob.vels)
    
    return v,b,l


#next do eq 11 from MD for <v> ====> DONE
#then go back and fiddle with dendrogram parameters, then feed the trunks into the ^ functions
#then take those l,b,v outputs & put into a sources.inp file to feed into reids distance calculator
#then we go from there


d_wocent_wonpx = Dendrogram.compute(wo_gal_cent,5*0.2,5*0.2,10)


pfar = 0.5


# colnames = ['M ['+Ms.unit.to_string()+']', 'A ['+As.unit.to_string()+']', 'd ['+ds.unit.to_string()+']']
# df = pd.DataFrame({colnames[0]:Ms, colnames[1]:As,colnames[2]:ds})


colnames=['ID','l [deg]','b [deg]','v_lsr [km/s]','pfar']
sources = pd.DataFrame(columns=colnames)


N = len(d_wocent_wonpx.trunk)
IDs = np.zeros(N)
ls = np.zeros(N)
bs = np.zeros(N)
vs = np.zeros(N)
ps = np.full(N, pfar)


sources['pfar'] = ps


# d_wocent_wonpx.save_to('d_wocent_10px.fits') #deprecated keyword


sources


sources['ID'] = IDs
sources['l [deg]'] = ls
sources['b [deg]'] = bs
sources['v_lsr [km/s]'] = vs


sources_with_ds = pd.read_csv('summary.prt', delim_whitespace=True,skiprows=2)


orig_sources = pd.read_csv('sources_info.inp',delim_whitespace=True)


orig_sources['l(deg)'] = round(orig_sources['l(deg)'],2)
orig_sources['b(deg)'] = round(orig_sources['b(deg)'],2)
orig_sources['vlsr(km/s)'] = round(orig_sources['vlsr(km/s)'],1)


orig_sources


sources_with_ds


merged = orig_sources.merge(sources_with_ds,on=['l(deg)','b(deg)','vlsr(km/s)'])


merged['d']


cwdf = merged[['ID','l(deg)','b(deg)','vlsr(km/s)','d1(kpc)','d1_prob','d1_unc','d2(kpc)','d2_prob','d2_unc']]


d100_sig = Dendrogram.compute(wo_gal_cent,0.2,0.2,100)


len(d100_sig.trunk)


for i in range(153):
    tr = d100_sig.trunk[int(i)]
    m = tr.get_mask()
    cb = m*DHT_CO.cube
    clob = DHT_CO.clone()
    clob.cube = cb
    
    
    lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
    plt.figure(figsize=(15,10))
    if tr.get_npix() > 300:
        plt.imshow(np.nansum(cb*clob.dv.value,axis=0),**kw_imshow)#,extent=lims)
        v,b,l = weighted_vbl(clob)
    #     plt.xlim(l-5,l+5)
    #     plt.ylim(b-5,b+5)
        plt.title(str(i)+'th trunk')
plt.show()


for i in (d100_sig.trunk):
    if i.get_npix() > 300:
        print(i.get_npix())


for i in range(153):
    tr = d100_sig.trunk[int(i)]
    m = tr.get_mask()
    cb = m*DHT_CO.cube
    clob = DHT_CO.clone()
    clob.cube = cb
    
    
    lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
    
    plt.imshow(np.nansum(cb*clob.dv.value,axis=0),**kw_imshow,extent=lims)
    v,b,l = weighted_vbl(clob)
    plt.xlim(l-5,l+5)
    plt.ylim(b-5,b+5)
    plt.title(str(i)+'th trunk')
    plt.show()


for i in cwdf['ID']:
    tr = d_wocent_wonpx.trunk[int(i)]
    m = tr.get_mask()
    cb = m*DHT_CO.cube
    clob = DHT_CO.clone()
    clob.cube = cb
    
    
#     lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
    
    plt.imshow(np.nansum(cb*clob.dv.value,axis=0),**kw_imshow)#,extent=lims)
    v,b,l = weighted_vbl(clob)
#     plt.xlim(l-5,l+5)
#     plt.ylim(b-5,b+5)
#     plt.title(str(i)+'th trunk')
plt.show()


for i in cwdf['ID']:
    tr = d_wocent_wonpx.trunk[int(i)]
    m = tr.get_mask()
    cb = m*DHT_CO.cube
    clob = DHT_CO.clone()
    clob.cube = cb
    
    
    lims = [clob.lons[0],clob.lons[-1],clob.lats[0],clob.lats[-1]]
    
    plt.imshow(np.nansum(cb*clob.dv.value,axis=0),**kw_imshow,extent=lims)
    v,b,l = weighted_vbl(clob)
    plt.xlim(l-5,l+5)
    plt.ylim(b-5,b+5)
    plt.title(str(i)+'th trunk')
    plt.show()


plt.scatter(cwdf['l(deg)'],cwdf['b(deg)'],marker='+')
# plt.scatter(orig_sources['l(deg)'],orig_sources['b(deg)'],marker='x')


np.argmin(np.abs(cwdf['l(deg)']-211.))


cwdf['l(deg)'][134],cwdf['b(deg)'][134],cwdf['ID'][134]


cat
