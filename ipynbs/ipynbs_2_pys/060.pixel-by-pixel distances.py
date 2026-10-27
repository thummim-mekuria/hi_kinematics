#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from astropy.table import QTable


import datacube as dc


from astrodendro import Dendrogram


co_dend = Dendrogram.load_from('COGAL_dendrogram.fits')


oriona = co_dend[246]


perseus = co_dend[214]


ophiuchus = co_dend[260]


corona = co_dend[275]


plt.imshow(np.nansum(oriona.get_mask(),axis=0),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(2500,2800)
plt.ylim(50,200)


plt.imshow(np.nansum(oriona.get_mask(),axis=1),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(2500,2800)
# plt.ylim(50,200)


plt.imshow(np.nansum(perseus.get_mask(),axis=0),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(100,300)
plt.ylim(30,150)


plt.imshow(np.nansum(perseus.get_mask(),axis=1),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(100,300)
# plt.ylim(30,150)


plt.imshow(np.nansum(ophiuchus.get_mask(),axis=0),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(1300,1700)
plt.ylim(300,450)


plt.imshow(np.nansum(ophiuchus.get_mask(),axis=1),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(1300,1700)
# plt.ylim(300,450)


plt.imshow(np.nansum(corona.get_mask(),axis=0),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(1300,1500)
plt.ylim(0,150)


plt.imshow(np.nansum(corona.get_mask(),axis=1),\
           origin='lower',cmap=dc.cmap_min0('viridis'))
plt.xlim(1300,1500)
# plt.ylim(0,150)


cat = pd.read_csv('../not_rn_3/physcl_prprts_catalog.csv')


np.where(cat['_idx']==246)[0][0]


ora_info = cat.iloc[120]


per_info = cat.iloc[100]


oph_info =  cat.iloc[131]


cra_info = cat.iloc[145]


ora_info


per_info


oph_info


cra_info


ora_info['area_sky<deg2>']/(0.125**2)


per_info['area_sky<deg2>']/(0.125**2)


oph_info['area_sky<deg2>']/(0.125**2)


cra_info['area_sky<deg2>']/(0.125**2)


len(oriona.indices()[0]),len(corona.indices()[0]),\
len(ophiuchus.indices()[0]),len(perseus.indices()[0])


#values from dht header
def v_kms(arrays):
    output = []
    for ar in arrays:
        this = ar*1.3-49.4
        output.append(this)
    return output

def b_deg(arrays):
    output = []
    for ar in arrays:
        this = ar*0.125-30.125
        output.append(this)
    return output

def l_deg(arrays):
    output = []
    for ar in arrays:
        this = ar*(-0.125)+180.125
        output.append(this)
    return output


ora_v,ora_b,ora_l = np.array(oriona.indices())


per_v,per_b,per_l = np.array(perseus.indices())


oph_v,oph_b,oph_l = np.array(ophiuchus.indices())


cra_v,cra_b,cra_l = np.array(corona.indices())


vs2 = v_kms(np.concatenate((ora_v[::2],per_v[::2],oph_v[::2],cra_v[::2])))
bs2 = b_deg(np.concatenate((ora_b[::2],per_b[::2],oph_b[::2],cra_b[::2])))
ls2 = l_deg(np.concatenate((ora_l[::2],per_l[::2],oph_l[::2],cra_l[::2])))

names2 = np.concatenate((np.array(['oriona']*len(ora_l[::2])),\
                         np.array(['perseus']*len(per_l[::2])),\
                         np.array(['ophiuchus']*len(oph_l[::2])),\
                         np.array(['corona']*len(cra_l[::2]))))

sources2 = pd.DataFrame(np.array((names2,ls2,bs2,vs2)).transpose(),\
             columns=['name','lon','lat','vel'])

sources2['v_unc'] = [0.6]*len(sources2)
sources2['pfar'] = [0.5]*len(sources2)

sources.to_csv('sources_info2.inp',index=False)


vs = v_kms(np.concatenate((ora_v[::5],per_v[::5],oph_v[::5],cra_v[::5])))
bs = b_deg(np.concatenate((ora_b[::5],per_b[::5],oph_b[::5],cra_b[::5])))
ls = l_deg(np.concatenate((ora_l[::5],per_l[::5],oph_l[::5],cra_l[::5])))


names = np.concatenate((np.array(['oriona']*971),np.array(['perseus']*539),\
                np.array(['ophiuchus']*711),np.array(['corona']*646)))


sources = pd.DataFrame(np.array((names,ls,bs,vs)).transpose(),\
             columns=['name','lon','lat','vel'])

sources2['v_unc'] = [0.6]*len(sources2)
sources2['pfar'] = [0.5]*len(sources2)


sources['v_unc'] = [0.6]*2867


sources['pfar'] = [0.5]*2867


sources.to_csv('sources_info.inp',index=False)


sources2.to_csv('sources_info.inp',index=False)


out2 = pd.read_csv('summary.prt',skiprows=4,header=None,delim_whitespace=True)


len(out2)


out = pd.read_csv('summary.prt',skiprows=4,header=None,delim_whitespace=True)


heads = list(pd.read_csv('summary.prt',skiprows=2).columns)[0].split()[1:]


heads


nice_heads = ['Long.',
 'Lat.',
 'Vlsr',
 'Dist1',
 'Prob1',
 'Arm1',
 'Dist2',
 'Prob2',
 'Arm2']


nice_out = out[[0,1,2,4,6,7,8,10,11]]


nice_out['name'] = sources['name']


nice_out = nice_out[['name','Long.', 'Lat.', 'Vlsr', 'Dist1',\
                     'Prob1', 'Arm1', 'Dist2', 'Prob2', 'Arm2']]


nice_out2 = out2[[0,1,2,4,6,7,8,10,11]]
nice_out2.columns = nice_heads
nice_out2['name'] = sources2['name']

nice_out2 = nice_out2[['name','Long.', 'Lat.', 'Vlsr', 'Dist1',\
                     'Prob1', 'Arm1', 'Dist2', 'Prob2', 'Arm2']]


ora_df = nice_out.iloc[range(0,970)]


per_df = nice_out.iloc[range(971,1510)] 


oph_df = nice_out.iloc[range(1510,2221)]


cra_df = nice_out.iloc[range(2221,2867)]


[i for i in range(len(nice_out2)-1) if nice_out2['name'][i]!=nice_out2['name'][i+1]]


ora_df2 = nice_out2.iloc[range(2427)]
per_df2 = nice_out2.iloc[range(2427,3773)]
oph_df2 = nice_out2.iloc[range(3773,5549)]
cra_df2 = nice_out2.iloc[range(5549,7163)]


#look at the pdf for the pixels across the entire cloud,
#add up the pdfs (or average?) and see


dht = fits.open('DHT_CO.fits')[0]


wco = np.nansum(dht.data*1.3,axis=0)


lws_msk = wco > 1.5


fig,(ax1,ax2,ax3) =plt.subplots(3,1,figsize=(7,10),layout='constrained')

ax1.imshow(wco*lws_msk,origin='lower',cmap=dc.cmap_min0('binary'),extent=[180,-180,-30,30],norm='log')
ax1.set_xlim(-140,-158)
ax1.set_ylim(-22,-15)
ax1.set_xticks([])
ax1.set_ylabel('b [deg]')

ora_ds = ax2.scatter(ora_df2['Long.'],ora_df2['Lat.'],c=ora_df2['Dist1'],alpha=0.7,s=15)
fig.colorbar(ora_ds,ax=ax2,orientation='horizontal')
ax2.set_xlim(-140,-158)
ax2.set_ylim(-22,-15)
ax2.set_xticks([])
ax2.set_ylabel('b [deg]')

ora_ds_2 = ax3.scatter(ora_df2['Long.'],ora_df2['Lat.'],c=ora_df2['Dist2'],alpha=0.7,s=15)
fig.colorbar(ora_ds_2,ax=ax3,orientation='horizontal')
ax3.set_xlim(-140,-158)
ax3.set_ylim(-22,-15)
ax3.set_xlabel('l [deg]')
ax3.set_ylabel('b [deg]')

fig.suptitle('orion a')


fig,(ax1,ax2,ax3) =plt.subplots(1,3,figsize=(10,5),layout='constrained')

ax1.imshow(wco*lws_msk,origin='lower',cmap=dc.cmap_min0('binary'),extent=[180,-180,-30,30],norm='log')
ax1.set_xlim(162,154)
ax1.set_ylim(-26,-15)
ax1.set_xlabel('l [deg]')
ax1.set_ylabel('b [deg]')

per_ds = ax2.scatter(per_df2['Long.'],per_df2['Lat.'],c=per_df2['Dist1'],alpha=0.7,s=15)

fig.colorbar(per_ds,ax=ax2,orientation='horizontal')
ax2.set_xlim(162,154)
ax2.set_ylim(-26,-15)
ax2.set_xlabel('b [deg]')
ax2.set_yticks([])

per_ds_2 = ax3.scatter(per_df2['Long.'],per_df2['Lat.'],c=per_df2['Dist2'],alpha=0.7,s=15)

fig.colorbar(per_ds_2,ax=ax3,orientation='horizontal')
ax3.set_xlim(162,154)
ax3.set_ylim(-26,-15)
ax3.set_xlabel('l [deg]')
ax3.set_yticks([])

fig.suptitle('perseus')


fig,(ax1,ax2,ax3) =plt.subplots(3,1,figsize=(6,12),layout='constrained')

ax1.imshow(wco*lws_msk,origin='lower',cmap=dc.cmap_min0('binary'),\
           extent=[180,-180,-30,30],norm='log')
ax1.set_xlim(2,-12)
ax1.set_ylim(12,20)
ax1.set_xticks([])
ax1.set_ylabel('b [deg]')

oph_ds = ax2.scatter(oph_df2['Long.'],oph_df2['Lat.'],\
                     c=oph_df2['Dist1'],alpha=0.7,s=15)

fig.colorbar(oph_ds,ax=ax2,orientation='horizontal')
ax2.set_xlim(2,-12)
ax2.set_ylim(12,20)
ax2.set_xticks([])
ax2.set_ylabel('b [deg]')


oph_ds_2 = ax3.scatter(oph_df2['Long.'],oph_df2['Lat.'],\
                       c=oph_df2['Dist2'],alpha=0.7,s=15)

fig.colorbar(oph_ds_2,ax=ax3,orientation='horizontal')
ax3.set_xlim(2,-12)
ax3.set_ylim(12,20)
ax3.set_xlabel('l [deg]')
ax3.set_ylabel('b [deg]')

fig.suptitle('ophiuchus')
# ax2.set_title('dist1')
# ax3.set_title('dist2')


# fig.tight_layout()


fig,(ax1,ax2,ax3) =plt.subplots(1,3,figsize=(11,5),layout='constrained')

ax1.imshow(wco*lws_msk,origin='lower',cmap=dc.cmap_min0('binary'),extent=[180,-180,-30,30],norm='log')

ax1.set_xlim(4,-2)
ax1.set_ylim(-24,-16)
ax1.set_ylabel('b [deg]')
ax1.set_xlabel('l [deg]')

cra_ds = ax2.scatter(cra_df2['Long.'],cra_df2['Lat.'],c=cra_df2['Dist1'],alpha=0.7,s=15)

fig.colorbar(cra_ds,ax=ax2,orientation='horizontal')
ax2.set_xlim(4,-2)
ax2.set_ylim(-24,-16)
ax2.set_yticks([])
ax2.set_xlabel('l [deg]')

cra_ds_2 = ax3.scatter(cra_df2['Long.'],cra_df2['Lat.'],c=cra_df2['Dist2'],alpha=0.7,s=15)

fig.colorbar(cra_ds_2,ax=ax3,orientation='horizontal')
ax3.set_xlim(4,-2)
ax3.set_ylim(-24,-16)
ax3.set_yticks([])
ax3.set_xlabel('l [deg]')

fig.suptitle('corona australis')


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


cat.iloc[70]


len(np.where(df1['Arm1'] != 'AqR')[0])


len(np.where(df1['Arm1'] == 'AqR')[0])


1*u.deg.to(u.rad)*1e3*u.pc


181*51*16


481*0.125 -30.125


lons = np.arange(180,-181,-2)
lats = np.arange(-30,31,4)
vels = np.arange(-50,51,2)
names = np.arange(135000)

grid = np.array(np.meshgrid(lons,lats,vels))


grid.transpose().reshape(-1,3).transpose()


out = pd.read_csv('summary.prt',skiprows=4,delim_whitespace=True,header=None)
nice_out = out[[0,1,2,4,6,7,8,10,11]]
nice_heads = ['Long.','Lat.','Vlsr','Dist1','Prob1',
'Arm1','Dist2','Prob2','Arm2']
nice_out.columns = nice_heads


nice_out.to_csv('whole_gal_dists.csv',index=False)


plt.figure(figsize=(20,7))
plt.scatter(nice_out['Long.'],nice_out['Lat.'],c=1e3*nice_out['Dist1'],norm='log',s=100)
plt.xlim(180,-180)
plt.colorbar()


galvels = np.arange(76)*1.3 - 49.4


gallons = np.arange(2881)*(-0.125) + 180.125


w_lv = np.nansum(dht.data,axis=1)


plt.figure(figsize=(20,7))
plt.imshow(w_lv,extent=[180,-180,-50,50],origin='lower',vmin=-10,cmap=dc.cmap_min0('viridis'))
plt.contour(gallons,galvels,w_lv,levels=[1],colors=['red'])
# plt.colorbar()


plt.figure(figsize=(20,7))
plt.imshow(np.nansum(dht.data,axis=2),extent=[-30,30,-50,50],\
           origin='lower',vmin=0,cmap=dc.cmap_min0('viridis'))
# plt.contour(gallons,galvels,w_lv,levels=[1],colors=['red'])
# plt.colorbar()


plt.figure(figsize=(20,7))
plt.scatter(nice_out['Long.'],nice_out['Lat.'],c=1e3*nice_out['Dist1'],norm='log')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

plt.xlim(180,-180)


plt.figure(figsize=(7,7))
plt.scatter(nice_out['Vlsr'][::-1],nice_out['Dist1'][::-1],c=nice_out['Long.'],s=5,cmap='PiYG_r')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

plt.ylim(-2,35)
plt.xlim(-50,50)


plt.figure(figsize=(7,7))
plt.scatter(nice_out['Vlsr'],nice_out['Dist1'],c=nice_out['Long.'],s=5,cmap='PiYG_r')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

plt.ylim(-2,35)
plt.xlim(-50,50)


plt.figure(figsize=(20,7))
plt.scatter(nice_out['Long.'],nice_out['Dist1'],c=nice_out['Vlsr'],s=15,cmap='PiYG')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

plt.xlim(180,-180)
plt.ylim(0,35)


plt.figure(figsize=(20,7))
plt.scatter(nice_out['Long.'],nice_out['Dist2'],c=nice_out['Vlsr'],s=15,cmap='PiYG')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

plt.xlim(180,-180)
plt.ylim(-2,35)


plt.figure(figsize=(20,7))
plt.scatter(nice_out['Long.'],nice_out['Vlsr'],c=1e3*nice_out['Dist1'],norm='log')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

plt.xlim(180,-180)


plt.figure(figsize=(20,7))
plt.scatter(nice_out['Vlsr'],nice_out['Lat.'],c=1e3*nice_out['Dist1'],norm='log')
plt.colorbar()
# plt.contour(gallons,galvels,w_lv,levels=[1,10,100],colors=['red'],linewidths=[0.4,0.7,1])

# plt.xlim(180,-180)


plt.scatter(nice_out['Long.'],nice_out['Dist1'])


x = nice_out['Dist1']*np.sin(nice_out['Long.']*u.deg.to(u.rad))


y = nice_out['Dist1']*np.cos(nice_out['Long.']*u.deg.to(u.rad))-8


plt.figure(figsize=(10,10))
plt.scatter(x,y+8,c=nice_out['Vlsr'],s=10)
plt.scatter(0,0,marker='x',color='red',s=20)
plt.scatter(0,8,marker='x',color='red',s=20)

plt.xlim(1,-1)
plt.ylim(-1,1)


plt.figure(figsize=(10,10))
plt.scatter(x,y+8,c=nice_out['Vlsr'],s=10)
plt.scatter(0,0,marker='x',color='red',s=20)
plt.scatter(0,8,marker='x',color='red',s=20)

plt.xlim(15,-15)
plt.ylim(-10,25)


dht.header


#okay even more pixel-by-pixel analyses
df1 = pd.read_csv('idx167.csv')
df2 = pd.read_csv('idx198.csv')
df3 = pd.read_csv('idx217.csv')
df4 = pd.read_csv('idx238.csv')
df5 = pd.read_csv('idx248.csv')
df6 = pd.read_csv('idx249.csv')
df7 = pd.read_csv('idx271.csv')
df8 = pd.read_csv('idx281.csv')


plt.figure(figsize=(15,7))
for df in [df1]:#,df2,df3,df4,df5,df6,df7,df8]:
    plt.scatter(df['Long.'],df['Lat.'],s=5,alpha=0.7,c=1e3*df['Dist1'],vmax=500)

plt.colorbar(label='d [pc]')
plt.xlim(38,11)
plt.ylim(0,20)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.show()


plt.figure(figsize=(15,7))
for df in [df1]:#,df2,df3,df4,df5,df6,df7,df8]:
    plt.scatter(df['Long.'],df['Vlsr'],s=15,alpha=0.7,c=1e3*df['Dist1'],vmax=500)

plt.colorbar(label='d [pc]')
plt.xlim(38,11)
plt.ylim(-10,25)
plt.xlabel('l [deg]')
plt.ylabel('v [km/s]')
plt.show()


plt.figure(figsize=(15,8))
for df in [df1,df2,df3,df4,df5,df6,df7,df8]:
    plt.scatter(df['Long.'],df['Lat.'],s=5,c=df['Dist2'])

plt.colorbar()
plt.xlim(45,20)
plt.ylim(-25,-10)
plt.show()


plt.figure(figsize=(15,8))
for df in [df1,df2,df3,df4,df5,df6,df7,df8]:
    plt.scatter(df['Long.'],df['Lat.'],s=5,c=df['Dist2'])

plt.colorbar()
plt.xlim(-7,-30)
plt.ylim(0,25)
plt.show()


plt.figure(figsize=(15,8))
for df in [df1,df2,df3,df4,df5,df6,df7,df8]:
    plt.scatter(df['Long.'],df['Lat.'],s=0.05,c=df['Dist2'])

plt.colorbar()
plt.xlim(50,-50)
plt.show()


dend[167].indices()


dend = Dendrogram.load_from('COGAL_dendrogram.fits')


kw = {'origin':'lower','extent':[180,-180,-30,30]}


for idx in [167,198,217,238,248,249,271,281]:
    plt.imshow(np.nansum(dend[idx].get_mask(),axis=0),cmap=dc.cmap_min0('viridis'),**kw)

plt.xlim(50,-50)
# plt.ylim(250,350)


galdss = pd.read_csv('summary.prt',skiprows=4,delim_whitespace=True,header=None)[[0,1,2,4,5,6,7,8,9,10,11]]


galdss.columns = ['lon<deg>', 'lat<deg>', 'vel<km/s>','d1<kpc>', 'd1_unc<kpc>', 'd1_prob',
       'arm1', 'd2<kpc>', 'd2_unc<kpc>', 'd2_prob', 'arm2']


galdss.to_csv('whole_gal_dists.csv',index=False)


galds = pd.read_csv('whole_gal_dists.csv')


zcat = QTable(fits.open('final_cat.fits')[1].data)


galds


zls = np.nanmean(zcat['lon<deg>'],axis=1)
zbs = np.nanmean(zcat['lat<deg>'],axis=1)


zvsish = np.argmax(zcat['vspectrum_GMC<K>'],axis=1)*1.3 - 100.1


rmins = np.zeros(19,dtype='int')


nearest_ds = np.zeros(19)
nearest_d_unc = np.zeros(19)


for i in range(19):
    r = np.sqrt((zbs[i] - galds['lat<deg>'])**2 +\
                 (zls[i] - galds['lon<deg>'])**2 + \
                 (zvsish[i]-galds['vel<km/s>'])**2)
    
    la = np.nanstd(zcat['lon<deg>'],axis=1)[i]
    ba = np.nanstd(zcat['lat<deg>'],axis=1)[i]
    va = zcat['vdispersion_GMC<km/s>'][i]
    
    nearest_ds[i] = np.nanmean(galds[r < max(la,ba,va)]['d1<kpc>'])
    nearest_d_unc[i] = np.nanmean(galds[r < max(la,ba,va)]['d1_unc<kpc>'])


plt.errorbar(zcat['names'],nearest_ds* 1e3,\
             yerr=1e3*nearest_d_unc,marker='x',label='reid_avg',\
            elinewidth=0.35)

plt.plot(zcat['names'],np.array(zcat['distance<pc>']),marker='x',label='zucker')


plt.xticks(rotation='vertical')
plt.ylim(0,2100)
plt.legend()
plt.show()


deld =  (nearest_ds* 1e3 - np.array(zcat['distance<pc>']))/np.array(zcat['distance<pc>'])


plt.plot(zcat['names'],np.abs(deld),marker='x')

plt.xticks(rotation='vertical')
plt.show()


181*51*16


del l = 2 del b = 4 del v = 2


ds['lon<deg>'].to_numpy().reshape(53,16,-1).shape


rawr=pd.read_csv('sources_info.inp',skiprows=2,header=None)[[1,2,3]]


rawr.columns = ['lon<deg>', 'lat<deg>', 'vel<km/s>']


otw = rawr.merge(ds,on=['lon<deg>', 'lat<deg>', 'vel<km/s>'],how='outer')


otw['lon<deg>'].to_numpy().reshape(51,-1,16)[0]#.swapaxes(0,1).transpose().shape


d1s = otw['d1<kpc>'].to_numpy().reshape(51,-1,16).swapaxes(1,2)
d1uncs = otw['d1_unc<kpc>'].to_numpy().reshape(51,-1,16).swapaxes(1,2)
d1probs = otw['d1_prob'].to_numpy().reshape(51,-1,16).swapaxes(1,2)


hd1 = fits.HDUList([fits.PrimaryHDU(d1s),\
              fits.ImageHDU(d1uncs),\
              fits.ImageHDU(d1probs)])


d2s = otw['d2<kpc>'].to_numpy().reshape(51,-1,16).swapaxes(1,2)
d2uncs = otw['d2_unc<kpc>'].to_numpy().reshape(51,-1,16).swapaxes(1,2)
d2probs = otw['d2_prob'].to_numpy().reshape(51,-1,16).swapaxes(1,2)


hd2 = fits.HDUList([fits.PrimaryHDU(d2s),\
              fits.ImageHDU(d2uncs),\
              fits.ImageHDU(d2probs)])


hd1.writeto('reid_d1s_gal.fits',overwrite=True)
hd2.writeto('reid_d2s_gal.fits',overwrite=True)


ds = pd.read_csv('whole_gal_dists.csv')


((ds['lon<deg>'] - 180)/-2).astype('int')


((ds['lat<deg>'] + 30)/4).astype('int')


((ds['vel<km/s>'] + 50)/2).astype('int')
