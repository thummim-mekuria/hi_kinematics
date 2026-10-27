#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle,Ellipse

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import models,fitting

import fun as fns
import datacube as dc

his = fits.open('hi_clouds.fits')
cos = fits.open('h2_clouds.fits')
# dust = fits.open('mean_and_std_lbd.fits')[1].data #for later


pcat = QTable(fits.open('post_flash_cat.fits')[1].data)


vels = [np.arange(pcat['vo'][i],pcat['vi'][i])*1.3-100 for i in range(22)]


pxpc = 0.125*u.deg.to(u.rad)*np.array(pcat['distance<pc>'])*u.pc


cmom0s = [np.nansum(cos[i].data,axis=0)*1.3*1.97e20 for i in range(22)]


cmom1s = [fns.get_mom1(cos[i].data,vels[i]) for i in range(22)]


cmom2s = [np.sqrt(fns.get_momN(cos[i].data,vels[i],2)) for i in range(22)]


hmom0s = [np.nansum(his[i].data,axis=0)*1.3*1.82e18 for i in range(22)]
hmom1s = [fns.get_mom1(his[i].data,vels[i]) for i in range(22)]
hmom2s = [fns.get_momN(his[i].data,vels[i],2) for i in range(22)]


coords = [np.array(np.meshgrid(np.arange(pcat['loc'][i],pcat['lic'][i])*-0.125+360,
             np.arange(pcat['boc'][i],pcat['bic'][i])*0.125-50)) for i in range(22)]


#surface density of the gmc -> can be inferred from the ratio  of aread
#between the hi envelope & gmc
#higher ratio->compacter gmc->higher sigma
for i in range(22):
    hlb = np.nansum(his[i].data,axis=0)*1.3*1.82e18
    clb = np.nansum(cos[i].data,axis=0)*1.3*1.97e20
    
    
    print(pcat['name'][i],'\n',
          round(len(np.ravel(hlb))/len(clb[clb>3e20]),2),
          round(pcat['Sigma_GMC<Msun/pc2>'][i],2))
    print('\n')
#     plt.scatter(np.ravel(hlb),np.ravel(clb),s=1e-2)
#     plt.show()
    

from sklearn.decomposition import PCA


pca2 = PCA(n_components=2)


#something is going on at velocity dispersion ~ 0.6 for the moleculae clouds


#calculating pca across the diff types of maps
#for all of the clouds, 2 components can explain >80% of the variance
#and for 15 of them, >90% of the variance
#not sure what this means though or how to take that
#and turn it into a clustering along the diff maps
for i in range(22):
    surden = (np.array(cmom0s[i])*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value
    d0dy,d0dx = np.gradient(surden)
    d0 = np.sqrt(d0dy**2 + d0dx**2)

    d1dy,d1dx = np.gradient(cmom1s[i])
    d1 = np.sqrt(d1dy**2 + d1dx**2)

    d2dy,d2dx = np.gradient(cmom2s[i])
    d2 = np.sqrt(d2dy**2 + d2dx**2)
    
    
    df = pd.DataFrame(np.array([np.ravel(coords[i][0]),np.ravel(coords[i][1]),
                                np.ravel(surden),
    np.ravel(cmom1s[i]),np.ravel(cmom2s[i]),
    np.ravel(d0),np.ravel(d1),np.ravel(d2)]).transpose(),
    columns=['l','b','w','v','s','dw','dv','ds'])

    df = df.dropna()
    df = df[df['w']>3]
    
    clpca = pca2.fit_transform(df)
    
    print(pcat['name'][i],'\n',pca2.explained_variance_ratio_,
          round(np.nansum(pca2.explained_variance_ratio_),2),'\n')


#printing out covariant coefficients across the diff types of maps
#for most the strongest correltion is between w & s
#which ig makes sense bc the wider the spectrum, the more mass when you add it all up
#but also, for the clouds that have strong velocity gradients, (eg california, perseus, l-orionis)
#it shows up as a strong correlation between l/b and v
#also strong contenders: dv&ds w&dw

for i in range(22):
    surden = (np.array(cmom0s[i])*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value
    d0dy,d0dx = np.gradient(surden)
    d0 = np.sqrt(d0dy**2 + d0dx**2)

    d1dy,d1dx = np.gradient(cmom1s[i])
    d1 = np.sqrt(d1dy**2 + d1dx**2)

    d2dy,d2dx = np.gradient(cmom2s[i])
    d2 = np.sqrt(d2dy**2 + d2dx**2)
    
    
    df = pd.DataFrame(np.array([np.ravel(coords[i][0]),np.ravel(coords[i][1]),
                                np.ravel(surden),
    np.ravel(cmom1s[i]),np.ravel(cmom2s[i]),
    np.ravel(d0),np.ravel(d1),np.ravel(d2)]).transpose(),
    columns=['l','b','w','v','s','dw','dv','ds'])

    df = df.dropna()
    df = df[df['w']>3]

    print(pcat['name'][i],)
    mtx = np.around(np.corrcoef(np.array([i[1] for i in df.items()])),2)
    print(mtx)
    print(np.nanmax(np.abs(mtx-np.identity(8))),np.argmax(np.abs(mtx-np.identity(8))))
    print('\n\n')
    

2.7e6*u.lightyear.to(u.pc)/1e6


from copy import deepcopy


pdata = deepcopy(pcat)


pdata['lon_GMC<deg>'] = np.nanmean(pcat['lon_GMC<deg>'],axis=1)


pdata['lat_GMC<deg>'] = np.nanmean(pcat['lat_GMC<deg>'],axis=1)


pdata['vel<km/s>'] = np.nanmean(pcat['vel<km/s>'],axis=1)


pdata.remove_columns(['name','loh','lih','boh','bih','loc','lic','bic','boc','vo','vi'])


pdata.columns


plt.scatter(pcat['vel_grad_GMC<km/s/pc>'],pcat['vel_grad_HI<km/s/pc>'])
# plt.scatter(,pcat['mass_HI<Msun>'])

plt.loglog()
plt.show()


plt.scatter(pcat['spec_ang_mom_GMC<cm2/s>'],pcat['spec_ang_mom_HI<cm2/s>'])
# plt.scatter(,pcat['mass_HI<Msun>'])

plt.loglog()
plt.show()


plt.scatter(pcat['spec_ang_mom_HI<cm2/s>'],pcat['mass_GMC<Msun>'])
plt.scatter(pcat['spec_ang_mom_GMC<cm2/s>'],pcat['mass_HI<Msun>'])

plt.loglog()
plt.show()


plt.scatter(pcat['spec_ang_mom_GMC<cm2/s>'],pcat['mass_GMC<Msun>'])
plt.scatter(pcat['spec_ang_mom_HI<cm2/s>'],pcat['mass_HI<Msun>'])

plt.loglog()
plt.show()


pd.DataFrame(np.around(np.corrcoef(np.array([i[1] for i in pdata.items()])),2),columns =pcat.columns[1:26]).to_csv('metrix.csv')


from IPython.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))


plt.scatter(pdata['lon_GMC<deg>'],pdata['vel_grad_HI<km/s/pc>'])
# plt.scatter(pdata['lon_GMC<deg>'],pdata['vel_grad_GMC<km/s/pc>'])

plt.xlim(180,-180)


allpca = pca2.fit_transform(pdata.to_pandas())


pcadf = pd.DataFrame(allpca)


#gemini is the outlier along the first principal component -> prob more of a cloud complex than a single cloud
#get rid of (i think it was also the furthest away)
#i don't have to throw it away, can just say, had at first selected it, then saw that it had these specific
#idiosyncracies, probably means it's not the same Type of object as the other things that we care about

#narrate choosing the catalog as: here are zuckers clouds, and then here are the iterative cuts we are doing in order
#to create a refined sample


#if we zoom in other than gemini, california, then monr2, rosette, maddalena, cma ob1, orion a and b
plt.scatter(pcat['name'],pcadf[0]/1e64)
plt.xticks(rotation='vertical')

plt.ylim(-11,-5)
plt.show()


#second component
#maddalena stands out as the clear(?) outslier
#then rosette
#then it's a little more blurred
plt.scatter(pcat['name'],pcadf[1]/1e64)
plt.xticks(rotation='vertical')

plt.ylim(-0.2,0.4)
plt.show()


plt.figure(figsize=(10,10))
plt.scatter(pcadf[0]/1e64,pcadf[1]/1e64)
for i in range(22):
    plt.annotate(pcat['name'][i],(pcadf[0][i]/1e64,pcadf[1][i]/1e64))
    
    
plt.xlim(-11,-6)
plt.ylim(-0.2,0.1)


plt.figure(figsize=(10,10))
plt.scatter(pcadf[0]/1e64,pcadf[1]/1e64)
for i in range(22):
    plt.annotate(pcat['name'][i],(pcadf[0][i]/1e64,pcadf[1][i]/1e64))
    
    
plt.xlim(-10.5,-9)
plt.ylim(-0.16,-0.05)


plt.figure(figsize=(10,10))
plt.scatter(pcadf[0]/1e64,pcadf[1]/1e64)
for i in range(22):
    plt.annotate(pcat['name'][i],(pcadf[0][i]/1e64,pcadf[1][i]/1e64))
    
    
plt.xlim(-10.2,-10)
plt.ylim(-0.13,-0.12)


#hmm i wonder if running pca without gemini will still give the same results
#also find the thing that tells you about which columns contributed most to defining the principal component
#(or alternativel if it is an even mix of all the columns)


scores_pc1 = pd.DataFrame(np.array([pdata.colnames[:25],np.abs(pca2.components_[0])]).transpose())


#the orders of magnitude on the numbers are a=making a big difference (look at scores pc1 -> the ang mom (erg.s) column dominates just cos its in units
#that make it go to ~10^60), so change units (and mayb ealso put stuff in log and re-run pca
pdata['ang_mom_GMC<Msun.pc2/yr>'] = pdata['ang_mom_GMC<erg.s>']*(u.erg*u.s).to(u.solMass*u.pc**2/u.yr)
pdata['ang_mom_HI<Msun.pc2/yr>'] = pdata['ang_mom_HI<erg.s>']*(u.erg*u.s).to(u.solMass*u.pc**2/u.yr)


pdata.remove_columns(['ang_mom_GMC<erg.s>','ang_mom_HI<erg.s>'])


pdata['spec_ang_mom_GMC<pc.km/s>'] = (pdata['spec_ang_mom_GMC<cm2/s>']*u.cm**2/u.s).to(u.pc*u.km/u.s)
pdata['spec_ang_mom_HI<pc.km/s>'] = (pdata['spec_ang_mom_HI<cm2/s>']*u.cm**2/u.s).to(u.pc*u.km/u.s)


pdata.remove_columns(['spec_ang_mom_GMC<cm2/s>','spec_ang_mom_HI<cm2/s>'])


pdata.remove_columns(['log_mass_GMC<Msun>','log_mass_HI<Msun>'])


pdata['drot_ax<deg>'] = np.abs(pdata['drot_ax<deg>'])


pdata.remove_columns(['rot_ax_GMC<deg>','rot_ax_HI<deg>'])


logdata = np.log10(pdata[pdata.colnames[4:]].to_pandas())


pca3 = PCA(n_components=3)


allpca = pca3.fit_transform(logdata)


pcadf = pd.DataFrame(allpca)


np.nanmean(pcadf),np.nanstd(pcadf)


#nvm doing it w the log-ified & rescaled values shows taht yes gemini is an outlier, but not so far away

fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.scatter(pcadf[0],pcadf[1])
for i in range(22):
    ax.annotate(pcat['name'][i],(pcadf[0][i],pcadf[1][i]))
    
 
ax.add_patch(Circle((0,0),1.7,fill=False,color='red',linewidth=0.7,linestyle='--'))
plt.xlim(-10,10)
plt.ylim(-5,5)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.scatter(pcadf[1],pcadf[2])
for i in range(22):
    ax.annotate(pcat['name'][i],(pcadf[1][i],pcadf[2][i]))
    
 
ax.add_patch(Circle((0,0),1.7,fill=False,color='red',linewidth=0.7,linestyle='--'))
# plt.xlim(-10,10)
# plt.ylim(-5,5)


fig,ax = plt.subplots(1,1,figsize=(20,10))
ax.scatter(pcadf[0],pcadf[2])
for i in range(22):
    ax.annotate(pcat['name'][i],(pcadf[0][i],pcadf[2][i]))
    
 
ax.add_patch(Circle((0,0),1.7,fill=False,color='red',linewidth=0.7,linestyle='--'))
# plt.xlim(-10,10)
# plt.ylim(-5,5)


#maybe what i can do with this is choose which are the 'interesting' clouds based on: pick top 3 furthest away from the center, and top 3 closest to the center


pca3.explained_variance_ratio_ #adds up to 97%


pc_coords = pcadf.to_numpy()


pca_ds = pd.DataFrame([(pcat['name'][i],np.sqrt(pc_coords[i][0]**2 + pc_coords[i][1]**2 + pc_coords[i][2]**2)) for i in range(22)]).sort_values(1)


pca_ds['rank'] = np.arange(22)


plt.figure(figsize=(10,10))

plt.scatter(pcat['mass_HI<Msun>']/pcat['mass_GMC<Msun>'],
            pcat['ang_mom_GMC<erg.s>']/pcat['ang_mom_HI<erg.s>'],c=pca_ds['rank'])
plt.loglog()
plt.show


hmom0s


surden


#okay so the pixels that we would use for pca arent all gonna be tthe same ones
#as how the size was determined, but thats okay
for i in range(22):
    surden = (np.array(cmom0s[i])*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value
    d0dy,d0dx = np.gradient(surden)
    d0 = np.sqrt(d0dy**2 + d0dx**2)

    d1dy,d1dx = np.gradient(cmom1s[i])
    d1 = np.sqrt(d1dy**2 + d1dx**2)

    d2dy,d2dx = np.gradient(cmom2s[i])
    d2 = np.sqrt(d2dy**2 + d2dx**2)


    df = pd.DataFrame(np.array([np.ravel(coords[i][0]),np.ravel(coords[i][1]),
                                np.ravel(surden),
    np.ravel(cmom1s[i]),np.ravel(cmom2s[i]),
    np.ravel(d0),np.ravel(d1),np.ravel(d2)]).transpose(),
    columns=['l','b','S','v','s','dw','dv','ds'])

#     df = df.dropna() #getting rid of this makes the math closer between rnew and rold
    #but gonna need to drop na for pca

    df = df[df['S']>5]

    dist = pcat['distance<pc>'][i]
    rnew = np.sqrt((len(df)*(0.125*u.deg.to(u.rad)*dist*u.pc)**2)/np.pi).value
    rold = pcat['radius_GMC<pc>'][i]
#     print(rnew,rold)
    print(pcat['name'][i],np.around(np.abs(np.around(rnew,2)-np.around(rold,2))/np.around(rnew,2)*100,2))
    
    print('\n')


#okay so the pixels that we would use for pca arent all gonna be tthe same ones
#as how the size was determined, but thats okay
dfs=[]
for i in range(22):
    surden = (np.array(cmom0s[i])*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value
    d0dy,d0dx = np.gradient(surden)
    d0 = np.sqrt(d0dy**2 + d0dx**2)

    d1dy,d1dx = np.gradient(cmom1s[i])
    d1 = np.sqrt(d1dy**2 + d1dx**2)

    d2dy,d2dx = np.gradient(cmom2s[i])
    d2 = np.sqrt(d2dy**2 + d2dx**2)


    df = pd.DataFrame(np.array([np.ravel(coords[i][0]),np.ravel(coords[i][1]),
                                np.ravel(surden),
    np.ravel(cmom1s[i]),np.ravel(cmom2s[i]),
    np.ravel(d0),np.ravel(d1),np.ravel(d2)]).transpose(),
    columns=['l','b','S','v','s','dw','dv','ds'])

    df = df.dropna() #getting rid of this makes the math closer between rnew and rold
    #but gonna need to drop na for pca

    df = df[df['S']>2.5]
    df = df[df['s']>1e-3]
    
    df['cloud'] = i#pcat['name'][i]
    dfs.append(df)
#     print(df.head())
#     print(pca3.fit_transform(df).shape)


#     dist = pcat['distance<pc>'][i]
#     rnew = np.sqrt((len(df)*(0.125*u.deg.to(u.rad)*dist*u.pc)**2)/np.pi).value
#     rold = pcat['radius_GMC<pc>'][i]
# #     print(rnew,rold)
#     print(pcat['name'][i],np.around(np.abs(np.around(rnew,2)-np.around(rold,2))/np.around(rnew,2)*100,2))
    
#     print('\n')


alltg = pd.concat(dfs)


allpaca = pca3.fit_transform(alltg[alltg.columns[2:-1]])


comps = pd.DataFrame(allpaca,columns=['pc1','pc2','pc3'])


np.array((alltg.columns,np.nanstd(alltg,axis=0))).transpose()


np.array((pca3.components_[0],alltg.columns[2:-1])).transpose()


plt.scatter(comps['pc1'][::10],comps['pc2'][::10],c=alltg['cloud'][::10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc2'][::10],comps['pc3'][::10],c=alltg['cloud'][::10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc1'][::10],comps['pc3'][::10],c=alltg['cloud'][::10],cmap='Paired',s=5)
plt.colorbar()
#orange is maddalena


#the thing with the largest standard deviation (S) contributes the most to the 
#principal components, so need to scale the data first


from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

pipeline = Pipeline([('scaling', StandardScaler()), ('pca', PCA(n_components=3))])

# pipeline.fit_transform(data)


scaled_alpaca = pipeline.fit_transform(alltg[alltg.columns[2:-1]])


comps = pd.DataFrame(scaled_alpaca,columns=['pc1','pc2','pc3']) 


plt.scatter(comps['pc1'][::10],comps['pc2'][::10],c=alltg['cloud'][::10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc1'][::-10],comps['pc2'][::-10],c=alltg['cloud'][::-10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc2'][::10],comps['pc3'][::10],c=alltg['cloud'][::10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc2'][::-10],comps['pc3'][::-10],c=alltg['cloud'][::-10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc1'][::10],comps['pc3'][::10],c=alltg['cloud'][::10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc1'][::-10],comps['pc3'][::-10],c=alltg['cloud'][::-10],cmap='Paired',s=5)
plt.colorbar()


plt.scatter(comps['pc1'][::-10],comps['pc2'][::-10],c=alltg['cloud'][::-10],cmap='Paired',s=2)
plt.colorbar()
plt.xlim(-2.5,2.5)
plt.ylim(-5,5)


#hmm still not seeing much difference
#other than (some?all?) maddalena pxs are off doing their own thing along the third axis
#and maybe aquila is also clustered tg?


scaled_alpaca_cat = pipeline.fit_transform(logdata)


comps = pd.DataFrame(scaled_alpaca_cat,columns=['pc1','pc2','pc3']) 


plt.figure(figsize=(10,10))
plt.scatter(comps['pc2'],comps['pc3'])
for i in range(22):
    plt.annotate(pcat['name'][i],(comps['pc2'][i],comps['pc3'][i]))

# plt.colorbar()


#hmm okay the edge clouds from before are still the edge clouds ish
#but it looks like there is much less information retained in just two components
#also figure out how to gget the explained variance ratio out of the pipeline


pcat


for i in range(22):
    surden = (np.array(cmom0s[i])*1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value
    
    d0dy,d0dx = np.gradient(surden)
    d0 = np.sqrt(d0dy**2 + d0dx**2)

    d1dy,d1dx = np.gradient(cmom1s[i])
    d1 = np.sqrt(d1dy**2 + d1dx**2)

    d2dy,d2dx = np.gradient(cmom2s[i])
    d2 = np.sqrt(d2dy**2 + d2dx**2)
    
    
    df = pd.DataFrame(np.array([np.ravel(coords[i][0]),np.ravel(coords[i][1]),
                                np.ravel(surden),
    np.ravel(cmom1s[i]),np.ravel(cmom2s[i]),
    np.ravel(d0),np.ravel(d1),np.ravel(d2)]).transpose(),
    columns=['l','b','S','v','s','dw','dv','ds'])

    df = df.dropna()
    df = df[df['S']>3]

    
    l0,b0 = np.nanmean(df['l']),np.nanmean(df['b'])
    rdg = np.sqrt((np.array(df['l'])-l0)**2 + (np.array(df['b'])-b0)**2)
    rpc = rdg*u.deg.to(u.rad)*pcat['distance<pc>'][i]
    
    df['r'] = rpc
#     print(len(df))
    binned_by_r = df.groupby(pd.cut(df['r'],bins=np.arange(0,np.nanmax(df['r']),0.5*np.nanstd(df['r'])))).mean()
    
    binned_by_r = binned_by_r.dropna()
    
#     print((pcat['radius_GMC<pc>'][i]/np.nanmax(df['r']))**-1)
#     print((np.nanstd(pcat['lon_GMC<deg>'][i])/np.nanmax(rdg))**-1)
#     print((np.nanstd(pcat['lat_GMC<deg>'][i])/np.nanmax(rdg))**-1)
#     print('\n')

    brs = np.zeros(len(binned_by_r))
    
    for i,idx in zip(range(1,len(binned_by_r)),binned_by_r.index[1:]):
        brs[i] = idx.mid
        

#     plt.scatter(df['r'],df['S'],s=5)
    plt.plot(brs,binned_by_r['S']/np.nanmax(binned_by_r['S']),linewidth=0.8,marker='.',color='black',alpha=0.5)
    
    plt.loglog()
#     plt.colorbar()
    plt.show()


cmom0_list = fits.HDUList([fits.PrimaryHDU(cmom0s[0])]+[fits.ImageHDU(cmom0s[i]) for i in range(1,22)])
cmom1_list = fits.HDUList([fits.PrimaryHDU(cmom1s[0])]+[fits.ImageHDU(cmom1s[i]) for i in range(1,22)])
cmom2_list = fits.HDUList([fits.PrimaryHDU(cmom2s[0])]+[fits.ImageHDU(cmom2s[i]) for i in range(1,22)])


cmom0_list.writeto('mom0_co.fits')
cmom1_list.writeto('mom1_co.fits')
cmom2_list.writeto('mom2_co.fits')


hmom0_list = fits.HDUList([fits.PrimaryHDU(hmom0s[0])]+[fits.ImageHDU(hmom0s[i]) for i in range(1,22)])
hmom1_list = fits.HDUList([fits.PrimaryHDU(hmom1s[0])]+[fits.ImageHDU(hmom1s[i]) for i in range(1,22)])
hmom2_list = fits.HDUList([fits.PrimaryHDU(hmom2s[0])]+[fits.ImageHDU(hmom2s[i]) for i in range(1,22)])


hmom0_list.writeto('mom0_hi.fits')
hmom1_list.writeto('mom1_hi.fits')
hmom2_list.writeto('mom2_hi.fits')


QTable(fits.open('post_flash_cat.fits')[1].data)


pcat


#also start keeping track of uncertainty sources
#for example, the threshold used for defining "significant" emission -> especially for the littler clouds this makes 
#a big difference (eg ~17% for aquila when choosing Sigma = 5 vs Sigma = 3)
#uncertainties in distance
