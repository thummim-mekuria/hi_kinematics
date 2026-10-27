#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable

import fun as fns
import datacube as dc


hi4pi,dht = fits.open('both_same_grid_720d.fits')

h3sig35kms = (hi4pi.data*(hi4pi.data>0.043*3))[50:105]
d3sig35kms = (dht.data*(dht.data>0.2*3))[50:105]


zcat = QTable(fits.open('final_cat.fits')[1].data)


h0 = np.nansum(h3sig35kms,axis=0)
c0 = np.nansum(d3sig35kms,axis=0)

h1 = fns.get_mom1(h3sig35kms,np.arange(55)*1.3-35)
c1 = fns.get_mom1(d3sig35kms,np.arange(55)*1.3-35)

h2 = fns.get_momN(h3sig35kms,np.arange(55)*1.3-35,2)
c2 = fns.get_momN(d3sig35kms,np.arange(55)*1.3-35,2)


for i in range(17):
    plt.contourf(np.log10(h0*1.3*1.82e18),extent=[360,-360,-50,50],
                 levels=np.arange(20,22,0.2))
    
    plt.contourf(np.log10(c0*1.3*1.97e20),extent=[360,-360,-50,50],
                 levels=np.arange(21,22,0.2),cmap='binary')
    
    plt.contour(np.log10(c0*1.3*1.97e20),extent=[360,-360,-50,50],
                 levels=[20.6],colors='red',linewidths=0.6)
    
    
    plt.xlim(zcat['lon_GMC<deg>'][i]+np.array((5,-5)))
    plt.ylim(zcat['lat_GMC<deg>'][i]+np.array((-5,5)))
    
    plt.vlines(zcat['lon_GMC<deg>'][i],
               zcat['lat_GMC<deg>'][i][0],zcat['lat_GMC<deg>'][i][1])
    
    plt.vlines(fake_zcat['lon_GMC<deg>'][i],
               fake_zcat['lat_GMC<deg>'][i][0],
               fake_zcat['lat_GMC<deg>'][i][1],
              colors='orange')
    
    
    plt.hlines(zcat['lat_GMC<deg>'][i],
           zcat['lon_GMC<deg>'][i][0],zcat['lon_GMC<deg>'][i][1])
    
    plt.hlines(fake_zcat['lat_GMC<deg>'][i],
              fake_zcat['lon_GMC<deg>'][i][0],
              fake_zcat['lon_GMC<deg>'][i][1],
              colors='orange')
    
    plt.show()


zcat[['names','lon_GMC<deg>','lat_GMC<deg>']] 


zcat


fake_zcat = pd.DataFrame([('california', [ 166.5,  154.], [-13,  -6.]),
('chamaeleon', [ -54.,  -67.], [-20., -12.]),
('cma_ob1', [-132.5, -138.5], [ -3.,   0.5]),
('crossbones', [-138.5, -142.], [-12.,  -6.5]),
('gemini_ob1', [-164.5, -173.5], [ -3.5,   2.5]),
('hercules', [  49.,   41.], [  7.,  11.]),
('lacerta', [  97.5,   95], [-16.5,  -8.]),
('maddalena', [-140.5, -146.], [ -4.5,   -1.]),
('mon_ob1', [-155., -162.], [  -0.5,   3.5]),
('mon_r2', [-141.5, -149.5], [-14.5, -10.5]),
('orion_a', [-143., -155.], [-21., -17.]),
('orion_b', [-151., -157.5], [-17., -10.]),
('orion_l', [-159, -172.5], [-18.,  -7.]),
('perseus', [ 162.,  155.], [-26., -15.5]),
('polaris', [ 127.,  118.], [ 20.,  31.]),
('rosette', [-151., -155.], [ -3.5,  -0.5]),
('taurus', [ 177.,  162.], [-21., -11.])],
columns=['names','lon_GMC<deg>','lat_GMC<deg>'])


for i in range(17):
    lo,li = ((np.array(fake_zcat['lon_GMC<deg>'][i])-360)/-0.125).astype('int')
    bo,bi = ((np.array(fake_zcat['lat_GMC<deg>'][i])+50)/0.125).astype('int')
    
    mom2 = np.nanmean(np.sqrt(h2[bo:bi,lo:li]))
    mom1 = np.nanmean(h1[bo:bi,lo:li])
    
    v0,v1 = fake_zcat_v['vel<km/s>'][i]
    vo,vi = ((np.array([v0,v1])+100)/1.3).astype('int')
    
    hsp = np.nanmean(hi4pi.data[:,bo:bi,lo:li],axis=(1,2))
    hsp2 = np.nanmean(hi4pi.data[vo:vi,bo:bi,lo:li],axis=(1,2))
    
    
    csp = np.nanmean(dht.data[:,bo:bi,lo:li],axis=(1,2))
    csp2 = np.nanmean(dht.data[vo:vi,bo:bi,lo:li],axis=(1,2))
    
    plt.plot(np.arange(155)*1.3-100,hsp)
    
    
    plt.plot(np.arange(vi-vo)*1.3+v0,hsp2)
    
    
    plt.plot(np.arange(155)*1.3-100,30*csp)
    
    
    plt.plot(np.arange(vi-vo)*1.3+v0,30*csp2)
    
    
#     print('i')
    plt.hlines(np.nanmax(hsp)*0.5,
               mom1-1.5*np.sqrt(2*np.log(2))*mom2,
               mom1+1.5*np.sqrt(2*np.log(2))*mom2)
#     print('ii')
    plt.show()
    print(fake_zcat['names'][i],mom1,mom2,np.nanmax(sp)*0.5)
    
#     print(lo,li,bo,bi)


fake_zcat_v = pd.DataFrame([('california', [ -35,35],),
('chamaeleon', [ -18,35],),
('cma_ob1',[ -35,29],),
('crossbones', [ -35,29],),
('gemini_ob1', [ -35,40], ),
('hercules', [ -35,60], ),
('lacerta',[ -25,35], ),
('maddalena', [ -35,42],),
('mon_ob1', [ -35,15], ),
('mon_r2', [ -35,27],),
('orion_a',[ -35,40], ),
('orion_b',[ -35,30], ),
('orion_l',[ -35,35], ),
('perseus', [ -25,35], ),
('polaris', [ -35,35], ),
('rosette', [ -35,30], ),
('taurus', [ -12,35], )],
columns=['names','vel<km/s>'])


news = pd.DataFrame([('aquila_s',[41,28],[-20,-12]),
('camelopardalis',[150,144],[16.5,18.5]),
('corona_australis',[3.5,-2],[-23.5,-16]),
('ophiuchus',[11.5,-10.],[12.5,25]),
('lupus',[-8,-27],[3.5,25.5])],columns=['names','lon_GMC<deg>','lat_GMC<deg>'])


for i in range(5):
    plt.contourf(np.log10(h0*1.3*1.82e18),extent=[360,-360,-50,50],
                 levels=np.arange(20,22,0.2))
    
    plt.contourf(np.log10(c0*1.3*1.97e20),extent=[360,-360,-50,50],
                 levels=np.arange(19,22,0.2),cmap='binary')
    
    plt.contour(np.log10(c0*1.3*1.97e20),extent=[360,-360,-50,50],
                 levels=[20.6],colors='red',linewidths=0.6)
    
    
    plt.xlim(news['lon_GMC<deg>'][i]+np.array((5,-5)))
    plt.ylim(news['lat_GMC<deg>'][i]+np.array((-5,5)))
    
    plt.vlines(news['lon_GMC<deg>'][i],
               news['lat_GMC<deg>'][i][0],news['lat_GMC<deg>'][i][1])
    

    plt.hlines(news['lat_GMC<deg>'][i],
           news['lon_GMC<deg>'][i][0],news['lon_GMC<deg>'][i][1])
    

    plt.show()


news_v = pd.DataFrame([('aquila_s',[-35,15]),
('camelopardalis',[-25,35]),
('corona_australis',[-35,35]),
('ophiuchus',[-10,35]),
('lupus',[-10,35])],columns=['names','vel<km/s>'])


for i in range(5):
    lo,li = ((np.array(news['lon_GMC<deg>'][i])-360)/-0.125).astype('int')
    bo,bi = ((np.array(news['lat_GMC<deg>'][i])+50)/0.125).astype('int')
    
    mom2 = np.nanmean(np.sqrt(h2[bo:bi,lo:li]))
    mom1 = np.nanmean(h1[bo:bi,lo:li])
    
    v0,v1 = news_v['vel<km/s>'][i]
    vo,vi = ((np.array([v0,v1])+100)/1.3).astype('int')
    
    hsp = np.nanmean(hi4pi.data[:,bo:bi,lo:li],axis=(1,2))
    hsp2 = np.nanmean(hi4pi.data[vo:vi,bo:bi,lo:li],axis=(1,2))
#     hsp2 = np.nanmean(h3sig35kms[:,bo:bi,lo:li],axis=(1,2))

    
    csp = np.nanmean(dht.data[:,bo:bi,lo:li],axis=(1,2))
    csp2 = np.nanmean(dht.data[vo:vi,bo:bi,lo:li],axis=(1,2))
#     csp2 = np.nanmean(d3sig35kms[:,bo:bi,lo:li],axis=(1,2))
    
    plt.plot(np.arange(155)*1.3-100,hsp)
    
    
    plt.plot(np.arange(vi-vo)*1.3+v0,hsp2)
#     plt.plot(np.arange(55)*1.3-35,hsp2)

    
    plt.plot(np.arange(155)*1.3-100,100*csp)
    
    
    plt.plot(np.arange(vi-vo)*1.3+v0,30*csp2)
    
#     plt.plot(np.arange(55)*1.3-35,100*csp2)
    

#     print('i')
    plt.hlines(np.nanmax(hsp)*0.5,
               mom1-1.5*np.sqrt(2*np.log(2))*mom2,
               mom1+1.5*np.sqrt(2*np.log(2))*mom2)
#     print('ii')
    plt.show()
    print(news['names'][i],mom1,mom2,np.nanmax(sp)*0.5)
    
    
#     plt.imshow(np.nansum())
    
#     print(lo,li,bo,bi)


#next steps -> plot patch by patch spectra on top of moment maps
#re-include the galactic-center-wards clouds (corona,ophiuchus,aqquila?,lupus,camelopardalis)

^done

#mass etc
#big picture plots (the one from h-alpha is v appealing in term of 
#pointing out which clouds are moving at significantly different v than surrounding hi)
#decide whether to cluster clouds together by accumulation radii (or - do the test of the
#boundary of the sum of the parts, vs the acc boundary defined by the sum of masses)
#but for crowded fields -> risks including other co not defined in the catalog, which distance would be used
#but the justification for clustering them is -> different behavior might emerge
#from clouds that have lots of co going on vs clouds that are relatively isolated
#definitey seems to be a bimodal distribution in genre


oldcatnew = fake_zcat.merge(fake_zcat_v,on='names')


oldcatnew


newcatnew


newcatnew = news.merge(news_v,on='names')


newcatnew.columns = ['names','lon','lat','vel']
oldcatnew.columns = ['names','lon','lat','vel']


meh = newcatnew.merge(oldcatnew,how='outer',on=['names'])


meh


newzlons = [np.array([41, 28]), np.array([150, 144]), np.array([3.5, -2]),
       np.array([11.5, -11.5]), np.array([-8, -27]),
            np.array([166.5, 154.0]), np.array([-54.0, -67.0]),
       np.array([-132.5, -138.5]), np.array([-138.5, -142.0]),
       np.array([-164.5, -173.5]), np.array([49.0, 41.0]), np.array([97.5, 95]),
       np.array([-140.5, -146.0]), np.array([-155.0, -162.0]),
       np.array([-141.5, -149.5]), np.array([-143.0, -155.0]),
       np.array([-151.0, -157.5]), np.array([-159, -172.5]), np.array([162.0, 155.0]),
       np.array([127.0, 118.0]), np.array([-151.0, -155.0]), np.array([177.0, 162.0])]


newzlats = [np.array([-20, -12]), np.array([16.5, 18.5]), np.array([-23.5, -16]),
       np.array([12.5, 25]),np.array([3.5, 25.5]),
            np.array([-13, -6.0]), np.array([-20.0, -12.0]),
       np.array([-3.0, 0.5]), np.array([-12.0, -6.5]), np.array([-3.5, 2.5]),
       np.array([7.0, 11.0]), np.array([-16.5, -8.0]), np.array([-4.5, -1.0]),
       np.array([-0.5, 3.5]), np.array([-14.5, -10.5]), np.array([-21.0, -17.0]),
       np.array([-17.0, -10.0]), list([-18.0, -7.0]), np.array([-26.0, -15.5]),
       np.array([20.0, 31.0]), list([-3.5, -0.5]), np.array([-21.0, -11.0])]


newzlatsnames = meh['names']


newzvels = [np.array([-35, 15]), np.array([-25, 35]), np.array([-35, 35]), np.array([-10, 35]),
            np.array([-10, 35]),np.array([-35, 35]), np.array([-18, 35]),
       np.array([-35, 29]), np.array([-35, 29]), np.array([-35, 40]), np.array([-35, 60]),
       np.array([-25, 35]), np.array([-35, 42]), np.array([-35, 15]), np.array([-35, 27]),
       np.array([-35, 40]), np.array([-35, 30]), np.array([-35, 35]), np.array([-25, 35]),
       np.array([-35, 35]), np.array([-35, 30]), np.array([-12, 35])]


new_cat = QTable([newzlatsnames,newzlons,newzlats,newzvels],
       names=['name','lon_GMC<deg>','lat_GMC<deg>','vel<km/s>'])


new_cat['name'] = ['aquila_s',
 'camelopardalis',
 'corona_australis',
 'ophiuchus',
 'lupus',
 'california',
 'chamaeleon',
 'cma_ob1',
 'crossbones',
 'gemini_ob1',
 'hercules',
 'lacerta',
 'maddalena',
 'mon_ob1',
 'mon_r2',
 'orion_a',
 'orion_b',
 'orion_l',
 'perseus',
 'polaris',
 'rosette',
 'taurus']


new_cat.write('re_expanded_cat.fits',overwrite=True)
