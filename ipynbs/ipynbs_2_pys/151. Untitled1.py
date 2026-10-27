#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
import numpy as np
import fun as fns
import matplotlib.pyplot as plt
from astropy.table import QTable

import pandas as pd

import astropy.units as u
import astropy.constants as c


hi4pi,dht = fits.open('both_same_grid_720d.fits')
dust = fits.open('mean_and_std_lbd.fits')[1].data


names = [
    'Aquila South',
    'California',
    'Camelopardalis',
    'Chamaeleon',
    'Canis Major OB1',
    'Corona Australis',
    'Crossbones',
    'Gemini OB1',
    'Hercules',
    'Lacerta',
    'Lupus',
    'Maddalena',
    'Monoceros OB1',
    'Monoceros R2',
    'Ophiuchus',
    'Orion A',
    'Orion B',
    r'$\lambda$-Orionis',
    'Perseus',
    'Polaris',
    'Rosette',
    'Taurus'  
]


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


cat = QTable(fits.open('sighhplease.fits')[1].data)


bounds = pd.read_csv('../Untitled Folder/lbv_bounds_both.csv')


hi_rad = np.array([22.24767, 123.66868, 21.751354, 46.572394, 110.41117, 28.498626, 
 69.668984, 245.52689, 38.891018, 43.873859, 64.646119, 153.29056, 104.20971, 112.08244,
 62.394432, 108.40114, 94.221693, 71.788496, 76.857215, 
           63.525783, 109.73185, 57.715274]) #from 145


h2h2 = fits.open('../Untitled Folder/h2h2.fits')


hihi = fits.open('../Untitled Folder/hihi.fits')


h2hi = fits.open('../Untitled Folder/h2hi.fits')


hi_msks = [(fits.open('../Untitled Folder/hi_masks.fits')[i].data-1).astype('bool')
           for i in range(22)]
h2_msks = [(fits.open('../Untitled Folder/h2_masks.fits')[i].data-1).astype('bool') 
           for i in range(22)]


from copy import deepcopy

h2_msks_hisz = deepcopy(hi_msks)
# for i in range(22):
#     h2_msks_hisz[i][bounds['boc'][i]-bounds['boh'][i]:bounds['bic'][i]-bounds['boh'][i],
#                     bounds['loc'][i]-bounds['loh'][i]:bounds['lic'][i]-bounds['loh'][i]] = ~h2_msks[i] 

#     h2_msks_hisz[i] = ~h2_msks_hisz[i]*hi_msks[i]


do,di = (cat['distance<pc>'] - hi_rad).astype('int'),(cat['distance<pc>'] + hi_rad).astype('int')


dlo = (((bounds['loh']*-0.125+360) - 180)/-0.229).astype('int')
dli = (((bounds['lih']*-0.125+360) - 180)/-0.229).astype('int')


dbo = (((bounds['boh']*0.125-50) + 90)/0.229).astype('int')
dbi = (((bounds['bih']*0.125-50) + 90)/0.229).astype('int')


d0ths = [np.nansum(2.8*dust[do[i]-68:di[i]-68,
                        max(0,dbo[i]):dbi[i],
                        max(0,dlo[i]):dli[i]],axis=0)*1.87e21
                   
                   for i in range(22)]


hi0ths = [np.nansum(hihi[i].data[bounds['vo'][i]:bounds['vi'][i],
                               :,:],axis=0)*1.82e18
          for i in range(22)]

h2hi0ths = [np.where(~h2_msks_hisz[i],
                     0,
                     np.nansum(h2hi[i].data[bounds['vo'][i]:bounds['vi'][i],
                               :,:],axis=0)*1.97e20)
          for i in range(22)]


d0_msks = [(
np.linalg.norm(np.indices(d0ths[i].shape)-d0ths[i].shape[0]/2,axis=0)*\
    0.229*u.deg.to(u.rad)*cat['distance<pc>'][i] < hi_rad[i]
)
 for i in range(22)]


lo_dg,li_dg = np.array((bounds['loh'],bounds['lih']))*-0.125+360
bo_dg,bi_dg = np.array((bounds['boh'],bounds['bih']))*0.125 -50


lvls = np.linspace(20,22,30)


c_map = 'cubehelix'


for i in [8,10,18]:
    fig,axs = plt.subplots(1,2,dpi=200,sharey=True,figsize=(10,5))
    
    lvls = np.linspace(19.9,22,20)
    
    axs[0].contourf(np.where(~hi_msks[i],
                             np.nan,np.log10(hi0ths[i]+h2hi0ths[i])),
                    extent=[lo_dg[i],li_dg[i],bo_dg[i],bi_dg[i]],
                    cmap=c_map,
                    levels=lvls,extend='both'),
    axs[0].set_aspect(1)

    axs[0].set_title('HI + $\mathrm{H_2}$')
    
    
    axs[1].contourf(np.where(~d0_msks[i],
                             np.nan,np.log10(d0ths[i])),
                    extent=[lo_dg[i],li_dg[i],bo_dg[i],bi_dg[i]],
                    cmap=c_map,
                    levels=lvls,extend='both'),
#     plt.colorbar(),
    axs[1].set_aspect(1),
    axs[1].set_title('Dust [E24]')
    
    
    for ax in axs:
        ax.invert_xaxis()
        ax.set_xticks(np.around(np.linspace(lo_dg[i],li_dg[i],3),1),
                     [str(tck)+'$\degree$'
                         for tck in np.around(np.linspace(lo_dg[i],li_dg[i],3),1)])
        
        ax.set_yticks(np.around(np.linspace(bo_dg[i],bi_dg[i],3),1),
                     [str(tck)+'$\degree$'
                         for tck in np.around(np.linspace(bo_dg[i],bi_dg[i],3),1)])
        
    
    fig.suptitle(names[i])
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    
#     fig.colorbar(axs[1])

    fig.savefig(cat['name'][i]+'_hih2dust',bbox_inches='tight')
    plt.show()


columns = np.array('x,y,z,vx,vy,vz,a,e,inc,Omega,omega,f,m,r'.split(','))


info = np.array(([-0.06671738977689852,-0.06285577673870507,0.00012636375509346694,14.232301868666564,-15.10670545626774,-0.013249228255757499,0.0916628318123935,1.390440899958571e-16,0.001519199559917549,1.8927382427267643,4.825406386228014,3.4622373064994374,0.00022231456930986776,0.0005077481990176487],
[-0.06010428226135066,-0.10347979749600535,-0.0006269953560445405,15.72739139771849,-9.136458963582115,0.1652334329077973,0.11965079557883092,2.740598420349355e-19,0.010489502776539174,-1.5738416683017606,4.1892535910251745,1.5707963267948966,0.002074145112982361,0.0005256839527997373],
[-0.10249219772900953,-0.03455191267242178,-0.0003449642202240893,6.1305419310114475,-18.19606612984009,-0.08826887819965291,0.10795547320191465,3.0953031151115705e-16,0.005612814325503162,-0.2795540805744614,0.6320773550714227,3.1126389555188623,0.0025284632633447504,0.0005256839527997373],
[-0.5678872217020126,0.6471312410138681,0.00846430474579255,-5.0535213370136045,-4.5428331228099355,0.050776641725353615,0.8609895836317306,1.5674332016048333e-16,0.012344863844420538,1.3689044620563924,0.7407417349241623,0.18086954266295496,0.00025016709898023084,0.0005158605324469098],
[1.3643496156481947,0.6990599389600342,-0.022499049338249765,-2.2730297065697034,4.463998273417785,-0.03342635685842273,1.533728630607286,1.4577868844945375e-16,0.01608143514534634,2.4663504402218686,3.9400957573444444,0.35021571733773893,0.0029413858362958662,0.0005256839527997373],
[0.004129647258136294,-0.04165503962581716,-0.0003242347548905357,30.302627707538242,0.420282627327971,-0.03168534267956497,0.043548193133090926,1.253015421591147e-16,0.006029052895045235,0.19044499247827792,1.2591963202188499,3.278478995460171,0.0023111451828164375,0.0005256839527997373],
[6.542459689582004,0.3635818870095953,-0.029602710968548177,-0.02414683741137992,2.4101105631353077,0.0047915940529181566,6.549081120058594,2.892319998960693e-17,0.004916686146880367,1.2164908635656355,3.551165952229615,1.5707963267948974,0.0007001372732119486,0.0005256839527997373],
[2.2705688633946526,-2.15062667399239,0.02793813313206451,2.579116033793892,2.539664278841728,-0.018351590218281567,3.1231527627900224,3.2521157981549473e-16,0.010329863008583268,-2.8495756047133676,2.6383671484452957,5.733865168536008,0.002314922176023078,0.0005256839527997373],
[-0.0853400875665874,-0.10306155363340645,0.0013192668963583455,12.151406024053172,-11.877326079220374,-0.039378417414862865,0.14048912209074102,2.1905153086755096e-16,0.009832066131770307,2.1255996546592146,4.527097148026105,3.566017118286897,6.312004008560186e-05,0.00038317378703589555],
[-2.027932346654878,-6.025485942731953,0.03851875380969535,2.4946913931592394,-0.8529897383445576,0.03084785059380855,6.359179289782496,1.7187303488630144e-16,0.013737272984408305,-2.3544523376313076,6.13224018659343,0.6078847784094572,0.0006975394546974033,0.0005256839527997373]))


info.shape


df = pd.DataFrame(info,columns=columns)


df2 = pd.DataFrame(info*1e10,columns=columns)


plt.hist(df2['m']/1e10)


df2['m']/1e10


df['m']


bighimom2 = np.sqrt(fns.get_momN(hi4pi.data,np.arange(155)*1.3-100,2))


bighimom1 = fns.get_mom1(hi4pi.data,np.arange(155)*1.3-100)


plt.figure(dpi=200,figsize=(13,7))
plt.contourf(bighimom2,extent=[360,-360,-50,50],cmap='Paired',levels=np.linspace(0,20,10),extend='both')
plt.xlim(180,-180)
plt.colorbar()


plt.figure(dpi=200,figsize=(13,7))
plt.contourf(bighimom2,extent=[360,-360,-50,50],levels=np.linspace(0,40,10),extend='both')
plt.contour(bighimom1,extent=[360,-360,-50,50],levels=20,extend='both',cmap='RdGy')

plt.xlim(180,-180)
# plt.colorbar()


import seaborn as sns


sns.kdeplot(x=np.ravel(bighimom2)[::1387], y=np.ravel(np.abs(bighimom1))[::1387])


plt.scatter([::875],[::875],s=1e-1)


bounds = pd.read_csv('../Untitled Folder/lbv_bounds_both.csv')


h2hi = [dht.data[:,bounds['boh'][i]:bounds['bih'][i],bounds['loh'][i]:bounds['lih'][i]] for i in range(22)]


fns.make_fits_from_list(h2hi,'h2hi.fits')


hlv = np.nansum(hi4pi.data,axis=1)*0.125*1.82e18
hbv = np.nansum(hi4pi.data,axis=2)*0.5*0.125*1.82e18


def make_fits_from_list(array,filename):
    hdulist = fits.HDUList([fits.PrimaryHDU(array[0])]+[fits.ImageHDU(array[i]) for i in range(1,len(array))])
    hdulist.writeto(filename,overwrite=True)


make_fits_from_list([hlv,hbv],'hi4pi_pv_maps.fits')


cat = QTable(fits.open('sighhplease.fits')[1].data)


cos_hisz = [dht.data[:,cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]] for i in range(22)]


fns.make_fits_from_list(cos_hisz,'h2_hisized.fits')


h0 = np.nansum(hi4pi.data,axis=0)*1.3*1.82e18


fits.PrimaryHDU(h0).writeto('bigh0.fits')


plt.figure(dpi=200)
plt.contourf(np.log10(h0),levels=np.arange(20,23,0.25),
             extent=[360,-360,-50,50])
plt.xlim(90,-90)


h1 = fns.get_mom1(hi4pi.data,np.arange(155)*1.3-100)


fits.PrimaryHDU(h1).writeto('bigh1.fits')


plt.figure(dpi=200)
plt.contourf(h1,levels=np.arange(-40,40,5),
             extent=[360,-360,-50,50],cmap='RdGy')
plt.xlim(90,-90)


plt.figure(dpi=200)
plt.contourf(h1,levels=np.arange(-40,40,5),
             extent=[360,-360,-50,50],cmap='RdGy')
plt.xlim(270,90)


c0 = np.nansum(dht.data,axis=0)*1.3*1.97e20


fits.PrimaryHDU(c0).writeto('big_zeroth_co.fits')


plt.contourf(np.log10(c0),levels=np.arange(20,23,0.1))


cat = QTable(fits.open('sighhplease.fits')[1].data)


lc = np.nanmean(cat['l'])


np.nanmean((cat['boh'],cat['bih']),axis=0) == np.nanmean((cat['boc'],cat['bic']),axis=0)


np.nanmean((cat['loh'],cat['lih']),axis=0) == np.nanmean((cat['loc'],cat['lic']),axis=0)


np.nanmean((cat['loc'],cat['lic']),axis=0)


cohi = [dht.data[:,
                   bc[i]-newr_px[i]:bc[i]+newr_px[i],
                   lc[i]-newr_px[i]:lc[i]+newr_px[i]] 
        for i in range(22)]


hi2hi = [hi4pi.data[:,
                   bc[i]-2*newr_px[i]:bc[i]+2*newr_px[i],
                   lc[i]-2*newr_px[i]:lc[i]+2*newr_px[i]] 
        for i in range(22)]

co2hi = [dht.data[:,
                   bc[i]-2*newr_px[i]:bc[i]+2*newr_px[i],
                   lc[i]-2*newr_px[i]:lc[i]+2*newr_px[i]] 
        for i in range(22)]


make_fits_from_list(coco,'h2_h2sized.fits')
make_fits_from_list(hico,'hi_h2sized.fits')
make_fits_from_list(hihi,'hi_hisized.fits')
make_fits_from_list(cohi,'co_hisized.fits')
make_fits_from_list(hi2hi,'hi_2xhisized.fits')
make_fits_from_list(co2hi,'co_2xhisized.fits')


cat['distance<pc>'][10]


dld = np.nansum(dust*2.8,axis=1)


plt.contourf(dld,extent=[180,-180,68,1250],levels=np.arange(3e-2,2,5e-2))
plt.ylim(70,200)
plt.xlim(50,20)


d01 = np.nansum(dust[80-68:100-68,:,:]*2.8,axis=0)


plt.contourf(d01,levels=np.arange(1e-2,0.3,1e-2),extent=[180,-180,-90,90])
plt.ylim(-40,60)
plt.xlim(100,-50)


d02 = np.nansum(dust[110-68:160-68,:,:]*2.8,axis=0)


plt.contourf(d02,levels=np.arange(1e-2,0.7,1e-2),extent=[180,-180,-90,90])
plt.ylim(-10,40)
plt.xlim(40,-30)
#hi ophiuchus!


plt.contourf(d02,levels=np.arange(1e-2,1.3,5e-2),extent=[180,-180,-90,90])
plt.contour(np.log10(c0),extent=[360,-360,-50,50],levels=np.arange(19.5,22,0.5),
            colors='red',linewidths=0.5)
plt.ylim(-25,-5)
plt.xlim(45,15)
#hi ophiuchus!


fits.open('mean_and_std_lbd.fits')[1].header


(np.arange(1572)-787)*-0.23


((np.arange(787)-394)*0.23)[230:380]


aq_dld = np.nansum(dust[:,230:380,:]*2.8,axis=1)


plt.contourf(aq_dld,extent=[180,-180,68,1250],levels=np.arange(3e-2,2,1e-1))
plt.ylim(90,170)
plt.xlim(45,15)
