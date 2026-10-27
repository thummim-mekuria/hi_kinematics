#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
from joblib import Parallel,delayed


cat = np.array((
['AqS',  41,     28,    -20,   -12,   -5.1,  9.2, 133],
['Cal',  174,    154,   -14,   -1,    -18.1, 18.3, 470],
['Cam',  150,    144,    16,    19,   -6.4,  4.0,  213],
['CMj', -132.5, -138.5, -4,     1.5,   9.2,  23.5, 1209],
['Chm', -54,    -67,    -20,   -12,   -3.8,  9.2,  183],
['CrA',  4,      -2,    -23.5, -16,    0.1,  10.5, 151],
['Xbn',  -138,   -142,   -12.5, -6.5,   4.0,  19.6, 886],
['Gem',  -164,   -174,   -4,     6,    -6.4,  14.4, 1786],
['Her',  49,     41,     7,     11.5, -10.3, 23.5, 227],
['Lac',  98,     94,    -17,   -7,    -6.4,  14.4, 503],
['Lup', -10,    -27,     5,     25.5, -7.7,  17.0, 189],
['Mad', -140,   -147,   -5,    -0.5,   13.1, 39.1, 2072],
['MOB', -155,   -162,   -0.5,   3.5,  -2.5,  15.7, 745],
['MR2',  -141,   -150,  -15,   -10,    2.7,  18.3, 778],
['Oph',  12,     -12,    3,     28,   -6.4,  13.1, 144],
['OrA', -140,   -156,   -24,   -16.5,  0.1,  14.4, 432],
['OrB', -149,   -158,   -18,   -9,    -1.2,  17,   423],
['LOr', -158.5, -173,   -19,   -6,    -11.6, 19.6, 402],
['Per',  162,    154,   -27.5, -15,   -1.2,  10.5, 294],
['Pol',  128,    117,    20,    30,   -12.9, 6.6,  352],
['Ros', -151,   -155,   -3.5,  -0.5,   1.4,  24.8, 1304],
['Tau',  177,    163,   -21,   -10,   -5.1,  17.0, 141],))


names, los, lis, bos, bis, vos, vis, ds = cat.transpose()
ds = ds.astype('float')


tab2 = np.array("AqS 0.88 0.87 5.45 22.25 1.3 3.2 87.6 80.8 105.43 108.46 0.80 1.23 \
Cal 151.78 98.49 38.66 123.67 5.3 8.5 51.0 58.1 4.94 -72.41 23.53 27.44 \
Cam 0.83 0.74 5.70 21.75 1.4 2.7 65.5 21.9 -125.43 -64.26 0.66 0.32 \
Can 108.02 54.36 32.43 110.41 2.6 3.8 40.7 12.0 144.76 152.0 2.31 1.02 \
Chm 8.11 4.05 8.92 46.57 1.8 3.2 94.1 15.2 115.04 72.57 13.21 4.51 \
CrA 1.86 1.08 4.74 28.5 1.1 2.7 111.2 2.2 -159.65 -61.01 0.77 0.05 \
Crs 27.14 20.78 19.72 69.67 1.9 4.2 65.8 31.8 20.21 82.82 7.90 4.76 \
Gem 11187.8 510.52 98.36 245.53 4.3 5.5 19.3 7.0 -46.43 -29.07 57.50 12.96 \
Her 4.72 6.71 7.25 38.89 3.1 7.5 323.1 68.5 -140.77 -93.36 5.24 3.20 \
Lac 6.78 4.85 13.34 43.87 2.5 3.7 176.8 15.6 27.39 132.85 9.71 0.92 \
Lup 21.68 13.98 16.32 64.65 2.3 5.8 67.0 15.8 -158.68 -66.35 5.51 2.03 \
Mad 289.06 174.49 57.88 153.29 4.2 7.3 50.6 19.9 44.86 103.88 52.31 14.42 \
MOB 90.82 52.48 24.66 104.21 3.0 4.7 107.7 7.8 97.29 -69.30 20.20 2.61 \
MR2 112.99 52.06 29.72 112.08 2.2 4.1 50.6 17.7 -147.28 92.80 13.79 6.88 \
Oph 19.49 14.48 15.04 62.39 2.4 4.8 79.2 29.6 -32.0 -61.49 5.53 3.55 \
OrA 102.22 28.20 24.40 108.40 2.9 3.7 70.6 4.8 110.03 130.47 12.96 1.75 \
OrB 67.13 37.99 17.91 94.22 2.9 4.6 44.3 9.5 -114.54 86.36 4.38 2.61 \
λOr 29.69 32.24 25.99 71.79 5.3 6.8 160.6 22.2 79.26 4.43 33.47 3.52 \
Per 36.43 11.49 15.50 76.86 3.0 3.1 159.7 18.8 31.59 12.41 11.83 3.42 \
Pol 20.57 5.64 17.15 63.53 2.2 5.1 63.8 22.8 31.14 54.87 5.79 2.84 \
Ros 106.03 96.62 30.71 109.73 4.2 6.2 87.6 10.5 -109.08 86.85 25.5 3.89 \
Tau 15.43 13.17 10.16 57.72 2.1 5.2 120.3 62.3 -105.49 -41.53 3.83 6.40".split(' '))

# mh2,mhi,rh2,rhi,sigh2,sighi,omgh2,omghi,thh2,thhi,jh2,jhi =  tab2.reshape(22,-1).transpose()[1:].astype('float')


tab2 = pd.DataFrame(tab2.reshape(22,-1),
                    columns = ['name','mh2','mhi','rh2','rhi',
                    'sigh2','sighi','omgh2','omghi','thh2','thhi','jh2','jhi'])


tab2['omgh2'].to_numpy().astype('float')*1e-3


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


# # refit for updated omega, theta


from scipy.optimize import curve_fit
def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y

def line(x,m,b):
    return m*x+b
    
def fit_gradient(mom1_map, dist, mask=0, weights=0, full=False):
    
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
        
    if type(weights)==int:#all the data is usable and/or a second moment map wasnt supplied
        weights = np.ones(mom1_map.shape).astype('bool')


    #the plane is pinned down to the geometric center of the map
    #the slope of the fit doesnt change with or without xo
    xdeg = (np.indices(mom1_map.shape)[1]-(mom1_map.shape[1]/2))*-0.125 
    ydeg = (np.indices(mom1_map.shape)[0]-(mom1_map.shape[0]/2))*0.125

    
    xpc = xdeg*u.deg.to(u.rad)*dist
    ypc = ydeg*u.deg.to(u.rad)*dist

    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map.data)
    
    
    raveldf = pd.DataFrame(np.array((
            np.ravel(xpc),np.ravel(ypc),
            np.ravel(mskd_mom1),
            np.ravel(weights),
            )).transpose(),
            columns=['x','y','v','unc'])
    
    
    noisy = raveldf
    raveldf = noisy[~(np.isnan(noisy['v']))] #get rid of the rows with masked data

    
    xslope,yslope,xyintercept = curve_fit(plane_rvld,
                                    np.array((raveldf['x'],raveldf['y'])),
                                    raveldf['v'],
                                    sigma=raveldf['unc'])[0]

    
    xmsig, ymsig, xybsig = np.sqrt(np.diag(
                                            curve_fit(plane_rvld,
                                                       np.array((raveldf['x'],raveldf['y'])),
                                                       raveldf['v'],
                                                       sigma=raveldf['unc'])[1]
                                    ))

    # print(xslope/xmsig, yslope/ymsig, xyintercept/xybsig)
    
    if not full:
        return xslope,yslope,xyintercept
    else:
        return {'fit_results':[xslope,yslope,xyintercept],
                '1_sigma_errors':[xmsig, ymsig, xybsig]}


def switch34(array):
    return np.concatenate((array[:3],array[4:5],array[3:4],array[5:]),axis=0)


hi_0s = switch34(fits.open('hi_0s.fits'))
h2_0s = switch34(fits.open('h2_0s.fits'))

hi_1s = switch34(fits.open('hi_1s.fits'))
h2_1s = switch34(fits.open('h2_1s.fits'))

hi_2s = switch34(fits.open('hi_2s.fits'))
h2_2s = switch34(fits.open('h2_2s.fits'))

hi_rs = [np.linalg.norm(np.indices(hi_1s[i].shape) - hi_1s[i].shape[0]/2,axis=0) for i in range(22)]
h2_rs = [np.linalg.norm(np.indices(h2_1s[i].shape) - np.array(h2_1s[i].shape)[:,np.newaxis,np.newaxis]/2,axis=0) for i in range(22)]

rhi_px = np.array([hi_1s[i].shape[0]/2 for i in range(22)]).astype('int')
hi_msk = [(hi_rs[i]<rhi_px[i])*(hi_0s[i].data>1e20) for i in range(22)]

h2_msk = [(switch34(fits.open('h2_0thmom_masks.fits'))[i].data - 1).astype('bool') for i in range(22)]


h2_msk[10] = h2_msk[10]* (h2_1s[10].data>2.)#lupus has a few extreme velocities
h2_msk[12] = h2_msk[12]* (h2_1s[12].data>-5.)#and mob2
h2_msk[14] = h2_msk[14]* (h2_1s[14].data<10.)#and and oph
h2_msk[19] = h2_msk[19]* (h2_1s[19].data<2.)#and pol
h2_msk[17] = h2_msk[17]* (h2_1s[17].data<15.)* (h2_1s[17].data>-5.)#and orl


h2_fitpars = np.array([fit_gradient(h2_1s[i].data,ds[i],h2_msk[i],h2_2s[i].data) for i in range(22)])
hi_fitpars = np.array([fit_gradient(hi_1s[i].data,ds[i],hi_msk[i],hi_2s[i].data) for i in range(22)])


# # stats on v vs rperp


h2dvdx = 1e-3*omgh2*np.cos(thh2*u.deg).value
h2dvdy = 1e-3*omgh2*np.sin(thh2*u.deg).value

hidvdx = 1e-3*omghi*np.cos(thhi*u.deg).value
hidvdy = 1e-3*omghi*np.sin(thhi*u.deg).value


hi_fitpars


h2_fitpars


h2_vos = h2_fitpars[:,2]
hi_vos = hi_fitpars[:,2]

hidvdxnew, hidvdynew = hi_fitpars[:,:2].transpose()
h2dvdxnew, h2dvdynew = h2_fitpars[:,:2].transpose()


magnt = (0.82, 0.27, 0.82)
green = (0.31, 0.57, 0.30)


names_full = np.array(['Aquila South', 'California', 'Camelopardalis', 'Canis Major', 'Chamaeleon', 'Corona Australis',
         'Crossbone', 'Gemini OB1', 'Hercules','Lacerta', 'Lupus', 'Maddalena', 
         'Monoceros OB1', 'Monoceros R2', 'Ophiuchus', 'Orion A', 'Orion B', r'$\lambda-$ Orionis',
         'Perseus', 'Polaris', 'Rosette', 'Taurus'])


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den


def plot_rvsv(i, show=False):  
    l0c,l1c = los.astype('float')[i], lis.astype('float')[i] 
    b0c,b1c = bos.astype('float')[i], bis.astype('float')[i] 

    lc,bc = np.nanmean((los.astype('float')[i],lis.astype('float')[i])), np.nanmean((bos.astype('float')[i],bis.astype('float')[i]))

    l0h,l1h = lc + rhi_px[i]*0.125, lc - rhi_px[i]*0.125
    b0h,b1h = bc - rhi_px[i]*0.125, bc + rhi_px[i]*0.125
    
    hlons = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0]
    hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1]
    
    
    clons = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0]
    clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1]
    

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)

    
    hbperp = (hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)/np.linalg.norm((hidvdx[i],hidvdy[i]))
    cbperp = (h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)/np.linalg.norm((h2dvdx[i],h2dvdy[i]))

    hbperp2 = (hidvdxnew[i]*hxs_pc + hidvdynew[i]*hys_pc)/np.linalg.norm((hidvdxnew[i],hidvdynew[i]))
    cbperp2 = (h2dvdxnew[i]*cxs_pc + h2dvdynew[i]*cys_pc)/np.linalg.norm((h2dvdxnew[i],h2dvdynew[i]))

    
    hv_model =  np.where(~hi_msk[i],np.nan, hi_vos[i] + hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)
    cv_model =  np.where(~h2_msk[i],np.nan, h2_vos[i] + h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)


    hv_model2 =  np.where(~hi_msk[i],np.nan, hi_vos[i] + hidvdxnew[i]*hxs_pc + hidvdynew[i]*hys_pc)
    cv_model2 =  np.where(~h2_msk[i],np.nan, h2_vos[i] + h2dvdxnew[i]*cxs_pc + h2dvdynew[i]*cys_pc)


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0s[i].data),
            np.ravel(hi_1s[i].data),
            np.ravel(hi_2s[i].data),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hbperp2),
            np.ravel(hv_model),
            np.ravel(hv_model2)
            )).transpose(),
            columns=['w','v','s','rper','rper2','vmod','vmod2'])


    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]

    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    meanhdf2 = hdf.groupby(pd.cut(hdf['rper2'],
                        bins=np.arange(np.nanmin(hdf['rper2'])-hR/15,
                        np.nanmax(hdf['rper2']),
                        hR/15))).mean().dropna()

    sighdf2 = hdf.groupby(pd.cut(hdf['rper2'],
                        bins=np.arange(np.nanmin(hdf['rper2'])-hR/15,
                        np.nanmax(hdf['rper2']),
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0s[i].data),
            np.ravel(h2_1s[i].data),
            np.ravel(cbperp),
            np.ravel(cbperp2),
            np.ravel(cv_model),
            np.ravel(cv_model2),
            )).transpose(),
            columns=['w','v','rper','rper2','vmod','vmod2'])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    meancdf2 = cdf.groupby(pd.cut(cdf['rper2'],
                        bins=np.arange(np.nanmin(cdf['rper2'])-cR/15,
                        np.nanmax(cdf['rper2']),
                        cR/15))).mean().dropna()


    sigcdf2 = cdf.groupby(pd.cut(cdf['rper2'],
                        bins=np.arange(np.nanmin(cdf['rper2'])-cR/15,
                        np.nanmax(cdf['rper2']),
                        cR/15))).std().dropna()

    
    if show:
        fig,axs = plt.subplots(1,2,figsize=(12,5))
        # fig.suptitle(names_full[i])
    
        cm,cb = curve_fit(line,meancdf2['rper2'],meancdf2['v'])[0]
        csm   = np.sqrt(curve_fit(line,meancdf2['rper2'],meancdf2['v'])[1][0,0])
    
        hm,hb = curve_fit(line,meanhdf2['rper2'],meanhdf2['v'])[0]
        hs = np.sqrt(curve_fit(line,meanhdf2['rper2'],meanhdf2['v'])[1][0,0])
    
        axs[0].fill_between(meancdf2['rper2'],
                            meancdf2['v']+sigcdf2['v'],
                            meancdf2['v']-sigcdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
        
        axs[0].fill_between(meancdf2['rper2'],
                            meancdf2['v']+2*sigcdf2['v'],
                            meancdf2['v']-2*sigcdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
    
     
        axs[0].plot(meancdf2['rper2'],meancdf2['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[0].plot(meancdf2['rper2'],meancdf2['vmod2'],color='red',linewidth=0.7,linestyle='--')
        
    
        axs[1].fill_between(meanhdf2['rper2'],
                            meanhdf2['v']+sighdf2['v'],
                            meanhdf2['v']-sighdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
        
        axs[1].fill_between(meanhdf2['rper2'],
                            meanhdf2['v']+2*sighdf2['v'],
                            meanhdf2['v']-2*sighdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
    
        
        axs[1].plot(meanhdf2['rper2'],meanhdf2['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[1].plot(meanhdf2['rper2'],meanhdf2['vmod2'],color='red',linewidth=0.7,linestyle='--')
    
    
        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
            # ax.grid()
        axs[0].set_ylabel(r'velocity [$\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names_full[i])
    
        # fig.savefig('v_vs_r_'+names[i],bbox_inches='tight')
        plt.show()

    return pearson(meancdf2['rper2'], meancdf2['v']),                   pearson(meanhdf2['rper2'], meanhdf2['v']),\
           pearson(sigcdf2['v'],      np.abs(meancdf2['v']-h2_vos[i])), pearson(sighdf2['v'],      np.abs(meanhdf2['v']-hi_vos[i])),\
           pearson(meancdf2['w'],     np.abs(meancdf2['v']-h2_vos[i])), pearson(meanhdf2['w'],     np.abs(meanhdf2['v']-hi_vos[i]))


ps = np.array([plot_rvsv(i) for i in range(22)]).transpose()


plt.scatter(names,ps[2],zorder=-20,alpha=0.5,linewidth=0.2,s=40, 
            label = r'$\sigma \leftrightarrow \Delta v$',color=green,edgecolor='black')
plt.xticks(rotation='vertical')

plt.show()


fig = plt.figure(dpi=200)#,figsize=(6,6))
ax = fig.add_subplot()

ax.scatter(ps[0],ps[1],zorder=10,linewidth=0.2,s=70,
            label = r'r$_\perp \leftrightarrow \overline{v_\mathrm{LOS}}$',color=magnt,edgecolor='black')
ax.scatter(ps[2],ps[3],zorder=-20,alpha=0.5,linewidth=0.2,s=40, 
            label = r'$\sigma \leftrightarrow \Delta v$',color=green,edgecolor='black')

# plt.scatter(ps[4],ps[5],zorder=-20,alpha=0.5,linewidth=0.2,s=40, 
#             label = r'$\sigma \leftrightarrow \Delta v$',color=magnt,edgecolor='black')

ax.set_xlabel(r'p (molecular clouds)')#$\mathrm{H_2}$)')
ax.set_ylabel('p (HI envelopes)')

ax.vlines(0.5,-1.1,1.5,color=magnt,linewidth=0.7)
ax.hlines(0.5,-1.1,1.1,color=magnt,linewidth=0.7)


lim = 1.05
ax.set_xlim(-lim,lim)
ax.set_ylim(-lim,lim)

ax.set_aspect(1)

ax.grid(zorder=-10,alpha=0.5,color='black',linewidth=0.3)

ax.legend(loc='upper left',fontsize='x-small')


xhist = ax.inset_axes([0, 1.05, 1, 0.25],)

# xhist.hist(ps[0],bins=np.arange(-1.,1.1,0.125),color=magnt,alpha=0.5,edgecolor='black',linewidth=0.5)#,histtype='barstacked')
xhist.hist(ps[2],bins=np.arange(-1.,1.1,0.125),color=green,alpha=0.5,edgecolor='black',linewidth=0.5)

xhist.set_xticks(np.arange(-1.,1.,0.5),[])
xhist.set_xlim(-lim,lim)
xhist.set_yticks([0,2,4,6,8])

# xhist.grid(zorder=-10,alpha=0.5,color='black',linewidth=0.5)


yhist = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
# yhist.hist(ps[1],bins=np.arange(-1.,1.1,0.125),color=magnt,alpha=0.5,orientation='horizontal',edgecolor='black',linewidth=0.5)
yhist.hist(ps[3],bins=np.arange(-1.,1.1,0.125),color=green,alpha=0.5,orientation='horizontal',edgecolor='black',linewidth=0.5)

yhist.set_yticks(np.arange(-1.,1.,0.5),[])
yhist.set_ylim(-lim,lim)
# yhist.grid(zorder=-10,alpha=0.5,color='black',linewidth=0.5)
yhist.set_xticks([0,4,8,12])

# plt.savefig('pvalues',bbox_inches='tight')

plt.show()


#add the names of the clouds with low pearson coefficient


# # patch by patch  / spectra along theta hat


los = los.astype('float')
lis = lis.astype('float')

bos = bos.astype('float')
bis = bis.astype('float')


lcs = ((np.nanmean((los,lis),axis=0) - 360)/-0.125).astype('int')
bcs = ((np.nanmean((bos,bis),axis=0) + 50)/0.125).astype('int')


hi4pi, dht = fits.open('/mnt/ceph/users/tmekuria/dust/both_same_grid_720d.fits')


his = [hi4pi.data[:,
            bcs[i]-rhi_px[i]:bcs[i]+rhi_px[i],
            lcs[i]-rhi_px[i]:lcs[i]+rhi_px[i]]
            for i in range(22)]


coshsz = [dht.data[:,
            bcs[i]-rhi_px[i]:bcs[i]+rhi_px[i],
            lcs[i]-rhi_px[i]:lcs[i]+rhi_px[i]]
            for i in range(22)]


losp,lisp = ((np.array((los,lis))-360)/-0.125).astype('int')
bosp,bisp = ((np.array((los,lis))-360)/-0.125).astype('int')


cos = [dht.data[:,
            bosp[i]:bisp[i],
            losp[i]:lisp[i]]
            for i in range(22)]


vos = vos.astype('float')
vis = vis.astype('float')


plt.scatter(names,np.abs(hidvdynew/hidvdxnew))

plt.xticks(rotation='vertical')

plt.yscale('log')

plt.grid()

plt.show()


def plot_spec(ax,j):
    bo,bi = slcs[j][0]
    lo,li = slcs[j][1]
    sp = np.nanmean(his[i][:,bo:bi,lo:li],axis=(1,2))
    
    ax.plot(np.arange(155)*1.3-100,(sp-np.nanmin(sp))/np.nanmax(sp),linewidth=0.7,color='black')


def plot_pbp(i):
    n = his[i].data.shape[1]
    
    num_ptch = 11
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    
    coords = np.array((inds[:-1],inds[1:])).transpose()
    
    slcs = np.array([(i,j) for j in coords for i in coords])
    

    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                           dpi=200)
    
    map_axs = np.flip(axs,axis=0).transpose().flat
    
    # res = Parallel(n_jobs=-1, verbose=4)(delayed(plot_spec)(ax=map_ax,j=j) for j,map_ax in enumerate(map_axs))
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(his[i][:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(155)*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),
                        linewidth=0.7,color='black')


        map_axs[j].set_xticks([]),
        map_axs[j].set_xlim(vos[i]-10,vis[i]+10),
        map_axs[j].vlines(np.nanmean((vos[i],vis[i])),0,1,linestyle='--',linewidth=0.5,color='black',)
        map_axs[j].set_yticks([])
        # map_axs[j].set_axis_off()
        
        
    bkg_map = fig.add_subplot(111)
    bkg_map.set_axis_off()
    bkg_map.contourf(np.where(~hi_msk[i],np.nan,hi_1s[i].data),levels=30,
                     cmap='Spectral',extend='both',alpha=0.3,antialiased=True)
    
    bkg_map.quiver(rhi_px[i],rhi_px[i],-hidvdxnew[i],hidvdynew[i])
#                      extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
#                              cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
#                     linewidths=0.4,levels=np.arange(-1.3+mom1_vs[i][0],1.3+mom1_vs[i][1],mom1_vs[i][2]/2))
# #     bkg_map.contourf(c1s[i].data,alpha=0.3,levels=10,cmap='coolwarm_r',
# #                     extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
# #                             cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
    

# #     bkg_map.invert_xaxis()
#     bkg_map.set_xlim(cat['lon_HI<deg>'][i])
#     bkg_map.set_ylim(cat['lat_HI<deg>'][i])
# #     bkg_map.set_aspect('equal')
        
    plt.subplots_adjust(wspace=-0.07,hspace=-0.07)
        
    fig.suptitle(names[i])
    # plt.savefig('patch_spec_'+lbls[i])
    plt.show()


plot_pbp(8)
