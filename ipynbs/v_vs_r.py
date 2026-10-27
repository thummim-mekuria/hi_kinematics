import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from   astropy.io import fits
import astropy.units as u
import astropy.constants as c

from scipy.optimize import curve_fit

from joblib import Parallel, delayed

plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


def moment_maps(i,sig=2):
    coords = coordses[i] 
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')

    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = vses[i][sig-1]#0 for 1sigma, 1 for 2 sigma .. 

    h2lb = np.log10(np.nansum(get_subcube(dht.data,coords[:4]),axis=0)*1.3*1.97e20)
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)

    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    h22 = second_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi2 = second_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    
    h2maps = np.array([np.where(~h2_msks[i],np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = np.linalg.norm(np.indices(hi1.shape) - \
                   np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) 
    
    himaps = np.array([np.where(hi_msk > rhi_px[i],np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps



def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den


def get_rperp_vmod(ccoords, hcoords, omgs, thts):
    clons, clats = ccoords
    hlons, hlats = hcoords
    omgh2,omghi = omgs
    ththi,thth2 = thts

    h2dvdx, h2dvdy = omgh2*np.cos(thth2*u.deg).value, omgh2*np.sin(thth2*u.deg).value
    hidvdx, hidvdy = omghi*np.cos(ththi*u.deg).value, omghi*np.sin(ththi*u.deg).value

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)


    hbperp = (hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)/omghis[i,1]#np.linalg.norm((hidvdx[i],hidvdy[i]))
    cbperp = (h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)/np.linalg.norm((h2dvdx[i],h2dvdy[i]))

    hv_model =  hi_vos[i] + hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc
    cv_model =  h2_vos[i] + h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc


    return cbperp,hbperp,cv_model,hv_model



def rvsv_pearson(ccoords, rhi_px, moment_maps, omgs, thts, show=False):  

    l0c,l1c,b0c,b1c = ccoords
        
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))

    l0h,l1h = lc + rhi_px*0.125, lc - rhi_px*0.125
    b0h,b1h = bc - rhi_px*0.125, bc + rhi_px*0.125

    hlons,hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))
    clons,clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))

    cbperp,hbperp,cv_model,hv_model = get_rperp_vmod([clons,clats],[hlons,hlats],omgs, thts)

    h2_0,h2_1,h2_2, hi_0,hi_1,hi_2 = moment_maps

    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','rper','vmod',])

    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>hi_Nlims[i][0]]
    
    
    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    
    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0),
            np.ravel(h2_1),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','rper','vmod'])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>h2_Nlims[i][0]]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    h2p, hip = pearson(meancdf['rper'], meancdf['v']),pearson(meanhdf['rper'], meanhdf['v'])
    
    
    
    fig,axs = plt.subplots(1,2,figsize=(12,5))
    axs[0].fill_between(meancdf['rper'],
                        meancdf['v']+2*sigcdf['v'],
                        meancdf['v']-2*sigcdf['v'],
                        color=(0.8,0.8,0.8,1),linewidth=0.5)
    
    axs[0].fill_between(meancdf['rper'],
                        meancdf['v']+sigcdf['v'],
                        meancdf['v']-sigcdf['v'],
                        color=(0.5,0.5,0.5,0.8),linewidth=0.5)
    

 
    axs[0].plot(meancdf['rper'],meancdf['vmod'],color='red',linewidth=2,linestyle='--')
    axs[0].plot(meancdf['rper'],meancdf['v'],color='black',marker='o',linewidth=0.5,ms=2)

    
    axs[1].fill_between(meanhdf['rper'],
                        meanhdf['v']+2*sighdf['v'],
                        meanhdf['v']-2*sighdf['v'],
                        color=(0.8,0.8,0.8,1),linewidth=0.5)

    axs[1].fill_between(meanhdf['rper'],
                        meanhdf['v']+sighdf['v'],
                        meanhdf['v']-sighdf['v'],
                        color=(0.5,0.5,0.5,0.8),linewidth=0.5)
    
    

    
    axs[1].plot(meanhdf['rper'],meanhdf['v'],color='black',marker='o',linewidth=0.5,ms=2)
    axs[1].plot(meanhdf['rper'],meanhdf['vmod'],color='red',linewidth=2,linestyle='--')
    
    axs[0].annotate(str(np.around(h2p,2)),(0,np.nanmax(meancdf['v'])+np.nanmax(sigcdf['v'])))
    
    axs[1].annotate(str(np.around(hip,2)),(meanhdf['rper'].iloc[0],
                                           np.nanmax(meanhdf['v'])+np.nanmean(sighdf['v'])))

    
    for ax in axs.flat:
        ax.set_xlabel(r'$\Delta r_\perp$ [pc]')

    axs[0].set_ylabel(r'velocity [$\mathrm{km\ s^{-1}}$]')
    fig.suptitle(names[i])

    fig.savefig('v_vs_r_'+names[i][:3]+names[i][-1],bbox_inches='tight')
    
    if show:
        plt.show()
    
    return h2p,hip