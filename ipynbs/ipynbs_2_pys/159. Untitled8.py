#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import cm

from astropy.io import fits
import astropy.units as u
import astropy.constants as c

from scipy.optimize import curve_fit
from joblib import Parallel, delayed


from astropy.modeling import models,fitting
gaus = models.Gaussian1D(1.,5.,5.)
# gaus = models.Gaussian1D(5.,10.,2.)
fttr = fitting.LMLSQFitter()


hi4pi, dht = fits.open('../../data/hi4piDHT_720deg.fits')
table1 = pd.read_csv('../../data/molecular_clouds.csv')
# table1 = pd.read_csv('../../data/zucker_clouds.csv')

coordses = table1[['l0 (°)','l1 (°)',
                   'b0 (°)','b1 (°)',
                  'v0 (km/s)','v1 (km/s)']].to_numpy().astype('float')               

ds    = table1['d (pc)'].to_numpy().astype('float')
names = table1['Name']
vs    = table1[['v0 (km/s)','v1 (km/s)']].to_numpy().astype('float')
vs_px = ((vs+100)/1.3).astype('int') + np.array([-15,15])[np.newaxis,:]
vs_px[3]  = ((np.array([9.2,30.])+100)/1.3).astype('int')
vs_px[11] = ((np.array([20.9,39.1])+100)/1.3).astype('int')


lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')


def get_subcube(bigcube, coords, dx=-0.125, dy=0.125, x0=360, y0=-50):
    lo,li,bo,bi = coords # assume in degrees
    xo,xi = ((np.array([lo,li]) - x0)/dx).astype('int')
    yo,yi = ((np.array([bo,bi]) - y0)/dy).astype('int')
    
    return bigcube[:,yo:yi,xo:xi]


def first_moment(cube,arr,oned=False):
    zeroth = np.nansum(cube,axis=0)
    
    
    if oned:
        num = np.nansum(cube * arr)
    else:
        num = np.nansum(cube * arr[:,np.newaxis,np.newaxis],axis=0)
        
    return num/zeroth, np.log10(zeroth*1.3*1.97e20)

def second_moment(cube,arr):
    zeroth = np.nansum(cube,axis=0)
    first = first_moment(cube,arr)[0]
    
    num = np.nansum(cube * (arr[:,np.newaxis,np.newaxis] - first[np.newaxis,:])**2, axis=0)
    return num/zeroth, first


def make_lb_lv_bv(cube, show=False):
    lb = np.nansum(cube,axis=0)*1.3*1.97e20
    lv = np.nansum(cube,axis=1)*0.125*1.97e20
    bv = np.nansum(cube,axis=2)*0.125*1.97e20
    
    if show:
        fig,axs = plt.subplots(1,3,dpi=200)
        axs[0].contourf(np.log10(lb),levels=np.linspace(20,22,30))
        axs[1].contourf(np.log10(lv),levels=np.linspace(20,22,30))
        axs[2].contourf(np.log10(bv),levels=np.linspace(20,22,30))

        for ax in axs:
            ax.set_aspect(1)
    #     for ax in axs[1:]:   
    #         ax.set_ylim(60,100)

        plt.show()
    return lb,lv#*(1/u.cm**2 *1.36*c.m_p).to(u.solMass/u.pc**2)#,lv,bv


def moment_maps(i,thresh):
    
    h2cube,hicube = make_h2_hi_cubes(i)
    vo,vi = vsnew[i]
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            np.log10(hilb)>thresh
                
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps

def calc_mass(cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    zeroth = np.nansum(cube,axis=0)*1.3*conv/u.cm**2
    
    mask = np.log10(zeroth.value)>thresh
    
    
#     plt.contourf(np.log10(zeroth.value)*mask,levels=np.linspace(threshs[i],22.5,30))
#     plt.gca().set_aspect(1)
#     plt.show()
    
    masked = zeroth*mask
    dA = (d*u.pc*0.125*u.deg.to(u.rad))**2
    area = np.nansum(masked/masked)*dA
    
    mass = np.nansum(masked*dA)*c.m_p*1.36*2
    return mass.to(u.solMass).value, np.sqrt(area/np.pi).value


def calc_masses(i): #cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    h20 = moment_maps(i,threshs[i])[0][0]/u.cm**2
    hi0 = moment_maps(i,threshs[i])[1][0] /u.cm**2   
    
    dA = (ds[i]*u.pc*0.125*u.deg.to(u.rad))**2
    h2m = (np.nansum(h20*dA)*c.m_p*1.36*2).to(u.solMass)
    him = (np.nansum(hi0*dA)*c.m_p*1.36*2).to(u.solMass)
    
    return h2m.value,him.value
    
    
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

    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map.data)*1e3 #m/s/pc instead of km
    
    
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
    
    
def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


def get_mh2_rh2_rhi(i,thresh):
    coords = coordses[i]
#     print(coords)
    h2cube = get_subcube(dht.data,coords[:4])
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=thresh)

    rhi = (((3*mh2*u.solMass)/(4*np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
    return mh2,rh2,rhi


def make_h2_hi_cubes(i):
    coords = coordses[i]

    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vsnew[i]

    hicube = hi4pi.data[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords[:4]) #velocity axis length = 155
    
    return h2cube,hicube


def spctrum_fit(i,which='hi'):
    coords = coordses[i]
    
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vs_px[i]
    
    hicube = hi4pi.data[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords[:4]) #velocity axis length = 155
    
    
#     sph2 = np.nanmean(h2cube,axis=(1,2)) #as things stand rn, we dont care abthe h2 fit
    if which=='hi':
        spctr = np.nanmean(hicube,axis=(1,2))
        vx = (np.arange(vo,vi)*1.3-100)
        
    else:
        spctr = np.nanmean(h2cube,axis=(1,2))
        vx = (np.arange(155)*1.3-100)
    
    fit_sp = fttr(gaus,vx,spctr)
    
    A   = fit_sp.amplitude.value
    mu  = fit_sp.mean.value
    std = fit_sp.stddev.value
    
    return A,mu,std
    
    
def fit_plane2maps(i,j):
    m0,m1,m2 = moment_maps(i,threshs[i])[j]
    
    pars = fit_gradient(m1,ds[i],weights=m2)
    
    return pars


def get_theta(xy):
    x,y = xy
    thet = np.arctan(y/x)*u.rad.to(u.deg)
    
    if x < 0:
        thet += 180
        
    if thet > 180:
        thet -= 360
        
    return thet


threshs = [20.5,21,20.5,21,20.5,20.5,20.5,21,21,20.5,21,21,21,21,21,20.5,21,20.5,20.5,20.5,21,20.5]


mh2,rh2,rhi = np.array([get_mh2_rh2_rhi(i,threshs[i]) for i in range(22)]).transpose()
rhi_px  = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')

vmnhi,sigvhi = np.array([spctrum_fit(i)[1:] for i in range(22)]).transpose()
vsnew = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')


h2_0ths = [make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0] for coords in coordses]
h2_msks = [np.log10(zeroth)>threshs[i] for i,zeroth in enumerate(h2_0ths)]


omgh2 = [np.linalg.norm(fit_plane2maps(i,0)[:2]) for i in range(22)]
omghi = [np.linalg.norm(fit_plane2maps(i,1)[:2]) for i in range(22)]


jh2 = 0.4*(omgh2*u.m/u.s/u.pc*(rh2*u.pc)**2).cgs
jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs


mhi = np.array([calc_masses(i)[1] for i in range(22)])


thth2 = [get_theta(fit_plane2maps(i,0)[:2]) for i in range(22)]
ththi = [get_theta(fit_plane2maps(i,1)[:2]) for i in range(22)]


tab2_newnew = pd.DataFrame(np.around(np.array((mh2/1e4,rh2,np.array(omgh2)/10,thth2,jh2/1e24,\
                                  mhi/1e5,rhi,np.array(omghi)/10,ththi,jhi/1e24)),3).transpose(),
                        columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


# tab2_newnew = pd.DataFrame(np.array((mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi)).transpose(),
#                         columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


tab2_newnew


# ibjh2 = (np.array((3.92,25.5,42.2,24.2,25.8))*u.pc*u.km/u.s).cgs #these are with c=0.5
# ibjhi = (np.array((25.2,120,107,60.3,73.4))*u.pc*u.km/u.s).cgs
ibrh2 = np.array([14.6,30.6,21.9,16.2,26.5,])
ibrhi = np.array([51.9,83.6,76.0,42.4,78.8])
iboh2 = np.array([4.6,6.68,22.,23.1,9.18])
ibohi = np.array([1.87,3.48,3.71,6.7,2.36])
ibjh2 = (0.4*(iboh2*1e-2*u.km/u.s/u.pc)*(ibrh2*u.pc)**2).cgs
ibjhi = (0.4*(ibohi*1e-2*u.km/u.s/u.pc)*(ibrhi*u.pc)**2).cgs

ibmh2 = np.array([2.5,3.0,2.4,1.4,4.0])*1e4*u.solMass

ib = pd.DataFrame(np.array((ibmh2,ibrh2,iboh2,ibjh2,ibrhi,ibohi,ibjhi)).transpose(),
                  columns='mh2,rh2,omgh2,jh2,rhi,omghi,jhi'.split(','))


in_ib = [12,13,15,18,20]


tab2_newnew.iloc[in_ib]


ib


np.nanmean((tab2_newnew['jhi']/tab2_newnew['jh2']).iloc[in_ib])


def make_pbp_spectra(i,save=False):
    coords = coordses[i]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    hicube = get_subcube(hi4pi.data,[loh,lih,boh,bih])
#     h2cube = get_subcube(dht.data,coords[:4])
    
    
    h2maps, himaps = moment_maps(i)
    
    hi1 = himaps[1]
    hiox, hioy, hivo = fit_plane2maps(i,1)
        
    n = hicube.shape[1]
    num_ptch = 15
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    coords = np.array((inds[:-1],inds[1:])).transpose()
    slcs = np.array([(i,j) for j in coords for i in coords])
    

    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                           dpi=400,figsize=(20,20))
    
    map_axs = np.flip(axs,axis=0).transpose().flat
    
    for j,ax in enumerate(map_axs):
        lo,li = slcs[j][1]
        bo,bi = slcs[j][0]
        
        sp = np.nanmean(hicube[vo:vi,bo:bi,lo:li],axis=(1,2))
        ax.plot(np.arange(vo,vi)*1.3-100,(sp-np.nanmin(sp)),linewidth=1.3,color='black')
        ax.set_xticks([]),
#     #     ax.set_xlim(vos[i]-10,vis[i]+10),
#         ax.vlines(hivo,0,np.nanmax(sp),linestyle='--',linewidth=0.5,color='black',)
        ax.set_yticks([])
    
    
#         cnt = int(num_ptch**2/2)
#         lo,li = slcs[cnt][1]
#         bo,bi = slcs[cnt][0]
#         sp = np.nanmean(hicube[vo:vi,bo:bi,lo:li],axis=(1,2))
#         ax.plot(np.arange(vo,vi)*1.3-100,(sp-np.nanmin(sp)),linewidth=0.5,color='green')
    
            
    plt.subplots_adjust(wspace=-0.07,hspace=-0.07)
    fig.suptitle(names[i]+r' $\Omega_x$ :'+str(np.around(hiox,2))+\
                 r' $\Omega_y$ :'+str(np.around(hioy,2))+') $10^{-2}\ \mathrm{km\ s^{-1}\ pc^{-1}}$')
    
#     fig.tight_layout()
    
#     fig.savefig(names[i][:3]+names[i][-1]+'_hipbp',bbox_inches='tight')


make_pbp_spectra(1)


plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


rhi_px


[ for i in [1]


tab_N20 = pd.read_csv('tab_N20.csv')


tab2_newnew.to_csv('tab_N21.csv',index=False)


plt.plot((tab2_newnew['jhi']/tab2_newnew['jh2'])/(tab_N20['jhi']/tab_N20['jh2'])
    ,marker='.')
plt.hlines(1,0,22)
plt.yscale('log')


def plot_spec(ax,j):
    bo,bi = slcs[j][0]
    lo,li = slcs[j][1]
    sp = np.nanmean(his[i][:,bo:bi,lo:li],axis=(1,2))
    
    ax.plot(np.arange(155)*1.3-100,(sp-np.nanmin(sp))/np.nanmax(sp),linewidth=0.7,color='black')


def plot_pbp(i):    
    h2cube, hicube = make_h2_hi_cubes(i)
#     print(hicube.shape)
    h2maps, himaps = moment_maps(i,threshs[i])
    
    hi1 = himaps[1]
    hiox, hioy, hivo = fit_plane2maps(i,1)
        
    n = hicube.shape[1]
        
    num_ptch = 11
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    
    coords = np.array((inds[:-1],inds[1:])).transpose()
    
    slcs = np.array([(i,j) for j in coords for i in coords])
    

    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                           dpi=200,figsize=(12,12))
    
    bkg_map = fig.add_subplot(111)
    bkg_map.set_axis_off()
    bkg_map.contourf(hi1,levels=30,cmap='Spectral',extend='both',alpha=0.3,antialiased=True,zorder=-30)
    
    bkg_map.quiver(rhi_px[i],rhi_px[i],-hiox,hioy,width=5e-3,color='black',zorder=-20,pivot='middle')
    
    bkg_map.set_xticks(np.array([0,1,2])*rhi_px[i],
                       labels = (lc_p[i]+np.array([-1,0,1])*rhi_px[i])*-0.125+360)
    
    bkg_map.set_yticks(np.array([0,1,2])*rhi_px[i],
                       labels = (bc_p[i]+np.array([-1,0,1])*rhi_px[i])*0.125-50)
    
    
    map_axs = np.flip(axs,axis=0).transpose().flat
        
    # res = Parallel(n_jobs=-1, verbose=4)(delayed(plot_spec)(ax=map_ax,j=j) for j,map_ax in enumerate(map_axs))
    
    mx = np.nanmax(np.nanmean(hicube,axis=(1,2)))
    vo,vi = vsnew[i]
    vax = np.arange(vo,vi)*1.3-100
    
    for ax,j in zip(map_axs,range(len(slcs))):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(hicube[:,bo:bi,lo:li],axis=(1,2))
        
        ax.plot(vax,(sp-np.nanmin(sp))/mx,
                        linewidth=0.7,color='black')

        ax.vlines(vmnhi[i],-0.5,2.5,linestyle='--',linewidth=0.5,color='black',)
#         print(sp.shape,vax.shape)
        ax.scatter(first_moment(sp,vax,oned=True)[0],0.15,color='black',s=4.)
        
        
#     for ax in map_axs:
#         ax.set_xticks([]),
        ax.set_xlim((np.array(vsnew[i])+np.array([-3,3]))*1.3-100, )
        ax.set_ylim(-0.15,2.3)
        
        ax.set_xticks(np.linspace(vo,vi,3)*1.3-100,
        labels = np.around(np.linspace(vo,vi,3)*1.3-100,1),
                      fontsize=5)#'xx-small')
        ax.set_yticks([0,1,2],labels = [0,round(mx,1),round(2*mx,1)],fontsize=5)#'xx-small')
        
#         ax.tick_params(labelsize=5)
        # map_axs[j].set_axis_off()
        
        
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
#     plt.subplots_adjust(wspace=-1,hspace=-1)

        
    fig.suptitle(names[i])
#     fig.tight_layout()
    plt.savefig('pbp_'+names[i][:3]+names[i][-1],bbox_inches='tight')
    plt.show()


[plot_pbp(i) for i in range(22)]


[(
plt.contourf(np.log10(h0),cmap='magma',levels=30,extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
    plt.xlim(coordses[i][:2]),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [h2maps]]


[(
plt.contourf(np.log10(h0),cmap='magma',levels=30,extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
    plt.xlim(coordses[i][:2]),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [h2maps]]


[(
plt.contourf(np.log10(h0),cmap='magma',levels=30,extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
    plt.xlim(coordses[i][:2]),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [h2maps]]


[(
plt.contourf(h1,cmap='Spectral',levels=30,extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
        plt.xlim(coordses[i][:2]),

plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [h2maps]]


[(
plt.contourf(h2,cmap='YlGnBu_r',levels=30),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [h2maps]]


hils = np.around((lc_p[:,np.newaxis]+np.array([-1,1])*rhi_px[:,np.newaxis])*-0.125+360,1)
hibs = np.around((bc_p[:,np.newaxis]+np.array([-1,1])*rhi_px[:,np.newaxis])*0.125-50,1)
hicoordses = np.stack((hils,hibs),axis=1).reshape(22,-1)


[(
plt.contourf(np.log10(h0),cmap='magma',levels=30,
             extent=hicoordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.xlim(hicoordses[i][:2]),
    plt.colorbar(),#np.around(np.linspace(threshs[i],np.nanmax(np.log10(h0)),3),1)),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [himaps]]


[(
plt.contourf(h1,cmap='Spectral',levels=30,extent=hicoordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
    plt.xlim(hicoordses[i][:2]),

plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [himaps]]


[(
plt.contourf(h2,cmap='YlGnBu_r',levels=30,extent=hicoordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i)] for h0,h1,h2 in [himaps]]


def make_zxz(i,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i)
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=35,vmin=5,#np.linspace(15,150,35),
                   cmap='magma',extend='max')

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',)#extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(5,20,50),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.9
                       )  
        
    for ax in axs.flat:#.transpose()[1]:
        ax.contour(h2Sig.value, extent=coords[:4],
                   levels=[1,5,10,30],#np.linspace(15,30,4),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hi0.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],
                    levels=[rhi_px[i],2*rhi_px[i]],colors='white')
        
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2,im3,im4],
                           axs.flat,
                          ['',r'$\Sigma \ [\mathrm{M_\odot\ pc^{-2}}]$','',
                           r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()
