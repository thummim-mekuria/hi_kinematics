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
table1 = pd.read_csv('/Users/thummimm/Downloads/molecular_clouds.csv')

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


def first_moment(cube,arr):
    zeroth = np.nansum(cube,axis=0)
    
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


def moment_maps(i):
    
    h2cube,hicube = make_h2_hi_cubes(i)
    vo,vi = vsnew[i]
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(~h2_msks[i],np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            np.log10(hilb)>20
                
    
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
    h20 = moment_maps(i)[0][0]/u.cm**2
    hi0 = moment_maps(i)[1][0] /u.cm**2   
    
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


def get_mh2_rh2_rhi(i):
    coords = coordses[i]
    print(coords)
    h2cube = get_subcube(dht.data,coords[:4])
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=20.)

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
    
    
def get_theta(xy):
    x,y = xy
    thet = np.arctan(y/x)*u.rad.to(u.deg)
    
    if x < 0:
        thet += 180
        
    if thet > 180:
        thet -= 360
        
    return thet
    

mh2,rh2 = np.array([get_mh2_rh2_rhi(i)[:2] for i in range(22)]).transpose()
rhi     = [get_mh2_rh2_rhi(i)[2] for i in range(22)]
rhi_px  = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')


vmnhi,sigvhi = np.array([spctrum_fit(i)[1:] for i in range(22)]).transpose()


vsnew = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')


def fit_plane2maps(i,j):
    mapses = moment_maps(i)[j]
    m0,m1,m2 = mapses
    
    pars = fit_gradient(m1,ds[i],weights=m2)
    return pars


h2_0ths = [make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0] for coords in coordses]
h2_msks = [np.log10(zeroth)>20. for zeroth in h2_0ths]


omgh2 = [np.linalg.norm(fit_plane2maps(i,0)[:2]) for i in range(22)]
omghi = [np.linalg.norm(fit_plane2maps(i,1)[:2]) for i in range(22)]


jh2 = 0.4*(omgh2*u.m/u.s/u.pc*(rh2*u.pc)**2).cgs
jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs


mhi[:,np.newaxis]


vosh2 = [fit_plane2maps(i,0)[2] for i in range(22)]
voshi = [fit_plane2maps(i,1)[2] for i in range(22)]


h2omgyx = [(fit_plane2maps(i,0)[1]/fit_plane2maps(i,0)[0]) for i in range(22)]
hiomgyx = [(fit_plane2maps(i,1)[1]/fit_plane2maps(i,1)[0]) for i in range(22)]


h2omgx,hiomgx = [[(fit_plane2maps(i,j)[0]) for i in range(22)] for j in [0,1]]


mh2,mhi = np.array([calc_masses(i) for i in range(22)]).transpose()


thth2 = [get_theta(fit_plane2maps(i,0)[:2]) for i in range(22)]
ththi = [get_theta(fit_plane2maps(i,1)[:2]) for i in range(22)]


tab2_new = pd.DataFrame(np.around(np.array((mh2/1e4,rh2,np.array(omgh2)/10,thth2,jh2/1e24,\
                                  mhi/1e5,rhi,np.array(omghi)/10,ththi,jhi/1e24)),2).transpose(),
                        columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


tab2_new.to_csv('tab_N20.csv',index=False)


def get_tdis(dat):
    egrav_hi = 0.6*(c.G * dat['mhi'].to_numpy()*u.solMass / (dat['rhi'].to_numpy()*u.pc)).cgs
    erot_hi = (0.4*dat['jhi'].to_numpy()*u.cm**2/u.s*dat['omghi'].to_numpy()*u.km/u.s/u.pc).cgs
    ekin_hi = (0.5*(sigvhi*u.km/u.s)**2).cgs

    etot_hi = erot_hi - egrav_hi + ekin_hi
    
    tdis = (dat['jhi']/etot_hi.value)#*u.s.to(u.Myr)
    
    jinit = (dat['mh2']*dat['jh2'] + dat['mhi']*dat['jhi'])/(dat['mh2'] + dat['mhi'])

    
    return tdis,jinit


tdis,jinit = get_tdis(tab2_new)


plt.figure(dpi=200)
plt.scatter(tdis*u.s.to(u.Myr),jh2/jinit)
[plt.annotate(names[i][:3]+names[i][-1],
              (tdis[i]*u.s.to(u.Myr),jh2[i].value/jinit[i]),fontsize='x-small') 
 for i in range(22)]
#               
plt.xlabel('$\Delta t_\mathrm{dis}$ [Myr]')
plt.ylabel('$j_\mathrm{H_2}\ /\ j_\mathrm{init}$')

plt.hlines([0.5,1],0,0.13)


get_subcube(hi4pi,coo)


hiomgx


plt.figure(dpi=200)

# plt.scatter(np.abs(hiomgyx),bc)
plt.quiver(np.abs(hiomgyx),bc,np.array(hiomgx)/np.array(omghi),
           np.array(hiomgyx)*np.array(hiomgx)/np.array(omghi),width=5e-3)
# [plt.annotate(names[i][:3]+names[i][-1],
#               (np.abs(ththi[i]),np.abs(bc[i])),fontsize='x-small') 
#  for i in range(22)]

plt.xscale('log')

plt.xlabel('$|\Omega_y\ /\ \Omega_x|$')
plt.ylabel('$b_c\ [\degree]$')
plt.hlines(0,-200,200,zorder=-100)
plt.ylim(-35,35)


plt.figure(dpi=200)
plt.scatter(ththi,bc)
# plt.scatter(np.abs(ththi),np.abs(bc))

# [plt.annotate(names[i][:3]+names[i][-1],
#               (np.abs(ththi[i]),np.abs(bc[i])),fontsize='x-small') 
#  for i in range(22)]

[plt.annotate(names[i][:3]+names[i][-1],
              (ththi[i],bc[i]),fontsize='x-small') 
 for i in range(22)]

#               
plt.xlabel(r'$\theta_\mathrm{HI}\ \cos b_c\ [\degree]$')
plt.ylabel('$b_c\ [\degree]$')

plt.hlines(0,-200,200)
plt.vlines(0,-40,40)
plt.xlim(-190,190)
plt.ylim(-30,30)
# 
plt.show()


plt.scatter(her2,her1)


def make_pbp_spectra(i,save=False):
    coords = coordses[i]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    hicube = get_subcube(hi4pi.data,[loh,lih,boh,bih])
    h2cube = get_subcube(dht.data,coords[:4])
    
    hi1 = moment_maps(i)[1][1]
        
    
    hix,hiy = np.indices(hi1.shape)
    hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
    
#     h21 = h21*(h21>h2_vlims[i][0])*(h21<h2_vlims[i][1])*(h2lb>h2_Nlims[i][0])
#     hi1 = hi1*(hi1>hi_vlims[i][0])*(hi1<hi_vlims[i][1])*(hicrc<rhi_px[i])

#     h2ox,h2oy,h2vo = fit_gradient(h21,ds[i])
    hiox,hioy,hivo = fit_gradient(hi1,ds[i])
    
        
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
    fig.suptitle(names[i]+r' $\Omega_x$ :'+str(np.around(hiox,4))+r' $\Omega_y$ :'+str(np.around(hioy,4)))
    
#     fig.tight_layout()
    if save:
        fig.savefig(names[i][:3]+names[i][-1]+'_hipbp',bbox_inches='tight')


make_pbp_spectra(8)


make_pbp_spectra(0)


here


make_h2_hi_cubes(0)[1].shape


varrs = [np.arange(vsnew[i][0],vsnew[i][1])*1.3-100 for i in range(22)]


nums = [np.nansum(make_h2_hi_cubes(i)[1]*varrs[i][:,np.newaxis,np.newaxis],axis=0) for i in range(22)]


hi0s = [moment_maps(i)[1][0] for i in range(22)]


#something else to try: for each pixel (or maybe a coarser grid) calculate the


nums[0].shape


[(
plt.scatter(hi0s[i],nums[i]),
plt.show())
for i in [17]]


[(
plt.scatter(hi0s[i],nums[i]),
plt.show())
for i in [12]]


[sigv_vs_omg(i) for i in range(22)]


def sigv_vs_omg(i):
    hi0,hi1,hi2 = moment_maps(i)[1]
    vohi = voshi[i]
    
    plt.scatter(np.ravel(hi2),
               np.ravel(np.abs(hi1-vohi/1e3)),alpha=0.7,s=3e-2)
    plt.xlabel('$T_\mathrm{B} \ [\mathrm{km\ s^{-1}}]$')
    plt.title(names[i])
    plt.show()


def sig_vs_omg(i):
    hi0,hi1,hi2 = moment_maps(i)[1]
    vohi = voshi[i]
    
    plt.scatter((np.ravel(hi0*u.cm**-2*c.m_p*1.36).to(u.solMass/u.pc**2).value),
               np.ravel(np.abs(hi1-vohi/1e3)),alpha=0.7,s=3e-2)
    plt.xlabel('$T_\mathrm{B} \ [\mathrm{km\ s^{-1}}]$')
    plt.title(names[i])
    plt.show()


#START HERE FOR 2) on <---yagnaw doc - then send ''revisions_v.10.17'' and ''referee report'' to nia
#actually i think the thing i want to be plotting is a density gradient vs a velocity gradient
#ie do the fitting a plane, but this time on the 0th moment maps
#and the vertical ones in velocity will be corelated with vertical ones in density


def fit_plane2maps0(i,j):
    mapses = moment_maps(i)[j]
    m0,m1,m2 = mapses
    
    m0Sig = m0*c.m_p*1.36
    
    pars = fit_gradient(m0Sig*1e3,ds[i],weights=m2)
    return pars


hithth_Sig[8]


ththi[8]


plt.scatter(ththi,hithth_Sig)


hithth_Sig = [get_theta(fit_plane2maps0(i,1)[:2]) for i in range(22)]


hiomg_Sig = [np.linalg.norm(fit_plane2maps0(i,1)[:2]) for i in range(22)]


hiomgx_Sig = [fit_plane2maps(i,1)[0] for i in range(22)]
hiomgyx_Sig = [(fit_plane2maps0(i,1)[1]/fit_plane2maps0(i,1)[0]) for i in range(22)]


plt.figure(dpi=200)

# plt.scatter(np.abs(hiomgyx),bc)
plt.quiver(np.abs(hiomgyx),bc,
           np.array(hiomgx)/np.array(omghi),
           np.array(hiomgyx)*np.array(hiomgx)/np.array(omghi),width=5e-3)

plt.quiver(np.abs(hiomgyx),bc,
           np.array(hiomgx_Sig),#/np.array(hiomg_Sig),
           np.array(hiomgyx_Sig)*np.array(hiomgx_Sig),#/np.array(hiomg_Sig),
           width=5e-3,color='red')

# [plt.annotate(names[i][:3]+names[i][-1],
#               (np.abs(ththi[i]),np.abs(bc[i])),fontsize='x-small') 
#  for i in range(22)]

plt.xscale('log')

plt.xlabel('$|\Omega_y\ /\ \Omega_x|$')
plt.ylabel('$b_c\ [\degree]$')
plt.hlines(0,-200,200,zorder=-100)
plt.ylim(-35,35)


plt.figure(dpi=200)

# plt.scatter(np.abs(hiomgyx),bc)
plt.quiver(np.abs(hiomgyx),np.array(omghi)/1e3,
           np.array(hiomgx)/np.array(omghi),
           np.array(hiomgyx)*np.array(hiomgx)/np.array(omghi),width=5e-3)

# plt.quiver(np.abs(hiomgyx),bc,
#            np.array(hiomgx_Sig),#/np.array(hiomg_Sig),
#            np.array(hiomgyx_Sig)*np.array(hiomgx_Sig),#/np.array(hiomg_Sig),
#            width=5e-3,color='red')

# [plt.annotate(names[i][:3]+names[i][-1],
#               (np.abs(ththi[i]),np.abs(bc[i])),fontsize='x-small') 
#  for i in range(22)]

plt.xscale('log')
plt.vlines(1,0,0.25,zorder=-10)

plt.xlabel('$|\Omega_y\ /\ \Omega_x|$')
plt.ylabel('$\Omega_\mathrm{HI}\ [\mathrm{km\ s^{-1}\ pc^{-1}}]$')
# plt.hlines(0,-200,200,zorder=-100)
plt.ylim(0,0.23)


coordses


bins = 10**np.arange(-2,2.1,0.2)
plt.hist(h2omgyx,bins=bins,edgecolor='white',alpha=0.5)
plt.hist(hiomgyx,bins=bins,edgecolor='white',alpha=0.5)
plt.xscale('log')


plt.hist([(fit_plane2maps(i,1)[1]/fit_plane2maps(i,1)[0]) for i in range(22)])


plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


mh2,rh2 = np.array([get_mh2_rh2_rhi(i)[:2] for i in range(22)]).transpose()


jhix = 0.4*(np.abs(hiomgx)*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs


plt.hist(jhix/jh2,alpha=0.5,bins=bins)
plt.vlines(np.nanmedian(jhix/jh2),0,7)
plt.hist(jhi/jh2,alpha=0.5,bins=bins)
plt.vlines(np.nanmedian(jhi/jh2),0,7)
plt.xscale('log')


h2_0ths = [make_lb_lv_bv(get_subcube(dht.data,coords))[0] for coords in coordses]
h2_msks = [np.log10(zeroth)>20. for zeroth in h2_0ths]


h2_vlims = np.array(([1,4],[-10,6],[-2,1],[12,25],[1,5],[4.8,6.3],[10,15],[0,10],[5,12],[-4,4],
                     [3,8],[23,28],[3,10],[9,13],[0.5,5],
                     [2,10],[4,11],[-2,10],[1,10],[-6,-2],[9,16],[4,7]))
# hi_vlims = 
hi_Slims =  np.array(([2,20],[6,68],[4,15],[12,32],[2,20],[3,12],[13,29],[25,65],[2,60],[6,21],
                     [2,55],[26,60],[8,72],[9,32],[2,52],
                     [1,30],[6,48],[6,48],[4,22],[1,16],[33,73],[1,53]))


thth2 = [get_theta(fit_plane2maps(i,0)[:2]) for i in range(22)]
ththi = [get_theta(fit_plane2maps(i,1)[:2]) for i in range(22)]


galhi = np.arccos(np.cos(ththi*u.deg)*np.cos(bc*u.deg))*u.rad.to(u.deg)
galh2 = np.arccos(np.cos(thth2*u.deg)*np.cos(bc*u.deg))*u.rad.to(u.deg)


ibb = pd.read_csv('ibb_tab1.csv',sep=',',nrows=33)
ibb_jhi = ((ibb['omghi'].to_numpy()*0.01*u.km/u.s/u.pc)*(70*u.pc)**2).cgs
ibb_jh2 = ((ibb['omgh2'].to_numpy()*0.01*u.km/u.s/u.pc)*(ibb['rh2'].to_numpy()*u.pc)**2).cgs
ibb['jh2'] = ibb_jh2
ibb['jhi'] = ibb_jhi


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


ib


in_ib = [12,13,15,18,20]


ibjh2


((10*u.solMass/u.pc**2)/(c.m_p*1.36)).to(u.cm**-2)


np.nanmean((np.array(jhi)[in_ib]/ibjhi))


np.nanmean((np.array(jh2)[in_ib]/ibjh2))


np.nanmean((np.array(omghi)[in_ib]/ibohi))


np.nanmean((np.array(omgh2)[in_ib]/iboh2))


np.nanmean(np.array(mhi)[in_ib]/ib['mhi'])


np.nanmean(np.array(rh2)[in_ib]/ib['rh2'])


np.nanmean(np.array(rhi)[in_ib]/ib['rhi'])


(jhi/jh2)[in_ib]


ibjhi/ibjh2


ibthi = np.array([-173,-11.2,-102,99.4,1.12])
ibth2 = np.array([47.8,-72.4,-111,-20.2,-132])


np.array(ththi)[in_ib]


ibthi


ibthi-np.array(ththi)[in_ib]


np.array(thth2)[in_ib]


ibth2


ibth2-np.array(thth2)[in_ib]


plt.scatter(jh2,jhi)
plt.scatter(ibb_jh2,ibb_jhi)


plt.loglog()


names


plt.plot(np.arange(22),galh2,marker='.')
plt.plot(np.arange(22),np.abs(thth2),marker='.')

plt.hlines([60,120],0,21)
plt.xticks(np.arange(0,22,3))
plt.grid()


plt.plot(np.arange(22),galhi,marker='.')
plt.plot(np.arange(22),np.abs(ththi),marker='.')

plt.hlines([60,120],0,21)
plt.xticks(np.arange(0,22,3))
plt.grid()


np.nanmean(np.array(galhi)[[0,1,3,5,6,8,11,13,15,16]])


np.tan(ththi*u.deg)[[0,1,3,5,6,8,11,13,15,16]]


h2p = np.array((0.94,0.96,0.64,0.74,0.72,-0.11,0.98,0.88,0.66,0.92,0.68,0.93,0.60,
                0.53,0.85,0.06,0.45,0.83,0.98,0.53,0.71,0.83))
hip = np.array((0.98,0.93,0.70,0.89,0.84,0.47,0.97,0.99,0.97,0.95,0.95,0.98,0.97,
               0.97,0.96,0.82,0.87,0.93,0.97,0.87,0.97,0.98))


bins = np.arange(0,1.1,0.1)


plt.figure(dpi=200,figsize=(5,5))
plt.hist(np.abs(h2p),edgecolor='white',color='black',bins=bins,label='$\mathrm{H_2}$')
plt.hist(np.abs(hip),edgecolor='red',histtype='step',bins=bins,linewidth=2.5,label='HI')
plt.ylabel('N (clouds)')
plt.xlabel('Pearson correlation coefficient')
plt.legend()

plt.savefig('pearson',bbox_inches='tight')


plt.scatter(h2p,hip)


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
