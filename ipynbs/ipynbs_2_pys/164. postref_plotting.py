#!/usr/bin/env python
# coding: utf-8


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

    
    mskd_mom1 = np.where(~mask,np.nan,mom1_map.data)*1e3
    
    
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
    
    h22 = np.sqrt(second_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(~h2_msks[i],np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = np.linalg.norm(np.indices(hi1.shape) - \
                   np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) 
    
    himaps = np.array([np.where(hi_msk > rhi_px[i],np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


table1 = pd.read_csv('/Users/thummimm/Downloads/molecular_clouds.csv')

coordses = table1[['l0 (°)','l1 (°)',
                   'b0 (°)','b1 (°)']].to_numpy()
#                    

ds = table1['d (pc)'].to_numpy()


hi4pi, dht = fits.open('../fitses/both_same_grid_720d.fits')


# threshs = [20. , 20.7, 20. , 20.5, 20.5, 20.5, 20.5, 20.8, 20.5, 20. , 20.5,
#        20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20. , 20.5, 20.5, 20.7, 21.1]

hi_Nlims = np.array([
    [20.1,  21.03],
    [20.56, 21.72],
    [20.5, 20.85],
    [20.575, 21.35],
    [20.4, 21.],
    [20.1,  20.9],
    [20.8, 21.25],
    [21, 21.64],
    [20.68, 21.6],
    [20.64, 21.22],
    [20.3, 21.78],
    [21.04, 21.56],
    [20.6, 21.56],
    [20.7, 21.35],
    [20.5, 21.68],
    [20., 21.4],
    [20.5, 21.5],
    [20.75, 21.6],
    [20.22, 21.09],
    [19.65, 21.2 ],
    [21.1,  21.58],
    [20.4,  21.64],
])


hi_vlims = np.array([
[-1.,   4.4],
[-8.5,  7.],
[-2.4,  -0.24],
[13.95, 18.6],
[-1.,   4.2],
[3.5, 6.2],
[8.25, 12.6],
[0.8, 7.],
[2.4, 11.7],
[-0.75,  3.],
[-1.8,  7.5],
[17.5, 28.],
[3.5, 10.],
[6.5, 12.3],
[-1.5,   5.],
[4.5,  10.],
[4.5,  11.],
[2,  8.5],
[1.4, 6.4],
[-6.75, 0.3],
[10.05, 14.7],
[-0.8, 10.]
])


h2_Nlims = np.array(([20.0,
 20.7,
 20.0,
 20.7,
 20.5,
 20.3,
 20.5,
 20.8,
 20.5,
 20.0,
 20.5,
 20.5,
 20.8,
 20.5,
 20.5,
 20.3,
 20.3,
 20.3,
 20.7,
 20.3,
 20.7,
 20.7]
    
    ,[21.25, 21.9, 21.3, 21.7, 21.7, 21.7, 21.8, 22., 21.7,
       21.5, 21.8 , 22, 22, 21.8, 22.2, 22.3, 22.3, 22,
       21.9, 21.5, 22, 21.9])).transpose()

h2_vlims = np.array((
    [0.9,  4.95],
    [-10,   8.],
    [-2.5,  0.7],
    [12, 20],
    [0.6,  6.],
    [4.56, 6.75],
    [10, 15.6],
    [0., 8],
    [4., 14.7],
    [-4.75,  4.25],
    [3.,  9.8],
    [20.,  31],
    [0, 10],
    [8, 14],
    [0.5, 6.],
    [2.5, 12.9],
    [3.9, 12],
    [-7,  15. ],
    [-1.,  9.5],
    [-6.25,  0],
    [7.2, 17],
    [4, 8]))


h2_0ths = [make_lb_lv_bv(get_subcube(dht.data,coords))[0] for coords in coordses]


h2_msks = [np.log10(zeroth)>20. for zeroth in h2_0ths]


names = table1['Name'].to_numpy()


lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')


from astropy.modeling import models,fitting
gaus = models.Gaussian1D(1.,5.,5.)
# gaus2 = models.Gaussian1D(1.,5.,5)
fttr = fitting.LMLSQFitter()


def get_rhi(i):
    coords = coordses[i]
    h2cube = get_subcube(dht.data,coords)
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=20.)

    rhi = (((3*mh2*u.solMass)/(4*np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
    return mh2,rh2,rhi


rhi = [get_rhi(i)[2] for i in range(22)]
rhi_px = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')


# def fit_specs(i,coords):
#     loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
#     boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    
#     hicube = hi4pi.data[:, boh:bih, loh:lih]
#     h2cube = get_subcube(dht.data,coords)

#     sph2 = np.nanmean(h2cube,axis=(1,2))
#     sphi = np.nanmean(hicube,axis=(1,2))
    
#     vx = (np.arange(155)*1.3-100)
#     sph2fit = fttr(gaus,vx,sph2)
    
#     Afit = sph2fit.amplitude.value
#     mufit = sph2fit.mean.value
#     stdfit = sph2fit.stddev.value
    
#     return Afit,mufit,stdfit
    
    
# def plot_spec(A,mu,std):
    
#     plt.figure(dpi=200,figsize=(5,5))    
#     mody = models.Gaussian1D(A,mu,std)
#     plt.vlines(mody.mean+np.array([-3,3])*mody.stddev,
#                    0,np.nanmax(sphi)+3,
# #                    alpha=1/(1+j),
#                    linewidth=1.5,color='black',linestyle='--')
    
        
#     plt.plot(vx,sph2*60,linewidth=1.5,color=(0.,0.8,0.0))
# #     axs[0].set_ylabel(r'$\overline{T_\mathrm{B}}\ [\times$'+' {} K]'.format(round(np.nanmax(sph2),3)),)
#     plt.plot(vx,sphi,linewidth=5,color=(0.8,0.,0.))
    
    
# #     plt.xticks(np.arange(-70,70,10))
#     plt.xlim(mufit+np.array([-1,1])*10*stdfit)
#     plt.ylim(0,np.nanmax(sphi)+3,)
# #     plt.vlines(vs[i],0,1)
    
#     plt.title(names[i])
#     plt.xlabel(r'$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
#     plt.ylabel(r'$\overline{T_\mathrm{B}}\ [\mathrm{K}]$')
    
# #     plt.gca().set_aspect(1)
#     plt.savefig(names[i][:3]+names[i][-1]+'_spec_zx')
#     plt.show()


for i,coords in enumerate(coordses):
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    
    hicube = hi4pi.data[:, boh:bih, loh:lih]
    h2cube = get_subcube(dht.data,coords)

    sph2 = np.nanmean(h2cube,axis=(1,2))
    sphi = np.nanmean(hicube,axis=(1,2))
    
    vx = (np.arange(155)*1.3-100)
    sph2fit = fttr(gaus,vx,sph2)
    
    Afit = sph2fit.amplitude.value
    mufit = sph2fit.mean.value
    stdfit = sph2fit.stddev.value
    
    plt.figure(dpi=200,figsize=(5,5))
    
    mody = models.Gaussian1D(Afit,mufit,2*stdfit)
    plt.vlines(mody.mean+np.array([-3,3])*mody.stddev,
                   0,np.nanmax(sphi)+3,
#                    alpha=1/(1+j),
                   linewidth=1.5,color='black',linestyle='--')
    
        
    plt.plot(vx,sph2*60,linewidth=1.5,color=(0.,0.8,0.0))
#     axs[0].set_ylabel(r'$\overline{T_\mathrm{B}}\ [\times$'+' {} K]'.format(round(np.nanmax(sph2),3)),)
    plt.plot(vx,sphi,linewidth=5,color=(0.8,0.,0.))
    
    
#     plt.xticks(np.arange(-70,70,10))
    plt.xlim(mufit+np.array([-1,1])*10*stdfit)
    plt.ylim(0,np.nanmax(sphi)+3,)
#     plt.vlines(vs[i],0,1)
    
    plt.title(names[i])
    plt.xlabel(r'$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
    plt.ylabel(r'$\overline{T_\mathrm{B}}\ [\mathrm{K}]$')
    
#     plt.gca().set_aspect(1)
    plt.savefig(names[i][:3]+names[i][-1]+'_spec_zx')
    plt.show()


#modify this to return vmeanhi,sighi and take in i,coords,vo,vi 
def get_vshi (i):
    coords = coordses[i]
    h2cube = get_subcube(dht.data,coords)
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=20.)
    
    rhi = (((3*mh2*u.solMass)/(4*np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
    rhi_px = (np.arctan(rhi/ds[i])*u.rad.to(u.deg)/0.125).astype('int')
    
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px
    
    hicube = hi4pi.data[:, boh:bih, loh:lih]

    sph2 = np.nanmean(h2cube,axis=(1,2))
    sphi = np.nanmean(hicube,axis=(1,2))
    
    vx = (np.arange(155)*1.3-100)
    sph2fit = fttr(gaus,vx,sph2)
    
    Afit = sph2fit.amplitude.value
    mufit = sph2fit.mean.value
    stdfit = sph2fit.stddev.value
    
    vs1 = sph2fit.mean+np.array([-3,3])*sph2fit.stddev
    vs2 = sph2fit.mean+np.array([-3,3])*2*sph2fit.stddev
    vs3 = sph2fit.mean+np.array([-3,3])*3*sph2fit.stddev
    
    return vs1,vs2,vs3


vses = ((np.array([get_vshi(i) for i in range(22)])+100)/1.3).astype('int')


def get_mhis(coordses):
    mses = []
    for i,coords in enumerate(coordses):
        loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
        boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]

        hicubes = [hi4pi.data[vo:vi, boh:bih, loh:lih] for vo,vi in zip(vses[i][:,0],vses[i][:,1])]

        ms = [calc_mass(hicube,ds[i],thresh=20.,conv=1.82e18)[0] for hicube in hicubes]
        mses.append(ms)
        
    return np.array(mses)


mh2 = [get_rhi(i)[0] for i in range(22)]
rh2 = [get_rhi(i)[1] for i in range(22)]


mhis = get_mhis(coordses)


def get_theta(xy):
    x,y = xy
    thet = np.arctan(y/x)*u.rad.to(u.deg)
    
    if x < 0:
        thet += 180
        
    if thet > 180:
        thet -= 360
        
    return thet
    
    
def fit_plane_to_maps(i,sig=2):
    h2maps, himaps = moment_maps(i,sig)
        
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    h2fitpars = fit_gradient(h21,ds[i],weights=h22)
    hifitpars = fit_gradient(hi1,ds[i],weights=hi2)
        
    return np.array(h2fitpars),np.array(hifitpars)


allFits = [[fit_plane_to_maps(i,j) for i in range(22)] for j in [0,1,2]]


omgh2 = np.linalg.norm(np.array(allFits[0])[:,0,:2],axis=1)
thth2 = [get_theta(np.array(allFits[0])[i,0,:2]) for i in range(22)]
voh2 = np.array(allFits[0])[:,0,2]*1e-3


omghi = np.linalg.norm(np.array(allFits[1])[:,1,:2],axis=1)
ththi = [get_theta(np.array(allFits[1])[i,1,:2]) for i in range(22)]
vohi = np.array(allFits[1])[:,1,2]*1e-3


omghi3 = np.linalg.norm(np.array(allFits[2])[:,1,:2],axis=1)


jhi = (omghi*u.m/u.s/u.pc * (rhi * u.pc)**2).cgs.value


jhi3 = (omghi3*u.m/u.s/u.pc * (rhi * u.pc)**2).cgs.value


jh2 = (omgh2*u.m/u.s/u.pc * (rh2 * u.pc)**2).cgs.value


np.nanmedian(jhi3/jh2)


curve_fit(line,np.log10(rhi),np.log10(jhi))


curve_fit(line,np.log10(rhi),np.log10(jhi3))


np.sqrt(curve_fit(line,np.log10(rhi),np.log10(jhi3))[1][0,0])


jhi3/jhi


np.save('table2_new.npy',np.array((mh2,rh2,omgh2,thth2,voh2,mhis[:,1],rhi,omghi,ththi,vohi)))


plt.figure(dpi=200)
[plt.vlines(mh2[i],mhis[i][0],mhis[i][2],color='red',zorder=10,lw=0.75) for i in range(22)]
[plt.vlines(mh2[i],mhis[i][0],mhis[i][1],color='red',zorder=10,lw=0.75) for i in range(22)]
[plt.vlines(mh2[i],mhis[i][1],mhis[i][2],color='red',zorder=10,lw=0.75) for i in range(22)]

plt.scatter(mh2,mhis[:,1],color='black',zorder=20,s=7)
#     plt.plot(mhi)
plt.loglog()
plt.ylabel('$M_\mathrm{HI} \ [M_\odot]$')
plt.xlabel('$M_\mathrm{H_2} \ [M_\odot]$')


def polar_hist(ax,data,bins=np.arange(0,361,15)*u.deg.to(u.rad),density='False'):
    data = np.asarray(data)
    data = np.mod(data,2*np.pi)
    counts,binsedges = np.histogram(data,bins=bins)#,density=density)
    
    widths = 0.7*np.diff(binsedges)
    bincnts = binsedges[:-1] + widths/2

    
    ax.bar(bincnts,counts, width=widths,bottom=0.,color=(0.,0.,0.8))
    ax.set_yticks([1,2,3,4],[])
    ax.set_xticks(np.arange(0,360,30)*u.deg.to(u.rad))
    ax.grid(color='black')
#     plt.show()


fig,ax = plt.subplots(subplot_kw={'projection':'polar'})
[polar_hist(ax,np.array(ththis[:,i])*u.deg.to(u.rad)) for i in [1]]

# [polar_hist(ax,np.array(ththis[:,i])*u.deg.to(u.rad)) for i in range(3)]


fig,ax = plt.subplots(subplot_kw={'projection':'polar'})

polar_hist(ax,np.array(thth2)*u.deg.to(u.rad))


del_thet = np.mod((ththis[:,1]-thth2)*u.deg.to(u.rad),2*np.pi)


fig,ax = plt.subplots(subplot_kw={'projection':'polar'})

polar_hist(ax,del_thet)


all_maps = [moment_maps(i) for i in range(22)]


h2_1maps = [moment_maps(i,sig=1)[0] for i in range(22)]


hi_1maps = [moment_maps(i,sig=1)[1] for i in range(22)]
hi_3maps = [moment_maps(i,sig=3)[1] for i in range(22)]


fits.HDUList([fits.PrimaryHDU(all_maps[0][0])]+[fits.ImageHDU(maps[0]) for maps in all_maps[1:]]).writeto('h2_maps.fits')


fits.HDUList([fits.PrimaryHDU(all_maps[0][1])]+[fits.ImageHDU(maps[1]) for maps in all_maps[1:]]).writeto('hi_maps.fits')


fits.HDUList([fits.PrimaryHDU(hi_1maps[0])]+[fits.ImageHDU(maps) for maps in hi_1maps[1:]]).writeto('hi_1sigmaps.fits',overwrite=True)


fits.HDUList([fits.PrimaryHDU(hi_3maps[0])]+[fits.ImageHDU(maps) for maps in hi_3maps[1:]]).writeto('hi_3sigmaps.fits',overwrite=True)


fits.HDUList([fits.PrimaryHDU(h2_1maps[0])]+[fits.ImageHDU(maps) for maps in h2_1maps[1:]]).writeto('h2_1sigmaps.fits',overwrite=True)


all_maps


np.array([plot_rvsv(i,coordses[i],show=True) for i in range(22)]).transpose()


ps = np.array([plot_rvsv(i,coords) for i,coords in enumerate(coordses)]).transpose()


plt.figure(dpi=200)
plt.scatter(ps[0],ps[1])

plt.hlines(0.5,-1,1.1,color='black',linestyle='--')
plt.vlines([-0.5,0,0.5],0.4,1.1,linestyles=['--','solid','--'],color='black')

plt.xlabel('p ($r_{\perp},\Delta v)_\mathrm{H_2}$')
plt.ylabel('p ($r_{\perp},\Delta v)_\mathrm{HI}$')

plt.ylim(0.45,1.05)
plt.xlim(-0.55,1.05)

[plt.annotate(names[i][:3]+names[i][-1], (ps[0][i],ps[1][i])) for i in range(22)]

plt.show()


plt.figure(dpi=200)
plt.scatter(ps[0],ps[1])

plt.hlines(0.5,-1,1.1,color='black',linestyle='--')
plt.vlines([-0.5,0,0.5],0.4,1.1,linestyles=['--','solid','--'],color='black')

plt.xlabel('p ($r_{\perp},\Delta v)_\mathrm{H_2}$')
plt.ylabel('p ($r_{\perp},\Delta v)_\mathrm{HI}$')

plt.ylim(0.45,1.05)
plt.xlim(-0.55,1.05)


plt.show()


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


mrks = ['o','s','^']


plt.figure(dpi=200)
[plt.vlines(omgh2[i],omghis[i][0],omghis[i][2],color='red',zorder=10,lw=0.75) for i in range(22)]
[plt.vlines(omgh2[i],omghis[i][0],omghis[i][1],color='red',zorder=10,lw=0.75) for i in range(22)]


[plt.scatter(omgh2,omghis[:,i],zorder=20,color='black',
             marker=mrks[i],alpha=1/(i+1),linewidth=0.) for i in [0,1,2]]
#     plt.plot(mhi)
plt.loglog()
plt.ylabel('$\Omega_\mathrm{HI} \ [\mathrm{km\ s^{-1}\ pc^{-1}}]$')
plt.xlabel('$\Omega_\mathrm{H_2} \ [\mathrm{km\ s^{-1}\ pc^{-1}}]$')


jhis = (omghis*u.km/u.s/u.pc*(np.array(rhi)[:,np.newaxis]*u.pc)**2).cgs.value


jh2 = (omgh2*u.km/u.s/u.pc*(np.array(rh2)*u.pc)**2).cgs.value


np.save('wx_hi_props_moThj.npy',np.array([mhis,omghis,ththis,jhis]))


np.save('h2_props_mRroThj.npy',np.array([mh2,rh2,rhi,omgh2,thth2,jh2]))


clrs = ['red','green','blue']


curve_fit(line,np.log10(rhi),np.log10(jhis[:,1]))


plt.figure(dpi=200)

plt.scatter(rh2,jh2,zorder=20,color='black',marker=mrks[1])
plt.loglog()
plt.ylabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.xlabel('$R_\mathrm{H_2}\ [\mathrm{pc}$]')


# plt.xticks(np.arange(20,200,30))


plt.figure(dpi=200)

plt.scatter(rhi,jhis[:,1],zorder=20,color='black',edgecolor=clrs[1], marker=mrks[1])
plt.loglog()
plt.ylabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.xlabel('$R_\mathrm{HI}\ [\mathrm{pc}$]')


# plt.xticks(np.arange(20,200,30))


plt.figure(dpi=200)
# [plt.vlines(rhi[i],jhis[i][0],jhis[i][2],color='black',zorder=10,lw=0.75,linestyle='--') for i in range(22)]
# [plt.vlines(rhi[i],jhis[i][0],jhis[i][1],color='black',zorder=10,lw=0.75,linestyle='--') for i in range(22)]


[plt.scatter(rhi,jhis[:,i],zorder=20,color='white',edgecolor=clrs[i],
             marker=mrks[i]) for i in [0,1,2]]
#     plt.plot(mhi)
plt.loglog()
# plt.yscale('log')
plt.ylabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.xlabel('$R_\mathrm{HI}\ [\mathrm{pc}$]')


# plt.xticks(np.arange(20,200,30))


def make_zxz(i,save=False):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((vs2[i]+100)/1.3).astype('int')
    
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)
    
    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2lb, extent=coords[:4],
                   levels=np.linspace(h2_Nlims[i][0]+0.1,h2_Nlims[i][1],25),
                   cmap='YlGn_r',extend='max')

        
    im2 = axs[0,1].contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="YlGn_r",
                    extend='both',zorder=-10
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2lb<threshs[i]+0.1,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='bwr_r',extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='bwr_r',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.88,antialiased=True,
                       )  
        
    for ax in axs.transpose()[1]:
        ax.contour(h2lb, extent=coords[:4],
                   levels=np.linspace(threshs[i]+0.1,22,5),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hilb.shape)
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
                          ['',r'log (N / $\mathrm{cm^{-2}}$)','',
                           r'$\overline{v_\mathrm{LSR}}\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


[make_zxz(i) for i in range(22)]


v2o,v2i = ((vs2+100)/1.3).astype('int')
    
    print(hicube[v2o:v2i].shape)


#remaking maps with new vbounds (+/-3sigma of a gausshn w 2x co-linewidth)
#fix the Nlims/vlims for the other clouds & update the appendix figures


hi_Nlims[18] = [20.6,21.3]
h2_Nlims[18] = [20.6,22.]


hi_vlims[18] = [-3.5,6.3]
# h2_vlims[18] = [20.6,22.]


make_zxz(18,save=True)


hi_0th = np.nansum(hi4pi.data,axis=0)*1.3*1.82e18


hi_0th_Msunpc2 = hi_0th*(1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value


np.arange(360,-360,-0.125)[[1420,4340]]


[-coords[0],coords[1],coords[2],coords[3]]


[-coords[0]]+coords[1:4]


longitudes = [r'180$\degree$','','',r'150$\degree$','','',r'120$\degree$','','',\
                    r'90$\degree$','','',r'60$\degree$','','',r'30$\degree$','','',0,\
                    '','',r'-30$\degree$','','',r'-60$\degree$', '','',r'-90$\degree$',\
                    '','',r'-120$\degree$','','',r'-150$\degree$','','',r'-180$\degree$',]

latitudes = ['',r'-40$\degree$','',r'-20$\degree$','',0,'',r'20$\degree$','',r'40$\degree$','']


fig = plt.figure(dpi=200,figsize=(20,10))
ax = fig.add_subplot()#projection='mollweide')

ax.contourf(hi_0th_Msunpc2[:,1420:4340],extent=[182.5,-182.5,-50,50],
             levels=np.linspace(10,150,30),extend='both',cmap='binary_r',norm='asinh',zorder=10)
# plt.colorbar(shrink=0.5,label='$\Sigma_\mathrm{HI} \ [M_\odot \ \mathrm{pc^{-2}}]$')
ax.set_xlim(190,-190)
ax.set_aspect(1)

ax.set_frame_on(False)


[ax.add_patch(Circle((lc[i],bc[i]),rhi_px[i]/8,fill=False,color='red',zorder=20)) for i in range(22)]
[ax.add_patch(Circle((-360+lc[i],bc[i]),rhi_px[i]/8,fill=False,color='red',zorder=20)) for i in [1,21]]


for i in range(22):
    coords = coordses[i]
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])

    ax.contourf(h2lb, extent=coords[:4],
                   levels=np.linspace(threshs[i]+0.1,22,25),
                   cmap='YlGn_r',
                   linewidths=0.75,zorder=30)

    
ax.set_xlim(182,-182)
ax.set_xticks(np.arange(180,-181,-10),
              labels=longitudes)
ax.set_yticks(np.arange(-50,51,10),labels=latitudes)

# ax.grid(zorder=-30,color='black',alpha=0.5)
plt.savefig('fig1_cat',bbox_inches='tight')

plt.show()


for i in range(22):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)
    
    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    
    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    hix,hiy = np.indices(hilb.shape)
    hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
    
    im1 = axs[0].contourf(h2lb, extent=coords[:4],
                   levels=25,
                   cmap='YlGn_r',
                   linewidths=0.75,)
    print(im1.levels[[0,-1]])


def make_h2_maps(i):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)
    
    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    
    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    hix,hiy = np.indices(hilb.shape)
    hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
    
    
    axs[0].contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="YlGn_r",
                    extend='both',zorder=-10
                       )    
    
        
    axs[1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='bwr_r',extend='both',zorder=-10)
    
    
    axs[0].contourf(h2lb, extent=coords[:4],
                   levels=np.linspace(h2_Nlims[i][0]+0.1,h2_Nlims[i][1],25),
                   cmap='YlGn_r',extend='max')
    
    im21 = axs[1].contourf(np.where(h2lb<threshs[i]+0.1,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='bwr_r',extend='both')
    
    
#     print(im21.levels[[0,-1]])
    axs[1].annotate(str(np.around(h2_vlims[i][0],1))+\
                    r'$\rightarrow$'+\
                    str(np.around(h2_vlims[i][1],1)),
                    (lih+rhi_px[i]/17,bih-rhi_px[i]/70),
                    fontsize='small')
    
    axs[0].annotate(str(np.around(h2_Nlims[i][0],1))+\
                    r'$\rightarrow$'+\
                    str(np.around(h2_Nlims[i][1],1)),
                    (lih+rhi_px[i]/17,bih-rhi_px[i]/70),
                    fontsize='small')
    
    axs[0].annotate(names[i],(loh-rhi_px[i]/100,bih-rhi_px[i]/70))
    
    
#     print(im.levels[[0,-1]])
#     plt.title(i)
#     plt.colorbar()
    
    for ax in axs:
        ax.contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="binary_r",
                    extend='both',zorder=-9,alpha=0.88,antialiased=True,
                       )  
        ax.set_aspect(1)
        
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],levels=[rhi_px[i],2*rhi_px[i]],colors='white'
                       )
        
        ax.set_xlim(loh,lih)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
    
        ax.tick_params(axis='both', which='both', length=0)
        
#         ax.label_outer()
        
#         fig.suptitle(names[i])
        fig.tight_layout()

        
#         ax.spines.clear()
        
        
#     axs[2].vlines(20.,0,5000,color='red')
#     axs[2].set_yscale('log')
#     plt.savefig(names[i][:3]+names[i][-1]+'_h2_0th_1st',bbox_inches='tight')
    plt.show()


[make_zxz(i,save=True) for i in range(22)]


def make_pbp_spectra(i,save=False):
    coords = coordses[i]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    hicube = get_subcube(hi4pi.data,[loh,lih,boh,bih])
    h2cube = get_subcube(dht.data,coords[:4])
    
        
#     h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],
#                        np.arange(vo,vi)*1.3-100)[0]
    
    hi1 = first_moment(hicube[vo:vi],
                       np.arange(vo,vi)*1.3-100)[0]
    
    hix,hiy = np.indices(hi1.shape)
    hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
    
#     h21 = h21*(h21>h2_vlims[i][0])*(h21<h2_vlims[i][1])*(h2lb>h2_Nlims[i][0])
    hi1 = hi1*(hi1>hi_vlims[i][0])*(hi1<hi_vlims[i][1])*(hicrc<rhi_px[i])

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
    
    fig.savefig(names[i][:3]+names[i][-1]+'_hipbp',bbox_inches='tight')


def fitline(x,y):
    fit = curve_fit(line,x,y)[0]
    std = np.sqrt(np.diag(curve_fit(line,x,y)[1]))
    
    return fit,std


jh2 = h2o*(u.km/u.s/u.pc).to(1/u.s)*(rh2*u.pc.to(u.cm))**2
jhi = hio*(u.km/u.s/u.pc).to(1/u.s)*(rhi*u.pc.to(u.cm))**2
