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


table1 = pd.read_csv('/Users/thummimm/Downloads/molecular_clouds.csv')

coordses = table1[['l0 (°)','l1 (°)',
                   'b0 (°)','b1 (°)']].to_numpy().astype('float')               

ds = table1['d (pc)'].to_numpy().astype('float')
names = table1['Name']
vs = table1[['v0 (km/s)','v1 (km/s)']].to_numpy().astype('float')


hi4pi, dht = fits.open('../../../data/hi4piDHT_720deg.fits')


table1


h2_0ths = [make_lb_lv_bv(get_subcube(dht.data,coords))[0] for coords in coordses]
h2_msks = [np.log10(zeroth)>20. for zeroth in h2_0ths]


lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')


def get_mh2_rh2_rhi(i):
    coords = coordses[i]
    h2cube = get_subcube(dht.data,coords)
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=20.)

    rhi = (((3*mh2*u.solMass)/(4*np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
    return mh2,rh2,rhi


rhi = [get_mh2_rh2_rhi(i)[2] for i in range(22)]
rhi_px = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')


mh2,rh2 = np.array([get_mh2_rh2_rhi(i)[:2] for i in range(22)]).transpose()


def spctrum_fit(i,which='hi'):
    coords = coordses[i]
    
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vs_px[i]
    
    hicube = hi4pi.data[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords) #velocity axis length = 155
    
    
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
    
    
plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


vs_px = ((vs+100)/1.3).astype('int') + np.array([-15,15])[np.newaxis,:]


vs_px[3] = ((np.array([9.2,30.])+100)/1.3).astype('int')
vs_px[11] = ((np.array([20.9,39.1])+100)/1.3).astype('int')


def draw_spectra(i):
    coords = coordses[i]
    
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vs_px[i]
    
    hicube = hi4pi.data[:, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords) #velocity axis length = 155
    
    
    hi_spctr = np.nanmean(hicube[vo:vi],axis=(1,2))
    hisp_long = np.nanmean(hicube,axis=(1,2))
    h2_spctr = np.nanmean(h2cube,axis=(1,2))

    Ahi,muhi,stdhi = spctrum_fit(i)
    mody = models.Gaussian1D(Ahi,muhi,stdhi)
    
    
    fig,ax = plt.subplots(dpi=200,figsize=(5,5))
    
    ax.plot(np.arange(155)*1.3-100,hisp_long,linewidth=2,color='red')
#     ax.plot(np.arange(vo,vi)*1.3-100,mody(np.arange(vo,vi)*1.3-100))
    
    ax.vlines(muhi+np.array([-1,1])*2*stdhi,0,np.nanmax(hi_spctr)*1.5,linestyle='--',
             color='black',linewidth=0.8)
    
    ax.set_xlim(muhi+np.array([-1,1])*stdhi*4)
    ax2 = ax.twinx()
    ax2.plot(np.arange(155)*1.3-100,h2_spctr,color='black',linewidth=1)
    ax2.set_ylim(0,1.5*np.nanmax(h2_spctr),)
    
    ax.set_ylim(0,1.3*np.nanmax(hi_spctr))
    
    plt.title(names[i])
    plt.xlabel('$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
    ax.set_ylabel('$T_\mathrm{B,21}\ [K]$',color='red')
    ax2.set_ylabel('$T_\mathrm{B,CO}\ [K]$',)
    ax.set_xlabel('$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
#     plt.show()
    
#     plt.plot(np.arange(155)*1.3-100,h2_spctr)


[(draw_spectra(i),
 plt.savefig('spec_'+names[i][:3]+names[i][-1],bbox_inches='tight')
 
 )for i in range(22)]


[(draw_spectra(i),
 plt.savefig('spec_'+names[i][:3]+names[i][-1],bbox_inches='tight')
 
 )for i in [11]]


[draw_spectra(i) for i in range(22)]


# #write a function that calls draw_spectra over a 3x4 grid

# #mentally prepare to change the scaling factor so it is custom for each cloud, but for now we are going with 100x


# fig,axs=plt.subplots(4,3,dpi=200,figsize=(10,10))#,sharey=True)

# [(
#     axs.flat[i-12].plot(np.arange(155)*1.3-100,np.nanmean(hihi[i].data,axis=(1,2)),linestyle=':',color='black',alpha=0.8,linewidth=3),
#     axs.flat[i-12].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),color='black',linewidth=4),
    
#     axs.flat[i-12].plot(np.arange(155)*1.3-100,100*np.nanmean(coco[i].data,axis=(1,2)),linestyle=':',color=fns.segmented_spectral(20)),
#     axs.flat[i-12].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,100*np.nanmean(coco[i].data[cat['vo'][i]:cat['vi'][i],:,:],axis=(1,2)),color=fns.segmented_spectral(20)),
    
#     axs.flat[i-12].set_xlim((cat['vo'][i]-20)*1.3-100,(cat['vi'][i]+15)*1.3-100),
#     axs.flat[i-12].set_title(names[i],y=0.8,x=0.23,fontsize='small')
    
    
# ) for i in range(12,22)]


# for ax in axs.flat[-2:]:
#     ax.set_axis_off()
    
# fig.supxlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')
# fig.supylabel('$T\mathrm{_B\ [K]}$')
    
# fig.tight_layout()
# plt.savefig('postage_spectra_2',bbox_inches='tight')


vmnhi,sigvhi = np.array([spctrum_fit(i)[1:] for i in range(22)]).transpose()


vsnew = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')


def make_h2_hi_cubes(i):
    coords = coordses[i]

    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vsnew[i]

    hicube = hi4pi.data[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords) #velocity axis length = 155
    
    return h2cube,hicube


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


h2_vlims = np.array(([1,4],[-10,6],[-2,1],[12,25],[1,5],[4.8,6.3],[10,15],[0,10],[5,12],[-4,4],
                     [3,8],[23,28],[3,10],[9,13],[0.5,5],
                     [2,10],[4,11],[-2,10],[1,10],[-6,-2],[9,16],[4,7]))
# hi_vlims = 
hi_Slims =  np.array(([2,20],[6,68],[4,15],[12,32],[2,20],[3,12],[13,29],[25,65],[2,60],[6,21],
                     [2,55],[26,60],[8,72],[9,32],[2,52],
                     [1,30],[6,48],[6,48],[4,22],[1,16],[33,73],[1,53]))


[make_zxz(i,save=True) for i in range(22)]


tab2_new


np.around(tab2_new['jh2']/1e24,2)


sigvhi


make_zxz(8,save=True)


make_zxz(17,save=True)


hi_0th = np.nansum(hi4pi.data,axis=0)*1.3*1.82e18
hi_0th_Msunpc2 = hi_0th*(1.36*c.m_p/u.cm**2).to(u.solMass/u.pc**2).value


longitudes = [r'180$\degree$','','',r'150$\degree$','','',r'120$\degree$','','',\
                    r'90$\degree$','','',r'60$\degree$','','',r'30$\degree$','','',0,\
                    '','',r'-30$\degree$','','',r'-60$\degree$', '','',r'-90$\degree$',\
                    '','',r'-120$\degree$','','',r'-150$\degree$','','',r'-180$\degree$',]

latitudes = ['',r'-40$\degree$','',r'-20$\degree$','',0,'',r'20$\degree$','',r'40$\degree$','']


# fig = plt.figure(dpi=200,figsize=(20,10))
# ax = fig.add_subplot()#projection='mollweide')

# ax.contourf(hi_0th_Msunpc2[:,1420:4340],extent=[182.5,-182.5,-50,50],
#              levels=np.linspace(10,80,30),extend='both',cmap='binary_r',norm='asinh',zorder=10)
# # plt.colorbar(shrink=0.5,label='$\Sigma_\mathrm{HI} \ [M_\odot \ \mathrm{pc^{-2}}]$')
# ax.set_xlim(190,-190)
# ax.set_aspect(1)

# ax.set_frame_on(False)


# [ax.add_patch(Circle((lc[i],bc[i]),rhi_px[i]/8,fill=False,color='red',zorder=20)) for i in range(22)]
# [ax.add_patch(Circle((-360+lc[i],bc[i]),rhi_px[i]/8,fill=False,color='red',zorder=20)) for i in [1,21]]


# for i in range(22):
#     coords = coordses[i]
#     h2lb = ((moment_maps(i)[0][0]/u.cm**2)*c.m_p).to(u.solMass/u.pc**2).value

#     ax.contourf(h2lb, extent=coords[:4],
#                    levels=np.linspace(5,50,25),
#                    cmap='magma',norm='asinh',
#                    linewidths=0.75,zorder=30)

    
# ax.set_xlim(182,-182)
# ax.set_xticks(np.arange(180,-181,-10),
#               labels=longitudes)
# ax.set_yticks(np.arange(-50,51,10),labels=latitudes)

# # ax.grid(zorder=-30,color='black',alpha=0.5)
# # plt.savefig('fig1_cat',bbox_inches='tight')

# plt.show()


def calc_mass(i): #cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    h20 = moment_maps(i)[0][0]/u.cm**2
    hi0 = moment_maps(i)[1][0] /u.cm**2   
    
    dA = (ds[i]*u.pc*0.125*u.deg.to(u.rad))**2
    h2m = (np.nansum(h20*dA)*c.m_p*1.36*2).to(u.solMass)
    him = (np.nansum(hi0*dA)*c.m_p*1.36*2).to(u.solMass)
    
    return h2m.value,him.value


mh2,mhi = np.array([calc_mass(i) for i in range(22)]).transpose()


mh2


mhi


def fit_plane2maps(i,j):
    mapses = moment_maps(i)[j]
    m0,m1,m2 = mapses
    
    pars = fit_gradient(m1,ds[i],weights=m2)
    return pars


omgh2 = [np.linalg.norm(fit_plane2maps(i,0)[:2]) for i in range(22)]
omghi = [np.linalg.norm(fit_plane2maps(i,1)[:2]) for i in range(22)]


omgh2


omghi


jh2 = 0.4*(omgh2*u.m/u.s/u.pc*(rh2*u.pc)**2).cgs
jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs


jh2.value


jhi.value


def get_theta(xy):
    x,y = xy
    thet = np.arctan(y/x)*u.rad.to(u.deg)
    
    if x < 0:
        thet += 180
        
    if thet > 180:
        thet -= 360
        
    return thet


thth2 = [get_theta(fit_plane2maps(i,0)[:2]) for i in range(22)]
ththi = [get_theta(fit_plane2maps(i,1)[:2]) for i in range(22)]


thth2


ththi


ibb = pd.read_csv('ibb_tab1.csv',sep=',',nrows=33)
ibb_jhi = ((ibb['omghi'].to_numpy()*0.01*u.km/u.s/u.pc)*(70*u.pc)**2).cgs
ibb_jh2 = ((ibb['omgh2'].to_numpy()*0.01*u.km/u.s/u.pc)*(ibb['rh2'].to_numpy()*u.pc)**2).cgs
ibb['jh2'] = ibb_jh2
ibb['jhi'] = ibb_jhi


# ibjh2 = (np.array((3.92,25.5,42.2,24.2,25.8))*u.pc*u.km/u.s).cgs
# ibjhi = (np.array((25.2,120,107,60.3,73.4))*u.pc*u.km/u.s).cgs
# ibrh2 = np.array([14.6,30.6,21.9,16.2,26.5,])
# ibrhi = np.array([51.9,83.6,76.0,42.4,78.8])
# iboh2 = np.array([4.6,6.68,22.,23.1,9.18])/100
# ibohi = np.array([1.87,3.48,3.71,6.7,2.36])/100
# ibjh2 = 0.4*((iboh2*u.km/u.s/u.pc)*(ibrh2*u.pc)**2).cgs
# ibjhi = 0.4*((ibohi*u.km/u.s/u.pc)*(ibrhi*u.pc)**2).cgs

# ibmh2 = np.array([2.5,3.0,2.4,1.4,4.0])*1e4*u.solMass


np.nanmean(ibb_jhi/ibb_jh2)


tab2_new = pd.DataFrame(np.array((mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi)).transpose(),
                        columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


plt.figure(dpi=200,figsize=(5,5))
bns=10**np.arange(21,27,1/4)
plt.hist(ibb['jh2'],bins=bns,color='green',histtype='step',
         zorder=20,linewidth=1.5)
plt.hist(tab2_new['jh2'],bins=bns,color='black',zorder=15)

plt.vlines(np.nanmedian(ibb['jh2']),0,8.5,linestyle='--',color='green',alpha=0.5)
plt.vlines(np.nanmedian(tab2_new['jh2']),0,8.5,color='black',alpha=0.5,linestyle='--')

plt.xscale('log')


# [plt.annotate(wch+' '+str(round(np.nanmean(jhi),3)),(2e24,4-0.25*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(2),['mw','m33',],[tab2['jh2'],ibb['jh2']],
#      ['black','green'])]


# [plt.annotate(wch+' '+str(round(np.nanstd(jhi),3)),(2e24,3-0.25*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(2),['mw','m33',],[tab2['jh2'],ibb['jh2']],
#      ['black','green'])]


# [plt.annotate(wch+' '+str(round(np.nanmedian(jhi),3)),(2e24,2-0.25*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(2),['mw','m33',],[tab2['jh2'],ibb['jh2']],
#      ['black','green'])]


plt.xlabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.ylabel('N (clouds)')
plt.xlim(1e23,5e26)
plt.ylim(0,8.5)
# plt.yticks(np.arange(5))

plt.savefig('jh2s',bbox_inches='tight')

plt.show()


ibb_omhi_lh2 = np.array([1.76,3.3,3.21,2.41,2.32,1.25,1.11,3.14,0.67,9.54,3.67,2.46,2.01,2.61,\
                           0.91,2.08,1.21,1.70,0.83,3.49,3.06,7.96,4.44,6.11,5.22,3.19,2.82,1.79,\
                           5.07,6.26,3.58,6.03,1.58,5.15,2.67,1.26,8.61,4.25,3.42,6.90,3.56,3.55])*0.01

ibb_mhi_lh2 = np.array([4.6,4.0,4.1,4.0,3.7,3.8,3.1,3.3,4.0,3.4,5.2,4.2,5.0,3.1,2.6,3.2,3.8,4.8,6.1,\
                         5.3,4.9,4.3,4.0,5.9,6.2,5.4,3.9,5.4,4.0,3.5,3.6,3.8,3.5,3.6,3.4,3.0,4.1,\
                         4.3,2.3,2.5,4.0,4.1])*0.01

ibb_jhi_lh2 = 0.4*((ibb_omhi_lh2*u.km/u.s/u.pc)*(70*u.pc)**2).cgs


plt.figure(dpi=200,figsize=(5,5))
bns=10**np.arange(22,27,1/4)

plt.hist(ibb_jhi_lh2.value,bins=bns,color='red',
         histtype='step',zorder=30,linewidth=1.5)
plt.hist(ibb['jhi'],bins=bns,color='green',histtype='step',zorder=30,linewidth=1.5,)
plt.hist(tab2_new['jhi'],bins=bns,color='black',zorder=25)

plt.vlines(np.nanmedian(ibb_jhi_lh2.value),0,18,linestyle='--',color='red',alpha=0.5)
plt.vlines(np.nanmedian(ibb['jhi']),0,18,linestyle='--',color='green',alpha=0.5)
plt.vlines(np.nanmedian(tab2_new['jhi']),0,18,alpha=0.5,color='black',linestyle='--')


# [plt.annotate(wch+' '+str(round(np.nanmean(jhi),3)),(2e24,8-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],[tab2['jhi'],ibb['jhi'],ibb_jhi_lh2.value],
#      ['black','red','green'])]


# [plt.annotate(wch+' '+str(round(np.nanstd(jhi),3)),(2e24,5-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],[tab2['jhi'],ibb['jhi'],ibb_jhi_lh2.value],
#      ['black','red','green'])]


# [plt.annotate(wch+' '+str(round(np.nanmedian(jhi),3)),(2e24,2-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],[tab2['jhi'],ibb['jhi'],ibb_jhi_lh2.value],
#      ['black','red','green'])]

plt.xscale('log')

plt.xlabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.ylabel('N (clouds)')
plt.xlim(1e23,5e26)
plt.ylim(0,17)
plt.yticks(np.arange(0,18,2))
plt.savefig('jhis',bbox_inches='tight')

plt.show()


plt.figure(dpi=200,figsize=(5,5))

plt.plot(np.linspace(1e23,1e27,5),
            6e6*np.linspace(1e23,1e27,5)**-0.25,color='black',linestyle='--',
            label=r'$\frac{j_\mathrm{HI}}{j_\mathrm{H_2}} \propto j_\mathrm{H_2}^{-0.25}$')


plt.plot(np.linspace(1e23,1e27,5),
            2e27*np.linspace(1e23,1e27,5)**-1.07,color='green',linestyle='--',alpha=0.2,
            label=r'$\frac{j_\mathrm{HI}}{j_\mathrm{H_2}} \propto j_\mathrm{H_2}^{-1.07}$')


plt.scatter(tab2_new['jh2'],tab2_new['jhi']/tab2_new['jh2'],color='black',edgecolor='white',s=80)


# plt.scatter(ibb_jh2,ibb_jhi/ibb_jh2,color='green',edgecolor='white',s=80,zorder=-10,alpha=0.2)
# plt.scatter(ibb_jh2,ibb_jhi/ibb_jh2,color='green',edgecolor='white',s=80,zorder=-10,alpha=0.2)


# plt.scatter(ib['jh2'],ib['jhi']/ib['jh2'],color='grey',edgecolor='white',s=50,zorder=20)
# [plt.plot([ib['jh2'].iloc[j],in_ib['jh2'].iloc[j]],
#           [(ib['jhi']/ib['jh2']).iloc[j],(in_ib['jhi']/in_ib['jh2']).iloc[j]],
#           linewidth=0.75,zorder=21,label=j) for j in range(5)]


plt.xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel(r'$j_\mathrm{HI}\ /\ j_\mathrm{H_2}$')

# plt.yticks([1,10],labels=['1','10'])
plt.xlim(3e23,2e26)
plt.ylim(5e-1,5e1)

plt.loglog()
plt.legend(fontsize='small')
# plt.savefig('jjj',bbox_inches='tight')


curve_fit(line,np.log10(tab2_new['rh2'].to_numpy()),np.log10(tab2_new['jh2'].to_numpy()))


np.sqrt(curve_fit(line,np.log10(tab2_new['rh2'].to_numpy()),np.log10(tab2_new['jh2'].to_numpy()))[1][0,0])


curve_fit(line,np.log10(tab2_new['rhi'].to_numpy()),np.log10(tab2_new['jhi'].to_numpy()))


np.sqrt(curve_fit(line,np.log10(tab2_new['rhi'].to_numpy()),np.log10(tab2_new['jhi'].to_numpy()))[1][0,0])


10**(22.66)


plt.figure(dpi=200,figsize=(5,5))

xs = np.linspace(1,200,10)


plt.scatter(tab2_new['rh2'],tab2_new['jh2'],color=(0.,0.8,0.),
            s=100,zorder=30,edgecolor='white')


# plt.scatter(ibb['rh2'],ibb['jh2'],color='black',alpha=0.5,
#             s=50,edgecolor='white',zorder=-20)
# plt.plot(xs, 5e22*xs**1.66,color='black',linestyle='--',linewidth=0.75,zorder=-30,alpha=0.5)


plt.xlim(4,200)
plt.ylim(3e23,4e26)
plt.loglog()

plt.plot(xs, 45e21*xs**1.54,color='black',linestyle='--')
plt.fill_between(xs,45e21*xs**(1.54+0.20),45e21*xs**(1.54-0.20),
                 alpha=0.2,color='black',linewidth=0.)


plt.xlabel(r'$R_\mathrm{H_2}$ [pc]')
plt.ylabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')

plt.annotate('$j_\mathrm{H_2} \propto R^{1.54 \pm 0.20}$',(4.5,1e26))
# plt.annotate('$j_\mathrm{H_2} \propto R^{1.66 \pm 0.19}$',(4.5,8e25),alpha=0.2)
plt.ylim(1e23,2e26)


plt.savefig('jR_h2',bbox_inches='tight')


plt.figure(dpi=200,figsize=(5,5))


plt.scatter(tab2_new['rhi'],tab2_new['jhi'],color=(0.8,0.,0.),
            s=100,edgecolor='white',zorder=30)

# [plt.vlines(tab2['rhi'][i],jhi1[i]-jhi[i],jhi3[i]-jhi[i]) for i in range(22)]


plt.xlim(20,300)
plt.ylim(15e23,6e26)
plt.loglog()

xs = np.linspace(10,300,10)
plt.plot(xs, 1e22*xs**1.76,color='black',linestyle='--')
plt.fill_between(xs,1e22*xs**(1.76+0.26),1e22*xs**(1.76-0.26),
                 alpha=0.2,color='black',linewidth=0.,zorder=10)


plt.annotate('$j \propto R^{1.76 \pm 0.26}$',(22,2e26))

# plt.vlines(70,1e23,7e26,zorder=-10,linestyle=':',color='black',alpha=0.5)

plt.xlabel(r'$R_\mathrm{HI}\ [\mathrm{pc}]$')
plt.ylabel(r'$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
plt.ylim(1e24,3e26)
plt.savefig('jR_hi',bbox_inches='tight')


curve_fit(line,np.log10(tab2_new['jh2']),np.log10(tab2_new['jhi']/tab2_new['jh2']))


np.sqrt(curve_fit(line,np.log10(tab2_new['jh2']),np.log10(tab2_new['jhi']/tab2_new['jh2']))[1][0,0])


np.log10(7e7)


plt.figure(dpi=200,figsize=(5,5))


plt.scatter(tab2_new['jh2'],tab2_new['jhi']/tab2_new['jh2'],color='black',
            s=100,edgecolor='white',zorder=30)

# [plt.vlines(tab2['rhi'][i],jhi1[i]-jhi[i],jhi3[i]-jhi[i]) for i in range(22)]


# plt.xlim(20,300)
# plt.ylim(15e23,6e26)
plt.loglog()

xs = np.linspace(1e23,1e27,10)
plt.plot(xs, 7e7*xs**-0.29,color='black',linestyle='--')
plt.plot(xs, 1e27*xs**-1.07,color='green',linestyle='--')

# plt.fill_between(xs,7e7*xs**(-0.29+0.12),7e7*xs**(-0.29-0.12),
#                  alpha=0.2,color='black',linewidth=0.,zorder=10)


plt.annotate(r'$\frac{j_\mathrm{HI}}{j_\mathrm{H_2}} \propto j_\mathrm{H_2}^{-0.29 \pm 0.12}$',(1e25,30))

# plt.vlines(70,1e23,7e26,zorder=-10,linestyle=':',color='black',alpha=0.5)

plt.xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2 s^{-1}}]$')
plt.ylabel(r'$j_\mathrm{HI}\ / \ j_\mathrm{H_2}$')
plt.ylim(5e-1,50)
plt.xlim(1e23,2e26)
plt.savefig('jjj',bbox_inches='tight')


plt.figure(dpi=200,figsize=(5,5))
plt.hist(np.arccos(np.cos(np.mod((tab2_new['ththi']-tab2_new['thth2'])*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg),
        bins=np.arange(0,361,15),edgecolor='white',color='black')

# plt.vlines(np.nanmean(np.arccos(np.cos(np.mod(ibb['deltht']*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg)),
#           0,5.5,color=(0.8,0.,0.),linestyle='--',linewidth=0.75,zorder=-10)

# plt.vlines(np.nanmean(np.arccos(np.cos(np.mod((tab2_new['ththi']-tab2_new['thth2'])*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg)),
#           0,5.5,color='black',linestyle='--')
# plt.hist(np.arccos(np.cos(np.mod(ibb['deltht']*u.deg.to(u.rad),2*np.pi)))*u.rad.to(u.deg),
#         bins=np.arange(0,361,15),color='green',histtype='step',linewidth=2)

# plt.hist(-360+np.mod((old_thhi-old_thh2)*u.deg.to(u.rad),2*np.pi)*u.rad.to(u.deg),
#         bins=np.arange(-360,361,10),edgecolor='white',color='red')

plt.xticks(np.arange(-360,361,30),labels = [str(num)+r'$\degree$' for num in np.arange(-360,361,30)],
           )

plt.xlim(-10,190)
plt.ylim(0,5.5)
plt.ylabel('N (clouds)')
# plt.grid()
plt.xlabel(r'$| \theta_\mathrm{HI} - \theta_\mathrm{H_2} |$')


plt.savefig('theta_h2_hi',bbox_inches='tight')

plt.show()
# plt.xticks


(u.pc*u.km/u.s).cgs


galhi = np.arccos(np.cos(tab2_new['ththi']*u.deg)*np.cos(bc*u.deg))*u.rad.to(u.deg)
galh2 = np.arccos(np.cos(tab2_new['thth2']*u.deg)*np.cos(bc*u.deg))*u.rad.to(u.deg)


plt.figure(dpi=200,figsize=(6,6))
plt.hist(np.mod(np.array(galh2)*u.deg.to(u.rad),2*np.pi)*u.rad.to(u.deg),
        bins=np.arange(0,361,15),edgecolor='white',color='black')


plt.xticks(np.arange(-360,361,30),labels = [str(num)+r'$\degree$' for num in np.arange(-360,361,30)])

plt.xlim(-10,190)
plt.ylim(0,5.5)
plt.ylabel('N (clouds)')
# plt.grid()
plt.xlabel(r'$| \theta_\mathrm{MW} - \theta_\mathrm{H_2} |$')

plt.savefig('theta_mw_h2',bbox_inches='tight')

plt.show()
# plt.xticks


plt.figure(dpi=200,figsize=(6,6))
plt.hist(galhi,bins=np.arange(0,361,15),edgecolor='white',color='black')

plt.xticks(np.arange(-360,361,30))#,labels = [str(num)+r'$\degree$' for num in np.arange(-360,361,30)])

plt.xlim(-10,190)
plt.ylim(0,5.5)
plt.ylabel('N (clouds)')
# plt.grid() 
plt.xlabel(r'$| \theta_\mathrm{MW} - \theta_\mathrm{HI} |$')

plt.savefig('theta_mw_hi',bbox_inches='tight')
plt.show()
# plt.xticks


p99_tab1 = pd.read_csv('../J_A+AS_134_241_table1.dat.gz.txt',skiprows=9,delimiter='|')
p99_tab1.columns = ['Name','Note','Ah Am RAs Ed Em DEs','l_Dist','Dist','u_Dist','RV','DV DV2',\
                   'l_Size','Size','n_Size','Size2','Line','L_OMEGA','OMEGA_OMEGA2','l_GAMMA',\
                   'GAMMA','l_Mass','Mass_Mass2','MassEm','THETA_THETA2','OMEGAo','u_OMEGAo','Type','r_']


pOmgs = np.array([ float(n.split()[0]) 
                  for i, n in enumerate(p99_tab1['OMEGA_OMEGA2'][1:-1])]) *\
                    1e-14*(1/u.s).to(u.km/u.s/u.pc).value

pRs = p99_tab1['Size'][1:-1].astype('float').to_numpy()

pMsraw = [pM.split() for pM in p99_tab1['Mass_Mass2'][1:-1].to_numpy()]

pTyps = p99_tab1['Type'][1:-1].to_numpy()


pMs = np.zeros(len(pMsraw))
for i,n in enumerate(pMsraw):
    if len(n)==0: # ""
        pMs[i] = 0
        
    elif len(n)==1: # "num"
        pMs[i] = float(n[0])
        
    elif len(n)==3: #"num","->","num"
        pMs[i] = np.nanmean(np.array(n)[[0,2]].astype('float'))
        
    elif n[0][-1]=='>':
        pMs[i] = np.nanmean(np.array((n[0][:6],n[1])).astype('float'))
        
    elif n[1][:2]=='->':
        pMs[i] = np.nanmean(np.array((n[0],n[1][2:])).astype('float'))
        
    else:
        pMs[i] = np.nanmean(np.array(pMsraw[126]).astype('float'))
        
        
p99_dat = pd.DataFrame(np.array([pTyps,pMs,pRs,pOmgs]).transpose(),
                      columns=['type','M<Msun>','R<pc>','Omega<km/s/pc>'])

p99_dat = p99_dat[~np.isinf(np.log10(pMs))] #M > 0
p99_dat = p99_dat[~(p99_dat['type']=='D/R ')] #no disks/rings
p99_dat = p99_dat[~np.array([('F' in nm) for nm in p99_dat['type']])] #no filaments

p99_dat['j<cm2/s>'] = ((p99_dat['R<pc>'].to_numpy()*u.pc)**2 *\
                        p99_dat['Omega<km/s/pc>'].to_numpy()*u.km/u.s/u.pc).cgs.value.astype('float')


chenx_rs = np.array([1703,1232,3530,1440,1272,1763,1329,2762])*u.AU.to(u.pc)
chenx_js = (np.array([0.12,0.12,0.45,0.13,0.24,0.12,0.42,0.45])*\
            1e-3*u.km/u.s*u.pc).cgs.value


chenh_rs = np.array([0.038,0.033,0.025,0.032,0.026,0.026,
                     0.034,0.054,0.042,0.020,0.022,0.024,0.19,0.028,0.045,
                    0.044,0.050,0.048])

chenh_js = (np.array([1.3,1.3,0.36,1.3,0.8,0.52,0.59,2.5,1.,0.8,1.5,\
 1.6,6.,0.37,0.25,0.34,0.46,0.54])*1e-3 * u.km/u.s*u.pc).cgs


goodman = pd.read_csv('../csvs/goodman_1993_data.csv')
tatmtsu_rs = np.array([
    0.049,0.088,0.083,0.084,0.059,0.099,0.091,0.049,0.098,0.095,0.052,\
    0.104,0.094,0.083,0.111,0.105,0.074,0.057,0.072,0.043,0.062,0.074,0.115,\
    0.083,0.084,0.062,0.074,0.081,0.062,0.052,0.083,0.157,0.057,0.072
])
tatmtsu_js = np.array([
    7.4e20,1.7e21,7.6e20,7.2e20,0,6e20,1.2e21,1.6e21,9.3e21,3.2e21,3e20,\
    2.3e21,1.8e21,2.2e21,2.1e21,1.9e21,2.2e21,0,0,8.4e20,1.5e21,0,1.9e21,4.4e21,\
    0,7.7e20,0,1.1e21,1.2e21,1.1e21,3e21,3.5e21,8.9e20,0
])


tatmtsu_rs = tatmtsu_rs[~np.isinf(np.log10(tatmtsu_js))]
tatmtsu_js = tatmtsu_js[~np.isinf(np.log10(tatmtsu_js))]


all_rs = np.concatenate((tab2_new['rh2'].to_numpy(),tab2_new['rhi'].to_numpy(),\
         p99_dat['R<pc>'].to_numpy().astype('float'),chenh_rs,chenx_rs,\
                         tatmtsu_rs,ibb['rh2'],))

all_js = np.concatenate((tab2_new['jh2'].to_numpy(),tab2_new['jhi'].to_numpy(),\
         p99_dat['j<cm2/s>'].to_numpy().astype('float'),chenh_js.value,chenx_js,\
                         tatmtsu_js,ibb['jh2'].to_numpy(),))


curve_fit(line,np.log10(all_rs),np.log10(all_js))


np.sqrt(curve_fit(line,np.log10(all_rs),np.log10(all_js))[1][0,0])


plt.figure(dpi=200,figsize=(5,5))


plt.scatter(tab2_new['rh2'],tab2_new['jh2'],
            color=(0.,0.8,0.),s=40,edgecolor='white',label='$\mathrm{H_2}$ (this work)')
plt.scatter(tab2_new['rhi'],tab2_new['jhi'],
            color=(0.8,0.,0.),s=40,edgecolor='white',label='HI (this work)')


plt.scatter(ibb['rh2'],ibb['jh2'],color='black',alpha=0.5,
            s=40,edgecolor='green',zorder=-20,label='Imara +11')


plt.scatter(p99_dat['R<pc>'],p99_dat['j<cm2/s>'],label='Phillips 99',
            color='black',marker='*',s=90,edgecolor='white')

plt.scatter(chenx_rs,chenx_js,s=40, 
            color='black',marker='<',edgecolor='white',label='Chen +07')

plt.scatter(tatmtsu_rs,tatmtsu_js,s=40,
            color='black',marker='o',edgecolor='white',label='Tatematsu +16')


plt.scatter(chenh_rs,chenh_js,s=40,
            color='black',marker='s',edgecolor='white',label='Chen +19')


# plt.scatter(np.full(len(ibb),70),ibb['jhi'],color='black',alpha=0.5,
#             s=70,edgecolor='red',zorder=30,label='ibb hi')

# plt.scatter(np.full(len(ibb_jhi_lh2),70),ibb_jhi_lh2,color='black',alpha=0.5,
#             s=70,edgecolor='orange',zorder=20,label='ibb hi wo h2')
 
# plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],
#             color=(0.6,0.6,0.6),marker='^',s=30,)


plt.xlim(3e-3,4e2)
plt.ylim(2e19,8e26)
plt.loglog()


xs = np.linspace(1e-3,5e2,10)
plt.plot(xs, 1e23*xs**1.46,color='black',linestyle='--',zorder=-10)
plt.fill_between(xs,1e23*xs**(1.46+0.03),1e23*xs**(1.46-0.03),
                 alpha=0.2,color='black',linewidth=0.,zorder=-10)


plt.annotate('$j \propto R^{1.46 \pm 0.03}$',(5e-3,2e26))

plt.legend(fontsize='xx-small')


plt.xlabel('R [pc]')
plt.ylabel('j [$\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('jR_all',bbox_inches='tight')


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den

def rvsv_pearson(i, show=False): 
    
    ccoords = coordses[i]
    d = ds[i]
    omgs = tab2_new[['omgh2','omghi']].to_numpy()[i]/1e3
    thts = tab2_new[['thth2','ththi']].to_numpy()[i]
    
    l0c,l1c,b0c,b1c = ccoords
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2maps, himaps = moment_maps(i)
    h2_0,h2_1,h2_2 = h2maps
    hi_0,hi_1,hi_2 = himaps
    
    rhi_px = hi_0.shape[0]/2

    l0h,l1h = lc + rhi_px*0.125, lc - rhi_px*0.125
    b0h,b1h = bc - rhi_px*0.125, bc + rhi_px*0.125

    hlons,hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))
    clons,clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))
    
    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*d
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*d

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*d
    cys_pc = (clats-bc)*u.deg.to(u.rad)*d
    
    
    h2_dvdss, hi_dvdss = [(om*np.cos(th*u.deg).value, om*np.sin(th*u.deg).value) for om,th in zip(omgs,thts)]
    
    hbperp = (hi_dvdss[0]*hxs_pc + hi_dvdss[1]*hys_pc)/omgs[1]
    cbperp = (h2_dvdss[0]*cxs_pc + h2_dvdss[1]*cys_pc)/omgs[0]


    hv_model =  hi_vos[i] + hi_dvdss[0]*hxs_pc + hi_dvdss[1]*hys_pc
    cv_model =  h2_vos[i] + h2_dvdss[0]*cxs_pc + hi_dvdss[1]*cys_pc
    
#     print(hbperp.shape)
#     plt.contourf(cbperp,cmap='PiYG')
#     plt.contour(cv_model)
#     plt.contourf(cv_model)
#     plt.contourf(h2_1,cmap='RdGy')
#     plt.gca().set_aspect(1)
#     plt.show()
    
    
#     plt.contourf(hv_model,cmap="RdGy")
#     plt.gca().set_aspect(1)
#     plt.colorbar()
#     plt.show()
    
#     plt.contourf(hi_1,cmap="PiYG")
#     plt.gca().set_aspect(1)
#     plt.colorbar()
#     plt.show()

    
    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hi_2),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod',])

#     hdf_noisy = hdf
#     hdf = hdf_noisy[hdf_noisy['w']>hi_Nlims[i][0]]
    
    
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
            np.ravel(h2_2),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod'])
    

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()
    
#     print(cdf)


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    h2p, hip = pearson(meancdf['rper'], meancdf['v']),pearson(meanhdf['rper'], meanhdf['v'])

#     print(meancdf['v'],meancdf['vmod'])
    
    if show:
    
        fig,axs = plt.subplots(1,2,figsize=(12,5))
#         axs[0].scatter(cdf['rper'],cdf['v'])
#         axs[1].scatter(hdf['rper'],hdf['v'])
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

        axs[0].annotate(str(np.around(h2p,2)),(meancdf['rper'].iloc[0],
                                               np.nanmax(meancdf['v'])+np.nanmean(sigcdf['v'])))

        axs[1].annotate(str(np.around(hip,2)),(meanhdf['rper'].iloc[0],
                                               np.nanmax(meanhdf['v'])+np.nanmean(sighdf['v'])))


        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')

        axs[0].set_ylabel(r'velocity [$\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names[i])

        fig.savefig('v_vs_r_'+names[i][:3]+names[i][-1],bbox_inches='tight')


        plt.show()
    
    return h2p,hip


rvsv_pearson(15,show=True)


h2p,hip = np.array([rvsv_pearson(i) for i in range(22)]).transpose()


[get_stats(pval) for pval in [h2p,hip]]


find_which(get_stats(jh2.value/jinits)[:2],jh2.value/jinits)


plt.scatter(rh2,jh2.value/jinits)
plt.loglog()


hi_vos = [fit_plane2maps(i,1)[2]/1e3 for i in range(22)]


h2_vos = [fit_plane2maps(i,0)[2]/1e3 for i in range(22)]


[rvsv_pearson(i,show=True) for i in range(22)]


ph2,phi=np.array([rvsv_pearson(i) for i in range(22)]).transpose()


plt.figure(dpi=200,figsize=(5,5))

bns=np.arange(-0.5,1.1,1/10)

plt.hist(np.abs(ph2),color='black',edgecolor='white',bins=bns)
plt.hist(np.abs(phi),histtype='step',linewidth=2,color='red',bins=bns)

plt.ylabel('N (clouds)')
plt.xlabel('$p (\Delta v,r_\perp)$')
# plt.xlim(-0.1,1.1)
plt.savefig('pvals',bbox_inches='tight')
plt.ylim(0,16)


def t_ff(M,R):
    rho = M*u.solMass/(4/3*np.pi*(R*u.pc)**3)
    num = 3*np.pi
    den = 32*c.G*rho
    
    return np.sqrt(num/den).to(u.s)#.to(u.Myr)


tffh2 = t_ff(mh2,rh2)
tffhi = t_ff(mhi,rhi)


np.nanmean(trothi/tffhi.to(u.Myr))


np.nanmedian(trothi/tffhi.to(u.Myr))


np.nanmedian(troth2/tffh2.to(u.Myr))


plt.plot(trothi/tffhi.to(u.Myr))


tffh2.to(u.Myr)


#do the calculation for t diss
#calculate free fall times 


sigvhi


def get_tdis(dat):
    egrav_hi = (c.G * dat['mhi'].to_numpy()*u.solMass / (dat['rhi'].to_numpy()*u.pc)).cgs
    erot_hi = (0.4*dat['jhi'].to_numpy()*u.cm**2/u.s*dat['omghi'].to_numpy()*u.km/u.s/u.pc).cgs
    ekin_hi = (0.5*(sigvhi*u.km/u.s)**2).cgs

    etot_hi = erot_hi - egrav_hi + ekin_hi
    
    tdis = (dat['jhi']/etot_hi.value)#*u.s.to(u.Myr)
    
    jinit = (dat['mh2']*dat['jh2'] + dat['mhi']*dat['jhi'])/(dat['mh2'] + dat['mhi'])

    
    return tdis,jinit


tdises,jinits = get_tdis(tab2_new)


troth2


bns


np.nanmean(tdises*u.s.to(u.Myr))


plt.figure(dpi=200)
bns = 10**np.arange(11,16,1/7)#*u.s.to(u.Myr)

plt.hist(tdises,bins=bns,color='black')

plt.hist(tffhi.value,alpha=0.5,bins=bns,color='red')
plt.hist(tffh2.value,alpha=0.5,bins=bns,color='green')

plt.hist(troth2.value,bins=bns,color='red',histtype='step')
plt.hist(trothi.value,bins=bns,color='green',histtype='step')

plt.xscale('log')
# plt.xlim(2e11,1e16)


plt.hist(jinits,bins=bns)
plt.hist(jh2.value,alpha=0.5,bins=bns)
plt.hist(jhi.value,alpha=0.5,bins=bns)

plt.xscale('log')


bns = 10**np.arange(22,27,1/4)


plt.figure(dpi=200,figsize=(5,5))
plt.hist(jinits,color='black',zorder=20,
          bins=bns) 

plt.hist(jh2.value,color='green',linewidth=2,zorder=25,
          bins=bns,histtype='step') 


plt.hist(jhi.value,color='red',linewidth=2,zorder=25,
          bins=bns,histtype='step') 


plt.vlines(np.nanmedian(jh2.value),0,12,linestyle='--',color='green',alpha=0.5)
plt.vlines(np.nanmedian(jhi.value),0,12,linestyle='--',color='red',alpha=0.5)
plt.vlines(np.nanmedian(jinits),0,12,linestyle='--',color='black')


# [plt.annotate(wch+' '+str(round(np.nanmean(jhi),3)),(2e24,8-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],
#      [get_tdis(tab2)[1],get_tdis(ibb)[1],ibb_jhi_lh2.value],
#      ['black','green','red'])]


# [plt.annotate(wch+' '+str(round(np.nanstd(jhi),3)),(2e24,5-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],
#      [get_tdis(tab2)[1],get_tdis(ibb)[1],ibb_jhi_lh2.value],
#      ['black','green','red'])]


# [plt.annotate(wch+' '+str(round(np.nanmedian(jhi),3)),(2e24,2-0.5*i),fontsize='small',color=clr) 
#  for i,wch,jhi,clr in 
#  zip(range(3),['mw','m33','m33woh2'],
#      [get_tdis(tab2)[1],get_tdis(ibb)[1],ibb_jhi_lh2.value],
#      ['black','green','red'])]


plt.xlabel('$j_\mathrm{init}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel('N(clouds)')
plt.xscale('log')
plt.ylim(0,6.5)
plt.xlim(1e23,5e26)
# plt.loglog()
plt.savefig('jinit',bbox_inches='tight')


plt.scatter(tdises,tffhi)
plt.scatter(tdises,tffh2)

plt.loglog()


plt.scatter(tdises/tffhi, jh2.value/jinits)


plt.scatter(tdises/tffh2, jh2.value/jinits)


tffh2.value/tdises


(1/(omgh2*u.m/u.s/u.pc)).to(u.Myr)


(1/(omghi*u.m/u.s/u.pc)).to(u.Myr)


ththi[8]


omgh2_old = np.array([87.6,51.0,65.5,40.7,94.1,111.2,65.8,19.3,323.1,176.8,67.0,\
                      50.6,107.77,50.6,79.2,70.6,44.5,160.6,159.7,63.8,87.6,120.3])
omghi_old = np.array([80.8,58.1,21.9,12.0,15.2,2.2,31.8,7.0,68.5,15.6,15.8,19.9,\
                      7.8,17.7,29.6,4.8,9.5,22.2,18.8,22.8,10.5,62.3])

rh2_old = np.array([5.45,38.66,5.70,32.43,8.92,4.74,19.72,98.36,7.25,13.34,16.32,\
                   57.88,24.66,29.72,15.04,24.40,17.91,25.99,15.50,17.15,30.71,10.16])
rhi_old = np.array([22.25, 123.67, 21.75,110.41,46.57,28.5,69.67,245.53,38.89,43.87,\
                   64.65,153.29,104.21,112.08,62.39,108.40,94.22,71.79,76.86,63.53,109.73,57.72])

jh2_old = ((omgh2_old*u.m/u.s/u.pc)*(rh2_old*u.pc)**2).cgs
jhi_old = ((omghi_old*u.m/u.s/u.pc)*(rhi_old*u.pc)**2).cgs

mhi_old = np.array([8.67,984.9,7.43,543.64,40.47,10.80,208.68,5105.2,67.14,48.53,\
                    139.80,1774.90,524.75,520.58,144.76,281.99,379.94,\
                    322.44,114.95,56.39,699.24,131.71])*1e3

mh2_old = np.array([0.88,151.78,0.83,108.02,8.11,1.86,27.14,11187.8,\
                    4.72,6.78,21.68,289.06,90.82,112.99,19.49,102.22,\
                   67.13,29.69,36.43,20.57,106.03,15.43])*1e3

oldtab2 = pd.DataFrame(np.array((mh2_old,rh2_old,omgh2_old,jh2_old,\
                                 mhi_old,rhi_old,omghi_old,jhi_old)).transpose(),
                        columns = 'mh2,rh2,omgh2,jh2,mhi,rhi,omghi,jhi'.split(','))


from astropy.table import QTable


oldoldcat = QTable(fits.open('../fitses/re_expanded_cat.fits')[1].data)


plt.scatter(oldoldcat['SFR<Msun/yr>'],tdises)


oldjh2 = (oldoldcat['spec_ang_mom_GMC<pc.km/s>']*u.pc*u.km/u.s).cgs

oldjhi = (oldoldcat['spec_ang_mom_HI<pc.km/s>']*u.pc*u.km/u.s).cgs


plt.plot(tab2_new['jhi'],tab2_new['jhi'])
plt.scatter(tab2_new['jhi'],jhi_old)
plt.scatter(tab2_new['jhi'],oldjhi)

plt.loglog()


plt.plot(tab2_new['jh2'],tab2_new['jh2'])
plt.scatter(tab2_new['jh2'],jh2_old)
plt.scatter(tab2_new['jh2'],oldjh2)

plt.loglog()


(2*np.pi/(12.2*1e-2*(u.km/u.s/u.pc))).to(u.s).to(u.Myr)


troth2 = (2*np.pi/(np.array(omgh2)*1e-3*(u.km/u.s/u.pc))).to(u.s)#.to(u.Myr)
trothi = (2*np.pi/(np.array(omghi)*1e-3*(u.km/u.s/u.pc))).to(u.s)#.to(u.Myr)


np.nanmean(troth2.to(u.Myr))


np.nanmean(trothi)


get_ind_std(rh2,jh2.value)


get_ind_std(rhi,jhi.value)


plt.figure(dpi=200)
plt.hist(np.abs(ththi),bins=np.arange(-180,180,15),edgecolor='white')
plt.show()


def get_stats(arr):
    return np.around([np.nanmin(arr),np.nanmax(arr), np.nanmean(arr), np.nanmedian(arr)],2)


[get_stats(tab2_new[col]/fact) 
for col,fact in zip(tab2_new.columns,[1e4,1,1e1,1,1e24,1e5,1,1e1,1,1e24])]


[(col,get_stats(np.mod(tab2_new[col]*u.deg.to(u.rad),np.pi)*u.rad.to(u.deg)))
 
for col,fact in zip(tab2_new.columns,[1e4,1,1e1,1,1e24,1e5,1,1e1,1,1e24])]


np.nanmedian(hip)


[(hi+'/'+h2,
 find_which(get_stats(np.mod(np.abs(tab2_new[hi]-tab2_new[h2])\
                             *u.deg.to(u.rad),np.pi))[:2]*u.rad.to(u.deg),
           np.mod(np.abs(tab2_new[hi]-tab2_new[h2])\
                             *u.deg.to(u.rad),np.pi)*u.rad.to(u.deg))
 ) for hi,h2 in zip(tab2_new.columns[5:],tab2_new.columns[:5])]


mhi/1e4


mhi/1e4


tab2_new.columns[0],tab2_new.columns[5]


[(hi+'/'+h2,
  get_stats(tab2_new[hi]/tab2_new[h2])) for hi,h2 in zip(tab2_new.columns[5:],tab2_new.columns[:5])]


1/0.4


[get_stats(np.mod(np.abs(tab2_new[hi]-tab2_new[h2])*u.deg.to(u.rad),np.pi))*u.rad.to(u.deg) 
           for hi,h2 in zip(tab2_new.columns[5:],tab2_new.columns[:5])]


[(hi,h2)for hi,h2 in zip(tab2_new.columns[5:],tab2_new.columns[:5])]


get_stats(jinits/1e24)


def find_which(vals,arr):
    
    return [np.argmin(np.abs(arr - val)) for val in vals]
    
    
#     return np.around([np.nanmin(arr),np.nanmax(arr), np.nanmean(arr), np.nanmedian(arr)],2)


find_which([0.74e24,197.3e24],jhi.value)


[(col,find_which(get_stats(tab2_new[col]/fact)[:2],tab2_new[col]/fact))
for col,fact in zip(tab2_new.columns,[1e4,1,1e1,1,1e24,1e5,1,1e1,1,1e24])]


[(hi+'/'+h2,
 find_which(get_stats(tab2_new[hi]/tab2_new[h2])[:2],tab2_new[hi]/tab2_new[h2]))
 for hi,h2 in zip(tab2_new.columns[5:],tab2_new.columns[:5])]


find_which(get_stats(jinits/1e24)[:2],jinits/1e24)


[(hi+'/'+h2,
 find_which(get_stats(np.mod(np.abs(tab2_new[hi]-tab2_new[h2])\
                             *u.deg.to(u.rad),np.pi))[:2]*u.rad.to(u.deg),
           np.mod(np.abs(tab2_new[hi]-tab2_new[h2])\
                             *u.deg.to(u.rad),np.pi)*u.rad.to(u.deg))
 ) for hi,h2 in zip(tab2_new.columns[5:],tab2_new.columns[:5])]


def get_ind_std(x,y,mod=line):
    fit = curve_fit(mod,np.log10(x),np.log10(y))
    return np.around((fit[0][0], np.sqrt(fit[1][0,0])),2)


get_ind_std(mh2,jh2.value)


get_ind_std(rh2,jh2.value)


get_ind_std(mhi,jhi.value)


get_ind_std(rhi,jhi.value)


get_ind_std(rh2,omgh2)


get_ind_std(rhi,omghi)


get_ind_std(jhi/jh2,jh2.value)


ibb     = pd.read_csv('ibb_tab1.csv',sep=',',nrows=33)
ibb_jhi = ((ibb['omghi'].to_numpy()*0.01*u.km/u.s/u.pc)*(70*u.pc)**2).cgs
ibb_jh2 = ((ibb['omgh2'].to_numpy()*0.01*u.km/u.s/u.pc)*(ibb['rh2'].to_numpy()*u.pc)**2).cgs
ibb['jh2'] = np.around(ibb_jh2.value,2)
ibb['jhi'] = np.around(ibb_jhi.value,2)
ibb['rhi'] = np.full(len(ibb),70)
ibb['thth2'] = ibb['ththi'] - ibb['deltht']


ibb_new = ibb[tab2_new.columns]


# ibjh2 = (np.array((3.92,25.5,42.2,24.2,25.8))*u.pc*u.km/u.s).cgs #these are with c=0.5
# ibjhi = (np.array((25.2,120,107,60.3,73.4))*u.pc*u.km/u.s).cgs
ibrh2 = np.array([14.6,30.6,21.9,16.2,26.5,])
ibrhi = np.array([51.9,83.6,76.0,42.4,78.8])
iboh2 = np.array([4.6,6.68,22.,23.1,9.18])
ibohi = np.array([1.87,3.48,3.71,6.7,2.36])
ibjh22 = (0.4*(iboh2*1e-2*u.km/u.s/u.pc)*(ibrh2*u.pc)**2).cgs
ibjhi2 = (0.4*(ibohi*1e-2*u.km/u.s/u.pc)*(ibrhi*u.pc)**2).cgs

ibmh2 = np.array([2.5,3.0,2.4,1.4,4.0])*1e4*u.solMass

ib = pd.DataFrame(np.array((ibmh2,ibrh2,iboh2,ibjh2,ibrhi,ibohi,ibjhi)).transpose(),
                  columns='mh2,rh2,omgh2,jh2,rhi,omghi,jhi'.split(','))


cosb = np.cos(bc*u.deg.to(u.rad))


plt.figure(dpi=200,figsize=(5,5))

plt.hist(np.abs(cosb*thth2),bins=np.arange(0,200,15),edgecolor='white',color='black')
plt.xlim(-5,185)
plt.xticks(np.arange(0,181,30))
plt.yticks([1,2,3])
plt.xlabel(r'$\phi_\mathrm{H_2}\ [\degree]$')

plt.vlines([60,120],0,4.5,color='red',linestyle='--')
plt.ylim(0,3.5)
plt.savefig('phih2',bbox_inches='tight')
plt.show()


plt.figure(dpi=200,figsize=(5,5))

plt.hist(np.abs(cosb*ththi),bins=np.arange(0,200,15),edgecolor='white',color='black')
plt.xlim(-5,185)
plt.xticks(np.arange(0,181,30))
plt.yticks([1,2,3,4])
plt.xlabel(r'$\phi_\mathrm{HI}\ [\degree]$')

plt.vlines([60,120],0,4.5,color='red',linestyle='--')

plt.ylim(0,4.5)

plt.savefig('phihi',bbox_inches='tight')

plt.show()


delphi = np.abs(np.array(ththi)-np.array(thth2))*cosb


plt.hist(np.mod(delphi*u.deg.to(u.rad,np.pi))


delthet = np.abs(np.array(thth2)-np.array(ththi))
for i in range(22):
    if delthet[i] > 180:
        delthet[i] = np.abs(-360+delthet[i])
        
#     if delthet[i] < -180:
#         delthet[i] += 360
        

delthet


plt.hist(np.abs(ibb['deltht']))


plt.figure(dpi=200,figsize=(5,5))

plt.hist(np.abs(cosb*delthet),bins=np.arange(0,200,15),edgecolor='white',color='black')

plt.hist(np.abs(ibb['deltht']),bins=np.arange(0,200,15),edgecolor='green',histtype='step')

plt.xlim(-5,185)
plt.xticks(np.arange(0,181,30))
plt.yticks([1,2,3,4,5])
plt.xlabel(r'$\Delta \phi\ [\degree]$')

plt.vlines([60,120],0,5.5,color='red',linestyle='--')

plt.ylim(0,5.5)

plt.annotate('corotating',(15,5.2),fontsize='x-small')
plt.annotate('orthogonal',(72,5.2),fontsize='x-small')

plt.annotate('counterrotating',(128,5.2),fontsize='x-small')


plt.savefig('delphi',bbox_inches='tight')

plt.show()


delphi,np.mod(delphi,180)


help(np.mod)


plt.scatter(np.mod(np.abs(cosb*ththi),180),np.abs(cosb*ththi))


def make_stacked_spec(i,n,thet=90)

for n in range(10):
    xo, xi =  n*33,(n+1)*33
    yo, yi =  n*33,(n+1)*33


i = 21

plt.figure(dpi=200,figsize=(7,10))
for n in range(10):
    
    xo, xi =  n*33,(n+1)*33
    yo, yi =  n*33,(n+1)*33
    
    
    this_vo = first_moment(np.nanmean(
        hihi[i].data[cat['vo'][i]:cat['vi'][i],
                    yo:yi,xo:xi],axis=(1,2)),
        np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)
    
    
    # vo_clr = (this_vo - h1_lvls[i][0] - np.nanmean(hi_1s[i].data))/(h1_lvls[i][1]- h1_lvls[i][0] - np.nanmean(hi_1s[i].data))
    
    if this_vo - np.nanmean(hi_1s[i].data) < 0:
        vo_clr = 0.5*(this_vo - h1_lvls[i][0])/(np.nanmean(hi_1s[i].data)-h1_lvls[i][0])
    
    else:
        vo_clr = 0.5+ 0.5*(this_vo - np.nanmean(hi_1s[i].data))/\
             (h1_lvls[i][1] - np.nanmean(hi_1s[i].data))
            
#     print(round(vo_clr,3))
    
    plt.plot(np.arange(155)*1.3-100,n*80 + np.nanmean(hihi[i].data[:,yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,alpha=0.5,linestyle=':'),
    
    
    plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
             n*80 + np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
             color=seg_spec_all_color(vo_clr),
            linewidth=1.5,),
    
    plt.scatter(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],yo:yi,xo:xi],axis=(1,2)),
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100),n*80,
                color=seg_spec_all_color(vo_clr)),
    
    plt.hlines(n*80,-100,105,linewidth=0.2,color='black',linestyle='--'),
    plt.plot(np.arange(155)*1.3-100,n*80 + 
             hico_scale[i]*np.nanmean(coco[i].data,axis=(1,2)),
             color='black',linewidth=0.6,alpha=0.5)
    
    
    plt.annotate(str(round((np.nanmean([xo,xi])+cat['loh'][i])*-0.125+360))+'$\degree$, '+
                 str(round((np.nanmean([yo,yi])+cat['boh'][i])*0.125-50))+'$\degree$',(-70,n*80+10),fontsize='small')
    
    
    # print(fns.get_mom1(np.nanmean(hihi[i].data[cat['vo'][i]:cat['vi'][i],n*13:(n+1)*13,60:75],axis=(1,2)),np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)/(h1_lvls[i][1]-h1_lvls[i][0]))
    

plt.xlim(-75,45)
plt.xlabel('$\mathrm{v_{LSR}\ [km\ s^{-1}]}$')

# plt.gca().set_aspect(0.2)

# plt.xticks(np.arange(-75,60,15))

plt.vlines(np.nanmean(hi_1s[i].data),0,810,color='black',linestyle='--',zorder=-20,linewidth=0.7,alpha=0.5)

plt.yticks([])
plt.ylim(-10,806)


# plt.title('Hercules',fontsize='x-small',y = 0.96,x=0.1)
# plt.savefig(names[i][:3]+'_stacked_patch_spectra',bbox_inches='tight')
