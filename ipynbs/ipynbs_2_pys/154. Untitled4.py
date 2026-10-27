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
gaus = models.Gaussian1D(2.,5.,5.)
# gaus = models.Gaussian1D(5.,10.,2.)
fttr = fitting.LMLSQFitter()


# x,y,z = np.indices((101,101,101))-50


# sphere = (np.sqrt(x**2 + y**2 + z**2) <=50)
# cylinder = (np.sqrt(z**2 + x**2 ) <=50)


# sphere_2d = np.nansum(sphere,axis=0)
# cylndr_2d = np.nansum(cylinder,axis=0)


# fig,axs = plt.subplots(1,2,dpi=200)
# axs[0].contourf(cylndr_2d/50,cmap='magma')
# one=axs[1].contourf(sphere_2d/50,cmap='magma')

# for ax in axs:
#     ax.set_aspect(1)
#     ax.set_xticks(np.linspace(0,100,5),labels=['-2R','-R','0','+R','+2R'])
#     ax.set_yticks(np.linspace(0,100,5),labels=['-2R','-R','0','+R','+2R'])
#     ax.label_outer()
    
# # fig.colorbar(one,axs[2])
# # plt.colorbar()


def get_subcube(bigcube, coords, dx=-0.125, dy=0.125, x0=360, y0=-50):
    lo,li,bo,bi = coords # assume in degrees
    xo,xi = ((np.array([lo,li]) - x0)/dx).astype('int')
    yo,yi = ((np.array([bo,bi]) - y0)/dy).astype('int')
    
    return bigcube[:,yo:yi,xo:xi]


def make_h2_hi_cubes(i, rhi,vovi):
    coords = coordses[i]
    vo,vi  = vovi

    loh,lih = lc_p[i]+np.array([-1,1])*rhi
    boh,bih = bc_p[i]+np.array([-1,1])*rhi

    hicube = hi4pi.data[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords[:4]) #velocity axis length = 155
    
    return h2cube,hicube


def moment_maps(i,thresh, rhi, vovi):
    h2cube,hicube = make_h2_hi_cubes(i,rhi, vovi)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.29*1.82e18
    
    vo,vi = vovi

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.29-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.29-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.log10(hilb)>20) #*\
             # (np.linalg.norm(np.indices(hi1.shape) - \
             #                np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi)
                            
#     print(hi_msk)
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


#------------------------------------------------------------------------------------------------------------------


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


#------------------------------------------------------------------------------------------------------------------

def calc_masses(i,rhi,vovi): #cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    h20,hi0 = [moment_maps(i,threshs[i],rhi,vovi)[j][0]/u.cm**2 for j in [0,1]]
    
    
    dA = (ds[i]*u.pc*0.125*u.deg.to(u.rad))**2
    h2m = (np.nansum(h20*dA)*c.m_p*1.36*2).to(u.solMass)
    him = (np.nansum(hi0*dA)*c.m_p*1.36).to(u.solMass)
    
    return h2m.value,him.value
    

def get_mh2_rh2_rhi(i,thresh):
    coords = coordses[i]
#     print(coords)
    h2cube = get_subcube(dht.data,coords[:4])
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=thresh)

    rhi = (((3*mh2*u.solMass)/(4*np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
    return mh2,rh2,rhi


def calc_mass(cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    zeroth = np.nansum(cube[vo:vi],axis=0)*1.3*conv/u.cm**2
    
    mask = np.log10(zeroth.value)>thresh
    
    
#     plt.contourf(np.log10(zeroth.value)*mask,levels=np.linspace(threshs[i],22.5,30))
#     plt.gca().set_aspect(1)
#     plt.show()
    
    masked = zeroth*mask
    dA = (d*u.pc*0.125*u.deg.to(u.rad))**2
    area = np.nansum(masked/masked)*dA
    
    mass = np.nansum(masked*dA)*c.m_p*1.36*2
    return mass.to(u.solMass).value, np.sqrt(area/np.pi).value


#------------------------------------------------------------------------------------------------------------------
    
def line(x,m,b):
    return m*x+b

def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y
    
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
    
    
#------------------------------------------------------------------------------------------------------------------

def fit_plane2maps(i,j,rhi,vovi):
    m0,m1,m2 = moment_maps(i,threshs[i],rhi,vovi)[j]
    
    out = fit_gradient(m1,ds[i],weights=m2,full=True)
    
    pars = out['fit_results']
    uncs = out['1_sigma_errors']
    
    return pars, uncs

def spctrum_fit(i,which='hi'):
    coords = coordses[i]

#     sph2 = np.nanmean(h2cube,axis=(1,2)) #as things stand rn, we dont care abthe h2 fit
    if which=='hi':
        loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
        boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
        vo,vi = vs_px[i]
        
        hicube = hi4pi.data[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    
        spctr = np.nanmean(hicube,axis=(1,2))
        vx = (np.arange(vo,vi)*1.3-100)
        
    else:
        h2cube = get_subcube(dht.data,coords[:4]) #velocity axis length = 155
        spctr = np.nanmean(h2cube,axis=(1,2))[39:117] #-50km/s to 50km/s
        vx = (np.arange(39,117)*1.3-100)
    
    fit_sp = fttr(gaus,vx,spctr)
    
    A   = fit_sp.amplitude.value
    mu  = fit_sp.mean.value
    std = fit_sp.stddev.value
    
    return A,mu,std
    
    
#------------------------------------------------------------------------------------------------------------------

def get_theta(xy):
    x,y = xy
    thet = np.arctan(y/x)*u.rad.to(u.deg)
    
    if x < 0:
        thet += 180
        
    if thet > 180:
        thet -= 360
        
    return thet


def Omg_vec(i,rhi,vovi,unc=False,xy=False):
    fitpars_all = np.array([fit_plane2maps(i,j,rhi,vovi) for j in [0,1]])
    
    Omg_xys     = fitpars_all[:,0,:] #omgx, omgy values
    Omg_xys_unc = fitpars_all[:,1,:] #1sig uncertainties
    
    if xy:
        return Omg_xys

    omgh2,omghi = np.linalg.norm(Omg_xys[:,:2],axis=1)
    omgh2_unc, omghi_unc = np.linalg.norm(Omg_xys[:,:2]*Omg_xys_unc[:,:2],axis=1)/\
                            np.linalg.norm(Omg_xys[:,:2],axis=1)/np.sqrt(2)
    
    
    thth2 = get_theta(Omg_xys[0][:2])
    ththi = get_theta(Omg_xys[1][:2])
                      
    if unc:
        return omgh2,omghi,omgh2_unc,omghi_unc, thth2,ththi
                      
    else:
        return omgh2,omghi,thth2,ththi


hi4pi, dht = fits.open('../../data/hi4piDHT_720deg.fits')
table1 = pd.read_csv('../../data/molecular_clouds.csv')

coordses = table1[['l0 (°)','l1 (°)',
                   'b0 (°)','b1 (°)',
                  'v0 (km/s)','v1 (km/s)']].to_numpy().astype('float')               

ds    = table1['d (pc)'].to_numpy().astype('float')
duncs = np.linalg.norm(np.array([2,7,2,24,8,11,4,60,3,9,3,8,4,44,4,89,1,11,5,25,9,9,6,104,3,37,\
                  3,39,2,7,2,22,2,21,1,20,2,15,1,18,5,65,2,7]).reshape(22,2),axis=1) #from table 1 in Zucker+19
names = table1['Name']

lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')

# muh2_3sig     = [mu+np.array([-1,1])*2*std  for i in range(22) for mu,std in [spctrum_fit(i,which='h2')[1:]]]
# muh2_3sig[11] = [mu+np.array([-1,1])*std    for i in [11] for mu,std in [spctrum_fit(i,which='h2')[1:]]][0]
# muh2_3sig[19] = [mu+np.array([-1,1])*3*std  for i in [19] for mu,std in [spctrum_fit(i,which='h2')[1:]]][0]
# vs_px     = ((np.array(muh2_3sig)+100)/1.3).astype('int')


#----------------------------------------------------------------------------------------------------------

threshs = [20.5,21,20.5,21,20.5,20.5,20.5,21,20.5,20.5,21,21,21,21,21,21,21,20.5,21,20.5,21,21]

#----------------------------------------------------------------------------------------------------------


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


mh2,rh2,rhi = np.array([get_mh2_rh2_rhi(i,threshs[i]) for i in range(22)]).transpose()
rhi_px  = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')
mhi = np.array([calc_masses(i,rhi_px[i],vsnew[i])[1] for i in range(22)])


rh2unc = 2*rh2*duncs/ds
rhiunc = 2*rhi*duncs/ds


#-----------------------------------------------------------------------------------------------------------------------


#how i did it in untitled8-copy1 <- the notebook where the calculatons + figures happened for the resubmission
#but now doesnt work for some reason so just putting in values from table1
# vmnhi,sigvhi = np.array([spctrum_fit(i)[1:] for i in range(22)]).transpose()

vmnhi  = np.array([2.4,-0.5,0.2,17.2,1.2,5.9,11.8,6.4,5.0,0.9,2.7,27.6,8.8,9.5,2.0,8.3,8.8,6.7,2.5,-1.6,11.4,4.1])
sigvhi = np.array([3.0,11.6,3.9,6.3,3.9,4.4,7.1,9.9,5.6,5.3,5.6,12.5,7.3,6.5,4.4,6.5,6.2,7.3,5.8,9.4,9.6,8.4])
vsnew = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')
vsnew[11] = ((vmnhi[11]+np.array([-1,1])*sigvhi[11] + 100)/1.3).astype('int')


#-----------------------------------------------------------------------------------------------------------------------
omgh2,omghi,omgh2unc, omghiunc, thth2,ththi = np.array([Omg_vec(i,rhi_px[i],vsnew[i],unc=True) 
                                                        for i in range(22)]).transpose()


#-----------------------------------------------------------------------------------------------------------------------


jh2 = 0.4*(omgh2*u.m/u.s/u.pc*(rh2*u.pc)**2).cgs
jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs

jh2unc = np.sqrt((omgh2_unc/omgh2)**2 + (2*rh2unc/rh2)**2)*jh2
jhiunc = np.sqrt((omghi_unc/omghi)**2 + (2*rhiunc/rhi)**2)*jhi


#-----------------------------------------------------------------------------------------------------------------------

# tab2_newnew = pd.DataFrame(np.around(np.array((mh2/1e4,rh2,np.array(omgh2)/10,thth2,jh2/1e24,\
#                                   mhi/1e5,rhi,np.array(omghi)/10,ththi,jhi/1e24)),3).transpose(),
#                         columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


rhi_cyl = lambda mh2, mu=1.36, nhi=0.57/u.cm**3, mp=c.m_p : (mh2/(2*np.pi*mu*mp*nhi))**(1/3)
rhi_c = rhi_cyl(mh2*u.solMass).to(u.pc)

rhi_pxc = (np.arctan(rhi_c.value/ds)*u.rad.to(u.deg)/0.125).astype('int')
mhi_c = np.array([calc_masses(i,rhi_pxc[i],vsnew[i])[1] for i in range(22)])


def chebyshev_grid(arr, orgn=0): #Returns an grid of same dim as arr where each cell contains its Chebyshev distance to 'orgn'.
    if orgn == 0:
        orgn = (np.array(arr.shape)/2).astype('int')
        
    coords = np.ogrid[tuple(slice(0, s) for s in arr.shape)]
    dx,dy = [np.abs(coord - p) for coord, p in zip(coords, orgn)]
    
    chebishv = np.nanmax(np.broadcast_arrays(dx,dy),axis=0)
    return chebishv, np.unique(chebishv)


#---------------------------------------------------------------------------------------------------------------------------------

def add_mass(mhi,hi_dm,chbshv,mh2,n=0): #put in a seed value of mhi as whats in the central pixel
    if (mhi >= mh2) or (n == chbshv[1][-1]):#reached the end of :
        # plt.imshow(hi_dm*(chbshv[0] <= n))
        # plt.show()
        return n, mhi/mh2, chbshv[1][-1]#, n
    
    
    n   += 1
    mhi = np.nansum(hi_dm*(chbshv[0] <= n))
    
    return add_mass(mhi,hi_dm,chbshv,mh2,n)


#RESUME HERE 


[
[(plt.imshow(himaps[0]*(grid==n)),
 plt.show())
    for h2maps,himaps in [moment_maps(i,threshs[i],rhi_pxc[i],vsnew[i])]
    for grid,nums in [chebyshev_grid(himaps[0])]
   for n in nums]
# plt.show()
# )
    for i in [0]]#range(22)] 

# plt.show()


[(
plt.plot([np.nansum(himaps[0]/u.cm**2*(grid<=n)*(ds[i]*u.pc*0.125*u.deg.to(u.rad))**2*c.m_p*1.36).to(u.Msun).value
    for h2maps,himaps in [moment_maps(i,threshs[i],rhi_pxc[i],vsnew[i])]
    for grid,nums in [chebyshev_grid(himaps[0])]
   for n in nums]),
    

# plt.show()
) for i in range(22)] 

# plt.loglog(),
plt.yscale('log')

plt.show()


[(
plt.plot([np.nansum(h2maps[0]/u.cm**2*(grid<=n)*(ds[i]*u.pc*0.125*u.deg.to(u.rad))**2*c.m_p*2*1.36).to(u.Msun).value
    for h2maps,himaps in [moment_maps(i,threshs[i],rhi_pxc[i],vsnew[i])]
    for grid,nums in [chebyshev_grid(h2maps[0])]
   for n in nums]),
# plt.show()
) for i in range(22)] 

# plt.loglog(),
plt.yscale('log')

plt.show()


rhi_c


[(
plt.plot([np.nansum(himaps[0]/u.cm**2*(grid<=n)*(ds[i]*u.pc*0.125*u.deg.to(u.rad))**2*c.m_p*1.36).to(u.Msun).value#/mh2[i]
    for grid,nums in [chebyshev_grid(himaps[0])]
   for n in nums],label='hi'),
    
    
plt.plot([np.nansum(h2maps[0]/u.cm**2*(grid<=n)*(ds[i]*u.pc*0.125*u.deg.to(u.rad))**2*c.m_p*2*1.36).to(u.Msun).value#/mh2[i]
    for h2maps,himaps in [moment_maps(i,threshs[i],rhi_pxc[i],vsnew[i])]
    for grid,nums in [chebyshev_grid(h2maps[0])]
   for n in nums],label='h2'),
    
plt.title(names[i]),
# plt.loglog(),
plt.yscale('log'),
plt.xlabel('R [pc]'),
plt.ylabel('$\int M ds\ [M_\mathrm{H_2}]$'),
plt.xticks(np.linspace(0,chebyshev_grid(himaps[0])[1][-1],5),
           labels=np.around(np.linspace(0,chebyshev_grid(himaps[0])[1][-1],5)*0.125*u.deg.to(u.rad)*ds[i],1)),
    # chebyshev_grid(himaps[0])[1][::10],labels=np.around(chebyshev_grid(himaps[0])[1][::10]*0.125*u.deg.to(u.rad)*ds[i],1)),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i],rhi_pxc[i],vsnew[i])]] 


# plt.yscale('log')


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den


def plot_rvsv(i,rhi,vovi, show=False):  
    h2dvds, hidvds = Omg_vec(i,pixelize(rhi_sum[i],ds[i]),vsnew[i],xy=True)/1e3
    
    omgh2 = np.linalg.norm(h2dvds[:2])
    omghi = np.linalg.norm(hidvds[:2])
    thth2 = get_theta(h2dvds[:2])
    ththi = get_theta(hidvds[:2])
    
    l0c,l1c,b0c,b1c = coordses[i][:4]
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2maps, himaps = moment_maps(i,threshs[i],rhi,vovi)
    h2_0,h2_1,h2_2 = h2maps
    hi_0,hi_1,hi_2 = himaps
    
    l0h,l1h = lc + rhi*0.125, lc - rhi*0.125
    b0h,b1h = bc - rhi*0.125, bc + rhi*0.125
    
    
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

    
    hbperp = (hidvds[0]*hxs_pc + hidvds[1]*hys_pc)/np.linalg.norm(hidvds)
    cbperp = (h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc)/np.linalg.norm(h2dvds)

    
    hv_model =  hidvds[2] + hidvds[0]*hxs_pc + hidvds[1]*hys_pc
    cv_model =  h2dvds[2] + h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hi_2),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod',])


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


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0),
            np.ravel(h2_1),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','rper','vmod',])

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


    h2p, hip = pearson(meancdf['rper'],  meancdf['v']),pearson(meanhdf['rper'],  meanhdf['v'])
#     vmn = min(min(meancdf['v']),min(meancdf['vmod']),min(meanhdf['v']),min(meanhdf['vmod'])) 
#     vmx = max(max(meancdf['v']),max(meancdf['vmod']),max(meanhdf['v']),max(meanhdf['vmod'])) 
    
    
    if show:
        fig,axs = plt.subplots(1,2,figsize=(12,5))
        # fig.suptitle(names_full[i])
    
        cm,cb = curve_fit(line,meancdf['rper'],meancdf['v'])[0]
        csm   = np.sqrt(curve_fit(line,meancdf['rper'],meancdf['v'])[1][0,0])
    
        hm,hb = curve_fit(line,meanhdf['rper'],meanhdf['v'])[0]
        hs = np.sqrt(curve_fit(line,meanhdf['rper'],meanhdf['v'])[1][0,0])
    
        axs[0].fill_between(meancdf['rper'],
                            meancdf['v']+sigcdf['v'],
                            meancdf['v']-sigcdf['v'],
                            color=(0.,0.8,0.8,0.6),linewidth=0.5)
        
        axs[0].fill_between(meancdf['rper'],
                            meancdf['v']+2*sigcdf['v'],
                            meancdf['v']-2*sigcdf['v'],
                            color=(0.,0.8,0.8,0.3),linewidth=0.5)
    
     
        axs[0].scatter(meancdf['rper'],meancdf['v'],color='black',s=15)
        axs[0].plot(meancdf['rper'],meancdf['vmod'],color='black',linewidth=0.7,linestyle='--')
        
#         axs[0].scatter(cdf['rper'],cdf['v'])
        
    
        axs[1].fill_between(meanhdf['rper'],
                            meanhdf['v']+sighdf['v'],
                            meanhdf['v']-sighdf['v'],
                            color=(0.8,0.,0.,0.6),linewidth=0.5)
        
        axs[1].fill_between(meanhdf['rper'],
                            meanhdf['v']+2*sighdf['v'],
                            meanhdf['v']-2*sighdf['v'],
                            color=(0.8,0.,0.,0.3),linewidth=0.5)
    
        
        axs[1].scatter(meanhdf['rper'],meanhdf['v'],color='black',s=15)
        axs[1].plot(meanhdf['rper'],meanhdf['vmod'],color='black',linewidth=0.7,linestyle='--')
    
    
        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
            # ax.grid()
        axs[0].set_ylabel(r'$\langle v_\mathrm{LSR}\ \rangle\ [\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names[i])
        
        axs[0].set_title(np.around(np.abs(h2p),2),y=0.87,x=0.13)
        axs[1].set_title(np.around(hip,2),y=0.87,x=0.13)

#         axs[1].annotate(str(np.around(hip,2)),(meanhdf['rper'].iloc[0],
#                                                np.nanmax(meanhdf['v'])+np.nanmean(sighdf['v'])))
        
    
        # fig.savefig('v_vs_r_'+names[i],bbox_inches='tight')
        plt.show()

    return h2p, hip
#     return pearson(meancdf['rper'],  meancdf['v']),                   pearson(meanhdf['rper'],  meanhdf['v']),\
#            pearson(sigcdf['v'],      np.abs(meancdf['v']-h2_vos[i])), pearson(sighdf['v'],      np.abs(meanhdf['v']-hi_vos[i])),\
#            pearson(meancdf['w'],     np.abs(meancdf['v']-h2_vos[i])), pearson(meanhdf['w'],     np.abs(meanhdf['v']-hi_vos[i]))


[plot_rvsv(i,pixelize(rhi_sum[i],ds[i]),vsnew[i],show=True) for i in range(22)]


[plot_rvsv(i,pixelize(rhi_sum[i],ds[i]),vsnew[i],show=True) for i in [21]]


[(
plt.contourf(himaps[1],cmap='Spectral',levels=30),
plt.gca().set_aspect(1),
    plt.title(names[i]),
plt.show()
)for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i],pixelize(rhi[i],ds[i]),vsnew[i])]] 


[(
plt.contourf(himaps[1],cmap='Spectral',levels=30),
plt.gca().set_aspect(1),
    plt.title(names[i]),
plt.show()
)for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i],pixelize(rhi_sum[i],ds[i]),vsnew[i])]] 


for i in range(22):
    h2maps,himaps = moment_maps(i,20.5,rhi_pxc)
    himaps

    fig,axs=plt.subplots(1,3,dpi=200,figsize=(10,5))

    axs[0].imshow(h2maps[0],cmap='viridis',origin='lower')#,extent=coordses[i][:4]),
    
    oghi = axs[1].imshow(himaps[0],cmap='magma',origin='lower',)#extent=coordses[i][:4]),
    
    cmin,cmax = oghi.get_clim()
    hi0_h2lb = np.nansum(get_subcube(hi4pi.data[vsnew[i][0]:vsnew[i][1]],coordses[i][:4]),axis=0)
    axs[2].imshow(hi0_h2lb,cmap='magma',origin='lower')#,vmin=cmin,vmax=cmax)
    axs[2].contour(h2maps[0],levels=5)

    for ax in axs:
        ax.set_aspect(1),
        # ax.set_xlim(coordses[i][:2]),
        
    plt.title(names[i]),
    # plt.colorbar(),
    

    plt.show()


summ = np.array([add_mass(0,himaps[0],chebyshev_grid(himaps[0]),np.nansum(h2maps[0])) 
        for i in range(22) 
        for h2maps,himaps in [moment_maps(i,threshs[i],rhi_pxc[i],vsnew[i])]])

rhi_sum = np.array(summ[:,0])*0.125*u.deg.to(u.rad)*ds

mhi_sum = summ[:,1]*mh2


plt.figure(dpi=200)
ss = 50

plt.scatter(rh2,rhi,label='sphere, fixed n_wnm',edgecolor='white',s=ss)
plt.scatter(rh2,rhi_c,label='cylinder, fixed n_wnm',edgecolor='white',s=ss)
plt.scatter(rh2,rhi_cc,label='cylinder, fixed n_cnm',edgecolor='white',s=ss)

plt.scatter(rh2,rhi_sum,label='evolving sum',edgecolor='white',s=ss)


plt.plot(rh2,rh2,color='black',label='rhi = rh2')


# plt.vlines(rh2,rhi,rhi_c,color='black',zorder=-10,linewidth=0.5)


plt.legend(fontsize='x-small')
# plt.xticks(np.arange(22),labels=names,rotation='vertical')

plt.ylabel(r'$R\mathrm{_{HI}\ [pc}]$')
plt.xlabel(r'$R\mathrm{_{H_2}\ [pc}]$')

plt.loglog()

# plt.grid()
plt.show()


plt.figure(dpi=200)
ss = 50
plt.scatter(mh2,mhi,label='sphere, fixed n_wnm',s=ss,alpha=0.5)
plt.scatter(mh2,mhi_c,label='cylinder, fixed n_wnm',s=ss,alpha=0.5)
plt.scatter(mh2,mhi_cc,label='cylinder, fixed n_cnm',s=ss,alpha=0.5)

plt.scatter(mh2,mhi_sum,label='evolving sum',s=ss,alpha=0.5)

# plt.scat(mh2,mhi_h2lb,mhi,color='black',zorder=-10,linewidth=0.5)
plt.scatter(mh2,mhi_h2lb,label='inside co',s=ss,alpha=0.5)

plt.vlines(mh2,mhi_cc,mhi,color='black',zorder=-10,linewidth=0.5)
plt.plot(mh2,mh2,color='black',)

# plt.plot(np.arange(22),mh2,label='mh2',zorder=10,color='black')
# plt.plot(np.arange(22),10*mh2,label='10 x mh2',zorder=10,linestyle='--',color='black')


# plt.xticks(np.arange(22),labels=names,rotation='vertical')
# plt.grid()

plt.legend(fontsize='x-small')

plt.ylabel(r'$M_\mathrm{HI}\ [M_\odot]$')
plt.xlabel(r'$M_\mathrm{H_2}\ [M_\odot]$')

# [plt.annotate(names[i], (mh2[i],mhi_h2lb[i])) for i in range(22) if mh2[i] > mhi_h2lb[i]] 
# [plt.annotate(names[i], (mh2[i],mhi[i])) for i in range(22) if mhi_c[i] < mhi_h2lb[i]] 

# [plt.annotate(names[i], (mh2[i],mhi_h2lb[i])) for i in range(22)]

plt.loglog()

plt.show()


[(
plt.plot(np.nanmean(get_subcube(dht.data,coordses[i][:4]),axis=(1,2)),marker='.'),
# plt.plot(1e-2*np.nanmean(get_subcube(hi4pi.data,coordses[i][:4]),axis=(1,2))),

plt.vlines((coordses[i][4:]+100)/1.3,
           0,np.nanmax(np.nanmean(get_subcube(dht.data,coordses[i][:4]),axis=(1,2)))),
# plt.vlines(vsnew[i],0,np.nanmax(np.nanmean(get_subcube(dht.data,coordses[i][:4]),axis=(1,2)))),

    plt.show()
) for i in range(22)]


[(
plt.contourf(h1,cmap='Spectral',levels=30,#np.linspace(vsnew[i][0]*1.3-100,vsnew[i][1]*1.3-100,30),
                                                   extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
        plt.xlim(coordses[i][:2]),

plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,20.5,rhi_pxc[i],vsnew[i])] for h0,h1,h2 in [h2maps]]


# # post meeting notes on 12/19
# ### the hi masses within a small region (eg the gmc footprint) are much higher than expected,given the fact that NHI < NH2
# - (from grid of the npdfs) turns out that this is only true for ~half the clouds: aqs,cam,cma,chm,cra,her,mob1,mr2,oph,ora,orb,per,pol --- weirdly tho california doesnt follow the expectation even tho in ib16 ...
# - expectation that NHI always < NH2: possibly bc high densities -> brightest in co -> were the first to be observed, have had a lot of follow-up -> are assumed to be the standard/benchmark 
# - when summing over the whole velocity range, the clouds where its still NH2 > NHI are cra, oph, ora, orb, pers, pol
# - aqs,xbn,her,lup,mad,ros have NHI > NH2 bedenb -> blending?
# 
# 
# 
# 
# ### cumulative mass function (eg lada lombardi alves, fig1) + other independent tests to check that the method for stepping out in some dr and calculating mhi is actually doing the correct thing
# - cumulative mass profiles - the ones where the blue and orange cross early = small rhi is sufficient for mhi=mh2
# - the clouds that cross early / dont correspond to the ones with small rhi -> really dense hi
# 
# 
# 
# 
# ### also look at the velocity gradients within the different rs and see what effect that has on the js
# - there is variation, but not large enough to change the overall population mean (ie omghi < omgh2 for most of the clouds regardless of rhi) 
# - the times when omghi > omgh2, it is usualy for the smaller rhi, the smaller rhi are also the ones with larger uncertainties (comparatively), but still ~<1% of respective omega
# - taurus and orion b are the only ones with all their omghi > omgh2


np.array(Omg_vec(i,pixelize(rhi_sum[i],ds[i]),vsnew[i]))[[1,3]]


Omg_vec(i,pixelize(rhi_sum[i],ds[i]),vsnew[i],unc=True)


pixelize = lambda rhi, d: (np.arctan(rhi/d)*u.rad.to(u.deg)/0.125).astype('int') 

omghi_sum, omghi_sumunc, ththi_sum,  = np.array([np.array(Omg_vec(i,pixelize(rhi_sum[i],ds[i]),vsnew[i],unc=True))[[1,3,5]] for i in range(22)]).transpose()
omghi_c, omghi_cunc, ththi_c = np.array([np.array(Omg_vec(i,rhi_pxc[i],vsnew[i],unc=True))[[1,3,5]] for i in range(22)]).transpose()
omghi_cc,omghi_ccunc, ththicc  = np.array([np.array(Omg_vec(i,rhi_pxcc[i],vsnew[i],unc=True))[[1,3,5]] for i in range(22)]).transpose()


mom2h2,mom2hi = [[np.nanmean(moment_maps(i,20.,pixelize(rhi[i],ds[i]),vsnew[i])[j][2]) for i in range(22)] for j in [0,1]]


plt.scatter(mom2h2,omgh2/1000)
# plt.plot(omgh2/1000,omgh2/1000)


plt.scatter(mom2hi,omghi/1000)
# plt.plot(omghi/1000,omghi/1000)

plt.ylabel('$\Omega\ [\mathrm{km\ s^{-1}\ pc^{-1}}$]')
plt.xlabel('$\sigma$')


plt.scatter(mom2hi,sigvhi)
plt.plot(np.arange(2,10),np.arange(2,10))
plt.ylabel('mean then std')
plt.xlabel('std then mean')


help(plt.errorbar)


omghi_cunc


omghi_sumunc


plt.figure(dpi=200)
ss=20
plt.scatter(names,omghi_unc,label='sphere',s=ss)
plt.scatter(names,omghi_cunc,label='cylinder',s=ss)
plt.scatter(names,omghi_ccunc,label='cnm',s=ss)
plt.scatter(names,omghi_sumunc,label='sum',s=ss)
plt.legend()
plt.xticks(rotation='vertical')
plt.show()


plt.figure(dpi=200)
ss = 50

plt.scatter(omgh2,omghi,label='sphere, fixed n_wnm',edgecolor='white',s=ss*4)
plt.scatter(omgh2,omghi_c,label='cylinder, fixed n_wnm',edgecolor='white',s=ss*3)
# plt.errorbar(omgh2,omghi_cc,yerr=omghi_ccunc,linestyle='',zorder=30)

plt.scatter(omgh2,omghi_cc,label='cylinder, fixed n_cnm',edgecolor='white',s=ss*2)

plt.scatter(omgh2,omghi_sum,label='evolving sum',edgecolor='white',s=ss)
plt.plot(omgh2,omgh2,color='black',label='rhi = rh2')


# plt.scatter(names,omghi,label='sphere, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(names,omghi_c,label='cylinder, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(names,omghi_sum,label='evolving sum',edgecolor='white',s=ss)
# plt.scatter(names,omghi_cc,label='cylinder, fixed n_cnm',edgecolor='white',s=ss)


# plt.vlines(rh2,rhi,rhi_c,color='black',zorder=-10,linewidth=0.5)


plt.legend(fontsize='x-small')
# plt.xticks(np.arange(22),labels=names,rotation='vertical')

plt.ylabel(r'$\Omega\mathrm{_{HI}\ [m\ s^{-1}\ pc^{-1}}]$')
plt.xlabel(r'$\Omega\mathrm{_{H_2}\ [m\ s^{-1}\ pc^{-1}}]$')

[plt.annotate(names[i],(omgh2[i],omghi_sum[i])) for i in range(22) if omghi[i]>omgh2[i]]

plt.loglog()

# plt.grid()
plt.show()


plt.figure(dpi=200)
ss = 50

# plt.scatter(omgh2,omghi,label='sphere, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(omgh2,omghi_c,label='cylinder, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(omgh2,omghi_cc,label='cylinder, fixed n_cnm',edgecolor='white',s=ss)

# plt.scatter(omgh2,omghi_sum,label='evolving sum',edgecolor='white',s=ss)
# plt.plot(omgh2,omgh2,color='black',label='rhi = rh2')


plt.scatter(names,omghi,label='sphere, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(names,omghi_c,label='cylinder, fixed n_wnm',edgecolor='white',s=ss)
plt.scatter(names,omghi_sum,label='evolving sum',edgecolor='white',s=ss)
# plt.scatter(names,omghi_cc,label='cylinder, fixed n_cnm',edgecolor='white',s=ss)


# plt.vlines(rh2,rhi,rhi_c,color='black',zorder=-10,linewidth=0.5)


plt.legend(fontsize='x-small')
plt.xticks(np.arange(22),labels=names,rotation='vertical')

plt.ylabel(r'$\Omega\mathrm{_{HI}\ [m\ s^{-1}\ pc^{-1}}]$')
# plt.xlabel(r'$\Omega\mathrm{_{H_2}\ [m\ s^{-1}\ pc^{-1}}]$')

# plt.loglog()

plt.grid()
plt.show()


# make maps to see if the gradients within the smaller radii are linear (or make the rperp vs delta v plots)


plt.figure(dpi=200)
ss = 50

plt.scatter(names,ththi,label='sphere, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(names,ththi_c,label='cylinder, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(names,ththicc,label='cylinder, fixed n_cnm',edgecolor='white',s=ss)

plt.scatter(names,ththi_sum,label='evolving sum',edgecolor='white',s=ss)
# plt.plot(thth2,thth2,color='black',label='rhi = rh2')


# plt.vlines(rh2,rhi,rhi_c,color='black',zorder=-10,linewidth=0.5)


plt.legend(fontsize='x-small')
plt.xticks(np.arange(22),labels=names,rotation='vertical')

plt.ylabel(r'$\theta\mathrm{_{HI}\ [\degree}]$')
# plt.xlabel(r'$\Omega\mathrm{_{H_2}\ [m\ s^{-1}\ pc^{-1}}]$')

# plt.loglog()

plt.grid()
plt.show()


get_j = lambda omg,r : (0.4*(omg*u.m/u.s/u.pc)*r**2).cgs


(get_j(omghi_sum,rhi_sum*u.pc)/(omghi_sum*u.m/u.s/u.pc*rhi_sum*u.pc)**2).cgs


(omghi_c*u.m/u.s/u.pc*rhi_c)


plt.figure(dpi=200)#,figsize=(6,6))
plt.scatter(names,(get_j(omghi_c,rhi_c)/(omghi_c*u.m/u.s/u.pc*rhi_c)**2).cgs.to(u.Myr),label='just erot')
# plt.scatter(names,(get_j(omghi_c,rhi_c)/\
#                    (-c.G*mhi_c*u.Msun/rhi_c + \
#                     (sigvhi*u.km/u.s)**2 +\
#                     (omghi_c*u.m/u.s/u.pc*rhi_c)**2)).cgs.to(u.Myr),label='-egrav + erot + ekin')
plt.scatter(names,rhi/sigvhi,label='hi crossing time')
# plt.scatter(names,rh2/np.nanmean(mom2h2),label='h2 crossing time')


plt.scatter(names,(2*np.pi/(omghi_c*u.m/u.s/u.pc)).cgs.to(u.Myr),label='hi rotation time')
# plt.scatter(names,(2*np.pi/(omgh2*u.m/u.s/u.pc)).cgs.to(u.Myr),label='h2 rotation time')


plt.yscale('log')
# plt.scatter(names,(get_j(omghi_sum,rhi_sum*u.pc)/(omghi_sum*u.m/u.s/u.pc*rhi_sum*u.pc)**2).cgs.to(u.Myr),label='just erot')
# plt.scatter(names,(get_j(omghi_sum,rhi_sum*u.pc)/\
#                    (omghi_sum*u.m/u.s/u.pc*rhi_sum*u.pc)**2).cgs.to(u.Myr),label='just erot')


plt.xticks(rotation='vertical')
plt.grid()
plt.legend()
plt.show()


plt.figure(dpi=200)
ss = 50

plt.scatter(jh2.value,get_j(omghi,rhi*u.pc),label='sphere, fixed n_wnm',edgecolor='white',s=ss)
# plt.scatter(jh2.value,get_j(omghi_c,rhi_c),label='cylinder, fixed n_wnm',edgecolor='white',s=ss)
plt.scatter(jh2.value,get_j(omghi_sum,rhi_sum*u.pc),label='evolving sum',edgecolor='white',s=ss)
# plt.scatter(jh2.value,get_j(omghi_cc,rhi_cc),label='cylinder, fixed n_cnm',edgecolor='white',s=ss)

plt.plot(jh2.value,jh2.value,color='black',label='jhi = jh2')


# plt.vlines(rh2,rhi,rhi_c,color='black',zorder=-10,linewidth=0.5)


plt.legend(fontsize='x-small')
# plt.xticks(np.arange(22),labels=names,rotation='vertical')

plt.ylabel(r'$j\mathrm{_{HI}\ [cm^2\ s^{-1}}]$')
plt.xlabel(r'$j\mathrm{_{H_2}\ [cm^2\ s^{-1}}]$')

plt.loglog()

# plt.grid()
plt.show()


curve_fit(line,np.log10(mhi),np.log10(jhi.value))


np.sqrt(0.00590627)


curve_fit(line,np.log10(rhi),np.log10(jhi.value))


np.sqrt(0.05350325)


help(moment_maps)


fig,axs = plt.subplots(4,6,figsize=(20,15),dpi=200)

[(
    # plt.imshow(h2maps[0],origin='lower'),
 # plt.imshow(np.log10(himaps[0]),origin='lower'),
 axs.flat[i].hist(np.ravel(np.log10(himaps[0])),bins=70,density=True,alpha=0.5),#bins=np.arange(18,22,0.1)),
 axs.flat[i].hist(np.ravel(np.log10(h2maps[0])),bins=70,density=True,alpha=0.5),
axs.flat[i].set_title(names[i]),
axs.flat[i].set_yscale('log'),
 # plt.colorbar(),
 # plt.show()
) for i in range(22)
  for h2maps,himaps 
     in [moment_maps(i,20.,pixelize(rhi_c[i].value,ds[i]),vsnew[i])]]

plt.show()


fig,axs = plt.subplots(4,6,figsize=(20,15),dpi=200)

[(
    # plt.imshow(h2maps[0],origin='lower'),
 # plt.imshow(np.log10(himaps[0]),origin='lower'),
 axs.flat[i].hist(np.ravel(np.log10(himaps[0])),bins=70,density=True,alpha=0.5),#bins=np.arange(18,22,0.1)),
 axs.flat[i].hist(np.ravel(np.log10(h2maps[0])),bins=70,density=True,alpha=0.5),
axs.flat[i].set_title(names[i]),
axs.flat[i].set_yscale('log'),
 # plt.colorbar(),
 # plt.show()
) for i in range(22)
  for h2maps,himaps 
     in [moment_maps(i,20.,pixelize(rhi_c[i].value,ds[i]),[0,155])]]

plt.show()


gaia_all = pd.read_csv('../../data/gaia_stars_2000pc_pm50kms.csv')


gaia_all


plt.scatter(gaia_all['l'],gaia_all['b'])
plt.scatter(gaia_all['l']-360,gaia_all['b'])
plt.xlim(180,-180)


[[(
plt.plot(np.arange(155)*1.3-100,np.nanmean(h2,axis=(1,2))),
plt.plot(np.arange(155)*1.3-100,np.nanmean(hi,axis=(1,2))),
plt.show()
)for h2, hi in [make_h2_hi_cubes(i,rhi_pxc[i],[0,155])]]

for i in range(22)]


[(
plt.contourf(chmap,),
# plt.contour(np.nansum(h2,axis=0)),
    plt.colorbar(),
plt.show()

) 
for h2, hi in [make_h2_hi_cubes(20,int(2.5*rhi_pxc[20]),vsnew[20])]
for chmap in hi]


import tensorly as tly
from tensorly.decomposition import parafac


weights


component_3d = np.einsum('i,j,k->ijk',factors[0][:,0],factors[1][:,0],factors[2][:,0])


plt.contourf(np.nanmean(tautens,axis=2))


plt.contourf(np.nanmean(tautens,axis=1))


plt.contourf(np.nansum(component_3d,axis=0))


for i in [21]:
    tens = tly.tensor(make_h2_hi_cubes(i,int(1.5*rhi_pxc[i]),vsnew[i])[1])
    weights,factors = parafac(tens, rank=5)
    
    plt.plot(np.nanmean(tens,axis=(1,2)))
    plt.plot(factors[0])
    plt.show()


[(
plt.contourf(chmap,vmax=100),
# plt.contour(np.nansum(h2,axis=0)),
    plt.colorbar(),
plt.show()

) 
for h2, hi in [make_h2_hi_cubes(i,int(1.5*rhi_pxc[i]),vsnew[i])]
for chmap in hi]


plt.rcParams['font.size'] = 23
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10),sharey=True)

jbns=10**np.arange(22,27,1/5)
thbns=np.arange(-180,181,90/4)

axs[0,0].hist(jh2.value,bins=jbns,color=(0.,0.8,0.8),zorder=15,edgecolor='white')
axs[1,0].hist(jhi.value,bins=jbns,color=(0.8,0.,0.),zorder=15,edgecolor='white')

for ax in axs[:,0]:
    ax.set_xscale('log')
    ax.set_xlim(1e23,2e26)
    ax.set_ylabel('N (clouds)')
    ax.set_ylim(0,6.1)
    
for ax in axs[:,1]:
    ax.set_xticks(np.arange(-180,181,90))
    ax.set_xlim(-185,185)
    ax.vlines([-90,90],0,6,linestyle='--',color='black',zorder=-20)
    ax.set_ylim(0,6.1)
    # ax.set_yticks([1,3,5])
    

axs[0,1].hist(thth2,bins=thbns,edgecolor='white',color=(0.,0.8,0.8))
axs[1,1].hist(ththi,bins=thbns,edgecolor='white',color=(0.8,0.,0.))


axs[0,0].set_xlabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')
axs[1,0].set_xlabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
    
axs[0,1].set_xlabel(r'$\theta_\mathrm{H_2}\ [\degree$]')
axs[1,1].set_xlabel(r'$\theta_\mathrm{HI}\ [\degree$]')
    

fig.tight_layout()

# plt.ylabel('N (clouds)')
# plt.xlim(8e22,3e26)
# plt.ylim(0,5.5)
# # plt.yticks(np.arange(5))

# plt.savefig('fig7_js',bbox_inches='tight')

# plt.show()


curve_fit(line,np.log10(jh2.value),np.log10(jhi.value),sigma=np.log10(jhiunc.value))


np.sqrt(0.01538317)


plt.rcParams['font.size']=18
plt.figure(dpi=200,figsize=(5,5))

xs = np.linspace(5e22,5e26,5)

plt.fill_between(xs,xs,10*xs,color='grey',zorder=-10,alpha=0.2,linewidth=0.)
plt.plot(xs,10**10.6*xs**0.60,color='red', linestyle='--',zorder=-9,
         label=r'$j_\mathrm{HI}\ \propto\ j_\mathrm{H_2}^{0.60\pm0.12}$')

plt.errorbar(jh2.value,jhi.value,xerr=jh2unc.value,yerr=jhiunc.value,linestyle='',zorder=-5,color='black',capsize=2,alpha=0.5)

plt.legend(fontsize='x-small')

plt.scatter(jh2.value,jhi.value,
            color='black',edgecolor='white',s=100)


plt.xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel(r'$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}]$')

# plt.yticks([1,10],labels=['1','10'])
plt.xlim(1e23,3e25)
plt.ylim(1e24,2e26)

# plt.hlines([1,10],5e22,1e27,linestyle='-.',color='black')
# plt.ylim(5e-1,5e1)

plt.loglog()
plt.savefig('fig8_jj',bbox_inches='tight')


curve_fit(line,np.log10(jh2.value),np.log10(jhi.value))


np.sqrt(0.01551816)


help(curve_fit)


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))

xs = np.linspace(1,200,10)


axs[0].scatter(tab2_newnew['rh2'],tab2_newnew['jh2']*1e24,color=(0.,0.8,0.8),
            s=100,zorder=30,edgecolor='white')

axs[0].scatter(tab2_newnew['rh2'][h2p<0.75],
               tab2_newnew['jh2'][h2p<0.75]*1e24,
               color='white',
            s=100,zorder=30,edgecolor=(0.,0.8,0.8))

axs[0].set_xlim(3,50)
axs[0].set_ylim(1e23,1e26)
axs[0].set_xlabel(r'$R_\mathrm{H_2}$ [pc]')
axs[0].set_ylabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')

    
axs[0].plot(xs, 4e22*xs**1.67,color='black',linestyle='--',
            label='$j_\mathrm{H_2} \propto R^{1.67 \pm 0.23}$')
axs[0].fill_between(xs,4e22*xs**(1.67+0.23),3e22*xs**(1.67-0.23),
                 alpha=0.2,color='black',linewidth=0.)


axs[1].scatter(tab2_newnew['rhi'],tab2_newnew['jhi']*1e24,color=(0.8,0.0,0.),
            s=100,zorder=30,edgecolor='white')


axs[1].scatter(tab2_newnew['rhi'][hip<0.75],
               tab2_newnew['jhi'][hip<0.75]*1e24,
               color='white',
            s=100,zorder=30,edgecolor=(0.8,0.,0.))

axs[1].set_xlim(10,200)


axs[1].plot(xs, 1e23*xs**1.36,color='black',linestyle='--',
            label='$j_\mathrm{HI} \propto R^{1.36 \pm 0.39}$')
axs[1].fill_between(xs,1e23*xs**(1.36+0.39),1e23*xs**(1.36-0.39),
                 alpha=0.2,color='black',linewidth=0.)


axs[1].set_xlabel(r'$R_\mathrm{HI}$ [pc]')
axs[1].set_ylabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')

axs[0].legend(loc='upper left')
axs[1].legend(loc='lower right')
for ax in axs:
    ax.loglog()
    
    ax.set_yticks([1e23,1e24,1e25,1e26])
    ax.set_ylim(1e23,2e26)
    
fig.tight_layout()


plt.savefig('fig9_jRs',bbox_inches='tight')
