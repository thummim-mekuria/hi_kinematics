#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib import cm

from astropy.io import fits
import astropy.units as u
import astropy.constants as c

from scipy.optimize import curve_fit
from joblib import Parallel, delayed


from astropy.modeling import models,fitting
gaus = models.Gaussian1D(10.,5.,20.)
fttr = fitting.TRFLSQFitter()


def get_subcube(bigcube, coords, dx=-0.125, dy=0.125, x0=360, y0=-50):
    lo,li,bo,bi = coords # assume in degrees
    xo,xi = ((np.array([lo,li]) - x0)/dx).astype('int')
    yo,yi = ((np.array([bo,bi]) - y0)/dy).astype('int')
    
    return bigcube[:,yo:yi,xo:xi]


def make_h2_hi_cubes(i,rhi,vovi,hidat,h2dat):
    coords = coordses[i]
    vo,vi  = vovi

    loh,lih = lc_p[i]+np.array([-1,1])*rhi
    boh,bih = bc_p[i]+np.array([-1,1])*rhi

    hicube = hidat[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(h2dat,coords[:4]) #velocity axis length = 155
    
    return h2cube,hicube


def moment_maps(i,thresh, rhi, vovi, hidat, h2dat):
    h2cube,hicube = make_h2_hi_cubes(i,rhi,vovi,hidat,h2dat)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.29*1.82e18
    
    vo,vi = vovi

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.29-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.29-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.log10(hilb)>20)
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


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


def calc_masses(i,rhi,vovi,hidat,h2dat,dx): #cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    h20,hi0 = [moment_maps(i,threshs[i],rhi,vovi,hidat,h2dat)[j][0]/u.cm**2 for j in [0,1]]
    
    dA = (ds[i]*u.pc*dx*u.deg.to(u.rad))**2
    h2m = (np.nansum(h20*dA)*c.m_p*1.36*2).to(u.solMass)
    him = (np.nansum(hi0*dA)*c.m_p*1.36).to(u.solMass)
    
    return h2m.value,him.value


def get_mh2_rh2_rhi(i,thresh,norm=0.5):
    coords = coordses[i]
#     print(coords)
    h2cube = get_subcube(dht.data,coords[:4])
    mh2,rh2 = calc_h2_mass(h2cube,ds[i],thresh=thresh)

    rhi = norm*(((mh2*u.solMass)/(np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
    return mh2,rh2,rhi


def calc_h2_mass(cube,d,thresh=20.6,vo=0,vi=155,conv=1.97e20):
    zeroth = np.nansum(cube[vo:vi],axis=0)*1.3*conv/u.cm**2
    
    mask = np.log10(zeroth.value)>thresh

    masked = zeroth*mask
    dA = (d*u.pc*0.125*u.deg.to(u.rad))**2
    area = np.nansum(masked/masked)*dA
    
    mass = np.nansum(masked*dA)*c.m_p*1.36*2
    return mass.to(u.solMass).value, np.sqrt(area/np.pi).value


def line(x,m,b):
    return m*x+b

def plane_rvld(xy,a,b,vo):
    x,y = xy
    xo = np.nanmean(np.unique(x))
    yo = np.nanmean(np.unique(y))
#     return vo + a*(x-xo) + b*(y-yo)
    return vo + a*x + b*y


def fit_gradient(mom1_map, dist,dx, mask=0, weights=0, full=False):
    if type(mask)==int:#all the data is usable and/or a mask wasnt supplied
        mask = np.ones(mom1_map.shape).astype('bool')
        
    if type(weights)==int:#all the data is usable and/or a second moment map wasnt supplied
        weights = np.ones(mom1_map.shape).astype('bool')


    #the plane is pinned down to the geometric center of the map
    #the slope of the fit doesnt change with or without xo
    xdeg = (np.indices(mom1_map.shape)[1]-(mom1_map.shape[1]/2))*-dx
    ydeg = (np.indices(mom1_map.shape)[0]-(mom1_map.shape[0]/2))*dx

    
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
    
    
def fit_plane2maps(i,j,rhi,vovi,hidat,h2dat,dx):
    m0,m1,m2 = moment_maps(i,threshs[i],rhi,vovi,hidat,h2dat)[j]
    
    out = fit_gradient(m1,ds[i],dx,weights=m2,full=True)
    
    pars = out['fit_results']
    uncs = out['1_sigma_errors']
    
    return pars, uncs


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den


def spctrum_fit(i,hidat,rhi=0,which='hi'):
    coords = coordses[i]

#     sph2 = np.nanmean(h2cube,axis=(1,2)) #as things stand rn, we dont care abthe h2 fit
    if which=='hi':
        vmnh2, sigvh2 = np.array(spctrum_fit(i,hidat,which='h2')[1:]).transpose()
        vo,vi = ((vmnh2+3*np.array([-1,1])*sigvh2+100)/1.3).astype('int')
                
        loh,lih = lc_p[i]+np.array([-1,1])*rhi
        boh,bih = bc_p[i]+np.array([-1,1])*rhi
        
        hicube = hidat[vo:vi, boh:bih, loh:lih] #velocity axis length = vi-vo
    
        # plt.contourf(np.nansum(hicube,axis=0),levels=30)
        # plt.gca().set_aspect(1)
        # plt.show()
        
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


fit_plane2maps


def get_theta(xy):
    x,y = xy
    thet = np.arctan(y/x)*u.rad.to(u.deg)
    
    if x < 0:
        thet += 180
        
    if thet > 180:
        thet -= 360
        
    return thet


def Omg_vec(i,rhi,vovi,hidat,h2dat,dx,unc=False,xy=False):
    fitpars_all = np.array([fit_plane2maps(i,j,rhi,vovi,hidat,h2dat,dx) for j in [0,1]])
    
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


pixelize = lambda rhi, d, dx: (np.arctan(rhi/d)*u.rad.to(u.deg)/dx).astype('int') 


hi4pi, dht = fits.open('../../data/hi4piDHT_720deg.fits')
table1 = pd.read_csv('../../data/mi_tab1.csv')

coordses = table1[['l0(°)','l1(°)',
                   'b0(°)','b1(°)',
                  'v0(km/s)','v1(km/s)']].to_numpy().astype('float')               

ds    = table1['d(pc)'].to_numpy().astype('float')
duncs = np.linalg.norm(np.array([2,7,2,24,8,11,4,60,3,9,3,8,4,44,4,89,1,11,5,25,9,9,6,104,3,37,\
                  3,39,2,7,2,22,2,21,1,20,2,15,1,18,5,65,2,7]).reshape(22,2),axis=1) #from table 1 in Zucker+19
names = table1['Name']

lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')

lc_pp = ((lc-360)/-0.083).astype('int')
bc_pp = ((bc+50)/0.083).astype('int')

#----------------------------------------------------------------------------------------------------------

threshs = [20.5,21,20.5,21,20.5,20.5,20.5,21,20.5,20.5,21,21,21,21,21,21,21,20.5,21,20.5,21,21]

#----------------------------------------------------------------------------------------------------------

vmnh2, sigvh2 = np.array([spctrum_fit(i,hidat,which='h2')[1:] for i in range(22)]).transpose()
vsh2 = ((vmnh2[:,np.newaxis]+3*np.array([-1,1])[np.newaxis,:]*sigvh2[:,np.newaxis] +100)/1.3).astype('int')

mh2,rh2,rhi = np.array([get_mh2_rh2_rhi(i,threshs[i]) for i in range(22)]).transpose()

rh2unc = rh2*duncs/ds
rhiunc = (2/3)*rhi*duncs/ds


Omg_vec


hidat = hi4pi.data
h2dat = dht.data

lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')

dx = 0.125

rhi_px  = (np.arctan(rhi/ds)*u.rad.to(u.deg)/dx).astype('int')


vmnhi, sigvhi = np.array([spctrum_fit(i,hidat,rhi_px[i],which='hi')[1:] for i in range(22)]).transpose()
vshi_px = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')

vshi_px[3]  = ((vmnhi[3] +np.array([-1,1])*sigvhi[3]  + 100)/1.3).astype('int')
vshi_px[11] = ((vmnhi[11]+np.array([-1,1])*sigvhi[11] + 100)/1.3).astype('int')
vshi_px[12] = ((vmnhi[12]+np.array([-1,1])*sigvhi[12] + 100)/1.3).astype('int')
vshi_px[20] = ((vmnhi[20]+np.array([-1,1])*sigvhi[20] + 100)/1.3).astype('int')

omgh2,omghi,omgh2unc,omghiunc,thth2,ththi = np.array([Omg_vec(i,rhi_px[i],vshi_px[i],hidat,h2dat,dx,unc=True) 
                                                        for i in range(22)]).transpose()

mhi = np.array([calc_masses(i,rhi_px[i],vshi_px[i],hidat,h2dat,dx)[1] for i in range(22)])


HI4PI = fits.open('../../data/HI4PI.fits')[0].data
HI4PI = np.concatenate((HI4PI[:,:,2162:],HI4PI,HI4PI[:,:,:2162]),axis=2)[78:78+155,:,:]


(np.arange(310)*1.288 - 199.572)[78+155]


HI4PI


rHI_px/rhi_px


hidat = HI4PI
dx = 0.083

lc_p = ((lc-360)/-0.083).astype('int')
bc_p = ((bc+50)/0.083).astype('int')


rHI_px  = (np.arctan(rhi/ds)*u.rad.to(u.deg)/dx).astype('int')


vmnHI, sigvHI = np.array([spctrum_fit(i,hidat,rHI_px[i],which='hi')[1:] for i in range(22)]).transpose()
vsHI_px = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvHI[:,np.newaxis] + 100)/1.3).astype('int')

vsHI_px[3]  = ((vmnHI[3] +np.array([-1,1])*sigvHI[3]  + 100)/1.3).astype('int')
vsHI_px[11] = ((vmnHI[11]+np.array([-1,1])*sigvHI[11] + 100)/1.3).astype('int')
vsHI_px[12] = ((vmnHI[12]+np.array([-1,1])*sigvHI[12] + 100)/1.3).astype('int')
vsHI_px[20] = ((vmnHI[20]+np.array([-1,1])*sigvHI[20] + 100)/1.3).astype('int')


omgh2,omgHI,omgh2unc,omgHIunc,thth2,thtHI = np.array([Omg_vec(i,rHI_px[i],vsHI_px[i],hidat,h2dat,dx,unc=True) 
                                                        for i in range(22)]).transpose()

mHI = np.array([calc_masses(i,rHI_px[i],vsHI_px[i],hidat,h2dat,dx)[1] for i in range(22)])


jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs
jHI = 0.4*(omgHI*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs

jHIunc = np.sqrt((omgHIunc/omgHI)**2 + (2*rhiunc/rhi)**2)*jHI
jhiunc = np.sqrt((omghiunc/omghi)**2 + (2*rhiunc/rhi)**2)*jhi


mhi/mHI


plt.scatter(names,mhi)
plt.scatter(names,mHI)
plt.xticks(rotation='vertical')
plt.grid()
plt.yscale('log')
plt.show()


plt.scatter(names,jhi)
plt.scatter(names,jHI)
plt.xticks(rotation='vertical')
plt.grid()
plt.yscale('log')
plt.show()


plt.scatter(names,thtHI)
plt.scatter(names,ththi)
plt.xticks(rotation='vertical')
plt.grid()
plt.show()


plt.scatter(names,omgHI)
plt.scatter(names,omghi)
plt.xticks(rotation='vertical')
plt.grid()
plt.show()


plt.scatter(names,vmnHI)
plt.scatter(names,vmnhi)
plt.xticks(rotation='vertical')
plt.grid()
plt.show()
