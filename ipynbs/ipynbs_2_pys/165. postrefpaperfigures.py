#!/usr/bin/env python
# coding: utf-8


#make the versions of everyting with rhi = assuming vertical cylinder, with fixed n


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


plt.rcParams['font.size'] = 23
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


pixelize = lambda rhi, d: (np.arctan(rhi/d)*u.rad.to(u.deg)/0.125).astype('int') 

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
    

def get_mh2_rh2_rhi(i,thresh,norm=0.5):
    coords = coordses[i]
#     print(coords)
    h2cube = get_subcube(dht.data,coords[:4])
    mh2,rh2 = calc_mass(h2cube,ds[i],thresh=thresh)

    rhi = norm*(((mh2*u.solMass)/(np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
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

def spctrum_fit(i,rhi=0,which='hi'):
    coords = coordses[i]

#     sph2 = np.nanmean(h2cube,axis=(1,2)) #as things stand rn, we dont care abthe h2 fit
    if which=='hi':
        vmnh2, sigvh2 = np.array(spctrum_fit(i,which='h2')[1:]).transpose()
        vo,vi = ((vmnh2+3*np.array([-1,1])*sigvh2+100)/1.3).astype('int')
                
        loh,lih = lc_p[i]+np.array([-1,1])*rhi
        boh,bih = bc_p[i]+np.array([-1,1])*rhi
        
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
    
    
#---------------------------------------------------------------------------------------------------------------------------------


def chebyshev_grid(arr, orgn=0): #Returns an grid of same dim as arr where each cell contains its Chebyshev distance to 'orgn'.
    if orgn == 0:
        orgn = (np.array(arr.shape)/2).astype('int')
        
    coords = np.ogrid[tuple(slice(0, s) for s in arr.shape)]
    dx,dy = [np.abs(coord - p) for coord, p in zip(coords, orgn)]
    
    chebishv = np.nanmax(np.broadcast_arrays(dx,dy),axis=0)
    return chebishv, np.unique(chebishv)


def add_mass(mhi,hi_dm,chbshv,mh2,n=0): #put in a seed value of mhi as whats in the central pixel
    if (mhi >= mh2) or (n == chbshv[1][-1]):#reached the end of :
        # plt.imshow(hi_dm*(chbshv[0] <= n))
        # plt.show()
        return n, mhi/mh2, chbshv[1][-1]#, n
    
    
    n   += 1
    mhi = np.nansum(hi_dm*(chbshv[0] <= n))
    
    return add_mass(mhi,hi_dm,chbshv,mh2,n)


#------------------------------------------------------------------------------------------------------------------


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den

#------------------------------------------------------------------------------------------------------------------    


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

ibmh2 = np.array([2.2,9.0,6.9,1.2,7.7])*1e4*u.solMass

ib = pd.DataFrame(np.array((ibmh2,ibrh2,iboh2,ibjh2,ibrhi,ibohi,ibjhi)).transpose(),
                  columns='mh2,rh2,omgh2,jh2,rhi,omghi,jhi'.split(','))


in_ib = [12,13,15,18,20]


# vmnhi  = np.array([2.4,-0.5,0.2,17.2,1.2,5.9,11.8,6.4,5.0,0.9,2.7,27.6,8.8,9.5,2.0,8.3,8.8,6.7,2.5,-1.6,11.4,4.1])
# sigvhi = np.array([3.0,11.6,3.9,6.3,3.9,4.4,7.1,9.9,5.6,5.3,5.6,12.5,7.3,6.5,4.4,6.5,6.2,7.3,5.8,9.4,9.6,8.4])
vmnh2, sigvh2 = np.array([spctrum_fit(i,which='h2')[1:] for i in range(22)]).transpose()
vsh2 = ((vmnh2[:,np.newaxis]+3*np.array([-1,1])[np.newaxis,:]*sigvh2[:,np.newaxis] +100)/1.3).astype('int')

mh2,rh2,rhi = np.array([get_mh2_rh2_rhi(i,threshs[i]) for i in range(22)]).transpose()


vmnhi, sigvhi = np.array([spctrum_fit(i,pixelize(rhi[i],ds[i]),which='hi')[1:] for i in range(22)]).transpose()
vshi_px = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')

vshi_px[3] = ((vmnhi[3]+np.array([-1,1])*sigvhi[3] + 100)/1.3).astype('int')
vshi_px[11] = ((vmnhi[11]+np.array([-1,1])*sigvhi[11] + 100)/1.3).astype('int')
vshi_px[12] = ((vmnhi[12]+np.array([-1,1])*sigvhi[12] + 100)/1.3).astype('int')
vshi_px[20] = ((vmnhi[20]+np.array([-1,1])*sigvhi[20] + 100)/1.3).astype('int')


#-----------------------------------------------------------------------------------------------------------------------

rhi_px  = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')
mhi = np.array([calc_masses(i,rhi_px[i],vshi_px[i])[1] for i in range(22)])


rh2unc = rh2*duncs/ds
rhiunc = (2/3)*rhi*duncs/ds


#-----------------------------------------------------------------------------------------------------------------------

omgh2,omghi,omgh2unc,omghiunc,thth2,ththi = np.array([Omg_vec(i,rhi_px[i],vshi_px[i],unc=True) 
                                                        for i in range(22)]).transpose()


#-----------------------------------------------------------------------------------------------------------------------


jh2 = 0.4*(omgh2*u.m/u.s/u.pc*(rh2*u.pc)**2).cgs
jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs

jh2unc = np.sqrt((omgh2unc/omgh2)**2 + (2*rh2unc/rh2)**2)*jh2
jhiunc = np.sqrt((omghiunc/omghi)**2 + (2*rhiunc/rhi)**2)*jhi


#-----------------------------------------------------------------------------------------------------------------------

# tab2_newnew = pd.DataFrame(np.around(np.array((mh2/1e4,rh2,np.array(omgh2)/10,thth2,jh2/1e24,\
#                                   mhi/1e5,rhi,np.array(omghi)/10,ththi,jhi/1e24)),3).transpose(),
#                         columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


np.nanmedian(np.sqrt(3*np.pi/(32*c.G*(mh2*u.solMass)/(4/3*np.pi*(rh2*u.pc)**3))).to(u.Myr))


np.nanmedian(np.sqrt(3*np.pi/(32*c.G*(mhi*u.solMass)/(4/3*np.pi*(rhi*u.pc)**3))).to(u.Myr))


np.nanmedian(rhi*u.pc/(sigvhi*u.km/u.s)).to(u.Myr)


np.nanmedian(rh2*u.pc/(sigvh2*u.km/u.s)).to(u.Myr)


plt.scatter(omghi,omgh2)


np.nansum(np.ones(22)[np.abs(ththi)<90])


np.nanmedian((jhi/jh2)[in_ib])


np.nanmedian(ib['jhi']/ib['jh2'])


plt.scatter(mh2[in_ib],ib['mh2'])
plt.plot(mh2[in_ib],mh2[in_ib],)


plt.scatter(rhi[in_ib],ib['rhi'])
plt.plot(rhi[in_ib],rhi[in_ib],)


ib


plt.scatter(names[in_ib],omgh2[in_ib]/10,)
plt.scatter(names[in_ib],ib['omgh2'])
# plt.plot(omgh2,omgh2)
# plt.loglog()
plt.xticks(rotation='vertical')
plt.show()


mh2[in_ib]/ib['mh2']


plt.scatter(names[in_ib],mh2[in_ib],)
plt.scatter(names[in_ib],ib['mh2'])
# plt.plot(omgh2,omgh2)
# plt.loglog()
plt.xticks(rotation='vertical')
plt.show()


plt.scatter(names[in_ib],rhi[in_ib],)
plt.scatter(names[in_ib],ib['rhi'])
# plt.plot(omgh2,omgh2)

# plt.loglog()
plt.xticks(rotation='vertical')
plt.show()


np.nanmean((jhi/jh2)[in_ib])


ib['mh2']


ib['omgh2']


# # figure 1


hi_0th_Msunpc2 = (np.nansum(hi4pi.data,axis=0)*1.3*1.82e18/u.cm**2*c.m_p*1.36).to(u.Msun/u.pc**2).value


longitudes = [r'180$\degree$','','',r'150$\degree$','','',r'120$\degree$','','',\
                    r'90$\degree$','','',r'60$\degree$','','',r'30$\degree$','','',0,\
                    '','',r'-30$\degree$','','',r'-60$\degree$', '','',r'-90$\degree$',\
                    '','',r'-120$\degree$','','',r'-150$\degree$','','',r'-180$\degree$',]

latitudes = [r'-30$\degree$',r'-15$\degree$',0,r'15$\degree$',r'30$\degree$']


siglvs = [5,10,15,20,25,30,35,40,50,75,100]


plt.rcParams['font.size'] = 23


fig = plt.figure(dpi=200,figsize=(20,10))
ax = fig.add_subplot()#projection='mollweide')

imhi = ax.contourf(hi_0th_Msunpc2,extent=[360,-360,-50,50],
             levels=siglvs,extend='both',cmap='magma',norm='asinh',zorder=10)
# fig.colorbar(imhi,orientation='vertical',shrink=0.5,label='$\Sigma_\mathrm{HI}\ [M_\odot\ \mathrm{pc^{-2}}]$')
# plt.colorbar(shrink=0.5,label='$\Sigma_\mathrm{HI} \ [M_\odot \ \mathrm{pc^{-2}}]$')
ax.set_xlim(190,-190)
ax.set_aspect(1)

ax.set_frame_on(False)


# [ax.add_patch(Circle((lc[i],bc[i]),rhi_px[i]/8,fill=False,color='black',zorder=20)) for i in range(22)]
# [ax.add_patch(Circle((-360+lc[i],bc[i]),rhi_px[i]/8,fill=False,color='black',zorder=20)) for i in [1,21]]

[ax.add_patch(Rectangle((lc[i]-rhi_px[i]/8,bc[i]-rhi_px[i]/8),rhi_px[i]/4,rhi_px[i]/4,fill=False,color='white',zorder=20)) for i in range(22)]
[ax.add_patch(Rectangle((lc[i]-rhi_px[i]/8-360,bc[i]-rhi_px[i]/8),rhi_px[i]/4,rhi_px[i]/4,fill=False,color='white',zorder=20)) for i in [1,21]]


# [ax.add_patch(Circle((lc[i],bc[i]),pixelize(rh2[i],ds[i])/8,fill=False,color='blue',zorder=40)) for i in range(22)]
# [ax.add_patch(Circle((-360+lc[i],bc[i]),pixelize(rh2[i],ds[i])/8,fill=False,color='blue',zorder=40)) for i in [1,21]]


# [ax.add_patch(Circle((lc[i],bc[i]),pixelize(rhi_cnm[i],ds[i])/8,fill=False,color='green',zorder=40)) for i in range(22)]
# [ax.add_patch(Circle((-360+lc[i],bc[i]),pixelize(rhi_cnm[i],ds[i])/8,fill=False,color='green',zorder=40)) for i in [1,21]]


for i in range(22):
    coords = coordses[i]
    mapses = moment_maps(i,threshs[i],rhi_px[i],vshi_px[i])
    h2lb,hilb = mapses[0][0],mapses[1][0]

    ax.contourf(np.where(np.log10(h2lb)<threshs[i],np.nan,np.log10(h2lb)),extent=coords[:4],
                   levels=10,#np.linspace(threshs[i],threshs[i]*1.05,9),
                   cmap='YlGnBu_r',extend='max',
                   linewidths=0.75,zorder=30)
    
    
    # ax.contourf(hilb,extent=[lc[i]-rhi_px[i]/8,lc[i]+rhi_px[i]/8,bc[i]-rhi_px[i]/8,bc[i]+rhi_px[i]/8],
    #                levels=10,#np.linspace(threshs[i],threshs[i]*1.05,9),
    #                cmap='magma',extend='max',
    #                linewidths=0.75,zorder=25)

# imh2 = ax.contourf((h2lb*c.m_p/u.cm**2*2*1.36).to(u.solMass/u.pc**2).value, extent=coords[:4],
#                    levels=np.linspace(20.4,22,25),
#                    cmap='magma',
#                    linewidths=0.75,zorder=30)

# fig.colorbar(imh2,orientation='horizontal',ticks=np.linspace(20.4,22,25))
    

ax.set_xlim(181,-181)
ax.set_ylim(-33,33)
ax.set_xticks(np.arange(180,-181,-10),
              labels=longitudes)
ax.set_yticks(np.arange(-30,31,15),labels=latitudes)

ax.set_xlabel('Galactic Longitude')
ax.set_ylabel('Galactic Latitude')
ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)

# ax.grid(zorder=-30,color='black',alpha=0.5)
plt.savefig('fig1_cat',bbox_inches='tight')

plt.show()


# # figure 2


def draw_spectra(i,rhi,vovi):
    coords = coordses[i]
    
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vovi
    
    hicube = hi4pi.data[:, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords[:4]) #velocity axis length = 155
    
    
    hisp_long = np.nanmean(hicube,axis=(1,2))
    hi_spctr = np.nanmean(hicube[vo:vi],axis=(1,2))
    h2_spctr = np.nanmean(h2cube,axis=(1,2))

    Ahi,muhi,stdhi = spctrum_fit(i,rhi)
    mody = models.Gaussian1D(Ahi,muhi,stdhi)
    
    
    Ah2,muh2,stdh2 = spctrum_fit(i,which='h2')
    h2mody = models.Gaussian1D(Ah2,muh2,stdh2)
    
    fig,ax = plt.subplots(dpi=200,figsize=(6,6))
    
    ax.plot(np.arange(155)*1.3-100,hisp_long,linewidth=6,color=(0.8,0.,0.))
    ax.vlines(vovi*1.3-100,0,np.nanmax(hi_spctr)*1.5,linestyle='--', color='black',linewidth=0.8,zorder=-30)
    
    ax.set_xlim(muhi+np.array([-1,1])*stdhi*5)
    ax2 = ax.twinx()
    ax2.plot(np.arange(155)*1.3-100,h2_spctr,color=(0.,0.8,0.8),linewidth=2.5)
              
    ax2.set_ylim(0,1.1*np.nanmax(h2_spctr),)
    ax.set_ylim(0,1.1*np.nanmax(hi_spctr))
    
    ax.set_yticks(np.around(np.linspace(0,np.nanmax(hi_spctr),4)))
    ax2.set_yticks(np.around(np.linspace(0,np.nanmax(h2_spctr),4),2))
    ax.set_xticks(np.around(muhi+np.array([-4,-2,0,2,4])*stdhi,1))

    plt.title(names[i])
    plt.xlabel('$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
    ax.set_ylabel('$T_\mathrm{HI}\ [K]$',color=(0.8,0.,0.))
    ax2.set_ylabel('$T_\mathrm{CO}\ [K]$',color=(0.,0.8,0.8))
    ax.set_xlabel('$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
    
    # ax.annotate(str(np.around(muhi,1))+'$\pm$'+str(np.around(stdhi,1)),(muhi-4*stdhi,np.nanmax(hi_spctr)))
    # ax.annotate()
    
    
    plt.savefig(names[i][:3]+names[i][-1]+'_spec',bbox_inches='tight')
#     plt.show()


[(
    draw_spectra(i,rhi_px[i],vshi_px[i]),
plt.show()
)
for i in range(22)]


# # figure 3


h2_vlims = np.array(([2,4],[-10,6],[-1.5,0.5],[14,21],[1,5.5],[4.8,6.3],[10,15],[1,9],[5,13],[0,4],
                     [3,9],[23,32],[2,12],[9,13.5],[1.5,5],
                     [4,12],[6,11],[-3,11],[1,9.5],[-5,-1.5],[9,17],[4.,7.5]))
# hi_vlims = 
hi_Slims =  np.array(([2,6.5],[7,55],[4.5,7.5],[8,19],[2.5,10],[3,8.5],[15,24.5],[37,59],[5.5,15],[7,18],
                     [5,30],[10,35],[7,24],[9,21],[5,18],
                     [4,21],[7,24],[9,30],[4,16],[2,9],[16,37],[3,30]))

hi_vlims = np.array([[  -1. ,  2 ],
                       [-14. ,   2. ],
                           [ -2.7,   -0.7],
                               [ 12 ,  15.5 ],
                     [  -2. ,   2. ],
                       [  2.5,   5.5],
                           [ 6 ,  12 ],
                                [  3. ,   10. ],
                    [  3. ,  6. ],
                        [  -3.7 ,   0.3 ],
                           [  -1.5 ,   5 ],
                               [ 18. ,  25. ],
                    [  5.5 ,  9 ],
                        [  4.5 ,  9.5 ],
                           [  -1.5 ,   3. ],
                               [  3. ,  8.5 ],
                    [  4.5 ,  9.5 ],
                        [ 3. ,  8. ],
                           [  -1 ,  6 ],
                               [ -7.5 ,  -3.5 ],
                    [  9. ,  12.5 ],
                       [  -2. ,   6]])


def make_zxz(i,rhi,vovi,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i,threshs[i],rhi,vovi)
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi
    
    xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
    ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
    
    
    h2dvds, hidvds = Omg_vec(i,rhi,vovi,xy=True)/1e3
    
    h2ox, h2oy, h2vo = h2dvds
    hiox, hioy, hivo = hidvds
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(13,13))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=np.linspace(5,3*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),35),
                   cmap='YlGnBu_r',extend='max',)

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                    ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)
    
    
    axs[0,0].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),#np.linspace(5,20,50),
                cmap="binary_r",
                extend='both',zorder=-9,
                alpha=0.9
                   )
    axs[1,0].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='binary_r',extend='both',zorder=-10)
        
    for ax in axs.flat:#.transpose()[1]:
        ax.contour(h2Sig.value, extent=coords[:4],
                   levels=np.array([1,2,3])*0.25*np.nanmax(h2Sig.value),
                   # np.array([-1,0,1])*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),#np.linspace(2,50,5),
                   colors='black',
                   linestyles=[':','solid','--'],
                   linewidths=1.2,)
    
                
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        # ax.label_outer()
        ax.set_frame_on(False)
        ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)

    
    for im,ax,lbl in zip([im1,im2],
                           axs.flat[:2],
                          ['$\Sigma_\mathrm{H_2} \ [\mathrm{M_\odot\ pc^{-2}}]$',r'$\Sigma_\mathrm{HI} \ [\mathrm{M_\odot\ pc^{-2}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl,norm='asinh')
        
        
    for im,ax,lbl in zip([im3,im4],
                           axs.flat[2:],
                          ['',r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl,norm='asinh')
        
        
    dx = (loh-lih)*0.1
    # print(dx)
    
    axs[1,0].add_patch(Circle((lih+dx,boh+dx),dx*0.8,color='white'))
    axs[1,0].quiver(lih+dx,boh+dx,-h2ox,h2oy,color=(0.,0.8,0.8),headwidth=3,zorder=30,pivot='tail',alpha=0.8)
    axs[1,0].quiver(lih+dx,boh+dx,-hiox,hioy,color=(0.8,0.,0.),headwidth=3,zorder=30,pivot='tail',alpha=0.8)
    
    axs[0,1].annotate(str(int(np.nanmax(h2Sig.value))),(loh-0.5*dx,boh+0.5*dx),color=(0.,0.8,0.8),fontsize='x-small')
    axs[0,1].contour(h2Sig.value, extent=coords[:4],
                   levels=np.array([1,2,3])*0.25*np.nanmax(h2Sig.value),
                   # np.array([-1,0,1])*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),#np.linspace(2,50,5),
                   colors=(0.,0.8,0.8),
                   linestyles=[':','solid','--'],
                   linewidths=1.2,)
        
        
    fig.suptitle(names[i])
        
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


[make_zxz(i,rhi_px[i],vshi_px[i],save=True) for i in range(22)]


[make_zxz(i,rhi_px[i],vshi_px[i],save=True) for i in range(22)]


# 

# # figure 4


#delta theta


delthet = np.array(ththi)-np.array(thth2)
for i in range(22):
    if delthet[i] > 180:
        delthet[i] -= 360
    if delthet[i] < -180:
        delthet[i] += 360


plt.rcParams['font.size'] = 24


plt.figure(dpi=200,figsize=(6,6))
plt.hist(delthet,
        bins=np.arange(-180,181,90/4),edgecolor='white',color='black')
plt.xticks(np.arange(-180,181,90))
plt.xlabel(r'$\Delta \theta\ [\degree]$')
plt.ylabel('N (clouds)')
plt.yticks([1,2,3,4,5,6])
plt.ylim(0,5.5)
plt.xlim(-185,185)
plt.vlines([-45,45],0,7,linestyle='--',color='grey',)

plt.annotate('counter-rotating',(-165,5.2),fontsize='xx-small')
plt.annotate('counter-rotating',(65,5.2),fontsize='xx-small')
plt.annotate('co-rotating',(-32,5.2),fontsize='xx-small')
plt.gca().tick_params(which='both',top=True,bottom=True, left=True, right=True)

# 
plt.savefig('fig4_delthet',bbox_inches='tight')


# 


# # figure 5


plt.rcParams['font.size'] = 23


def plot_rvsv(i,rhi,vovi, show=False):  
    h2dvds, hidvds = Omg_vec(i,rhi,vovi,xy=True)/1e3
    
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

    
    hbperp = (hidvds[0]*hxs_pc + hidvds[1]*hys_pc)/omghi
    cbperp = (h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc)/omgh2

    
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
                        np.nanmax(hdf['rper'])+hR/15,
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper'])+hR/15,
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
                        np.nanmax(cdf['rper'])+cR/15,
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper'])+cR/15,
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
        axs[0].set_xticks(np.around([-cR,-cR/2,0,cR/2,cR]))
        axs[0].scatter(0,h2dvds[2],marker='+',color='black',)
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
        axs[1].set_xticks(np.around([-hR,-hR/2,0,hR/2,hR]))
        axs[1].scatter(0,hidvds[2],marker='+',color='black',)
    
        
        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
            # ax.grid()
            ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)

        axs[0].set_ylabel(r'$\langle v_\mathrm{LSR}\ \rangle\ [\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names[i])
        
        axs[0].set_title(np.around(np.abs(h2p),2),y=0.87,x=0.13)
        axs[1].set_title(np.around(hip,2),y=0.87,x=0.13)

#         axs[1].annotate(str(np.around(hip,2)),(meanhdf['rper'].iloc[0],
#                                                np.nanmax(meanhdf['v'])+np.nanmean(sighdf['v'])))
        
    
        fig.savefig('v_vs_r_'+names[i],bbox_inches='tight')
        plt.show()

    return h2p, hip
#     return pearson(meancdf['rper'],  meancdf['v']),                   pearson(meanhdf['rper'],  meanhdf['v']),\
#            pearson(sigcdf['v'],      np.abs(meancdf['v']-h2_vos[i])), pearson(sighdf['v'],      np.abs(meanhdf['v']-hi_vos[i])),\
#            pearson(meancdf['w'],     np.abs(meancdf['v']-h2_vos[i])), pearson(meanhdf['w'],     np.abs(meanhdf['v']-hi_vos[i]))


[plot_rvsv(i,pixelize(rhi[i],ds[i]),vshi_px[i],show=True) for i in range(22)]


np.array([omgh2,omghi]).transpose()/1e3


np.array([plot_rvsv(i,pixelize(rhi[i],ds[i]),vshi_px[i],show=True) for i in range(22)]).transpose()


np.nanmedian(hip)


np.nanmedian(h2p)


# # figure 6


h2p,hip = np.array([plot_rvsv(i,pixelize(rhi[i],ds[i]),vshi_px[i]) for i in range(22)]).transpose()


plt.rcParams['font.size'] = 17
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


fig = plt.figure(dpi=200)#,figsize=(6,6))
ax = fig.add_subplot()

[ax.scatter(np.abs(pval),np.abs(tht),marker='.', alpha=0.8,s=500,c=c,edgecolor='white',label=h) 
        for tht,pval,c,h
        in zip([thth2,thth2],[h2p,hip],[(0.,0.8,0.8),(0.8,0.,0.)],['$\mathrm{H_2}$','$\mathrm{HI}$'])]

ax.vlines([0.5,0.75],0,190,zorder=-30,color='black',linestyle='--')
ax.hlines([60,120],-0.1,1.1,zorder=-30,color='black',linestyle='--')


# plt.xticks(np.arange(0,22,3),labels=names[::3],rotation='vertical')
# plt.grid()

ax.set_ylabel(r'$\theta\ [\degree]$')
ax.set_xlabel(r'$\rho\ (r_\perp,v)$')
ax.set_ylim(0,185)
ax.set_xlim(-0.02/3,1.02)
ax.set_xticks(np.arange(0.,1.1,0.25),[0,0.25,0.5,0.75,1])
ax.set_yticks(np.arange(0,181,30))
ax.legend(loc='upper left',fontsize='x-small',)
# ax.tick_params(top=True,bottom=True, left=True, right=True)


xhist = ax.inset_axes([0, 1.05, 1, 0.25],)

# xhist.hist(ps[0],bins=np.arange(-1.,1.1,0.125),color=magnt,alpha=0.5,edgecolor='black',linewidth=0.5)#,histtype='barstacked')
[xhist.hist(np.abs(pval),bins=np.arange(-0.5,1.1,0.5/10),color=c,
            edgecolor='white',linewidth=1.5,alpha=0.5)
for pval,c
in zip([h2p,hip],[(0.,0.8,0.8),(0.8,0.,0.)],)]


xhist.set_xticks(np.arange(0.,1.1,0.25),[])
xhist.set_xlim(-0.02,1.02)
xhist.set_yticks([0,4,8])
xhist.set_ylabel('N')
# xhist.tick_params(top=True,left=True, right=True)


# xhist.hlines(10,0,1)
plt.savefig('fig6_theta_vs_p',bbox_inches='tight')

plt.show()


# # figure 7


plt.rcParams['font.size'] = 23


phihi = ththi-90
phih2 = thth2-90

for i in range(22):
    if phihi[i] < -180:
        phihi[i] += 360
        
    if phih2[i] < -180:
        phih2[i] += 360


plt.rcParams['font.size'] = 22


fig,axs = plt.subplots(2,2,dpi=200,figsize=(8,8),sharey=True)

jbns=10**np.arange(22,27,1/5)
thbns=np.arange(-180,181,90/4)

axs[0,0].hist(jh2.value,bins=jbns,color=(0.,0.8,0.8),zorder=15,edgecolor='white',label='$\mathrm{H_2}$')
axs[1,0].hist(jhi.value,bins=jbns,color=(0.8,0.,0.),zorder=15,edgecolor='white',label='$\mathrm{HI}$')

for ax in axs[:,0]:
    ax.set_xscale('log')
    ax.set_xlim(1e23,1e26)
    ax.set_ylabel('N (clouds)')
    ax.set_ylim(0,5.1)
    ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)
    ax.label_outer()

    
for ax in axs[:,1]:
    ax.set_xticks(np.arange(-180,181,90))
    ax.set_xlim(-185,185)
    # ax.vlines([-90,90],0,6,linestyle='--',color='black',zorder=-20)
    ax.set_ylim(0,5.1)
    ax.set_yticks([1,2,3,4,5])
    ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)
    ax.label_outer()
    

axs[0,1].hist(phih2,bins=thbns,edgecolor='white',color=(0.,0.8,0.8))
axs[1,1].hist(phihi,bins=thbns,edgecolor='white',color=(0.8,0.,0.))


# axs[0,0].set_xlabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')
axs[1,0].set_xlabel('$j\ [\mathrm{cm^2\ s^{-1}}$]')
    
# axs[0,1].set_xlabel(r'$\theta_\mathrm{H_2}\ [\degree$]')
axs[1,1].set_xlabel(r'$\phi\ [\degree$]')
    
    
axs[0,0].legend(loc='upper left',fontsize='x-small')
axs[1,0].legend(loc='upper left',fontsize='x-small')


fig.tight_layout()

# plt.ylabel('N (clouds)')
# plt.xlim(8e22,3e26)
# plt.ylim(0,5.5)
# # plt.yticks(np.arange(5))

plt.savefig('fig7_js',bbox_inches='tight')

# plt.show()


# # figure 8


curve_fit(line,np.log10(jh2.value),np.log10(jhi.value),sigma=np.log10(jhiunc.value))


np.sqrt(0.01508102)


np.log10(3.7e9)


plt.rcParams['font.size']=20
plt.figure(dpi=200,figsize=(5,5))

xs = np.linspace(5e22,5e26,5)

plt.plot(xs,xs,color='grey',zorder=-15,alpha=0.5)
# plt.fill_between(xs,3.7e9*xs**(0.62+0.12),3.7e9*xs**(0.62-0.12),alpha=0.5)
plt.plot(xs,3.7e9*xs**0.62,color='black', linestyle='--',zorder=-9,
         label=r'$j_\mathrm{HI}\ \propto\ j_\mathrm{H_2}^{0.62\pm0.12}$')

plt.errorbar(jh2.value,jhi.value,xerr=jh2unc.value,yerr=jhiunc.value,
             linestyle='',zorder=-5,color='black',capsize=2,alpha=0.3)

plt.legend(loc='upper left')

plt.scatter(jh2.value,jhi.value,
            color='black',edgecolor='white',s=100)
plt.gca().tick_params(which='both',top=True,bottom=True, left=True, right=True)


plt.xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel(r'$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}]$')

# plt.yticks([1,10],labels=['1','10'])
plt.xlim(9e22,4e25)
plt.ylim(2e23,8e25)
# plt.hlines([1,10],5e22,1e27,linestyle='-.',color='black')
# plt.ylim(5e-1,5e1)

plt.loglog()
plt.savefig('fig8_jj',bbox_inches='tight')


np.array((names,np.around(rh2unc,1),np.around(rhiunc,1))).transpose()


# # figure 9


curve_fit(line,np.log10(rh2),np.log10(jh2.value))


np.sqrt(0.05415663)


curve_fit(line,np.log10(rhi),np.log10(jhi.value))


np.sqrt(0.07194111)


10**22.45


np.log10(7e21)


fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5),sharey=True,)

xs = np.linspace(1,200,10)


axs[0].scatter(rh2,jh2,color=(0.,0.8,0.8),s=100,zorder=30,edgecolor='white')
axs[0].scatter(rh2[h2p<0.5],jh2[h2p<0.5],color='white',s=100,zorder=30,edgecolor=(0.,0.8,0.8))


axs[1].scatter(rhi,jhi,color=(0.8,0.0,0.),s=100,zorder=30,edgecolor='white')
axs[1].scatter(rhi[hip<0.75],jhi[hip<0.75],color='white',s=100,zorder=30,edgecolor=(0.8,0.,0.))

    
axs[0].plot(xs, 3e22*xs**1.67,color='black',linestyle='--',label='$j_\mathrm{H_2} \propto R^{1.67 \pm 0.22}$')
axs[1].plot(xs, 1e22*xs**1.71,color='black',linestyle='--',label='$j_\mathrm{HI} \propto R^{1.71 \pm 0.27}$')

axs[0].fill_between(xs,3e22*xs**(1.67+0.22),3e22*xs**(1.67-0.22),alpha=0.2,color='grey',linewidth=0.)
axs[1].fill_between(xs,1e22*xs**(1.71+0.27),7e21*xs**(1.71-0.27),alpha=0.2,color='grey',linewidth=0.)


axs[0].errorbar(rh2,jh2.value,xerr=rh2unc,yerr=jh2unc.value,linestyle='',zorder=-5,color='black',capsize=2,alpha=0.3)
axs[1].errorbar(rhi,jhi.value,xerr=rhiunc,yerr=jhiunc.value,linestyle='',zorder=-5,color='black',capsize=2,alpha=0.3)


axs[0].set_xlim(3,100)
axs[0].set_xlabel(r'$R_\mathrm{H_2}$ [pc]')
axs[0].set_ylabel('$j\ [\mathrm{cm^2\ s^{-1}}$]')

axs[1].set_xlim(3,100)
axs[1].set_xlabel(r'$R_\mathrm{HI}$ [pc]')


axs[0].legend(loc='upper left')
axs[1].legend(loc='upper left')

for ax in axs:
    ax.loglog()
    
    ax.set_yticks([1e23,1e24,1e25,1e26])
    ax.set_ylim(9e22,1e26)
    ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)

    
fig.tight_layout()


plt.savefig('fig9_jRs',bbox_inches='tight')


# # figure 10


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


#--------------------------------------------------------------------------------------------------------------------


p99_tab1 = pd.read_csv('../../data/J_A+AS_134_241_table1.dat.gz.txt',skiprows=9,delimiter='|')
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
# p99_dat = p99_dat[~(p99_dat['type']=='D/R ')] #no disks/rings
# p99_dat = p99_dat[~np.array([('F' in nm) for nm in p99_dat['type']])] #no filaments

p99_dat['j<cm2/s>'] = (0.4*(p99_dat['R<pc>'].to_numpy()*u.pc)**2 *p99_dat['Omega<km/s/pc>'].to_numpy()*u.km/u.s/u.pc).cgs.value.astype('float')


chenx_rs = np.array([1703,1232,3530,1440,1272,1763,1329,2762])*u.AU.to(u.pc)
chenx_js = (np.array([0.12,0.12,0.45,0.13,0.24,0.12,0.42,0.45])*1e-3*u.km/u.s*u.pc).cgs.value


chenh_rs = np.array([0.038,0.033,0.025,0.032,0.026,0.026,0.034,0.054,0.042,0.020,0.022,0.024,0.19,0.028,0.045,0.044,0.050,0.048])
chenh_js = (np.array([1.3,1.3,0.36,1.3,0.8,0.52,0.59,2.5,1.,0.8,1.5,1.6,6.,0.37,0.25,0.34,0.46,0.54])*1e-3 * u.km/u.s*u.pc).cgs.value


goodman = pd.read_csv('../csvs/goodman_1993_data.csv')

tatmtsu_rs = np.array([0.049,0.088,0.083,0.084,0.059,0.099,0.091,0.049,0.098,0.095,0.052,0.104,0.094,0.083,0.111,0.105,0.074,0.057,\
                       0.072,0.043,0.062,0.074,0.115,0.083,0.084,0.062,0.074,0.081,0.062,0.052,0.083,0.157,0.057,0.072])
tatmtsu_js = np.array([7.4e20,1.7e21,7.6e20,7.2e20,0,6e20,1.2e21,1.6e21,9.3e21,3.2e21,3e20,2.3e21,1.8e21,2.2e21,2.1e21,1.9e21,2.2e21,\
                       0,0,8.4e20,1.5e21,0,1.9e21,4.4e21,0,7.7e20,0,1.1e21,1.2e21,1.1e21,3e21,3.5e21,8.9e20,0
])


tatmtsu_rs = tatmtsu_rs[~np.isinf(np.log10(tatmtsu_js))]
tatmtsu_js = tatmtsu_js[~np.isinf(np.log10(tatmtsu_js))]


caselli_dat = pd.read_csv('caselli_tab3.csv').merge(pd.read_csv('caselli_tab5.csv'),on='Core')
caselli_js = 0.4*((caselli_dat['r (pc)'].to_numpy()*u.pc)**2 *caselli_dat['G (km s^-1 pc^-1)'].to_numpy()*u.km/u.s/u.pc).cgs.value


ga1985 = pd.read_csv('goldsmitharquila_1985.csv')
ga1985_rs = ga1985['Size (pc)'].to_numpy()
ga1985_js = 0.4*((ga1985_rs*u.pc)**2 * (ga1985['omega (km s^-1 pc^-1)'].to_numpy()*u.km/u.s/u.pc)).cgs


pirogov = pd.read_csv('pirogov_tab7.csv').merge(pd.read_csv('pirogov_tab5.csv'),on='Source')
pirogov_j = (pirogov['Gamma_km_s-1_pc-1'].to_numpy()*u.km/u.pc/u.s * (pirogov['d_pc'].to_numpy()*u.pc)**2).cgs.value


xu = pd.read_csv('xu_tab1.csv')
xu_rs = np.nanmean(xu[['Major_arcsec','Minor_arcsec']].to_numpy(),axis=1)*u.arcsec.to(u.rad)*432
xu_js = (0.4*(xu_rs*u.pc)**2 * xu['G_km_s-1_pc-1'].to_numpy()*u.km/u.s/u.pc).cgs.value

liu =  pd.read_csv('../../data/liu_tab1.csv')
liu_j = (0.4*(liu['Rc'].to_numpy()*u.pc)**2 *liu['omega_obs'].to_numpy()*u.km/u.s/u.pc).cgs.value


tobin_dat = pd.read_csv('tobin_tab1.csv').merge(pd.read_csv('tobin_tab6.csv'),on='Source').merge(pd.read_csv('tobin_tab9.csv'),on='Source')
tobin_rs = np.nanmean(tobin_dat[['Radius_arcsec','r2']],axis=1)*u.arcsec.to(u.rad) * tobin_dat['Distance (pc)'].to_numpy()
tobin_js = (0.4*(tobin_rs*u.pc)**2 * (tobin_dat['GradFit2D_km_s_pc'].to_numpy()*u.km/u.s/u.pc)).cgs.value


roslwsky07 = pd.read_csv('rosolowsky07.csv')
roslwsky07_js = (0.4*(roslwsky07['R_e (pc)'].to_numpy()*u.pc)**2 * roslwsky07['|grad_v| (km/s/pc)'].to_numpy()*u.km/u.s/u.pc).cgs


yen_dat = pd.read_csv('yen_tab1.csv').merge(pd.read_csv('yen_tab5.csv'),on='Source').merge(pd.read_csv('yen_tab8.csv'),on='Source')
yen_rs = np.nanmean(np.array([rs.split('x') for rs in yen_dat['Deconvolved Size (arcsec)']]).astype('float'),axis=1)*u.arcsec.to(u.rad)*yen_dat['Distance (pc)']
yen_js = np.nanmean(yen_dat[['jmin (km s^-1 pc)','jmax']],axis=1)*(u.km/u.s*u.pc).to(u.cm**2/u.s)

liu21 = pd.read_csv('liu2021_Table1.csv')
liu21_js = (0.4*(liu21['      R'][1:].to_numpy().astype('float')*u.pc)**2 *liu21['    Omega'][1:].to_numpy().astype('float')*u.km/u.s/u.pc).cgs

puvanova_dat = pd.read_csv('puvanova_dat.csv')


all_rs,all_js = np.concatenate(([rh2,jh2.value],
                         [rhi,jhi.value],
                         [p99_dat['R<pc>'],p99_dat['j<cm2/s>']],
                         [chenh_rs,chenh_js],
                         [chenx_rs,chenx_js],
                         [tatmtsu_rs,tatmtsu_js],
                         [ibb['rh2'],ibb['jh2']],
                         # [goodman['R_eff<pc>'],goodman['j<cm2/s>']],
                         [pirogov['d_pc'],pirogov_j],
                         [caselli_dat['r (pc)'],caselli_js],
                         [xu_rs,xu_js],
                         [liu['Rc'],liu_j],
                         [tobin_rs,tobin_js],
                         [puvanova_dat['R (pc)'],1e20*puvanova_dat['J/M (10²⁰ cm² s⁻¹)']],
                         [roslwsky07['R_e (pc)'],roslwsky07_js],
                         [liu21['      R'][1:].to_numpy().astype('float'),liu21_js],
                         [yen_rs,yen_js]
                        ),axis=1).astype('float')


cond = (all_rs>0.) & (all_js>0.)


curve_fit(line,np.log10(all_rs[cond]),np.log10(all_js[cond]))


np.sqrt(3.07078834e-04)


np.log10(8.5e22)


plt.rcParams['font.size'] = 24


plt.figure(dpi=200,figsize=(8,8))

ss = 45
lw = 0.75

xs = 10**np.linspace(-3,3,10)
plt.plot(xs, 8.5e22*xs**1.50,color='black',linewidth=0.75,zorder=-6,linestyle='--',label='$j \propto R^{1.50\pm0.02}$')
plt.fill_between(xs, 8.5e22*xs**(1.50+0.02),8.5e22*xs**(1.50-0.02),color='black',linewidth=0.,alpha=0.2,zorder=-25)

# [plt.plot(xs,(0.4*(omg*u.km/u.s/u.pc)*(xs*u.pc)**2).cgs) for omg in [1e-2,1e-1,1,1e1,1e2,1e3]]
# [plt.plot(xs, omg*(u.km/u.s/u.pc*(xs*u.pc**2)).cgs,color='grey',zorder=-10,alpha=0.3,linewidth=0.75) for omg in [1e-3,1e-2,1e-1,1,1e1,1e2]]


plt.scatter(rh2,jh2,                              label=r'$\mathrm{H_2}$',color=(0.,0.8,0.8),s=1.5*ss,linewidth=lw,edgecolor='white',zorder=15)
plt.scatter(rhi,jhi,                              label='HI',             color=(0.8,0.,0.), s=1.5*ss,linewidth=lw,edgecolor='white',zorder=15)

plt.scatter(liu21['      R'][1:].to_numpy().astype('float'),
                                            liu21_js,label='L+21', color='white', s=ss, marker='o',edgecolor='black',linewidth=lw,zorder=-10)
plt.scatter(liu['Rc'],liu_j,                        label='L+25', color='white', s=ss, marker='s',edgecolor='black',linewidth=lw,zorder=-10)
plt.scatter(xu_rs,xu_js,                            label='X+20', color='black', s=ss, marker='X',edgecolor='white',linewidth=lw,)
plt.scatter(chenh_rs,chenh_js,                      label="C+19", color='black',s=2*ss,marker='*',edgecolor='white',linewidth=lw,)
plt.scatter(puvanova_dat['R (pc)'],
            1e20*puvanova_dat['J/M (10²⁰ cm² s⁻¹)'],label='P+18', color='black', s=ss, marker='p',edgecolor='white',linewidth=lw,zorder=-10)
plt.scatter(tatmtsu_rs,tatmtsu_js,                  label='T+16', color='black', s=ss*0.75, marker='D',edgecolor='white',linewidth=lw,)
plt.scatter(yen_rs,yen_js,                          label='Y+15', color='black', s=ss, marker='1',edgecolor='white',linewidth=lw*2)
plt.scatter(ibb['rh2'],ibb['jh2'],                  label='I+11', color='white', s=ss, marker='s',edgecolor='black',linewidth=lw,zorder=10)
plt.scatter(tobin_rs,tobin_js,                      label='T+11', color='black', s=ss, marker='s',edgecolor='white',linewidth=lw)
plt.scatter(chenx_rs,chenx_js,                      label='C+07', color='black', s=ss, marker='^',edgecolor='white',linewidth=lw,)
plt.scatter(roslwsky07['R_e (pc)'],roslwsky07_js,    label='R07', color='white', s=ss, marker='^',edgecolor='black',linewidth=lw)
plt.scatter(pirogov['d_pc'],pirogov_j,              label='P+03', color='black', s=ss, marker='2',edgecolor='white',linewidth=lw*2,)
plt.scatter(caselli_dat['r (pc)'],caselli_js,       label='C+02', color='black', s=ss, marker='<',edgecolor='white',linewidth=lw,)
plt.scatter(p99_dat['R<pc>'],p99_dat['j<cm2/s>'],   label='P99',  color='black', s=ss, marker='o',edgecolor='white',linewidth=lw,zorder=10)
# plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],label='G+93',color='black', s=ss, marker='>',edgecolor='white',linewidth=lw,)
# plt.scatter(ga1985_rs,ga1985_js,                     label='G&A85', color='black',s=ss, marker='1',edgecolor='white',linewidth=0.5,)


# plt.scatter(all_rs,all_js)
# plt.scatter(yen_rs_d,yen_js)


plt.xlim(8e-4,1.5e2)
plt.ylim(7e18,5e26)
plt.yticks([1e20,1e22,1e24,1e26])
plt.loglog()


plt.gca().tick_params(which='both',top=True,bottom=True, left=True, right=True)

plt.legend(fontsize=15.5,loc='upper left',ncol=2)
# plt.legend(fontsize=14.5,loc='lower right',ncol=2)

# plt.annotate(,(5e-2,7e25),color='black')

plt.xlabel('$R$ [pc]')
plt.ylabel('$j\ [\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('fig10_jR_all',bbox_inches='tight')


np.nanmean(rhi/rh2)


# # figure 11


omgxyses = np.array([Omg_vec(i,rhi_px[i],vshi_px[i],xy=True)[:,:2]/1e3 for i in range(22)])


jinitxy = np.nansum((0.4*omgxyses*u.km/u.s/u.pc*\
                     np.stack([mh2,mhi]).transpose()[:,:,np.newaxis]*u.solMass*\
                     (np.stack([rh2,rhi]).transpose()[:,:,np.newaxis]*u.pc)**2)/\
                     ((mh2+mhi)[:,np.newaxis,np.newaxis]*u.solMass),axis=1).cgs.value


jinit_sc = (jh2*mh2 + jhi*mhi)/(mhi+mh2)


jinit_mag = np.linalg.norm(jinitxy,axis=1)


plt.figure(dpi=200)
plt.scatter(jh2,jinit_mag,s=60,alpha=0.5,label='vector sum')
plt.scatter(jh2,jinit_sc,s=60,alpha=0.5,label='scalar sum')

plt.loglog()
plt.xlabel('jh2')
plt.ylabel('jinit')
plt.legend()


def get_tdis(jhi,rhi,mhi,omghi,sigvhi):
    egrav = c.G*mhi*u.solMass/(rhi*u.pc)
    erot = (omghi*u.m/u.s/u.pc*rhi*u.pc)**2
    ekin = (sigvhi*u.km/u.s)**2
        
    return (jhi/(-egrav+erot+ekin).cgs).to(u.Myr)


tdis = get_tdis(jhi,rhi,mhi,omghi,sigvhi)


np.nanmean(tdis)


omgxyses


i


np.nanmedian(jinit_mag)


np.nanmin(jinit_mag),np.nanmax(jinit_mag)


np.nanmean(jinit_mag/jh2),np.nanmedian(jinit_mag/jh2)


bns=10**np.arange(22,26,0.2)

# plt.hist(jh2.value,alpha=0.5,bins=bns)
plt.hist(jhi.value,alpha=0.5,bins=bns)
plt.hist(jinit_mag,alpha=0.5,bins=bns)

plt.xscale('log')


np.median(jinit_mag/jhi)


np.nanmedian(jhi/jinit_mag)


plt.hist(tdis.value,bins=np.arange(0,3,0.1))


np.nanmean(mhi/mh2),np.nanmedian(mhi/mh2)


np.nanmean(omgh2),np.nanmean(omghi)


np.nanmedian(omgh2),np.nanmedian(omghi)


np.nanmean(jhi/jh2)


np.nansum(np.ones(22)[np.abs(ththi)<90])


thth2


16/22


np.nanmedian(tdis)


plt.hist(ththi,bins=np.arange(-180,181,15),edgecolor='white')
plt.xticks(np.arange(180,-180,-30))


15/22


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(lc,bc,ththi)
ax.view_init(azim=90,elev=40)


np.nanmean(rhi*u.pc/(sigvhi*u.km/u.s)).to(u.Myr)


np.nanmean(rh2*u.pc/(sigvh2*u.km/u.s)).to(u.Myr)


np.nanmean(rh2/sigvh2)


moment_maps


15/22


(jhi/(omghi*u.m/u.s/u.pc * rhi*u.pc)**2).to(u.Myr)[[12]]


plt.scatter(names,(jhi/(omghi*u.m/u.s/u.pc * rhi*u.pc)**2).to(u.Myr),)
plt.xticks(rotation='vertical')

plt.hlines([5,20],0,22)
plt.show()


200/13


np.nanmedian((jhi/(omghi*u.m/u.s/u.pc * rhi*u.pc)**2).to(u.Myr))


duncs


# # figure 12


plt.figure(dpi=200,figsize=(6,5))
plt.scatter(ththi,bc,c=lc,cmap='vanimo',edgecolor='white',linewidth=.5,s=100,vmin=-180,vmax=180)
plt.vlines([-90,90],-40,40,zorder=-30,color='black',linestyle='--')
plt.hlines([15,-15],-200,200,zorder=-30,color='black',linestyle='--')

plt.colorbar(label='$l\ [\degree]$',ticks=np.linspace(-180,180,5))#[-180,-90,0,90,180])
# plt.scatter(np.abs(ththi),np.abs(bc))

# [plt.annotate(names[i][:3]+names[i][-1],
#               (np.abs(ththi[i]),np.abs(bc[i])),fontsize='x-small') 
#  for i in range(22)]

# [plt.annotate(names[i],
#               (lc[i],ththi[i]),fontsize='xx-small') 
#  for i in range(22)]

#               
plt.ylabel(r'$b\ [\degree]$')
plt.xlabel(r'$\theta_\mathrm{HI}\ [\degree]$')

plt.ylim(-30,30)
plt.xticks(np.arange(-180,181,90))
plt.xlim(-185,185)
plt.yticks(np.arange(-30,31,15))
plt.gca().tick_params(which='both',top=True,bottom=True, left=True, right=True)


plt.savefig('fig12_thetas',bbox_inches='tight')
# 
plt.show()


# # tables


# # notes


# - fig10: add in data points with smaller radii (from acvs, the only one left is:
#     also look on my own) + rosolowsky/braine data points /others at larger end
# - - annotate this more with (contours for diff evolutionary stages), clean up lines of constant omega?)
# - might also be in teresting to keep the disks to show that there is a clear breaak at smaller radii, so the message of this plot becomes more complex ie we are probingt the top and bottom ends of the supersonic turbulence cascade in the ism
# (or at least even if we do the low end we can say this is when we go from omega goes as R to () to keplerian
# 
# 
# 
# --hmm also wondering whether its more straightforward to plot omega vs r, ie dont need to assume rotation
# 
# #remove individual citation labels (except phillips, 7 in arroyo just say is there, can cite them individually in the body where i say they combine the data from ...ppl) in the legend and replace with numbers then in the caption say who each is
# 
# 
# in the text, say how many different tracers, parts of the galaxy are jointly represented by the observations (!!find southern hemisphere stuff esp! make map of where these are? i think i have everyones lb) -- and also method of measurement (i think tatematsu is the only one that doesnt fit a plane?exclude?)
# 
# potential:
# 
# https://ui.adsabs.harvard.edu/abs/2019MNRAS.490..527C/abstract
# https://ui.adsabs.harvard.edu/abs/2019ApJ...872...30S/abstract
# https://ui.adsabs.harvard.edu/abs/2019ApJ...870...81T/abstract
# https://ui.adsabs.harvard.edu/abs/2018ApJ...853..169D/abstract
# https://ui.adsabs.harvard.edu/abs/2017ApJ...849...90S/abstract
# https://ui.adsabs.harvard.edu/abs/2017ApJ...840..119I/abstract
# https://ui.adsabs.harvard.edu/abs/2017ApJ...840..119I/abstract
# https://www.aanda.org/articles/aa/pdf/2015/03/aa22955-13.pdf
# 
# 
# 
# 
# 
# for moere, resume at page 4 of: https://ui.adsabs.harvard.edu/search/fq=%7B!type%3Daqp%20v%3D%24fq_database%7D&fq_database=(database%3Aastronomy%20OR%20database%3Aphysics)&q=molecular%20cloud%20core%20velocity%20gradient&sort=date%20desc%2C%20bibcode%20desc&p_=3 
# + theres a window with tabs open of relevant pprs to add to ^ list
# 
# extragalactic;
# braine+roslowsky
# see if eric / friends have stuff (esp if hi, which would help with filling in the top end)
# https://ui.adsabs.harvard.edu/abs/2025MNRAS.541.3081L/abstract


np.log10(60*u.AU.to(u.pc))
