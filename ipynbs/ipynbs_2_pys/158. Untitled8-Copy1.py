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
    vo,vi = vsnew[i]
    h2cube,hicube = make_h2_hi_cubes(i)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18
    

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            (np.log10(hilb)>20)
#     print(hi_msk)
    
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
    h20,hi0 = [moment_maps(i,threshs[i])[j][0]/u.cm**2 for j in [0,1]]
    
    
    dA = (ds[i]*u.pc*0.125*u.deg.to(u.rad))**2
    h2m = (np.nansum(h20*dA)*c.m_p*1.36*2).to(u.solMass)
    him = (np.nansum(hi0*dA)*c.m_p*1.36).to(u.solMass)
    
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


hi4pi, dht = fits.open('../../data/hi4piDHT_720deg.fits')
table1 = pd.read_csv('../../data/molecular_clouds.csv')

coordses = table1[['l0 (°)','l1 (°)',
                   'b0 (°)','b1 (°)',
                  'v0 (km/s)','v1 (km/s)']].to_numpy().astype('float')               

ds    = table1['d (pc)'].to_numpy().astype('float')
names = table1['Name']

lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)
lc_p = ((lc-360)/-0.125).astype('int')
bc_p = ((bc+50)/0.125).astype('int')

muh2_3sig     = [mu+np.array([-1,1])*2*std  for i in range(22) for mu,std in [spctrum_fit(i,which='h2')[1:]]]
muh2_3sig[11] = [mu+np.array([-1,1])*std    for i in [11] for mu,std in [spctrum_fit(i,which='h2')[1:]]][0]
muh2_3sig[19] = [mu+np.array([-1,1])*3*std  for i in [19] for mu,std in [spctrum_fit(i,which='h2')[1:]]][0]
vs_px     = ((np.array(muh2_3sig)+100)/1.3).astype('int')


threshs = [20.5,21,20.5,21,20.5,20.5,20.5,21,20.5,20.5,21,21,21,21,21,21,21,20.5,21,20.5,21,21]


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

vmnhi,sigvhi = np.array([spctrum_fit(i)[1:] for i in range(22)]).transpose()
vsnew = ((vmnhi[:,np.newaxis]+np.array([-2,2])[np.newaxis,:]*sigvhi[:,np.newaxis] + 100)/1.3).astype('int')
vsnew[11] = ((vmnhi[11]+np.array([-1,1])*sigvhi[11] + 100)/1.3).astype('int')


mhi = np.array([calc_masses(i)[1] for i in range(22)])

omgh2 = [np.linalg.norm(fit_plane2maps(i,0)[:2]) for i in range(22)]
omghi = [(np.linalg.norm(fit_plane2maps(i,1)[:2])) for i in range(22)]

jh2 = 0.4*(omgh2*u.m/u.s/u.pc*(rh2*u.pc)**2).cgs
jhi = 0.4*(omghi*u.m/u.s/u.pc*(rhi*u.pc)**2).cgs

thth2 = [get_theta(fit_plane2maps(i,0)[:2]) for i in range(22)]
ththi = [get_theta(fit_plane2maps(i,1)[:2]) for i in range(22)]

tab2_newnew = pd.DataFrame(np.around(np.array((mh2/1e4,rh2,np.array(omgh2)/10,thth2,jh2/1e24,\
                                  mhi/1e5,rhi,np.array(omghi)/10,ththi,jhi/1e24)),3).transpose(),
                        columns = 'mh2,rh2,omgh2,thth2,jh2,mhi,rhi,omghi,ththi,jhi'.split(','))


np.around(tab2_newnew,2)


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


from astropy.table import QTable


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


all_rs = np.concatenate((tab2_newnew['rh2'].to_numpy(),tab2_newnew['rhi'].to_numpy(),\
         p99_dat['R<pc>'].to_numpy().astype('float'),chenh_rs,chenx_rs,\
                         goodman['R_eff<pc>'].to_numpy(),tatmtsu_rs,ibb['rh2'],))

all_js = np.concatenate((tab2_newnew['jh2'].to_numpy()*1e24,tab2_newnew['jhi'].to_numpy()*1e24,\
         p99_dat['j<cm2/s>'].to_numpy().astype('float'),chenh_js.value,chenx_js,\
                         goodman['j<cm2/s>'].to_numpy(),tatmtsu_js,ibb['jh2'].to_numpy(),))


plt.rcParams['font.size'] = 23
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10),sharey=True)

jbns=10**np.arange(22,27,1/5)
thbns=np.arange(-180,181,30)

axs[0,0].hist(tab2_newnew['jh2']*1e24,bins=jbns,color=(0.,0.8,0.8),zorder=15,edgecolor='white')
axs[1,0].hist(tab2_newnew['jhi']*1e24,bins=jbns,color=(0.8,0.,0.),zorder=15,edgecolor='white')

for ax in axs[:,0]:
    ax.set_xscale('log')
    ax.set_xlim(1e23,2e26)
    ax.set_ylabel('N (clouds)')
    
for ax in axs[:,1]:
    ax.set_xticks(np.arange(-180,181,90))
    ax.set_xlim(-185,185)
    ax.vlines([-90,90],0,6,linestyle='--',color='black',zorder=-20)
    ax.set_ylim(0,5.1)
    ax.set_yticks([1,3,5])
    

axs[0,1].hist(tab2_newnew['thth2'],bins=thbns,edgecolor='white',color=(0.,0.8,0.8))
axs[1,1].hist(tab2_newnew['ththi'],bins=thbns,edgecolor='white',color=(0.8,0.,0.))


axs[0,0].set_xlabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}$]')
axs[1,0].set_xlabel('$j_\mathrm{HI}\ [\mathrm{cm^2\ s^{-1}}$]')
    
axs[0,1].set_xlabel(r'$\theta_\mathrm{H_2}\ [\degree$]')
axs[1,1].set_xlabel(r'$\theta_\mathrm{HI}\ [\degree$]')
    

fig.tight_layout()

# plt.ylabel('N (clouds)')
# plt.xlim(8e22,3e26)
# plt.ylim(0,5.5)
# # plt.yticks(np.arange(5))

plt.savefig('fig7_js',bbox_inches='tight')

# plt.show()


plt.rcParams['font.size'] = 15


curve_fit(line,np.log10(jh2.value),np.log10(jhi/jh2))


np.sqrt(0.01962425)


np.log10(3e8)


plt.figure(dpi=200,figsize=(5,5))

xs = np.linspace(5e22,1e27,5)

plt.plot(xs,3e8*xs**-0.32,color='black',linestyle='--',
            label=r'$j_\mathrm{HI}\ /\ j_\mathrm{H_2} \propto\ j_\mathrm{H_2}^{-0.32\pm0.14}$')

plt.legend(fontsize='x-small',)

plt.scatter(tab2_newnew['jh2']*1e24,tab2_newnew['jhi']/tab2_newnew['jh2'],
            color='black',edgecolor='white',s=100)


plt.xlabel(r'$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
plt.ylabel(r'$j_\mathrm{HI}\ /\ j_\mathrm{H2}$')

# plt.yticks([1,10],labels=['1','10'])
plt.xlim(1e23,1e26)
plt.hlines([1,10],5e22,1e27,linestyle='-.',color='black')
# plt.ylim(5e-1,5e1)

plt.loglog()
plt.savefig('fig8_jj',bbox_inches='tight')


curve_fit(line,np.log10(tab2_newnew['rh2']),np.log10(jh2.value))


np.sqrt(0.05415548)


curve_fit(line,np.log10(tab2_newnew['rh2'][h2p>0.75]),np.log10(jh2[h2p>0.75].value))


np.sqrt(0.04783445)


curve_fit(line,np.log10(tab2_newnew['rhi'][hip>0.75]),np.log10(jhi[hip>0.75].value))


np.sqrt(0.04885612)


curve_fit(line,np.log10(tab2_newnew['rhi']),np.log10(tab2_newnew['jhi']*1e24))


np.sqrt(0.15357292)


10**22.46


tab2_newnew


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


np.around(omghi,)


plt.rcParams['font.size'] = 18#20


plt.figure(dpi=200,figsize=(5,5))

xs = np.linspace(1e-3,5e2,10)
plt.plot(xs, 1e23*xs**1.42,color='black',linestyle='--',zorder=-10,label='$j \propto R^{1.42 \pm 0.03}$')


plt.scatter(tab2_newnew['rh2'],tab2_newnew['jh2']*1e24,
            color=(0.,0.8,0.8),s=40,edgecolor='white',label='$\mathrm{H_2}$')
plt.scatter(tab2_newnew['rhi'],tab2_newnew['jhi']*1e24,
            color=(0.8,0.,0.),s=40,edgecolor='white',label='HI')


plt.scatter(ibb['rh2'],ibb['jh2'],color='white',
            s=40,edgecolor='black',zorder=-20,label='Imara+ 11')


plt.scatter(p99_dat['R<pc>'],p99_dat['j<cm2/s>'],label='Phillips 99',
            color='black',marker='*',s=90,edgecolor='white')

plt.scatter(chenx_rs,chenx_js,s=40, 
            color='black',marker='<',edgecolor='white',label='Chen+ 07')

plt.scatter(tatmtsu_rs,tatmtsu_js,s=40,
            color='black',marker='o',edgecolor='white',label='Tatematsu+ 16')


plt.scatter(chenh_rs,chenh_js,s=40,
            color='black',marker='s',edgecolor='white',label="Chen+ 19")


# plt.scatter(np.full(len(ibb),70),ibb['jhi'],color='black',alpha=0.5,
#             s=70,edgecolor='red',zorder=30,label='ibb hi')

# plt.scatter(np.full(len(ibb_jhi_lh2),70),ibb_jhi_lh2,color='black',alpha=0.5,
#             s=70,edgecolor='orange',zorder=20,label='ibb hi wo h2')
 
# plt.scatter(goodman['R_eff<pc>'],goodman['j<cm2/s>'],
#             color=(0.6,0.6,0.6),marker='^',s=30,)


plt.xlim(3e-3,4e2)
plt.ylim(2e19,8e26)
plt.loglog()


plt.fill_between(xs,1e23*xs**(1.42+0.03),1e23*xs**(1.42-0.03),
                 alpha=0.2,color='black',linewidth=0.,zorder=-10)


plt.legend(fontsize='xx-small',loc='lower right')#,ncol=2)


plt.xlabel('R [pc]')
plt.ylabel('j [$\mathrm{cm^2\ s^{-1}}$]')

plt.savefig('fig9_jR_all',bbox_inches='tight')


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


def get_tdis(dat):
    egrav_hi = (c.G * dat['mhi'].to_numpy()*1e5*u.solMass / (dat['rhi'].to_numpy()*u.pc)).cgs
    erot_hi = (0.4*dat['jhi'].to_numpy()*1e24*u.cm**2/u.s*dat['omghi'].to_numpy()*u.km/u.s/u.pc).cgs
    ekin_hi = (0.5*(sigvhi*u.km/u.s)**2).cgs

    etot_hi = erot_hi - egrav_hi + ekin_hi
    
    tdis = (dat['jhi']*1e24/etot_hi.value)#*u.s.to(u.Myr)
    
    jinit = (dat['mh2']*1e4*dat['jh2']*1e24 + dat['mhi']*1e5*dat['jhi']*1e24)/(dat['mh2']*1e4 + dat['mhi']*1e5)

    
    return tdis,jinit


tdises,jinits = get_tdis(tab2_newnew)


curve_fit(line,np.log10(tdises*u.s.to(u.Myr)),np.log10(jinits/jh2.value))


np.nanmedian(tdises*u.s.to(u.Myr))


np.nanmean(t_ff(mhi,rhi)*u.s.to(u.Myr))


np.nanmedian(jinits/jh2.value)


5.998/0.68


np.nanmedian(jinits)


1e14*u.s.to(u.Myr)


plt.figure(dpi=200,figsize=(5,5))
plt.scatter(tdises.to_numpy()*u.s.to(u.Myr),jinits/jh2,s=100,color='black',edgecolor='white')

plt.ylabel('$j_\mathrm{init}\ /\ j_\mathrm{H_2}$')
plt.xlabel('$\Delta t_\mathrm{dis}$ [Myr]')
# plt.loglog()

plt.xticks([1,2,3,])
plt.yticks(np.arange(0,41,10))

# plt.show()
# plt.yscale('log')
plt.savefig('fig11_tdis',bbox_inches='tight')


ib['jh2']


tab2_newnew['jh2'][in_ib]*1e24


3.3/18991


3.3/59.8


plt.figure(dpi=200,figsize=(5,5))
plt.scatter(jh2,jhi/jh2,s=100,color='black',edgecolor='white')

plt.plot(np.linspace(7e22,4e26),2e11*np.linspace(8e22,3e26)**-0.43,
         linestyle='--',color='black',
         label='$j_\mathrm{HI}\ /\ j_\mathrm{H_2} \propto j_\mathrm{H_2}^{-0.43}$')

plt.ylabel('$j_\mathrm{HI}\ /\ j_\mathrm{H_2}$')
plt.xlabel('$j_\mathrm{H_2}\ [\mathrm{cm^2\ s^{-1}}]$')
# plt.xlabel('$\Delta t_\mathrm{dis}$ [Myr]')
# plt.yscale('log')
plt.loglog()
plt.xlim(8e22,2e26)
plt.legend(fontsize='x-small')
plt.savefig('jj',bbox_inches='tight')


curve_fit(line,np.log10(jh2.value),np.log10(jhi/jh2))


curve_fit(line,np.log10(jh2.value),np.log10(jhi.value))


delthet = np.array(ththi)-np.array(thth2)
for i in range(22):
    if delthet[i] > 180:
        delthet[i] -= 360
    if delthet[i] < -180:
        delthet[i] += 360


plt.rcParams['font.size'] = 25


plt.figure(dpi=200,figsize=(5,5))
plt.hist(delthet,
        bins=np.arange(-180,181,30),edgecolor='white',color='black')
plt.xticks(np.arange(-180,181,90))
plt.xlabel(r'$\Delta \theta\ [\degree]$')
plt.ylabel('N (clouds)')
plt.yticks([0,1,2,3,4,5])
plt.ylim(0,4.5)
plt.xlim(-185,185)
plt.vlines([-90,90],0,5,linestyle='--',color='black',)

plt.annotate('counter-\nrotating',(-165,3.95),fontsize='xx-small')
plt.annotate('counter-\nrotating',(110,3.95),fontsize='xx-small')
plt.annotate('co-rotating',(-40,4.2),fontsize='xx-small')

# 
plt.savefig('fig4_delthet',bbox_inches='tight')


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

# plt.savefig('jh2s',bbox_inches='tight')

plt.show()


plt.vlines([-90,90],-40,40,zorder=-30,color='black',linestyle='--')
plt.hlines([15,-15],-100,100,zorder=-30,color='black',linestyle='--')


plt.rcParams['font.size'] = 22


plt.figure(dpi=200,figsize=(6,5))
plt.scatter(ththi,bc,c=lc,cmap='PiYG',edgecolor='black',linewidth=.5,s=100)
plt.vlines([-90,90],-40,40,zorder=-30,color='black',linestyle='--')
# plt.hlines([15,-15],-200,200,zorder=-30,color='black',linestyle='--')

plt.colorbar(label='$l\ [\degree]$',ticks=[-180,-90,0,90,180])
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

# plt.hlines(0,-200,200,linestyle='--',zorder=-10)
# plt.vlines(0,-40,40,linestyle='--',color='black',zorder=-10)
plt.ylim(-30,30)
plt.xticks(np.arange(-180,181,90))
plt.xlim(-185,185)
plt.yticks(np.arange(-30,31,15))

plt.savefig('fig12_thetas',bbox_inches='tight')
# 
plt.show()


plt.figure(dpi=200,figsize=(5,5))
plt.scatter(lc,ththi,c=bc,cmap='PiYG',edgecolor='black',linewidth=.5,s=100)


plt.hlines([-90,90],-190,190,zorder=-30,color='black',linestyle='--')
plt.vlines([90,-90],-200,200,zorder=-30,color='black',linestyle='--')

#               
plt.xlabel(r'$l_c\ [\degree]$')
plt.ylabel(r'$\theta_\mathrm{HI}\ [\degree]$')

# plt.hlines(0,-200,200,linestyle='--',zorder=-10)
# plt.vlines(0,-40,40,linestyle='--',color='black',zorder=-10)
plt.xlim(180,-180)
plt.xticks(np.arange(-180,181,90))
plt.yticks(np.arange(-180,181,90))
plt.ylim(-185,185)
plt.xlim(185,-185)

plt.savefig('theta_l',bbox_inches='tight')
# 
plt.show()


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


hi_vos = [fit_plane2maps(i,1)[2]/1e3 for i in range(22)]
h2_vos = [fit_plane2maps(i,0)[2]/1e3 for i in range(22)]


hidvdx = np.array(omghi)*np.cos(ththi*u.deg).value/1e3
hidvdy = np.array(omghi)*np.sin(ththi*u.deg).value/1e3

h2dvdx = np.array(omgh2)*np.cos(thth2*u.deg).value/1e3
h2dvdy = np.array(omgh2)*np.sin(thth2*u.deg).value/1e3


def pearson(xs,ys):
    N = len(xs)
    xmn = np.nanmean(xs)
    ymn = np.nanmean(ys)

    num = np.nansum(xs*ys) - N*xmn*ymn
    den = np.sqrt(np.nansum(xs*xs) - N*xmn**2) * np.sqrt(np.nansum(ys*ys) - N*ymn**2)

    return num/den


def plot_rvsv(i, show=False):  
    ccoords = coordses[i][:4]
    d = ds[i]
    omgs = tab2_newnew[['omgh2','omghi']].to_numpy()[i]/1e3
    thts = tab2_newnew[['thth2','ththi']].to_numpy()[i]
    
    l0c,l1c,b0c,b1c = ccoords
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2maps, himaps = moment_maps(i,threshs[i])
    h2_0,h2_1,h2_2 = h2maps
    hi_0,hi_1,hi_2 = himaps
    
    rhi_px = hi_0.shape[0]/2

    l0h,l1h = lc + rhi_px*0.125, lc - rhi_px*0.125
    b0h,b1h = bc - rhi_px*0.125, bc + rhi_px*0.125
    
    
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

    
    hv_model =  hi_vos[i] + hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc
    cv_model =  h2_vos[i] + h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc


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
        
    
        fig.savefig('v_vs_r_'+names[i],bbox_inches='tight')
        plt.show()

    return h2p, hip
#     return pearson(meancdf['rper'],  meancdf['v']),                   pearson(meanhdf['rper'],  meanhdf['v']),\
#            pearson(sigcdf['v'],      np.abs(meancdf['v']-h2_vos[i])), pearson(sighdf['v'],      np.abs(meanhdf['v']-hi_vos[i])),\
#            pearson(meancdf['w'],     np.abs(meancdf['v']-h2_vos[i])), pearson(meanhdf['w'],     np.abs(meanhdf['v']-hi_vos[i]))


plt.rcParams['font.size'] = 22


hi_vos[16] = 8


[plot_rvsv(i,show=True) for i in range(22)]


h2p,hip = np.array([plot_rvsv(i)[:2] for i in range(22)]).transpose()


curve_fit(line,np.log10(rhi[hip>0.75]),np.log10(jhi.value[hip>0.75]))


np.sqrt(0.0488552 )


curve_fit(line,np.log10(rh2[h2p>0.75]),np.log10(jh2.value[h2p>0.75]))


np.sqrt(0.04783331)


plt.rcParams['font.size'] = 19


fig = plt.figure(dpi=200)#,figsize=(6,6))
ax = fig.add_subplot()

[ax.scatter(np.abs(pval),np.abs(tht),marker='.', s=500,c=c,edgecolor='white',alpha=0.8) 
        for tht,pval,c
        in zip([thth2,thth2],[h2p,hip],[(0.,0.8,0.8),(0.8,0.,0.)])]

ax.vlines([0.5,0.75],0,190,zorder=-30,color='black',linestyle='--')
ax.hlines([60,120],-0.1,1.1,zorder=-30,color='black',linestyle='--')
# plt.xticks(np.arange(0,22,3),labels=names[::3],rotation='vertical')
# plt.grid()

ax.set_ylabel(r'$\theta\ [\degree]$')
ax.set_xlabel(r'$p\ (r_\perp,v)$')
ax.set_ylim(0,185)
ax.set_xlim(-0.02/3,1.02)
ax.set_xticks(np.arange(0.,1.1,0.25),[0,0.25,0.5,0.75,1])
ax.set_yticks(np.arange(0,181,30))


xhist = ax.inset_axes([0, 1.05, 1, 0.25],)

# xhist.hist(ps[0],bins=np.arange(-1.,1.1,0.125),color=magnt,alpha=0.5,edgecolor='black',linewidth=0.5)#,histtype='barstacked')
[xhist.hist(np.abs(pval),bins=np.arange(-0.5,1.1,0.5/8),color=c,label=h,
            edgecolor='white',linewidth=1.5,alpha=0.35)
for pval,c,h
in zip([h2p,hip],[(0.,0.8,0.8),(0.8,0.,0.)],['$\mathrm{H_2}$','$\mathrm{HI}$'])]


xhist.set_xticks(np.arange(0.,1.1,0.25),[])
xhist.set_xlim(-0.02,1.02)
xhist.set_yticks([0,5,10])
xhist.legend(loc='upper left',fontsize='xx-small',)
xhist.set_ylabel('N (clouds)')

# xhist.hlines(10,0,1)
plt.savefig('fig6_theta_vs_p',bbox_inches='tight')

plt.show()


fig = plt.figure(dpi=200)#,figsize=(6,6))
ax = fig.add_subplot()

[ax.scatter(pval,np.abs(tht),marker='.', s=500,c=c,edgecolor='white') 
        for tht,pval,c
        in zip([thth2,thth2],[h2p,hip],[(0.,0.8,0.8),(0.8,0.,0.)])]

ax.vlines([0.5,0.75],0,190,zorder=-30,color='black',linestyle='--')
ax.hlines([60,120],-0.1,1.1,zorder=-30,color='black',linestyle='--')
# plt.xticks(np.arange(0,22,3),labels=names[::3],rotation='vertical')
# plt.grid()

ax.set_ylabel(r'$\theta\ [\degree]$')
ax.set_xlabel(r'$p\ (r_\perp,v)$')
ax.set_ylim(0,185)
ax.set_xlim(-0.02/3,1.02)
ax.set_xticks(np.arange(0.,1.1,0.25),[0,0.25,0.5,0.75,1])
ax.set_yticks(np.arange(0,181,30))


xhist = ax.inset_axes([0, 1.05, 1, 0.25],)

# xhist.hist(ps[0],bins=np.arange(-1.,1.1,0.125),color=magnt,alpha=0.5,edgecolor='black',linewidth=0.5)#,histtype='barstacked')
[xhist.hist(pval,bins=np.arange(-0.5,1.1,0.5/8),color=c,label=h,
            edgecolor='white',linewidth=1.5,alpha=0.35)
for pval,c,h
in zip([h2p,hip],[(0.,0.8,0.8),(0.8,0.,0.)],['$\mathrm{H_2}$','$\mathrm{HI}$'])]


xhist.set_xticks(np.arange(0.,1.1,0.25),[])
xhist.set_xlim(-0.02,1.02)
xhist.set_yticks([0,5,10])
xhist.legend(loc='upper left',fontsize='xx-small',)
xhist.set_ylabel('N (clouds)')

# xhist.hlines(10,0,1)
# plt.savefig('fig6_theta_vs_p',bbox_inches='tight')

plt.show()


np.nanmean(mhi/mh2)


curve_fit(line,np.log10(rh2),np.log10(jh2.value))


np.sqrt(0.05415662)


curve_fit(line,np.log10(rhi),np.log10(jhi.value))


np.sqrt(0.15357622)


plt.scatter(rhi,jhi)
plt.loglog()
plt.xlim(1e1,2e2)
plt.ylim(1e23,2e26)


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
    mo
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


(1e20/u.cm**2*c.m_p*2*1.36).to(u.solMass/u.pc**2)


(10**(20.97)/u.cm**2*c.m_p*1.36).to(u.solMass/u.pc**2)


np.log10(9e20)


10**(20.97)


[(
    plt.figure(figsize=(5,5)),
plt.hist(np.ravel(np.log10(h0)),bins=np.arange(19.5,22.5,0.1),density=True),
    plt.vlines(threshs[i],0,1,color='red'),
# plt.gca().set_aspect(1),
#     plt.yscale('log'),
plt.title(names[i]),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,19)] for h0,h1,h2 in [himaps]]


_here_


threshs = [20.5,21,20.5,21,20.5,20.5,20.5,21,20.5,20.5,21,21,21,21,21,21,21,20.5,21,20.5,21,21]


[(
plt.contourf(np.log10(h0),cmap='magma',levels=np.linspace(threshs[i],np.nanmax(np.log10(h0)),20),
             extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
    plt.xlim(coordses[i][:2]),
plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [h2maps]]


np.nanmedian(jhi)


np.nanmedian(h2p)


[(
plt.contourf(h1,cmap='Spectral',levels=np.linspace(vs_px[i][0]*1.3-100,vs_px[i][1]*1.3-100,30),
                                                   extent=coordses[i][:4]),
plt.gca().set_aspect(1),
plt.title(names[i]),
    plt.colorbar(),
        plt.xlim(coordses[i][:2]),

plt.show()
) for i in range(22) for h2maps,himaps in [moment_maps(i,20.5)] for h0,h1,h2 in [h2maps]]


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
) for i in range(22) for h2maps,himaps in [moment_maps(i,threshs[i])] for h0,h1,h2 in [himaps]]


vs


def draw_spectra(i):
    coords = coordses[i]
    
    loh,lih = lc_p[i]+np.array([-1,1])*rhi_px[i]
    boh,bih = bc_p[i]+np.array([-1,1])*rhi_px[i]
    vo,vi = vs_px[i]
    
    hicube = hi4pi.data[:, boh:bih, loh:lih] #velocity axis length = vi-vo
    h2cube = get_subcube(dht.data,coords[:4]) #velocity axis length = 155
    
    
    hisp_long = np.nanmean(hicube,axis=(1,2))
    hi_spctr = np.nanmean(hicube[vo:vi],axis=(1,2))
    h2_spctr = np.nanmean(h2cube,axis=(1,2))

    Ahi,muhi,stdhi = spctrum_fit(i)
    mody = models.Gaussian1D(Ahi,muhi,stdhi)
    
    
    Ah2,muh2,stdh2 = spctrum_fit(i,which='h2')
    h2mody = models.Gaussian1D(Ah2,muh2,stdh2)
    
    fig,ax = plt.subplots(dpi=200,figsize=(5,5))
    
    ax.plot(np.arange(155)*1.3-100,hisp_long,linewidth=6,color=(0.8,0.,0.))
    ax.vlines(vsnew[i]*1.3-100,0,np.nanmax(hi_spctr)*1.5,linestyle='--', color='black',linewidth=0.8,zorder=-30)
    
    ax.set_xlim(muhi+np.array([-1,1])*stdhi*5)
    ax2 = ax.twinx()
    ax2.plot(np.arange(155)*1.3-100,h2_spctr,color=(0.,0.8,0.8),linewidth=2.5)
              
    ax2.set_ylim(0,1.5*np.nanmax(h2_spctr),)
    
    ax.set_ylim(0,1.3*np.nanmax(hi_spctr))
    
    plt.title(names[i])
    plt.xlabel('$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
    ax.set_ylabel('$T_\mathrm{HI}\ [K]$',color=(0.8,0.,0.))
    ax2.set_ylabel('$T_\mathrm{CO}\ [K]$',color=(0.,0.8,0.8))
    ax.set_xlabel('$v_\mathrm{LSR}\ [\mathrm{km\ s^{-1}}]$')
    
    plt.savefig(names[i][:3]+names[i][-1]+'_spec',bbox_inches='tight')
#     plt.show()
    
#     plt.plot(np.arange(155)*1.3-100,h2_spctr)


plt.rcParams['font.size'] = 23


[(
    draw_spectra(i),
plt.show()
)
for i in range(22)]


h2_vlims = np.array(([2,4],[-10,6],[-1.5,0.5],[14,21],[1,6],[4.8,6.3],[10,15],[1,9],[5,12],[0,4],
                     [3,9],[23,32],[2,12],[9,14],[1,5],
                     [4,12],[6,11],[-3,11],[1,10],[-5,-2],[9,17],[5,7.5]))
# hi_vlims = 
hi_Slims =  np.array(([1,10],[2,60],[4,10],[12,30],[2,15],[2,11],[12,26],[30,63],[4,22],[5,18],
                     [2,40],[10,38],[6,54],[10,28],[1,30],
                     [2,27],[2,30],[4,44],[3,21],[1,15],[20,46],[2,50]))

hi_vlims = np.array([[  0. ,  4 ],
                       [-12. ,   5. ],
                           [ -1.5,   1],
                               [ 14. ,  19.5 ],
                     [  -3. ,   3. ],
                       [  3,   7],
                           [ 6.5 ,  13.5 ],
                                [  2. ,   11. ],
                    [  3. ,  8. ],
                        [  -2.5 ,   1.5 ],
                           [  -1.5 ,   6. ],
                               [ 21. ,  28. ],
                    [  7. ,  12. ],
                        [  5. ,  13. ],
                           [  -1.5 ,   4. ],
                               [  5. ,  12. ],
                    [  5. ,  12. ],
                        [ 3. ,  9. ],
                           [  -3. ,  6. ],
                               [ -8. ,  0.5 ],
                    [  12. ,  17. ],
                       [  -4. ,   7.5]])


def make_zxz(i,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i,threshs[i])
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(13,13))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=np.linspace(10,np.nanmax(h2Sig.value),35),
                   cmap='magma',extend='max',norm='asinh')

        
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
    
    
    h2ox, h2oy, h2vo = fit_plane2maps(i,0)
    hiox, hioy, hivo = fit_plane2maps(i,1)
    
#     axs[1,0].quiver((coords[0]+coords[1])/2,(coords[2]+coords[3])/2,-h2ox,h2oy,scale=1100,headwidth=6,
#                     color='magenta',zorder=30,pivot='middle')
#     axs[1,1].quiver((loh+lih)/2,(boh+bih)/2,-hiox,hioy,color='magenta',
#                     scale=150,headwidth=6,zorder=30,pivot='middle')
    
    
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
                   levels=np.array([-1,0,1])*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),#np.linspace(2,50,5),
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
        
#         ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2],
                           axs.flat[:2],
                          ['',r'$\Sigma \ [\mathrm{M_\odot\ pc^{-2}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3))
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl,norm='asinh')
        
        
    for im,ax,lbl in zip([im3,im4],
                           axs.flat[2:],
                          ['',r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl,norm='asinh')
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


[make_zxz(i,save=True) for i in range(22)]


[make_zxz(i,save=True) for i in range(22)]


longitudes = [r'180$\degree$','','',r'150$\degree$','','',r'120$\degree$','','',\
                    r'90$\degree$','','',r'60$\degree$','','',r'30$\degree$','','',0,\
                    '','',r'-30$\degree$','','',r'-60$\degree$', '','',r'-90$\degree$',\
                    '','',r'-120$\degree$','','',r'-150$\degree$','','',r'-180$\degree$',]

latitudes = ['',r'-40$\degree$','',r'-20$\degree$','',0,'',r'20$\degree$','',r'40$\degree$','']


hi_0th_Msunpc2 = (np.nansum(hi4pi.data,axis=0)*1.3*1.82e18/u.cm**2*c.m_p*1.36).to(u.solMass/u.pc**2).value


plt.rcParams['font.size'] = 30


fig = plt.figure(dpi=200,figsize=(20,10))
ax = fig.add_subplot()#projection='mollweide')

imhi = ax.contourf(hi_0th_Msunpc2[:,1420:4340],extent=[182.5,-182.5,-50,50],
             levels=np.arange(5,110,5),extend='both',cmap='binary',norm='asinh',zorder=10)
# fig.colorbar(imhi,orientation='horizontal',label='$\Sigma_\mathrm{HI}\ [M_\odot\ \mathrm{pc^{-2}}]$')
# plt.colorbar(shrink=0.5,label='$\Sigma_\mathrm{HI} \ [M_\odot \ \mathrm{pc^{-2}}]$')
ax.set_xlim(190,-190)
ax.set_aspect(1)

ax.set_frame_on(False)


[ax.add_patch(Circle((lc[i],bc[i]),rhi_px[i]/8,fill=False,color='black',zorder=20)) for i in range(22)]
[ax.add_patch(Circle((-360+lc[i],bc[i]),rhi_px[i]/8,fill=False,color='black',zorder=20)) for i in [1,21]]


for i in range(22):
    coords = coordses[i]
    h2lb = make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0]

    ax.contourf(np.where(np.log10(h2lb)<threshs[i],np.nan,np.log10(h2lb)),extent=coords[:4],
                   levels=10,#np.linspace(threshs[i],threshs[i]*1.05,9),
                   cmap='magma',extend='max',
                   linewidths=0.75,zorder=30)

# imh2 = ax.contourf((h2lb*c.m_p/u.cm**2*2*1.36).to(u.solMass/u.pc**2).value, extent=coords[:4],
#                    levels=np.linspace(20.4,22,25),
#                    cmap='magma',
#                    linewidths=0.75,zorder=30)

# fig.colorbar(imh2,orientation='horizontal',ticks=np.linspace(20.4,22,25))
    

ax.set_xlim(182,-182)
ax.set_xticks(np.arange(180,-181,-10),
              labels=longitudes)
ax.set_yticks(np.arange(-50,51,10),labels=latitudes)

ax.set_xlabel('Galactic Longitude')
ax.set_ylabel('Galactic Latitude')


# ax.grid(zorder=-30,color='black',alpha=0.5)
plt.savefig('fig1_cat',bbox_inches='tight')

plt.show()


rhoh2 = ((mh2*u.solMass)/(4/3*np.pi*(rh2*u.pc)**3)).cgs
rhohi = ((mhi*u.solMass)/(4/3*np.pi*(rhi*u.pc)**3)).cgs


tffh2 = np.sqrt(3*np.pi/(32*c.G*rhoh2)).to(u.Myr)
tffhi = np.sqrt(3*np.pi/(32*c.G*rhohi)).to(u.Myr)


troth2 = (2*np.pi/(np.array(omgh2)*1e-3*u.km/u.s/u.pc)).to(u.Myr)
trothi = (2*np.pi/(np.array(omghi)*1e-3*u.km/u.s/u.pc)).to(u.Myr)
