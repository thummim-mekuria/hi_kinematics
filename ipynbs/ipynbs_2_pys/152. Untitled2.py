#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from   astropy.io import fits
import astropy.units as u
import astropy.constants as c

from scipy.optimize import curve_fit

from joblib import Parallel, delayed


table1 = pd.read_csv('/Users/thummimm/Downloads/molecular_clouds.csv')


ds = table1['d (pc)'].to_numpy()


coordses = table1[['l0 (°)','l1 (°)',
                   'b0 (°)','b1 (°)',
                   'v0 (km/s)','v1 (km/s)']].to_numpy()


hi4pi, dht = fits.open('../fitses/both_same_grid_720d.fits')


def first_moment(cube,arr):
    zeroth = np.nansum(cube,axis=0)
    
    num = np.nansum(cube * arr[:,np.newaxis,np.newaxis],axis=0)
    return num/zeroth, np.log10(zeroth*1.3*1.97e21)


def get_subcube(bigcube, coords, dx=-0.125, dy=0.125, x0=360, y0=-50):
    lo,li,bo,bi = coords # assume in degrees
    xo,xi = ((np.array([lo,li]) - x0)/dx).astype('int')
    yo,yi = ((np.array([bo,bi]) - y0)/dy).astype('int')
    
    return bigcube[:,yo:yi,xo:xi]


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


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
    return lb,lv,bv


for coords in coordses:
    make_lb_lv_bv(get_subcube(dht.data,coords[:4]))


threshs = [20. , 20.7, 20. , 20.5, 20.5, 20.5, 20.5, 20.8, 20.5, 20. , 20.5,
       20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20. , 20.5, 20.5, 20.7, 21.1]


for i,coords in enumerate(coordses):
    
    first,zeroth = first_moment(get_subcube(dht.data,coords[:4]),np.arange(155)*1.3-100)

    fig,axs = plt.subplots(1,2,dpi=200)
    axs[0].contourf(np.where(zeroth<21.5,np.nan,first),extent=coords[:4],
                 cmap='Spectral',levels=np.linspace(coords[4]+2.6,-2.6+coords[5],30))
    
    axs[1].contourf(zeroth,extent=coords[:4],
                 cmap='magma',levels=np.linspace(threshs[i],22.5,30))
#     axs[1].contour(zeroth,extent=coords[:4],
#                  colors='black',levels=[21.5])
    
#     plt.colorbar(shrink=0.7)
    for ax in axs:
        ax.set_xlim(coords[:2])
        ax.set_aspect(1)
    plt.show()


def calc_mass(cube,d,thresh=20.6,vo=0,vi=155):
    zeroth = np.nansum(cube,axis=0)*1.3*1.97e20/u.cm**2
    
    mask = np.log10(zeroth.value)>thresh
    
    
#     plt.contourf(np.log10(zeroth.value)*mask,levels=np.linspace(threshs[i],22.5,30))
#     plt.gca().set_aspect(1)
#     plt.show()
    
    masked = zeroth*mask
    dA = (d*u.pc*0.125*u.deg.to(u.rad))**2
    area = np.nansum(masked/masked)*dA
    
    mass = np.nansum(masked*dA)*c.m_p*1.36*2
    return mass.to(u.solMass).value, np.sqrt(area/np.pi).value
    
    
mh2, rh2 = pd.read_csv('/Users/thummimm/Downloads/coordinate_table.csv')[['X','Y']].to_numpy().transpose()
rhi = (((3*mh2*u.solMass*1e3)/(4*np.pi*1.36*c.m_p*0.57/u.cm**3))**(1/3)).to(u.pc).value
rhi_px = (np.arctan(rhi/ds)*u.rad.to(u.deg)/0.125).astype('int')


lc,bc = np.nanmean(coordses[:,:2],axis=1),np.nanmean(coordses[:,2:4],axis=1)


names = table1['Name']


#this is in place of threshs (errors above)
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


def make_hi_maps(i):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])
    
     
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]

    
    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    axs[0].contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="YlGn_r",
                    extend='both',
                       )
    
    
    axs[0].contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="YlGn_r",
                    extend='both',
                       )    
    
        
    axs[1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='bwr_r',extend='both')
    
    
    axs[1].annotate(str(np.around(hi_vlims[i][0],1))+\
                    r'$\rightarrow$'+\
                    str(np.around(hi_vlims[i][1],1)),
                    (lih+rhi_px[i]/17,bih-rhi_px[i]/70),
                    fontsize='small')
    
    axs[0].annotate(str(np.around(hi_Nlims[i][0],1))+\
                    r'$\rightarrow$'+\
                    str(np.around(hi_Nlims[i][1],1)),
                    (lih+rhi_px[i]/17,bih-rhi_px[i]/70),
                    fontsize='small')
    
    axs[0].annotate(names[i],(loh-rhi_px[i]/100,bih-rhi_px[i]/70))
    
    
    for ax in axs:
        ax.set_aspect(1)
        ax.contour(h2lb, extent=coords[:4],
                   levels=np.linspace(threshs[i]+0.1,22,5),
                   colors='black',
                   linewidths=0.75,)
        
        
        hix,hiy = np.indices(hilb.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],levels=[rhi_px[i],2*rhi_px[i]],colors='white'
                       )
        
        ax.set_xlim(loh,lih)
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
    
        ax.tick_params(axis='both', which='both', length=0)
                
        fig.suptitle(names[i])
        fig.tight_layout()

        
#         ax.spines.clear()
        
        
#     axs[2].vlines(20.,0,5000,color='red')
#     axs[2].set_yscale('log')
    plt.savefig(names[i][:3]+names[i][-1]+'_0th_1st',bbox_inches='tight')
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


h2_Nlims = np.array((threshs,[21.25, 21.9, 21.4, 22., 21.7, 21.7, 21.7, 22., 21.7,
       21.75, 21.8 , 21.7, 22.2, 22., 22.2, 22.5, 22.5, 22.,
       22.2, 21.5, 21.8, 21.8])).transpose()


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


def make_zxz(i,save=False):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
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


[make_zxz(i,save=True) for i in range(22)]


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


# h2_fitpars = np.array([fit_gradient(h2_1s[i].data,ds[i],h2_msk[i],h2_2s[i].data) for i in range(22)])
# hi_fitpars = np.array([fit_gradient(hi_1s[i].data,ds[i],hi_msk[i],hi_2s[i].data) for i in range(22)])


def plot_spec(ax,cube,coords):
    lo,li,bo,bi = coords
    
    sp = np.nanmean(cube[:,bo:bi,lo:li],axis=(1,2))
    
    ax.plot(np.arange(155)*1.3-100,(sp-np.nanmin(sp))/np.nanmax(sp),linewidth=0.7,color='black')
    
    ax.set_xticks([]),
    ax.set_xlim(vos[i]-10,vis[i]+10),
    ax.vlines(np.nanmean((vos[i],vis[i])),0,1,linestyle='--',linewidth=0.5,color='black',)
    ax.set_yticks([])


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


[make_pbp_spectra(i) for i in range(22)]


make_pbp_spectra(8)


def calc_os_js(i,save=False):
    coords = coordses[i]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
 
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    hicube = get_subcube(hi4pi.data,[loh,lih,boh,bih])
    h2cube = get_subcube(dht.data,coords[:4])
    
    h2lb = make_lb_lv_bv(h2cube)[0]
    hilb = make_lb_lv_bv(hicube)[0]
        
    h21 = first_moment(h2cube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    hix,hiy = np.indices(hi1.shape)
    hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
    
    h21 = h21*(h21>h2_vlims[i][0])*(h21<h2_vlims[i][1])*(h2lb>h2_Nlims[i][0])
    hi1 = hi1*(hi1>hi_vlims[i][0])*(hi1<hi_vlims[i][1])*(hicrc<rhi_px[i])

    h2ox,h2oy,h2vo = fit_gradient(h21,ds[i])
    hiox,hioy,hivo = fit_gradient(hi1,ds[i])
    
    h2o = np.linalg.norm((h2ox,h2oy))
    hio = np.linalg.norm((hiox,hioy))
    
    mhi = calc_mass(hicube[vo:vi],ds[i],thresh=hi_Nlims[i][0])[0] 
    
    
    return mh2[i]*1e3,rh2[i],h2o,mhi,rhi[i],hio

#want to also have ^ calculate + return theta, but thats for later


mh2,rh2,h2o,mhi,rhi,hio = np.array([calc_os_js(i) for i in range(22)]).reshape(22,6).transpose()


def fitline(x,y):
    fit = curve_fit(line,x,y)[0]
    std = np.sqrt(np.diag(curve_fit(line,x,y)[1]))
    
    return fit,std


jh2 = h2o*(u.km/u.s/u.pc).to(1/u.s)*(rh2*u.pc.to(u.cm))**2
jhi = hio*(u.km/u.s/u.pc).to(1/u.s)*(rhi*u.pc.to(u.cm))**2


def make_jR_jM ():
    jRh2_fit = np.array(fitline(np.log10(rh2),np.log10(jh2))).transpose()
    jRhi_fit = np.array(fitline(np.log10(rhi),np.log10(jhi))).transpose()
    
#     jMh2_fit = np.array(fitline(np.log10(mh2),np.log10(jh2))).transpose()
#     jMhi_fit = np.array(fitline(np.log10(mhi),np.log10(jhi))).transpose()


    fig,axs = plt.subplots(1,2,dpi=200,figsize=(10,5))
    for ax, x, y, fitpars, clr, sub in\
    zip(axs,[rh2,rhi],[jh2,jhi],[jRh2_fit,jRhi_fit], ['magenta','green']):
        ax.scatter(x,y,color='black')
        ax.plot(x,(x**fitpars[0,0])*(10**fitpars[1,0]),color='red')
        
        ax.plot(x,(x**(fitpars[0,0]+fitpars[0,1]))*(10**fitpars[1,0]),color='pink',alpha=0.5)
        ax.plot(x,(x**(fitpars[0,0]-fitpars[0,1]))*(10**fitpars[1,0]),color='pink',alpha=0.5)
        ax.loglog()
        
    fig.supxlabel('R [pc]')
    fig.supylabel(r'j [$\mathrm{cm^2\ s^{-1}}$]')
    
    fig.tight_layout()


make_jR_jM()


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
