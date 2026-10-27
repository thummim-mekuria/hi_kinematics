#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
from scipy.optimize import curve_fit

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datacube as dc

import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
import fun as fns


hi4pi,dht = fits.open('both_same_grid_720d.fits')


dht_3sig = dht.data*(dht.data>0.6)


hmom1 = fns.get_mom1(hi4pi.data,[0,155])


hmom0 = np.nansum(hi4pi.data,axis=0)
cmom0 = np.nansum(dht.data,axis=0)


cmom02 = np.nansum(dht_3sig,axis=0)


hmom1 = fns.get_mom1(hi4pi.data,[0,155])


cmom1 = fns.get_mom1(dht.data,[0,155])


cmom12 = fns.get_mom1(dht_3sig,[0,155])


hmom2 = fns.get_momN(hi4pi.data,[0,155],2)


cmom2 = fns.get_momN(dht.data,[0,155],2)


cmom22 = fns.get_momN(dht_3sig,[0,155],2)


zcat = QTable(fits.open('final_cat.fits')[1].data)


for i in range(17):
   
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm1 = cmom1[bo-10:bi+10,lo-10:li+10]
    cm2 = cmom2[bo-10:bi+10,lo-10:li+10]
    cm0 = cmom0[bo-10:bi+10,lo-10:li+10]

    
    fig,axs = plt.subplots(1,3,figsize=(20,5))
    axs[0].contourf(cm1,levels=20,cmap=dc.cmap_min0('RdYlGn'))
    axs[1].contourf(cm2,levels=20,cmap=dc.cmap_min0('viridis'),norm='asinh')
    axs[2].contourf(cm0,levels=20,cmap=dc.cmap_min0('viridis'))
    axs[2].contour(cm0,levels=[4e20/(1.3*1.97e20)],colors='white')

    print(zcat['names'][i])
    plt.show()


init = models.Linear1D(1.,0.)


fttr = fitting.LinearLSQFitter()


inti = models.Planar2D(1.,1.)


ftrt = fitting.LMLSQFitter() 


for i in range(17):
    if i!=1 and i!=6 and i!=14:
        lo,li = zcat['l_pix_GMC'][i]
        bo,bi = zcat['b_pix_GMC'][i]
        cm0 = cmom0[bo:bi,lo:li]
        cm0 = cm0
        cm1 = cmom1[bo:bi,lo:li] 
        cm2 = cmom2[bo:bi,lo:li]

        cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
        df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
                                   np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
                         columns=['l','b','w','v','s'])
        dffl = df[df['w']>1.56]
#         aw = fttr(init,np.ravel(cl)-np.nanmean(cl),np.ravel(cm1/(cm2**2)))
        a = fttr(init,dffl['l']*-0.125,dffl['v'])

        b = fttr(init,dffl['b']*0.125,dffl['v'])
#         bw = fttr(init,np.ravel(cb)-np.nanmean(cb),np.ravel(cm1/(cm2**2)))

        pfit = ftrt(inti,cl,cb,cm1)

        vnew_l = np.nanmean(dffl['v']) + -0.125*a.slope*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
                              + 0.125*b.slope*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])
        
        vnew_p = np.nanmean(dffl['v']) + pfit.slope_x*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
                              + pfit.slope_y*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])
        
#         print(zcat['names'][i],'\n',a.slope,pfit.slope_x,'\n',b.slope,pfit.slope_y,'\n\n')
        
        plt.title(zcat['names'][i])
#         plt.imshow(cm1*(cm0>1.56),extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50],
#                    origin='lower',cmap=dc.cmap_min0('RdYlGn'),
#                    vmin=100,vmax=101)
# #         plt.show()
#         plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50,
#                 c=dffl['v'],cmap=dc.cmap_min0('RdYlGn')
#                ,vmin=np.nanmin(vnew_p),vmax=np.nanmax(vnew_p),s=5)
# #         plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-(bi-bo)*0.125,
# #                     c=vnew_l,cmap=dc.cmap_min0('RdYlGn')
# #                    ,vmin=np.nanmin(vnew_l),vmax=np.nanmax(vnew_l),s=5)
#         plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-2*(bi-bo)*0.125,
#                     c=vnew_p,cmap=dc.cmap_min0('RdYlGn'),s=5
#                    ,vmin=np.nanmin(vnew_p),vmax=np.nanmax(vnew_p))
#         plt.colorbar()
#         plt.show()
#         plt.scatter(dffl['l']*-0.125+360,(vnew_l-dffl['v'])/dffl['v'],s=3)
#         plt.scatter(dffl['l']*-0.125+360,(vnew_p-dffl['v'])/dffl['v'],s=3)
#         plt.ylim(-1,1)
#         plt.show()
#         plt.scatter(dffl['b']*-0.125+360,(vnew_l-dffl['v'])/dffl['v'],s=3)
#         plt.scatter(dffl['b']*-0.125+360,(vnew_p-dffl['v'])/dffl['v'],s=3)

#         plt.ylim(-1,1)
#         plt.show()
        plt.hist((vnew_l-dffl['v'])/dffl['v'],bins=50,density=True,alpha=0.5)
        plt.hist((vnew_p-dffl['v'])/dffl['v'],bins=50,density=True,alpha=0.5)
        plt.yscale('log')
#         plt.xscale('log')
        plt.show()


#final verdict: go with fitting planar2d


for i in range(17):
    if i!=1 and i!=6 and i!=14:
        lo,li = zcat['l_pix_GMC'][i]
        bo,bi = zcat['b_pix_GMC'][i]
        cm0 = cmom02[bo:bi,lo:li]
        cm1 = cmom12[bo:bi,lo:li] 
        cm2 = cmom22[bo:bi,lo:li]

        cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
        df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
                                   np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
                         columns=['l','b','w','v','s'])
        dffl = df[df['w']>1.56]
#         aw = fttr(init,np.ravel(cl)-np.nanmean(cl),np.ravel(cm1/(cm2**2)))
        a = fttr(init,dffl['l']*-0.125,dffl['v'])

        b = fttr(init,dffl['b']*0.125,dffl['v'])
#         bw = fttr(init,np.ravel(cb)-np.nanmean(cb),np.ravel(cm1/(cm2**2)))

        vnew = np.nanmean(dffl['v']) + -0.125*a.slope*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
                              + 0.125*b.slope*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])
        
        fig,axs=plt.subplots(1,4,figsize=(15,5))
        
        fig.suptitle(zcat['names'][i])
        axs[0].scatter(dffl['b']*0.125-50,dffl['v'],s=5)
        axs[0].plot(dffl['b']*0.125-50,b(dffl['b']*0.125),color='red')
#         axs[0].plot(np.ravel(cb),bw(np.ravel(cb)-np.nanmean(cb)),color='orange')

        axs[1].scatter(dffl['l']*-0.125+360,dffl['v'],s=5)
        axs[1].plot(dffl['l']*-0.125+360,a(dffl['l']*-0.125),color='red')
#         axs[1].plot(np.ravel(cl),aw(np.ravel(cl)-np.nanmean(cl)),color='orange')
        
        
        axs[2].imshow(cm1,origin='lower',cmap=dc.cmap_min0('RdYlGn'))
        axs[3].scatter(dffl['l'],dffl['b'],c=vnew,cmap=dc.cmap_min0('RdYlGn'))
        
        
        Omega = np.sqrt(a.slope**2 + \
                        b.slope**2)/(1*u.deg.to(u.rad)*(zcat['distance<pc>'][i])**2) #turn 0.125 into rad
        #& two factors of one d is to convert slope from deg-> pc, other one is in omega eq 
        theta = np.arctan(b.slope/a.slope)
        
        print(Omega,theta)
        print(a.parameters,b.parameters)
        
        
        plt.show()


for i in range(17):
    if i!=1 and i!=6 and i!=14:
        lo,li = zcat['l_pix_GMC'][i]
        bo,bi = zcat['b_pix_GMC'][i]
        cm0 = cmom0[bo:bi,lo:li]
        cm1 = cmom1[bo:bi,lo:li] 
        cm2 = cmom2[bo:bi,lo:li]

        cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
        df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
                                   np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
                         columns=['l','b','w','v','s'])
        dffl = df[df['w']>1.56]
        a = fttr(init,dffl['l']*-0.125,dffl['v'])

        b = fttr(init,dffl['b']*0.125,dffl['v'])

        vnew = np.nanmean(dffl['v']) + -0.125*a.slope*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
                              + 0.125*b.slope*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])
        
        fig,axs=plt.subplots(1,4,figsize=(15,5))
        
        fig.suptitle(zcat['names'][i])
        axs[0].scatter(dffl['b']*0.125-50,dffl['v'],s=5)
        axs[0].plot(dffl['b']*0.125-50,b(dffl['b']*0.125),color='red')

        axs[1].scatter(dffl['l']*-0.125+360,dffl['v'],s=5)
        axs[1].plot(dffl['l']*-0.125+360,a(dffl['l']*-0.125),color='red')
        
        
        axs[2].imshow(cm1,origin='lower',cmap=dc.cmap_min0('RdYlGn'))
        axs[3].scatter(dffl['l'],dffl['b'],c=vnew,cmap=dc.cmap_min0('RdYlGn'))
        
        
        Omega = np.sqrt(a.slope**2 + \
                        b.slope**2)/(1*u.deg.to(u.rad)*(zcat['distance<pc>'][i])**2) #turn 0.125 into rad
        #& two factors of one d is to convert slope from deg-> pc, other one is in omega eq 
        theta = np.arctan(b.slope/a.slope)
        
        print(Omega,theta)
        print(a.parameters,b.parameters)
        
        
        plt.show()


for i in range(17):
#     if i!=1 and i!=6 and i!=14:
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm1 = cmom1[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]

    cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
                               np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
                     columns=['l','b','w','v','s'])
    dffl = df[df['w']>1.56]
#         aw = fttr(init,np.ravel(cl)-np.nanmean(cl),np.ravel(cm1/(cm2**2)))
    a = fttr(init,dffl['l']*-0.125,dffl['v'])

    b = fttr(init,dffl['b']*0.125,dffl['v'])
#         bw = fttr(init,np.ravel(cb)-np.nanmean(cb),np.ravel(cm1/(cm2**2)))

    if i in [1,6,14]:
        pfit = ftrt(inti,dffl['l'],dffl['b'],dffl['v'])
    else:
        pfit = ftrt(inti,cl,cb,cm1)

    vnew_l = np.nanmean(dffl['v']) + -0.125*a.slope*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
                          + 0.125*b.slope*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])

    vnew_p = np.nanmean(dffl['v']) + pfit.slope_x*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
                          + pfit.slope_y*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])

#         print(zcat['names'][i],'\n',a.slope,pfit.slope_x,'\n',b.slope,pfit.slope_y,'\n\n')

    plt.title(zcat['names'][i])
    plt.imshow(cm1*(cm0>1.56),extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50],
               origin='lower',cmap=dc.cmap_min0('RdYlGn'),
               vmin=100,vmax=101)
#         plt.show()
    plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50,
            c=dffl['v'],cmap=dc.cmap_min0('RdYlGn')
           ,vmin=np.nanmin(vnew_p),vmax=np.nanmax(vnew_p),s=5)
    
    plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-(bi-bo)*0.125,
                c=vnew_l,cmap=dc.cmap_min0('RdYlGn')
               ,vmin=np.nanmin(vnew_l),vmax=np.nanmax(vnew_l),s=5)
    plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-2*(bi-bo)*0.125,
                c=vnew_p,cmap=dc.cmap_min0('RdYlGn'),s=5
               ,vmin=np.nanmin(vnew_p),vmax=np.nanmax(vnew_p))
    plt.colorbar()
    plt.show()


#         plt.scatter(dffl['l']*-0.125+360,(vnew_l-dffl['v'])/dffl['v'],s=3)
#         plt.scatter(dffl['l']*-0.125+360,(vnew_p-dffl['v'])/dffl['v'],s=3)
#         plt.ylim(-1,1)
#         plt.show()
#         plt.scatter(dffl['b']*-0.125+360,(vnew_l-dffl['v'])/dffl['v'],s=3)
#         plt.scatter(dffl['b']*-0.125+360,(vnew_p-dffl['v'])/dffl['v'],s=3)
#         plt.ylim(-1,1)
#         plt.show()


#     plt.hist((vnew_l-dffl['v'])/dffl['v'],bins=50,density=True,alpha=0.5)
#     plt.hist((vnew_p-dffl['v'])/dffl['v'],bins=50,density=True,alpha=0.5)
#     plt.yscale('log')
# #         plt.xscale('log')
#     plt.show()


for i in range(17):
#     if i!=1 and i!=6 and i!=14:
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm0 = cm0
    cm1 = cmom1[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]

    cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
            np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    delvdell = fttr(init,dffl['l']*-0.125,dffl['v'])
    delvdelb = fttr(init,dffl['b']*0.125,dffl['v'])

    if i in [1,6,14]:
        pfit = ftrt(inti,dffl['l']*-0.125,dffl['b']*0.125,dffl['v'])
    else:
        pfit = ftrt(inti,cl*-0.125,cb*0.125,cm1)
        
    dvdl = pfit.slope_x
    dvdb = pfit.slope_y
    
    
    vnew_l = np.nanmean(dffl['v']) +\
            -0.125*delvdell.slope*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])+\
             0.125*delvdelb.slope*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])
    #the 0.125 is for dffl['l'] etc to go from pixel -> degree units and match the slopes

#     vnew_p = np.nanmean(dffl['v']) + \
#             dvdl*(dffl['l']-dffl['l'].iloc[int(len(dffl)/2)])\
#             + dvdb*(dffl['b']-dffl['b'].iloc[int(len(dffl)/2)])
    
    vnew_p = np.nanmean(dffl['v']) + \
            -0.125*dvdl*(dffl['l']-dffl['l'].iloc[0])\
            + 0.125*dvdb*(dffl['b']-dffl['b'].iloc[0])
    
    plt.title(zcat['names'][i])
    plt.imshow(cm1*(cm0>1.56),extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50],
               origin='lower',cmap=dc.cmap_min0('RdYlGn'),
               vmin=100,vmax=101)
#         plt.show()
    plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50,
            c=dffl['v'],cmap=dc.cmap_min0('RdYlGn'),s=5)
#     plt.colorbar()
#     plt.show()
    
    plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-(bi-bo)*0.125,
                c=vnew_l,cmap=dc.cmap_min0('RdYlGn'),s=5)
#     plt.colorbar()
#     plt.show()
    plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-2*(bi-bo)*0.125,
                c=vnew_p,cmap=dc.cmap_min0('RdYlGn'),s=5)
    plt.colorbar()
    plt.show()


#     print(zcat['names'][i],':',delvdell.slope.value,dvdl.value,'\n',
#           delvdelb.slope.value,dvdb.value,'\n\n\n')
    
    
    theta_p = np.arctan(dvdb/dvdl)*u.rad.to(u.deg)
    theta_l = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    rotax_p = theta_p-90
    rotax_l = theta_l-90
    
    
    print(zcat['names'][i],dvdl.value,dvdb.value,theta_p,'\n',\
          delvdell.slope.value,delvdelb.slope.value,theta_l,'\n\n',)


#gonna go with fitting l & b separately
for i in range(17):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm0 = cm0
    cm1 = cmom1[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]

    cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
            np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    lpc = ((dffl['l']-np.nanmean(dffl['l']))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((dffl['b']-np.nanmean(dffl['b']))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    delvdell = fttr(init,lpc,dffl['v'])
    delvdelb = fttr(init,bpc,dffl['v'])

    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    #account for similar angles, based on the signs of dx and dy
    if delvdell.slope > 0 : #first & fourth quadrants
        theta = theta
    elif delvdell.slope < 0: #second quadrant & third
        theta = theta + 180
        
    if theta > 180:
        theta -=360
        
        
#     plt.title(zcat['names'][i])
#     plt.imshow(cm1*(cm0>1.56),
#                extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50],
#                origin='lower',cmap=dc.cmap_min0('RdYlGn'))
#     plt.xlim(li*-0.125+360,lo*-0.125+360,)
#     plt.show()
#     plt.plot(dffl['l']*-0.125+360,(dffl['l']*-0.125+360)/np.cos(theta_p*u.deg))
#     plt.show()
    print(zcat['names'][i], Omega)
    print(delvdell.slope.value,delvdelb.slope.value,theta,'\n\n')


#gonna go with fitting l & b separately
Omegas = np.zeros(17)
thetas = np.zeros(17)
for i in range(17):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm1 = cmom1[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]

    cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
            np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    lpc = ((dffl['l']-np.nanmean(dffl['l']))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((dffl['b']-np.nanmean(dffl['b']))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    delvdell = fttr(init,lpc,dffl['v'])
    delvdelb = fttr(init,bpc,dffl['v'])
    
    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
#     #account for similar angles, based on the signs of dx and dy
#     if delvdell.slope > 0 : #first & fourth quadrants
#         theta = theta
#     elif delvdell.slope < 0: #second quadrant & third
#         theta = theta + 180
        
#     if theta > 180:
#         theta -=360
        
#     Omegas[i] = Omega
#     thetas[i] = theta
        
    plt.scatter(((dffl['l']-np.nanmean(np.unique(dffl['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i], 
            ((dffl['b']-np.nanmean(np.unique(dffl['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i],
            c=dffl['v'],alpha=0.6,cmap='RdYlGn')
    
    plt.plot(70*np.arange(-Omega,Omega,Omega/10)*np.cos(-90*u.deg+theta*u.deg),
         70*np.arange(-Omega,Omega,Omega/10)*np.sin(-90*u.deg+theta*u.deg),
             color='black',linestyle='--') #rotation axis, 90 degree away from the velocity gradient vector 
    
    plt.show()
  

#gonna go with fitting l & b separately
for i in range(17):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom02[bo:bi,lo:li]
    cm1 = cmom12[bo:bi,lo:li] 
    cm2 = cmom22[bo:bi,lo:li]

    cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
            np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    lpc = ((dffl['l']-np.nanmean(dffl['l']))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((dffl['b']-np.nanmean(dffl['b']))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    delvdell = fttr(init,lpc,dffl['v'])
    delvdelb = fttr(init,bpc,dffl['v'])
    
    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    #account for similar angles, based on the signs of dx and dy
    if delvdell.slope > 0 : #first & fourth quadrants
        theta = theta
    elif delvdell.slope < 0: #second quadrant & third
        theta = theta + 180
        
    if theta > 180:
        theta -=360
        
    Omegas[i] = Omega
    thetas[i] = theta
        
    plt.scatter(((dffl['l']-np.nanmean(np.unique(dffl['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i], 
            ((dffl['b']-np.nanmean(np.unique(dffl['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i],
            c=dffl['v'],alpha=0.6,cmap='RdYlGn')
    
    plt.plot(100*np.arange(-Omega,Omega,Omega/10)*np.cos(-90*u.deg+theta*u.deg),
         100*np.arange(-Omega,Omega,Omega/10)*np.sin(-90*u.deg+theta*u.deg),
             color='black',linestyle='--') #rotation axis, 90 degree away from the velocity gradient vector 
    
    plt.show()
  

#gonna go with fitting l & b separately
for i in range(17):
    loh,lih = zcat['l_pix_HI'][i]
    boh,bih = zcat['b_pix_HI'][i]
    hm0 = hmom0[bo:bi,lo:li]
    hm1 = hmom1[bo:bi,lo:li] 
    hm2 = hmom2[bo:bi,lo:li]
    
    lo,li = zcat['l_pix_HI'][i]
    bo,bi = zcat['b_pix_HI'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm1 = cmom12[bo:bi,lo:li]
    
    
    hl,hb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(hl),np.ravel(hb),
            np.ravel(hm0),np.ravel(hm1),np.ravel(hm2))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    lpc = ((dffl['l']-np.nanmean(dffl['l']))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((dffl['b']-np.nanmean(dffl['b']))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    delvdell = fttr(init,lpc,dffl['v'])
    delvdelb = fttr(init,bpc,dffl['v'])
    
    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    #account for similar angles, based on the signs of dx and dy
    if delvdell.slope > 0 : #first & fourth quadrants
        theta = theta
    elif delvdell.slope < 0: #second quadrant & third
        theta = theta + 180
        
    if theta > 180:
        theta -=360
        
#     Omegas[i] = Omega
#     thetas[i] = theta
        
        
    plt.title(zcat['names'][i])
    plt.imshow(hm1*(hm0>1.56),
               extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50],
               origin='lower',cmap=dc.cmap_min0('RdYlGn'))
    plt.colorbar()
    plt.contourf(np.arange(lo,li)*-0.125+360,np.arange(bo,bi)*0.125-50,
                cm1,cmap='binary_r',alpha=0.9)
    plt.colorbar()

#     plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50,
#             c=dffl['v'],cmap=dc.cmap_min0('RdYlGn'),s=5)
    
#     plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-(bi-bo)*0.125,
#                 c=vnew_l,cmap=dc.cmap_min0('RdYlGn'),s=5)

    
#     plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-2*(bi-bo)*0.125,
#                 c=vnew_p,cmap=dc.cmap_min0('RdYlGn'),s=5)
   
    plt.show()  
#     plt.scatter(((dffl['l']-np.nanmean(np.unique(dffl['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i], 
#             ((dffl['b']-np.nanmean(np.unique(dffl['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i],
#             c=dffl['v'],alpha=0.6,cmap='RdYlGn')
    
#     plt.plot(100*np.arange(-Omega,Omega,Omega/10)*np.cos(-90*u.deg+theta*u.deg),
#          100*np.arange(-Omega,Omega,Omega/10)*np.sin(-90*u.deg+theta*u.deg),
#              color='black',linestyle='--') #rotation axis, 90 degree away from the velocity gradient vector 
    
#     plt.show()
  

#gonna go with fitting l & b separately
for i in range(17):
    lo,li = zcat['l_pix_HI'][i]
    bo,bi = zcat['b_pix_HI'][i]
    hm0 = hmom0[bo:bi,lo:li]
    cm0 = cmom0[bo:bi,lo:li]

    hm1 = hmom1[bo:bi,lo:li] 
    cm1 = cmom12[bo:bi,lo:li]
    
    hm2 = hmom2[bo:bi,lo:li]

    hl,hb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(hl),np.ravel(hb),
            np.ravel(hm0),np.ravel(hm1),np.ravel(hm2))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    lpc = ((dffl['l']-np.nanmean(dffl['l']))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((dffl['b']-np.nanmean(dffl['b']))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    delvdell = fttr(init,lpc,dffl['v'])
    delvdelb = fttr(init,bpc,dffl['v'])
    
    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    #account for similar angles, based on the signs of dx and dy
    if delvdell.slope > 0 : #first & fourth quadrants
        theta = theta
    elif delvdell.slope < 0: #second quadrant & third
        theta = theta + 180
        
    if theta > 180:
        theta -=360
        
#     Omegas[i] = Omega
#     thetas[i] = theta
        
        
    plt.title(zcat['names'][i])
    plt.imshow(hm1*(hm0>1.56),
               extent=[lo*-0.125+360,li*-0.125+360,bo*0.125-50,bi*0.125-50],
               origin='lower',cmap=dc.cmap_min0('RdYlGn'))
    plt.colorbar()
    plt.contourf(np.arange(lo,li)*-0.125+360,np.arange(bo,bi)*0.125-50,
                cm1,cmap='binary_r',alpha=0.9)
    plt.colorbar()

#     plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50,
#             c=dffl['v'],cmap=dc.cmap_min0('RdYlGn'),s=5)
    
#     plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-(bi-bo)*0.125,
#                 c=vnew_l,cmap=dc.cmap_min0('RdYlGn'),s=5)

    
#     plt.scatter(dffl['l']*-0.125+360,dffl['b']*0.125-50-2*(bi-bo)*0.125,
#                 c=vnew_p,cmap=dc.cmap_min0('RdYlGn'),s=5)
   
    plt.show()  
#     plt.scatter(((dffl['l']-np.nanmean(np.unique(dffl['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i], 
#             ((dffl['b']-np.nanmean(np.unique(dffl['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i],
#             c=dffl['v'],alpha=0.6,cmap='RdYlGn')
    
#     plt.plot(100*np.arange(-Omega,Omega,Omega/10)*np.cos(-90*u.deg+theta*u.deg),
#          100*np.arange(-Omega,Omega,Omega/10)*np.sin(-90*u.deg+theta*u.deg),
#              color='black',linestyle='--') #rotation axis, 90 degree away from the velocity gradient vector 
    
#     plt.show()
  

zcat.write('final_cat.fits',overwrite=True)


for i in range(16):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm1 = cmom1[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]


#gonna go with fitting l & b separately
for i in range(17):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm1 = cmom1[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]

    cl,cb = np.meshgrid(np.arange(lo,li),np.arange(bo,bi))
    df = pd.DataFrame(np.array((np.ravel(cl),np.ravel(cb),
            np.ravel(cm0),np.ravel(cm1),np.ravel(cm1))).transpose(),
            columns=['l','b','w','v','s'])
    
    dffl = df[df['w']>1.56]
    
    lpc = ((dffl['l']-np.nanmean(dffl['l']))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]
    bpc = ((dffl['b']-np.nanmean(dffl['b']))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i]

    delvdell = fttr(init,lpc,dffl['v'])
    delvdelb = fttr(init,bpc,dffl['v'])
    
#     plt.title(zcat['names'][i])
    plt.scatter(lpc,dffl['v'],s=1)
    plt.plot(lpc,delvdell(lpc),color='red')
    plt.show()

    Omega = np.sqrt(delvdell.slope**2 + delvdelb.slope**2) 
    
    theta = np.arctan(delvdelb.slope/delvdell.slope)*u.rad.to(u.deg)
    
    #account for similar angles, based on the signs of dx and dy
    if delvdell.slope > 0 : #first & fourth quadrants
        theta = theta
    elif delvdell.slope < 0: #second quadrant & third
        theta = theta + 180
        
    if theta > 180:
        theta -=360
        
    Omegas[i] = Omega
    thetas[i] = theta
        
    plt.scatter(((dffl['l']-np.nanmean(np.unique(dffl['l'])))*-0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i], 
            ((dffl['b']-np.nanmean(np.unique(dffl['b'])))*0.125)*u.deg.to(u.rad)*zcat['distance<pc>'][i],
            c=dffl['v'],alpha=0.6,cmap='RdYlGn')
    
    plt.plot(70*np.arange(-Omega,Omega,Omega/10)*np.cos(-90*u.deg+theta*u.deg),
         70*np.arange(-Omega,Omega,Omega/10)*np.sin(-90*u.deg+theta*u.deg),
             color='black',linestyle='--') #rotation axis, 90 degree away from the velocity gradient vector 
    
    plt.show()
  

for i in range(17):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom02[bo:bi,lo:li]
    plt.scatter(zcat['names'][i],len(cm0[cm0>1.56]))
    
plt.xticks(rotation='vertical')
plt.yscale('log')
plt.hlines(100,'california','taurus')
plt.show()


plt.scatter(zcat['names'],(thetas-np.array(zcat['vel_grad_GMC<deg>']\
                )))
plt.xticks(rotation='vertical')
plt.hlines(0,'california','taurus')
plt.show()


zcat['vel_grad_GMC<km/s/pc>'] = Omegas
zcat['vel_grad_GMC<deg>'] = thetas


zcat[13]


zcat[10]


plt.contourf(cm0)
plt.colorbar()


plt.contourf(cm1)
plt.colorbar()


plt.contourf(cm2,levels=np.arange(0,10,2))
plt.colorbar()
plt.contour(cm0,linewidths=0.8,levels=[2*1.56],colors='white')


for i in range(17):
    lo,li = zcat['l_pix_GMC'][i]
    bo,bi = zcat['b_pix_GMC'][i]
    cm0 = cmom0[bo:bi,lo:li]
    cm02 = cmom02[bo:bi,lo:li] 
    cm2 = cmom2[bo:bi,lo:li]
    cm22 = cmom22[bo:bi,lo:li]
    
    fig,axs=plt.subplots(1,2,figsize=(10,5))
    axs[0].contourf(cm2,levels=np.arange(0,10,1))
#     plt.colorbar()
    axs[0].contour(cm0,linewidths=0.8,levels=[2*1.56],colors='white')
#     plt.scatter(95,10,color='red',marker='x')
#     plt.scatter(44,30,color='red',marker='x')
#     plt.scatter(75,30,color='red',marker='x')


    axs[1].contourf(cm22,levels=np.arange(0,10,2))
#     plt.colorbar()
    axs[1].contour(cm02,linewidths=0.8,levels=[2*1.56],colors='white')
    plt.show()

    
ros = dht.data[:,bo:bi,lo:li]
ros2 = ros*(ros>0.6)
cm12 = cmom12[bo:bi,lo:li]
cm22 = cmom22[bo:bi,lo:li]
dvdy = np.gradient(cm12)[0]
dvdx = np.gradient(cm12)[1]


plt.imshow(cm12,origin='lower',norm='asinh')


plt.contourf(dvdx,levels=np.arange(-5,5,0.5),cmap='RdYlGn')
plt.colorbar()


plt.contourf(dvdy,levels=np.arange(-5,5,0.5),cmap='RdYlGn')
plt.colorbar()


plt.contourf(np.arctan(dvdy/dvdx)*u.rad.to(u.deg))
plt.colorbar()


plt.contourf(np.sqrt(dvdx**2 + dvdy**2),levels=np.arange(-5,5,0.5))
plt.colorbar()


cm12 = cmom12[bo:bi,lo:li]
cm22 = cmom22[bo:bi,lo:li]
dvdy = np.gradient(cm12)[0]
dvdx = np.gradient(cm12)[1]


plt.imshow(cm12,origin='lower')
plt.colorbar()


x,y=np.meshgrid(np.arange(96),np.arange(32))


plt.quiver(x,y,dvdx,dvdy,np.sqrt(dvdx**2 + dvdy**2),angles='xy')


plt.scatter(np.ravel(x),np.ravel(cm12),s=1)
plt.scatter(np.ravel(x),np.ravel(dvdx),s=1)


plt.contourf(dvdx,levels=np.arange(-5,5,0.5),cmap='RdYlGn')
plt.colorbar()


plt.contourf(dvdy,levels=np.arange(-5,5,0.5),cmap='RdYlGn')
plt.colorbar()


plt.contourf(np.sqrt(dvdx**2 + dvdy**2),levels=np.arange(0,3,0.05))
plt.colorbar()


plt.figure(figsize=(10,5))
plt.contourf(np.arctan(dvdy/dvdx)*u.rad.to(u.deg))
plt.colorbar()


cm12 = cmom12[bo:bi,lo:li]
cm22 = cmom22[bo:bi,lo:li]
dvdy = np.gradient(cm12)[0]
dvdx = np.gradient(cm12)[1]


plt.imshow(cm12,origin='lower')
plt.colorbar()


x,y=np.meshgrid(np.arange(56),np.arange(56))


plt.quiver(x,y,dvdx,dvdy,np.sqrt(dvdx**2 + dvdy**2),angles='xy')


plt.figure(figsize=(10,5))
plt.contourf(np.arctan(dvdy/dvdx)*u.rad.to(u.deg))
plt.colorbar()


plt.contourf(np.sqrt(dvdx**2 + dvdy**2),levels=np.arange(0,3,0.05))
plt.colorbar()


cm12 = cmom12[bo:bi,lo:li]
cm22 = cmom22[bo:bi,lo:li]
dvdy = np.gradient(cm12)[0]
dvdx = np.gradient(cm12)[1]


plt.imshow(cm12,origin='lower')
plt.colorbar()


x,y=np.meshgrid(np.arange(104),np.arange(48))


plt.quiver(x,y,dvdx,dvdy,np.arctan(dvdy/dvdx)*u.rad.to(u.deg),angles='xy',scale=30)
plt.colorbar()


plt.figure(figsize=(10,5))
plt.contourf(np.arcsin(dvdy/np.sqrt(dvdx**2 + dvdy**2))*u.rad.to(u.deg))
plt.colorbar()
