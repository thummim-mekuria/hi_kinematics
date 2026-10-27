#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
from astropy.table import QTable

import numpy as np
import matplotlib.pyplot as plt


hi4pi,dht = fits.open('both_same_grid_720d.fits')
zcat = QTable(fits.open('final_cat.fits')[1].data)


d_unc = np.sqrt((np.array([24,9,60,44,89,11,25,104,37,39,22,21,20,15,18,65,7]))**2 +\
(np.array([2,3,4,4,4,1,5,6,3,3,2,2,1,2,1,5,2]))**2)
zcat['d_uncertainty<percent>'] = 100*d_unc/np.array(zcat['distance<pc>'])


plt.plot(zcat['names'],np.array(zcat['vel<km/s>']),marker='+')

plt.vlines(np.arange(17),\
           np.array(zcat['vel<km/s>']).transpose()[0],\
           np.array(zcat['vel<km/s>']).transpose()[1],alpha=0.5)

plt.xticks(rotation='vertical')
# plt.hlines(np.array(zcat['vel<km/s>']).transpose()[1],\
#            0,np.arange(17),alpha=0.5)

plt.hlines([0,15],'california','taurus')

plt.show()


from matplotlib.patches import Circle,Ellipse
import mpl_toolkits.mplot3d.art3d as art3d


lcs,bcs,vcs = np.nanmean(zcat['lon<deg>'],axis=1),\
              np.nanmean(zcat['lat<deg>'],axis=1),\
              np.nanmean(zcat['vel<km/s>'],axis=1)


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')


for i in range(17):
    ax.scatter(lcs[i],bcs[i],vcs[i])
    ax.errorbar(lcs[i],bcs[i],vcs[i],\
                xerr=2*zcat['radius_GMC<deg>'][i],\
                yerr=zcat['radius_GMC<deg>'][i],\
                zerr=zcat['vdispersion_GMC<km/s>'][i],\
               fmt='o',label=zcat['names'][i])
    
    
ax.set(xlabel='lon [deg]',ylabel='lat [deg]',zlabel='vel [km/s]')

ax.view_init(azim=0,elev=20)
plt.legend()


nicks = ['cal','chm','cma','xbn','gem','her',\
         'lac','mad','mob','mr2','ora','orb',\
         'orL','per','pol','ros','tau']


ds = zcat['distance<pc>']


fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(projection='3d')

# ax.set_xlim(-130,-190)
ax.set_ylim(-50,50)
# ax.set_zlim(-5,20)

for i in range(17):
    r = 2*zcat['radius_GMC<deg>'][i]
    l = 2*zcat['radius_GMC<pc>'][i]
#     s = 1.5*zcat['vdispersion_GMC<km/s>'][i]
    ax.scatter(lcs[i],bcs[i],ds[i],label=zcat['names'][i])
    lb = Circle((lcs[i],bcs[i]),r,alpha=0.2)
    
    lz = Ellipse((lcs[i],ds[i]),2*r,2*l,alpha=0.2)

    bz = Ellipse((bcs[i],ds[i]),2*r,2*l,alpha=0.2)

    
    ax.add_patch(lb)
    ax.add_patch(lz)
    ax.add_patch(bz)
    
    art3d.pathpatch_2d_to_3d(lb,z=ds[i],zdir='z')
    art3d.pathpatch_2d_to_3d(lz,z=bcs[i],zdir='y')
    art3d.pathpatch_2d_to_3d(bz,z=lcs[i],zdir='x')

    
    ax.text(lcs[i],bcs[i],ds[i],nicks[i])
    
ax.set(xlabel='lon [deg]',ylabel='lat [deg]',zlabel='distance [pc]')

ax.view_init(azim=210,elev=0,roll=0)


fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(projection='3d')

# ax.set_xlim(-130,-190)
# ax.set_ylim(-30,5)
# ax.set_zlim(-5,20)

for i in range(17):
    r = 2*zcat['radius_GMC<deg>'][i]
    s = 1.5*zcat['vdispersion_GMC<km/s>'][i]
    ax.scatter(lcs[i],bcs[i],vcs[i],label=zcat['names'][i])
    lb = Circle((lcs[i],bcs[i]),r,alpha=0.2)
    
    lv = Ellipse((lcs[i],vcs[i]),2*r,s,alpha=0.2)

    bv = Ellipse((bcs[i],vcs[i]),2*r,s,alpha=0.2)

    
    ax.add_patch(lb)
    ax.add_patch(lv)
    ax.add_patch(bv)
    
    art3d.pathpatch_2d_to_3d(lb,z=vcs[i],zdir='z')
    art3d.pathpatch_2d_to_3d(lv,z=bcs[i],zdir='y')
    art3d.pathpatch_2d_to_3d(bv,z=lcs[i],zdir='x')

    
    ax.text(lcs[i],bcs[i],vcs[i],nicks[i])
    
ax.set(xlabel='lon [deg]',ylabel='lat [deg]',zlabel='vel [km/s]')

ax.view_init(azim=60,elev=0,roll=0)


california, perseus+taurus, polaris, lacerta, hercules, chamaeleon, maddalena, cma ob1, rosette,
crossbones+mr2, ora+orb+orl, gemini, mob1


xbn_mr2 = zcat[[3,9]]


per_tau = zcat[[13,16]]


orions = zcat[[10,11,12]]


seprt = zcat[[0,1,4,5,6,7,8,14,15]]


cats = fits.HDUList([fits.PrimaryHDU(),
              fits.BinTableHDU(xbn_mr2),
              fits.BinTableHDU(per_tau),
              fits.BinTableHDU(orions),
              fits.BinTableHDU(seprt)])

cats.writeto('segmented_cats.fits',overwrite=True)
