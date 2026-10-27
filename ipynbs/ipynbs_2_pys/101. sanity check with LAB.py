#!/usr/bin/env python
# coding: utf-8


next steps -> do the sanity check  on a model velocity 
-> moment2 map, compare the dispersion of HI as a function of location relative to the CO cloud


create my own velocity gradients, withh known magnitude and direction
then see if the fitting routine is accurately recovering those numbers


make sure to indicate that column density is being calculated under teh optically thin
assumption


make plots binning the data perpendicular to the rotation axis - a check for how good  plane isb


what are other possible explanations for the velocity gradient? what would the signature be 
of those physical processes (could turbulence explain this?)


look at plt.quiver again


look at figure 4 of contreras et al (https://arxiv.org/pdf/1805.01802.pdf) infall signatures 
-> want to do something similar showing
how the spectra (overlapped hi + co) to see how linewiddths change as a cuntion of location

also see if diff los through the cloud  have different peak velocities that match the velocity field
^ says somethign about whether the width of cloud spectrum is bc of turbulence or bc of a steep gradient


model a velocity gradient due to rotation, then on top of that velocity field, add in randomly
moving / turbulent eddies (and that gas would have v within the fwhm of the spectrum) -> disentangling
turbulent motion from the ordered rotational motion


reasons for velocity -> inflows, outflows, rotation, turbulence, shear


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import fun as fns
import datacube as dc

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
from astropy.table import QTable
from astropy.modeling import models,fitting


labh = fits.open('labh.fits')[0].data


zcat = QTable(fits.open('final_cat.fits')[1].data)


ora_ra_nia = 76*u.pc #recalculated from m_gmc, is the same 
ora_d_nia = 414*u.pc
ora_lc_nia = 210*u.deg
ora_bc_nia = -19*u.deg
ora_m_nia = 6.9e4*u.solMass
ora_v0h_nia,ora_v1h_nia = -7*u.km/u.s, 22*u.km/u.s


ora_dl_nia = int(((ora_ra_nia/ora_d_nia)*u.rad.to(u.deg)/(0.5)).value)


ora_dl_nia #pixels


ora_lc_px = int(((ora_lc_nia.value-360) - 180)/-0.5)


ora_bc_px = int((ora_bc_nia.value+90)/0.5)


ora_lc_px,ora_bc_px


ora_v0_px, ora_v1_px = int((ora_v0h_nia.value+458)/1.03), int((ora_v1h_nia.value+458)/1.03)


ora_v0_px,ora_v1_px


ora_lab = labh[ora_v0_px:ora_v1_px,
               ora_bc_px-ora_dl_nia:ora_bc_px+ora_dl_nia,
               ora_lc_px-ora_dl_nia:ora_lc_px+ora_dl_nia]


mom0_ora_lab = ((np.nansum(ora_lab,0)*u.K*1.03*u.km/u.s*1.82e18*(
    u.cm**-2/u.K/(u.km/u.s))*1.36*c.m_p).to(u.solMass/u.pc**2)).value


# plt.figure(figsize=(5,5))
plt.contourf(mom0_ora_lab,
            levels=np.arange(6,21),cmap='inferno',
            extent=[220,200,-29,-9])

plt.ylim(-24,-14)
plt.xlim(216,205)
plt.colorbar()


mom1_ora_lab = fns.get_mom1(ora_lab,np.arange(-7,22,1.03))


plt.contourf(mom1_ora_lab,
             extent=[220,200,-29,-9],
            levels=np.arange(5,11.5,0.5))
plt.colorbar()
plt.ylim(-26,-14)
plt.xlim(216,205)


oral,orb = np.meshgrid(np.arange(ora_bc_px-ora_dl_nia,ora_bc_px+ora_dl_nia),
                       np.arange(ora_lc_px-ora_dl_nia,ora_lc_px+ora_dl_nia))


ora_df = pd.DataFrame(np.array((
            np.ravel(oral),np.ravel(orb),
            np.ravel(mom0_ora_lab),
            np.ravel(mom1_ora_lab),
            )).transpose(),
            columns=['l','b','w','v'])


ora_ls_pc = ((ora_df['l']-np.nanmean(np.unique(ora_df['l'])))*-0.5)*u.deg.to(u.rad)*414
ora_bs_pc = ((ora_df['b']-np.nanmean(np.unique(ora_df['b'])))*0.5)*u.deg.to(u.rad)*414


init = models.Linear1D(1.,0.)
fttr = fitting.LinearLSQFitter()


ora_delvdell = fttr(init,ora_ls_pc,ora_df['v'])
ora_delvdelb = fttr(init,ora_bs_pc,ora_df['v'])


ora_Omega = np.sqrt(ora_delvdell.slope**2 + ora_delvdelb.slope**2) 


ora_theta = np.arctan(ora_delvdelb.slope/ora_delvdell.slope)*u.rad.to(u.deg)


ora_plane = 7.5 + \
    ora_delvdell.slope*(((ora_df['l']-np.nanmean(np.unique(ora_df['l'])))*-0.5)*u.deg.to(u.rad)*414)\
    + ora_delvdelb.slope*(((ora_df['b']-np.nanmean(np.unique(ora_df['b'])))*0.5)*u.deg.to(u.rad)*414)


init_P = models.Planar2D(1.,1.)
fttr_P = fitting.LMLSQFitter()


plt.imshow(np.array(ora_plane).reshape(42,42),origin='lower',
             )
plt.colorbar()
plt.xlim(41,0)


ora_delvdell.slope,ora_delvdelb.slope


ora_Pfit


ora_Pfit = fttr_P(init_P,(oral*-0.5+180)*u.deg.to(u.rad)*414,(orb*0.5-90)*u.deg.to(u.rad)*414,mom1_ora_lab)


plt.scatter((ora_df['l']*-0.5+180)*u.deg.to(u.rad)*414,ora_df['v'],s=1)
plt.plot((ora_df['l']*-0.5+180)*u.deg.to(u.rad)*414,
         -52+ora_delvdell((ora_df['l']*-0.5+180))*u.deg.to(u.rad)*414,color='red')


plt.scatter((ora_df['b']*0.5-90)*u.deg.to(u.rad)*414,ora_df['v'],s=1)
plt.plot((ora_df['b']*0.5-90)*u.deg.to(u.rad)*414,
         -76+ora_delvdelb((ora_df['b']*0.5-90))*u.deg.to(u.rad)*414,color='red')


ora_Omega/0.01


ora_theta


if ora_delvdell.slope > 0 : #first & fourth quadrants
    ora_theta = ora_theta
elif ora_delvdell.slope < 0: #second quadrant & third
    ora_theta = ora_theta + 180

if ora_theta > 180:
    ora_theta -=360


ora_theta


per_ra_nia = 42.4*u.pc #recalculated from m_gmc, is the same 
per_d_nia = 320*u.pc
per_lc_nia = 159*u.deg
per_bc_nia = -20*u.deg
per_m_nia = 1.4e4*u.solMass
per_v0h_nia,per_v1h_nia = -8*u.km/u.s, 14*u.km/u.s


per_dl_nia = int((per_ra_nia/per_d_nia)*u.rad.to(u.deg)/(0.5))


per_dl_nia


per_lc_px = int((per_lc_nia.value - 180)/-0.5)


per_bc_px = int((per_bc_nia.value + 180)/0.5)


per_lc_px,per_bc_px


per_v0_px, per_v1_px = ((np.array([per_v0h_nia.value,per_v1h_nia.value])+458)/1.03).astype('int')


per_v0_px,per_v1_px


(-8+458.6)/1.03


per_lab = labh[per_v0_px:per_v1_px,
               per_bc_px-per_dl_nia:per_bc_px+per_dl_nia,
               per_lc_px-per_dl_nia:per_lc_px+per_dl_nia]


per_lab = labh[437:459,140-15:140+15,42-15:42+15]


mom0_per_lab = ((np.nansum(per_lab,0)*u.K*1.03*u.km/u.s*1.82e18*(
    u.cm**-2/u.K/(u.km/u.s))*1.36*c.m_p).to(u.solMass/u.pc**2)).value


mom0_per_lab.shape


# plt.figure(figsize=(5,5))
plt.contourf(mom0_per_lab,
            cmap='inferno',levels=np.arange(4,16),
            extent=[166,152,-27,-12])

# plt.ylim(-2,-14)
plt.xlim(166,152)
plt.colorbar()


mom1_per_lab = fns.get_mom1(per_lab,np.arange(-8,14,1.03))


plt.contourf(mom1_per_lab,
             extent=[166,152,-27,-12],
            levels=np.arange(0,8,0.5))
plt.colorbar()
plt.xlim(166,152)


per_plane = 7.5 + \
    per_delvdell.slope*(((per_df['l']-np.nanmean(np.unique(per_df['l'])))*-0.5)*u.deg.to(u.rad)*320)\
    + per_delvdelb.slope*(((per_df['b']-np.nanmean(np.unique(per_df['b'])))*0.5)*u.deg.to(u.rad)*320)


plt.imshow(np.array(per_plane).reshape(30,30),origin='lower',
             )
plt.colorbar()


plt.scatter((per_df['l']*-0.5+180)*u.deg.to(u.rad)*320,per_df['v'],s=1)
plt.plot((per_df['l']*-0.5+180)*u.deg.to(u.rad)*320,
         -37+per_delvdell((per_df['l']*-0.5+180))*u.deg.to(u.rad)*320,color='red')


plt.scatter((per_df['b']*0.5-90)*u.deg.to(u.rad)*320,per_df['v'],s=1)
plt.plot((per_df['b']*0.5-90)*u.deg.to(u.rad)*320,
         -13+per_delvdelb((per_df['b']*0.5-90))*u.deg.to(u.rad)*320,color='red')


perl,perb = np.meshgrid(np.arange(42-15,42+15),np.arange(140-15,140+15))


per_df = pd.DataFrame(np.array((
            np.ravel(perl),np.ravel(perb),
            np.ravel(mom0_per_lab),
            np.ravel(mom1_per_lab),
            )).transpose(),
            columns=['l','b','w','v'])


per_ls_pc = ((per_df['l']-np.nanmean(np.unique(per_df['l'])))*-0.5)*u.deg.to(u.rad)*320
per_bs_pc = ((per_df['b']-np.nanmean(np.unique(per_df['b'])))*0.5)*u.deg.to(u.rad)*320


per_delvdell = fttr(init,per_ls_pc,per_df['v'])
per_delvdelb = fttr(init,per_bs_pc,per_df['v'])


per_Omega = np.sqrt(per_delvdell.slope**2 + per_delvdelb.slope**2) 
per_theta = np.arctan(per_delvdelb.slope/per_delvdell.slope)*u.rad.to(u.deg)


per_P = fttr_P(init_P,perl,perb,mom1_per_lab)


per_P


per_delvdell.slope,per_delvdelb.slope


per_Omega/0.01


per_theta


per_Omega_P = np.sqrt(per_P.slope_x**2 + per_P.slope_y**2) 
per_theta_P = np.arctan(per_P.slope_y/per_P.slope_x)*u.rad.to(u.deg)


per_Omega_P/0.01


per_theta_P


if per_delvdell.slope > 0 : #first & fourth quadrants
    per_theta = per_theta
elif per_delvdell.slope < 0: #second quadrant & third
    per_theta = per_theta + 180

if per_theta > 180:
    per_theta -=360


per_theta


zcat[13]


ora_As = np.nanmax(ora_lab,0)


mom2_ora_lab = fns.get_momN(ora_lab,np.arange(-7,22,1.03),2)


[(plt.plot(np.arange(-7,22,1.03),
           ora_lab.transpose().swapaxes(0,1).reshape(-1,29)[i]),
  plt.plot(np.arange(-7,22,1.03),
           models.Gaussian1D(np.ravel(ora_As)[i],
                            np.ravel(mom1_ora_lab)[i],
                            np.ravel(np.sqrt(mom2_ora_lab))[i])(np.arange(-7,22,1.03))),
 plt.show())
 for i in range(0,len(np.ravel(mom0_ora_lab)),10)]


per_As = np.nanmax(per_lab,0)


mom2_per_lab = fns.get_momN(per_lab,np.arange(-8,14,1.03),2)


per_lab.shape


[(plt.plot(np.arange(-8,14,1.03),
           per_lab.transpose().swapaxes(0,1).reshape(-1,22)[i]),
  plt.plot(np.arange(-8,14,1.03),
           models.Gaussian1D(np.ravel(per_As)[i],
                            np.ravel(mom1_per_lab)[i],
                            np.ravel(np.sqrt(mom2_per_lab))[i])(np.arange(-8,14,1.03))),
 plt.show())
 for i in range(0,len(np.ravel(mom0_per_lab)),10)]


zcat[[10,13]]


omegas = np.array([[0.068/0.01,0.089/0.01],[1.36,9.27],[3.71,6.7]])
thetas = np.array([[68,-64.7],[98,-44],[-102,99.4]])


[(plt.plot(['hi4pi','thummim lab','nia lab'],omegas.transpose()[i],marker='+'),
) for i in range(2)]
plt.show()


[(plt.plot(['hi4pi','thummim lab','nia lab'],thetas.transpose()[i],marker='+'),
) for i in range(2)]
plt.show()
