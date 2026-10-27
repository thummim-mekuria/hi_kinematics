#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as c
from astropy.io import fits
from astropy.table import QTable
from astropy.modeling import fitting,models
fttr = fitting.LinearLSQFitter()
init = models.Linear1D(1.,0.)

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit as fit

def line(x,m,b):
    return m*x+b

import datacube as dc
import fun as fns
import interpolate as ip


c0s = fits.open('h2_mom0.fits')
h0s = fits.open('hi_mom0.fits')

c1s = fits.open('h2_mom1.fits')
h1s = fits.open('hi_mom1.fits')

c2s = fits.open('h2_mom2.fits')
h2s = fits.open('hi_mom2.fits')


cat = QTable(fits.open('lastlastLAST_cat.fits')[1].data)


names = ['Aquila South', 'California', 'Camelopardalis', 'Chamaeleon',
       'Canis Major OB1', 'Corona Australis', 'Crossbones', 'Gemini OB1',
       'Hercules', 'Lacerta', 'Lupus', 'Maddalena', 'Monoceros OB1', 'Monoceros R2',
       r'Ophiuchus', 'Orion A', 'Orion B', r'$\lambda$-Orionis', 'Perseus', 'Polaris',
       'Rosette', 'Taurus']


lbls = ['AqS', 'Cal', 'Cam', 'Chm',
       'CMa', 'CrA', 'Xbn', 'Gem',
       'Her', 'Lac', 'Lup', 'Mad', 'MOB', 'MR2',
       r'Oph', 'OrA', 'OrB', r'$\lambda$Or', 'Per', 'Pol',
       'Ros', 'Tau']


plt.rcParams['font.size'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'


his_2r = fits.open('hi_clouds_2r.fits')
his_3r = fits.open('hi_clouds_3r.fits')


h0s_2r = [np.nansum(his_2r[i].data,axis=0)*1.3*1.82e18 for i in range(22)]
h0s_3r = [np.nansum(his_3r[i].data,axis=0)*1.3*1.82e18 for i in range(22)]

h1s_2r = [fns.get_mom1(his_2r[i].data,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100) for i in range(22)]
h1s_3r = [fns.get_mom1(his_3r[i].data,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100) for i in range(22)]


h2s_2r = [np.sqrt(fns.get_momN(his_2r[i].data,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2)) for i in range(22)]
h2s_3r = [np.sqrt(fns.get_momN(his_3r[i].data,np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,2)) for i in range(22)]


rh = (cat['radius_HI<pc>']/cat['distance<pc>'])*u.rad.to(u.deg)


plt.hist(np.ravel(h1s_3r[0].data))
plt.hist(np.ravel(h1s_2r[0].data))
plt.hist(np.ravel(h1s[0].data))


for i in [5]:

    plt.figure(dpi=200)
    plt.contourf(h1s_3r[i].data,extent=[cat['lon_HI<deg>'][i][0]+2*rh[i],cat['lon_HI<deg>'][i][1]-2*rh[i],
                                    cat['lat_HI<deg>'][i][0]-2*rh[i],cat['lat_HI<deg>'][i][1]+2*rh[i]],
                levels=np.arange(np.nanmin(h1s_3r[i].data),np.nanmax(h1s_3r[i].data),0.65/2),cmap='Spectral')
    plt.colorbar(label='v [km/s]')

    plt.vlines([cat['lon_HI<deg>'][i][0]+rh[i],cat['lon_HI<deg>'][i][1]-rh[i]],
               cat['lat_HI<deg>'][i][0]-rh[i],cat['lat_HI<deg>'][i][1]+rh[i],linestyle='--',
              colors='white')

    plt.hlines([cat['lat_HI<deg>'][i][0]-rh[i],cat['lat_HI<deg>'][i][1]+rh[i]],
              cat['lon_HI<deg>'][i][0]+rh[i],cat['lon_HI<deg>'][i][1]-rh[i],
               linestyle='--',
              colors='white')


    plt.vlines([cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1]],
               cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],linestyle='--',
              colors='white')

    plt.hlines([cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
              cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
               linestyle='--',
              colors='white')


    plt.contour(np.log10(c0s[i].data),extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                                    cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
                levels=[20.3,20.4,20.5],colors='black',linewidths=0.5)

    plt.xlim(cat['lon_HI<deg>'][i][0]+2*rh[i],cat['lon_HI<deg>'][i][1]-2*rh[i])

    plt.gca().set_aspect('equal')
#     plt.savefig('cal_mom1_r2r3',bbox_inches='tight')
    plt.show()


cat['rot_ax_HI<deg>'][1]+90,Th2[1], Th3[1]


cat['rot_ax_HI<deg>'][6]+90,Th2[6], Th3[6]


plt.hist(np.ravel(np.log10(c0s[0].data)),bins=np.arange(15,25,0.25))
plt.show()


# at some point in talking about mass, mention that we are only considering emission above 1e20 for
# both the mass and velocity gradient calculations in the case of the molecular clouds


[(plt.contourf(np.log10(h0s_3r[i].data)*(np.log10(h0s_3r[i].data)<20)*(np.log10(h0s_3r[i].data)>17)),
  plt.title(cat['name'][i]),
plt.show()) for i in range(22)]


Omh2 = np.zeros(22)
Omh3 = np.zeros(22)

Th2 = np.zeros(22)
Th3 = np.zeros(22)


hi_dvds_2r


hxfit


hi_dvds = []
hi_dvds_2r = []
hi_dvds_3r = []
co_dvds = []


for i in range(22):
   
    xc,yc = np.meshgrid(np.arange(cat['loc'][i],cat['lic'][i])*-0.125+360,
                        np.arange(cat['boc'][i],cat['bic'][i])*0.125-50)
    xh,yh = np.meshgrid(np.arange(cat['loh'][i],cat['lih'][i])*-0.125+360,
                        np.arange(cat['boh'][i],cat['bih'][i])*0.125-50)
    
    rh = (cat['bih'][i]-cat['boh'][i])/2
    
    
    xh2,yh2 = np.meshgrid(np.arange(cat['loh'][i]-int(rh),cat['lih'][i]+int(rh))*-0.125+360,
                        np.arange(cat['boh'][i]-int(rh),cat['bih'][i]+int(rh))*0.125-50)
    
    
    xh3,yh3 = np.meshgrid(np.arange(cat['loh'][i]-2*rh,cat['lih'][i]+2*rh)*-0.125+360,
                        np.arange(cat['boh'][i]-2*rh,cat['bih'][i]+2*rh)*0.125-50)
    
    
    #---------------------------------------------------------------------------------
    
    hdf = pd.DataFrame(np.array((
            np.ravel(xh),np.ravel(yh),
            np.ravel(h0s[i].data),
            np.ravel(h1s[i].data),
            np.ravel(h2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy = hdf
    
    hdf = hdf_noisy[#not adding the w>1e20 req here cos all the data is > 20
                    (hdf_noisy['s']>0) & \
                   ~(np.isinf(hdf_noisy['s'])) &\
                   ~(np.isnan(hdf_noisy['s']))] 
    #make sure that the linewidths are a valid number 
    
    
    cdf = pd.DataFrame(np.array((
            np.ravel(xc),np.ravel(yc),
            np.ravel(c0s[i].data),
            np.ravel(c1s[i].data),
            np.ravel(c2s[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    cdf_noisy = cdf
    cdf = cdf_noisy[(cdf_noisy['w']>1e20) &\
                    (cdf_noisy['s']>0) &\
                   ~(np.isinf(cdf_noisy['s'])) &\
                   ~(np.isnan(cdf_noisy['s']))]
    
    
    hdf2 = pd.DataFrame(np.array((
            np.ravel(xh2),np.ravel(yh2),
            np.ravel(h0s_2r[i].data),
            np.ravel(h1s_2r[i].data),
            np.ravel(h2s_2r[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy2 = hdf2
    
    #this and r3 do have emission < 1e20 (polaris, ora)
    
    hdf2 = hdf_noisy2[(hdf_noisy2['w']>1e20) &\
                    (hdf_noisy2['s']>0) &\
                   ~(np.isinf(hdf_noisy2['s'])) &\
                   ~(np.isnan(hdf_noisy2['s']))] 
    
    
    hdf3 = pd.DataFrame(np.array((
            np.ravel(xh3),np.ravel(yh3),
            np.ravel(h0s_3r[i].data),
            np.ravel(h1s_3r[i].data),
            np.ravel(h2s_3r[i].data)
            )).transpose(),
            columns=['l','b','w','v','s'])
    
    hdf_noisy3 = hdf3
    
    #this and r3 do have emission < 1e20 (polaris, ora)
    
    hdf3 = hdf_noisy3[(hdf_noisy3['w']>1e20) &\
                    (hdf_noisy3['s']>0) &\
                   ~(np.isinf(hdf_noisy3['s'])) &\
                   ~(np.isnan(hdf_noisy3['s']))] 
    
    
    #------------------------------------------------------------------------------------
    
    #converting to pc now makes life wayyy easier
    hxpc = hdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc = hdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    cxpc = cdf['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    cypc = cdf['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    hxpc2 = hdf2['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc2 = hdf2['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    
    
    hxpc3 = hdf3['l']*u.deg.to(u.rad)*cat['distance<pc>'][i]
    hypc3 = hdf3['b']*u.deg.to(u.rad)*cat['distance<pc>'][i]

    #------------------------------------------------------------------------------------

    
#     dvdxh = fttr(init,hxpc,hdf['v']) 
#     dvdyh = fttr(init,hypc,hdf['v'])
    
#     hi_dvds[i] = [dvdxh.slope.value,dvdyh.slope.value]
    
    
#     dvdxc = fttr(init,cxpc,cdf['v'])
#     dvdyc = fttr(init,cypc,cdf['v'])
    
#     co_dvds[i] = [dvdxc.slope.value,dvdyc.slope.value]
      
    
    hxfit = fit(line,hxpc,hdf['v'])
    hyfit = fit(line,hypc,hdf['v'])
 
    
    cxfit = fit(line,cxpc,cdf['v'])
    cyfit = fit(line,cypc,cdf['v'])
    
    
    hxfit2 = fit(line,hxpc2,hdf2['v'])
    hyfit2 = fit(line,hypc2,hdf2['v'])
    
    
    hxfit3 = fit(line,hxpc3,hdf3['v'])
    hyfit3 = fit(line,hypc3,hdf3['v'])
    
    hi_dvds.append([hxfit,hyfit])
    hi_dvds_2r.append([hxfit2,hyfit2])
    hi_dvds_3r.append([hxfit3,hyfit3])
    co_dvds.append([cxfit,cyfit])


#     print(names[i],np.around(np.array((hxfit[0][0], hxfit2[0][0], hxfit3[0][0])),5),'\n')
#     plt.plot([0,hxfit[0][0]],[0,hyfit[0][0]])
#     plt.plot([0,hxfit2[0][0]],[0,hyfit2[0][0]])
#     plt.plot([0,hxfit3[0][0]],[0,hyfit3[0][0]])
    
#     plt.plot([0,cxfit[0][0]],[0,cyfit[0][0]])
    
#     plt.show()
    
    #looking at gradient directions as the aperture gets bigger
    #definitely changes a lot for the hi
    #doesnt seem to have a pattern though
    
    
    #and now looking at the errors, for most of the clouds, the standard deviation goes down
    #uniformly with increasing radius
    
#     plt.plot(np.abs(np.around(np.array((hxfit[0][0]/np.sqrt(hxfit[1][0,0]),
#                               hxfit2[0][0]/np.sqrt(hxfit2[1][0,0]),
#                               hxfit3[0][0]/np.sqrt(hxfit3[1][0,0]))),5)))
    
    
#     plt.plot(np.abs(np.around(np.array((hyfit[0][0]/np.sqrt(hyfit[1][0,0]),
#                               hyfit2[0][0]/np.sqrt(hyfit2[1][0,0]),
#                               hyfit3[0][0]/np.sqrt(hyfit3[1][0,0]))),5)))
#     plt.title(cat['name'][i])
#     plt.show()
    
    
    #------------------------------------------------------------------------------------
    
    
#     omh2 = np.sqrt(hxfit2[0][0]**2 + hyfit2[0][0]**2)
#     omh3 = np.sqrt(hxfit3[0][0]**2 + hyfit3[0][0]**2)
    
#     Omh2[i] = omh2
#     Omh3[i] = omh3

#     th2 = np.arctan(hyfit2[0][0]/hxfit2[0][0])*u.rad.to(u.deg)
#     th3 = np.arctan(hyfit3[0][0]/hxfit3[0][0])*u.rad.to(u.deg)
    
    
#     if hxfit2[0][0] > 0 : #first & fourth quadrants
#         th2 = th2
#     elif hxfit2[0][0] > 0: #second quadrant & third
#         th2 = th2 + 180
        
#     if hxfit3[0][0] > 0 : #first & fourth quadrants
#         th3 = th3  
#     elif hxfit3[0][0] < 0: #second quadrant & third
#         th3 = th3 + 180

    
#     Th2[i] = th2
#     Th3[i] = th3
    
#     thc = np.arctan(cyfit[0][0]/cxfit[0][0])*u.rad.to(u.deg) - 90
    
#     cat['rot_ax_GMC<deg>'][i] = thc #fixing since before the 2nd & 3rd qdrnt angles were calc wrong


pd.DataFrame()


hixslopes = [hi_dvds[i][0][0][0] for i in range(22)]


mom1_fit_params = pd.DataFrame(hixslopes,columns=['hi_x_m'])


mom1_fit_params['hi_x_b'] = [hi_dvds[i][0][0][1] for i in range(22)]


mom1_fit_params['co_x_m'] = [co_dvds[i][0][0][0] for i in range(22)]
mom1_fit_params['co_x_b'] = [co_dvds[i][0][0][1] for i in range(22)]

mom1_fit_params['hi2_x_m'] = [hi_dvds_2r[i][0][0][0] for i in range(22)]
mom1_fit_params['hi2_x_b'] = [hi_dvds_2r[i][0][0][1] for i in range(22)]

mom1_fit_params['hi3_x_m'] = [hi_dvds_3r[i][0][0][0] for i in range(22)]
mom1_fit_params['hi3_x_b'] = [hi_dvds_3r[i][0][0][1] for i in range(22)]


mom1_fit_params['hi_y_m'] = [hi_dvds[i][1][0][0] for i in range(22)]
mom1_fit_params['hi_y_b'] = [hi_dvds[i][1][0][1] for i in range(22)]

mom1_fit_params['co_y_m'] = [co_dvds[i][1][0][0] for i in range(22)]
mom1_fit_params['co_y_b'] = [co_dvds[i][1][0][1] for i in range(22)]

mom1_fit_params['hi2_y_m'] = [hi_dvds_2r[i][1][0][0] for i in range(22)]
mom1_fit_params['hi2_y_b'] = [hi_dvds_2r[i][1][0][1] for i in range(22)]

mom1_fit_params['hi3_y_m'] = [hi_dvds_3r[i][1][0][0] for i in range(22)]
mom1_fit_params['hi3_y_b'] = [hi_dvds_3r[i][1][0][1] for i in range(22)]


mom1_fit_params['hi_x_msig'] = [np.sqrt(hi_dvds[i][0][1][0,0]) for i in range(22)]
mom1_fit_params['hi_x_bsig'] = [np.sqrt(hi_dvds[i][0][1][1,1]) for i in range(22)]

mom1_fit_params['co_x_msig'] = [np.sqrt(co_dvds[i][0][1][0,0]) for i in range(22)]
mom1_fit_params['co_x_bsig'] = [np.sqrt(co_dvds[i][0][1][1,1]) for i in range(22)]

mom1_fit_params['hi2_x_msig'] = [np.sqrt(hi_dvds_2r[i][0][1][0,0]) for i in range(22)]
mom1_fit_params['hi2_x_bsig'] = [np.sqrt(hi_dvds_2r[i][0][1][1,1]) for i in range(22)]

mom1_fit_params['hi3_x_msig'] = [np.sqrt(hi_dvds_3r[i][0][1][0,0]) for i in range(22)]
mom1_fit_params['hi3_x_bsig'] = [np.sqrt(hi_dvds_3r[i][0][1][1,1]) for i in range(22)]


mom1_fit_params['hi_y_msig'] = [np.sqrt(hi_dvds[i][1][1][0,0]) for i in range(22)]
mom1_fit_params['hi_y_bsig'] = [np.sqrt(hi_dvds[i][1][1][1,1]) for i in range(22)]

mom1_fit_params['co_y_msig'] = [np.sqrt(co_dvds[i][1][1][0,0]) for i in range(22)]
mom1_fit_params['co_y_bsig'] = [np.sqrt(co_dvds[i][1][1][1,1]) for i in range(22)]

mom1_fit_params['hi2_y_msig'] = [np.sqrt(hi_dvds_2r[i][1][1][0,0]) for i in range(22)]
mom1_fit_params['hi2_y_bsig'] = [np.sqrt(hi_dvds_2r[i][1][1][1,1]) for i in range(22)]

mom1_fit_params['hi3_y_msig'] = [np.sqrt(hi_dvds_3r[i][1][1][0,0]) for i in range(22)]
mom1_fit_params['hi3_y_bsig'] = [np.sqrt(hi_dvds_3r[i][1][1][1,1]) for i in range(22)]


mom1_fit_params['name'] = cat['name']


mom1_fit_params.to_csv('fit_params_mom1planefitting.csv',index=False)


delthet = cat['rot_ax_GMC<deg>']-cat['rot_ax_HI<deg>']


for i in range(22):
    if delthet[i] < -180:
        delthet[i] += 360
        
    if delthet[i] > 180:
        delthet[i] -= 360


plt.hist(delthet,bins=np.arange(-180,210,20))
plt.show()


#new delphi histogram ^


plt.scatter(cat['vel_grad_HI<km/s/pc>'],Omh2)
plt.scatter(cat['vel_grad_HI<km/s/pc>'],Omh3)
plt.plot(np.arange(0,0.3,0.05),np.arange(0,0.3,0.05))


plt.rcParams['font.size']=10


(cat['rot_ax_HI<deg>'] - (Th3-90))[10]


diff2 = cat['rot_ax_HI<deg>'] - (Th2-90)
diff3 = cat['rot_ax_HI<deg>'] - (Th3-90)

for i in range(22):
    if diff2[i] > 180:
        diff2[i] -=360
    if diff2[i] < -180:
        diff2[i] +=360
        
        
    if diff3[i] < -180:
        diff3[i] +=360
        
    if diff3[i] > 180:
        diff3[i] -=360


plt.figure(dpi=200)
plt.plot(cat['name'],diff2,marker='+',label=r'$\phi(R_H) - \phi(2R_H)$')
plt.plot(cat['name'],diff3,marker='+',label=r'$\phi(R_H) - \phi(3R_H)$')

# plt.plot(cat['name'],Th3-90,marker='+',label=r'$R = 2R_H$')
# plt.plot(cat['name'],Th2-90,marker='+',label=r'$R = 3R_H$')

# plt.plot(cat['name'],cat['rot_ax_GMC<deg>'],marker='+')

plt.xticks(rotation='vertical')
plt.legend()
plt.ylabel('$\Delta \phi\ [\degree]$')
# plt.ylim(180,-180)

plt.savefig('vel_grad_diff',bbox_inches='tight')
plt.show()


cos = fits.open('h2_clouds.fits')


plt.figure(dpi=200)
plt.vlines(np.nanmean(h1s[1].data),0,40,linestyle='--')
plt.plot(np.arange(cat['vo'][1],cat['vi'][1])*1.3-100,np.nanmean(his[1].data,axis=(1,2)),label='HI')
plt.plot(np.arange(cat['vo'][1],cat['vi'][1])*1.3-100,50*np.nanmean(cos[1].data,axis=(1,2)),label=r'$H_2$')

plt.xlabel('v [km/s]')
plt.legend()
plt.savefig('cal_spec')


plt.scatter(cat['ang_mom_HI<erg.s>'],(cat['mass_HI<Msun>']*u.solMass*
            Omh2*u.km/u.s/u.pc*(cat['radius_HI<pc>']*u.pc)**2).cgs)

plt.scatter(cat['ang_mom_HI<erg.s>'],(cat['mass_HI<Msun>']*u.solMass*
            Omh3*u.km/u.s/u.pc*(cat['radius_HI<pc>']*u.pc)**2).cgs)
plt.loglog()

plt.plot(np.arange(1e60,1e66,1e60),np.arange(1e60,1e66,1e60))
plt.show()


#hmm not sure whether to include the momentum stuff for 2r and 3r
#ie whether it says anything
#but maybe just state it as here is also another thing i did and what came out of it

#then then then
#do the patch by patch spectra / sigma vs r analysis
#and make the mom0/1/2 maps

#also, make for clouds near each other, show tg (look at lewis)
#and ghosts of co nearby that wasnt included as just one line at 21 (or 20.7 or sth)


#look at belloche (https://arxiv.org/pdf/1305.0627) for review on ang mom 
#and tobin (https://iopscience.iop.org/article/10.1088/0004-637X/748/1/16) for vel grad x flow geometry x rotation 


his = fits.open('hi_clouds.fits')


[(
plt.plot(np.nanmedian(his[i].data,axis=(1,2))),
plt.plot(np.nanmedian(his_2r[i].data,axis=(1,2))),
plt.plot(np.nanmedian(his_3r[i].data,axis=(1,2))),
plt.title(names[i]),
plt.show()
)

for i in range(22)]


cos = fits.open('h2_clouds.fits')


coshsz = fits.open('h2_hisz_clouds.fits')


mom1_vs = np.array((
[0.5,2,0.65/4],
[-15,0,1.3],
[-2.5,-1,0.65/4],
[-0.5,3,0.65/2],
[14,18,0.65/2],
[4,7,0.65/4],
[8.5,13.5,0.65/2],
[4,13,0.65],
[10.5,18.5,0.65],
[-4.5,0,0.65/2],
[-1,7,0.65],
[27,29.5,0.65/2],
[5,8,0.65/2],
[6.5,13.5,0.65],
[0,5,0.65],
[4,10,0.65],
[5,12.5,0.65],
[3,8.5,0.65],
[0,6.5,0.65],
[-4,-1,0.65/2],
[11,14,0.65/2],
[0,7,0.65]))


for i in range(22):
    n = his[i].data.shape[1]
    
    num_ptch = 6
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    
    coords = np.array((inds[:-1],inds[1:])).transpose()
    
    slcs = np.array([(i,j) for j in coords for i in coords])
    

    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=200)
    
    map_axs = np.flip(axs,axis=0).transpose().flat
    
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(his[i].data[:,bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(coshsz[i].data[:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.3,
                        linewidth=2,color='black')
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
            map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                            (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')
            
            map_axs[j].scatter(np.nanmean(fns.get_mom1(coshsz[i].data[:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,s=1,color='black')
            
            
        map_axs[j].vlines(np.nanmedian(h1s[i].data),0,1.5,color='black',
                          linestyle='--',linewidth=1.5,alpha=0.2)
        map_axs[j].vlines(np.nanmedian(c1s[i].data),0,1.5,color='black',linestyle='--',linewidth=0.4)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(his[i].data[:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
                           s=7,color='black',alpha=0.3,linewidth=0)
        
        
#         map_axs[j].set_xlim(np.array((mom1_vs[i][0],mom1_vs[i][1]))+7*np.array((-1.3,1.3,)))             
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_axis_off()
        
        
    bkg_map = fig.add_subplot(111)
    bkg_map.set_axis_off()
    bkg_map.contour(h1s[i].data,alpha=0.9,cmap='coolwarm_r',extend='both',
                     extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                             cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
                    linewidths=0.4,levels=np.arange(-1.3+mom1_vs[i][0],1.3+mom1_vs[i][1],mom1_vs[i][2]/2))
#     bkg_map.contourf(c1s[i].data,alpha=0.3,levels=10,cmap='coolwarm_r',
#                     extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                             cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
    

#     bkg_map.invert_xaxis()
    bkg_map.set_xlim(cat['lon_HI<deg>'][i])
    bkg_map.set_ylim(cat['lat_HI<deg>'][i])
#     bkg_map.set_aspect('equal')
        
    plt.subplots_adjust(wspace=-0.07,hspace=-0.07)
        
    fig.suptitle(names[i])
    plt.savefig('patch_spec_'+lbls[i])
    plt.show()


linx = np.array((3,3,3,3,3,3,3,4,4,2,4,3,3,3,4,4,3,3,3,3,3,4))
liny = np.array((4,4,2,4,3,4,4,3,3,4,5,4,3,4,4,3,4,3,4,4,3,4))


mom0_range = np.array(([20.2,21.2], #0
                    [20.6,22.],
                    [20,21.3],
                    [20.7,21.7], #3
                    [20.3,21.8],
                    [20.2,21.8],
                    [20.3,21.5], #6
                    [20.7,22.],
                    [20.3,21.6],
                    [20.3,21.6], #9
                    [20.5,21.8],
                    [20.4,21.6],
                    [20.4,22.3],#12
                    [20.3,21.8],
                    [20.3,22.1],
                    [20.6,22.3], #15
                    [20.6,22.5],
                    [20.5,22.],
                    [20.3,22], #18
                    [20.3,21.5],
                    [20.5,21.7],
                    [20.6,21.8],#21
                ))


[( 
    plt.figure(dpi=200),
        
    plt.contourf(np.log10(c0s[i].data),levels=np.linspace(mom0_range[i][0],mom0_range[i][1],15),
                                              cmap='YlGnBu_r',extend='both',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.colorbar(label=r'$N_{H_2}\ \ [cm^{-2}]$',pad=0.05,
                 shrink=0.8,ticks=np.around(np.linspace(mom0_range[i][0],mom0_range[i][1],3),1)),

    
    plt.contourf(np.log10(c0s[i].data),levels=[18,mom0_range[i][0]],colors='white',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
#     plt.contour(np.log10(c0s[i].data),levels=np.arange(20,21.5,0.3),colors='black',
#                 linewidths=[0.3,0.3,0.3,0.9],
#                extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                        cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                ),
    
    
    plt.contour(np.log10(c0s[i].data),levels=[21],colors='black',
                linewidths=[0.3],
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
               ),

 plt.xlim(cat['lon_GMC<deg>'][i]+np.array((0.1,-0.1))),
 plt.ylim(cat['lat_GMC<deg>'][i]+np.array((-0.1,0.1))),
    
 plt.xticks(np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                        cat['lon_GMC<deg>'][i][1],linx[i]),1),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lon_GMC<deg>'][i][0],
                              cat['lon_GMC<deg>'][i][1],linx[i]),1)]),
    
 plt.yticks(np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                        cat['lat_GMC<deg>'][i][1],liny[i]),1),
           labels=[str(tck)+r'$\degree$' for tck in
                  np.around(np.linspace(cat['lat_GMC<deg>'][i][0],
                              cat['lat_GMC<deg>'][i][1],liny[i]),1)]),
    
 plt.gca().set_aspect('equal'),
    
 plt.xlabel('Galactic Longitude'),
 plt.ylabel('Galactic Latitude'),
    
 plt.title(names[i]),
 plt.savefig('mom0_'+lbls[i],bbox_inches='tight'),

 plt.show())
              
              
 for i in range(22)]


bigc0 = fits.open('allsky_0th_co.fits')[0].data


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [1,18,21]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.7,22.2,0.1),
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title('California, Taurus, Perseus'),
# plt.colorbar(label=r'$N_{H_2}\ \ [cm^{-2}]$',pad=0.05,
#                  shrink=0.8,ticks=np.linspace(20.9,22,3)),


plt.xlim(180,150),
plt.ylim(-27,-5),


plt.xticks(np.linspace(180,150,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(180,150,4).astype('int')]),

plt.yticks(np.linspace(-25,-5,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-25,-5,3).astype('int')]),


plt.colorbar(shrink=0.79,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[21,21.5,22])
plt.savefig('mom0_calpertau',bbox_inches='tight')


help(plt.contourf)


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [6,13,15,16,17]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.5,22.2,0.15),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Crossbones, Monoceros R2, Orion'),


plt.xlim(-135,-175),
plt.ylim(-23,-3),


plt.xticks(np.linspace(-140,-170,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-140,-170,4).astype('int')]),

plt.yticks(np.linspace(-20,-5,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-20,-5,4).astype('int')]),


plt.colorbar(shrink=0.54,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.5,21.3,22.1])
plt.savefig('mom0_xbnmr2ors',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.8,22.21,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Rosette, Monoceros OB1, Gemini OB1'),


plt.xlim(-150,-173),
plt.ylim(-5,5),


plt.xticks(np.linspace(-150,-170,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-150,-170,3).astype('int')]),

plt.yticks(np.linspace(-5,5,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-5,5,3).astype('int')]),


plt.colorbar(shrink=0.48,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.8,21.5,22.2])
plt.savefig('mom0_rosmobgem',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.8,22.01,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Canis Major OB1, Maddalena'),


plt.xlim(-132,-147),
plt.ylim(-6,1),


plt.xticks(np.linspace(-135,-145,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-135,-145,3).astype('int')]),

plt.yticks(np.linspace(-6,0,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-6,0,3).astype('int')]),


plt.colorbar(shrink=0.5,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.8,21.4,22.])
plt.savefig('mom0_cmamad',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [10,14]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.8,22.1,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Ophiuchus, Lupus'),


plt.xlim(15,-30),
plt.ylim(2,28),


plt.xticks(np.linspace(10,-20,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(10,-20,4).astype('int')]),

plt.yticks(np.linspace(5,25,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(5,25,3).astype('int')]),


plt.colorbar(shrink=0.62,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.8,21.4,22.])
plt.savefig('mom0_ophlup',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [2,19]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.2,21.3,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0]-0.1,cat['lat_GMC<deg>'][i][1]+0.2,
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i]+np.array([-0.1,0.2]),cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Camelopardalis'),


plt.xlim(153,141),
plt.ylim(15.5,26),


plt.xticks(np.linspace(152,142,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(152,142,4).astype('int')]),

plt.yticks(np.linspace(16,25,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(16,25,4).astype('int')]),


plt.colorbar(shrink=0.9,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.2,20.7,21.2,])
plt.savefig('mom0_cam',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [2,19]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.5,21.7,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Polaris Flare'),


plt.xlim(129,116),
plt.ylim(18.5,32),


plt.xticks(np.linspace(128,116,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(128,116,4).astype('int')]),

plt.yticks(np.linspace(19,31,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(19,31,4).astype('int')]),


plt.colorbar(shrink=0.9,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.5,21,21.5,])
plt.savefig('mom0_pol',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [8]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.5,21.71,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Hercules'),


plt.xlim(50,40),
plt.ylim(6,12),


plt.xticks(np.linspace(50,40,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(50,40,3).astype('int')]),

plt.yticks(np.linspace(6,12,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(6,12,3).astype('int')]),


plt.colorbar(shrink=0.66,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.5,21.1,21.7,])
plt.savefig('mom0_her',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [3]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.5,21.71,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Chamaeleon'),


plt.xlim(-53,-68),
plt.ylim(-21,-11),


plt.yticks(np.linspace(-20,-12,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-20,-12,3).astype('int')]),

plt.xticks(np.linspace(-53,-68,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-53,-68,4).astype('int')]),


plt.colorbar(shrink=0.69,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.5,21.1,21.7,])
plt.savefig('mom0_chm',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [9]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.2,21.7,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Lacerta'),


plt.xlim(99,93),
plt.ylim(-17.5,-7),


plt.yticks(np.linspace(-17,-8,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-17,-8,4).astype('int')]),

plt.xticks(np.linspace(99,93,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(99,93,3).astype('int')]),


plt.colorbar(shrink=0.8,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.2,20.9,21.6,])
plt.savefig('mom0_lac',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.9,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [5]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.6,22.01,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Corona Australis'),


plt.xlim(5,-3.5),
plt.ylim(-25,-15),


plt.yticks(np.linspace(-25,-16,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-25,-16,4).astype('int')]),

plt.xticks(np.linspace(4,-2,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(4,-2,4).astype('int')]),


plt.colorbar(shrink=0.9,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.6,21.3,22,])
plt.savefig('mom0_cra',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(bigc0),
             levels=np.arange(19.7,22,0.1),
             extent=[360,-360,-50,50],cmap='binary_r')

for i in [0]:
    
        
    plt.contourf(np.log10(c0s[i].data),levels=np.arange(20.2,21.3,0.1),#extend='max',
                                              cmap='YlGnBu_r',
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_GMC<deg>'][i],cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    plt.hlines(cat['lat_GMC<deg>'][i],cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
           linewidths=0.6,colors='red',linestyle='--'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
plt.title(r'Aquila South'),


plt.xlim(42,22),
plt.ylim(-22,-10),


plt.yticks(np.linspace(-22,-10,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-22,-10,4).astype('int')]),

plt.xticks(np.linspace(42,22,5).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(42,22,5).astype('int')]),


plt.colorbar(shrink=0.65,label=r'$\mathrm{N_{H_2}\ [cm^{-2}]}$',ticks=[20.2,20.7,21.2,])
plt.savefig('mom0_aqs',bbox_inches='tight')


np.linspace(20.8,22.,3)


mom1_vs = np.array(([1,5,0.65/2],
[-10,5,1.3],
[-3,1,0.65/2],
[1,6,0.65],
[12,20,0.65],
[5,7,0.65/2],
[10,15,0.65],
[-3,11,1.3],
[5,15,1.3],
[-5,5,1.3],
[1,10,0.65],
[20,34,1.3],
[0,11,1.3],
[8,15,0.65],
[-1,7,0.65],
[3,15,1.3],
[6,13,0.65],
[-7,16,2.6],
[-1,10,1.3],
[-6,1,0.65],
[8,20,1.3],
[4,11,0.65]))


n=25
mom0_range_hi = [
    np.linspace(20.6,20.8,n),#0
    np.linspace(21,21.7,n),
    np.linspace(20.75,20.95,n),#2
    np.linspace(20.7,21.,n),
    np.linspace(21.1,21.5,n),#4
    np.linspace(20.6,21.,n),
    np.linspace(21.2,21.5,n), #6
    np.linspace(21.5,21.8,n),
    np.linspace(21,21.5,n),#8
    np.linspace(20.9,21.3,n),
    np.linspace(20.6,21.4,n),#10
    np.linspace(20.7,21.5,n),
    np.linspace(20.8,21.4,n),#12
    np.linspace(21.05,21.5,n),
    np.linspace(20.9,21.3,n),#14
    np.linspace(20.7,21.4,n),
    np.linspace(21,21.4,n),#16
    np.linspace(21,21.5,n),
    np.linspace(20.7,21.2,n),#18
    np.linspace(20.3,20.8,n),
    np.linspace(21.3,21.6,n),#20
    np.linspace(20.6,21.4,n),
]


bigh0 = fits.open('allsky_0th_hi.fits')[0].data


cat['vi'][i]-cat['vo'][i]


for i in [1,18,21]:
    plt.hist(np.ravel(np.log10(0.8*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=250),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.6,21.6,18),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.7,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.6,21.1,21.6],pad=0.1)


plt.contour(np.log10(bigc0),
             levels=np.arange(20.5,22,0.25),linewidths=0.3,alpha=0.2,
             extent=[360,-360,-50,50],cmap='binary_r')


for i in [1,18,21]:
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.7,22.2,0.2),alpha=0.9,
                                              cmap='magma',linewidths=0.5,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title('California, Taurus, Perseus'),
# plt.colorbar(label=r'$N_{H_2}\ \ [cm^{-2}]$',pad=0.05,
#                  shrink=0.8,ticks=np.linspace(20.9,22,3)),


plt.xlim(190,140),
plt.ylim(-30,0),


plt.xticks(np.linspace(190,140,6).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(190,140,6).astype('int')]),

plt.yticks(np.linspace(-25,-5,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-25,-5,3).astype('int')]),


plt.savefig('mom0_calpertau',bbox_inches='tight')


for i in [4,11]:
    plt.hist(np.ravel(np.log10(0.3*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.3*bigh0),
             levels=np.linspace(21,21.4,18),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')

plt.colorbar(shrink=0.6,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[21,21.2,21.4])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,21.2,0.25),linewidths=0.3,alpha=0.2,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.8,21.8,0.15),extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Canis Major OB1, Maddalena'),


plt.xlim(-131,-150),
plt.ylim(-7,3),


plt.xticks(np.linspace(-130,-150,5).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-130,-150,5).astype('int')]),

plt.yticks(np.linspace(-6,3,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-6,3,4).astype('int')]),


plt.savefig('mom0_cmamad',bbox_inches='tight')


for i in [7,12,20]:
    plt.hist(np.ravel(np.log10(0.5*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.5*bigh0),
             levels=np.linspace(20.9,21.6,18),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.53,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.9,21.2,21.5])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')


for i in [7,12,20]:
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.8,22.21,0.2),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Rosette, Monoceros OB1, Gemini OB1'),


plt.xlim(-148,-175),
plt.ylim(-7,7),


plt.xticks(np.linspace(-150,-174,5).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-150,-174,5).astype('int')]),

plt.yticks(np.linspace(-5,5,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-5,5,3).astype('int')]),


plt.savefig('mom0_rosmobgem',bbox_inches='tight')


for i in [6,13,15,16,17]:
    plt.hist(np.ravel(np.log10(0.8*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.4,21.6,18),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.63,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.4,21.,21.6,])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.8,22.5,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Crossbones, Monoceros R2, Orion'),


plt.xlim(-135,-175),
plt.ylim(-28,-5),


plt.xticks(np.linspace(-140,-170,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-140,-170,4).astype('int')]),

plt.yticks(np.linspace(-25,-5,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-25,-5,3).astype('int')]),


plt.savefig('mom0_xbnmr2ors',bbox_inches='tight')


plt.figure(dpi=200),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.5,21.7,18),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.63,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.5,21.1,21.7,])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.8,22.1,0.2),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Ophiuchus, Lupus'),


plt.xlim(15,-30),
plt.ylim(2,28),


plt.xticks(np.linspace(10,-20,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(10,-20,4).astype('int')]),

plt.yticks(np.linspace(5,25,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(5,25,3).astype('int')]),


plt.savefig('mom0_ophlup',bbox_inches='tight')


for i in [2]:
    plt.hist(np.ravel(np.log10(0.8*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.7,21.1,23),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.63,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.7,20.9,21.1,])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.4,21.3,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Camelopardalis'),


plt.xlim(151,143),
plt.ylim(15,20),


plt.xticks(np.linspace(151,143,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(151,143,3).astype('int')]),

plt.yticks(np.linspace(16,20,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(16,20,3).astype('int')]),


plt.savefig('mom0_cam',bbox_inches='tight')


for i in [8]:
    plt.hist(np.ravel(np.log10(0.8*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.9,21.6,21),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.87,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.9,21.2,21.5])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.5,21.71,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Hercules'),


plt.xlim(50,40),
plt.ylim(4,14),


plt.xticks(np.linspace(50,40,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(50,40,3).astype('int')]),

plt.yticks(np.linspace(5,13,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(5,13,3).astype('int')]),


plt.savefig('mom0_her',bbox_inches='tight')


for i in [9]:
    plt.hist(np.ravel(np.log10(0.8*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.8,21.4,21),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.87,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.8,21.1,21.4])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.2,21.7,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Lacerta'),


plt.xlim(99,93),
plt.ylim(-17.5,-7),


plt.yticks(np.linspace(-17,-8,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-17,-8,4).astype('int')]),

plt.xticks(np.linspace(99,93,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(99,93,3).astype('int')]),


plt.savefig('mom0_lac',bbox_inches='tight')


for i in [3]:
    plt.hist(np.ravel(np.log10(0.8*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.8*bigh0),
             levels=np.linspace(20.6,21.2,21),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.87,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.6,20.9,21.2])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.5,21.71,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Chamaeleon'),


plt.xlim(-53,-68),
plt.ylim(-23,-9),


plt.yticks(np.linspace(-23,-9,3).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-23,-9,3).astype('int')]),

plt.xticks(np.linspace(-53,-68,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-53,-68,4).astype('int')]),


plt.savefig('mom0_chm',bbox_inches='tight')


for i in [0]:
    plt.hist(np.ravel(np.log10(0.65*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.65*bigh0),
             levels=np.linspace(20.4,21.2,21),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.7,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.4,20.8,21.2])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.2,21.3,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Aquila South'),


plt.xlim(42,27),
plt.ylim(-22,-10),


plt.yticks(np.linspace(-22,-10,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-22,-10,4).astype('int')]),

plt.xticks(np.linspace(42,27,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(42,27,4).astype('int')]),

plt.savefig('mom0_aqs',bbox_inches='tight')


for i in [19]:
    plt.hist(np.ravel(np.log10(0.6*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.6*bigh0),
             levels=np.linspace(20.4,21.,21),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.7,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.4,20.7,21.])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.25),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.6,21.7,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Polaris Flare'),


plt.xlim(129,116),
plt.ylim(18.5,32),


plt.xticks(np.linspace(128,116,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(128,116,4).astype('int')]),

plt.yticks(np.linspace(19,31,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(19,31,4).astype('int')]),

plt.savefig('mom0_pol',bbox_inches='tight')


for i in [5]:
    plt.hist(np.ravel(np.log10(0.9*bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5)
#     plt.hist(np.log10((cat['vi'][i]-cat['vo'][i])*1.3/155*np.ravel(
#        bigh0[cat['boh'][i]:cat['bih'][i],cat['loh'][i]:cat['lih'][i]])),alpha=0.5),
    plt.hist(np.ravel(np.log10(h0s[i].data)),alpha=0.5)
    plt.title(names[i])
    plt.show()


plt.figure(dpi=200),
plt.contourf(np.log10(0.9*bigh0),
             levels=np.linspace(20.6,21.,21),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')


plt.colorbar(shrink=0.83,label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$',ticks=[20.6,20.8,21.])


plt.contour(np.log10(bigc0),
             levels=np.arange(19.9,22,0.2),alpha=0.2,linewidths=0.3,
             extent=[360,-360,-50,50],cmap='binary_r')

for i in range(22):
    
        
    plt.contour(np.log10(c0s[i].data),levels=np.arange(20.6,21.7,0.15),#extend='max',
                                              cmap='magma',linewidths=0.5,alpha=0.9,
               extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                       cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
    
    plt.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    plt.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
           linewidths=0.7,colors='white',linestyle=':'),
    
plt.gca().set_aspect('equal'),
    
plt.xlabel('Galactic Longitude'),
plt.ylabel('Galactic Latitude'),
                  
# plt.title(r'Corona Australis'),


plt.xlim(5,-3.5),
plt.ylim(-25,-15),


plt.yticks(np.linspace(-25,-16,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-25,-16,4).astype('int')]),

plt.xticks(np.linspace(4,-2,4).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(4,-2,4).astype('int')]),

plt.savefig('mom0_cra',bbox_inches='tight')


from matplotlib.patches import Circle


plt.rcParams['font.size'] = 15


fig,ax = plt.subplots(1,1,dpi=200,figsize=(25,7))

h0im = ax.contourf(np.log10(bigh0),
             levels=np.linspace(20.,22.,18),extend='both',
             extent=[360,-360,-50,50],cmap='YlGnBu_r')

fig.colorbar(h0im,ax=ax,ticks=[20,21,22],label=r'$\mathrm{log\ (N_{HI}\ [cm^{-2}])}$')

ax.contourf(np.log10(bigc0),
             levels=np.linspace(20.4,22.6,18),alpha=0.8,
             extent=[360,-360,-50,50],cmap='magma')

# plt.gca().set_aspect('equal'),

ax.set_xlim(180,-180)

for i in range(22):
    ax.add_patch(Circle((np.nanmean(cat['lon_GMC<deg>'][i]),np.nanmean(cat['lat_GMC<deg>'][i])),
                       (cat['radius_HI<pc>'][i]/cat['distance<pc>'][i])*u.rad.to(u.deg),fill=False,
                       edgecolor='white',linestyle='--',linewidth=0.7))
    
    
ax.set_xticks(np.linspace(180,-180,7).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(180,-180,7).astype('int')])


ax.set_yticks(np.linspace(-40,40,5).astype('int'),
           
           labels=[str(tck)+r'$\degree$' for tck in
                  np.linspace(-40,40,5).astype('int')])


ax.set_xlabel('Galactic Longitude')
ax.set_ylabel('Galactic Latitude')

plt.savefig('catalog_0th',bbox_inches='tight')


for i in range(22):

    plt.figure(dpi=200)
    plt.contourf(h1s_3r[i].data,extent=[cat['lon_HI<deg>'][i][0]+2*rh[i],
                                        cat['lon_HI<deg>'][i][1]-2*rh[i],
                                    cat['lat_HI<deg>'][i][0]-2*rh[i],
                                        cat['lat_HI<deg>'][i][1]+2*rh[i]],
                levels=np.arange(np.nanmin(h1s_2r[i].data),
                                 np.nanmax(h1s_2r[i].data),0.65/2),
                 cmap='Spectral',extend='both')
    plt.colorbar(label='v [km/s]')

    plt.vlines([cat['lon_HI<deg>'][i][0]+rh[i],cat['lon_HI<deg>'][i][1]-rh[i]],
               cat['lat_HI<deg>'][i][0]-rh[i],cat['lat_HI<deg>'][i][1]+rh[i],linestyle='--',
              colors='white')

    plt.hlines([cat['lat_HI<deg>'][i][0]-rh[i],cat['lat_HI<deg>'][i][1]+rh[i]],
              cat['lon_HI<deg>'][i][0]+rh[i],cat['lon_HI<deg>'][i][1]-rh[i],
               linestyle='--',
              colors='white')


    plt.vlines([cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1]],
               cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],linestyle='--',
              colors='white')

    plt.hlines([cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
              cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
               linestyle='--',
              colors='white')


    plt.contour(np.log10(c0s[i].data),extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                                    cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
                levels=[20.3,20.4,20.5],colors='black',linewidths=0.5)

    plt.xlim(cat['lon_HI<deg>'][i][0]+2*rh[i],cat['lon_HI<deg>'][i][1]-2*rh[i])

    plt.gca().set_aspect('equal')
#     plt.savefig('cal_mom1_r2r3',bbox_inches='tight')
    plt.show()


from photutils.centroids import centroid_com as com


[
for i in range(22)]


# [(cat['name'][i],
#  com(np.log10(h0s[i].data)).astype('int')-(np.array(h0s[i].data.shape)*0.5))
    
    
#     for i in range(22)
# ]
# #these are at most 1.5 pix off from each other, so its okay to just use the 
# #geometric mean as the center of mass


def calc_jsum(mom0_maps,mom1_maps,cl_typ,cat):
    n = len(mom1_maps)
    
    if n == 1:
        dist = [cat['distance<pc>']]
        if cl_typ=='CO':
            lo = [cat['loc']]
            li = [cat['lic']]
            bo = [cat['boc']]
            bi = [cat['bic']]
            mu = 1.36*2

        elif cl_typ=='21cm':
            lo = [cat['loh']]
            li = [cat['lih']]
            bo = [cat['boh']]
            bi = [cat['bih']]
            mu = 1.36
    
    else:
        dist = cat['distance<pc>']
        if cl_typ=='CO':
            lo = cat['loc']
            li = cat['lic']
            bo = cat['boc']
            bi = cat['bic']
            mu = 1.36*2

        elif cl_typ=='21cm':
            lo = cat['loh']
            li = cat['lih']
            bo = cat['boh']
            bi = cat['bih']
            mu = 1.36
        
        
    com_coords = [(np.array(mom1_maps[i].data.shape)*0.5).astype('int') 
                  for i in range(n)]
    
#     print(com_coords)
    
    coords = [np.array((
                 np.meshgrid(
                 np.arange(lo[i],li[i])*-0.125+360,
                 np.arange(bo[i],bi[i])*0.125-50)
                 ))*u.deg.to(u.rad)*dist[i]
              
                    for i in range(n)]
    
    normed_h1s = [mom1_maps[i].data-\
                  mom1_maps[i].data[com_coords[i][0],com_coords[i][1]]
                     for i in range(n)]
    
    for i in range(22): #if the geometric center is somewhere with no velocity data
        if np.all(np.isnan(normed_h1s[i])):
            normed_h1s[i] = mom1_maps[i].data-\
                      np.nanmedian(mom1_maps[i].data)
                         
    
    mass_maps = [((mom0_maps[i].data*(mu*c.m_p)/u.cm**2\
                  ).to(u.solMass/u.pc**2))*\
                 (0.125*u.deg.to(u.rad)*dist[i]*u.pc)**2
                 for i in range(n)]
    
    r_vector = [np.array((
                 coords[i][0] - np.nanmean(coords[i][0]),
                 coords[i][1] - np.nanmean(coords[i][1]),
                 np.zeros(coords[i][0].shape)
                 ))   for i in range(n)]
    
#     print(r_vector[0].reshape(3,-1).transpose())
    
    v_vector = [np.array((
                np.zeros(normed_h1s[i].shape),
                np.zeros(normed_h1s[i].shape),
                normed_h1s[i]
                ))  for i in range(n)]
    
#     print(v_vector[0].reshape(3,-1).transpose())


    jsum_kmpc_s = [np.nansum(np.cross(
                 r_vector[i].reshape(3,-1).transpose()*u.pc,
                 v_vector[i].reshape(3,-1).transpose()*u.km/u.s)*\
                 mass_maps[i].ravel()[:,np.newaxis],
                 axis=0)/np.nansum(mass_maps[i])
                   for i in range(n)]
    
    
    return jsum_kmpc_s,mass_maps,r_vector,v_vector


#hmm okay something that might turn out to be relevant is, i am calculating the 
#summed up j assuming by calculating relative to the center of mass
#and also relative to a point, rather than a rotational axis
#and in most cases, the calculated angular momentum when summed up is smaller 
#compared to the ones calculated with the velocity gradient

#can i say, if rotation was a good description of the motion / dynamics of the objects
#then the angular momentum from j=omega r2 should be colse to summed up
#but maybe a more robust comparison requires assuming the same rotation axis?

#and for radius, just say here is the ways in which the radius we calculate 
#with the geometric meaning is different than with the area


sumjs_co = np.linalg.norm(np.array(
    [calc_jsum(c0s,c1s,'CO',cat)[0][i].to(u.cm**2/u.s).value 
          for i in range(22)]),axis=1)


sumjs_co


sumjs_hi = np.linalg.norm(np.array(
    [calc_jsum(h0s,h1s,'21cm',cat)[0][i].to(u.cm**2/u.s).value 
          for i in range(22)]),axis=1)


cat['sum_j_HI<cm2.s>'] = sumjs_hi
cat['sum_j_GMC<cm2.s>'] = sumjs_co


jspec_gmc = (cat['ang_mom_GMC<erg.s>']*u.erg*u.s/\
 (cat['mass_GMC<Msun>']*u.solMass)).cgs.value

jspec_hi = (cat['ang_mom_HI<erg.s>']*u.erg*u.s/\
 (cat['mass_HI<Msun>']*u.solMass)).cgs.value


areas_npix = np.array(([len(c0s[i].data[c0s[i].data>1e20])*\
                 (0.125*u.deg.to(u.rad)*cat['distance<pc>'][i])**2 
                        for i in range(22)]))


cat['area_GMC<pc2>'] = areas_npix


plt.plot(cat['name'],
          np.abs(np.sqrt(areas_npix/np.pi)-
         cat['radius_GMC<pc>'])/cat['radius_GMC<pc>'],
         marker='+')


plt.xticks(rotation='vertical')
# plt.hlines(0.8,'aquila_s','taurus')
plt.show()


plt.scatter(cat['name'],
         np.abs(sumjs_co-(0.4*jspec_gmc))/sumjs_co,marker='+')

plt.xticks(rotation='vertical')
# plt.hlines(1,'aquila_s','taurus')
plt.show()


plt.plot(sumjs_hi/(0.4*jspec_hi))
plt.plot(np.arange(22),np.full(22,1))


np.array((sumjs_hi/jspec_hi)).transpose()


fttr(init,np.log10(jspec_gmc),np.log10(jspec_hi))


fttr(init,np.log10(sumjs_co),np.log10(sumjs_hi))


plt.scatter(jspec_gmc,jspec_hi)
plt.scatter(sumjs_co,sumjs_hi)

plt.loglog()
plt.show()


plt.scatter(cat['radius_GMC<pc>'],jspec_gmc)
plt.scatter(cat['radius_GMC<pc>'],sumjs_co)


plt.loglog()
plt.show()


plt.scatter(cat['radius_HI<pc>'],jspec_hi)
plt.scatter(cat['radius_HI<pc>'],sumjs_hi)


plt.loglog()
plt.show()


plt.plot(1/3*jspec_hi/sumjs_hi)
# plt.plot(1/3*jspec_gmc/sumjs_co)


plt.scatter(cat['vel_grad_HI<km/s/pc>'],sumjs_hi)
plt.scatter(cat['vel_grad_HI<km/s/pc>'],jspec_hi)

# plt.loglog()
plt.yscale('log')
plt.show()


plt.scatter(cat['vel_grad_GMC<km/s/pc>'],sumjs_co)
plt.scatter(cat['vel_grad_GMC<km/s/pc>'],jspec_gmc)

# plt.loglog()
plt.yscale('log')
plt.show()


plt.contourf(h1s[1].data)
plt.colorbar()


[(
plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
         np.nanmean(his[i].data,axis=(1,2))),
plt.plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
         100*np.nanmean(cos[i].data,axis=(1,2))),

plt.show()

)for i in range(22)]


i = 1
plt.plot(np.nanmean(his[i].data[17:,:,:],axis=(1,2))),
plt.plot(100*np.nanmean(cos[i].data[17:,:,:],axis=(1,2))),

plt.show()


calh1_short = fns.get_mom1(his[i].data[17:,:,:],
             np.arange(cat['vo'][1]+17,cat['vi'][1])*1.3-100)


plt.contourf(h1s[1].data,levels=np.arange(-8,5,0.65),cmap='RdGy')
plt.colorbar()


plt.contourf(calh1_short,levels=20,cmap='RdGy')
plt.colorbar()


for i in [1]:
    n = his[i].data.shape[1]
    
    num_ptch = 6
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    
    coords = np.array((inds[:-1],inds[1:])).transpose()
    
    slcs = np.array([(i,j) for j in coords for i in coords])
    

    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=200)
    
    map_axs = np.flip(axs,axis=0).transpose().flat
    
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(his[i].data[17:,bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(coshsz[i].data[17:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(cat['vo'][i]+17,cat['vi'][i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.3,
                        linewidth=2,color='black')
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
            map_axs[j].plot(np.arange(cat['vo'][i]+17,cat['vi'][i])*1.3-100,
                            (sp2-np.nanmin(sp2))/np.nanmax(sp2),
                            linewidth=0.4,color='black')
            
            map_axs[j].scatter(np.nanmean(fns.get_mom1(coshsz[i].data[17:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i]+17,cat['vi'][i])*1.3-100)),0,
                               s=1,color='black')
            
            
        map_axs[j].vlines(np.nanmedian(calh1_short),0,1.5,color='black',
                          linestyle='--',linewidth=1.5,alpha=0.2)
        map_axs[j].vlines(np.nanmedian(h1s[i].data),0,1.5,color='green',
                          linestyle='--',linewidth=1.5,alpha=0.2)
        map_axs[j].vlines(np.nanmedian(c1s[i].data),0,1.5,color='black',linestyle='--',linewidth=0.4)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(his[i].data[17:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i]+17,cat['vi'][i])*1.3-100)),0,
                           s=7,color='black',alpha=0.3,linewidth=0)
        
        
#         map_axs[j].set_xlim(np.array((mom1_vs[i][0],mom1_vs[i][1]))+7*np.array((-1.3,1.3,)))             
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_axis_off()
        
        
    bkg_map = fig.add_subplot(111)
    bkg_map.set_axis_off()
    bkg_map.contour(calh1_short,alpha=0.9,cmap='coolwarm_r',extend='both',
                     extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                             cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
                    linewidths=0.4,levels=np.arange(-8,5,0.65))
#     bkg_map.contourf(c1s[i].data,alpha=0.3,levels=10,cmap='coolwarm_r',
#                     extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                             cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
    

#     bkg_map.invert_xaxis()
    bkg_map.set_xlim(cat['lon_HI<deg>'][i])
    bkg_map.set_ylim(cat['lat_HI<deg>'][i])
#     bkg_map.set_aspect('equal')
        
    plt.subplots_adjust(wspace=-0.07,hspace=-0.07)
        
    fig.suptitle(names[i])
#     plt.savefig('patch_spec_'+lbls[i])
    plt.show()


i = 8
plt.plot(np.nanmean(his[i].data[:36,:,:],axis=(1,2))),
plt.plot(100*np.nanmean(cos[i].data[:36,:,:],axis=(1,2))),

plt.show()


herh1_short = fns.get_mom1(his[i].data[:36,:,:],
             np.arange(cat['vo'][i],cat['vo'][i]+36)*1.3-100)


plt.contourf(herh1_short,cmap='RdGy',levels=20)
plt.colorbar()


plt.contourf(h1s[8].data,cmap='RdGy',levels=20)#,levels=np.arange(0,12,0.65))
plt.colorbar()


for i in range(22):
    n = his[i].data.shape[1]
    
    num_ptch = 8
    
    inds = np.linspace(0,n,num_ptch+1).astype('int')
    
    coords = np.array((inds[:-1],inds[1:])).transpose()
    
    slcs = np.array([(i,j) for j in coords for i in coords])
    

    fig,axs = plt.subplots(num_ptch,num_ptch,
                           sharex=True,sharey=True,
                          dpi=200)
    
    map_axs = np.flip(axs,axis=0).transpose().flat
    
    
    for j in range(len(slcs)):
        bo,bi = slcs[j][0]
        lo,li = slcs[j][1]
        sp = np.nanmean(his[i].data[:,bo:bi,lo:li],axis=(1,2))
        sp2 = np.nanmean(coshsz[i].data[:,bo:bi,lo:li],axis=(1,2))
        
        map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                        (sp-np.nanmin(sp))/np.nanmax(sp),alpha=0.3,
                        linewidth=2,color='black')
        #the minus min normalizes it to start at 0 and the divide by max normalizes the top to be at 1
        if np.nansum(np.abs(sp2))>3e-1: #ie theres data here
            map_axs[j].plot(np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                            (sp2-np.nanmin(sp2))/np.nanmax(sp2),linewidth=0.4,color='black')
            
            map_axs[j].scatter(np.nanmean(fns.get_mom1(coshsz[i].data[:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,s=1,color='black')
            
            
        map_axs[j].vlines(np.nanmedian(herh1_short),0,1.5,color='black',
                          linestyle='--',linewidth=1.5,alpha=0.2)
        map_axs[j].vlines(np.nanmedian(c1s[i].data),0,1.5,color='black',linestyle='--',linewidth=0.4)
        
        map_axs[j].scatter(np.nanmean(fns.get_mom1(his[i].data[:,bo:bi,lo:li],
                                   np.arange(cat['vo'][i],cat['vi'][i])*1.3-100)),0,
                           s=7,color='black',alpha=0.3,linewidth=0)
        
#         map_axs[j].annotate(str(np.around(np.nanmean(
#                             np.sqrt(fns.get_momN(sp, 
#                             np.arange(cat['vo'][i],
#                             cat['vi'][i])*1.3-100,2))),2)),
#                             (cat['vo'][i]*1.3-100,0.2),fontsize=5)
        
#         map_axs[j].set_xlim(np.array((mom1_vs[i][0],mom1_vs[i][1]))+7*np.array((-1.3,1.3,)))             
        
        map_axs[j].set_xticks([]),
        map_axs[j].set_yticks([])
        map_axs[j].set_axis_off()
        
        
    bkg_map = fig.add_subplot(111)
    bkg_map.set_axis_off()
    bkg_map.contour(herh1_short,cmap='coolwarm_r',extend='both',
                     extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                             cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
                    linewidths=0.4,levels=np.arange(5,11,0.65/2))
#     bkg_map.contourf(c1s[i].data,alpha=0.3,levels=10,cmap='coolwarm_r',
#                     extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                             cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]])
    
    
#     bkg_map.contour(np.log10(c0s[i].data),levels=[20.5],colors='purple',
#                     linestyles='--',zorder=-30,linewdiths=0.7,alpha=0.7,
#                     extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                             cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                    )
    

#     bkg_map.invert_xaxis()
    bkg_map.set_xlim(cat['lon_HI<deg>'][i])
    bkg_map.set_ylim(cat['lat_HI<deg>'][i])
#     bkg_map.set_aspect('equal')
        
    plt.subplots_adjust(wspace=-0.07,hspace=-0.07)
        
    fig.suptitle(names[i])
#     plt.savefig('crowded_patch_spec_'+lbls[i])
    plt.savefig('patch_spec_'+lbls[i])

    plt.show()


#also seems that the linewidth is smaller -> in this case, that means there is just one emission peak, rather 
#than two or three blended together ->in the places with co emission


c2s_hsz = [np.sqrt(fns.get_momN(coshsz[i].data,
                               np.arange(cat['vo'][i],cat['vi'][i])*1.3-100,
                               2))
    
    for i in range(22)]


c0s_hsz = [np.nansum(coshsz[i].data,axis=0)*1.3*1.97e20
    
    for i in range(22)]


plt.contourf(h2s[i].data,levels=np.arange(0.5,15))
plt.colorbar()


plt.contourf((c0s_hsz[i]>1e20)*h2s[i].data,levels=np.arange(0.5,15))
plt.colorbar()


c0s_hsz[i].size/40


i


[(
plt.hist(h1s[i].data[c0s_hsz[i]>1e20],alpha=0.5,bins=25,density=True),
plt.hist(h1s[i].data[c0s_hsz[i]<1e20],alpha=0.5,bins=25,density=True),
plt.hist(np.ravel(h1s[i].data),alpha=0.5,bins=25,density=True,histtype='step'),
    plt.title(names[i]),
#     plt.xlim(3,),
plt.show()
) for i in range(22)]


[(
plt.hist(h2s[i].data[c0s_hsz[i]>1e20],alpha=0.5,bins=25,density=True),
plt.hist(h2s[i].data[c0s_hsz[i]<1e20],alpha=0.5,bins=25,density=True),
plt.hist(np.ravel(h2s[i].data),alpha=0.5,bins=25,density=True,histtype='step'),
    plt.title(names[i]),
#     plt.xlim(3,),
plt.show()
) for i in range(22)]


#okay moving forward, formalize the results from the patch by patch spectra
#maybe do the sigma vs v thing
#formalize results from two diff methods of calculating j


bigc1 = fits.open('allsky_1st_co.fits')[0].data
bigc1 = np.where(bigc0<1e20,np.nan,bigc1)
fits.PrimaryHDU(bigc1).writeto('allsky_1st_co.fits',overwrite=True)


for i in [1]:
    plt.contourf(bigc1,
                 levels=np.linspace(np.nanmean(c1s[i].data)-np.nanstd(c1s[i].data),
                                    np.nanmean(c1s[i].data)+np.nanstd(c1s[i].data),
                                  20),
                 extent=[360,-360,-50,50],extend='both',alpha=0.5,
                cmap='PiYG')
    plt.colorbar()
    plt.xlim(cat['lon_GMC<deg>'][i])
    plt.ylim(cat['lat_GMC<deg>'][i])
    plt.show()


for i in range(22):
    fig, axs = plt.subplots(1,3,dpi=200,figsize=(15,5))
    
    h1im = axs[1].contourf(h1s_3r[i].data,
         levels=np.linspace(np.nanmean(h1s_2r[i].data)-3*np.nanstd(h1s_2r[i].data),
                               np.nanmean(h1s_2r[i].data)+3*np.nanstd(h1s_2r[i].data),
                                  40),extend='both',
        extent=[cat['lon_HI<deg>'][i][0]+2*rh[i],cat['lon_HI<deg>'][i][1]-2*rh[i],
                cat['lat_HI<deg>'][i][0]-2*rh[i],cat['lat_HI<deg>'][i][1]+2*rh[i]],
         cmap='RdGy')
    
    
    axs[1].contourf(h1s[i].data,
         levels=np.linspace(np.nanmean(h1s[i].data)-3*np.nanstd(h1s[i].data),
                               np.nanmean(h1s[i].data)+3*np.nanstd(h1s[i].data),
                                  40),extend='both',
        extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
                cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]],
         cmap='RdGy')


#     fig.colorbar(h1im,ax=axs[1],shrink=0.87,label=r'$\mathrm{V_{LOS}\ [km\ s^{-1}]}$')

#     axs[1].contour(np.log10(bigc0),
#              levels=np.arange(19.9,22.,0.2),
#                    alpha=0.2,linewidths=0.3,
#              extent=[360,-360,-50,50],cmap='binary_r')
    
    
    axs[1].contour(np.log10(c0s[i].data),levels=[20.5,21],
                   alpha=0.8,linewidths=0.6,
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],colors='black')
    
    axs[0].contour(np.log10(c0s[i].data),levels=[20.5,21],
                   alpha=0.8,linewidths=0.6,
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],colors='black',
                  zorder=20)
    
    
#     axs[1].contourf(c1s[i].data,
#              levels=np.linspace(np.nanmean(c1s[i].data)-2.2*np.nanstd(c1s[i].data),
#                                np.nanmean(c1s[i].data)+2.2*np.nanstd(c1s[i].data),
#                                   15),
#                    linewidths=0.4,alpha=0.3,
#                    extend='both',
#              extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                 cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
#                     cmap='RdGy')
    
    
    axs[0].contourf(bigc1,
             levels=np.linspace(np.nanmean(c1s[i].data)-2.2*np.nanstd(c1s[i].data),
                               np.nanmean(c1s[i].data)+2.2*np.nanstd(c1s[i].data),
                                  10),
                   alpha=0.2,antialiased=True,
                   extend='both',
             extent=[360,-360,-50,50],cmap='PiYG')
    
    axs[0].contour(c1s[i].data,
             levels=np.linspace(np.nanmean(c1s[i].data)-2.2*np.nanstd(c1s[i].data),
                               np.nanmean(c1s[i].data)+2.2*np.nanstd(c1s[i].data),
                                  40),
                   linewidths=0.5,
                   extend='both',
             extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
                cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]],
                    cmap='PiYG')
    
    
#     c1im = axs[0].contour(c1s[i].data,levels=20,#extend='max',
#                cmap='Spectral',linewidths=0.5,alpha=0.9,
#                extent=[cat['lon_GMC<deg>'][i][0],cat['lon_GMC<deg>'][i][1],
#                        cat['lat_GMC<deg>'][i][0],cat['lat_GMC<deg>'][i][1]]),
    
#     fig.colorbar(c1im,ax=axs[0],label='vh2')
    
    
    for ax in axs.flat[:2]:
        ax.vlines(cat['lon_HI<deg>'][i],cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1],
               linewidths=0.7,colors='black',linestyle=':'),
        ax.hlines(cat['lat_HI<deg>'][i],cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
               linewidths=0.7,colors='black',linestyle=':'),
        ax.set_aspect('equal'),

        ax.set_xlabel('Galactic Longitude'),
#         ax.set_ylabel('Galactic Latitude'),


        ax.set_xlim(cat['lon_HI<deg>'][i][0]+1.2*rh[i],cat['lon_HI<deg>'][i][1]-1.2*rh[i]),
        ax.set_ylim(cat['lat_HI<deg>'][i][0]-1.2*rh[i],cat['lat_HI<deg>'][i][1]+1.2*rh[i]),


        xtick_locs = np.linspace(cat['lon_HI<deg>'][i][0]+1.*rh[i],
                              cat['lon_HI<deg>'][i][1]-1.*rh[i],3).astype('int')
        ax.set_xticks(xtick_locs,
               labels=[str(tck)+r'$\degree$' for tck in
                      xtick_locs]),
        
        
    ytick_locs = np.linspace(cat['lat_HI<deg>'][i][0]-0.9*rh[i],
                              cat['lat_HI<deg>'][i][1]+0.9*rh[i],3).astype('int')
    axs[0].set_yticks(ytick_locs,

               labels=[str(tck)+r'$\degree$' for tck in
                      ytick_locs]),

    axs[0].set_ylabel('Galactic Latitude')
    axs[1].set_yticks([])

    axs[2].hist(np.ravel(h1s_3r[i].data),
                bins=np.linspace(np.nanmean(h1s_2r[i].data)-3*np.nanstd(h1s_2r[i].data),
                np.nanmean(h1s_2r[i].data)+3*np.nanstd(h1s_2r[i].data),40),
                density=True,histtype='step',color='black',zorder=-10)
    
    nh, bins, patches  = axs[2].hist(np.ravel(h1s[i].data),
                bins=np.linspace(np.nanmean(h1s[i].data)-3*np.nanstd(h1s[i].data),
                np.nanmean(h1s[i].data)+3*np.nanstd(h1s[i].data),40),
                density=True)
    
    for num, patch in enumerate(patches):
        patch.set_facecolor(plt.cm.RdGy(num/len(bins)))
        
    
    nc, bins, patches = axs[2].hist(np.ravel(c1s[i].data),
                bins=np.linspace(np.nanmean(c1s[i].data)-2.2*np.nanstd(c1s[i].data),
                np.nanmean(c1s[i].data)+2.2*np.nanstd(c1s[i].data),40),
                density=True,)
        
    for num, patch in enumerate(patches):
        patch.set_facecolor(plt.cm.PiYG(num/len(bins)))
        patch.set_height(-patch.get_height())
    
    
    axs[2].axhline(0, color="k",linewidth=0.7)
    
    
#     nc, bins, patches = axs[2].hist(np.ravel(c1s[i].data),
#                 bins=np.linspace(np.nanmean(c1s[i].data)-2.2*np.nanstd(c1s[i].data),
#                 np.nanmean(c1s[i].data)+2.2*np.nanstd(c1s[i].data),40),
#                 density=True,histtype='step',color='black')
    
#     for num, patch in enumerate(patches):
#         patch.set_height(-patch.get_height())
    
    
    axs[2].set_xlabel(r'$\mathrm{\overline{V_{LOS}}\ [km\ s^{-1}]}$')
    axs[2].yaxis.tick_right()
    axs[2].set_ylim(1.2*min(-nc),1.2*max(nh))
    axs[2].set_xticks(np.around(np.linspace(
                 min(np.nanmean(c1s[i].data)-2.2*np.nanstd(c1s[i].data),
                     np.nanmean(h1s_2r[i].data)-3*np.nanstd(h1s_2r[i].data)),
                 max(np.nanmean(c1s[i].data)+2.2*np.nanstd(c1s[i].data),
                     np.nanmean(h1s_2r[i].data)+3*np.nanstd(h1s_2r[i].data)),5),1))
#     axs[2].set_yticks([])
    axs[2].set_yticks(np.around(np.array([min(-nc),0.5*min(-nc),
                       0,max(nh)*0.5,max(nh)]),1),
                     labels=np.abs(np.around(np.array([min(-nc),0.5*min(-nc),
                       0,max(nh)*0.5,max(nh)]),1)))
    
#     axs[2].set_yticks(np.around(np.linspace(min(-nc),max(nh),5),1))

    
    plt.subplots_adjust(wspace=0.05)
    fig.tight_layout()
    plt.savefig('mom1_'+lbls[i],bbox_inches='tight')
    plt.show()


cat.write('lastlastLAST_cat.fits',overwrite=True)


#look back through and check if i want to save/export anything
#then kill this notebook, make the scaling relations-y ones in a new one
#then send nia email with agenda: 20min ra stuff, 20min velocity gradient is from
#additional emission components, 20min looking over draft (or maybe ask for 1.5hr meeting
#and do 30, 30, 30? or propose 30/30, and written feedback on the draft ?
#also ! look back at outline in jane st notebook, other notebooks ive written out 
#'heres what i should say about blah'


rad_gmc_circ = np.sqrt(areas_npix/np.pi)


cat[cat['radius_GMC<pc>']>cat['radius_HI<pc>']]


la_r_gmc = cat[(np.nanstd(cat['lat_GMC<deg>'],axis=1)>np.nanstd(cat['lat_HI<deg>'],axis=1))|
   (np.nanstd(cat['lon_GMC<deg>'],axis=1)>np.nanstd(cat['lon_HI<deg>'],axis=1))]


la_r_ind = [0,2,5,6,9,10,14,17,19]


fig,axs=plt.subplots(3,3,figsize=(10,10),dpi=200)
for i,ind in enumerate(la_r_ind):
#     plt.contourf(np.log10(h0s[i].data),extent=[cat['lon_HI<deg>'][i][0],cat['lon_HI<deg>'][i][1],
#                                     cat['lat_HI<deg>'][i][0],cat['lat_HI<deg>'][i][1]])
    
    
    axs.flat[i].vlines(cat['lon_HI<deg>'][ind],
                       cat['lat_HI<deg>'][ind][0],
                       cat['lat_HI<deg>'][ind][1],
              linestyle='--',color='black')
    axs.flat[i].hlines(cat['lat_HI<deg>'][ind],
                       cat['lon_HI<deg>'][ind][0],cat['lon_HI<deg>'][ind][1],
               linestyle='--',color='black')
    
    
    axs.flat[i].vlines(cat['lon_GMC<deg>'][ind],
                       cat['lat_GMC<deg>'][ind][0],
                       cat['lat_GMC<deg>'][ind][1],
              linestyle='--',color='red')
    axs.flat[i].hlines(cat['lat_GMC<deg>'][ind],
                       cat['lon_GMC<deg>'][ind][0],
                       cat['lon_GMC<deg>'][ind][1],
               linestyle='--',color='red')
    
    
    axs.flat[i].contourf(np.log10(c0s[ind].data),
                         extent=[cat['lon_GMC<deg>'][ind][0],cat['lon_GMC<deg>'][ind][1],
                                 cat['lat_GMC<deg>'][ind][0],cat['lat_GMC<deg>'][ind][1]])
    
    axs.flat[i].set_aspect('equal')
    
    axs.flat[i].set_xlim(cat['lon_GMC<deg>'][ind]+np.array([3,-3]))
    axs.flat[i].set_ylim(cat['lat_GMC<deg>'][ind]+np.array([-3,3]))
    axs.flat[i].set_title(r'$R_H$ = '+str(round(cat['radius_HI<pc>'][ind],2))+\
                      r'pc, $R_{GMC}$ = '+str(round(cat['radius_GMC<pc>'][ind],2)),
                         fontsize='small')
                      

fig.tight_layout()
fig.supxlabel('Galactic Longitude',y=-0.05)
fig.supylabel('Galactic Latitude',x=-0.05)


fig.savefig('small_ra',bbox_inches='tight')

#     plt.show()


r'$R_H$ = '+str(round(cat['radius_HI<pc>'][ind],2))+\
                      r'pc, $R_{GMC}$ = '+str(round(cat['radius_GMC<pc>'][ind],2))


str(round(cat['radius_HI<pc>'][ind],2))


help(fig.supxlabel)


cat['lat_GMC<deg>'][19] = [20.,31.]


cat['lat_GMC<deg>'][19]


np.nanstd(cat['lon_HI<deg>'][19])


np.nanstd((30.25,20.625))


cat['lat_HI<deg>'][19]


np.nanmean((cat['loc'][19],cat['lic'][19]))


np.nanmean((cat['loh'][19],cat['lih'][19]))


np.nanmean((cat['boc'][19],cat['bic'][19]))


np.nanmean((cat['boh'][19],cat['bih'][19]))


cat['boh'][19],cat['bih'][19]


cat['boc'][19],cat['bic'][19]


642-565


648-560


cat['lat_GMC<deg>'][19]
