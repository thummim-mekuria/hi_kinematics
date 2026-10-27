#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datacube as dc


ct = pd.read_csv('physcl_prprts_catalog.csv')


weird_r = ct[ct['effective_radius<pc>']>ct['accumulation_radius<pc>']]


mekcat = pd.read_csv('Mekuria_catalog.csv',skiprows=2)


mekcat


plt.figure(figsize=(13,7))
plt.scatter(mekcat['l_c<deg>'],mekcat['b_c<deg>'],marker='s',s=50,alpha=0.7)
plt.scatter(weird_r['lon<deg>'],weird_r['lat<deg>'],marker='x',c=weird_r['dist<kpc>'],cmap='magma')
plt.colorbar(label='d [kpc]')
plt.xlim(185,-185)
plt.hlines(0,180,-180,color='grey',alpha=0.4)
plt.vlines(0,25,-25,color='grey',alpha=0.4)
plt.show()


weird_r.sort_values('dist<kpc>')


# 275 = corona
# 242,260 = ophiuchus
# 


from astrodendro import Dendrogram
dnd = Dendrogram.load_from('COGAL_dendrogram.fits')


plt.imshow(np.nansum(dnd[275].get_mask() * dnd.data,axis=0),\
           cmap=dc.cmap_grey0('viridis'),origin='lower',extent=[180,-180,-30,30])
plt.xlim(10,-10)
plt.ylim(-25,-10)


plt.imshow(np.nansum(dnd[242].get_mask() * dnd.data,axis=0),\
           cmap=dc.cmap_grey0('viridis'),origin='lower',extent=[180,-180,-30,30])
plt.imshow(np.nansum(dnd[260].get_mask() * dnd.data,axis=0),\
           cmap=dc.cmap_grey0('viridis'),origin='lower',extent=[180,-180,-30,30])
plt.scatter(-2,17,color='red')

plt.xlim(5,-12)
plt.ylim(10,25)


oph_bs = np.arange(10,25,0.125)


oph_ls = np.arange(15,-15,-0.125)


oph_ls_i = ((oph_ls- (180.125))/(-0.125)).astype('int')


oph_bs_i = ((oph_bs - (-29.875))/(0.125)).astype('int')


oph_bs_i,oph_ls_i


oph_field = dnd.index_map[:,319:438,1321:1560]


have_data_vbl = np.where(oph_field>-1)


bs = have_data_vbl[1] + 319
ls = have_data_vbl[2] + 1321
vs = have_data_vbl[0]


have_data_vbl = np.array((vs,bs,ls))


vs_kms = (vs*1.3)+(-49.4)


bs_deg = (bs*0.125)+(-30.125)


ls_deg = (ls*-0.125)+(180.125)


oph_sourcs = pd.DataFrame(np.array((ls_deg,bs_deg,vs_kms)).transpose(),columns=['lon','lon','lat'])


oph_sourcs['v_unc'] = np.full(len(oph_sourcs),0.6)


oph_sourcs['p_far'] = np.full(len(oph_sourcs),0.5)


oph_sourcs.to_csv('sources_info.inp')


from astropy.io import fits
dht = fits.open('DHT_CO.fits')[0]


plt.imshow(np.nansum(dnd[260].get_mask()*dnd.data,axis=0),\
           cmap=dc.cmap_grey0('viridis'),origin='lower',extent=[180,-180,-30,30])
plt.xlim(5,-10)
plt.ylim(10,20)
plt.colorbar()


oph_260 = dnd[260].get_mask()*dnd.data


ophbs,ophls = np.where(oph_260>0)[1:]


num_px_data = len(np.unique(np.where(oph_260>0)[1:],axis=1)[0])


num_px_data


import astropy.constants as c


# import astropy.units as u
ang_ar = ((u.deg*0.125)**2)*num_px_data


r_l_oph = 19.7*u.pc
d_l_oph = 144*u.pc
d_r_oph = 18470*u.pc
m_l_oph = 40.5*1e3*c.M_sun
m_r_oph = 1.594765e38*u.kg


np.sqrt((((ang_ar.to(u.rad**2)).value)*(d_l_oph**2))/np.pi)


np.sqrt((((ang_ar.to(u.rad**2)).value)*(d_r_oph**2))/np.pi)


876.37251/6.832574


d_r_oph/d_l_oph


m_r_oph/m_l_oph


(d_r_oph/d_l_oph)**2


plt.imshow(np.nansum(dnd[275].get_mask()*dnd.data,axis=0),\
           cmap=dc.cmap_grey0('viridis'),origin='lower',extent=[180,-180,-30,30])
plt.xlim(5,-5)
plt.ylim(-25,-15)
plt.colorbar()


r_l_cra = 3.2*u.pc
d_l_cra = 130*u.pc
d_r_cra = 8340*u.pc
m_l_cra = 1.2*1e3*c.M_sun
m_r_cra = 1.214388e37*u.kg


m_r_cra/m_l_cra


d_r_cra/d_l_cra


(d_r_cra/d_l_cra)**2


m_r_ora/m_l_ora


d_r_ora/d_l_ora


(d_r_ora/d_l_ora)**2


r_l_ora = 23*u.pc
d_l_ora = 432*u.pc
d_r_ora = 730*u.pc
m_l_ora = 85.2*1e3*c.M_sun
m_r_ora = 4.302619e35*u.kg


ct.iloc[120]


weird_r


pdf16026 = pd.read_csv('16026_final_distance_pdf.dat',delim_whitespace=True,skiprows=2,header=None)


d16026 =  pdf16026[0]


pd16026 =  pdf16026[1]


plt.plot(d16026,pd16026)


oph_dive = pd.read_csv('summary.prt',skiprows=4,delim_whitespace=True,header=None)


clmns = pd.read_csv('summary.prt',nrows=3).iloc[1]


cl_names = clmns.to_list()[0].split()[1:]


cl_names[4:8] = ['Dist1.','+/-1','Prob1','Arm1']
cl_names[8:] = ['Dist2.','+/-2','Prob2','Arm2']


oph_dive.columns=cl_names


plt.scatter(oph_dive['Long.'],oph_dive['Lat.'],marker='.',s=5,c=oph_dive['Dist1.'])
plt.xlim(15,-15)
plt.colorbar()


plt.show()


oph_dive.query("Arm1 != '...?' and Arm1 != '...??'")


oph_dive[oph_dive['Long.']== -5.25]


oph_dive


oph_dive[oph_dive['Arm1']=='...?']['Arm1'] = 0


oph_dive[oph_dive['Arm1']=='...?']


oph_dive


plt.scatter(oph_dive['Long.'],oph_dive['Lat.'],marker='.',s=1,c=oph_dive['Arm1'])
plt.xlim(15,-15)
plt.colorbar()

plt.scatter(-5.25,16.25,color='red')

plt.show()


pdf1852 = pd.read_csv('Reid_BDC/1852_final_distance_pdf.dat',delim_whitespace=True,skiprows=2,header=None)


pdf1887 = pd.read_csv('Reid_BDC/1887_final_distance_pdf.dat',delim_whitespace=True,skiprows=2,header=None)


pdf1623 = pd.read_csv('Reid_BDC/1623_final_distance_pdf.dat',delim_whitespace=True,skiprows=2,header=None)


plt.scatter(pdf1887[0],pdf1887[1]/np.max(pdf1887[1]),marker='+')
plt.scatter(pdf1852[0],pdf1852[1]/np.max(pdf1852[1]),marker='+')
plt.scatter(pdf1623[0],pdf1623[1]/np.max(pdf1623[1]),marker='+')


# plt.hlines(0.06,-1,25)
# plt.xlim(0,1)


pdf1887_kd = pd.read_csv('Reid_BDC/18870_kinematic_distance_pdf.dat',delim_whitespace=True,skiprows=2,header=None)


pd.read_csv('Reid_BDC/1887',\
            delim_whitespace=True)#,skiprows=2,header=None)


weird_r
