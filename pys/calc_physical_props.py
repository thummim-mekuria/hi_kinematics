#===================================================================================
import pandas as pd
from astrodendro import Dendrogram
import astropy.constants as c
import astropy.units as u
import numpy as np
import clouds as cld
import datacube as dc
#===================================================================================


#===================================================================================
cat = pd.read_csv('src_sm_ppv_merged.csv')
dnd = Dendrogram.load_from('COGAL_dendrogram.fits')
#===================================================================================


#===================================================================================
# idxs = cat['_idx'].to_numpy()
# msks = [dnd[idxs[i]].get_mask() for i in range(len(cat))]

# cube = dnd.data
# clouds = [(cube*msks[i]*u.K) for i in range(len(cat))]
# ws = [np.nansum(1.3*u.km/u.s * clouds[i],axis=0) for i in range(len(cat))]

# x_co = 1.97e20 * u.s/(u.K*u.km*u.cm**2)
# ns = [ws[i]*x_co for i in range(len(cat))]

# dSs = [ns[i]*c.m_p*2*1.36 for i in range(len(cat))] 
# Ss = [np.nansum(dSs[i]) for i in range(len(cat))]

# ds = cat['dist<kpc>'].to_numpy()* u.kpc
# px_area = (0.125)**2 * (u.deg**2).to(u.rad**2)

# As = [ds[i]*px_area for i in range(len(cat))]

# Ms = [Ss[i]*As[i] for i in range(len(cat))]
#===================================================================================


#===================================================================================
#forloop is quicker than ^
Ms = np.zeros(len(cat))
As = np.zeros(len(cat))
RAs = np.zeros(len(cat))
RA_ss = np.zeros(len(cat))
sigmas = np.zeros(len(cat))
#===================================================================================


#===================================================================================
for i,idx in enumerate(cat['_idx']):
	msk = dnd[idx].get_mask()
	cloud = msk*dnd.data
	#-------------------------------------------------------------------------------
	d = cat['dist<kpc>'].loc[i] * u.kpc
	b = cat['lat<deg>'].loc[i] * u.deg
	#-------------------------------------------------------------------------------
	dS, dA_px = cld.surface_density(cloud,'CO')

	M = cld.mass(dS,d,dA_px)
	Ms[i] = M.value

	A = np.nansum(d**2 * dA_px * (dS!=0)).to(u.pc**2)
	As[i] = A.value
	#-------------------------------------------------------------------------------
	ra, s = cld.accumulation_radius(M,d,b)
	RAs[i] = ra.value
	RA_ss[i] = s.value
	#-------------------------------------------------------------------------------
	spec = dc.spectrum1D(cloud)
	sig = dc.dispersion_1sigma(spec)

	if len(sig) == 0 or (np.argmax(spec)<5*sig) or (np.argmax(spec)>len(spec)-5*sig):
		sigmas[i] = np.nan #the peak is invalid /too close to the edges of the v axis 
	else:
		sigmas[i] = sig * 1.3 #to put in km/s instead of pixel units
#===================================================================================


#===================================================================================
cat['mass<kg>'] = Ms
cat['area<pc2>'] = As
cat['effective_radius<pc>'] = np.sqrt(As/np.pi)
cat['effective_radius<deg>'] = np.sqrt(cat['area_sky<deg2>']/np.pi)
cat['accumulation_radius<pc>'] = RAs
cat['accumulation_radius<deg>'] = RA_ss
cat['velocity_dispersion<km/s>'] = sigmas
#===================================================================================


#===================================================================================
cat.to_csv('physcl_prprts_catalog.csv',index=False)
#===================================================================================





