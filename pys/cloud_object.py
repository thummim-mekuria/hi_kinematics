#might not need all these imports, but they're there for now
import numpy as np
import pandas as pd
import astropy.units as u 
import astropy.constants as c
import matplotlib.pyplot as plt
import datacube as dc
import clouds as cld
from scipy.signal import peak_widths, find_peaks
from copy import deepcopy

from astropy.io import fits
from astrodendro import Dendrogram


from matplotlib import cm
from matplotlib.colors import ListedColormap



"""
might change how i do this later because this is all quickly becoming a tall jenga tower
for now, instantiate a Cloud object by inputting the dictionary that is output by get_parent_HI
which is in the form of {comp: ___, 0th: ___, spectrum: ___} 
and comp is in the form of {cube: ___, vels: ___, lats: ___, lons: ___}
"""


class Cube:
	def __init__(self, unpacked_dict):
		"""
		unpacked_dict is the thing that *used to be* output by dc.unpack_fits_cube & cld.extrct_cloud
		after defining this class, what instead happens is that the dictionaries get created in those functions, 
		THEN that dictionary gets used to instantiate a Cube object (tbh this is prob v similar to spectral_cube but
		i wanna figure this out for myself so it makes) from that dictionary
		"""
		self.units = unpacked_dict['units']

		self.cube = unpacked_dict['cube']
		self.vels = unpacked_dict['vels']
		self.lats = unpacked_dict['lats']
		self.lons = unpacked_dict['lons']


		# cube = unpacked_dict['cube']*self.units['data']
		# self.cube = cube

		# vels = unpacked_dict['vels']*self.units['vels']
		# self.vels = vels

		# lats = unpacked_dict['lats']*self.units['lats']
		# self.lats = lats

		# lons = unpacked_dict['lons']*self.units['lons']
		# self.lons = lons

		self.data_header = unpacked_dict['header']


		#then calculate the properties that depend on ^ that come up a lot
		v1 = self.vels[1]
		v0 = self.vels[0]
		self.dv = (v1 - v0)*self.units['vels']

		b1 = self.lats[1]
		b0 = self.lats[0]
		self.db = (b1 - b0)*self.units['lats']

		l1 = self.lats[1]
		l0 = self.lats[0]
		self.dl = (l1 - l0)*self.units['lons']


		# self.zeroth_moment =  np.nansum(self.cube*self.dv.value, axis=0)


	

	############ a series of calculation methods ##################
	
	def calculate_0th_moment(self):
		zeroth_moment =  np.nansum(self.cube*self.units['data']*self.dv, axis=0)
		return zeroth_moment


	def calculate_rms_noise(self):
		rms_noise = np.nanstd(self.cube)
		return rms_noise



	def clone(self):	
		new = deepcopy(self)
		return new

	######### BEGIN: subcube (used to be extract cloud in clouds.py) #########

	def create_subcube (self, coords, v_guess = None, trim_vel = True, n_sigma_vchannels = 10):
		#----------------- step 1: extract based on spatial coords ---------------

		b0, b1, l0, l1 = coords
		#for now, this method is using upper & lower l and b bounds
		#will be nice to eventually get to a point where i take in a central l, b coordinate, then
		#do some dendogram decomposition of that general area to figure out what the size of the most coherent structure is in tht
		#vicinity, then define the cube bounds that way (or at least do this for CO & GMC-less clouds, and use accumulation radius)
		#for GMC-full HI clouds

		       
		#convert physical coords -> index locations within the larger cube
		lat_inds = cld.ind_ang('get index', self.lats, (b0,b1))
		lon_inds = cld.ind_ang('get index', self.lons, (l0,l1))

		#fixes index orders in case data is in decreasing order
		lat_inds.sort()
		lon_inds.sort()
		    
		b0, b1 = lat_inds #get each element separately
		l0, l1 = lon_inds
		    
		scube = self.cube[:, b0:b1, l0:l1]
		sublats = self.lats[b0:b1]
		sublons = self.lons[l0:l1]


		#------------------------------------------------------------------





		# #---------------------- step 2: extract based on signal --------------
		d = Dendrogram.compute(scube, )





		spec = dc.spectrum1D(scube)
		v_dispersion = dc.dispersion_1sigma(spec)

		if v_guess:
			v_cent = v_guess
		else:
			v_cent = np.nanargmax(spec)

		v0 = v_cent - int(0.5*n_sigma_vchannels*v_dispersion)
		v1 = v_cent + int(0.5*n_sigma_vchannels*v_dispersion)



		#save the cube across the full range of v & redefine the bounds of the cloud cube
		scube_long = scube
		scube = scube_long[v0:v1, :,:]

		subvels = self.vels[v0:v1]


		if trim_vel:
			output_dict = {'cube':scube, 'vels':subvels, 'lats':sublats, 'lons':sublons,'units':self.units, 'header': self.data_header}

		else:
			output_dict = {'cube':scube_long, 'vels':self.vels, 'lats':sublats, 'lons':sublons,'units':self.units, 'header': self.data_header}
		return Cube(output_dict)

	######### END: subcube #########
		



class Cloud:
	def __init__(self, full_dict):
		"""
		might change how i do this later because this is quickly becoming a tall jenga tower
		but for now, instantiate a Cloud object by inputting the dictionary that is output by get_parent_HI
		which is in the form of {comp: ___, 0th: ___, spectrum: ___} 
		and comp is in the form of {cube: ___, vels: ___, lats: ___, lons: ___}
		"""

