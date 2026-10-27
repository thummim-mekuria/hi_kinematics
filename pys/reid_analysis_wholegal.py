import numpy as np
import pandas as pd

from astropy import units as u
from astropy import constants as c

def make_source(names,lons,lats,vels):
	sources = pd.DataFrame(np.array((names,lons,lats,vels)).transpose(),\
		columns = ['name','lon','lat','vel'])

	sources['v_unc'] = [0.6]*len(sources)
	sources['p_far'] = [0.5]*len(sources)
	print('made dataframe',len(sources))
	sources.to_csv('sources_info.inp',index=False)

	return sources

def read_summaryprt(names):
	out = pd.read_csv('summary.prt',skiprows=4,
		header=None,delim_whitespace=True)
	

	nice_out = out[[0,1,2,4,6,7,8,10,11]]
	nice_heads = ['Long.','Lat.','Vlsr','Dist1','Prob1',
	'Arm1','Dist2','Prob2','Arm2']

	nice_out.columns = nice_heads
	nice_out['name'] = names

	nice_out = nice_out[['name','Long.','Lat.','Vlsr','Dist1',\
	'Prob1','Arm1','Dist2','Prob2','Arm2']]
	
	nice_out.to_csv('nice_out.csv',index=False)

	return nice_out




if __name__=="__main__":
	# this is the first part
	lons = np.arange(180,-181,-2)
	lats = np.arange(-30,31,4)
	vels = np.arange(-50,51,2)
	grid = np.array(np.meshgrid(lons,lats,vels))

	l,b,v = grid.transpose().reshape(-1,3).transpose()
	names = np.arange(147696)

	srcs = make_source(names,l,b,v)
	# read_and_segmentsummaryprt(srcs['name'])

	#now distance calculations have been done
	#this is the second part





