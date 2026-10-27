import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from astropy.table import QTable

import datacube as dc

from astrodendro import Dendrogram

dend = Dendrogram.load_from('COGAL_dendrogram.fits')

def get_leaves(array):
	out = [dend[i] for i in range(len(dend)) if i in array]
	out = np.array(out)
	return out

# cat = pd.read_csv('physcl_prprts_catalog.csv')
def find_in_cat(array):
	out = [i for i in range(len(cat)) if cat['_idx'][i] in array]
	out_df = cat.iloc[out]
	return out_df

def v_kms(arrays):
    output = []
    for ar in arrays:
        this = ar*1.3-49.4
        output.append(this)
    return output

def b_deg(arrays):
    output = []
    for ar in arrays:
        this = ar*0.125-30.125
        output.append(this)
    return output

def l_deg(arrays):
    output = []
    for ar in arrays:
        this = ar*(-0.125)+180.125
        output.append(this)
    return output

def vbl_arrays(array):
	vs = []
	bs = []
	ls = []
	ns = []
	for ar in array:
		v,b,l = dend[ar].indices()
		vs.append(v)
		bs.append(b)
		ls.append(l)
		ns.append([str(ar)]*len(v))

	print('made vs bs etc')
	vslong = np.concatenate(vs)
	bslong = np.concatenate(bs)
	lslong = np.concatenate(ls)

	vkms = v_kms(vslong)
	bdeg = b_deg(bslong)
	ldeg = l_deg(lslong)
	print('put them in physical units')
	names = np.concatenate(ns)
	return names,ldeg,bdeg,vkms

def make_source(names,lons,lats,vels):
	sources = pd.DataFrame(np.array((names,lons,lats,vels)).transpose(),\
		columns = ['name','lon','lat','vel'])

	sources['v_unc'] = [0.6]*len(sources)
	sources['p_far'] = [0.5]*len(sources)
	print('made dataframe',len(sources))
	sources.to_csv('sources_info.inp',index=False)

	return sources

def read_and_segmentsummaryprt(names):
	out = pd.read_csv('summary.prt',skiprows=4,
		header=None,delim_whitespace=True)
	

	nice_out = out[[0,1,2,4,6,7,8,10,11]]
	nice_heads = ['Long.','Lat.','Vlsr','Dist1','Prob1',
	'Arm1','Dist2','Prob2','Arm2']

	nice_out.columns = nice_heads
	nice_out['name'] = names

	nice_out = nice_out[['name','Long.','Lat.','Vlsr','Dist1',\
	'Prob1','Arm1','Dist2','Prob2','Arm2']]
	
	# nice_out.to_csv('nice_out.csv',index=False)


	demarc = [i for i in range(len(nice_out)-1) \
	if nice_out['name'][i]!=nice_out['name'][i+1]]

	
	# print('\n\n\n')
	# for i in range(len(nice_out)-1):
		# print(i,nice_out['name'][i],nice_out['name'][i+1])
		# print('\n\n')

	inn = np.array(demarc) + 1
	demarc = np.concatenate(([0],inn,[len(nice_out)]))

	
	for i in range(1,len(demarc)):
		stt = demarc[i-1]
		end = demarc[i]

		sect = nice_out.iloc[range(stt,end)]
		print(sect)
		sect.to_csv('idx'+str(sect['name'].iloc[0])+'.csv',index=False)






if __name__=="__main__":
	# this is the first part
	dangr = [167,198,217,238,248,249,271,281]
	n,l,b,v = vbl_arrays(dangr)
	srcs = make_source(n,l,b,v)
	read_and_segmentsummaryprt(srcs['name'])

	#now distance calculations have been done
	#this is the second part





