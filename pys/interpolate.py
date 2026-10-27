
import numpy as np
from numba import njit




@njit
def interpolate1D(line,max_gap=1):
	ln = len(line)

	output = np.zeros(ln)

	output[0] = line[0]
	output[-1] = line[-1]
	for n in range(ln-1):
		num = line[n]
		if num <= 0. or np.isnan(num):
			x1 = False
			x2 = False
			for k in range(1,max_gap+1):
				# print(k)
				if not np.isnan(line[n-k]) and line[n-k] > 0:
					x1 = line[n-k]
					k1 = n-k
					# print(n,k,n-k,x1)

			for k in range(1,max_gap+1):
				if line[n+k] > 0 and not np.isnan(line[n+k]):
					x2 = line[n+k]
					k2 = n+k
				
			if x1 and x2:
				output[n] = (x1+x2)/(k2-k1)

		else:
			output[n] = num

	return output



@njit
def interpolate_vlb(cube,max_gap=1):
	#order of interpolation:
	#across spectral axis, then across l, then across b
	#assume the cube axes are 0=v,1=b,2=l
	lbv = np.transpose(cube)
	spec_smooth = np.zeros_like(lbv)

	for l in range(len(lbv)):
		for b in range(len(lbv[l])):
			spec = lbv[l][b]
			spec_smooth[l][b] = interpolate1D(spec,max_gap)

	vbl = np.transpose(spec_smooth)
	l_smooth = np.zeros_like(vbl)
	for v in range(len(vbl)):
		for b in range(len(vbl[v])):
			lon = vbl[v][b]
			l_smooth[v][b] = interpolate1D(lon,max_gap)


	vlb = np.swapaxes(l_smooth,1,2)
	b_smooth = np.zeros_like(vlb)
	for v in range(len(vlb)):
		for l in range(len(vlb[v])):
			lat = vlb[v][l]
			b_smooth[v][l] = interpolate1D(lat,max_gap)


	interp_cube = np.swapaxes(b_smooth,1,2)

	return interp_cube







