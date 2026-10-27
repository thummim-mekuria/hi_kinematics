#===================================================================================
from astropy.io import fits
import astropy.units as u
import numpy as np
from scipy.signal import peak_widths, find_peaks

from matplotlib import cm
from matplotlib.colors import ListedColormap

import cloud_object as clob
from numba import njit
import pandas as pd
#===================================================================================

#===================================================================================
def unpack_fits_cube(name, cube_format):
    
    """
    if cube_format = 'DHT', then unpacks data cubes that are in the form of 
    Dame's CO survey outputs:
    https://lweb.cfa.harvard.edu/rtdc/CO/CompositeSurveys/
    
    if cube_format = 'HI4PI', then unpacks datacubes that are in the form of 
    the HI4PI survey output:
    https://cdsarc.cds.unistra.fr/ftp/J/A+A/594/A116/
    """
    
    
    #----- open the image and the header from fits ---------------------------
    image = fits.open(name)[0].data
    header = fits.open(name)[0].header
    #-------------------------------------------------------------------------
       
    
    #--- extract information labelling velocity, longitude, latitude axes ----
    
    #figure out which axes are which
    for i in range(3):
        num = i+1 #axes are 1, 2, or 3
        if header['CTYPE'+str(num)] == 'GLON-CAR':
            l_num = num
    
        if header['CTYPE'+str(num)] == 'GLAT-CAR':
            b_num = num
            
    
    
    if cube_format == 'DHT':
        for i in range(3):
            num = i+1 #axes are 1, 2, or 3
            if header['CTYPE'+str(num)] == 'VELO-LSR':
                v_num = num
                
        #hard-coding because this information is in the comments        
        velUnit = u.km/u.s 
        latUnit = u.deg
        lonUnit = u.deg
        dataUnit = u.K
                
    elif cube_format == 'HI4PI':
        for i in range(3):
            num = i+1 #axes are 1, 2, or 3
            if header['CTYPE'+str(num)] == 'VRAD':
                v_num = num
                
        velUnit = u.Unit(header['CUNIT'+str(v_num)])
        latUnit = u.Unit(header['CUNIT'+str(b_num)])
        lonUnit = u.Unit(header['CUNIT'+str(l_num)])
        dataUnit = u.Unit(header['BUNIT'])


    elif cube_format == 'LAB':
        for i in range(3):
            num = i+1 #axes are 1, 2, or 3
            if header['CTYPE'+str(num)] == 'VELO-LSR':
                v_num = num
                
        #units are put in in all caps, so need to first convert to lowercase
        vunit_str = header['CUNIT'+str(v_num)].lower()
        velUnit = u.Unit(vunit_str)

        latunit_str = header['CUNIT'+str(b_num)].lower()
        latUnit = u.Unit(latunit_str)

        lonunit_str = header['CUNIT'+str(l_num)].lower()
        lonUnit = u.Unit(lonunit_str)
        
        dataUnit = u.Unit(header['BUNIT'])
    #------------------------------------------------------------------------

                
                    

    #----- pull out vel-related, lat-related,lon-related values from header 
    # then create the arrays holding the real/physical axis values -----------
    velLen = header['NAXIS'+str(v_num)]
    velMin = header['CRVAL'+str(v_num)]
    velDel = header['CDELT'+str(v_num)]
    velRef = header['CRPIX'+str(v_num)]    
    

                
    latLen = header['NAXIS'+str(b_num)]
    latMin = header['CRVAL'+str(b_num)]
    latDel = header['CDELT'+str(b_num)]
    latRef = header['CRPIX'+str(b_num)] 
    
        
    lonLen = header['NAXIS'+str(l_num)]
    lonMin = header['CRVAL'+str(l_num)]
    lonDel = header['CDELT'+str(l_num)]
    lonRef = header['CRPIX'+str(l_num)]

    vels = np.arange(velLen)
    vels = (vels - velRef)*velDel + velMin

    glat = np.arange(latLen)
    glat = (glat - latRef)*latDel + latMin
    
    glon = np.arange(lonLen)
    glon = (glon - lonRef)*lonDel + lonMin
        
    
    units = {'vels':velUnit, 'lats':latUnit, 'lons':lonUnit, 'data': dataUnit} 
    #------------------------------------------------------------------------
     
    
    #----- rearrange fits cube axes ------------------------------------------
    
    #for how numpy displays / manipulates arrays, we want 
    #axis0=velocity, axis1=latitude, & axis2=longitude
    
    #!!!!! this whole thing is assuming that all three axes have different lengths, 
    #will run into a problem if that is not the case
    
    oldArrangement = image
    x, y, z = oldArrangement.shape
    oldAxes = np.array((x,y,z))
    
    image = np.zeros((velLen,latLen,lonLen))
    
    
    #first check if it is already in the right dimensions
    if ((x == velLen) and (y==latLen) and (z==lonLen)):
        image = oldArrangement   
    
    #check conditions that would only require one swap
    elif ((x == velLen) or (y==latLen) or (z==lonLen)):
        if x == velLen:
            image = np.swapaxes(oldArrangement, 1, 2)
            
        elif y == latLen:
            image = np.swapaxes(oldArrangement, 0, 2)
        
        elif z == lonLen:
            image = np.swapaxes(oldArrangement, 0, 1)
   
    else: #two swaps need to be made
        #indices are to get the actual number
        shuffle = np.where(oldAxes == velLen)[0][0] 
        
        #if shuffle = 1 that means axes need to shuffle <- way
        #if shuffle = 2 that means axes need to shuffle -> that way
        image = np.swapaxes(oldArrangement, 0, 1)
        
        if shuffle == 1: 
            image = np.swapaxes(image, 1, 2)
        elif shuffle == 2:
            image = np.swapaxes(image, 0, 2)
            
    #---------------------------------------------------------------------------
    
    

    output = clob.Cube({'cube': image, 'vels': vels, 'lats':glat, 'lons':glon, 'units':units, 'header':header})

    return output
#===================================================================================


#===================================================================================
@njit
def spectrum1D (array):
    """
    calculates 1D spectrum from a cube or 2d array
    assume axis = 0 contains the information to be spectrum-ed
    ie axis = 0 is the spectral axis
    
    works for both 2D and 3D arrays
    """
    l = len(array)
    spectrum = np.zeros(l)

    for i, ii in enumerate(array):
        spectrum[i] = np.nanmean(ii)



    # dim = len(array.shape) 
    
    # if dim == 2:
    #     l = len(array)

    #     spectrum = np.zeros(l)
    #     for i,ii in enumerate(array):
    #         spectrum[i] = np.nanmean(ll)
      
    # elif dim == 3:
    #     l, j = array.shape[:-1]
    #     spectrum = np.zeros(l,)
    #     for i,ii in enumerate(array):
    #         for k,kk in enumerate(ii):
    #             spectrum = np.nanmean(array, axis = (1,2))

         
    return spectrum
#===================================================================================


#===================================================================================
#might eventually have to rethink how i do this in the cases with multiple peaks
def dispersion_1sigma (spectrum, threshold=0.75):
    '''
    calculates the 1 sigma dispersion from a spectrum
    threshold how tall peaks must be (factor of max value) in order to count
    '''
    sqrt2ln2x2 = 2*(np.sqrt(2*np.log(2)))
    
    maxim = np.nanmax(spectrum)
    
    p = find_peaks(spectrum, height = threshold * maxim)
    #first unpack what find_peaks returns
    peak = p[0] #array of peak locations
    peak_info = p[1] #array of info about each peak
    
    if len(peak) > 1:#ie if there are multiple peaks found that are above the threshold
        #look at the peak heights & find the location of the tallest one
        i = np.argmax(peak_info['peak_heights']) 
        #then use that tallest peak for fwhm purposes
        peak = [peak[i]] #use a list because peak_widths wants an iterable
        
    fwhm = peak_widths(spectrum, peak)[0] #the first from peak_width output is fwhm
    sigma = fwhm/sqrt2ln2x2 

    #returns the 1-sigma dispersion in pixel/index units
    return sigma  
    # return fwhm  
#===================================================================================


#===================================================================================
def cmap_grey0 (default):
    """
    define a color-map where the minimum value is transparent 
    #& the rest is a matplotlib default map
    default = string corresponding to any valid matplotlib maps
    from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html

    following along example from:
    https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html
    not sure why get_cmap has 256 passed in
    """

    mpl_map = cm.get_cmap(default, 256)
    newcolors = mpl_map(np.linspace(0, 1, 256))
    
    min_val = np.array([240/256,240/256,240/256, 0]) # a lightish grey -- actually transparent now
    newcolors[:1, :] = min_val #set the colormap so that the min value is the grey defined above
    
    
    new_cmap = ListedColormap(newcolors)
    
    return new_cmap
#===================================================================================


#===================================================================================
@njit
def numba_nansum_ax0(cube):
    x,y,z = cube.shape
    trans = np.transpose(cube)
    listsum = np.zeros((y*z))
    
    for zz in range(z):
        this = trans[zz]
        for yy in range(y):
            thisthis = this[yy]
            
            listsum[y*zz + yy] = np.nansum(thisthis)
    
    cube0sum = np.transpose(listsum.reshape((z,y)))
            
    
    return cube0sum
#===================================================================================


#===================================================================================
#eq 6 from MD
@njit
def weighted_l(cube,lons,dv):
    wco = numba_nansum_ax0(cube*dv) #include data units at somepoint
    wco_tot = np.nansum(wco)
    lfront = np.transpose(cube)
    num = 0
    for l,ll in enumerate(lons):
        this_slice = lfront[l] 
        wcoi = np.nansum(this_slice*dv*ll) #include data units at somepoint
        num += wcoi
        
    lc = num/wco_tot
    return lc
#-----------------------------------------------------------------------------------
#eq 7 from MD
@njit
def weighted_b(cube,lats,dv):
    wco = numba_nansum_ax0(cube*dv) #include data units at somepoint
    wco_tot = np.nansum(wco)
    bfront = np.swapaxes(cube,0,1)
    num = 0
    for b,bb in enumerate(lats):
        this_slice = bfront[b] 
        wcoi = np.nansum(this_slice*dv*bb) #include data units at somepoint
        num += wcoi
        
    bc = num/wco_tot
    return bc
#-----------------------------------------------------------------------------------
#eq 11 from MD
@njit
def weighted_v(cube,vels,dv):
    #except instead of the form in 11, it looks more like the forms in 6 and 7 
    #since T_B is not from gaussian decomp
    wco = numba_nansum_ax0(cube*dv) #include data units at somepoint
    wco_tot = np.nansum(wco)
    
    num = 0
    for v,vv in enumerate(vels):
        this_slice = cube[v] 
        wcoi = np.nansum(this_slice*dv*vv) #include data units at somepoint
        num += wcoi
        
    vc = num/wco_tot
    return vc
#-----------------------------------------------------------------------------------
@njit
def weighted_vbl(cube,axes,ress):
    vels, lats, lons = axes
    dv, db, dl = ress
    l = weighted_l(cube, lons, dv)
    b = weighted_b(cube, lats, dv)
    v = weighted_v(cube, vels, dv)
    
    return v,b,l
#-----------------------------------------------------------------------------------
@njit
def weighted_coords(data, vbl_axes, ax_ress): 
    #prob very imprecise, since the clusters are spread out in l, b, v space
    #might need to do multiple distance calculations???

    N = len(data)
    IDs, ls, bs, vs = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    for i in range(N):
        cb = data[i]
        v,b,l = weighted_vbl(cb,vbl_axes,ax_ress)
        IDs[i] = i
        ls[i] = l
        bs[i] = b
        vs[i] = v

    return IDs, vs, bs, ls
#===================================================================================


#===================================================================================
def generate_sources(colnames, data):
    #assume data is a list of arrays& each corresponds to a column from colnames
    
    if colnames[0][0] != '!':
        colnames[0] = '! '+ colnames[0]
        
    if len(colnames) < 6: #assume that only id, lon, lat, and vlsr info is given
        N = len(data[0])
        pfar = np.full(N, 0.5)
        v_unc = np.full(N, 0.7)
        
        data = data+[v_unc, pfar]
        colnames = colnames + ['v_unc', 'pfar']
    
    df = pd.DataFrame(np.transpose(data),columns = colnames)
    #add second row thats commented out 
    df.columns = [df.columns.to_list(),['!'] + [' ']*(len(df.columns)-1)]
    
    df.to_csv('sources_info.inp',sep=' ',index=False)
    
    return df
#===================================================================================





