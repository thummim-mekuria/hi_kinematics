#!/usr/bin/env python
# coding: utf-8

def surface_density (clob, data_type, W_thresh = 1.5, v_channel_width = 25, show_Sigma = True):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from extract_cloud/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    and whether to assume that the spectrum is noisy / will need to isolate velocity peak or not (this is mostly
    for HI clouds)
    3) W_thresh (using value from lewis, lada, dame, but think about whether i wanna change it)

    '''
    spec = dc.spectrum1D(clob.cube)
    plt.plot(spec)
    plt.show()
    
    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/2203.07480.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        thresh = np.nanmin(spec)
        

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        thresh = np.nanmean(spec) + 3 * np.nanstd(spec)
    
        
    p = find_peaks (spec, height = thresh) 
    #something counts as a peak only if its above the mean-3sigma value 
    #this is because the assumption here is that the majority of the velocity channels
    #are noisy (ie a smaller fraction of the channels have data vs noise. this means that 
    #the mean of the values across the spectral axis is telling you about noise, and 
    #any value that is within the 3sigma of the noise mean is also noise. ie what counts as
    #a velocity channel containing signal is one whose value is > mean + 3sigma)
    #but then we run into a problem when every pixel contains data / the noise has already been clipped away
    print(p)
    peak = p[0] #array of peaks locations
    peak_info = p[1] #array of peak heights
    
    if len(peak) > 1:
        i = np.argmax(peak_info['peak_heights']) #if there are multiple peaks, pick the tallest one #this is problematic though if there is blending
        peak = [peak[i]]
    vel_center = peak[0]

    #this assumes that vel_center is at a place along the spectral axis that there is
    # still data +/- v_channel_width/2 around it ... seems like an ok assumption for now but might 
    # run into problems at some point
    v_channel_start =  vel_center - int(v_channel_width/2)
    v_channel_end = vel_center + int(v_channel_width/2)

    subcube = (clob.cube[v_channel_start:v_channel_end+1, :, :])*clob.units['data']

    #calculate integrated intensity
    W = np.nansum(subcube*clob.dv, axis=0)
    # W = np.nansum(cube*dv, axis=0)

    
    #use W_thresh to create a mask:
    W_mask = W.value > W_thresh
    W_old = W
    W = W_old * W_mask
    
    #calculate column density from the integrated intensity & conversion factor
    N = W * conversion_factor
    
    #define surface density here
    dSig = mmw * (c.m_p.to(u.g)) * N
    
    
    #calculate on sky area & cloud area at each pixel
    dAsky = np.abs(clob.dl * clob.db)
    

    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        kw_imshow = {'cmap': cmap0, 'origin':'lower'}

    
        dSig_Msunpc2 = dSig.to(c.M_sun/u.pc**2)

        fig, ax = plt.subplots (figsize = (15, 10), dpi = 200)
        Sig_plot = ax.imshow(dSig_Msunpc2.value, extent = [clob.lons[0], clob.lons[-1], \
                                                clob.lats[0], clob.lats[-1]], **kw_imshow)


        cbar = fig.colorbar(Sig_plot, ax=ax, label = Sig_label +r' [$M_\odot \ pc^{-2}$]', orientation = 'vertical')

        if data_type == '21cm':
            KMT_thresh = plt.contour(clob.lons, clob.lats, dSig_Msunpc2.value,  linewidths = 0.7, \
            levels=[10, 100], colors = ['xkcd:red', 'xkcd:orange']) #atomic-to-molecular transition for HI maps 
            cbar.add_lines(KMT_thresh)


        ax.set_xlabel(r'Galactic Longitude [$\degree$]')
        ax.set_ylabel(r'Galactic Latitude [$\degree$]')
        ax.set_title(Sig_label)


        plt.show()

    
    return {'surface density per pixel': dSig, 'pixel area': dAsky}


def surface_density (clob, data_type, W_thresh = 1.5, v_channel_width = 25, show_Sigma = True):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from extract_cloud/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    and whether to assume that the spectrum is noisy / will need to isolate velocity peak or not (this is mostly
    for HI clouds)
    3) W_thresh (using value from lewis, lada, dame, but think about whether i wanna change it)

    '''
    spec = dc.spectrum1D(clob.cube)
    plt.plot(spec)
    plt.show()
    
    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/2203.07480.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        thresh = np.nanmin(spec)
        

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        thresh = np.nanmean(spec) + 3 * np.nanstd(spec)
    
        
    p = find_peaks (spec, height = thresh) 
    #something counts as a peak only if its above the mean-3sigma value 
    #this is because the assumption here is that the majority of the velocity channels
    #are noisy (ie a smaller fraction of the channels have data vs noise. this means that 
    #the mean of the values across the spectral axis is telling you about noise, and 
    #any value that is within the 3sigma of the noise mean is also noise. ie what counts as
    #a velocity channel containing signal is one whose value is > mean + 3sigma)
    #but then we run into a problem when every pixel contains data / the noise has already been clipped away
    print(p)
    peak = p[0] #array of peaks locations
    peak_info = p[1] #array of peak heights
    
    if len(peak) > 1:
        i = np.argmax(peak_info['peak_heights']) #if there are multiple peaks, pick the tallest one #this is problematic though if there is blending
        peak = [peak[i]]
    vel_center = peak[0]

    #this assumes that vel_center is at a place along the spectral axis that there is
    # still data +/- v_channel_width/2 around it ... seems like an ok assumption for now but might 
    # run into problems at some point
    v_channel_start =  vel_center - int(v_channel_width/2)
    v_channel_end = vel_center + int(v_channel_width/2)

    subcube = (clob.cube[v_channel_start:v_channel_end+1, :, :])*clob.units['data']

    #calculate integrated intensity
    W = np.nansum(subcube*clob.dv, axis=0)
    # W = np.nansum(cube*dv, axis=0)

    
    #use W_thresh to create a mask:
    W_mask = W.value > W_thresh
    W_old = W
    W = W_old * W_mask
    
    #calculate column density from the integrated intensity & conversion factor
    N = W * conversion_factor
    
    #define surface density here
    dSig = mmw * (c.m_p.to(u.g)) * N
    
    
    #calculate on sky area & cloud area at each pixel
    dAsky = np.abs(clob.dl * clob.db)
    

    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        kw_imshow = {'cmap': cmap0, 'origin':'lower'}

    
        dSig_Msunpc2 = dSig.to(c.M_sun/u.pc**2)

        fig, ax = plt.subplots (figsize = (15, 10), dpi = 200)
        Sig_plot = ax.imshow(dSig_Msunpc2.value, extent = [clob.lons[0], clob.lons[-1], \
                                                clob.lats[0], clob.lats[-1]], **kw_imshow)


        cbar = fig.colorbar(Sig_plot, ax=ax, label = Sig_label +r' [$M_\odot \ pc^{-2}$]', orientation = 'vertical')

        if data_type == '21cm':
            KMT_thresh = plt.contour(clob.lons, clob.lats, dSig_Msunpc2.value,  linewidths = 0.7, \
            levels=[10, 100], colors = ['xkcd:red', 'xkcd:orange']) #atomic-to-molecular transition for HI maps 
            cbar.add_lines(KMT_thresh)


        ax.set_xlabel(r'Galactic Longitude [$\degree$]')
        ax.set_ylabel(r'Galactic Latitude [$\degree$]')
        ax.set_title(Sig_label)


        plt.show()

    
    return {'surface density per pixel': dSig, 'pixel area': dAsky}


def surface_density (clob, data_type, W_thresh = 1.5, v_channel_width = 25, show_Sigma = True):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from extract_cloud/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    and whether to assume that the spectrum is noisy / will need to isolate velocity peak or not (this is mostly
    for HI clouds)
    3) W_thresh (using value from lewis, lada, dame, but think about whether i wanna change it)

    '''
    spec = dc.spectrum1D(clob.cube)
    plt.plot(spec)
    plt.show()
    
    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/2203.07480.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        thresh = np.nanmin(spec)
        

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        thresh = np.nanmean(spec) + 3 * np.nanstd(spec)
    
        
    p = find_peaks (spec, height = thresh) 
    #something counts as a peak only if its above the mean-3sigma value 
    #this is because the assumption here is that the majority of the velocity channels
    #are noisy (ie a smaller fraction of the channels have data vs noise. this means that 
    #the mean of the values across the spectral axis is telling you about noise, and 
    #any value that is within the 3sigma of the noise mean is also noise. ie what counts as
    #a velocity channel containing signal is one whose value is > mean + 3sigma)
    #but then we run into a problem when every pixel contains data / the noise has already been clipped away
    print(p)
    peak = p[0] #array of peaks locations
    peak_info = p[1] #array of peak heights
    
    if len(peak) > 1:
        i = np.argmax(peak_info['peak_heights']) #if there are multiple peaks, pick the tallest one #this is problematic though if there is blending
        peak = [peak[i]]
    vel_center = peak[0]

    #this assumes that vel_center is at a place along the spectral axis that there is
    # still data +/- v_channel_width/2 around it ... seems like an ok assumption for now but might 
    # run into problems at some point
    v_channel_start =  vel_center - int(v_channel_width/2)
    v_channel_end = vel_center + int(v_channel_width/2)

    subcube = (clob.cube[v_channel_start:v_channel_end+1, :, :])*clob.units['data']

    #calculate integrated intensity
    W = np.nansum(subcube*clob.dv, axis=0)
    # W = np.nansum(cube*dv, axis=0)

    
    #use W_thresh to create a mask:
    W_mask = W.value > W_thresh
    W_old = W
    W = W_old * W_mask
    
    #calculate column density from the integrated intensity & conversion factor
    N = W * conversion_factor
    
    #define surface density here
    dSig = mmw * (c.m_p.to(u.g)) * N
    
    
    #calculate on sky area & cloud area at each pixel
    dAsky = np.abs(clob.dl * clob.db)
    

    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        kw_imshow = {'cmap': cmap0, 'origin':'lower'}

    
        dSig_Msunpc2 = dSig.to(c.M_sun/u.pc**2)

        fig, ax = plt.subplots (figsize = (15, 10), dpi = 200)
        Sig_plot = ax.imshow(dSig_Msunpc2.value, extent = [clob.lons[0], clob.lons[-1], \
                                                clob.lats[0], clob.lats[-1]], **kw_imshow)


        cbar = fig.colorbar(Sig_plot, ax=ax, label = Sig_label +r' [$M_\odot \ pc^{-2}$]', orientation = 'vertical')

        if data_type == '21cm':
            KMT_thresh = plt.contour(clob.lons, clob.lats, dSig_Msunpc2.value,  linewidths = 0.7, \
            levels=[10, 100], colors = ['xkcd:red', 'xkcd:orange']) #atomic-to-molecular transition for HI maps 
            cbar.add_lines(KMT_thresh)


        ax.set_xlabel(r'Galactic Longitude [$\degree$]')
        ax.set_ylabel(r'Galactic Latitude [$\degree$]')
        ax.set_title(Sig_label)


        plt.show()

    
    return {'surface density per pixel': dSig, 'pixel area': dAsky}


def surface_density (clob, data_type, n_sigma = 3, show_Sigma = False):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from create_subcube/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    '''

    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/astro-ph/0411089.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        rms_noise = 0.25 *u.K #mean calculated from values in Table 1 of the paper announcing the full DHT dataset (https://iopscience.iop.org/article/10.1086/318388/pdf)

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        rms_noise = 0.043 *u.K #from Table 1 of the paper announcing the full HI4PI dataset (https://arxiv.org/pdf/1610.06175.pdf)
    
   #assume that the clob has already been processed such that 
   #the noise (values < 3 sigma_rms) in each channel has been suppressed & there is only *one* peak along the velocity axis
   #this works fine for now because working with CO-bright clouds & we can specify what velocity range
   #might need to change for later if CO-dark clouds need this function
   
    #lets assume (at least for now) that the clob 
    #has already been processed through dendrogram.trunk[i].get_mask()

    W = np.nansum(clob.cube * clob.units['data'] * clob.dv, axis=0) #calculate integrated intensity
    N = W * conversion_factor #calculate column density from the integrated intensity & conversion factor
    dSig = mmw * (c.m_p.to(u.g)) * N #define surface density here
    
    mask = (W!=0)
    dAsky = mask * np.ones_like(mask) 
    A_tot_ang = np.nansum(dAsky) * (clob.dl * clob.db)#.to(u.rad**2)
    print(np.nansum(dAsky))
    
    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        
    # return {'surface density per pixel': dSig, 'angular area per pixel': dAsky_signal}
    return dSig, A_tot_ang


def surface_density (clob, data_type, n_sigma = 3, show_Sigma = False):
    '''
    define a function to take in 
    1) a clob (=cloud object)  (the output from create_subcube/unpack_fits_cube)
    2) string that indicates whether the data is CO or 21cm (to decide different which conversion factor to use)
    '''

    if data_type == 'CO':
        conversion_factor = 1.97e20 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #X_CO value from lewis, lada, dame 2022 paper
        mmw = 2 * 1.36 #2 bc molecule, 1.36 = atomic mmw value from https://arxiv.org/pdf/astro-ph/0411089.pdf 
        label = r'H_2'
        Sig_label = r'$\Sigma_{H_2}$'
        rms_noise = 0.25 *u.K #mean calculated from values in Table 1 of the paper announcing the full DHT dataset (https://iopscience.iop.org/article/10.1086/318388/pdf)

    elif data_type == '21cm':
        conversion_factor =  1.82e18 * ((u.cm)**-2)/(u.K*(u.km/u.s)) #from rearranging eq 11 in https://www.mso.anu.edu.au/~krumholz/teaching/astr8007_s1_2021/notes5.pdf
        mmw = 1.36
        Sig_label = r'$\Sigma_{HI}$'
        rms_noise = 0.043 *u.K #from Table 1 of the paper announcing the full HI4PI dataset (https://arxiv.org/pdf/1610.06175.pdf)
    
   #assume that the clob has already been processed such that 
   #the noise (values < 3 sigma_rms) in each channel has been suppressed & there is only *one* peak along the velocity axis
   #this works fine for now because working with CO-bright clouds & we can specify what velocity range
   #might need to change for later if CO-dark clouds need this function
   
    #lets assume (at least for now) that the clob 
    #has already been processed through dendrogram.trunk[i].get_mask()

    W = np.nansum(clob.cube * clob.units['data'] * clob.dv, axis=0) #calculate integrated intensity
    N = W * conversion_factor #calculate column density from the integrated intensity & conversion factor
    dSig = mmw * (c.m_p.to(u.g)) * N #define surface density here
    
    mask = (W!=0)
    dAsky = mask * np.ones_like(W.value) 
    A_tot_ang = np.nansum(dAsky) * (clob.dl * clob.db).to(u.rad**2)

    
    #might also want to figure out a way to more efficiently plot things (ie instead of each function having its own plot routine, 
    #create a function to make the plots, then call it when show_____ is True, .... if possible to generalize the plotting tings)
    #plotting the surface density map
    if show_Sigma:
        cmap0 = dc.cmap_grey0('viridis')
        kw_imshow = {'cmap': cmap0, 'origin':'lower'}

        dSig_Msunpc2 = dSig.to(c.M_sun/u.pc**2)

        fig, ax = plt.subplots (figsize = (15, 10), dpi = 200)
        Sig_plot = ax.imshow(dSig_Msunpc2.value, extent = [clob.lons[0], clob.lons[-1], \
                                                clob.lats[0], clob.lats[-1]], **kw_imshow)


        cbar = fig.colorbar(Sig_plot, ax=ax, label = Sig_label +r' [$M_\odot \ pc^{-2}$]', orientation = 'vertical')

        if data_type == '21cm':
            KMT_thresh = plt.contour(clob.lons, clob.lats, dSig_Msunpc2.value,  linewidths = 0.7, \
            levels=[10, 100], colors = ['xkcd:red', 'xkcd:orange']) #atomic-to-molecular transition for HI maps 
            cbar.add_lines(KMT_thresh)


        ax.set_xlabel(r'Galactic Longitude [$\degree$]')
        ax.set_ylabel(r'Galactic Latitude [$\degree$]')
        ax.set_title(Sig_label)


        plt.show()

    # return {'surface density per pixel': dSig, 'angular area per pixel': dAsky_signal}
    return dSig, A_tot_ang
