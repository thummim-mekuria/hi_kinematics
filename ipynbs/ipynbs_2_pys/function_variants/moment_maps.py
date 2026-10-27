#!/usr/bin/env python
# coding: utf-8

def moment_maps(i,thresh, rhi, vovi, hidat, h2dat):
    h2cube,hicube = make_h2_hi_cubes(i,rhi,vovi,hidat,h2dat)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.29*1.82e18
    
    vo,vi = vovi

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.29-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.29-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.log10(hilb)>20)
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


def moment_maps(i,thresh, rhi, vovi):
    h2cube,hicube = make_h2_hi_cubes(i,rhi, vovi)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.29*1.82e18
    
    vo,vi = vovi

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.29-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.29-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.log10(hilb)>20) #*\
             # (np.linalg.norm(np.indices(hi1.shape) - \
             #                np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi)


def moment_maps(i):
    
    h2cube,hicube = make_h2_hi_cubes(i)
    vo,vi = vsnew[i]
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(~h2_msks[i],np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            np.log10(hilb)>20
                
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


def moment_maps(i):
    
    h2cube,hicube = make_h2_hi_cubes(i)
    vo,vi = vsnew[i]
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(~h2_msks[i],np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            np.log10(hilb)>20
                
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


def moment_maps(i,thresh):
    vo,vi = vsnew[i]
    h2cube,hicube = make_h2_hi_cubes(i)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18
    

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            (np.log10(hilb)>20)


def moment_maps(i,thresh):
    
    h2cube,hicube = make_h2_hi_cubes(i)
    vo,vi = vsnew[i]
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.3*1.82e18

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.3-100)[0]
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.linalg.norm(np.indices(hi1.shape) - \
                            np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi_px[i])*\
                            np.log10(hilb)>thresh
                
    
    himaps = np.array([np.where(~hi_msk,np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


def moment_maps(i,sig=2):
    coords = coordses[i] 
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')

    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = vses[i][sig-1]#0 for 1sigma, 1 for 2 sigma .. 

    h2lb = np.log10(np.nansum(get_subcube(dht.data,coords[:4]),axis=0)*1.3*1.97e20)
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)

    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    h22 = np.sqrt(second_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0])
    
    
    h2maps = np.array([np.where(~h2_msks[i],np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = np.linalg.norm(np.indices(hi1.shape) - \
                   np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) 
    
    himaps = np.array([np.where(hi_msk > rhi_px[i],np.nan,himap) for himap in [hilb,hi1,hi2]])
    
    return h2maps, himaps


def moment_maps(i,thresh, rhi, vovi):
    h2cube,hicube = make_h2_hi_cubes(i,rhi, vovi)
    
    h2lb = np.nansum(h2cube,axis=0)*1.3*1.97e20
    hilb = np.nansum(hicube,axis=0)*1.29*1.82e18
    
    vo,vi = vovi

    h21 = first_moment(h2cube,np.arange(155)*1.3-100)[0]
    hi1 = first_moment(hicube,np.arange(vo,vi)*1.29-100)[0]
    
    
    h22 = np.sqrt(second_moment(h2cube,np.arange(155)*1.3-100)[0])
    hi2 = np.sqrt(second_moment(hicube,np.arange(vo,vi)*1.29-100)[0])
    
    
    h2maps = np.array([np.where(np.log10(h2lb)<thresh,np.nan,h2map) for h2map in [h2lb,h21,h22]])
    
    hi_msk = (np.log10(hilb)>20) #*\
             # (np.linalg.norm(np.indices(hi1.shape) - \
             #                np.array(hi1.shape)[:,np.newaxis,np.newaxis]/2,axis=0) < rhi)
