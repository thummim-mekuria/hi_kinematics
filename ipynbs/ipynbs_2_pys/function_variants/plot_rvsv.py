#!/usr/bin/env python
# coding: utf-8

def plot_rvsv(i, show=False):  
    l0c,l1c = los.astype('float')[i], lis.astype('float')[i] 
    b0c,b1c = bos.astype('float')[i], bis.astype('float')[i] 

    lc,bc = np.nanmean((los.astype('float')[i],lis.astype('float')[i])), np.nanmean((bos.astype('float')[i],bis.astype('float')[i]))

    l0h,l1h = lc + rhi_px[i]*0.125, lc - rhi_px[i]*0.125
    b0h,b1h = bc - rhi_px[i]*0.125, bc + rhi_px[i]*0.125
    
    hlons = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0]
    hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1]
    
    
    clons = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0]
    clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1]
    

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)

    
    hbperp = (hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)/np.linalg.norm((hidvdx[i],hidvdy[i]))
    cbperp = (h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)/np.linalg.norm((h2dvdx[i],h2dvdy[i]))

    hbperp2 = (hidvdxnew[i]*hxs_pc + hidvdynew[i]*hys_pc)/np.linalg.norm((hidvdxnew[i],hidvdynew[i]))
    cbperp2 = (h2dvdxnew[i]*cxs_pc + h2dvdynew[i]*cys_pc)/np.linalg.norm((h2dvdxnew[i],h2dvdynew[i]))

    
    hv_model =  np.where(~hi_msk[i],np.nan, hi_vos[i] + hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)
    cv_model =  np.where(~h2_msk[i],np.nan, h2_vos[i] + h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)


    hv_model2 =  np.where(~hi_msk[i],np.nan, hi_vos[i] + hidvdxnew[i]*hxs_pc + hidvdynew[i]*hys_pc)
    cv_model2 =  np.where(~h2_msk[i],np.nan, h2_vos[i] + h2dvdxnew[i]*cxs_pc + h2dvdynew[i]*cys_pc)


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0s[i].data),
            np.ravel(hi_1s[i].data),
            np.ravel(hi_2s[i].data),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hbperp2),
            np.ravel(hv_model),
            np.ravel(hv_model2)
            )).transpose(),
            columns=['w','v','s','rper','rper2','vmod','vmod2'])


    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]

    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    meanhdf2 = hdf.groupby(pd.cut(hdf['rper2'],
                        bins=np.arange(np.nanmin(hdf['rper2'])-hR/15,
                        np.nanmax(hdf['rper2']),
                        hR/15))).mean().dropna()

    sighdf2 = hdf.groupby(pd.cut(hdf['rper2'],
                        bins=np.arange(np.nanmin(hdf['rper2'])-hR/15,
                        np.nanmax(hdf['rper2']),
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0s[i].data),
            np.ravel(h2_1s[i].data),
            np.ravel(cbperp),
            np.ravel(cbperp2),
            np.ravel(cv_model),
            np.ravel(cv_model2),
            )).transpose(),
            columns=['w','v','rper','rper2','vmod','vmod2'])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    meancdf2 = cdf.groupby(pd.cut(cdf['rper2'],
                        bins=np.arange(np.nanmin(cdf['rper2'])-cR/15,
                        np.nanmax(cdf['rper2']),
                        cR/15))).mean().dropna()


    sigcdf2 = cdf.groupby(pd.cut(cdf['rper2'],
                        bins=np.arange(np.nanmin(cdf['rper2'])-cR/15,
                        np.nanmax(cdf['rper2']),
                        cR/15))).std().dropna()

    
    if show:
        fig,axs = plt.subplots(1,2,figsize=(12,5))
        # fig.suptitle(names_full[i])
    
        cm,cb = curve_fit(line,meancdf2['rper2'],meancdf2['v'])[0]
        csm   = np.sqrt(curve_fit(line,meancdf2['rper2'],meancdf2['v'])[1][0,0])
    
        hm,hb = curve_fit(line,meanhdf2['rper2'],meanhdf2['v'])[0]
        hs = np.sqrt(curve_fit(line,meanhdf2['rper2'],meanhdf2['v'])[1][0,0])
    
        axs[0].fill_between(meancdf2['rper2'],
                            meancdf2['v']+sigcdf2['v'],
                            meancdf2['v']-sigcdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
        
        axs[0].fill_between(meancdf2['rper2'],
                            meancdf2['v']+2*sigcdf2['v'],
                            meancdf2['v']-2*sigcdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
    
     
        axs[0].plot(meancdf2['rper2'],meancdf2['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[0].plot(meancdf2['rper2'],meancdf2['vmod2'],color='red',linewidth=0.7,linestyle='--')
        
    
        axs[1].fill_between(meanhdf2['rper2'],
                            meanhdf2['v']+sighdf2['v'],
                            meanhdf2['v']-sighdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
        
        axs[1].fill_between(meanhdf2['rper2'],
                            meanhdf2['v']+2*sighdf2['v'],
                            meanhdf2['v']-2*sighdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
    
        
        axs[1].plot(meanhdf2['rper2'],meanhdf2['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[1].plot(meanhdf2['rper2'],meanhdf2['vmod2'],color='red',linewidth=0.7,linestyle='--')
    
    
        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
            # ax.grid()
        axs[0].set_ylabel(r'velocity [$\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names_full[i])
    
        # fig.savefig('v_vs_r_'+names[i],bbox_inches='tight')
        plt.show()

    return pearson(meancdf2['rper2'], meancdf2['v']),                   pearson(meanhdf2['rper2'], meanhdf2['v']),\
           pearson(sigcdf2['v'],      np.abs(meancdf2['v']-h2_vos[i])), pearson(sighdf2['v'],      np.abs(meanhdf2['v']-hi_vos[i])),\
           pearson(meancdf2['w'],     np.abs(meancdf2['v']-h2_vos[i])), pearson(meanhdf2['w'],     np.abs(meanhdf2['v']-hi_vos[i]))


def plot_rvsv(i,rhi,vovi, show=False):  
    h2dvds, hidvds = Omg_vec(i,pixelize(rhi_sum[i],ds[i]),vsnew[i],xy=True)/1e3
    
    omgh2 = np.linalg.norm(h2dvds[:2])
    omghi = np.linalg.norm(hidvds[:2])
    thth2 = get_theta(h2dvds[:2])
    ththi = get_theta(hidvds[:2])
    
    l0c,l1c,b0c,b1c = coordses[i][:4]
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2maps, himaps = moment_maps(i,threshs[i],rhi,vovi)
    h2_0,h2_1,h2_2 = h2maps
    hi_0,hi_1,hi_2 = himaps
    
    l0h,l1h = lc + rhi*0.125, lc - rhi*0.125
    b0h,b1h = bc - rhi*0.125, bc + rhi*0.125
    
    
    hlons = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0]
    hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1]
    
    
    clons = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0]
    clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1]
    

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)

    
    hbperp = (hidvds[0]*hxs_pc + hidvds[1]*hys_pc)/np.linalg.norm(hidvds)
    cbperp = (h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc)/np.linalg.norm(h2dvds)

    
    hv_model =  hidvds[2] + hidvds[0]*hxs_pc + hidvds[1]*hys_pc
    cv_model =  h2dvds[2] + h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hi_2),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod',])


    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]

    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0),
            np.ravel(h2_1),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','rper','vmod',])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    h2p, hip = pearson(meancdf['rper'],  meancdf['v']),pearson(meanhdf['rper'],  meanhdf['v'])


def plot_rvsv(i, show=False):  
    ccoords = coordses[i][:4]
    d = ds[i]
    omgs = tab2_newnew[['omgh2','omghi']].to_numpy()[i]/1e3
    thts = tab2_newnew[['thth2','ththi']].to_numpy()[i]
    
    l0c,l1c,b0c,b1c = ccoords
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2maps, himaps = moment_maps(i,threshs[i])
    h2_0,h2_1,h2_2 = h2maps
    hi_0,hi_1,hi_2 = himaps
    
    rhi_px = hi_0.shape[0]/2

    l0h,l1h = lc + rhi_px*0.125, lc - rhi_px*0.125
    b0h,b1h = bc - rhi_px*0.125, bc + rhi_px*0.125
    
    
    hlons = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0]
    hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1]
    
    
    clons = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0]
    clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1]
    

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)

    
    hbperp = (hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)/np.linalg.norm((hidvdx[i],hidvdy[i]))
    cbperp = (h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)/np.linalg.norm((h2dvdx[i],h2dvdy[i]))

    
    hv_model =  hi_vos[i] + hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc
    cv_model =  h2_vos[i] + h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hi_2),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod',])


    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]

    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0),
            np.ravel(h2_1),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','rper','vmod',])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    h2p, hip = pearson(meancdf['rper'],  meancdf['v']),pearson(meanhdf['rper'],  meanhdf['v'])


def plot_rvsv(i, show=False):  
    l0c,l1c = los.astype('float')[i], lis.astype('float')[i] 
    b0c,b1c = bos.astype('float')[i], bis.astype('float')[i] 

    lc,bc = np.nanmean((los.astype('float')[i],lis.astype('float')[i])), np.nanmean((bos.astype('float')[i],bis.astype('float')[i]))

    l0h,l1h = lc + rhi_px[i]*0.125, lc - rhi_px[i]*0.125
    b0h,b1h = bc - rhi_px[i]*0.125, bc + rhi_px[i]*0.125
    
    hlons = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0]
    hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1]
    
    
    clons = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0]
    clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1]
    

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)

    
    hbperp = (hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)/np.linalg.norm((hidvdx[i],hidvdy[i]))
    cbperp = (h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)/np.linalg.norm((h2dvdx[i],h2dvdy[i]))

    hbperp2 = (hidvdxnew[i]*hxs_pc + hidvdynew[i]*hys_pc)/np.linalg.norm((hidvdxnew[i],hidvdynew[i]))
    cbperp2 = (h2dvdxnew[i]*cxs_pc + h2dvdynew[i]*cys_pc)/np.linalg.norm((h2dvdxnew[i],h2dvdynew[i]))

    
    hv_model =  np.where(~hi_msk[i],np.nan, hi_vos[i] + hidvdx[i]*hxs_pc + hidvdy[i]*hys_pc)
    cv_model =  np.where(~h2_msk[i],np.nan, h2_vos[i] + h2dvdx[i]*cxs_pc + h2dvdy[i]*cys_pc)


    hv_model2 =  np.where(~hi_msk[i],np.nan, hi_vos[i] + hidvdxnew[i]*hxs_pc + hidvdynew[i]*hys_pc)
    cv_model2 =  np.where(~h2_msk[i],np.nan, h2_vos[i] + h2dvdxnew[i]*cxs_pc + h2dvdynew[i]*cys_pc)


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0s[i].data),
            np.ravel(hi_1s[i].data),
            np.ravel(hi_2s[i].data),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hbperp2),
            np.ravel(hv_model),
            np.ravel(hv_model2)
            )).transpose(),
            columns=['w','v','s','rper','rper2','vmod','vmod2'])


    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]

    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper']),
                        hR/15))).std().dropna()


    meanhdf2 = hdf.groupby(pd.cut(hdf['rper2'],
                        bins=np.arange(np.nanmin(hdf['rper2'])-hR/15,
                        np.nanmax(hdf['rper2']),
                        hR/15))).mean().dropna()

    sighdf2 = hdf.groupby(pd.cut(hdf['rper2'],
                        bins=np.arange(np.nanmin(hdf['rper2'])-hR/15,
                        np.nanmax(hdf['rper2']),
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0s[i].data),
            np.ravel(h2_1s[i].data),
            np.ravel(cbperp),
            np.ravel(cbperp2),
            np.ravel(cv_model),
            np.ravel(cv_model2),
            )).transpose(),
            columns=['w','v','rper','rper2','vmod','vmod2'])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper']),
                        cR/15))).std().dropna()


    meancdf2 = cdf.groupby(pd.cut(cdf['rper2'],
                        bins=np.arange(np.nanmin(cdf['rper2'])-cR/15,
                        np.nanmax(cdf['rper2']),
                        cR/15))).mean().dropna()


    sigcdf2 = cdf.groupby(pd.cut(cdf['rper2'],
                        bins=np.arange(np.nanmin(cdf['rper2'])-cR/15,
                        np.nanmax(cdf['rper2']),
                        cR/15))).std().dropna()

    
    if show:
        fig,axs = plt.subplots(1,2,figsize=(12,5))
        # fig.suptitle(names_full[i])
    
        cm,cb = curve_fit(line,meancdf2['rper2'],meancdf2['v'])[0]
        csm   = np.sqrt(curve_fit(line,meancdf2['rper2'],meancdf2['v'])[1][0,0])
    
        hm,hb = curve_fit(line,meanhdf2['rper2'],meanhdf2['v'])[0]
        hs = np.sqrt(curve_fit(line,meanhdf2['rper2'],meanhdf2['v'])[1][0,0])
    
        axs[0].fill_between(meancdf2['rper2'],
                            meancdf2['v']+sigcdf2['v'],
                            meancdf2['v']-sigcdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
        
        axs[0].fill_between(meancdf2['rper2'],
                            meancdf2['v']+2*sigcdf2['v'],
                            meancdf2['v']-2*sigcdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
    
     
        axs[0].plot(meancdf2['rper2'],meancdf2['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[0].plot(meancdf2['rper2'],meancdf2['vmod2'],color='red',linewidth=0.7,linestyle='--')
        
    
        axs[1].fill_between(meanhdf2['rper2'],
                            meanhdf2['v']+sighdf2['v'],
                            meanhdf2['v']-sighdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
        
        axs[1].fill_between(meanhdf2['rper2'],
                            meanhdf2['v']+2*sighdf2['v'],
                            meanhdf2['v']-2*sighdf2['v'],
                            color=(0.7,0.7,0.7,0.1),linewidth=0.5)
    
        
        axs[1].plot(meanhdf2['rper2'],meanhdf2['v'],color='black',marker='o',linewidth=0.5,ms=2)
        axs[1].plot(meanhdf2['rper2'],meanhdf2['vmod2'],color='red',linewidth=0.7,linestyle='--')
    
    
        for ax in axs.flat:
            ax.set_xlabel(r'$\Delta r_\perp$ [pc]')
            # ax.grid()
        axs[0].set_ylabel(r'velocity [$\mathrm{km\ s^{-1}}$]')
        fig.suptitle(names_full[i])
    
        # fig.savefig('v_vs_r_'+names[i],bbox_inches='tight')
        plt.show()

    return pearson(meancdf2['rper2'], meancdf2['v']),                   pearson(meanhdf2['rper2'], meanhdf2['v']),\
           pearson(sigcdf2['v'],      np.abs(meancdf2['v']-h2_vos[i])), pearson(sighdf2['v'],      np.abs(meanhdf2['v']-hi_vos[i])),\
           pearson(meancdf2['w'],     np.abs(meancdf2['v']-h2_vos[i])), pearson(meanhdf2['w'],     np.abs(meanhdf2['v']-hi_vos[i]))


def plot_rvsv(i,rhi,vovi, show=False):  
    h2dvds, hidvds = Omg_vec(i,rhi,vovi,xy=True)/1e3
    
    omgh2 = np.linalg.norm(h2dvds[:2])
    omghi = np.linalg.norm(hidvds[:2])
    
    
    thth2 = get_theta(h2dvds[:2])
    ththi = get_theta(hidvds[:2])
    
    l0c,l1c,b0c,b1c = coordses[i][:4]
    lc,bc = np.nanmean((l0c,l1c)), np.nanmean((b0c,b1c))
    
    h2maps, himaps = moment_maps(i,threshs[i],rhi,vovi)
    h2_0,h2_1,h2_2 = h2maps
    hi_0,hi_1,hi_2 = himaps
    
    l0h,l1h = lc + rhi*0.125, lc - rhi*0.125
    b0h,b1h = bc - rhi*0.125, bc + rhi*0.125
    
    
    hlons = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[0]
    hlats = np.meshgrid(np.arange(l0h,l1h,-0.125),np.arange(b0h,b1h,0.125))[1]
    
    
    clons = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[0]
    clats = np.meshgrid(np.arange(l0c,l1c,-0.125),np.arange(b0c,b1c,0.125))[1]
    

    hxs_pc = (hlons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    hys_pc = (hlats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    cxs_pc = (clons-lc)*u.deg.to(u.rad)*ds.astype('float')[i]
    cys_pc = (clats-bc)*u.deg.to(u.rad)*ds.astype('float')[i]

    hrs_pc = np.linalg.norm(np.meshgrid(hxs_pc[i],hys_pc[i]),axis=0)
    crs_pc = np.linalg.norm(np.meshgrid(cxs_pc[i],cys_pc[i]),axis=0)

    
    hbperp = (hidvds[0]*hxs_pc + hidvds[1]*hys_pc)/omghi
    cbperp = (h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc)/omgh2

    
    hv_model =  hidvds[2] + hidvds[0]*hxs_pc + hidvds[1]*hys_pc
    cv_model =  h2dvds[2] + h2dvds[0]*cxs_pc + h2dvds[1]*cys_pc


    hdf = pd.DataFrame(np.array((
            np.ravel(hi_0),
            np.ravel(hi_1),
            np.ravel(hi_2),
            # np.ravel(hrs_pc),
            np.ravel(hbperp),
            np.ravel(hv_model),
            )).transpose(),
            columns=['w','v','s','rper','vmod',])


    hdf_noisy = hdf
    hdf = hdf_noisy[hdf_noisy['w']>4e20]

    hR = np.nanmax(hdf['rper'])
    meanhdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper'])+hR/15,
                        hR/15))).mean().dropna()

    sighdf = hdf.groupby(pd.cut(hdf['rper'],
                        bins=np.arange(np.nanmin(hdf['rper'])-hR/15,
                        np.nanmax(hdf['rper'])+hR/15,
                        hR/15))).std().dropna()


    cdf = pd.DataFrame(np.array((
            np.ravel(h2_0),
            np.ravel(h2_1),
            np.ravel(cbperp),
            np.ravel(cv_model),
            )).transpose(),
            columns=['w','v','rper','vmod',])

    cdf_noisy = cdf
    cdf = cdf_noisy[cdf_noisy['w']>4e20]

    cR = np.nanmax(cdf['rper'])
    meancdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper'])+cR/15,
                        cR/15))).mean().dropna()


    sigcdf = cdf.groupby(pd.cut(cdf['rper'],
                        bins=np.arange(np.nanmin(cdf['rper'])-cR/15,
                        np.nanmax(cdf['rper'])+cR/15,
                        cR/15))).std().dropna()


    h2p, hip = pearson(meancdf['rper'],  meancdf['v']),pearson(meanhdf['rper'],  meanhdf['v'])
