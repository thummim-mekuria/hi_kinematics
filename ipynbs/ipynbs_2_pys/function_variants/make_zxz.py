#!/usr/bin/env python
# coding: utf-8

def make_zxz(i,save=False):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((coordses[i][4:]+100)/1.3).astype('int')
    
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)
    
    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2lb, extent=coords[:4],
                   levels=np.linspace(h2_Nlims[i][0]+0.1,h2_Nlims[i][1],25),
                   cmap='YlGn_r',extend='max')

        
    im2 = axs[0,1].contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="YlGn_r",
                    extend='both',zorder=-10
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2lb<threshs[i]+0.1,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='bwr_r',extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='bwr_r',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.88,antialiased=True,
                       )  
        
    for ax in axs.transpose()[1]:
        ax.contour(h2lb, extent=coords[:4],
                   levels=np.linspace(threshs[i]+0.1,22,5),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hilb.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],
                    levels=[rhi_px[i],2*rhi_px[i]],colors='white')
        
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2,im3,im4],
                           axs.flat,
                          ['',r'log (N / $\mathrm{cm^{-2}}$)','',
                           r'$\overline{v_\mathrm{LSR}}\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


def make_zxz(i,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i)
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=35,vmin=5,#np.linspace(15,150,35),
                   cmap='magma',extend='max')

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',)#extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(5,20,50),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.9
                       )  
        
    for ax in axs.flat:#.transpose()[1]:
        ax.contour(h2Sig.value, extent=coords[:4],
                   levels=[1,5,10,30],#np.linspace(15,30,4),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hi0.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],
                    levels=[rhi_px[i],2*rhi_px[i]],colors='white')
        
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2,im3,im4],
                           axs.flat,
                          ['',r'$\Sigma \ [\mathrm{M_\odot\ pc^{-2}}]$','',
                           r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


def make_zxz(i,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i)
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=35,vmin=5,#np.linspace(15,150,35),
                   cmap='magma',extend='max')

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',)#extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(5,20,50),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.9
                       )  
        
    for ax in axs.flat:#.transpose()[1]:
        ax.contour(h2Sig.value, extent=coords[:4],
                   levels=[1,5,10,30],#np.linspace(15,30,4),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hi0.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],
                    levels=[rhi_px[i],2*rhi_px[i]],colors='white')
        
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2,im3,im4],
                           axs.flat,
                          ['',r'$\Sigma \ [\mathrm{M_\odot\ pc^{-2}}]$','',
                           r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


def make_zxz(i,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i,threshs[i])
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(13,13))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=np.linspace(10,np.nanmax(h2Sig.value),35),
                   cmap='magma',extend='max',norm='asinh')

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                    ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)
    
    
    h2ox, h2oy, h2vo = fit_plane2maps(i,0)
    hiox, hioy, hivo = fit_plane2maps(i,1)


def make_zxz(i,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i)
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=35,vmin=5,#np.linspace(15,150,35),
                   cmap='magma',extend='max')

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',)#extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=35,#np.linspace(5,20,50),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.9
                       )  
        
    for ax in axs.flat:#.transpose()[1]:
        ax.contour(h2Sig.value, extent=coords[:4],
                   levels=[1,5,10,30],#np.linspace(15,30,4),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hi0.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],
                    levels=[rhi_px[i],2*rhi_px[i]],colors='white')
        
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2,im3,im4],
                           axs.flat,
                          ['',r'$\Sigma \ [\mathrm{M_\odot\ pc^{-2}}]$','',
                           r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


def make_zxz(i,save=False):
    coords = coordses[i]
    
    lc_p = ((lc[i]-360)/-0.125).astype('int')
    bc_p = ((bc[i]+50)/0.125).astype('int')
    
    hicube = hi4pi.data[:, bc_p-rhi_px[i]:bc_p+rhi_px[i], lc_p-rhi_px[i]:lc_p+rhi_px[i],]
    vo,vi = ((vs2[i]+100)/1.3).astype('int')
    
    h2lb = np.log10(make_lb_lv_bv(get_subcube(dht.data,coords[:4]))[0])
    hilb = np.log10(np.nansum(hicube[vo:vi],axis=0)*1.3*1.82e18)
    
    h21 = first_moment(get_subcube(dht.data,coords[:4])[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    hi1 = first_moment(hicube[vo:vi],np.arange(vo,vi)*1.3-100)[0]
    
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi_px[i]
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi_px[i]
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(10,10))
    im1 = axs[0,0].contourf(h2lb, extent=coords[:4],
                   levels=np.linspace(h2_Nlims[i][0]+0.1,h2_Nlims[i][1],25),
                   cmap='YlGn_r',extend='max')

        
    im2 = axs[0,1].contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="YlGn_r",
                    extend='both',zorder=-10
                       ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2lb<threshs[i]+0.1,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='bwr_r',extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='bwr_r',extend='both',zorder=-10)


    for ax in axs.transpose()[0]:
        ax.contourf(hilb,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Nlims[i][0],hi_Nlims[i][1],40),
                    cmap="binary_r",
                    extend='both',zorder=-9,
                    alpha=0.88,antialiased=True,
                       )  
        
    for ax in axs.transpose()[1]:
        ax.contour(h2lb, extent=coords[:4],
                   levels=np.linspace(threshs[i]+0.1,22,5),
                   colors='black',
                   linewidths=0.75,)
    
    for ax in axs.flat:
        hix,hiy = np.indices(hilb.shape)
        hicrc = np.sqrt((hix-rhi_px[i])**2 + (hiy-rhi_px[i])**2)
        
        xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
        ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
        
        ax.contourf(hicrc,extent=[loh,lih,boh,bih],
                    levels=[rhi_px[i],2*rhi_px[i]],colors='white')
        
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        ax.label_outer()
        ax.set_frame_on(False)
    
    fig.suptitle(names[i])
    
    
    for im,ax,lbl in zip([im1,im2,im3,im4],
                           axs.flat,
                          ['',r'log (N / $\mathrm{cm^{-2}}$)','',
                           r'$\overline{v_\mathrm{LSR}}\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl)
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()


def make_zxz(i,rhi,vovi,save=False):
    coords = coordses[i]
    h2maps,himaps = moment_maps(i,threshs[i],rhi,vovi)
    h20,h21,h22 = h2maps
    hi0,hi1,hi2 = himaps
    
    loh,lih = lc[i]+np.array([1,-1])*0.125*rhi
    boh,bih = bc[i]+np.array([-1,1])*0.125*rhi
    
    xtcks = np.around(np.linspace(loh-0.6,lih+0.6,3),1)
    ytcks = np.around(np.linspace(boh+0.6,bih-0.6,3),1)
    
    
    h2dvds, hidvds = Omg_vec(i,rhi,vovi,xy=True)/1e3
    
    h2ox, h2oy, h2vo = h2dvds
    hiox, hioy, hivo = hidvds
    
    
    h2Sig = ((h20/u.cm**2)*c.m_p*2*1.36).to(u.solMass/u.pc**2)
    hiSig = ((hi0/u.cm**2)*c.m_p*1.36).to(u.solMass/u.pc**2)
    
    
    fig,axs = plt.subplots(2,2,dpi=200,figsize=(13,13))
    im1 = axs[0,0].contourf(h2Sig.value, extent=coords[:4],
                   levels=np.linspace(5,3*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),35),
                   cmap='YlGnBu_r',extend='max',)

        
    im2 = axs[0,1].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),
                    cmap="magma",
                    extend='both',zorder=-10,
                    ) 
    
    
    im3 = axs[1,0].contourf(np.where(h2Sig.value<10,np.nan,h21), extent=coords[:4],
                    levels=np.linspace(h2_vlims[i][0],h2_vlims[i][1],40),
                    cmap='Spectral',extend='both')
    
    im4 = axs[1,1].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),#np.linspace(-3,5,50),
                    cmap='Spectral',extend='both',zorder=-10)
    
    
    axs[0,0].contourf(hiSig.value,extent=[loh,lih,boh,bih],
                levels=np.linspace(hi_Slims[i][0],hi_Slims[i][1],40),#np.linspace(5,20,50),
                cmap="binary_r",
                extend='both',zorder=-9,
                alpha=0.9
                   )
    axs[1,0].contourf(hi1,extent=[loh,lih,boh,bih],
                    levels=np.linspace(hi_vlims[i][0],hi_vlims[i][1],40),
                    cmap='binary_r',extend='both',zorder=-10)
        
    for ax in axs.flat:#.transpose()[1]:
        ax.contour(h2Sig.value, extent=coords[:4],
                   levels=np.array([1,2,3])*0.25*np.nanmax(h2Sig.value),
                   # np.array([-1,0,1])*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),#np.linspace(2,50,5),
                   colors='black',
                   linestyles=[':','solid','--'],
                   linewidths=1.2,)
    
                
        ax.set_aspect(1)
        ax.set_xlim(loh,lih)
    
        ax.set_xticks(xtcks,[str(xtc)[:-2]+r"$^\degree$"+str(xtc)[-2:] for xtc in xtcks])
        ax.set_yticks(ytcks,[str(ytc)[:-2]+r"$^\degree$"+str(ytc)[-2:] for ytc in ytcks])
        ax.tick_params(axis='both', which='both', length=7.5)
        
        # ax.label_outer()
        ax.set_frame_on(False)
        ax.tick_params(which='both',top=True,bottom=True, left=True, right=True)

    
    for im,ax,lbl in zip([im1,im2],
                           axs.flat[:2],
                          ['$\Sigma_\mathrm{H_2} \ [\mathrm{M_\odot\ pc^{-2}}]$',r'$\Sigma_\mathrm{HI} \ [\mathrm{M_\odot\ pc^{-2}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl,norm='asinh')
        
        
    for im,ax,lbl in zip([im3,im4],
                           axs.flat[2:],
                          ['',r'$\langle v_\mathrm{LSR} \rangle\ [\mathrm{km\ s^{-1}}]$']):
        
        tcks = np.around(np.linspace(im.levels[0],im.levels[-1],3),1)
        
        fig.colorbar(im,ax=ax,shrink=0.5,ticks=tcks,label=lbl,norm='asinh')
        
        
    dx = (loh-lih)*0.1
    # print(dx)
    
    axs[1,0].add_patch(Circle((lih+dx,boh+dx),dx*0.8,color='white'))
    axs[1,0].quiver(lih+dx,boh+dx,-h2ox,h2oy,color=(0.,0.8,0.8),headwidth=3,zorder=30,pivot='tail',alpha=0.8)
    axs[1,0].quiver(lih+dx,boh+dx,-hiox,hioy,color=(0.8,0.,0.),headwidth=3,zorder=30,pivot='tail',alpha=0.8)
    
    axs[0,1].annotate(str(int(np.nanmax(h2Sig.value))),(loh-0.5*dx,boh+0.5*dx),color=(0.,0.8,0.8),fontsize='x-small')
    axs[0,1].contour(h2Sig.value, extent=coords[:4],
                   levels=np.array([1,2,3])*0.25*np.nanmax(h2Sig.value),
                   # np.array([-1,0,1])*np.nanstd(h2Sig.value)+np.nanmean(h2Sig.value),#np.linspace(2,50,5),
                   colors=(0.,0.8,0.8),
                   linestyles=[':','solid','--'],
                   linewidths=1.2,)
        
        
    fig.suptitle(names[i])
        
    
    fig.supxlabel('Galactic Longitude')
    fig.supylabel('Galactic Latitude')
    
    fig.tight_layout()
    if save:
        plt.savefig(names[i][:3]+names[i][-1]+'_zxz',bbox_inches='tight')
    plt.show()
