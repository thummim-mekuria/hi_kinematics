#!/usr/bin/env python
# coding: utf-8


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import datacube as dc
import pandas as pd
from astropy.table import Table
from astropy.modeling import models
from scipy.signal import find_peaks


hi4pi,dht = fits.open('both_same_grid_720d.fits')


vels = np.arange(155)*1.3 - 100


dht18030 = dht.data[:,160:640,1440:4320]


plt.imshow(np.nansum(dht18030,0),cmap=dc.cmap_min0('viridis'))


mom0 = np.nansum(dht18030,0)


np.nansum(mom0 > 1)


cospecs = dht18030.T.reshape(-1,155)


lb_pxs = np.array(np.meshgrid(np.arange(2880)*-0.125 + 180,np.arange(480)*0.125 - 30)).T.reshape(-1,2)


goodspecs = [np.nansum(cosp) > 1 for cosp in cospecs]


co_spectra = Table([lb_pxs[goodspecs],cospecs[goodspecs]],names=['coord','spec'])


raw_ffts  = [np.fft.rfftn(spec) for spec in cospecs[goodspecs]]
fft_mags  = [np.abs(fft)**2 for fft in raw_ffts]
fft_norm0 = [fft_mag[0] for fft_mag in fft_mags]


co_spectra['raw_fft'] = raw_ffts
co_spectra['fft_amp'] = fft_mags
co_spectra['k0_T'] = fft_norm0


ks = np.fft.rfftfreq(155)


co_spectra['coord'].T[0]


plt.figure(figsize=(24,8))
plt.scatter(co_spectra['coord'].T[0],co_spectra['coord'].T[1],s=1,
            c=np.nansum(co_spectra['spec'],axis=1),norm='asinh')
plt.colorbar()
plt.xlim(180,-180)


gaus2 = models.Gaussian1D(1.,0.,2.)
gaus5 = models.Gaussian1D(1.,0.,5.)
gaus10 = models.Gaussian1D(1.,0.,10.)
gaus20 = models.Gaussian1D(1.,0.,20.)
gaus50 = models.Gaussian1D(1.,0.,50.)


# plt.plot(gaus20(vels))
plt.plot(ks,np.fft.fftshift(np.abs(np.fft.fftn(gaus2(vels)+1))**2)[77:]/14) #normalizing values calculated by hand
plt.plot(ks,np.fft.fftshift(np.abs(np.fft.fftn(gaus5(vels)+1))**2)[77:]/90)
plt.plot(ks,np.fft.fftshift(np.abs(np.fft.fftn(gaus10(vels)+1))**2)[77:]/337)
plt.plot(ks,np.fft.fftshift(np.abs(np.fft.fftn(gaus20(vels)+1))**2)[77:]/1007)
plt.plot(ks,np.fft.fftshift(np.abs(np.fft.fftn(gaus50(vels)+1))**2)[77:]/1009)

plt.vlines(0.12,1e-5,1e5,linestyle='--',color='grey')
plt.hlines(0.05,0,0.5,linestyle='--',color='grey')
plt.hlines(0.26,0,0.5,linestyle='--',color='grey')

# plt.plot(np.abs(np.fft.rfftn(gaus10(vels)+1))**2)
# plt.yscale('log')
plt.loglog()
plt.ylim(1e-3,1e0)


gaus10_fft = np.fft.fftshift(np.abs(np.fft.fftn(gaus5(vels)+1))**2)


plt.plot(np.fft.fftshift(np.abs(np.fft.fftn(gaus10_fft))**2)/1e8 - 7)
plt.plot(gaus10(vels))


plt.plot(gaus10_fft)


plt.plot(cospecs[258022])


cmplx_fft = np.fft.fftshift(np.abs(np.fft.fftn(cospecs[258022]))**2)


cmplx_fft2 = np.fft.fftshift(np.abs(np.fft.fftn(cmplx_fft))**2)


plt.plot((np.abs(np.fft.fftn(cmplx_fft))**2))
# plt.plot(cmplx_fft[77:])
# plt.loglog()


vcmplx = cospecs[157200]


def iter_fft(array):
    ft = np.abs(
           np.fft.fftn(
                       array))
    
#     non_redund_ft = ft[:int(0.5*len(ft))]
      
    return ft


from scipy.ndimage import gaussian_filter1d


gaussian_filter1d(iter_fft(vcmplx),sigma=10)


fft = iter_fft(vcmplx)


fft = np.fft.fftn(vcmplx)


smd = gaussian_filter1d(np.abs(fft),sigma=7)


len(np.abs(np.fft.ifftn(smd)))


plt.plot(vcmplx)
plt.plot(np.fft.ifft(fft-smd))


plt.plot(np.fft.fftshift(np.fft.fftn(gaus10(vels))))
plt.plot(np.fft.fftshift(np.fft.fftn(vcmplx)))
# plt.xlim(10,80)
# plt.yscale


plt.plot(fft)
plt.plot(np.fft.ifft(fft-smd))
plt.plot(smd)


plt.plot(iter_fft(vcmplx))
plt.plot(gaussian_filter1d(iter_fft(vcmplx),sigma=7))
# plt.loglog()


def count_ffts(array)
    n=1
    ft1 = iter_fft(array)
        
    arrmin = np.nanmin(ft1)
    arrstd = np.nanstd(ft1)
    peak_num = find_peaks(ft1,arrmin+arrstd)[0]
        
    while len(peak_num) > 1:
        n += 1
        ft = iter_fft(ft)
        
#     else:
    print(n)


iter_fft(vcmplx)


np.nanmin(iter_fft(vcmplx)),np.nanstd(iter_fft(vcmplx))


msk = np.full(len(iter_fft(vcmplx)),0)
msk[72:82] = np.full(10,1)


plt.plot(np.fft.fftshift(np.fft.fftn(vcmplx))*msk)
# plt.plot(np.fft.fftshift(iter_fft(vcmplx)))


co_spectra


mskd_cmplz = np.abs(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(vcmplx))*msk)))


np.fft.fftshift(np.abs(np.fft.ifftn(vcmplx)))*155


np.fft.fftshift(np.abs(np.fft.fftn(vcmplx)))


plt.plot(np.fft.fftshift(np.abs(np.fft.ifftn(vcmplx)))*155)
plt.plot(np.fft.fftshift(np.abs(np.fft.fftn(vcmplx))))

# plt.plot(vcmplx)


plt.plot(vcmplx-mskd_cmplz)
plt.plot(vcmplx)
# plt.ylim(0)


len(find_peaks(iter_fft(vcmplx))[0])


find_peaks(iter_fft(vcmplx),height=np.nanmin(iter_fft(vcmplx))+
          np.nanstd(iter_fft(vcmplx)))[0]


np.nanmin(iter_fft2(gaus10_fft)) + np.nanstd(iter_fft2(gaus10_fft))


0.1*np.max(iter_fft2(gaus10_fft))


iter_fft2(gaus10_fft)


help(find_peaks)


find_peaks(iter_fft2(gaus10_fft))[0]


np.nanstd(iter_fft2(gaus10_fft))


plt.plot(iter_fft2(gaus10_fft))


plt.plot(iter_fft2(gaus10_fft))
plt.yscale('log')


1382400/2


[(plt.plot(ks,co_spectra['fft_amp'][i]/co_spectra['k0_T'][i],alpha=0.8),
  print(i),plt.loglog()) 
 for i in range(0,len(co_spectra),68217)]


plt.plot(vels,(models.Gaussian1D(1.,-1.8,2.1) + models.Gaussian1D(1.,1.5,2.))(vels))

# plt.plot(vels,(models.Gaussian1D(1.4,-1.,1.) + models.Gaussian1D(1.,1.,1.))(vels))
plt.plot(vels,models.Gaussian1D(1.46,0.,2.3)(vels))
# plt.plot(vels,models.Lorentz1D(1.2,0.,4.)(vels))


plt.xlim(-15,15)
