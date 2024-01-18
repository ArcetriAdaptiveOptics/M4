def PSD_2D(filename, image, DIA, var_name):
    # the function compute the 2D PSD of an SQUARE image
    # based on:
    # https://bertvandenbroucke.netlify.app/2019/05/24/computing-a-power-spectrum-in-python/    
    #------------------------------------------------------------------------------    
    import numpy as np
    import scipy.stats as stats
    import matplotlib.pyplot as plt
    import scipy.io as sio
    #------------------------------------------------------------------------------
    
    filename_fig = filename[:-4]+'_PSD.png'
    
    mask = np.isnan(image)!=1
    media = np.mean(image[mask])
    image [np.isnan(image)==1] = media
    
    plt.figure()
    plt.imshow(image)
    plt.colorbar()
    
    npix = image.shape[0]
    
    fourier_image = np.fft.fftn(image)
    fourier_amplitudes = np.abs(fourier_image)**2 #ampiezze al quadrato -> circa quasi PSD
    fourier_amplitudes = fourier_amplitudes * (DIA/npix**2)**2 # PSD
    
    kfreq = np.fft.fftfreq(npix) * npix *1/DIA # *1/DIA = from normalized freq to spatial freq
    kfreq2D = np.meshgrid(kfreq, kfreq)
    knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
    
    knrm = knrm.flatten()
    fourier_amplitudes = fourier_amplitudes.flatten()
    
    kbins = np.arange(0.5, npix//2+1, 1.)*1/DIA # *1/DIA = from normalized freq to spatial freq
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    
    Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes,
                                         statistic = "mean",
                                         bins = kbins)
    
    Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2) # 
    
    freq = kvals*1.
    PSD = Abins*1.
    
    RMS_image = np.sqrt(np.mean(image[mask]**2))
    RMS_PSD = np.sqrt(np.sum(PSD*(freq[-1]-freq[-2])))
    print('RMS image:___'+str(np.round(RMS_image*1e9,2))+' nm')
    print('RMS psd:_____'+str(np.round(RMS_PSD*1e9,2))+' nm')
    
    plt.figure()
    plt.loglog(freq, PSD,'k',linewidth=3)
    plt.xlabel("1/m")
    plt.ylabel("m$^3$")
    plt.title('RMS from PSD: '+str(np.round(RMS_PSD*1e9,2))+' nm \nRMS from Image: '+str(np.round(RMS_image*1e9,2))+' nm')
    plt.tight_layout()
    plt.grid(b=True, which='major', color='grey', linestyle='-')
    plt.grid(b=True, which='minor', color='grey', linestyle=':')
    # ISO limit A/f^B: A=1e-15, B=2
    plt.loglog(freq,1e-15/freq**2,'r:')
    plt.savefig(filename_fig, dpi = 300, bbox_inches = "tight")
       
    return(freq,PSD)