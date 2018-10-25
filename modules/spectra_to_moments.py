'''
# spectra_to_moments
# translated from Heike's Matlab function
# determination of radar moments of Dopler spectrum over range of Doppler velocity bins

# input:
# spectra_lin_norm: time-height-FFT points of Doppler spectra ([mm^6 / m^3 ] / (m/s))
# with noise removed(!?)
# velbins         : FFTpoint-long spectral velocity bins (m/s)
# left_edge       : left edge of Doppler spectrum (v1) for moment determination
# right_edge      : right edge of Doppler spectrum (v2) for moment determination

# left_edge and right_edge can be determined using findEdges.py

# output: 
# Ze              : 0. moment = reflectivity over range of Doppler velocity bins v1 to v2 [dBZ]
# mdv             : 1. moment = mean Doppler velocity over range of Doppler velocity bins v1 to v2 [m/s]
# spectrum_width  : 2. moment = spectrum width over range of Doppler velocity bins v1 to v2  [m/s]
# skew            : 3. moment = skewness over range of Doppler velocity bins v1 to v2 
# kurt            : 4. moment = kurtosis over range of Doppler velocity bins v1 to v2 
# pwr_nrm_out     : normalized power (the sum of pwr_norm is 1)
'''
import numpy as np

# velbins = the_bins
# spectra_lin_norm=np.ma.masked_less_equal(the_hspec, -999.)+np.ma.masked_less_equal(the_vspec, -999.)
# left_edge = 0
# right_edge = 128

def spectra_to_moments(velbins, spectra_lin_norm, left_edge,right_edge):
   
    # initialize variables:
    # create empty arrays for output    
    Ze = np.full((spectra_lin_norm.shape[0],spectra_lin_norm.shape[1]),-999.0000)
    mdv= np.full((spectra_lin_norm.shape[0],spectra_lin_norm.shape[1]),-999.0000)
    spectrum_width = np.full((spectra_lin_norm.shape[0],spectra_lin_norm.shape[1]),-999.0000)
    pwr_nrm    = np.full((spectra_lin_norm.shape[0],spectra_lin_norm.shape[1],velbins.shape[0]),-999.0000)
    skew = np.full((spectra_lin_norm.shape[0],spectra_lin_norm.shape[1]),-999.0000)
    kurt = np.full((spectra_lin_norm.shape[0],spectra_lin_norm.shape[1]),-999.0000)
   
    # velocity difference between velocity bins:
    delta_vel=np.nanmean([j-i for i, j in zip(velbins[:-1], velbins[1:])])
    # Python code explanation:
    # [:-1] returns all elements [:] except the last one -1
    # [1:] returns all elements [:] except the first one which would be "0"
    
    print('Determine Doppler spectrum moments...')
    
    for i in range(0,spectra_lin_norm.shape[0]):
        for j in range(0,spectra_lin_norm.shape[1]):
            # extract power spectra and velocity bins in chosen range:
            spectra_lin_norm_extr=spectra_lin_norm[i,j,left_edge:right_edge]
            velbins_extr=velbins[left_edge:right_edge]
            ze_lin_int = np.nansum(spectra_lin_norm_extr * delta_vel)
            ze_lin_int_no_dV = np.nansum(spectra_lin_norm_extr )
   
            if not (ze_lin_int is np.ma.masked):
                Ze[i,j]=(10*np.log10(ze_lin_int))
                pwr_nrm[i,j,left_edge:right_edge]=(spectra_lin_norm_extr/ze_lin_int_no_dV)
                mdv[i,j]=(np.nansum(velbins_extr*pwr_nrm[i,j,left_edge:right_edge]))
                #spectrum_width.append(np.sqrt(np.abs(np.nansum(pwr_nrm * (velbins_extr - mdv))^2)))
                #skew.append(np.nansum(pwr_nrm * (velbins_extr . mdv)^3) / spectrum_width^3)
                #kurt.append(np.nansum (pwr_nrm * velbins_extr - mdv ^4) / spectrum_width^4)
         
    return Ze, mdv, spectrum_width, skew, kurt, pwr_nrm


    