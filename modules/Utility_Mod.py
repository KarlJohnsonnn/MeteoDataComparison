import calendar
import datetime
import os, sys
import numpy as np
import pandas as pd
from numba import jit

from modules.Parameter_Mod import *


def Print_Head():
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '          LIMRAD94 - MIRA35  Comparison          ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)


# Unixtime conversion functions (from Johannes):
def seconds_since_epoch(dtm):
    """ """
    return calendar.timegm(dtm.utctimetuple())


def datetime_from_seconds(seconds):
    """ """
    return datetime.datetime.utcfromtimestamp(seconds)


def findBasesTops(dbz_m, range_v):
    """
    % FINDBASESSTOPS
    %
    % functionality:
    % find cloud bases and tops from radar reflectivity profiles for up to 10 cloud layers
    % no cloud = NaN
    %
    % input:
    %   dbz       ... reflectivity matrix         [dbz] (range x time)
    %   range     ... radar height vector         [m or km] (range)
    %
    % output:
    %   bases     ... matrix of indices (idx) of cloud bases  (10 x time)
    %   tops      ... matrix of indices (idx) of cloud tops   (10 x time)
    %   base_m    ... matrix of heights of cloud bases  [m or km, same unit as range] (10 x time), 1st base = -1 means no cloud detected
    %   top_m     ... matrix of heights of cloud tops   [m or km, same unit as range] (10 x time), 1st top  = -1 means no cloud detected
    %   thickness ... matrix of cloud thickness         [m or km, same unit as range] (10 x time)
    """

    shape_dbz = dbz_m.shape
    len_time = shape_dbz[1]
    len_range = len(range_v)

    bases = np_NaN(10, len_time)
    tops = np_NaN(10, len_time)
    thickness = np_NaN(10, len_time)

    top_m = np_NaN(10, len_time)  # max. 10 cloud layers detected
    base_m = np_NaN(10, len_time)  # max. 10 cloud layers detected

    if pts:
        print('')
        print(' dBZ(i,:) = ', [dbz_m[i, :] for i in range(shape_dbz[0])])

    for i in range(0, len_time):

        in_cloud = 0
        layer_idx = 0
        current_base = np.nan

        if pts: print("    Searching for cloud bottom and top ({} of {}) time steps".format(i + 1, len_time), end="\r")

        # found the first base in first bin.
        if not np.isnan(dbz_m[0, i]):
            layer_idx = 1
            current_base = 1
            in_cloud = 1

        for j in range(1, len_range):

            if in_cloud == 1:  # if in cloud

                # cloud top found at (j-1)
                if np.isnan(dbz_m[j, i]):
                    current_top = j - 1
                    thickness[layer_idx, i] = range_v[current_top] - range_v[current_base]
                    bases[layer_idx, i] = current_base  # bases is an idx
                    tops[layer_idx, i] = current_top  # tops is an idx

                    base_m[layer_idx, i] = range_v[current_base]  # cloud base in m or km
                    top_m[layer_idx, i] = range_v[current_top]  # cloud top in m or km

                    print(str(i) + ': found ' + str(layer_idx) + '. cloud [' + str(bases[layer_idx, i]) + ', ' +
                          str(tops[layer_idx, i]) + '], thickness: ' + str(thickness) + 'km')

                    in_cloud = 0

            else:  # if not in cloud

                # cloud_base found at j
                if not np.isnan(dbz_m[j, i]):
                    layer_idx += 1
                    current_base = j
                    in_cloud = 1

        # at top height but still in cloud, force top
        if in_cloud == 1:
            tops[layer_idx, i] = len(range_v)
            top_m[layer_idx, i] = max(range_v)  # cloud top in m or km

    ###
    # keep only first 10 cloud layers
    bases = bases[:10, :]
    tops = tops[:10, :]
    base_m = base_m[:10, :]
    top_m = top_m[:10, :]
    thickness = thickness[:10, :]
    # give clear sky flag when first base_m ==NaN (no dbz detected over all heights),
    # problem: what if radar wasn't working, then dbz would still be NaN!
    loc_nan = np.where(np.isnan(base_m[0, :]))

    if pts:
        print('vor bases = ', [ba for ba in base_m])

    base_m[0, np.where(np.isnan(base_m[0, :]))] = -1
    top_m[0, np.where(np.isnan(top_m[0, :]))] = -1

    if pts:
        print('')
        print('npisnan = ', [loc_nan[i] for i in range(len(loc_nan))])

        if pts: print('')
        for ba in base_m:
            print(' bases = ', ba)

        for to in top_m:
            print(' tops = ', to)

    return bases, tops, base_m, top_m, thickness


def dim(a):
    if not type(a) == list:  return []
    return [len(a)] + dim(a[0])


def string_to_datetime(ds, string):
    hour, minu, sec = string.split(':')
    return datetime.datetime(ds.year, ds.month, ds.day, hour=int(hour), minute=int(minu), second=int(sec))


def np_NaN(n, m):
    mat = np.zeros((n, m))
    mat[:, :] = np.nan
    return mat


def lookupNearest(x0, y0, x, y, data):
    xi = np.abs(x - x0).argmin()
    yi = np.abs(y - y0).argmin()
    return data[yi, xi]


def correlation(v1, v2):
    rho = np.array([[]])
    for i in range(len(v1[0, :])):
        df_v1 = pd.DataFrame(v1[:, i])
        df_v2 = pd.DataFrame(v2[:, i])
        rho = np.append(rho, df_v1.corrwith(df_v2))

    return np.ma.masked_invalid(rho)


def gather_user_input(LV0):
    str_time = input('     Enter time in (HH:MM:SS) ::  ')
    str_height = input('     Enter height in (km)     ::  ')

    bsp_time = string_to_datetime(LV0, str_time)
    bsp_height = float(str_height)

    bsp_height0 = min(LV0.height_all, key=lambda x: abs(x - bsp_height))
    bsp_time0 = min(LV0.t_plt, key=lambda x: abs(x - bsp_time))

    itime = LV0.t_plt.index(bsp_time0)

    for ic in range(LV0.no_c):
        try:
            idx_height = list(LV0.height[ic]).index(bsp_height0)
            if idx_height > 0:
                ichirp = ic
                iheight = idx_height
                break
        except Exception as e:
            print('Something went wrong X-(    ' + str(e))
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, ' at Line ', exc_tb.tb_lineno)
            sys.exit(0)

    return ichirp, itime, iheight


@jit(nopython=True, fastmath=True)
def estimate_noise_hs74(spectrum, navg=1, std_div=0.0):
    """
    Estimate noise parameters of a Doppler spectrum.
    Use the method of estimating the noise level in Doppler spectra outlined
    by Hildebrand and Sehkon, 1974.
    Parameters
    ----------
    spectrum : array like
        Doppler spectrum in linear units.
    navg : int, optional
        The number of spectral bins over which a moving average has been
        taken. Corresponds to the **p** variable from equation 9 of the
        article.  The default value of 1 is appropiate when no moving
        average has been applied to the spectrum.
    Returns
    -------
    mean : float-like
        Mean of points in the spectrum identified as noise.
    threshold : float-like
        Threshold separating noise from signal.  The point in the spectrum with
        this value or below should be considered as noise, above this value
        signal. It is possible that all points in the spectrum are identified
        as noise.  If a peak is required for moment calculation then the point
        with this value should be considered as signal.
    var : float-like
        Variance of the points in the spectrum identified as noise.
    nnoise : int
        Number of noise points in the spectrum.
    References
    ----------
    P. H. Hildebrand and R. S. Sekhon, Objective Determination of the Noise
    Level in Doppler Spectra. Journal of Applied Meteorology, 1974, 13,
    808-811.
    """
    n_spec = len(spectrum)
    sorted_spectrum = np.sort(spectrum)
    nnoise = n_spec  # default to all points in the spectrum as noise
    for npts in range(1, n_spec + 1):
        partial = sorted_spectrum[:npts]
        mean = np.mean(partial)
        var = np.var(partial)
        if var * navg < mean * mean:
            nnoise = npts
        else:
            # partial spectrum no longer has characteristics of white noise
            break

    noise_spectrum = sorted_spectrum[:nnoise]
    mean = np.mean(noise_spectrum)
    var = np.var(noise_spectrum)

    if std_div == -1.0:
        threshold = sorted_spectrum[nnoise - 1]
    else:
        threshold = mean + np.sqrt(var) * std_div
        # print(' normthres = ', sorted_spectrum[nnoise - 1])
        # print(' mean+std  = ', mean + np.sqrt(var) * std_div)

    left_intersec = -1
    right_intersec = -1

    if nnoise < n_spec:
        idxMaxSignal = np.argmax(spectrum)

        for ispec in range(idxMaxSignal, n_spec):
            if spectrum[ispec] <= threshold:
                right_intersec = ispec
                break

        for ispec in range(idxMaxSignal, -1, -1):
            if spectrum[ispec] <= threshold:
                left_intersec = ispec
                break

    #    left_intersec  = -1
    #    right_intersec = -1
    #    if nnoise < n_spec:
    #        for ispec in range(n_spec):
    #            if spectrum[ispec] > threshold:
    #                left_intersec  = ispec-1
    #                break
    #
    #        for ispec in range(n_spec-1, -1, -1):
    #            if spectrum[ispec] > threshold:
    #                right_intersec = ispec+1
    #                break

    return mean, threshold, var, nnoise, left_intersec, right_intersec


# @jit(nopython=True, fastmath=True)
def remove_noise(ds, std_div=0.0):
    """
    :param std_div: how many standart diviations from mean noise will be added to mean noise for threshold
    :param ds:  data set from LIMRAD94 LV0
    :return:    noise floor estimation for all time and range points
    """

    n_t = ds.n_time

    # allocate numpy arrays
    mean_noise = []
    threshold = []
    variance = []
    numnoise = []
    integration_bounds = []

    for ic in range(ds.no_c):
        n_r = ds.n_height[ic]
        mean_noise.append(np.zeros((n_t, n_r)))
        threshold.append(np.zeros((n_t, n_r)))
        variance.append(np.zeros((n_t, n_r)))
        numnoise.append(np.zeros((n_t, n_r)))
        integration_bounds.append(np.zeros((n_t, n_r, 2)))

    # gather noise level etc. for all chirps, range gates and times
    for ic in range(ds.no_c):
        for iR in range(ds.n_height[ic]):
            for iT in range(ds.n_time):
                mean, thresh, var, nnoise, left_intersec, right_intersec = \
                    estimate_noise_hs74(ds.VHSpec[ic][iT, iR, :], navg=ds.no_av[ic], std_div=std_div)

                mean_noise[ic][iT, iR] = mean
                variance[ic][iT, iR] = var
                numnoise[ic][iT, iR] = nnoise
                threshold[ic][iT, iR] = thresh
                integration_bounds[ic][iT, iR, :] = [left_intersec, right_intersec]

    return mean_noise, threshold, variance, numnoise, integration_bounds


def spectra_to_moments(spectra_linear_units, velocity_bins, bounds, DoppRes):
    """
    # spectra_to_moments
    # translated from Heike's Matlab function
    # determination of radar moments of Doppler spectrum over range of Doppler velocity bins

    # input:
    # spectra_linear_units: time-height-FFT points of Doppler spectra ([mm^6 / m^3 ] / (m/s))
    # with noise removed(!?)
    # velocity_bins         : FFTpoint-long spectral velocity bins (m/s)
    # left_edge       : left edge of Doppler spectrum (v1) for moment determination
    # right_edge      : right edge of Doppler spectrum (v2) for moment determination

    # left_edge and right_edge can be determined using findEdges.py

    # output:
    # Ze              : 0. moment = reflectivity over range of Doppler velocity bins v1 to v2 [mm6/m3]
    # mdv             : 1. moment = mean Doppler velocity over range of Doppler velocity bins v1 to v2 [m/s]
    # sw              : 2. moment = spectrum width over range of Doppler velocity bins v1 to v2  [m/s]
    # skew            : 3. moment = skewness over range of Doppler velocity bins v1 to v2
    # kurt            : 4. moment = kurtosis over range of Doppler velocity bins v1 to v2
    # pwr_nrm_out     : normalized power (the sum of pwr_norm is 1)
    """

    # contains the dimensionality of the Doppler spectrum, (nTime, nRange, nDopplerbins)
    no_chirps = len(spectra_linear_units)
    no_times = spectra_linear_units[0].shape[0]
    no_Dbins = [velocity_bins[ic].size for ic in range(no_chirps)]
    no_ranges_chrip = [spectra_linear_units[ic].shape[1] for ic in range(no_chirps)]
    no_ranges = sum(no_ranges_chrip)

    # initialize variables:
    # create empty arrays for output
    Ze = np.full((no_times, no_ranges), np.nan)
    mdv = np.full((no_times, no_ranges), np.nan)
    sw = np.full((no_times, no_ranges), np.nan)
    skew = np.full((no_times, no_ranges), np.nan)
    kurt = np.full((no_times, no_ranges), np.nan)

    pwr_nrm_out = [None] * no_chirps
    delta_vel = [None] * no_chirps

    for ic in range(no_chirps):  # ith chirp, depends on chirp table configuration

        # add new list element, mean velocity difference between velocity bins:
        pwr_nrm_out[ic] = np.full((no_times, no_ranges, no_Dbins[ic]), np.nan)
        delta_vel[ic] = np.nanmean(np.diff(velocity_bins[ic]))

        for iR in range(no_ranges_chrip[ic]):  # range dimension
            for iT in range(no_times):  # time dimension

                if bounds[ic][iT, iR, 0] > -1:  # check if signal was detected by estimate_noise routine

                    lb = int(bounds[ic][iT, iR, 0])

                    if bounds[ic][iT, iR, 1] < 0:
                        ub = None
                    else:
                        ub = int(bounds[ic][iT, iR, 1])

                    if ic > 0:
                        iR_out = iR + sum(no_ranges_chrip[:ic])
                    else:
                        iR_out = iR

                    # hier signal durch 2 teilen .... warum?? weil VH = Vspec+Hspec???? ask Alexander Myagkov why
                    signal = spectra_linear_units[ic][iT, iR, lb:ub]  # extract power spectra in chosen range

                    #if ic == 0:
                    #    if np.sum(signal) < 5.e-6:
                    #        hydro_meteor = False
                    #    elif np.sum(np.ma.diff(np.ma.masked_less_equal(signal, -999.))) < 5.e-6:
                    #        hydro_meteor = False
                    #    else:
                    #        hydro_meteor =  True
                    #else:
                    hydro_meteor = True

                    velocity_bins_extr = velocity_bins[ic][lb:ub]  # extract velocity bins in chosen Vdop bin range
                    signal_sum = np.nansum(signal)  # linear full spectrum Ze [mm^6/m^3], scalar

                    if np.isfinite(signal_sum) and hydro_meteor:  # check if Ze_linear is not NaN

                        # Ze_linear = signal_sum
                        Ze_linear = signal_sum / 2.0
                        Ze[iT, iR_out] = Ze_linear  # copy temporary Ze_linear variable to output variable

                        pwr_nrm = signal / signal_sum  # determine normalized power (NOT normalized by Vdop bins)
                        pwr_nrm_out[ic][iT, iR, lb:ub] = pwr_nrm  # create output matrix of normalized power

                        mdv[iT, iR_out] = np.nansum(velocity_bins_extr * pwr_nrm)
                        sw[iT, iR_out] = np.sqrt(
                            np.abs(np.nansum(np.multiply(pwr_nrm, np.square(velocity_bins_extr - mdv[iT, iR_out])))))

                        skew[iT, iR_out] = np.nansum(
                            pwr_nrm * np.power(velocity_bins_extr - mdv[iT, iR_out], 3.0)) / np.power(sw[iT, iR_out],
                                                                                                      3.0)
                        kurt[iT, iR_out] = np.nansum(
                            pwr_nrm * np.power(velocity_bins_extr - mdv[iT, iR_out], 4.0)) / np.power(sw[iT, iR_out],
                                                                                                      4.0)

                        mdv[iT, iR_out] = mdv[iT, iR_out] - DoppRes[ic] / 2.0  # ask Alexander Myagkov why

    Ze = np.ma.masked_invalid(Ze)
    mdv = np.ma.masked_invalid(mdv)
    sw = np.ma.masked_invalid(sw)
    skew = np.ma.masked_invalid(skew)
    kurt = np.ma.masked_invalid(kurt)

    return Ze, mdv, sw, skew, kurt, pwr_nrm_out


def compare_datasets(lv0, lv1):
    # Z_norm = 10.0 * np.log10(np.linalg.norm(np.ma.subtract(Ze1, Ze2), ord='fro'))
    # mdv_norm = np.linalg.norm(np.subtract(lv0.mdv, lv1.mdv), ord='fro')
    # sw_norm  = np.linalg.norm(np.subtract(lv0.sw, lv1.sw), ord='fro')
    Z_norm = np.mean(np.ma.subtract(lv0.ZeLin, lv1.ZeLin))
    mdv_norm = np.mean(np.subtract(lv0.mdv, lv1.mdv))
    sw_norm = np.mean(np.subtract(lv0.sw, lv1.sw))

    # convert to dBZ
    print()
    print(f'    ||Ze_lv0  -  Ze_lv1|| = {Z_norm:.12f} [mm6/m3]')
    print(f'    ||mdv_lv0 - mdv_lv1|| = {mdv_norm:.12f} [m/s]')
    print(f'    ||sw_lv0  -  sw_lv1|| = {sw_norm:.12f} [m/s]')
    print()

    pass


# compare elements of two lists
def Diff(li1, li2):
    return (list(set(li1) - set(li2)))


def Create_NeuralNet_Input(ds):
    from sklearn.preprocessing import MinMaxScaler

    Z = ds.variables['Z']
    v = ds.variables['v']
    width = ds.variables['width']
    ldr = ds.variables['ldr']

    beta = ds.variables['beta']
    lidar_depolarisation = ds.variables['beta']

    Times = []
    Heights = []
    radar_data = []
    lidar_label = []

    for iT in range(ds.dimensions['time']):
        for iH in range(ds.dimensions['height']):

            # if beta[iT, iH] > -999.0 and Z[iT, iH] > -999.0:
            if beta[iT, iH] and Z[iT, iH] > -999.0:
                Times.append(iT)
                Heights.append(iH)

                # assign radar moments reflectivity, mean doppler velocity, spectral width, linear deplo ratio
                Zij = Z[iT, iH]
                vij = v[iT, iH]
                widthij = width[iT, iH]
                ldrij = ldr[iT, iH]
                radar_data.append(np.array([iT, iH, Zij, vij, widthij, ldrij]))

                # assign lidar moments beta, depol
                betaij = beta[iT, iH]
                lidar_depolarisationij = lidar_depolarisation[iT, iH]
                lidar_label.append(np.array([betaij, lidar_depolarisationij]))

    Times = np.array(Times)
    Heights = np.array(Heights)

    scaler = MinMaxScaler(feature_range=(0, 1))
    scaled_radar_data = scaler.fit_transform(radar_data)
    scaled_lidar_label = scaler.fit_transform(lidar_label)

    return scaled_radar_data, scaled_lidar_label
