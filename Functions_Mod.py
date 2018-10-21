import datetime

import numpy as np
import pandas as pd
from numba import jit

from Parameter_Mod import *


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


@jit(nopython=True, fastmath=True)
def estimate_noise_hs74(spectrum, navg=1):
    sorted_spectrum = np.sort(spectrum)
    nnoise = len(spectrum)  # default to all points in the spectrum as noise
    for npts in range(1, len(sorted_spectrum) + 1):
        partial = sorted_spectrum[:npts]
        mean = np.mean(partial)
        var = np.var(partial)
        if var * navg < mean ** 2.:
            nnoise = npts
        else:
            # partial spectrum no longer has characteristics of white noise
            break

    noise_spectrum = sorted_spectrum[:nnoise]
    mean = np.mean(noise_spectrum)
    threshold = sorted_spectrum[nnoise - 1]
    var = np.var(noise_spectrum)
    return mean, threshold, var, nnoise
