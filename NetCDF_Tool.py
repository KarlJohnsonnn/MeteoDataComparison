import numpy as np

import netCDF4
import datetime
from datetime import timezone
import glob
import os

import sys

from Parameter_Mod import *



# class LIMRad94_LV1 contains the radar data for a given time intervall within one day
# reflectivity is converted into log unit [dBZ], chirps and file transitions will
# concatinate automaticly
#

class LIMRad94_LV1():


    frequency = 94.0            # [GHz]
    radar_wavelength = 0.00319  # [m]
    beamwidth = 0.48            # []

    #def __init__(self, date, time_int, h_bounds):
    def __init__(self, *args):

        os.chdir(LIMRad_path)  # path to data needs to be fit to the devices file structure

        if len(args) < 1:
            print('You need to specify a date at least!')
        elif len(args) == 1:
            date     = args[0]
            time_int = '0000-2400'
            h_bounds = [0.0, 12.0]
        else:
            date     = args[0]
            time_int = args[1]
            h_bounds = args[2]

        comp_hours = [int(time_int[0:2]), int(time_int[5:7])]
        comp_minutes = [int(time_int[2:4]), int(time_int[7:9])]

        clock = np.array(comp_hours) + np.divide(comp_minutes, 60.)  # [hours] + [minutes]/60#

        # -- gathering year, month, day for convertion to UTC time
        plotyear = int('20' + date[:2])
        plotmonth = int(date[2:4])
        plotday = int(date[4:6])

        time = [0, 0, 0, 0]
        time[3] = datetime.datetime(plotyear, plotmonth, plotday, hour=int(comp_hours[1]), minute=int(comp_minutes[1]))
        time[2] = time[3] - datetime.timedelta(seconds=15)

        first_file = int(clock[0])
        if clock[1] - int(clock[1]) > 0.0:
            last_file = int(clock[1]) + 1
        else:
            last_file = int(clock[1])

        range_file_list = list(range(first_file, last_file, 1))

        self.ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob(date + '_' + str(il).zfill(2) + '*.LV1.NC'))
            self.ncfiles.append(file_name[2:-2])

            if file_name[2:-2] == '':
                print('   Error!  File: "' + file_name + '" not found --> exit!')
                print('   Check LIMRAD folder!')
                exit(0)

        file = self.ncfiles[0]
        n_nc_files = len(self.ncfiles)
        no_c = get_nc_dimension(file, 'Chirp')

        # find the number of range gates per chirp sequence,
        # also find the resolution of each chirp and
        # calculate the vector containing the height-steps
        self.range_gates = np.zeros((no_c,), dtype='int')
        self.cum_range_gates = np.zeros((no_c + 1,), dtype='int')
        self.range_res = np.array(get_nc_data(file, 'RangeRes'))
        self.height = [chirpTable_min_height]

        for ichirp in range(0, no_c):
            dummy = get_nc_data(file, 'C' + str(ichirp + 1) + 'MeanVel')
            self.range_gates[ichirp] = len(dummy[0, :])
            self.cum_range_gates[ichirp + 1] = self.cum_range_gates[ichirp] + self.range_gates[ichirp]
            n_height = len(self.height)
            for i in range(n_height, self.range_gates[ichirp] + n_height):
                self.height = np.append(self.height, self.height[i - 1] + self.range_res[ichirp])

        self.cum_time_gates = np.zeros((n_nc_files + 1,), dtype='int')
        time_samp = []
        self.CBH  = []
        self.DDTb = []
        self.LWP  = []
        self.Rain = []
        self.SurfPres   = []
        self.SurfRelHum = []
        self.SurfTemp   = []
        self.SurfWD     = []
        self.SurfWS     = []

        i_nc_file = 0
        for file in self.ncfiles:
            nc_data_set = netCDF4.Dataset(file, 'r')

            time_chirp = np.array(nc_data_set.variables['Time'])                        # Number of seconds since 1/1/2001 00:00:00 [UTC]
            time_samp = np.append(time_samp, time_chirp)
            self.cum_time_gates[i_nc_file + 1] = self.cum_time_gates[i_nc_file] + len(time_chirp)

            self.CBH  = np.append(self.CBH,  np.array(nc_data_set.variables['CBH']))    # Cloud Bottom Height [m]
            self.DDTb = np.append(self.DDTb, np.array(nc_data_set.variables['DDTb']))   # Direct detection brightness temperature [K]
            self.LWP  = np.append(self.LWP,  np.array(nc_data_set.variables['LWP']))    # Liquid Water Path [g/m2]
            self.Rain = np.append(self.Rain, np.array(nc_data_set.variables['Rain']))   # Rain rate from weather station [mm/h]

            self.SurfPres   = np.append(self.SurfPres,   np.array(nc_data_set.variables['SurfPres']))   # Surface atmospheric pressure from weather station [hPa]
            self.SurfRelHum = np.append(self.SurfRelHum, np.array(nc_data_set.variables['SurfRelHum'])) # Relative humidity from weather station [%]
            self.SurfTemp   = np.append(self.SurfTemp,   np.array(nc_data_set.variables['SurfTemp']))   # Surface temperature from weather station [K]
            self.SurfWD     = np.append(self.SurfWD,     np.array(nc_data_set.variables['SurfWD']))     # Surface wind direction from weather station [deg]
            self.SurfWS     = np.append(self.SurfWS,     np.array(nc_data_set.variables['SurfWS']))     # Surface wind speed from weather station [m/s]

            i_nc_file += 1
            nc_data_set.close()

        self.n_height = len(self.height)
        self.n_time = len(time_samp)

        # gahter radar data values and stack them together
        i_nc_file = 0

        Ze_chirps   = np.zeros((self.n_time, self.n_height))        # Equivalent radar reflectivity factor [mm6/m3]
        Ze45_chirps = np.zeros((self.n_time, self.n_height))        # Slanted equivalent radar reflectivity factor[mm6/m3]
        ZDR_chirps  = np.zeros((self.n_time, self.n_height))        # Differential reflectivity [dB]
        mdv_chirps  = np.zeros((self.n_time, self.n_height))        # Mean Doppler velocity [m/s]
        sw_chirps   = np.zeros((self.n_time, self.n_height))        # Spectrum width [m/s]
        ldr_chirps  = np.zeros((self.n_time, self.n_height))        # Slanted linear depolarization ratio [dB]
        kurt_chirps = np.zeros((self.n_time, self.n_height))        # Kurtosis [linear]
        DiffAtt_chirps = np.zeros((self.n_time, self.n_height))     # Differential attenuation [dBZ/km]
        Skew_chirps = np.zeros((self.n_time, self.n_height))        # Skewness [linear]

        for file in self.ncfiles:
            #if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

            nc_data_set = netCDF4.Dataset(file, 'r')

            lb_t = self.cum_time_gates[i_nc_file]
            ub_t = self.cum_time_gates[i_nc_file + 1]

            for ichirp in range(no_c):
                lb_h = self.cum_range_gates[ichirp]
                ub_h = self.cum_range_gates[ichirp + 1]
                Ze_chirps[lb_t:ub_t, lb_h:ub_h]   = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZE'])
                Ze45_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZE45'])
                ZDR_chirps[lb_t:ub_t, lb_h:ub_h]  = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZDR'])
                mdv_chirps[lb_t:ub_t, lb_h:ub_h]  = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'MeanVel'])
                sw_chirps[lb_t:ub_t, lb_h:ub_h]   = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SpecWidth'])
                ldr_chirps[lb_t:ub_t, lb_h:ub_h]  = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLDR'])
                kurt_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Kurt'])
                DiffAtt_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'DiffAtt'])
                Skew_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Skew'])

            nc_data_set.close()
            i_nc_file += 1

        #if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format

        time_plot = [datetime.datetime(2001, 1, 1, 0, 0, 0)
                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

        i = 0
        for zeit in time_plot:
            if time[0] <= zeit <= time[1]: min_time = i
            if time[2] <= zeit <= time[3]: max_time = i
            i += 1

        self.t_plt = time_plot[min_time:max_time]
        self.t_unix = [ts.replace(tzinfo=timezone.utc).timestamp() for ts in self.t_plt]

        imin_h, imax_h = get_height_boundary(self.height, 1000 * h_bounds[0], 1000 * h_bounds[1])
        self.height = np.divide(self.height[imin_h:imax_h], 1000)

        # build stacked chirps
        self.Ze = np.array(Ze_chirps[min_time:max_time, imin_h:imax_h])
        self.Ze = np.transpose(self.Ze)
        self.Ze = np.ma.masked_less_equal(self.Ze, 0.)
        self.Ze = np.ma.log10(self.Ze) * 10

        self.Ze45 = np.array(Ze45_chirps[min_time:max_time, imin_h:imax_h])
        self.Ze45 = np.transpose(self.Ze45)
        self.Ze45 = np.ma.masked_less_equal(self.Ze, 0.)
        self.Ze45 = np.ma.log10(self.Ze45) * 10

        self.ZDR = np.array(ZDR_chirps[min_time:max_time, imin_h:imax_h])
        self.ZDR = np.transpose(self.ZDR)
        self.ZDR = np.ma.masked_less_equal(self.ZDR, -999.)

        self.mdv = np.array(mdv_chirps[min_time:max_time, imin_h:imax_h])
        self.mdv = np.transpose(self.mdv)
        self.mdv = np.ma.masked_less_equal(self.mdv, -999.)

        self.sw = np.array(sw_chirps[min_time:max_time, imin_h:imax_h])
        self.sw = np.transpose(self.sw)
        self.sw = np.ma.masked_less_equal(self.sw, -999.)

        self.ldr = np.array(ldr_chirps[min_time:max_time, imin_h:imax_h])
        self.ldr = np.transpose(self.ldr)
        self.ldr = np.ma.masked_less_equal(self.ldr, -999.)

        self.kurt = np.array(kurt_chirps[min_time:max_time, imin_h:imax_h])
        self.kurt = np.transpose(self.kurt)
        self.kurt = np.ma.masked_less_equal(self.kurt, -999.)   # fill value correct?

        self.DiffAtt = np.array(DiffAtt_chirps[min_time:max_time, imin_h:imax_h])
        self.DiffAtt = np.transpose(self.DiffAtt)
        self.DiffAtt = np.ma.masked_less_equal(self.DiffAtt, -999.) # fill value correct?

        self.Skew = np.array(Skew_chirps[min_time:max_time, imin_h:imax_h])
        self.Skew = np.transpose(self.Skew)
        self.Skew = np.ma.masked_less_equal(self.Skew, -999.)

        self.kurt = np.array(kurt_chirps[min_time:max_time, imin_h:imax_h])
        self.kurt = np.transpose(self.kurt)
        self.kurt = np.ma.masked_less_equal(self.kurt, -999.)


def get_nc_data(thisfile, varname):
    # if pts: print('loading variable '+varname +' from ' + thisfile)
    ncfile = netCDF4.Dataset(thisfile, 'r')
    var = ncfile.variables[varname]

    if ncfile.isopen == 1: ncfile.close()
    return var


def get_nc_date(thisfile):
    # if pts: print('loading variable '+varname +' from ' + thisfile)
    ncfile = netCDF4.Dataset(thisfile, 'r')
    year = ncfile.year
    month = ncfile.month
    day = ncfile.day

    if ncfile.isopen == 1: ncfile.close()
    return year, month, day


def get_nc_dimension(thisfile, dim_name):
    ncfile = netCDF4.Dataset(thisfile, 'r')
    dim = ncfile.dimensions[dim_name].size
    if ncfile.isopen == 1: ncfile.close()
    return dim


def get_height_boundary(Array, hmin, hmax):
    i = 0
    imin = 0
    for h in Array:
        # print(' (min) h = ',h,hmin)
        if (hmin <= h):
            imin = i
            break
        i += 1

    i = 0
    imax = 0
    for h in reversed(Array):
        # print(' (max) h = ',h,hmax)
        if (hmax >= h):
            imax = len(Array) - i
            break
        i += 1

    return imin, imax


def extract_dataset(date, time, clock, h_bounds, fext, kind):
    '''
    Extract data from NetCDF files
    :param date: string of desired date 'YYMMDD'
    :param time: list of datetimes
    :param clock: deciaml hours [begin, end]
    :param fext:  file extention of dataset
    :return:
    '''
    if fext == '*mira.nc':
        os.chdir(meteo_path + 'MIRA/calibrated')  # path to data needs to be fit to the devices file structure
        ncfiles = glob.glob('20' + date + '*mira.nc')

        if pts: print("    Loading MIRA35 (mira.nc) NC-files ({} of {})".format(0, 1), end="\r")

        file = ncfiles[0]

        if file == '':
            print('   Error!  File: "' + file + '" not found --> exit!')
            print('   Check MIRA folder!')
            exit(0)

        height = np.array(get_nc_data(file, 'range'))
        imin_h, imax_h = get_height_boundary(height, h_bounds[0], h_bounds[1])

        # conversion from deciaml hour to datetime
        time_samp = np.array(get_nc_data(file, 'time'))

        mira_y, mira_m, mira_d = get_nc_date(ncfiles[0])
        dt_midnight = datetime.datetime(mira_y, mira_m, mira_d)

        time_plot = [dt_midnight + datetime.timedelta(seconds=int(t * 3600.)) for t in time_samp]

        i = 0
        for zeit in time_plot:
            if (time[0] <= zeit <= time[1]): min_time = i
            if (time[2] <= zeit <= time[3]): max_time = i
            i += 1

        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        height = height[imin_h:imax_h]

        # get data from .nc files
        Ze = np.array(get_nc_data(file, 'Zh'))
        mdv = np.array(get_nc_data(file, 'v'))
        sw = np.array(get_nc_data(file, 'width'))

        # np.warnings.filterwarnings('ignore')

        Ze = Ze[min_time:max_time, imin_h:imax_h]
        mdv = mdv[min_time:max_time, imin_h:imax_h]
        sw = sw[min_time:max_time, imin_h:imax_h]

        # stack variables of individual chirps
        Ze = np.transpose(Ze)
        Ze = np.ma.masked_less_equal(Ze, -999.)

        # conversion to numpy array for truncation
        mdv = np.transpose(mdv)
        mdv = np.ma.masked_less_equal(mdv, -999.)

        sw = np.transpose(sw)
        sw = np.ma.masked_invalid(sw)
        sw = np.ma.masked_less_equal(sw, -999.)

        os.chdir('../../')  # path to data needs to be fit to the devices file structure
        if pts: print("    Loading MIRA35 (mira.nc) NC-files ({} of {})".format(1, 1))

    elif fext == '*.mmclx':

        os.chdir(meteo_path + 'MIRA/mmclx')  # path to data needs to be fit to the devices file structure

        # ncfiles = glob.glob('20' + date + '*.mmclx')     # 20180727_000013.mmclx

        first_file = int(clock[0]) - np.remainder(int(clock[0]), 3)
        if clock[1] - int(clock[1]) > 0.0:
            last_file = int(clock[1]) + 1
        else:
            last_file = int(clock[1])

        range_file_list = list(range(first_file, last_file, 3))

        ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob('20' + date + '_' + str(il).zfill(2) + '*.mmclx'))
            ncfiles.append(file_name[2:-2])

        if file_name[2:-2] == '':
            print('   Error!  File not found --> exit!')
            print('   Check LIMRAD folder!')
            exit(0)

        file = ncfiles[0]

        if file == '':
            print('   Error!  File: "' + file + '" not found --> exit!')
            print('   Check MIRA folder!')
            exit(0)

        # extract date (year, month, day)
        # mira_y, mira_m, mira_d = get_nc_date(ncfiles[0])

        # extract range array
        height = np.array(get_nc_data(file, 'range'))
        imin_h, imax_h = get_height_boundary(height, h_bounds[0] * 1000, h_bounds[1] * 1000)
        # 20180728_060014.mmclx

        if kind == 'cl':
            kind1 = ''
        elif kind == 'g':
            kind1 = 'g'
        else:
            kind1 = 'e'

        # conversion from deciaml hour to datetime
        time_samp = np.array(get_nc_data(file, 'time'))
        Ze = np.array(get_nc_data(file, 'Z' + kind1))
        mdv = np.array(get_nc_data(file, 'VEL' + kind))
        sw = np.array(get_nc_data(file, 'RMS' + kind))

        i_nc_file = 0
        n_nc_file = len(ncfiles)

        for file in ncfiles[1:]:
            i_nc_file += 1
            if pts: print("    Loading MIRA35 (mmclx) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

            time_samp = np.append(time_samp, get_nc_data(file, 'time'))
            Ze = np.append(Ze, get_nc_data(file, 'Z' + kind1), axis=0)
            mdv = np.append(mdv, get_nc_data(file, 'VEL' + kind), axis=0)
            sw = np.append(sw, get_nc_data(file, 'RMS' + kind), axis=0)

        if pts: print("    Loading MIRA35 (mmclx) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        time_plot = [datetime.datetime(1970, 1, 1, 0, 0, 0)
                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

        i = 0
        for zeit in time_plot:
            if (time[0] <= zeit <= time[1]): min_time = i
            if (time[2] <= zeit <= time[3]): max_time = i
            i += 1

        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        height = np.divide(height[imin_h:imax_h], 1000.0)

        Ze = Ze[min_time:max_time, imin_h:imax_h]
        mdv = mdv[min_time:max_time, imin_h:imax_h]
        sw = sw[min_time:max_time, imin_h:imax_h]

        Ze = np.ma.masked_invalid(Ze).T
        Ze = np.ma.log10(Ze) * 10
        mdv = np.ma.masked_invalid(mdv).T
        sw = np.ma.masked_invalid(sw).T

        os.chdir('../../')  # path to data needs to be fit to the devices file structure


    elif fext == '*.LV1.NC':

        os.chdir(meteo_path + 'LIMRad94/')  # path to data needs to be fit to the devices file structure

        first_file = int(clock[0])
        if clock[1] - int(clock[1]) > 0.0:
            last_file = int(clock[1]) + 1
        else:
            last_file = int(clock[1])

        range_file_list = list(range(first_file, last_file, 1))

        ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob(date + '_' + str(il).zfill(2) + '*.LV1.NC'))
            ncfiles.append(file_name[2:-2])

            if file_name[2:-2] == '':
                print('   Error!  File: "' + file_name + '" not found --> exit!')
                print('   Check LIMRAD folder!')
                exit(0)

        # initialize variables
        height = []
        range_gates = []
        dummy = None

        file = ncfiles[0]

        no_c = get_nc_dimension(file, 'Chirp')

        # find the number of range gates per chirp sequence
        range_res = get_nc_data(file, 'RangeRes')
        for i in range(0, no_c):
            dummy = get_nc_data(file, 'C' + str(i + 1) + 'MeanVel')
            range_gates = np.append(range_gates, len(dummy[0, :]))

        # calculate LR_height levels
        for i in range(0, int(range_gates[0])):
            height = np.append(height, (i + 1) * range_res[0])
        for i in range(0, int(range_gates[1])):
            i_temp = int(range_gates[0] - 1 + i)
            height = np.append(height, height[i_temp] + range_res[1])
        for i in range(0, int(range_gates[2])):
            i_temp = int(range_gates[0] + range_gates[1] - 1 + i)
            height = np.append(height, height[i_temp] + range_res[2])

        # get data from .nc files
        time_samp = get_nc_data(file, 'Time')
        Ze1 = get_nc_data(file, 'C1ZE')
        Ze2 = get_nc_data(file, 'C2ZE')
        Ze3 = get_nc_data(file, 'C3ZE')
        mdv1 = get_nc_data(file, 'C1MeanVel')
        mdv2 = get_nc_data(file, 'C2MeanVel')
        mdv3 = get_nc_data(file, 'C3MeanVel')
        sw1 = get_nc_data(file, 'C1SpecWidth')
        sw2 = get_nc_data(file, 'C2SpecWidth')
        sw3 = get_nc_data(file, 'C3SpecWidth')

        i_nc_file = 0
        n_nc_file = len(ncfiles)

        for file in ncfiles[1:]:
            i_nc_file += 1
            if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

            time_samp = np.append(time_samp, get_nc_data(file, 'Time'))
            Ze1 = np.append(Ze1, get_nc_data(file, 'C1ZE'), axis=0)
            Ze2 = np.append(Ze2, get_nc_data(file, 'C2ZE'), axis=0)
            Ze3 = np.append(Ze3, get_nc_data(file, 'C3ZE'), axis=0)
            mdv1 = np.append(mdv1, get_nc_data(file, 'C1MeanVel'), axis=0)
            mdv2 = np.append(mdv2, get_nc_data(file, 'C2MeanVel'), axis=0)
            mdv3 = np.append(mdv3, get_nc_data(file, 'C3MeanVel'), axis=0)
            sw1 = np.append(sw1, get_nc_data(file, 'C1SpecWidth'), axis=0)
            sw2 = np.append(sw2, get_nc_data(file, 'C2SpecWidth'), axis=0)
            sw3 = np.append(sw3, get_nc_data(file, 'C3SpecWidth'), axis=0)

        if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format

        time_plot = [datetime.datetime(2001, 1, 1, 0, 0, 0)
                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

        i = 0
        for zeit in time_plot:
            if (time[0] <= zeit <= time[1]): min_time = i
            if (time[2] <= zeit <= time[3]): max_time = i
            i += 1

        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        # for zeit in UTC_time_LR:
        #    print('zeit = ', zeit)

        height = height + chirpTable_min_height * 1000.  # get LR_height in km

        imin_h, imax_h = get_height_boundary(height, 1000 * h_bounds[0], 1000 * h_bounds[1])
        height = np.divide(height[imin_h:imax_h], 1000)

        # stack variables of individual chirps
        Ze = np.hstack((Ze1, Ze2, Ze3))
        Ze = Ze[min_time:max_time, imin_h:imax_h]
        Ze = np.transpose(Ze)
        Ze = np.ma.masked_less_equal(Ze, 0.)
        Ze = np.ma.log10(Ze) * 10

        mdv = np.hstack((mdv1, mdv2, mdv3))
        mdv = mdv[min_time:max_time, imin_h:imax_h]
        mdv = np.transpose(mdv)
        mdv = np.ma.masked_less_equal(mdv, -999.)

        sw = np.hstack((sw1, sw2, sw3))
        sw = sw[min_time:max_time, imin_h:imax_h]
        sw = np.transpose(sw)
        sw = np.ma.masked_less_equal(sw, -999.)

        os.chdir('../')  # path to data needs to be fit to the devices file structure


    else:
        sys.exit('Error!  Unknown data file extension.')

    # convert datetime to unix time (seconds since 1.1.1970)
    time_samp = [(ts.replace(tzinfo=timezone.utc).timestamp()) for ts in time_plot]

    return time_samp, time_plot, height, Ze, mdv, sw


def save_log_data(filename, meth, hmin, hmax, comp_date, comp_time_int):
    file = open(filename + '.log', 'w')

    file.write('')
    file.write(' This is the log file for LIMRad 94GHz - MIRA 35GHz Radar Data' + '\n' * 2)
    file.write(' # User Inputs\n')
    file.write('    minimum height = ' + str(hmin) + '\n')
    file.write('    maximum height = ' + str(hmax) + '\n' * 2)

    file.write('    date = ' + str(comp_date) + ' # in "YYMMDD" \n')
    file.write('    time intervall = ' + str(comp_time_int) + ' # in "HHMM-HHMM"' + '\n' * 2)

    file.write(' # Logicals \n')
    file.write('    pts = ' + str(pts) + ' # print to screen\n')
    file.write('    dbg = ' + str(dbg) + ' # debugging flag, show some parameter\n')
    file.write(
        '    plot_RectBivariateSpline   = ' + str(plot_RectBivariateSpline) + ' # bivariate interpolation plot\n')
    file.write('    plot_interpolation_scatter = ' + str(
        plot_interpolation_scatter) + ' # phase diagram of eqv. radar reflectivity\n')
    file.write('    plot_compare_mira_mmclx    = ' + str(
        plot_compare_mira_mmclx) + ' # plot Ze, mdv, sw of LIMRad and MIRA (all targes, only hydrometeors)\n')
    file.write(
        '    plot_radar_results = ' + str(plot_radar_results) + ' # plot radar data of LIMRad and MIRA (Ze, mdv, sw)\n')
    file.write('    plot_comparisons   = ' + str(plot_comparisons) + ' # plot time/height-averaged data\n')
    file.write('    plot_interp2d      = ' + str(
        plot_interp2d) + ' # plot 2D interpolation; LIMRad data onto MIRA grid and vice versa' + '\n' * 2)

    if plot_interpolation_scatter:
        file.write(' # phase-diagram interpolation parameter\n')
        file.write('    interpolation method = ' + meth + '\n')
        # file.write('    resolution of interpolated points = ' + str(res_interp) + '\n')

    file.close()
