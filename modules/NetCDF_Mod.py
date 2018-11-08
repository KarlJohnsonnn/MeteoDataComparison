import datetime
import glob
import os
import time
from datetime import timezone

import netCDF4
import numpy as np

from modules.Interpolation_Mod import interpolate2d
from modules.Parameter_Mod import *


class LIMRAD94_LV0():
    frequency = 94.0  # [GHz]
    radar_wavelength = 0.00319  # [m]
    beamwidth = 0.48  # []

    def __init__(self, *args):

        os.chdir(LIMRAD_path)  # path to data needs to be fit to the devices file structure

        if len(args) < 1:
            print('You need to specify a date at least!')
            exit(0)

        elif len(args) == 1:
            date = args[0]
            time_int = '0000-2400'
            h_bounds = [0.0, 12.0]

        else:
            date = args[0]
            time_int = args[1]
            h_bounds = args[2]

        comp_hours = [int(time_int[0:2]), int(time_int[5:7])]
        comp_minutes = [int(time_int[2:4]), int(time_int[7:9])]

        clock = np.array(comp_hours) + np.divide(comp_minutes, 60.)  # [hours] + [minutes]/60#

        # -- gathering self.year, self.month, self.day for convertion to UTC time
        self.time_int = time_int
        self.year = int('20' + date[:2])
        self.month = int(date[2:4])
        self.day = int(date[4:6])

        time = [0, 0, 0, 0]
        time[0] = datetime.datetime(self.year, self.month, self.day, hour=int(comp_hours[0]),
                                    minute=int(comp_minutes[0]))
        time[1] = time[0] + datetime.timedelta(seconds=15)
        time[3] = datetime.datetime(self.year, self.month, self.day, hour=int(comp_hours[1]),
                                    minute=int(comp_minutes[1]))
        time[2] = time[3] - datetime.timedelta(seconds=15)

        first_file = int(clock[0])
        if clock[1] - int(clock[1]) > 0.0:
            last_file = int(clock[1]) + 1
        else:
            last_file = int(clock[1])

        range_file_list = list(range(first_file, last_file, 1))

        self.ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob('*' + date + '_' + str(il).zfill(2) + '*.LV0.NC'))
            self.ncfiles.append(file_name[2:-2])

            if file_name[2:-2] == '':
                print('   Error!  File: "' + file_name + '" not found --> exit!')
                print('   Check LIMRAD folder!')
                exit(0)

        #        self.ncfiles = []
        #        flat_list = []
        #        for il in range_file_list:
        #            file_name = str(glob.glob('*' + date + '_' + str(il).zfill(2) + '*.LV0.NC'))
        #
        #            flat_list = file_name.split(",")
        #
        #            for iFile in flat_list:
        #                self.ncfiles.append(iFile[2:-1])
        #
        #            if file_name[2:-2] == '':
        #                print('   Error!  File: "' + file_name + '" not found --> exit!')
        #                print('   Check LIMRAD folder!')
        #                exit(0)
        #
        #        n_nc_files = len(self.ncfiles)

        n_nc_files = len(self.ncfiles)

        file = self.ncfiles[0]
        nc_data_set = netCDF4.Dataset(file, 'r')

        # find the number of range gates per chirp sequence,
        # also find the resolution of each chirp and
        # calculate the vector containing the height-steps

        self.nc_VariableList = list(nc_data_set.variables.keys())

        self.TAlt = nc_data_set.dimensions['TAlt'].size
        self.HAlt = nc_data_set.dimensions['HAlt'].size
        self.no_c = nc_data_set.dimensions['Chirp'].size
        self.Time = nc_data_set.dimensions['Time'].size

        self.cum_range_gates = np.zeros((self.no_c + 1,), dtype='int')
        self.range_gates = np.zeros((self.no_c,), dtype='int')
        self.vel_gates = np.zeros((self.no_c,), dtype='int')
        self.range_res = np.array(nc_data_set.variables['RangeRes'])
        self.maxvel = np.array(nc_data_set.variables['MaxVel'])
        self.doplen = np.array(nc_data_set.variables['DoppLen'])
        self.height = [None] * self.no_c  # [chirpTable_min_height]
        self.height_all = [chirpTable_min_height]
        self.n_height = [None] * self.no_c  # [chirpTable_min_height]
        self.latitude = float(nc_data_set.variables['GPSLat'][:])
        self.longitude = float(nc_data_set.variables['GPSLon'][:])
        self.DopplerBins = [None] * self.no_c
        self.CRanges = []

        # "Number of chirps averaged (coherently and non-coherently) for a single time sample"
        self.no_avg = np.array(nc_data_set.variables['AvgNum'][:])

        # From Nils' Code:
        # nAvg = data.SeqAvg./data.DoppLen; % number of spectral averages
        self.no_av = np.divide(self.no_avg, self.doplen)
        # Alexander Myagkov:
        # As you certainly know the standard deviation of the noise power density goes down by sqrt(Nav)
        # in the case of non-coherent averaging and this fact the Hildebrand algorithm takes into account.
        # So you need to use Nav which is AvgNum / DoppLen for the corresponding chirp sequence.

        self.DoppMax = np.array(nc_data_set.variables['MaxVel'][:])
        cnt = 0
        for ichirp in range(self.no_c):
            self.CRanges.append(np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Range'][:]))
            self.range_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Range'].size
            self.vel_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Vel'].size
            self.cum_range_gates[ichirp + 1] = self.cum_range_gates[ichirp] + self.range_gates[ichirp]

            # And the number of Doppler bins
            # Evenly space between the minimum and maximum unambiguous Doppler velocity
            # The number of Doppler bins for each of the chirps
            self.DopplerBins[ichirp] = np.linspace(-self.maxvel[ichirp], self.maxvel[ichirp], self.vel_gates[ichirp])
            self.height[ichirp] = np.zeros(self.range_gates[ichirp])

            for ih in range(self.range_gates[ichirp]):
                self.height_all.append(self.height_all[cnt] + self.range_res[ichirp])
                cnt += 1

                self.height[ichirp][ih] = self.height_all[cnt]/1000.0

            self.n_height[ichirp] = len(self.height[ichirp])

            # for each chirp, get the range (height) dimension
#
#            self.height[ichirp] = np.divide(np.array(range_dummy), 1000.0)
#            self.n_height[ichirp] = len(range_dummy)
        self.height_all.pop(0)
        self.height_all = np.divide(np.array(self.height_all), 1000.)
        self.CRanges = [iR for Cchirp in self.CRanges for iR in Cchirp]
        self.CRanges = np.divide(self.CRanges, 1000.)

        nc_data_set.close()

        self.cum_time_gates = np.zeros((n_nc_files + 1,), dtype='int')
        time_samp = []
        self.CBH  = []
        self.DDTb = []
        self.LWP  = []
        self.Rain = []
        self.SurfPres   = []
        self.SurfRelHum = []
        self.SurfTemp   = []
        self.SurfWD = []
        self.SurfWS = []

        i_nc_file = 0
        for file in self.ncfiles:
            nc_data_set = netCDF4.Dataset(file, 'r')

            time_chirp = np.array(nc_data_set.variables['Time'])  # Number of seconds since 1/1/2001 00:00:00 [UTC]
            time_samp = np.append(time_samp, time_chirp)
            self.cum_time_gates[i_nc_file + 1] = self.cum_time_gates[i_nc_file] + len(time_chirp)

            self.CBH = np.append(self.CBH, np.array(nc_data_set.variables['CBH']))  # Cloud Bottom Height [m]
            self.DDTb = np.append(self.DDTb, np.array(nc_data_set.variables['DDTb']))  # Direct detection brightness temperature [K]
            self.LWP = np.append(self.LWP, np.array(nc_data_set.variables['LWP']))  # Liquid Water Path [g/m2]
            self.Rain = np.append(self.Rain, np.array(nc_data_set.variables['Rain']))  # Rain rate from weather station [mm/h]
            self.SurfPres = np.append(self.SurfPres, np.array(
                nc_data_set.variables['SurfPres']))  # Surface atmospheric pressure from weather station [hPa]
            self.SurfRelHum = np.append(self.SurfRelHum, np.array(
                nc_data_set.variables['SurfRelHum']))  # Relative humidity from weather station [%]
            self.SurfTemp = np.append(self.SurfTemp, np.array(
                nc_data_set.variables['SurfTemp']))  # Surface temperature from weather station [K]
            self.SurfWD = np.append(self.SurfWD, np.array(
                nc_data_set.variables['SurfWD']))  # Surface wind direction from weather station [deg]
            self.SurfWS = np.append(self.SurfWS, np.array(
                nc_data_set.variables['SurfWS']))  # Surface wind speed from weather station [m/s]

            i_nc_file += 1
            nc_data_set.close()

        self.n_time = len(time_samp)

        # gahter radar data values and stack them together
        # Fr_chirps = np.zeros((self.n_time, self.n_height))  # Range factor

        VHSpec_chirps = [None] * self.no_c  # Doppler spectrum at vertical+horizontal polarization [mm6/m3]
        ReVHSpec_chirps = [None] * self.no_c  # Real part of covariance spectrum [mm6/m3]
        ImVHSpec_chirps = [None] * self.no_c  # Imaginary part of covariance spectrum [mm6/m3]
        HSpec_chirps = [None] * self.no_c  # Doppler spectrum at horizontal polarization [mm6/m3 ]
        SLh_chirps = [None] * self.no_c  # Sensitivity limit for horizontal polarization [mm6/m3]
        SLv_chirps = [None] * self.no_c  # Sensitivity limit for vertical polarization [mm6/m3]

        for ichirp in range(self.no_c):
            VHSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
            ReVHSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
            ImVHSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
            HSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
            SLh_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp]))
            SLv_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp]))

        i_nc_file = 0
        n_nc_file = len(self.ncfiles)

        for file in self.ncfiles:
            if pts: print("    Loading LIMRAD94 (LV0.NC) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

            nc_data_set = netCDF4.Dataset(file, 'r')

            lb_t = self.cum_time_gates[i_nc_file]
            ub_t = self.cum_time_gates[i_nc_file + 1]

            for ichirp in range(self.no_c):
                VHSpec_chirps[ichirp][lb_t:ub_t, :, :]   = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'VSpec'])
                ReVHSpec_chirps[ichirp][lb_t:ub_t, :, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ReVHSpec'])
                ImVHSpec_chirps[ichirp][lb_t:ub_t, :, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ImVHSpec'])
                HSpec_chirps[ichirp][lb_t:ub_t, :, :]    = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'HSpec'])
                SLh_chirps[ichirp][lb_t:ub_t, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLh'])
                SLv_chirps[ichirp][lb_t:ub_t, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLv'])

            nc_data_set.close()
            i_nc_file += 1

        if pts: print("    Loading LIMRAD94 (LV0.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format

        time_plot = [datetime.datetime(2001, 1, 1, 0, 0, 0)
                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

        #min_t, max_t = get_time_boundary(time_plot, time)
        min_t, max_t = [0, None]
        min_h, max_h = [0, None]

        self.t_plt = time_plot
        self.t_unix = [ts.replace(tzinfo=timezone.utc).timestamp() for ts in self.t_plt]
        self.n_time = len(self.t_unix)

        # build stacked chirps and prune arrays
        self.CBH = self.CBH
        self.DDTb = self.DDTb
        self.LWP = self.LWP
        self.Rain = self.Rain

        self.SurfPres = self.SurfPres
        self.SurfRelHum = self.SurfRelHum
        self.SurfTemp = self.SurfTemp
        self.SurfWD = self.SurfWD
        self.SurfWS = self.SurfWS

        self.VHSpec = list(block for block in VHSpec_chirps)
        self.ReVHSpec = list(block for block in ReVHSpec_chirps)
        self.ImVHSpec = list(block for block in ImVHSpec_chirps)
        self.HSpec = list(block for block in HSpec_chirps)
        self.SLh = list(block for block in SLh_chirps)
        self.SLv = list(block for block in SLv_chirps)

        # self.VarDict = {'CBH': False, 'DDTb': False, 'LWP': False, 'Rain': False,
        #                'SurfPres': False, 'SurfRelHum': False, 'SurfTemp': False, 'SurfWD': False, 'SurfWS': False,
        #                'Ze': True, 'ZDR': False, 'mdv': True, 'sw': True, 'ldr': True, 'kurt': True, 'Skew': True,
        #                'DiffAtt': False}


    def save(self, path, Ze, mdv, sw, skew, kurt):

        # copy new values to LIMRAD94_lv0 class

        self.Ze   = np.ma.log10(Ze.T)*10.0
        self.mdv  = mdv.T
        self.sw   = sw.T
        self.skew = skew.T
        self.kurt = kurt.T

        ds_name = path + 'concatinated/' + str(self.year) + str(self.month).zfill(2) \
                  + str(self.day).zfill(2) + '_' + self.time_int + '_LIMRAD94_moments.nc'

        ds = netCDF4.Dataset(ds_name, "w", format="NETCDF4")

        ds.description = 'Concatenated data files of LIMRAD 94GHz - FMCW Radar, moments from spectra'
        ds.history = 'Created ' + time.ctime(time.time())
        ds.source = ''

        n_tot_height = sum(self.n_height)
        tot_height   = [iheight for ic_height in self.height for iheight in ic_height]

        # caution ?
        self.height  = tot_height

        ds.createDimension('TAlt',  self.TAlt)
        ds.createDimension('HAlt',  self.HAlt)
        ds.createDimension('Chirp', self.no_c)
        ds.createDimension('time',  self.n_time)
        ds.createDimension('range', n_tot_height)

        for ic in range(self.no_c):
            ds.createDimension('C' + str(ic + 1) + 'Range', self.range_gates[ic])
            ds.createDimension('C' + str(ic + 1) + 'Vel', self.vel_gates[ic])

        self.nc_add_variable(ds, 'time', np.int, ('time',), 'Seconds since 01.01.1970 00:00 UTC', '[sec]', self.t_unix)
        self.nc_add_variable(ds, 'range', np.float32, ('range',), 'range', '[m]', np.copy(np.multiply(1000.0, tot_height)))

        self.nc_add_variable(ds, 'Ze', np.float32, ('time', 'range',), 'Equivalent radar reflectivity factor', '[dBZ]', self.Ze, -999.0)
        self.nc_add_variable(ds, 'vm', np.float32, ('time', 'range',), 'Mean Doppler velocity', '[m/s]', mdv.T, -999.)
        self.nc_add_variable(ds, 'sigma', np.float32, ('time', 'range',), 'Spectrum width', '[m/s]', sw.T, -999.)
        self.nc_add_variable(ds, 'kurt', np.float32, ('time', 'range',), 'Kurtosis', '[linear]', kurt.T, -999.)
        self.nc_add_variable(ds, 'Skew', np.float32, ('time', 'range',), 'Skewness', '[linear]', skew.T, -999.)

        self.nc_add_variable(ds, 'latitude', np.float32, (), 'GPS latitude', '[deg]', self.latitude)
        self.nc_add_variable(ds, 'longitude', np.float32, (), 'GPS longitude', '[deg]', self.longitude)
        self.nc_add_variable(ds, 'DoppMax', np.float32, ('Chirp',), 'Unambiguous Doppler velocity (+/-)', '[m/s]',
                             self.DoppMax)

        self.nc_add_variable(ds, 'cbh', np.float32, ('time',), 'Cloud Bottom Height', '[m]', self.CBH)
        self.nc_add_variable(ds, 'bt', np.float32, ('time',), 'Direct detection brightness temperature', '[m]',
                             self.DDTb)
        self.nc_add_variable(ds, 'lwp', np.float32, ('time',), 'Liquid water path', '[g/m^2]', self.LWP)
        self.nc_add_variable(ds, 'rain', np.float32, ('time',), 'Rain rate from weather station', '[mm/h]',
                             self.Rain)

        self.nc_add_variable(ds, 'SurfPres', np.float32, ('time',),
                             'Surface atmospheric pressure from weather station', '[hPa]', self.SurfPres)
        self.nc_add_variable(ds, 'SurfRelHum', np.float32, ('time',), 'Relative humidity from weather station',
                             '[%]', self.SurfRelHum)
        self.nc_add_variable(ds, 'SurfTemp', np.float32, ('time',), 'Surface temperature from weather station',
                             '[K]', self.SurfTemp)
        self.nc_add_variable(ds, 'SurfWD', np.float32, ('time',), 'Surface wind direction from weather station',
                             '[deg]', self.SurfWD)
        self.nc_add_variable(ds, 'SurfWS', np.float32, ('time',), 'Surface wind speed from weather station',
                             '[deg]', self.SurfWS)

        ds.close()

        print('')
        print('    Concatenated nc file written: ', ds_name)

    @staticmethod
    def nc_add_variable(*args):

        if len(args) < 7:
            print(' check arguments for adding a netCDF variable')

        elif len(args) >= 7:
            datastruct = args[0]
            var_name = args[1]
            type = args[2]
            dim = args[3]
            long_name = args[4]
            unit = args[5]
            data = np.copy(args[6]).T
            fillval = None

            if len(args) == 8:
                fillval = args[7]
                data[data == np.ma.masked] = fillval

        var = datastruct.createVariable(var_name, type, dim, fill_value=fillval)
        var.long_name = long_name
        var.unit = unit
        var[:] = data


class LIMRAD94_LV1():


    frequency = 94.0            # [GHz]
    radar_wavelength = 0.00319  # [m]
    beamwidth = 0.48            # []

    def __init__(self, *args):

        os.chdir(LIMRAD_path)  # path to data needs to be fit to the devices file structure

        if len(args) < 1:
            print('You need to specify a date at least!')
            exit(0)

        elif len(args) == 1:
            date     = args[0]
            time_int = '0000-2400'
            h_bounds = [0.0, 12.0]

        else:
            date     = args[0]
            time_int = args[1]
            h_bounds = args[2]

        comp_hours   = [int(time_int[0:2]), int(time_int[5:7])]
        comp_minutes = [int(time_int[2:4]), int(time_int[7:9])]

        clock = np.array(comp_hours) + np.divide(comp_minutes, 60.)  # [hours] + [minutes]/60#

        # -- gathering self.year, self.month, self.day for convertion to UTC time
        self.time_int = time_int
        self.year = int('20' + date[:2])
        self.month = int(date[2:4])
        self.day = int(date[4:6])

        time = [0, 0, 0, 0]
        time[0] = datetime.datetime(self.year, self.month, self.day, hour=int(comp_hours[0]), minute=int(comp_minutes[0]))
        time[1] = time[0] + datetime.timedelta(seconds=15)
        time[3] = datetime.datetime(self.year, self.month, self.day, hour=int(comp_hours[1]), minute=int(comp_minutes[1]))
        time[2] = time[3] - datetime.timedelta(seconds=15)

        first_file = int(clock[0])
        if clock[1] - int(clock[1]) > 0.0:
            last_file = int(clock[1]) + 1
        else:
            last_file = int(clock[1])

        range_file_list = list(range(first_file, last_file, 1))

        self.ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob('*' + date + '_' + str(il).zfill(2) + '*.LV1.NC'))
            self.ncfiles.append(file_name[2:-2])

            if file_name[2:-2] == '':
                print('   Error!  File: "' + file_name + '" not found --> exit!')
                print('   Check LIMRAD folder!')
                exit(0)

        n_nc_files = len(self.ncfiles)

        file = self.ncfiles[0]
        nc_data_set = netCDF4.Dataset(file, 'r')

        # find the number of range gates per chirp sequence,
        # also find the resolution of each chirp and
        # calculate the vector containing the height-steps

        self.nc_VariableList = list(nc_data_set.variables.keys())

        self.TAlt = nc_data_set.dimensions['TAlt'].size
        self.HAlt = nc_data_set.dimensions['HAlt'].size
        self.no_c = nc_data_set.dimensions['Chirp'].size
        self.Time = nc_data_set.dimensions['Time'].size
        

        self.range_gates = np.zeros((self.no_c,), dtype='int')
        self.vel_gates   = np.zeros((self.no_c,), dtype='int')
        self.cum_range_gates = np.zeros((self.no_c + 1,), dtype='int')
        self.range_res = np.array(nc_data_set.variables['RangeRes'])
        self.height    = [chirpTable_min_height]

        self.latitude  = float(nc_data_set.variables['GPSLat'][:])
        self.longitude = float(nc_data_set.variables['GPSLon'][:])
        #self.azimuth   = nc_data_set.variables['Azm'][:]

        self.DoppMax = np.array(nc_data_set.variables['MaxVel'][:])

        cnt = 0
        for ichirp in range(self.no_c):
            self.range_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Range'].size
            self.vel_gates[ichirp]   = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Vel'].size
            self.cum_range_gates[ichirp + 1] = self.cum_range_gates[ichirp] + self.range_gates[ichirp]

            for _ in range(self.range_gates[ichirp]):
                self.height.append(self.height[cnt] + self.range_res[ichirp])
                cnt += 1

        self.height.pop(0)
        self.height = np.divide(np.array(self.height), 1000.0)
        self.n_height = len(self.height)
        nc_data_set.close()

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
            time_samp  = np.append(time_samp, time_chirp)
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

        Ze_chirps   = np.zeros((self.n_time, self.n_height))        # Equivalent radar reflectivity factor [mm6/m3]
        ZDR_chirps  = np.zeros((self.n_time, self.n_height))        # Differential reflectivity [dB]
        mdv_chirps  = np.zeros((self.n_time, self.n_height))        # Mean Doppler velocity [m/s]
        sw_chirps   = np.zeros((self.n_time, self.n_height))        # Spectrum width [m/s]
        ldr_chirps  = np.zeros((self.n_time, self.n_height))        # Slanted linear depolarization ratio [dB]
        kurt_chirps = np.zeros((self.n_time, self.n_height))        # Kurtosis [linear]
        DiffAtt_chirps = np.zeros((self.n_time, self.n_height))     # Differential attenuation [dBZ/km]
        Skew_chirps = np.zeros((self.n_time, self.n_height))        # Skewness [linear]

        i_nc_file = 0
        n_nc_file = len(self.ncfiles)

        for file in self.ncfiles:
            if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

            nc_data_set = netCDF4.Dataset(file, 'r')

            lb_t = self.cum_time_gates[i_nc_file]
            ub_t = self.cum_time_gates[i_nc_file + 1]

            for ichirp in range(self.no_c):
                lb_h = self.cum_range_gates[ichirp]
                ub_h = self.cum_range_gates[ichirp + 1]
                Ze_chirps[lb_t:ub_t, lb_h:ub_h]   = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZE'])
                ZDR_chirps[lb_t:ub_t, lb_h:ub_h]  = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZDR'])
                mdv_chirps[lb_t:ub_t, lb_h:ub_h]  = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'MeanVel'])
                sw_chirps[lb_t:ub_t, lb_h:ub_h]   = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SpecWidth'])
                ldr_chirps[lb_t:ub_t, lb_h:ub_h]  = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLDR'])
                kurt_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Kurt'])
                DiffAtt_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'DiffAtt'])
                Skew_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Skew'])

            nc_data_set.close()
            i_nc_file += 1

        if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format

        time_plot = [datetime.datetime(2001, 1, 1, 0, 0, 0)
                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

        min_t, max_t = get_time_boundary(time_plot, time)
        min_h, max_h = get_height_boundary(self.height, h_bounds[:])

        # print('min_H, max_H :: ',min_h, max_h, self.height[:5])
        self.t_plt = time_plot[min_t:max_t]
        self.t_unix = [ts.replace(tzinfo=timezone.utc).timestamp() for ts in self.t_plt]
        self.n_time = len(self.t_unix)

        self.height = self.height[min_h:max_h]
        self.n_height = len(self.height)

        # build stacked chirps and prune arrays
        self.CBH  = self.CBH[min_t:max_t]
        self.DDTb = self.DDTb[min_t:max_t]
        self.LWP  = self.LWP[min_t:max_t]
        self.Rain = self.Rain[min_t:max_t]

        self.SurfPres   = self.SurfPres[min_t:max_t]
        self.SurfRelHum = self.SurfRelHum[min_t:max_t]
        self.SurfTemp   = self.SurfTemp[min_t:max_t]
        self.SurfWD = self.SurfWD[min_t:max_t]
        self.SurfWS = self.SurfWS[min_t:max_t]

        self.Ze  = np.ma.log10(np.ma.masked_less_equal(Ze_chirps[min_t:max_t, min_h:max_h].T, 0.)) * 10.0
        self.ZDR = np.ma.masked_less_equal(ZDR_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.mdv = np.ma.masked_less_equal(mdv_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.sw  = np.ma.masked_less_equal(sw_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.ldr = np.ma.masked_less_equal(ldr_chirps[min_t:max_t, min_h:max_h].T, - 999.)
        self.kurt = np.ma.masked_less_equal(kurt_chirps[min_t:max_t, min_h:max_h].T, -999.)   # fill value correct?
        self.Skew = np.ma.masked_less_equal(Skew_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.kurt = np.ma.masked_less_equal(kurt_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.DiffAtt = np.ma.masked_less_equal(DiffAtt_chirps[min_t:max_t, min_h:max_h].T, -999.) # fill value correct?

        self.VarDict = {'CBH': False, 'DDTb': False, 'LWP': False, 'Rain': False,
                        'SurfPres': False, 'SurfRelHum': False, 'SurfTemp': False, 'SurfWD': False, 'SurfWS': False,
                        'Ze': True, 'ZDR': False, 'mdv': True, 'sw': True, 'ldr': True, 'kurt': True, 'Skew': True,
                        'DiffAtt': False}

    def avg_time(self):

        self.timeavg_Ze  = np.average(self.Ze,  axis=1)
        self.timeavg_mdv = np.average(self.mdv, axis=1)
        self.timeavg_sw  = np.average(self.sw,  axis=1)


    def avg_height(self):

        self.heightavg_Ze  = np.average(self.Ze,  axis=0)
        self.heightavg_mdv = np.average(self.mdv, axis=0)
        self.heightavg_sw  = np.average(self.sw,  axis=0)

    def save(self, path):

        ds_name = path + 'concatinated/' + str(self.year) + str(self.month).zfill(2) \
                  + str(self.day).zfill(2) + '_' + self.time_int + '_LIMRAD94.nc'

        ds = netCDF4.Dataset(ds_name, "w", format="NETCDF4")

        ds.description = 'Concatenated data files of LIMRAD 94GHz - FMCW Radar'
        ds.history  = 'Created ' + time.ctime(time.time())
        ds.source   = ''


        ds.createDimension('TAlt', self.TAlt)
        ds.createDimension('HAlt', self.HAlt)
        ds.createDimension('Chirp', self.no_c)
        ds.createDimension('time', len(self.t_unix))
        ds.createDimension('range', len(self.height))

        for ic in range(self.no_c):
            ds.createDimension('C'+str(ic+1)+'Range', self.range_gates[ic])
            ds.createDimension('C'+str(ic+1)+'Vel',   self.vel_gates[ic])


        self.nc_add_variable(ds, 'time', np.int, ('time',),   'Seconds since 01.01.1970 00:00 UTC', '[sec]', self.t_unix)
        self.nc_add_variable(ds, 'range', np.float32, ('range',), 'range', '[m]',
                             np.copy(np.multiply(1000.0, self.height)))

        self.nc_add_variable(ds, 'Ze', np.float32, ('time', 'range',), 'Equivalent radar reflectivity factor', '[dBZ]',
                             self.Ze, -999.)
        self.nc_add_variable(ds, 'vm', np.float32, ('time', 'range',), 'Mean Doppler velocity', '[m/s]', self.mdv,
                             -999.)
        self.nc_add_variable(ds, 'sigma', np.float32, ('time', 'range',), 'Spectrum width', '[m/s]', self.sw, -999.)
        self.nc_add_variable(ds, 'ldr', np.float32, ('time', 'range',), 'Slanted linear depolarization ratio', '[dB]',
                             self.ldr, -999.)
        self.nc_add_variable(ds, 'kurt', np.float32, ('time', 'range',), 'Kurtosis', '[linear]', self.kurt, -999.)
        self.nc_add_variable(ds, 'Skew', np.float32, ('time', 'range',), 'Skewness', '[linear]', self.Skew, -999.)
        self.nc_add_variable(ds, 'DiffAtt', np.float32, ('time', 'range',), 'Differential attenuation', '[dB/km]',
                             self.DiffAtt, -999.)

        self.nc_add_variable(ds, 'latitude',  np.float32, (), 'GPS latitude',  '[deg]', self.latitude)
        self.nc_add_variable(ds, 'longitude', np.float32, (), 'GPS longitude', '[deg]', self.longitude)
        self.nc_add_variable(ds, 'DoppMax', np.float32, ('Chirp',), 'Unambiguous Doppler velocity (+/-)', '[m/s]', self.DoppMax)

        self.nc_add_variable(ds, 'cbh', np.float32, ('time',), 'Cloud Bottom Height', '[m]', self.CBH)
        self.nc_add_variable(ds, 'bt', np.float32, ('time',), 'Direct detection brightness temperature', '[K]',
                             self.DDTb)
        self.nc_add_variable(ds, 'lwp', np.float32, ('time',), 'Liquid water path', '[g/m^2]', self.LWP)
        self.nc_add_variable(ds, 'rain', np.float32, ('time',), 'Rain rate from weather station', '[mm/h]', self.Rain)

        self.nc_add_variable(ds, 'SurfPres',   np.float32, ('time',), 'Surface atmospheric pressure from weather station', '[hPa]', self.SurfPres)
        self.nc_add_variable(ds, 'SurfRelHum', np.float32, ('time',), 'Relative humidity from weather station', '[%]', self.SurfRelHum)
        self.nc_add_variable(ds, 'SurfTemp',   np.float32, ('time',), 'Surface temperature from weather station', '[K]', self.SurfTemp)
        self.nc_add_variable(ds, 'SurfWD',     np.float32, ('time',), 'Surface wind direction from weather station', '[deg]', self.SurfWD)
        self.nc_add_variable(ds, 'SurfWS',     np.float32, ('time',), 'Surface wind speed from weather station', '[deg]', self.SurfWS)

        if interpolate_cn:
            ds.createDimension('time_interp', len(self.t_unix_interp))
            ds.createDimension('range_interp', len(self.height_interp))

            self.nc_add_variable(ds, 'time_interp_res', np.float32, (), 'Time resolution of interpolated data', '[sec]',
                                 interp_time_res)
            self.nc_add_variable(ds, 'range_interp_res', np.float32, (), 'Range resolution of interpolated data', '[m]',
                                 interp_range_res)

            self.nc_add_variable(ds, 'time_interp', np.int, ('time_interp',),
                                 'Seconds since 01.01.1970 00:00 UTC (interpolated)', '[sec]', self.t_unix_interp)
            self.nc_add_variable(ds, 'range_interp', np.float32, ('range_interp',), 'height (interpolated)', '[m]',
                                 np.copy(np.multiply(1000.0, self.height_interp)))

            self.nc_add_variable(ds, 'Ze_interp', np.float32, ('time_interp', 'range_interp',),
                                 'Equivalent radar reflectivity factor (interpolated)', '[dBZ]', self.Ze_interp, -999.)
            self.nc_add_variable(ds, 'vm_interp', np.float32, ('time_interp', 'range_interp',),
                                 'Mean Doppler velocity (interpolated)', '[m/s]', self.mdv_interp, -999.)
            self.nc_add_variable(ds, 'sigma_interp', np.float32, ('time_interp', 'range_interp',),
                                 'Spectrum width (interpolated)', '[m/s]', self.sw_interp, -999.)
            self.nc_add_variable(ds, 'ldr_interp', np.float32, ('time_interp', 'range_interp',),
                                 'Slanted linear depolarization ratio (interpolated)', '[dB]', self.ldr_interp, -999.)

        ds.close()

        print('')
        print('    Concatenated nc file written: ', ds_name)

    @staticmethod
    def nc_add_variable(*args):

        if len(args) < 7:
            print(' check arguments for adding a netCDF variable')

        elif len(args) >= 7:
            datastruct = args[0]
            var_name = args[1]
            type = args[2]
            dim = args[3]
            long_name = args[4]
            unit = args[5]
            data = np.copy(args[6]).T
            fillval = None

            if len(args) == 8:
                fillval = args[7]
                data[data == np.ma.masked] = fillval


        var = datastruct.createVariable(var_name, type, dim, fill_value=fillval)
        var.long_name = long_name
        var.unit = unit
        var[:] = data

    def interpolate_cn(self, t_res, r_res, method):

        self.t_unix_interp = np.arange(self.t_unix[0], self.t_unix[-1], t_res)
        self.height_interp = np.arange(self.height[0], self.height[-1], np.divide(r_res, 1000.0))
        len_t = len(self.t_unix_interp)
        len_h = len(self.height_interp)

        coordinates = np.empty((len_t * len_h, 2))

        cnt = 0
        for iTime in self.t_unix_interp:
            for iRange in self.height_interp:
                coordinates[cnt, 0] = iTime
                coordinates[cnt, 1] = iRange
                cnt += 1

        if method == 'biliniear':
            mth = 'linear'
        else:
            mth = 'constant'

        # for var in self.VarDict:
        #    if var:  interp_z = interpolate2d(self.t_unix, self.height, z1.T, coordinates, mode=mth, bounds_error=False)

        if self.VarDict['Ze']:  interp_Ze = interpolate2d(self.t_unix, self.height, self.Ze.T, coordinates, mode=mth,
                                                          bounds_error=False)
        if self.VarDict['mdv']: interp_mdv = interpolate2d(self.t_unix, self.height, self.mdv.T, coordinates, mode=mth,
                                                           bounds_error=False)
        if self.VarDict['sw']:  interp_sw = interpolate2d(self.t_unix, self.height, self.sw.T, coordinates, mode=mth,
                                                          bounds_error=False)
        if self.VarDict['ldr']: interp_ldr = interpolate2d(self.t_unix, self.height, self.ldr.T, coordinates, mode=mth,
                                                           bounds_error=False)

        # interp_z = np.ma.masked_equal(interp_Ze, 0.0)

        interp_Ze = np.ma.masked_less_equal(interp_Ze, -80.0)
        interp_Ze = np.ma.masked_invalid(interp_Ze)
        interp_Ze = np.reshape(interp_Ze, (len_t, len_h)).T

        interp_mdv = np.ma.masked_less_equal(interp_mdv, -30.0)
        interp_mdv = np.ma.masked_invalid(interp_mdv)
        interp_mdv = np.reshape(interp_mdv, (len_t, len_h)).T

        interp_sw = np.ma.masked_less_equal(interp_sw, -30.0)
        interp_sw = np.ma.masked_invalid(interp_sw)
        interp_sw = np.reshape(interp_sw, (len_t, len_h)).T

        interp_ldr = np.ma.masked_less_equal(interp_ldr, -30.0)
        interp_ldr = np.ma.masked_invalid(interp_ldr)
        interp_ldr = np.reshape(interp_ldr, (len_t, len_h)).T

        self.Ze_interp = interp_Ze
        self.mdv_interp = interp_mdv
        self.sw_interp = interp_sw
        self.ldr_interp = interp_ldr



class MIRA35_LV1():

    def __init__(self, *args):

        if len(args) < 1:
            print('You need to specify a date at least!')
            exit(0)

        elif len(args) == 1:
            date     = args[0]
            time_int = '0000-2400'
            h_bounds = [0.0, 12.0]
            fext     = '*mira.nc'
        elif len(args) < 4:
            date     = args[0]
            time_int = args[1]
            h_bounds = args[2]
            fext     = '*mira.nc'
        else:
            date     = args[0]
            time_int = args[1]
            h_bounds = args[2]
            fext     = args[3]


        comp_hours = [int(time_int[0:2]), int(time_int[5:7])]
        comp_minutes = [int(time_int[2:4]), int(time_int[7:9])]

        clock = np.array(comp_hours) + np.divide(comp_minutes, 60.)  # [hours] + [minutes]/60#

        # -- gathering self.year, self.month, self.day for convertion to UTC time
        self.year  = int('20' + date[:2])
        self.month = int(date[2:4])
        self.day   = int(date[4:6])

        time = [0, 0, 0, 0]
        time[0] = datetime.datetime(self.year, self.month, self.day, hour=int(comp_hours[0]), minute=int(comp_minutes[0]))
        time[1] = time[0] + datetime.timedelta(seconds=15)
        time[3] = datetime.datetime(self.year, self.month, self.day, hour=int(comp_hours[1]), minute=int(comp_minutes[1]))
        time[2] = time[3] - datetime.timedelta(seconds=15)



        if fext == '*mira.nc':
            os.chdir(MIRA_path + 'calibrated/')  # path to data needs to be fit to the devices file structure
            self.ncfiles = glob.glob('20' + date + '*mira.nc')

            #if pts: print("    Loading MIRA35 (mira.nc) NC-files ({} of {})".format(0, 1), end="\r")

            file = self.ncfiles[0]

            nc_data_set = netCDF4.Dataset(file, 'r')

            if file == '':
                print('   Error!  File: "' + file + '" not found --> exit!')
                print('   Check MIRA folder!')
                exit(0)

            height = np.array(nc_data_set.variables['range'])

            # conversion from decimal hour to datetime
            time_samp = np.array(nc_data_set.variables['time'])

            dt_midnight = datetime.datetime(self.year, self.month, self.day)
            time_plot   = [dt_midnight + datetime.timedelta(seconds=int(t * 3600.)) for t in time_samp]

            min_h, max_h = get_height_boundary(height, h_bounds)
            min_t, max_t = get_time_boundary(time_plot, time)

            self.t_unix = time_samp[min_t:max_t]
            self.t_plt  = time_plot[min_t:max_t]
            self.height = height[min_h:max_h]

            self.n_time   = len(self.t_plt)
            self.n_height = len(self.height)

            # --- get data from .nc files ---


            # Calibration convention: in the absence of attenuation, a cloud at 273 K
            # containing one million 100-micron droplets per cubic metre will have a reflectivity of 0 dBZ at all frequencies
            self.Ze  = np.array(nc_data_set.variables['Zh'])        # Radar reflectivity factor [dBZ], Calibrated reflectivity.
            self.mdv = np.array(nc_data_set.variables['v'])         # Doppler velocity [m/s], radial component of the velocity, with positive velocities are away from the radar
            self.sw  = np.array(nc_data_set.variables['width'])     # Spectral width [m/s], standard deviation of the reflectivity-weighted velocities in the radar pulse volume
            self.ldr = np.array(nc_data_set.variables['ldr'])       # Linear depolarisation ratio [dB], ratio of cross-polar to co-polar reflectivity
            self.prf = np.array(nc_data_set.variables['prf'])       # Pulse repetition frequency [Hz]
            self.PulseWidth = np.array(nc_data_set.variables['PulseWidth'])  # Pulse width [s]
            self.frequency  = np.array(nc_data_set.variables['frequency'])   # Radar frequency [GHz]


            # np.warnings.filterwarnings('ignore')

            self.Ze  = self.Ze[min_t:max_t, min_h:max_h]
            self.mdv = self.mdv[min_t:max_t, min_h:max_h]
            self.sw  = self.sw[min_t:max_t, min_h:max_h]
            self.ldr  = self.ldr[min_t:max_t, min_h:max_h]

            # stack variables of individual chirps
            self.Ze = np.ma.masked_less_equal(self.Ze, -999.).T

            # conversion to numpy array for truncation
            self.mdv = np.ma.masked_less_equal(self.mdv, -999.).T

            self.sw = np.ma.masked_invalid(self.sw)
            self.sw = np.ma.masked_less_equal(self.sw, -999.).T

            self.ldr = np.ma.masked_invalid(self.ldr)
            self.ldr = np.ma.masked_less_equal(self.ldr, -999.).T

            self.VarDict = {'CBH': False, 'DDTb': False, 'LWP': False, 'Rain': False,
                            'SurfPres': False, 'SurfRelHum': False, 'SurfTemp': False, 'SurfWD': False, 'SurfWS': False,
                            'Ze': True, 'ZDR': False, 'mdv': True, 'sw': True, 'ldr': True, 'kurt': False,
                            'Skew': False,
                            'DiffAtt': False}

            nc_data_set.close()
            if pts: print("    Loading MIRA35 (mira.nc) NC-files ({} of {})".format(1, 1))

        elif fext == '*.mmclx':

            os.chdir(MIRA_path + 'mmclx/')  # path to data needs to be fit to the devices file structure


            first_file = int(clock[0]) - np.remainder(int(clock[0]), 3)
            if clock[1] - int(clock[1]) > 0.0:
                last_file = int(clock[1]) + 1
            else:
                last_file = int(clock[1])

            range_file_list = list(range(first_file, last_file, 3))

            self.ncfiles = []
            for il in range_file_list:
                file_name = str(glob.glob('20' + date + '_' + str(il).zfill(2) + '*.mmclx'))
                self.ncfiles.append(file_name[2:-2])

            if file_name[2:-2] == '':
                print('   Error!  File not found --> exit!')
                print('   Check LIMRAD folder!')
                exit(0)

            file = self.ncfiles[0]

            if file == '':
                print('   Error!  File: "' + file + '" not found --> exit!')
                print('   Check MIRA folder!')
                exit(0)


            # conversion from decimal hour to datetime
            i_nc_file  = 0
            n_nc_files = len(self.ncfiles)

            time_samp = []
            self.cum_time_gates = np.zeros((n_nc_files + 1,), dtype='int')

            for file in self.ncfiles:

                nc_data_set = netCDF4.Dataset(file, 'r')
                tmp = np.array(nc_data_set.variables['time'])
                time_samp = np.append(time_samp, tmp)                        # Number of seconds since 1/1/2001 00:00:00 [UTC]
                self.cum_time_gates[i_nc_file + 1] = self.cum_time_gates[i_nc_file] + len(tmp)
                i_nc_file += 1

                nc_data_set.close()

            # convert unix to datetime format
            time_plot = [datetime.datetime(1970, 1, 1, 0, 0, 0)
                         + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

            height = np.array(get_nc_data(file, 'range'))
            min_h, max_h = get_height_boundary(height, np.multiply(1000.0, h_bounds[:]))
            min_t, max_t = get_time_boundary(time_plot, time)

            # gahter radar data values and stack them together

            nc_data_set = netCDF4.Dataset(file, 'r')
            self.drg = nc_data_set.variables['drg'][:]
            nc_data_set.close()
            self.t_unix = time_samp
            self.t_plt  = time_plot
            self.height = np.divide(height, 1000.0)

            self.n_time   = len(self.t_unix)
            self.n_height = len(self.height)

            self.Ze        = np.zeros((self.n_time, self.n_height))  # Equivalent radar reflectivity factor Ze of all targets [dBZ]
            self.Ze_hydro  = np.zeros((self.n_time, self.n_height))  # Equivalent Radar Reflectivity Factor Ze of hydrometeors [dBZ]
            self.mdv       = np.zeros((self.n_time, self.n_height))  # Mean Doppler velocity of all targets [m/s]
            self.mdv_hydro = np.zeros((self.n_time, self.n_height))  # Mean Doppler velocity of hydrometeors [m/s]
            self.sw        = np.zeros((self.n_time, self.n_height))  # Spectrum width of all targets [m/s]
            self.sw_hydro  = np.zeros((self.n_time, self.n_height))  # Spectrum width of hydrometeors [m/s]
            self.ldr       = np.zeros((self.n_time, self.n_height))  # linear depolarization ratio of all targets [dB]
            self.ldr_hydro = np.zeros((self.n_time, self.n_height))  # linear depolarization ratio of hydrometeors [dB]
            self.LWC       = np.zeros((self.n_time, self.n_height))  # Liquid Water Content [g/m3]
            self.RainRate  = np.zeros((self.n_time, self.n_height))  # Liquid Water Content [g/m3]


            i_nc_file = 0
            n_nc_file = len(self.ncfiles)

            for file in self.ncfiles:
                if pts: print("    Loading MIRA35 (mmclx) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

                nc_data_set = netCDF4.Dataset(file, 'r')

                lb_t = self.cum_time_gates[i_nc_file]
                ub_t = self.cum_time_gates[i_nc_file + 1]

                self.Ze[lb_t:ub_t, :]  = np.array(nc_data_set.variables['Zg'])
                self.mdv[lb_t:ub_t, :] = np.array(nc_data_set.variables['VELg'])
                self.sw[lb_t:ub_t, :]  = np.array(nc_data_set.variables['RMSg'])
                self.ldr[lb_t:ub_t, :] = np.array(nc_data_set.variables['LDRg'])
                self.LWC[lb_t:ub_t, :] = np.array(nc_data_set.variables['LWC'])

                self.Ze_hydro[lb_t:ub_t, :]  = np.array(nc_data_set.variables['Ze'])
                self.mdv_hydro[lb_t:ub_t, :] = np.array(nc_data_set.variables['VEL'])
                self.sw_hydro[lb_t:ub_t, :]  = np.array(nc_data_set.variables['RMS'])
                self.ldr_hydro[lb_t:ub_t, :] = np.array(nc_data_set.variables['LDR'])
                self.RainRate[lb_t:ub_t, :]  = np.array(nc_data_set.variables['RR'])


                i_nc_file += 1
                nc_data_set.close()

            if pts: print("    Loading MIRA35 (mmclx) NC-files ({} of {})".format(n_nc_file, n_nc_file))


            self.t_unix = time_samp[min_t:max_t]
            self.t_plt  = time_plot[min_t:max_t]
            self.height = np.divide(height[min_h:max_h], 1000.0)

            self.n_time   = len(self.t_plt)
            self.n_height = len(self.height)


            self.Ze  = np.ma.log10(np.ma.masked_invalid(self.Ze[min_t:max_t, min_h:max_h]).T) * 10.
            self.mdv = np.ma.masked_invalid(self.mdv[min_t:max_t, min_h:max_h]).T
            self.sw  = np.ma.masked_invalid(self.sw[min_t:max_t, min_h:max_h]).T
            self.ldr = np.ma.log10(np.ma.masked_invalid(self.ldr[min_t:max_t, min_h:max_h]).T) * 10.

            self.Ze_hydro  = np.ma.log10(np.ma.masked_invalid(self.Ze_hydro[min_t:max_t, min_h:max_h]).T) * 10
            self.mdv_hydro = np.ma.masked_invalid(self.mdv_hydro[min_t:max_t, min_h:max_h]).T
            self.sw_hydro  = np.ma.masked_invalid(self.sw_hydro[min_t:max_t, min_h:max_h]).T
            self.ldr_hydro = np.ma.log10(np.ma.masked_invalid(self.ldr_hydro[min_t:max_t, min_h:max_h]).T) * 10

            self.VarDict = {'CBH': False, 'DDTb': False, 'LWP': False, 'Rain': False,
                            'SurfPres': False, 'SurfRelHum': False, 'SurfTemp': False, 'SurfWD': False, 'SurfWS': False,
                            'Ze': True, 'ZDR': False, 'mdv': True, 'sw': True, 'ldr': True, 'kurt': False,
                            'Skew': False,
                            'DiffAtt': False}

            #self.Ze = np.ma.log10(self.Ze) * 10


    def avg_time(self):

        self.timeavg_Ze  = np.average(self.Ze,  axis=1)
        self.timeavg_mdv = np.average(self.mdv, axis=1)
        self.timeavg_sw  = np.average(self.sw,  axis=1)


    def avg_height(self):

        self.heightavg_Ze  = np.average(self.Ze,  axis=0)
        self.heightavg_mdv = np.average(self.mdv, axis=0)
        self.heightavg_sw  = np.average(self.sw,  axis=0)

    def interpolate_cn(self, t_res, r_res, method):

        self.t_unix_interp = np.arange(self.t_unix[0], self.t_unix[-1], t_res)
        self.height_interp = np.arange(self.height[0], self.height[-1], np.divide(r_res, 1000.0))
        len_t = len(self.t_unix_interp)
        len_h = len(self.height_interp)

        coordinates = np.empty((len_t * len_h, 2))

        cnt = 0
        for iTime in self.t_unix_interp:
            for iRange in self.height_interp:
                coordinates[cnt, 0] = iTime
                coordinates[cnt, 1] = iRange
                cnt += 1

        if method == 'biliniear':
            mth = 'linear'
        else:
            mth = 'constant'

        # for var in self.VarDict:
        #    if var:  interp_z = interpolate2d(self.t_unix, self.height, z1.T, coordinates, mode=mth, bounds_error=False)

        if self.VarDict['Ze']:
            interp_Ze = interpolate2d(self.t_unix, self.height, self.Ze.T, coordinates, mode=mth, bounds_error=False)
            interp_Ze = np.ma.masked_less_equal(interp_Ze, -80.0)
            interp_Ze = np.ma.masked_invalid(interp_Ze)
            interp_Ze = np.reshape(interp_Ze, (len_t, len_h)).T
            self.Ze_interp = interp_Ze

        if self.VarDict['mdv']:
            interp_mdv = interpolate2d(self.t_unix, self.height, self.mdv.T, coordinates, mode=mth, bounds_error=False)
            interp_mdv = np.ma.masked_less_equal(interp_mdv, -30.0)
            interp_mdv = np.ma.masked_invalid(interp_mdv)
            interp_mdv = np.reshape(interp_mdv, (len_t, len_h)).T
            self.mdv_interp = interp_mdv

        if self.VarDict['sw']:
            interp_sw = interpolate2d(self.t_unix, self.height, self.sw.T, coordinates, mode=mth, bounds_error=False)
            interp_sw = np.ma.masked_less_equal(interp_sw, -30.0)
            interp_sw = np.ma.masked_invalid(interp_sw)
            interp_sw = np.reshape(interp_sw, (len_t, len_h)).T
            self.sw_interp = interp_sw

        if self.VarDict['ldr']:
            interp_ldr = interpolate2d(self.t_unix, self.height, self.ldr.T, coordinates, mode=mth, bounds_error=False)
            interp_ldr = np.ma.masked_less_equal(interp_ldr, -30.0)
            interp_ldr = np.ma.masked_invalid(interp_ldr)
            interp_ldr = np.reshape(interp_ldr, (len_t, len_h)).T
            self.ldr_interp = interp_ldr

        # interp_z = np.ma.masked_equal(interp_Ze, 0.0)







def get_nc_data(thisfile, varname):
    # #if pts: print('loading variable '+varname +' from ' + thisfile)
    ncfile = netCDF4.Dataset(thisfile, 'r')
    var = ncfile.variables[varname][:]

    if ncfile.isopen == 1: ncfile.close()
    return var


def get_nc_date(thisfile):
    # #if pts: print('loading variable '+varname +' from ' + thisfile)
    ncfile = netCDF4.Dataset(thisfile, 'r')
    year = ncfile.self.year
    month = ncfile.self.month
    day = ncfile.self.day

    if ncfile.isopen == 1: ncfile.close()
    return year, month, day


def get_nc_dimension(thisfile, dim_name):
    ncfile = netCDF4.Dataset(thisfile, 'r')
    dim = ncfile.dimensions[dim_name].size
    if ncfile.isopen == 1: ncfile.close()
    return dim


def get_height_boundary(Array, h_int):
    i = 0
    imin = 0
    for h in Array:
        # print(' (min) h = ',h,hmin)
        if (h_int[0] <= h):
            imin = i
            break
        i += 1

    i = 0
    imax = None
    for h in reversed(Array):
        # print(' (max) h = ',h,hmax)
        if h_int[1] >= h:
            imax = len(Array) - i
            if i == 0: imax = None
            break
        i += 1

    return imin, imax

def get_time_boundary(Array, t_int):
    i = 0
    t_min = 0
    t_max = None
    for zeit in Array:
        if t_int[0] <= zeit <= t_int[1]: t_min = i
        if t_int[2] <= zeit <= t_int[3]: t_max = i
        i += 1

    return t_min, t_max


