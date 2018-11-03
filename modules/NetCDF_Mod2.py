import datetime
import glob
import os
import os.path
import sys
import time
from datetime import timezone

import netCDF4
import numpy as np

from modules.Interpolation_Mod import interpolate2d
from modules.Parameter_Mod import *

# fixed paramters
max_MDF_files = 10


class LIMRAD94():
    frequency = 94.0  # [GHz]
    beamwidth = 0.48  # []
    radar_wavelength = 0.00319  # [m]

    def __init__(self, *args):

        os.chdir(LIMRAD_path)  # path to data needs to be fit to the devices file structure

        try:
            # check input parameter
            if len(args) < 1:
                print('You need to specify a date at least!')
                exit(0)

            # if one argument is given it contains the path to one specific file
            elif len(args) == 1:
                file_path = args[0]
                file_str = file_path[file_path.rfind('/') + 1:]

                self.ncfiles = file_path
                self.n_files = 1

                # search for the position where the program number appears
                pos_prog = file_str.find('_P')

                # gathering self.year, self.month, self.day for convertion to UTC time
                self.year = int('20' + file_str[pos_prog - 13:pos_prog - 11])
                self.month = int(file_str[pos_prog - 11:pos_prog - 9])
                self.day = int(file_str[pos_prog - 9:pos_prog - 7])

                self.time_int = file_str[pos_prog - 6:pos_prog] + '-xxxxxx'
                self.lvl = file_str[-6:-3]

                self.file_MDF = [[file_path]]
                self.name_MDF = file_str[pos_prog + 1:pos_prog + 4]
                self.num_MDF = [1]

            # if there are four values given, the args contain the:
            #   - path to the data folder (string):                     -> args[0]
            #   - date in format (string):              YYMMDD          -> args[1]
            #   - time interval of the date (string):   HHMMSS-HHMMSS   -> args[2]
            #   - height range from/to in km (2*float): [h_min, h_max]  -> args[3]
            elif len(args) > 3:

                folder_path = args[0]
                date_str = args[1]
                time_str = args[2]
                height_str = args[3]

                # gathering self.year, self.month, self.day for convertion to UTC time
                self.time_int = time_str
                self.year = int('20' + date_str[:2])
                self.month = int(date_str[2:4])
                self.day = int(date_str[4:6])

                # check path for LV0 or LV1 files, if subfolder /LVx/ (x=0 or x=1) not existent raise error
                try:
                    pos_lvl = folder_path.find('LV')
                    if pos_lvl < 0: raise Exception('Folder path does not contain LVx information (x=0 or x=1)!')
                except Exception as e:
                    print('Something went wrong:', e)
                    print('Change to: [path_to_data]/YYMMDD/LVx/')
                    sys.exit()
                else:
                    self.lvl = folder_path[pos_lvl:pos_lvl + 3]

                # count LVx files in the given folder
                files_path = folder_path + '*' + self.lvl + '.NC'
                all_ncfiles = [name for name in glob.glob(files_path)]

                # Save only the files which are in between the time_int boundaries
                ncfiles = []
                for iFile in all_ncfiles:
                    pos_time = iFile.find('_P') - 6  # search for the position where the program number appears
                    if int(self.time_int[:6]) <= int(iFile[pos_time:pos_time + 6]) <= int(self.time_int[7:]):
                        ncfiles.append(iFile)

                self.n_files = len(ncfiles)

                # sort the list of netcdf files
                pos_time = all_ncfiles[0].find('_P') - 6
                only_times = [None] * len(all_ncfiles)
                for itime in range(len(all_ncfiles)):
                    only_times[itime] = int(all_ncfiles[itime][pos_time:pos_time + 6])

                permutation = np.argsort(np.array(only_times))

                self.ncfiles = []
                for ifile in range(self.n_files):
                    self.ncfiles.append(ncfiles[permutation[ifile]])

                # Check data for different measurement definition files (MDF), count different MDFs
                # max number of definition files = max_MDF_files. If the user provides different chirp
                # tables for one day, one should know that time and height
                # resolution of each chip program (MDF) is likely to be different.

                file_MDF = [None] * max_MDF_files
                num_MDF = [None] * max_MDF_files
                name_MDF = [None] * max_MDF_files
                for iMDF in range(max_MDF_files):

                    file_this_MDF = glob.glob(folder_path + '*_P' + str(iMDF + 1).zfill(2) + '_*')
                    num_this_MDF = len(file_this_MDF)

                    if num_this_MDF > 0:
                        file_MDF[iMDF] = file_this_MDF
                        num_MDF[iMDF] = num_this_MDF
                        name_MDF[iMDF] = 'P' + str(iMDF + 1).zfill(2)

                # remove entries with None value
                self.file_MDF = [fMDF for fMDF in file_MDF if fMDF is not None]
                self.name_MDF = [MDF for MDF in name_MDF if MDF is not None]
                self.num_MDF = [iMDF for iMDF in num_MDF if iMDF is not None]

                # setting flag if multiple MDFs are given
                if len(self.num_MDF) > 1:
                    self.multiple_MDFs = True
                else:
                    self.multiple_MDFs = False

            else:
                raise Exception(' Input routine')

        except Exception as e:
            print('Something went wrong :(    ' + str(e))
            sys.exit(0)

        # extract the dimensionality information from all different MDFs and store it in a list of dictionaries
        self.dimensions = [None] * len(self.num_MDF)
        self.num_chirps = [None] * len(self.num_MDF)
        self.range_gates = [None] * len(self.num_MDF)

        for iMDF in range(len(self.num_MDF)):
            nc_data_set = netCDF4.Dataset(self.file_MDF[iMDF][0], 'r')
            self.num_chirps[iMDF] = nc_data_set.dimensions['Chirp'].size

            Nav = np.divide(np.array(nc_data_set.variables['AvgNum'][:]),
                            np.array(nc_data_set.variables['DoppLen'][:]))

            self.dimensions[iMDF] = {
                'Program': self.name_MDF[iMDF],
                'NumChirps': self.num_chirps[iMDF],
                'RangeRes': list(nc_data_set.variables['RangeRes'][:]),
                'MaxVel': list(nc_data_set.variables['MaxVel'][:]),
                'AvgNum': list(nc_data_set.variables['AvgNum'][:]),
                'DoppLen': list(nc_data_set.variables['DoppLen'][:]),
                'Nav': list(Nav),
                'Range': [nc_data_set.dimensions['C' + str(iC + 1) + 'Range'].size for iC in
                          range(self.num_chirps[iMDF])],
                'Vel': [nc_data_set.dimensions['C' + str(iC + 1) + 'Vel'].size for iC in range(self.num_chirps[iMDF])]
            }

            self.latitude = float(nc_data_set.variables['GPSLat'][:])
            self.longitude = float(nc_data_set.variables['GPSLon'][:])

            # Alexander Myagkov: (regarding Nav values above)
            # As you certainly know the standard deviation of the noise power density goes down by sqrt(Nav)
            # in the case of non-coherent averaging and this fact the Hildebrand algorithm takes into account.
            # So you need to use Nav which is AvgNum / DoppLen for the corresponding chirp sequence.

            # range bins and velocity bins are collected a view line above, time bins will be extracted later
            for idim in nc_data_set.dimensions:
                if str(idim).find('C') < 0 or str(idim).find('Time') < 0:
                    self.dimensions[iMDF].update({str(idim): nc_data_set.dimensions[idim].size})

            nc_data_set.close()

        # separate different variables by dimension
        # use an arbitrary file because all files contain the same variables (LV0 OR LV1)

        self.var_names_1D = [None] * len(self.num_MDF)
        self.var_names_2D = [None] * len(self.num_MDF)
        self.var_names_3D = [None] * len(self.num_MDF)

        max_chirp_num = max(self.num_chirps)
        i_nc_file = 0
        n_nc_file = len(self.ncfiles)

        for iMDF in range(len(self.num_MDF)):
            self.var_names_1D[iMDF] = dict()
            self.var_names_2D[iMDF] = dict()
            self.var_names_3D[iMDF] = dict()

            for iFile in range(self.num_MDF[iMDF]):

                if pts: print("    Loading LIMRAD94 NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")
                i_nc_file += 1

                nc_data_set = netCDF4.Dataset(self.file_MDF[iMDF][iFile], 'r')
                for ivar in nc_data_set.variables.keys():
                    var = nc_data_set.variables[ivar]
                    if 'Units' in var.ncattrs():  # these variables have Units

                        if iFile == 0:
                            if len(var.shape) == 1:
                                if var.shape[0] > max_chirp_num:
                                    self.var_names_1D[iMDF].update({ivar:
                                                                        {'Dim': [var.shape], 'LongName': var.Name,
                                                                         'Unit': var.Units, 'Val': [var[:]]}})
                            elif len(var.shape) == 2:
                                self.var_names_2D[iMDF].update({ivar:
                                                                    {'Dim': [var.shape], 'LongName': var.Name,
                                                                     'Unit': var.Units, 'Val': [var[:, :]]}})
                            elif len(var.shape) == 3:
                                self.var_names_3D[iMDF].update({ivar:
                                                                    {'Dim': [var.shape], 'LongName': var.Name,
                                                                     'Unit': var.Units, 'Val': [var[:, :, :]]}})
                        else:
                            if len(var.shape) == 1:
                                if var.shape[0] > max_chirp_num:
                                    self.var_names_1D[iMDF][ivar]['Dim'].append(var.shape)
                                    self.var_names_1D[iMDF][ivar]['Val'].append(var[:])
                            elif len(var.shape) == 2:
                                self.var_names_2D[iMDF][ivar]['Dim'].append(var.shape)
                                self.var_names_2D[iMDF][ivar]['Val'].append(var[:, :])
                            elif len(var.shape) == 3:
                                self.var_names_3D[iMDF][ivar]['Dim'].append(var.shape)
                                self.var_names_3D[iMDF][ivar]['Val'].append(var[:, :, :])

                nc_data_set.close()

        if pts: print("    Loading LIMRAD94 NC-files ({} of {})".format(n_nc_file, n_nc_file))

        #
        #        # time_MDF
        ranges = [None] * self.num_chirps[iMDF]

    #        range_all = [chirpTable_min_height]
    #        num_ranges = [None] * self.num_chirps[iMDF]
    #        DopplerBins = [None] * self.num_chirps[iMDF]
    #
    #        self.CRanges = []
    #
    #        # "Number of chirps averaged (coherently and non-coherently) for a single time sample"
    #        self.no_avg = np.array(nc_data_set.variables['AvgNum'][:])
    #
    #        # From Nils' Code:
    #        # nAvg = data.SeqAvg./data.DoppLen; % number of spectral averages
    #
    #        self.DoppMax = np.array(nc_data_set.variables['MaxVel'][:])
    #        cnt = 0
    #        for ichirp in range(self.no_c):
    #            self.CRanges.append(np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Range'][:]))
    #            self.range_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Range'].size
    #            self.vel_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Vel'].size
    #            self.cum_range_gates[ichirp + 1] = self.cum_range_gates[ichirp] + self.range_gates[ichirp]
    #
    #            # And the number of Doppler bins
    #            # Evenly space between the minimum and maximum unambiguous Doppler velocity
    #            # The number of Doppler bins for each of the chirps
    #            self.DopplerBins[ichirp] = np.linspace(-self.maxvel[ichirp], self.maxvel[ichirp], self.vel_gates[ichirp])
    #            self.height[ichirp] = np.zeros(self.range_gates[ichirp])
    #
    #            for ih in range(self.range_gates[ichirp]):
    #                self.height_all.append(self.height_all[cnt] + self.range_res[ichirp])
    #                cnt += 1
    #
    #                self.height[ichirp][ih] = self.height_all[cnt] / 1000.0
    #
    #            self.n_height[ichirp] = len(self.height[ichirp])
    #
    #            # for each chirp, get the range (height) dimension
    #        #
    #        #            self.height[ichirp] = np.divide(np.array(range_dummy), 1000.0)
    #        #            self.n_height[ichirp] = len(range_dummy)
    #        self.height_all.pop(0)
    #        self.height_all = np.divide(np.array(self.height_all), 1000.)
    #        self.CRanges = [iR for Cchirp in self.CRanges for iR in Cchirp]
    #        self.CRanges = np.divide(self.CRanges, 1000.)
    #
    #        nc_data_set.close()
    #
    #        self.cum_time_gates = np.zeros((n_nc_files + 1,), dtype='int')
    #
    #        self.n_time = len(time_samp)
    #
    #        # gahter radar data values and stack them together
    #        # Fr_chirps = np.zeros((self.n_time, self.n_height))  # Range factor
    #
    #        VHSpec_chirps = [None] * self.no_c  # Doppler spectrum at vertical+horizontal polarization [mm6/m3]
    #        ReVHSpec_chirps = [None] * self.no_c  # Real part of covariance spectrum [mm6/m3]
    #        ImVHSpec_chirps = [None] * self.no_c  # Imaginary part of covariance spectrum [mm6/m3]
    #        HSpec_chirps = [None] * self.no_c  # Doppler spectrum at horizontal polarization [mm6/m3 ]
    #        SLh_chirps = [None] * self.no_c  # Sensitivity limit for horizontal polarization [mm6/m3]
    #        SLv_chirps = [None] * self.no_c  # Sensitivity limit for vertical polarization [mm6/m3]
    #
    #        for ichirp in range(self.no_c):
    #            VHSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
    #            ReVHSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
    #            ImVHSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
    #            HSpec_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp], self.vel_gates[ichirp]))
    #            SLh_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp]))
    #            SLv_chirps[ichirp] = np.zeros((self.n_time, self.n_height[ichirp]))
    #
    #        i_nc_file = 0
    #        n_nc_file = len(self.ncfiles)
    #
    #        for file in self.ncfiles:
    #            if pts: print("    Loading LIMRAD94 (LV0.NC) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")
    #
    #            nc_data_set = netCDF4.Dataset(file, 'r')
    #
    #            lb_t = self.cum_time_gates[i_nc_file]
    #            ub_t = self.cum_time_gates[i_nc_file + 1]
    #
    #            for ichirp in range(self.no_c):
    #                VHSpec_chirps[ichirp][lb_t:ub_t, :, :] = np.array(
    #                    nc_data_set.variables['C' + str(ichirp + 1) + 'VSpec'])
    #                ReVHSpec_chirps[ichirp][lb_t:ub_t, :, :] = np.array(
    #                    nc_data_set.variables['C' + str(ichirp + 1) + 'ReVHSpec'])
    #                ImVHSpec_chirps[ichirp][lb_t:ub_t, :, :] = np.array(
    #                    nc_data_set.variables['C' + str(ichirp + 1) + 'ImVHSpec'])
    #                HSpec_chirps[ichirp][lb_t:ub_t, :, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'HSpec'])
    #                SLh_chirps[ichirp][lb_t:ub_t, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLh'])
    #                SLv_chirps[ichirp][lb_t:ub_t, :] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLv'])
    #
    #            nc_data_set.close()
    #            i_nc_file += 1
    #
    #        if pts: print("    Loading LIMRAD94 (LV0.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))
    #
    #        # convert times in datetime format
    #
    #        time_plot = [datetime.datetime(2001, 1, 1, 0, 0, 0)
    #                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]
    #
    #        # min_t, max_t = get_time_boundary(time_plot, time)
    #        min_t, max_t = [0, None]
    #        min_h, max_h = [0, None]
    #
    #        self.t_plt = time_plot
    #        self.t_unix = [ts.replace(tzinfo=timezone.utc).timestamp() for ts in self.t_plt]
    #        self.n_time = len(self.t_unix)
    #
    #        # build stacked chirps and prune arrays
    #        self.CBH = self.CBH
    #        self.DDTb = self.DDTb
    #        self.LWP = self.LWP
    #        self.Rain = self.Rain
    #
    #        self.SurfPres = self.SurfPres
    #        self.SurfRelHum = self.SurfRelHum
    #        self.SurfTemp = self.SurfTemp
    #        self.SurfWD = self.SurfWD
    #        self.SurfWS = self.SurfWS
    #
    #        self.VHSpec = list(block for block in VHSpec_chirps)
    #        self.ReVHSpec = list(block for block in ReVHSpec_chirps)
    #        self.ImVHSpec = list(block for block in ImVHSpec_chirps)
    #        self.HSpec = list(block for block in HSpec_chirps)
    #        self.SLh = list(block for block in SLh_chirps)
    #        self.SLv = list(block for block in SLv_chirps)

    def save(self, path, Ze, mdv, sw, skew, kurt):

        # copy new values to LIMRAD94_lv0 class

        self.Ze = np.ma.log10(Ze.T) * 10.0
        self.mdv = mdv.T
        self.sw = sw.T
        self.skew = skew.T
        self.kurt = kurt.T

        ds_name = path + 'concatinated/' + str(self.year) + str(self.month).zfill(2) \
                  + str(self.day).zfill(2) + '_' + self.time_int + '_LIMRAD94_moments.nc'

        ds = netCDF4.Dataset(ds_name, "w", format="NETCDF4")

        ds.description = 'Concatenated data files of LIMRAD 94GHz - FMCW Radar, moments from spectra'
        ds.history = 'Created ' + time.ctime(time.time())
        ds.source = ''

        n_tot_height = sum(self.n_height)
        tot_height = [iheight for ic_height in self.height for iheight in ic_height]

        # caution ?
        self.height = tot_height

        ds.createDimension('TAlt', self.TAlt)
        ds.createDimension('HAlt', self.HAlt)
        ds.createDimension('Chirp', self.no_c)
        ds.createDimension('time', self.n_time)
        ds.createDimension('range', n_tot_height)

        for ic in range(self.no_c):
            ds.createDimension('C' + str(ic + 1) + 'Range', self.range_gates[ic])
            ds.createDimension('C' + str(ic + 1) + 'Vel', self.vel_gates[ic])

        self.nc_add_variable(ds, 'time', np.int, ('time',), 'Seconds since 01.01.1970 00:00 UTC', '[sec]', self.t_unix)
        self.nc_add_variable(ds, 'range', np.float32, ('range',), 'range', '[m]',
                             np.copy(np.multiply(1000.0, tot_height)))

        self.nc_add_variable(ds, 'Ze', np.float32, ('time', 'range',), 'Equivalent radar reflectivity factor', '[dBZ]',
                             self.Ze, -999.0)
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

        n_files = list(range(first_file, last_file, 1))

        self.ncfiles = []
        for il in n_files:
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
        self.vel_gates = np.zeros((self.no_c,), dtype='int')
        self.cum_range_gates = np.zeros((self.no_c + 1,), dtype='int')
        self.range_res = np.array(nc_data_set.variables['RangeRes'])
        self.height = [chirpTable_min_height]

        self.latitude = float(nc_data_set.variables['GPSLat'][:])
        self.longitude = float(nc_data_set.variables['GPSLon'][:])
        # self.azimuth   = nc_data_set.variables['Azm'][:]

        self.DoppMax = np.array(nc_data_set.variables['MaxVel'][:])

        cnt = 0
        for ichirp in range(self.no_c):
            self.range_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Range'].size
            self.vel_gates[ichirp] = nc_data_set.dimensions['C' + str(ichirp + 1) + 'Vel'].size
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
        self.CBH = []
        self.DDTb = []
        self.LWP = []
        self.Rain = []
        self.SurfPres = []
        self.SurfRelHum = []
        self.SurfTemp = []
        self.SurfWD = []
        self.SurfWS = []

        i_nc_file = 0
        for file in self.ncfiles:
            nc_data_set = netCDF4.Dataset(file, 'r')

            time_chirp = np.array(nc_data_set.variables['Time'])  # Number of seconds since 1/1/2001 00:00:00 [UTC]
            time_samp = np.append(time_samp, time_chirp)
            self.cum_time_gates[i_nc_file + 1] = self.cum_time_gates[i_nc_file] + len(time_chirp)

            self.CBH = np.append(self.CBH, np.array(nc_data_set.variables['CBH']))  # Cloud Bottom Height [m]
            self.DDTb = np.append(self.DDTb, np.array(
                nc_data_set.variables['DDTb']))  # Direct detection brightness temperature [K]
            self.LWP = np.append(self.LWP, np.array(nc_data_set.variables['LWP']))  # Liquid Water Path [g/m2]
            self.Rain = np.append(self.Rain,
                                  np.array(nc_data_set.variables['Rain']))  # Rain rate from weather station [mm/h]

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

        self.n_height = len(self.height)
        self.n_time = len(time_samp)

        # gahter radar data values and stack them together

        Ze_chirps = np.zeros((self.n_time, self.n_height))  # Equivalent radar reflectivity factor [mm6/m3]
        ZDR_chirps = np.zeros((self.n_time, self.n_height))  # Differential reflectivity [dB]
        mdv_chirps = np.zeros((self.n_time, self.n_height))  # Mean Doppler velocity [m/s]
        sw_chirps = np.zeros((self.n_time, self.n_height))  # Spectrum width [m/s]
        ldr_chirps = np.zeros((self.n_time, self.n_height))  # Slanted linear depolarization ratio [dB]
        kurt_chirps = np.zeros((self.n_time, self.n_height))  # Kurtosis [linear]
        DiffAtt_chirps = np.zeros((self.n_time, self.n_height))  # Differential attenuation [dBZ/km]
        Skew_chirps = np.zeros((self.n_time, self.n_height))  # Skewness [linear]

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
                Ze_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZE'])
                ZDR_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'ZDR'])
                mdv_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'MeanVel'])
                sw_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SpecWidth'])
                ldr_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'SLDR'])
                kurt_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Kurt'])
                DiffAtt_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(
                    nc_data_set.variables['C' + str(ichirp + 1) + 'DiffAtt'])
                Skew_chirps[lb_t:ub_t, lb_h:ub_h] = np.array(nc_data_set.variables['C' + str(ichirp + 1) + 'Skew'])

            nc_data_set.close()
            i_nc_file += 1

        if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format

        time_plot = [datetime.datetime(2001, 1, 1, 0, 0, 0)
                     + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp))]

        min_t, max_t = get_time_boundary(time_plot, time)
        min_h, max_h = get_height_boundary(self.height, np.multiply(1000.0, h_bounds[:]))

        self.t_plt = time_plot[min_t:max_t]
        self.t_unix = [ts.replace(tzinfo=timezone.utc).timestamp() for ts in self.t_plt]
        self.n_time = len(self.t_unix)

        self.height = self.height[min_h:max_h]
        self.n_height = len(self.height)

        # build stacked chirps and prune arrays
        self.CBH = self.CBH[min_t:max_t]
        self.DDTb = self.DDTb[min_t:max_t]
        self.LWP = self.LWP[min_t:max_t]
        self.Rain = self.Rain[min_t:max_t]

        self.SurfPres = self.SurfPres[min_t:max_t]
        self.SurfRelHum = self.SurfRelHum[min_t:max_t]
        self.SurfTemp = self.SurfTemp[min_t:max_t]
        self.SurfWD = self.SurfWD[min_t:max_t]
        self.SurfWS = self.SurfWS[min_t:max_t]

        self.Ze = np.ma.log10(np.ma.masked_less_equal(Ze_chirps[min_t:max_t, min_h:max_h].T, 0.)) * 10.0
        self.ZDR = np.ma.masked_less_equal(ZDR_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.mdv = np.ma.masked_less_equal(mdv_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.sw = np.ma.masked_less_equal(sw_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.ldr = np.ma.masked_less_equal(ldr_chirps[min_t:max_t, min_h:max_h].T, - 999.)
        self.kurt = np.ma.masked_less_equal(kurt_chirps[min_t:max_t, min_h:max_h].T, -999.)  # fill value correct?
        self.Skew = np.ma.masked_less_equal(Skew_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.kurt = np.ma.masked_less_equal(kurt_chirps[min_t:max_t, min_h:max_h].T, -999.)
        self.DiffAtt = np.ma.masked_less_equal(DiffAtt_chirps[min_t:max_t, min_h:max_h].T, -999.)  # fill value correct?

        self.VarDict = {'CBH': False, 'DDTb': False, 'LWP': False, 'Rain': False,
                        'SurfPres': False, 'SurfRelHum': False, 'SurfTemp': False, 'SurfWD': False, 'SurfWS': False,
                        'Ze': True, 'ZDR': False, 'mdv': True, 'sw': True, 'ldr': True, 'kurt': True, 'Skew': True,
                        'DiffAtt': False}

    def avg_time(self):

        self.timeavg_Ze = np.average(self.Ze, axis=1)
        self.timeavg_mdv = np.average(self.mdv, axis=1)
        self.timeavg_sw = np.average(self.sw, axis=1)

    def avg_height(self):

        self.heightavg_Ze = np.average(self.Ze, axis=0)
        self.heightavg_mdv = np.average(self.mdv, axis=0)
        self.heightavg_sw = np.average(self.sw, axis=0)

    def save(self, path):

        ds_name = path + 'concatinated/' + str(self.year) + str(self.month).zfill(2) \
                  + str(self.day).zfill(2) + '_' + self.time_int + '_LIMRAD94.nc'

        ds = netCDF4.Dataset(ds_name, "w", format="NETCDF4")

        ds.description = 'Concatenated data files of LIMRAD 94GHz - FMCW Radar'
        ds.history = 'Created ' + time.ctime(time.time())
        ds.source = ''

        ds.createDimension('TAlt', self.TAlt)
        ds.createDimension('HAlt', self.HAlt)
        ds.createDimension('Chirp', self.no_c)
        ds.createDimension('time', len(self.t_unix))
        ds.createDimension('range', len(self.height))

        for ic in range(self.no_c):
            ds.createDimension('C' + str(ic + 1) + 'Range', self.range_gates[ic])
            ds.createDimension('C' + str(ic + 1) + 'Vel', self.vel_gates[ic])

        self.nc_add_variable(ds, 'time', np.int, ('time',), 'Seconds since 01.01.1970 00:00 UTC', '[sec]', self.t_unix)
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

        self.nc_add_variable(ds, 'latitude', np.float32, (), 'GPS latitude', '[deg]', self.latitude)
        self.nc_add_variable(ds, 'longitude', np.float32, (), 'GPS longitude', '[deg]', self.longitude)
        self.nc_add_variable(ds, 'DoppMax', np.float32, ('Chirp',), 'Unambiguous Doppler velocity (+/-)', '[m/s]',
                             self.DoppMax)

        self.nc_add_variable(ds, 'cbh', np.float32, ('time',), 'Cloud Bottom Height', '[m]', self.CBH)
        self.nc_add_variable(ds, 'bt', np.float32, ('time',), 'Direct detection brightness temperature', '[K]',
                             self.DDTb)
        self.nc_add_variable(ds, 'lwp', np.float32, ('time',), 'Liquid water path', '[g/m^2]', self.LWP)
        self.nc_add_variable(ds, 'rain', np.float32, ('time',), 'Rain rate from weather station', '[mm/h]', self.Rain)

        self.nc_add_variable(ds, 'SurfPres', np.float32, ('time',), 'Surface atmospheric pressure from weather station',
                             '[hPa]', self.SurfPres)
        self.nc_add_variable(ds, 'SurfRelHum', np.float32, ('time',), 'Relative humidity from weather station', '[%]',
                             self.SurfRelHum)
        self.nc_add_variable(ds, 'SurfTemp', np.float32, ('time',), 'Surface temperature from weather station', '[K]',
                             self.SurfTemp)
        self.nc_add_variable(ds, 'SurfWD', np.float32, ('time',), 'Surface wind direction from weather station',
                             '[deg]', self.SurfWD)
        self.nc_add_variable(ds, 'SurfWS', np.float32, ('time',), 'Surface wind speed from weather station', '[deg]',
                             self.SurfWS)

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
