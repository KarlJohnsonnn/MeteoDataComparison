import glob
import os
import os.path
import re
import sys

import netCDF4
import numpy as np

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
                heightminmax = args[3]

                # gathering self.year, self.month, self.day for convertion to UTC time
                self.time_int = time_str
                self.year = int('20' + date_str[:2])
                self.month = int(date_str[2:4])
                self.day = int(date_str[4:6])

                # set minimum and maximum height (for plotting)
                self.h_min = heightminmax[0]
                self.h_max = heightminmax[1]

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
                files_path = folder_path + '*' + date_str + '*' + self.lvl + '.NC'
                all_ncfiles = [name for name in glob.glob(files_path)]

                # Save only the files which are in between the time_int boundaries
                ncfiles = []
                for iFile in all_ncfiles:
                    pos_time = iFile.find('_P') - 6  # search for the position where the program number appears
                    if int(self.time_int[:6]) <= int(iFile[pos_time:pos_time + 6]) <= int(self.time_int[7:]):
                        ncfiles.append(iFile)

                self.n_files = len(ncfiles)

                if self.n_files == 0: print('   No Files found.')

                # sort the list of netcdf files
                pos_time = ncfiles[0].find('_P') - 6
                only_times = [None] * len(ncfiles)
                for itime in range(len(ncfiles)):
                    only_times[itime] = int(ncfiles[itime][pos_time:pos_time + 6])

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

                    tmp_files = []

                    for File in self.ncfiles:
                        if File.find('_P' + str(iMDF + 1).zfill(2)) > 0:
                            tmp_files.append(File)

                    # file_this_MDF = glob.glob(folder_path + '*' + date_str + '*' + '*_P' + str(iMDF + 1).zfill(2) + '_*')
                    file_this_MDF = tmp_files
                    num_this_MDF = len(file_this_MDF)

                    if num_this_MDF > 0:
                        num_MDF[iMDF] = num_this_MDF
                        name_MDF[iMDF] = 'P' + str(iMDF + 1).zfill(2)

                        # sort the list of netcdf files
                        pos_time = file_this_MDF[0].find('_P') - 6
                        only_times = [None] * num_this_MDF
                        for itime in range(len(file_this_MDF)):
                            # only_times[itime] = int(all_ncfiles[itime][pos_time:pos_time + 6])
                            only_times[itime] = int(file_this_MDF[itime][pos_time:pos_time + 6])

                        permutation = np.argsort(np.array(only_times))
                        file_MDF[iMDF] = []
                        for ifile in range(num_MDF[iMDF]):
                            file_MDF[iMDF].append(file_this_MDF[permutation[ifile]])

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
            print('Something went wrong X-(    ' + str(e))
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, ' at Line ', exc_tb.tb_lineno)
            sys.exit(0)

        # extract the dimensionality information from all different MDFs and store it in a list of dictionaries
        self.dimensions = [None] * len(self.num_MDF)
        self.num_chirps = [None] * len(self.num_MDF)

        for iMDF in range(len(self.num_MDF)):
            nc_data_set = netCDF4.Dataset(self.file_MDF[iMDF][0], 'r')
            self.num_chirps[iMDF] = nc_data_set.dimensions['Chirp'].size

            Nav = np.divide(np.array(nc_data_set.variables['AvgNum'][:]),
                            np.array(nc_data_set.variables['DoppLen'][:]))

            self.dimensions[iMDF] = {
                'Program': self.name_MDF[iMDF],
                'NumChirps': self.num_chirps[iMDF],
                'TAlt': nc_data_set.dimensions['TAlt'].size,
                'HAlt': nc_data_set.dimensions['HAlt'].size,
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
            for idim in nc_data_set.dimensions:
                if str(idim).find('C') == 0:
                    self.dimensions[iMDF].update({str(idim) + str(iFile): nc_data_set.dimensions[idim].size})

            # extract first time dimension, if more than one file see loop below
            self.dimensions[iMDF].update({'Time': [nc_data_set.dimensions['Time'].size]})
            nc_data_set.close()

            # start at 1 instead of 0 because the first time dimension was already extrated see lines above
            for iFile in range(1, self.num_MDF[iMDF]):
                nc_data_set = netCDF4.Dataset(self.file_MDF[iMDF][iFile], 'r')
                # range bins and velocity bins are collected a view line above, time bins will be extracted later
                self.dimensions[iMDF]['Time'].append(nc_data_set.dimensions['Time'].size)
                nc_data_set.close()

        # separate different variables by dimension
        # use an arbitrary file because all files contain the same variables (LV0 OR LV1)
        self.time_series_1D = [None] * len(self.num_MDF)
        self.time_series_2D = [None] * len(self.num_MDF)
        self.time_series_3D = [None] * len(self.num_MDF)

        i_nc_file = 0
        n_nc_file = len(self.ncfiles)

        var_list = nc_data_set.variables.keys()

        for iMDF in range(len(self.num_MDF)):
            self.time_series_1D[iMDF] = dict()
            self.time_series_2D[iMDF] = dict()
            self.time_series_3D[iMDF] = dict()

            for iFile in range(self.num_MDF[iMDF]):

                if pts: print("    Loading LIMRAD94 NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")
                i_nc_file += 1

                nc_data_set = netCDF4.Dataset(self.file_MDF[iMDF][iFile], 'r')
                for ivar in var_list:
                    var = nc_data_set.variables[ivar]
                    if 'Units' in var.ncattrs():  # these variables have Units

                        try:
                            if iFile == 0:
                                if len(var.shape) == 1:
                                    self.time_series_1D[iMDF].update({ivar: {
                                        'Dim': list(var.shape), 'LongName': var.Name,
                                        'Unit': var.Units, 'Val': np.array(var[:])}})
                                elif len(var.shape) == 2:
                                    self.time_series_2D[iMDF].update({ivar: {
                                        'Dim': list(var.shape), 'LongName': var.Name,
                                        'Unit': var.Units, 'Val': np.array(var[:, :])}})
                                elif len(var.shape) == 3:
                                    self.time_series_3D[iMDF].update({ivar: {
                                        'Dim': list(var.shape), 'LongName': var.Name,
                                        'Unit': var.Units, 'Val': np.array(var[:, :, :])}})
                            else:
                                if len(var.shape) == 1:
                                    if ivar in ['TAlts', 'HAtls']:
                                        self.time_series_1D[iMDF][ivar]['Val'] = np.append(
                                            self.time_series_1D[iMDF][ivar]['Val'], var[:], axis=0)
                                    else:
                                        self.time_series_1D[iMDF][ivar]['Dim'][0] = \
                                            self.time_series_1D[iMDF][ivar]['Dim'][0] + var.shape[0]

                                        self.time_series_1D[iMDF][ivar]['Val'] = np.append(
                                            self.time_series_1D[iMDF][ivar]['Val'], var[:], axis=0)

                                elif len(var.shape) == 2:
                                    self.time_series_2D[iMDF][ivar]['Dim'][0] = \
                                        self.time_series_2D[iMDF][ivar]['Dim'][0] + var.shape[0]
                                    self.time_series_2D[iMDF][ivar]['Val'] = np.append(
                                        self.time_series_2D[iMDF][ivar]['Val'], var[:, :], axis=0)

                                elif len(var.shape) == 3:
                                    self.time_series_3D[iMDF][ivar]['Dim'][0] = \
                                        self.time_series_3D[iMDF][ivar]['Dim'][0] + var.shape[0]
                                    self.time_series_3D[iMDF][ivar]['Val'] = np.append(
                                        self.time_series_3D[iMDF][ivar]['Val'], var[:, :, :], axis=0)
                        except Exception as e:
                            print('Something went wrong during data type construction: ', e)
                            print('Variable :: ', ivar, '\n')

                            exc_type, exc_obj, exc_tb = sys.exc_info()
                            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                            print(exc_type, fname, ' at Line ', exc_tb.tb_lineno)

                nc_data_set.close()

            # extract the range (height) values
            cnt = 0
            height = [chirpTable_min_height]
            for iC in range(self.num_chirps[iMDF]):
                for _ in range(self.dimensions[iMDF]['Range'][iC]):
                    height.append(height[cnt] + self.dimensions[iMDF]['RangeRes'][iC])
                    cnt += 1

            height.pop(0)
            height = np.array(height) / 1000.0  # convert to km

            self.time_series_1D[iMDF].update({'Height': {'Dim': len(height), 'LongName': 'Range',
                                                         'Unit': 'km', 'Val': np.array(height)}})

            # concatenate different chirps in level 1 data
            if self.lvl == 'LV1':

                for iC in range(1, self.num_chirps[iMDF] + 1):

                    for ivar in var_list:

                        # find different chirp variables using regex
                        regex = re.compile('C' + str(iC))
                        match = re.match(regex, ivar)

                        if match is not None:
                            try:
                                if iC == 1:
                                    self.time_series_2D[iMDF].update({ivar[2:]: {
                                        'Dim': self.time_series_2D[iMDF][ivar]['Dim'],
                                        'LongName': self.time_series_2D[iMDF][ivar]['LongName'][:-9],
                                        'Unit': self.time_series_2D[iMDF][ivar]['Unit'],
                                        'Val': self.time_series_2D[iMDF][ivar]['Val']}})

                                else:
                                    self.time_series_2D[iMDF][ivar[2:]]['Dim'][1] = \
                                        self.time_series_2D[iMDF][ivar[2:]]['Dim'][1] + \
                                        self.time_series_2D[iMDF][ivar]['Dim'][1]
                                    self.time_series_2D[iMDF][ivar[2:]]['Val'] = np.concatenate((
                                        self.time_series_2D[iMDF][ivar[2:]]['Val'],
                                        self.time_series_2D[iMDF][ivar]['Val']), axis=1)

                                # delete data for individual chirps
                                del self.time_series_2D[iMDF][ivar]

                            except Exception as e:
                                print('Something went wrong during data type construction: ', e)
                                print('Variable :: ', ivar, '\n')

                                exc_type, exc_obj, exc_tb = sys.exc_info()
                                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                                print(exc_type, fname, ' at Line ', exc_tb.tb_lineno)

        if pts: print("    Loading LIMRAD94 NC-files ({} of {})".format(n_nc_file, n_nc_file))
