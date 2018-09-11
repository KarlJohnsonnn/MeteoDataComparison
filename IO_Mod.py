
import numpy as np

import netCDF4
import datetime
from datetime import timezone
import glob
import os

import sys

import time

from Parameter_Mod import *

def get_nc_data(thisfile, varname):
    # if pts: print('loading variable '+varname +' from ' + thisfile)
    ncfile = netCDF4.Dataset(thisfile,'r')
    var    = ncfile.variables[varname]

    if ncfile.isopen==1: ncfile.close()
    return var

def get_nc_date(thisfile):
    # if pts: print('loading variable '+varname +' from ' + thisfile)
    ncfile = netCDF4.Dataset(thisfile,'r')
    year  = ncfile.year
    month = ncfile.month
    day   = ncfile.day

    if ncfile.isopen==1: ncfile.close()
    return year,month,day

def get_nc_dimension(thisfile,dim_name):
    ncfile = netCDF4.Dataset(thisfile,'r')
    dim    = ncfile.dimensions[dim_name].size
    if ncfile.isopen==1: ncfile.close()
    return dim

def get_height_boundary(Array,hmin,hmax):
    i = 0
    imin = 0
    for h in Array:
        #print(' (min) h = ',h,hmin)
        if (hmin <= h):
            imin = i
            break
        i += 1

    i = 0
    imax = 0
    for h in reversed(Array):
        #print(' (max) h = ',h,hmax)
        if (hmax >= h):
            imax = len(Array)-i
            break
        i += 1

    return imin, imax


def extract_dataset(date,time,clock,h_bounds,fext,kind):
    '''
    Extract data from NetCDF files
    :param date: string of desired date 'YYMMDD'
    :param time: list of datetimes
    :param clock: deciaml hours [begin, end]
    :param fext:  file extention of dataset
    :return:
    '''
    if fext == '*mira.nc':
        os.chdir(meteo_path+'MIRA/calibrated')  # path to data needs to be fit to the devices file structure
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
            if (time[0] <= zeit <= time[1] ): min_time = i
            if (time[2] <= zeit <= time[3] ): max_time = i
            i += 1

        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        height    = height[imin_h:imax_h]

        # get data from .nc files
        Ze  = np.array(get_nc_data(file, 'Zh'))
        mdv = np.array(get_nc_data(file, 'v'))
        sw  = np.array(get_nc_data(file, 'width'))

        # np.warnings.filterwarnings('ignore')

        Ze  = Ze [min_time:max_time, imin_h:imax_h]
        mdv = mdv[min_time:max_time, imin_h:imax_h]
        sw  = sw [min_time:max_time, imin_h:imax_h]

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

        os.chdir(meteo_path+'MIRA/mmclx')  # path to data needs to be fit to the devices file structure

        #ncfiles = glob.glob('20' + date + '*.mmclx')     # 20180727_000013.mmclx

        first_file = int(clock[0]) - np.remainder(int(clock[0]), 3)
        if clock[1] - int(clock[1]) > 0.0 :
            last_file  = int(clock[1]) + 1
        else:
            last_file  = int(clock[1])

        range_file_list = list(range(first_file, last_file, 3))


        ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob( '20' + date + '_' + str(il).zfill(2) + '*.mmclx' ))
            ncfiles.append(file_name[2:-2])

        if file_name[2:-2] == '':
            print('   Error!  File not found --> exit!')
            print('   Check LIMRAD folder!')
            exit(0)

        file = ncfiles[0]

        if file == '':
            print('   Error!  File: "'+file+'" not found --> exit!')
            print('   Check MIRA folder!')
            exit(0)

        # extract date (year, month, day)
        #mira_y, mira_m, mira_d = get_nc_date(ncfiles[0])

        # extract range array
        height = np.array(get_nc_data(file, 'range'))
        imin_h, imax_h = get_height_boundary(height,h_bounds[0]*1000,h_bounds[1]*1000)
        #20180728_060014.mmclx

        if   kind == 'cl':
            kind1 = ''
        elif kind == 'g':
            kind1 ='g'
        else:
            kind1 = 'e'


        # conversion from deciaml hour to datetime
        time_samp = np.array(get_nc_data(file, 'time'))
        Ze  = np.array(get_nc_data(file, 'Z'+kind1))
        mdv = np.array(get_nc_data(file, 'VEL'+kind))
        sw  = np.array(get_nc_data(file, 'RMS'+kind))


        i_nc_file = 0
        n_nc_file = len(ncfiles)

        for file in ncfiles[1:]:
            i_nc_file += 1
            if pts: print("    Loading MIRA35 (mmclx) NC-files ({} of {})".format(i_nc_file,n_nc_file), end="\r")

            time_samp = np.append(time_samp, get_nc_data(file,'time'))
            Ze  = np.append(Ze,  get_nc_data(file, 'Z'+kind1),  axis=0)
            mdv = np.append(mdv, get_nc_data(file, 'VEL'+kind), axis=0)
            sw  = np.append(sw,  get_nc_data(file, 'RMS'+kind), axis=0)

        if pts: print("    Loading MIRA35 (mmclx) NC-files ({} of {})".format(n_nc_file,n_nc_file))

        time_plot = [ datetime.datetime(1970, 1, 1, 0, 0, 0)
                    + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp)) ]


        i = 0
        for zeit in time_plot:
            if ( time[0]  <= zeit <= time[1] ): min_time = i
            if ( time[2]  <= zeit <= time[3] ): max_time = i
            i += 1


        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        height    = np.divide(height[imin_h:imax_h],1000.0)

        Ze  = Ze [min_time:max_time, imin_h:imax_h]
        mdv = mdv[min_time:max_time, imin_h:imax_h]
        sw  = sw [min_time:max_time, imin_h:imax_h]

        Ze  = np.ma.masked_invalid(Ze).T
        Ze  = np.ma.log10(Ze) * 10
        mdv = np.ma.masked_invalid(mdv).T
        sw  = np.ma.masked_invalid(sw).T


        os.chdir('../../')  # path to data needs to be fit to the devices file structure


    elif fext == '*.LV1.NC':

        os.chdir(meteo_path+'LIMRad94/')  # path to data needs to be fit to the devices file structure

        first_file = int(clock[0])
        if clock[1]-int(clock[1]) > 0.0:
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
        hgt_gates = []
        dummy = None

        file = ncfiles[0]

        no_c = get_nc_dimension(file, 'Chirp')

        # find the number of range gates per chirp sequence
        hgt_res = get_nc_data(file, 'RangeRes')
        for i in range(0, no_c):
            dummy = get_nc_data(file, 'C' + str(i + 1) + 'MeanVel')
            hgt_gates = np.append(hgt_gates, len(dummy[0, :]))


        # calculate LR_height levels
        for i in range(0, int(hgt_gates[0])):
            height = np.append(height, (i + 1) * hgt_res[0])
        for i in range(0, int(hgt_gates[1])):
            i_temp = int(hgt_gates[0] - 1 + i)
            height = np.append(height, height[i_temp] + hgt_res[1])
        for i in range(0, int(hgt_gates[2])):
            i_temp = int(hgt_gates[0] + hgt_gates[1] - 1 + i)
            height = np.append(height, height[i_temp] + hgt_res[2])


        # get data from .nc files
        time_samp = get_nc_data(file, 'Time')
        Ze1  = get_nc_data(file, 'C1ZE')
        Ze2  = get_nc_data(file, 'C2ZE')
        Ze3  = get_nc_data(file, 'C3ZE')
        mdv1 = get_nc_data(file, 'C1MeanVel')
        mdv2 = get_nc_data(file, 'C2MeanVel')
        mdv3 = get_nc_data(file, 'C3MeanVel')
        sw1  = get_nc_data(file, 'C1SpecWidth')
        sw2  = get_nc_data(file, 'C2SpecWidth')
        sw3  = get_nc_data(file, 'C3SpecWidth')

        i_nc_file = 0
        n_nc_file = len(ncfiles)

        for file in ncfiles[1:]:
            i_nc_file += 1
            if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

            time_samp = np.append(time_samp, get_nc_data(file, 'Time'))
            Ze1  = np.append(Ze1, get_nc_data(file, 'C1ZE'), axis=0)
            Ze2  = np.append(Ze2, get_nc_data(file, 'C2ZE'), axis=0)
            Ze3  = np.append(Ze3, get_nc_data(file, 'C3ZE'), axis=0)
            mdv1 = np.append(mdv1, get_nc_data(file, 'C1MeanVel'), axis=0)
            mdv2 = np.append(mdv2, get_nc_data(file, 'C2MeanVel'), axis=0)
            mdv3 = np.append(mdv3, get_nc_data(file, 'C3MeanVel'), axis=0)
            sw1  = np.append(sw1, get_nc_data(file, 'C1SpecWidth'), axis=0)
            sw2  = np.append(sw2, get_nc_data(file, 'C2SpecWidth'), axis=0)
            sw3  = np.append(sw3, get_nc_data(file, 'C3SpecWidth'), axis=0)

        if pts: print("    Loading LIMRAD94 (LV1.NC) NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format


        time_plot = [ datetime.datetime(2001, 1, 1, 0, 0, 0)
                      + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp)) ]

        i = 0
        for zeit in time_plot:
            if (time[0] <= zeit <= time[1]): min_time = i
            if (time[2] <= zeit <= time[3]): max_time = i
            i += 1

        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        # for zeit in UTC_time_LR:
        #    print('zeit = ', zeit)

        height = height + chirpTable_min_height*1000.  # get LR_height in km

        imin_h, imax_h = get_height_boundary(height, 1000*h_bounds[0], 1000*h_bounds[1])
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



def save_log_data(filename,meth,hmin,hmax,comp_date,comp_time_int):
    file = open(filename + '.log', 'w')

    file.write('')
    file.write(' This is the log file for LIMRad 94GHz - MIRA 35GHz Radar Data'+'\n'*2)
    file.write(' # User Inputs\n')
    file.write('    minimum height = ' + str(hmin) + '\n')
    file.write('    maximum height = ' + str(hmax) + '\n'*2)

    file.write('    date = ' + str(comp_date) + ' # in "YYMMDD" \n')
    file.write('    time intervall = ' + str(comp_time_int) + ' # in "HHMM-HHMM"'+'\n'*2)


    file.write(' # Logicals \n')
    file.write('    pts = ' + str(pts) + ' # print to screen\n')
    file.write('    dbg = ' + str(dbg) + ' # debugging flag, show some parameter\n')
    file.write('    plot_RectBivariateSpline   = ' + str(plot_RectBivariateSpline) + ' # bivariate interpolation plot\n')
    file.write('    plot_interpolation_scatter = ' + str(plot_interpolation_scatter) + ' # phase diagram of eqv. radar reflectivity\n')
    file.write('    plot_compare_mira_mmclx    = ' + str(plot_compare_mira_mmclx) + ' # plot Ze, mdv, sw of LIMRad and MIRA (all targes, only hydrometeors)\n')
    file.write('    plot_radar_results = ' + str(plot_radar_results) + ' # plot radar data of LIMRad and MIRA (Ze, mdv, sw)\n')
    file.write('    plot_comparisons   = ' + str(plot_comparisons) + ' # plot time/height-averaged data\n')
    file.write('    plot_interp2d      = ' + str(plot_interp2d) + ' # plot 2D interpolation; LIMRad data onto MIRA grid and vice versa'+'\n'*2)

    if plot_interpolation_scatter:
        file.write(' # phase-diagram interpolation parameter\n')
        file.write('    interpolation method = ' + meth + '\n')
        #file.write('    resolution of interpolated points = ' + str(res_interp) + '\n')

    file.close()