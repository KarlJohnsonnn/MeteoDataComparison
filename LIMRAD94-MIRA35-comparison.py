
import numpy as np

import netCDF4

import matplotlib.pyplot as plt
import matplotlib        as mpl

from matplotlib import rc
from matplotlib import colors

rc('font',family='serif')
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from matplotlib              import dates
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.colors as colors

import datetime
from datetime import timezone
import glob
import os

import sys
import warnings

import time

import pandas as pd
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline


# colorbar sw


##################################################################################################
#
#       ##     ##  ######  ######## ########     #### ##    ## ########  ##     ## ########
#       ##     ## ##    ## ##       ##     ##     ##  ###   ## ##     ## ##     ##    ##
#       ##     ## ##       ##       ##     ##     ##  ####  ## ##     ## ##     ##    ##
#       ##     ##  ######  ######   ########      ##  ## ## ## ########  ##     ##    ##
#       ##     ##       ## ##       ##   ##       ##  ##  #### ##        ##     ##    ##
#       ##     ## ##    ## ##       ##    ##      ##  ##   ### ##        ##     ##    ##
#        #######   ######  ######## ##     ##    #### ##    ## ##         #######     ##
#
##################################################################################################

# Logicals for different tasks
plot_interp2d              = True
plot_RectBivariateSpline   = False           #interp2 testen
plot_radar_results         = False
plot_comparisons           = False
plot_interpolation_scatter = False
plot_compare_mira_mmclx    = False

pts = True # print to screen
dbg = False

# file type
LIMRad_file_extension = '*.LV1.NC'
#mira_file_extension = '*mira.nc'
mira_file_extension = '*.mmclx'

#constants
chirpTable_min_height = 0.1 # (km)

#dots per inch parameter for .png
dpi_val = 200




#################################################################################################
#
#       ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
#       ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
#       ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
#       ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
#       ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
#       ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
#       ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######
#
#################################################################################################

# Subroutine definition and start of the program
def dim(a):
    if not type(a) == list: return []
    return [len(a)] + dim(a[0])

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

def save_log_data(filename,meth,res_interp):
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
    file.write('    plot_RectBivariateSpline = ' + str(plot_RectBivariateSpline) + ' # bivariate interpolation plot\n')
    file.write('    plot_interpolation_scatter    = ' + str(plot_interpolation_scatter) + ' # phase diagram of eqv. radar reflectivity\n')
    file.write('    plot_radar_results = ' + str(plot_radar_results) + ' # plot radar data of LIMRad and MIRA (Ze, mdv, sw)\n')
    file.write('    plot_comparisons   = ' + str(plot_comparisons) + ' # plot time/height-averaged data'+'\n'*2)

    if plot_interpolation_scatter:
        file.write(' # phase-diagram interpolation parameter\n')
        file.write('    interpolation method = ' + meth + '\n')
        file.write('    resolution of interpolated points = ' + str(res_interp) + '\n')

    file.close()

def lookupNearest(x0, y0):
   xi = numpy.abs(x-x0).argmin()
   yi = numpy.abs(y-y0).argmin()
   return data[yi,xi]

def extract_dataset(date,time,clock,fext,kind):
    '''
    Extract data from NetCDF files
    :param date: string of desired date 'YYMMDD'
    :param time: list of datetimes
    :param clock: deciaml hours [begin, end]
    :param fext:  file extention of dataset
    :return:
    '''
    if fext == '*mira.nc':
        os.chdir('MIRA')  # path to data needs to be fit to the devices file structure
        ncfiles = glob.glob('20' + date + '*mira.nc')

        if pts: print("    Loading MIRA35 NC-files ({} of {})".format(0, 1), end="\r")

        file = ncfiles[0]

        if file == '':
            print('   Error!  File: "' + file + '" not found --> exit!')
            print('   Check MIRA folder!')
            exit(0)

        height = np.array(get_nc_data(file, 'range'))
        imin_h, imax_h = get_height_boundary(height, hmin, hmax)

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

        os.chdir('../')  # path to data needs to be fit to the devices file structure
        if pts: print("    Loading MIRA35 NC-files ({} of {})".format(1, 1))

    elif fext == '*.mmclx':

        os.chdir('MIRA')  # path to data needs to be fit to the devices file structure

        #ncfiles = glob.glob('20' + date + '*.mmclx')     # 20180727_000013.mmclx

        first_file = int(clock[0]) - np.remainder(int(clock[0]), 3)
        if comp_minutes[1] > 0.0 :
            last_file  = int(clock[1]) + 1
        else:
            last_file  = int(clock[1])

        range_file_list = list(range(first_file, last_file, 3))


        ncfiles = []
        for il in range_file_list:
            file_name = str(glob.glob( '20' + comp_date + '_' + str(il).zfill(2) + mira_file_extension ))
            ncfiles.append(file_name[2:-2])

        if file_name[2:-2] == '':
            print('   Error!  File: "'+file+'" not found --> exit!')
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
        imin_h, imax_h = get_height_boundary(height,hmin*1000,hmax*1000)
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
            if pts: print("    Loading MIRA35 NC-files ({} of {})".format(i_nc_file,n_nc_file), end="\r")

            time_samp = np.append(time_samp, get_nc_data(file,'time'))
            Ze  = np.append(Ze,  get_nc_data(file, 'Z'+kind1),  axis=0)
            mdv = np.append(mdv, get_nc_data(file, 'VEL'+kind), axis=0)
            sw  = np.append(sw,  get_nc_data(file, 'RMS'+kind), axis=0)

        if pts: print("    Loading MIRA35 NC-files ({} of {})".format(n_nc_file,n_nc_file))

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


        os.chdir('../')  # path to data needs to be fit to the devices file structure


    elif fext == '*.LV1.NC':

        os.chdir('LIMRAD')  # path to data needs to be fit to the devices file structure

        first_file = int(clock[0])
        if comp_minutes[1] > 0.0:
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


        ##flat_list = [item for sublist in ncfiles for item in sublist]
        #flat_list = []
        #for sublist in ncfiles:
        #    print('asdfa = ',sublist)
        #    for item in sublist:
        #        flat_list.append(item)
        #        print('über = ',item)
#
 #       ncfiles = flat_list
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
            if pts: print("    Loading LIMRAD94 NC-files ({} of {})".format(i_nc_file, n_nc_file), end="\r")

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

        if pts: print("    Loading LIMRAD94 NC-files ({} of {})".format(n_nc_file, n_nc_file))

        # convert times in datetime format


        time_plot = [ datetime.datetime(2001, 1, 1, 0, 0, 0)
                      + datetime.timedelta(seconds=int(time_samp[i])) for i in range(len(time_samp)) ]

        i = 0
        for zeit in time_plot:
            if (time_int[0] <= zeit <= time_int[1]): min_time = i
            if (time_int[2] <= zeit <= time_int[3]): max_time = i
            i += 1

        time_samp = time_samp[min_time:max_time]
        time_plot = time_plot[min_time:max_time]
        # for zeit in UTC_time_LR:
        #    print('zeit = ', zeit)

        height = height + chirpTable_min_height*1000.  # get LR_height in km

        imin_h, imax_h = get_height_boundary(height, 1000*hmin, 1000*hmax)
        height = np.divide(height[imin_h:imax_h], 1000)

        # stack variables of individual chirps
        Ze = np.hstack((Ze1, Ze2, Ze3))
        Ze = Ze[min_time:max_time, imin_h:imax_h]
        Ze = np.transpose(Ze)
        Ze = np.ma.masked_less_equal(Ze, -999.)
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

    return time_samp, time_plot, height, Ze, mdv, sw

def place_text(plot,pos,text):
    plot.text( pos[0], pos[1], text,
               fontweight = 'bold',
               horizontalalignment = 'left',
               transform = plot.transAxes,
               bbox = dict( facecolor = 'white',
                            edgecolor ='black',
                            pad = 5.)
               )

def place_statistics(plot,pos,stat,vn):
    if vn=='Ze':
        text  = r'$\overline{ Z_e^{94} - Z_e^{35} }=$' + \
                 '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'$\rho(Z_e^{94}, Z_e^{35}) = $' + '{:6.2f}'.format(stat[1])
    if vn=='mdv':
        text = r'$\mathrm{mean}_h(\mathrm{mdv}_{\mathrm{lr}} - \mathrm{mdv}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])
    if vn=='sw':
        text = r'$\mathrm{mean}_h(\mathrm{sw}_{\mathrm{lr}} - \mathrm{sw}_{\mathrm{mi}} ) =$' \
               +'{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])

    plot.text( pos[0], pos[1], text, fontweight = 'bold',
               horizontalalignment = 'center',
               transform = plot.transAxes,
               bbox=dict(facecolor='none',
                         edgecolor='black',
                         pad=5.))
    plot.text(pos[0]+0.6, pos[1], text2, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))

def place_statistics_mean_corcoef(plot,pos,stat,vn):
    if vn=='Ze':
        text  = r'$\overline{Z_e^{94} - Z_e^{35} } =$' + \
                 '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'corr$(Z_e^{94}, Z_e^{35}) = $'+ r'$=$' + '{:6.2f}'.format(stat[1])
    if vn=='mdv':
        text = r'$\mathrm{mean}_h(\mathrm{mdv}_{\mathrm{lr}} - \mathrm{mdv}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])
    if vn=='sw':
        text = r'$\mathrm{mean}_h(\mathrm{sw}_{\mathrm{lr}} - \mathrm{sw}_{\mathrm{mi}} ) =$' \
               +'{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])

    plot.text( pos[0], pos[1], text, fontweight = 'bold',
               horizontalalignment = 'center',
               transform = plot.transAxes,
               bbox=dict(facecolor='none',
                         edgecolor='black',
                         pad=5.))
    plot.text(pos[0]+0.6, pos[1], text2, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))

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

def np_NaN(n,m):
    mat = np.zeros((n,m))
    mat[:,:] = np.nan
    return mat

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
    len_time  = shape_dbz[1]
    len_range = len(range_v)

    bases = np_NaN(10,len_time)
    tops = np_NaN(10,len_time)
    thickness = np_NaN(10,len_time)

    top_m  = np_NaN(10,len_time) # max. 10 cloud layers detected
    base_m = np_NaN(10,len_time) # max. 10 cloud layers detected

    if pts:
        print('')
        print(' dBZ(i,:) = ', [dbz_m[i, :] for i in range(shape_dbz[0])])

    for i in range(0,len_time):

        in_cloud  = 0
        layer_idx = 0
        current_base = np.nan

        if pts: print("    Searching for cloud bottom and top ({} of {}) time steps".format(i+1,len_time), end="\r")

        # found the first base in first bin.
        if ( not np.isnan(dbz_m[0,i]) ):
            layer_idx = 1
            current_base = 1
            in_cloud = 1

        for j in range(1,len_range):

            if (in_cloud == 1): # if in cloud

                # cloud top found at (j-1)
                if np.isnan(dbz_m[j,i]):

                    current_top = j-1
                    thickness[layer_idx,i] = range_v[current_top] - range_v[current_base]
                    bases[layer_idx,i]     = current_base  # bases is an idx
                    tops[layer_idx,i]      = current_top    # tops is an idx

                    base_m[layer_idx,i]    = range_v[current_base]  # cloud base in m or km
                    top_m[layer_idx,i]     = range_v[current_top]   # cloud top in m or km

                    print(str(i)+': found '+str(layer_idx)+'. cloud ['+str(bases[layer_idx,i])+', '+\
                        str(tops[layer_idx,i])+'], thickness: '+str(thickness)+'km')

                    in_cloud = 0

            else: # if not in cloud

                # cloud_base found at j
                if ( not np.isnan(dbz_m[j,i]) ):
                    layer_idx += 1
                    current_base = j
                    in_cloud = 1

        # at top height but still in cloud, force top
        if ( in_cloud == 1 ):
            tops[layer_idx,i]  = len(range_v)
            top_m[layer_idx,i] = max(range_v)    #  cloud top in m or km



    ###
    # keep only first 10 cloud layers
    bases     = bases[:10,:]
    tops      = tops [:10,:]
    base_m    = base_m[:10,:]
    top_m     = top_m [:10,:]
    thickness = thickness[:10,:]
    # give clear sky flag when first base_m ==NaN (no dbz detected over all heights),
    # problem: what if radar wasn't working, then dbz would still be NaN!
    loc_nan = np.where(np.isnan(base_m[0,:]))

    if pts:
        print('vor bases = ',[ba for ba in base_m])

    base_m[0,np.where(np.isnan(base_m[0,:]))] = -1
    top_m[0,np.where(np.isnan(top_m[0,:]))]   = -1

    if pts:
        print('')
        print('npisnan = ',[loc_nan[i] for i in range(len(loc_nan))])


        if pts: print('')
        for ba in base_m:
            print(' bases = ',ba)

        for to in top_m:
            print(' tops = ', to)


    return bases, tops, base_m, top_m, thickness


def plot_data_set(axh,text,x,y,z,vmi,vma,x_min,x_max,y_min,y_max,x_lab,y_lab,z_lab):
    text = r'\textbf{'+text+'}'
    if text.find('sw')>0:
        cp = axh.pcolormesh(x, y, z,  norm=colors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
    else:
        place_text(axh, [.02, 1.05], text )
        cp = axh.pcolormesh(x, y, z, vmin=vmi, vmax=vma, cmap='jet')
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax0 = divider1.append_axes("right", size="3%", pad=0.1)
    cbar= fig.colorbar(cp, cax=cax0, ax=axh)
    cbar.set_label(z_lab)
    axh.axes.tick_params(axis='both', direction='inout', length=10, width=1.5)
    axh.set_ylabel(y_lab)
    axh.set_xlim( left = x_min, right = x_max )
    axh.set_ylim( bottom = y_min, top = y_max, )

    # exceptions
    if x_lab == '':
        axh.axes.xaxis.set_ticklabels([])
    else:
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
        axh.set_xlabel(x_lab)

def plot_correlation_matrix(axh,text,z,vmi,vma,x_lab,y_lab,z_lab):
    from pylab import pcolor
    #text = r'\textbf{'+text+'}'
    #place_text(axh, [.02, 1.05], text )
    #z = z[1500:2000,1500:2000]
    cp = axh.pcolor(z, norm=colors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax0 = divider1.append_axes("right", size="3%", pad=0.1)
    cbar= fig.colorbar(cp, cax=cax0, ax=axh)
    cbar.set_label(z_lab)
    axh.axes.tick_params(axis='both', direction='inout', length=10, width=1.5)

    # exceptions
    #axh.axes.xaxis.set_ticklabels([]) if x_lab == '' else axh.set_xlabel(x_lab)
    #axh.axes.yaxis.set_ticklabels([]) if y_lab == '' else axh.set_ylabel(y_lab)


    #axh.set_aspect('equal', 'box')

def plot_correlation(axh,text,x,y,label,marker,x_min,x_max,y_min,y_max,x_lab,y_lab):
    place_text(axh, [.15, 1.1], text )
    axh.plot(x,y,marker,label=label)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    divider1 = make_axes_locatable(axh)
    cax = divider1.append_axes("right", size="3%", pad=0.1)
    axh.legend(loc="upper right")

    cax.set_facecolor('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        cax.spines[axis].set_linewidth(0)
    cax.set_xticks([])
    cax.set_yticks([])

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))

def plot_avg_data_set(axh,text,x1,y1,x2,y2,label1,marker1,label2,marker2,x_min,x_max,y_min,y_max,x_lab,y_lab,ax):
    text = r'\textbf{' + text + '}'
    place_text(axh, [.02, 1.05], text )
    axh.scatter(x1,y1,marker=marker1,label=label1)
    axh.scatter(x2,y2,marker=marker2,label=label2)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    if ax == 'y':
        divider1 = make_axes_locatable(axh)
        cax = divider1.append_axes("right", size="3%", pad=0.1)
        axh.legend(loc="upper right")

        cax.set_facecolor('none')
        for axis in ['top','bottom','left','right']:
            cax.spines[axis].set_linewidth(0)
        cax.set_xticks([])
        cax.set_yticks([])

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))

def plot_phase_data_set(axh,text,x,y,marker,x_min,x_max,y_min,y_max,x_lab,y_lab):
    text = r'\textbf{' + text + '}'
    place_text(axh, [.02, 1.05], text )
    axh.scatter(x,y,marker=marker)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax = divider1.append_axes("right", size="3%", pad=0.25)
    axh.legend(loc="upper right")

    cax.set_facecolor('none')
    for axis in ['top','bottom','left','right']:
        cax.spines[axis].set_linewidth(0)
    cax.set_xticks([])
    cax.set_yticks([])

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))

def plot_interpol_data_set(axh,text,x1,y1,x2,y2,label1,marker1,label2,marker2,x_min,x_max,y_min,y_max,x_lab,y_lab):
    place_text(axh, [.15, 1.1], text )
    axh.plot(x1,y1,marker1,label=label1)
    axh.plot(x2,y2,marker2,label=label2)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))



def plot_scatter(axh,text,x,y,marker,x_min,x_max,y_min,y_max,x_lab,y_lab):
    import matplotlib.ticker as ticker
    from matplotlib.patches     import Polygon
    from matplotlib.collections import PatchCollection

    #place_text(axh, [.05, 1.1], text )
    axh.plot(x,y,marker)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min,   right = x_max )
    axh.set_ylim( bottom = y_min, top = y_max)
    axh.xaxis.set_ticks((np.arange(x_min, x_max, (x_max-x_min)/4.0)))
    axh.yaxis.set_ticks((np.arange(y_min, y_max, (y_max-y_min)/4.0)))
    axh.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    axh.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")
    axh.set_aspect('equal', 'box')

    # plot 1:1 line
    N = 100
    X = np.linspace(x_min, x_max, N)
    axh.plot(X, X, 'k--')

    # add patches
    colors = [np.divide([31,119,180],255.) , np.divide([255,127,14],255.)]
    patches = [ Polygon( [[x_min,y_min],[x_max,y_min],[x_max,y_max]], facecolor='C0', fill=True),
                Polygon( [[x_min,y_min],[x_min,y_max],[x_max,y_max]], facecolor='C1', fill=True)]

    p = PatchCollection(patches, alpha=0.25)
    p.set_color(colors)
    axh.add_collection(p)

def get_plot_ybounds(y1,y2,pm):
    import math
    return max(y1.min(),y2.min())-pm, min(y1.max(),y2.max())+pm


def gather_time_lists(dtime,dt_max):

    dbg1 = False
    i_bin = 0
    index = 0
    last_time = dtime[0]
    sutime_samp = 0.0
    idx_list = [[]]
    times_list = [[]]
    if dbg1: print(' type = ', dt_max, dt_max)

    if isinstance(dt_max,list):
        if dbg1: print('')
        for curr_time in dtime[1:]:

            if dbg1: print(' curr_time = ', curr_time, sutime_samp)
            if ( sutime_samp < dt_max[i_bin] ):
                seconds = (curr_time-last_time).total_seconds()
                sutime_samp += seconds
                idx_list[i_bin].append(index)
                times_list[i_bin].append(seconds)
                if dbg1:  print(' time = ',curr_time,'     -     ',curr_time-last_time,' time to float = ', seconds)

                if ( sutime_samp >= dt_max[i_bin] ):
                    i_bin += 1
                    times_list.append([])
                    idx_list.append([])
                    idx_list[i_bin].append(index)
                    times_list[i_bin].append(seconds)
                    sutime_samp = 0.0

            last_time = curr_time
            index += 1

    else:
        if dbg1: print('')
        for curr_time in dtime[1:]:

            if dbg1: print(' curr_time = ', curr_time, sutime_samp)
            if ( sutime_samp < dt_max ):
                seconds = (curr_time-last_time).total_seconds()
                sutime_samp += seconds
                idx_list[i_bin].append(index)
                times_list[i_bin].append(seconds)
                if dbg1: print(' time = ',curr_time,'     -     ',curr_time-last_time,' time to float = ', seconds)

                if ( sutime_samp >= dt_max ):
                    i_bin += 1
                    times_list.append([])
                    idx_list.append([])
                    idx_list[i_bin].append(index)
                    times_list[i_bin].append(seconds)
                    sutime_samp = 0.0

            last_time = curr_time
            index += 1


    return idx_list, times_list


def interpolate_data(x,y,xnew,method):

    # create a callable function from the actual data
    fnew = interpolate.interp1d(x, y, kind=method)

    # calculate the interpolation dataset
    ynew = fnew(xnew)

    # mask values ( minInterpol=minAcutalData, etc. )
    ynew = np.ma.masked_greater_equal(ynew, y.max())
    ynew = np.ma.masked_less_equal(ynew,    y.min())

    return ynew

def Interpolate_NearestND_2D(x1,y1,z1,x2,y2):
    len_x1 = len(x1);   len_x2 = len(x2)
    len_y1 = len(y1);   len_y2 = len(y2)
    coord  = np.empty(( len_x1*len_y1, 2))
    values = np.empty(( len_x1*len_y1, 1))
    cnt = 0
    for i in range(len_x1):
        for j in range(len_y1):
            coord[cnt,0] = x1[i]
            coord[cnt,1] = y1[j]
            values[cnt] = z1[j,i]
            cnt += 1

    values = np.ma.masked_invalid(values)
    fcn = interpolate.NearestNDInterpolator(coord, values)

    Interp_z = np.empty((len_x2*len_y2))
    cnt = 0
    for xi in x2:
        for yi in y2:
            Interp_z[cnt] = fcn(xi, yi)
            cnt += 1

    Interp_z = np.ma.masked_equal(Interp_z, 0.0)
    Interp_z = np.ma.masked_less_equal(Interp_z, -100.0)
    Interp_z = np.ma.masked_invalid(Interp_z)
    Interp_z = np.reshape(Interp_z, (len_x2, len_y2)).T

    return Interp_z


####################################################################################################################
#
#   ##     ##    ###    #### ##    ##         ########  ########   #######   ######   ########     ###    ##     ##
#   ###   ###   ## ##    ##  ###   ##         ##     ## ##     ## ##     ## ##    ##  ##     ##   ## ##   ###   ###
#   #### ####  ##   ##   ##  ####  ##         ##     ## ##     ## ##     ## ##        ##     ##  ##   ##  #### ####
#   ## ### ## ##     ##  ##  ## ## ## ####### ########  ########  ##     ## ##   #### ########  ##     ## ## ### ##
#   ##     ## #########  ##  ##  ####         ##        ##   ##   ##     ## ##    ##  ##   ##   ######### ##     ##
#   ##     ## ##     ##  ##  ##   ###         ##        ##    ##  ##     ## ##    ##  ##    ##  ##     ## ##     ##
#   ##     ## ##     ## #### ##    ##         ##        ##     ##  #######   ######   ##     ## ##     ## ##     ##
#
####################################################################################################################

start_time = time.clock()

# Print Head
print(' ')
print('  \u250F' + 49*'\u2501' + '\u2513')
print('  \u2503' + '          LIMRAD94 - MIRA35  Comparison          ' + '\u2503')
print('  \u2517' + 49*'\u2501' + '\u251B' + '\n')
print('')
print('')


# gather arguments
if len(sys.argv) == 6:

    hmin, hmax = float(sys.argv[1]), float(sys.argv[2])
    comp_date  = str(sys.argv[3])
    comp_time_int = str(sys.argv[4])+'-'+str(sys.argv[5])

else:

    ## cirrus
    hmin = 8.50 #(km)  - lower y-axis limit
    hmax = 10.0 #(km) - upper y-axis limit, highest range gate may be higher
    comp_date     = '180728'     # in YYMMDD
    comp_time_int = '0740-0810'  # in HHMM-HHMM

    ##cummulis
    #hmin = 5.0 #(km)  - lower y-axis limit
    #hmax = 12.0 #(km) - upper y-axis limit, highest range gate may be higher
    #comp_date     = '180802'     # in YYMMDD
    #comp_time_int = '0330-1200'  # in HHMM-HHMM

    ##nimbus
    #hmin = 0.0 #(km)  - lower y-axis limit
    #hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    #comp_date     = '180805'     # in YYMMDD
    #comp_time_int = '0510-0620'  # in HHMM-HHMM

    #hmin = 0.0 #(km)  - lower y-axis limit
    #hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    #comp_date     = '180805'     # in YYMMDD
    #comp_time_int = '1030-1200'  # in HHMM-HHMM

    #hmin = 1.0 #(km)  - lower y-axis limit
    #hmax = 2.5 #(km) - upper y-axis limit, highest range gate may be higher
    #comp_date     = '180805'     # in YYMMDD
    #comp_time_int = '0700-1210'  # in HHMM-HHMM

    ##nimbus
    #hmin = 0.0 #(km)  - lower y-axis limit
    #hmax = 12.00 #(km) - upper y-axis limit, highest range gate may be higher
    #comp_date     = '180808'     # in YYMMDD
    #comp_time_int = '1330-1700'  # in HHMM-HHMM

warnings.filterwarnings("ignore")

# calculate the time in decimal hours
comp_hours    = [int(comp_time_int[0:2]) , int(comp_time_int[5:7])]
comp_minutes  = [int(comp_time_int[2:4]) , int(comp_time_int[7:9])]

clock_time  = np.array(comp_hours) + np.divide(comp_minutes, 60.) # [hours] + [minutes]/60#


# -- gathering year, month, day for convertion to UTC time
plotyear  = int('20'+comp_date[:2])
plotmonth = int(comp_date[2:4])
plotday   = int(comp_date[4:6])


time_int    = [0, 0, 0, 0]
time_int[0] = datetime.datetime(plotyear, plotmonth, plotday,
                               hour=int(comp_hours[0]), minute=int(comp_minutes[0]))
time_int[1] = time_int[0] + datetime.timedelta(seconds=15)
time_int[3] = datetime.datetime(plotyear, plotmonth, plotday,
                               hour=int(comp_hours[1]), minute=int(comp_minutes[1]))
time_int[2] = time_int[3] - datetime.timedelta(seconds=15)




######################################################################################################
#
#   ########     ###    ########    ###                  #### ##    ## ########  ##     ## ########
#   ##     ##   ## ##      ##      ## ##                  ##  ###   ## ##     ## ##     ##    ##
#   ##     ##  ##   ##     ##     ##   ##                 ##  ####  ## ##     ## ##     ##    ##
#   ##     ## ##     ##    ##    ##     ##    #######     ##  ## ## ## ########  ##     ##    ##
#   ##     ## #########    ##    #########                ##  ##  #### ##        ##     ##    ##
#   ##     ## ##     ##    ##    ##     ##                ##  ##   ### ##        ##     ##    ##
#   ########  ##     ##    ##    ##     ##               #### ##    ## ##         #######     ##
#
######################################################################################################


# LIMRad 94GHz Radar data extraction

LR_time, UTC_time_LR, LR_height, \
LR_Ze,   LR_mdv,       LR_sw      = extract_dataset(comp_date, time_int, clock_time, LIMRad_file_extension, '')


# MIRA 35GHz Radar data extraction

mira_time, UTC_time_mira, mira_height, \
mira_Ze,   mira_mdv,       mira_sw      = extract_dataset(comp_date, time_int, clock_time, mira_file_extension , 'g')

print('')

if dbg:
    print('')
    print('   dim(LR_height,UTC_time_LR) = ',LR_height.shape,len(UTC_time_LR))
    print('          LR_Ze (dim,min,max)  = ',LR_Ze.shape, LR_Ze.min(), LR_Ze.max())
    print('          LR_mdv (dim,min,max) = ',LR_mdv.shape, LR_mdv.min(), LR_mdv.max())
    print('          LR_sw  (dim,min,max) = ',LR_sw.shape,  LR_sw.min(),  LR_sw.max())


    print('')
    print('   dim(mira_height,UTC_time_mira) = ',mira_height.shape,len(UTC_time_mira))
    print('            mira_Ze (dim,min,max)  = ',mira_Ze.shape, mira_Ze.min(), mira_Ze.max())
    print('            mira_mdv (dim,min,max) = ',mira_mdv.shape, mira_mdv.min(), mira_mdv.max())
    print('            mira_sw  (dim,min,max) = ',mira_sw.shape,  mira_sw.min(),  mira_sw.max())



####################################################################################################################
#
#       ######  ########    ###    ######## ####  ######  ######## ####  ######   ######
#       ##    ##    ##      ## ##      ##     ##  ##    ##    ##     ##  ##    ## ##    ##
#       ##          ##     ##   ##     ##     ##  ##          ##     ##  ##       ##
#        ######     ##    ##     ##    ##     ##   ######     ##     ##  ##        ######
#             ##    ##    #########    ##     ##        ##    ##     ##  ##             ##
#       ##    ##    ##    ##     ##    ##     ##  ##    ##    ##     ##  ##    ## ##    ##
#        ######     ##    ##     ##    ##    ####  ######     ##    ####  ######   ######
#
####################################################################################################################

# time averaged values
LR_timeavg_Ze = np.average(LR_Ze, axis=1)
LR_timeavg_mdv = np.average(LR_mdv, axis=1)
LR_timeavg_sw  = np.average(LR_sw,  axis=1)
mira_timeavg_Ze = np.average(mira_Ze, axis=1)
mira_timeavg_mdv = np.average(mira_mdv, axis=1)
mira_timeavg_sw  = np.average(mira_sw,  axis=1)

# height averaged values
LR_heightavg_Ze = np.average(LR_Ze, axis=0)
LR_heightavg_mdv = np.average(LR_mdv, axis=0)
LR_heightavg_sw  = np.average(LR_sw,  axis=0)
mira_heightavg_Ze = np.average(mira_Ze, axis=0)
mira_heightavg_mdv = np.average(mira_mdv, axis=0)
mira_heightavg_sw  = np.average(mira_sw,  axis=0)

# adapted means (under construction)

#t_mira, tp_mira, h_mira, \
#Z_mira, VEL_mira, RMS_mira = extract_dataset(comp_date, time_int, clock_time, '*mira.nc', '')
#
#t_mmclx, tp_mmclx, h_mmclx, \
#Ze_mmclx, VEL_mmclx, RMS_mmclx = extract_dataset(comp_date, time_int, clock_time, '*.mmclx', '')
#
#t_mmclx, tp_mmclx, h_mmclx, \
#Zg_mmclx, VELg_mmclx, RMSg_mmclx = extract_dataset(comp_date, time_int, clock_time, '*.mmclx', 'g')
#
#t_lr, tp_lr, h_lr, \
#Ze_lr, mdv_lr, sw_lr = extract_dataset(comp_date, time_int, clock_time, '*.LV1.NC', '')
#
#
#print('         MIRA                MMCLX               LIMRAD')
#for t1, t2, t3 in zip(tp_mira,tp_mmclx,tp_lr):
#    print(t1,' | ',t2,' | ',t3)
#
#print('         MIRA                MMCLX               LIMRAD')
#for h1, h2, h3 in zip(h_mira, h_mmclx, h_lr):
#    print('      ',h1, '     |     ', h2, '       |     ', h3)
#
#
#print('mira.nc\n')
#print('   dim(mira_height,UTC_time_mira) = ', h_mira.shape, len(tp_mira))
#print('            mira_Z (dim,min,max)   = ', Z_mira.shape, Z_mira.min(), Z_mira.max())
#print('            mira_VEL (dim,min,max) = ', VEL_mira.shape, VEL_mira.min(), VEL_mira.max())
#print('            mira_RMS (dim,min,max) = ', RMS_mira.shape, RMS_mira.min(), RMS_mira.max())
#
#
#print('.mmclx\n')
#print('   dim(mira_height,UTC_time_mira) = ', h_mmclx.shape, len(tp_mmclx))
#print('            mira_Ze (dim,min,max)  = ', Ze_mmclx.shape, Ze_mmclx.min(), Ze_mmclx.max())
#print('            mira_VEL (dim,min,max) = ', VEL_mmclx.shape, VEL_mmclx.min(), VEL_mmclx.max())
#print('            mira_RMS  (dim,min,max) = ', RMS_mmclx.shape, RMS_mmclx.min(), RMS_mmclx.max())
#
#
#print('.mmclx\n')
#print('   dim(mira_height,UTC_time_mira) = ', h_mmclx.shape, len(tp_mmclx))
#print('            mira_Ze (dim,min,max)  = ', Zg_mmclx.shape, Zg_mmclx.min(), Zg_mmclx.max())
#print('            mira_VEL (dim,min,max) = ', VELg_mmclx.shape, VELg_mmclx.min(), VELg_mmclx.max())
#print('            mira_RMS  (dim,min,max) = ', RMSg_mmclx.shape, RMSg_mmclx.min(), RMSg_mmclx.max())
#
#
#print('.LV1.NC\n')
#print('   dim(LR_height,UTC_time_LR) = ', h_lr.shape, len(tp_lr))
#print('          LR_Ze (dim,min,max)  = ', Ze_lr.shape, Ze_lr.min(), Ze_lr.max())
#print('          LR_mdv (dim,min,max) = ', mdv_lr.shape, mdv_lr.min(), mdv_lr.max())
#print('          LR_sw  (dim,min,max) = ', sw_lr.shape, sw_lr.min(), sw_lr.max())
#
#
#
#dummy = 5
#for ti in UTC_time_mira:
#    print('ti = ', ti)

#LR_idx_bins,   LR_time_bins   = gather_time_lists( UTC_time_LR , 40.0 )

#LR_dt_sums = []
#for sec in LR_time_bins:
#    LR_dt_sums.append(np.sum(sec))



#mira_idx_bins, mira_time_bins = gather_time_lists( UTC_time_mira , LR_dt_sums )

#for LR_s,MR_s in zip(LR_dt_sums,mira_time_bins):
#    print('MIRA indizes = ', LR_s, np.sum(MR_s))
#
#print( ' length of time series LR/Mira  =  ',len(LR_idx_bins),len(mira_idx_bins))

if plot_interp2d:

    npts =  500
    stat_pos    = [0.2 , -0.4]

    method = 'linear'

    # MIRA 35GHz Radar data extraction
    mira_timeNC, UTC_time_miraNC, mira_heightNC, \
    miraNC_Z, miraNC_VEL, miraNC_RMS = extract_dataset(comp_date, time_int, clock_time, '*mira.nc', '')

    _, _, _, \
    mira_Ze, mira_VEL, mira_RMS = extract_dataset(comp_date, time_int, clock_time, '*.mmclx', '')

    mira_time, UTC_time_mira, mira_height, \
    mira_Zg, mira_VELg, mira_RMSg = extract_dataset(comp_date, time_int, clock_time, '*.mmclx', 'g')

    LR_time, UTC_time_LR, LR_height, \
    LR_Ze, LR_mdv, LR_sw = extract_dataset(comp_date, time_int, clock_time, '*.LV1.NC', '')

    LR_unix_t   = [(ts.replace(tzinfo=timezone.utc).timestamp()) for ts in UTC_time_LR]
    mira_unix_t = [(ts.replace(tzinfo=timezone.utc).timestamp()) for ts in UTC_time_mira]


    #fcn_LR_to_MIRA_grid = interpolate.interp2d(LR_unix_t,   LR_height,   LR_Ze,   kind=method)
    #fcn_MIRA_to_LR_grid = interpolate.interp2d(mira_unix_t, mira_height, mira_Zg, kind=method)

    LRtoMIRA_Ze = Interpolate_NearestND_2D(LR_unix_t,   LR_height,   LR_Ze,   mira_unix_t, mira_height)
    MIRAtoLR_Zg = Interpolate_NearestND_2D(mira_unix_t, mira_height, mira_Ze, LR_unix_t,   LR_height)

    LR_Ze_mm6m3_I   = np.power(np.divide(LRtoMIRA_Ze,10.0), 10.0)
    LR_Ze_mm6m3     = np.power(np.divide(LR_Ze,10.0), 10.0)
    mira_Zg_mm6m3_I = np.power(np.divide(MIRAtoLR_Zg, 10.0), 10.0)
    mira_Zg_mm6m3   = np.power(np.divide(mira_Zg, 10.0), 10.0)

    mean_diff_LRtoM = np.mean(LRtoMIRA_Ze - mira_Zg)
    mean_diff_MtoLR = np.mean(LR_Ze - MIRAtoLR_Zg)

    #mean_diff_LRtoM = 10*np.log10( np.mean(LR_Ze_mm6m3_I - mira_Zg_mm6m3) )
    #mean_diff_MtoLR = 10*np.log10( np.mean(LR_Ze_mm6m3 - mira_Zg_mm6m3_I) )

    cor_coef = np.array([[]])

    for i in range(len(mira_unix_t)):
        df_LR   = pd.DataFrame(LRtoMIRA_Ze[:,i])
        df_mira = pd.DataFrame(mira_Zg[:,i])
        R = df_LR.corrwith(df_mira)
        cor_coef = np.append(cor_coef, R)


    correlation_LR_MIRA  = np.ma.masked_invalid(cor_coef)

    cor_coef = np.array([[]])
    for i in range(len(LR_unix_t)):
        df_LR   = pd.DataFrame(LR_Ze[:, i])
        df_mira = pd.DataFrame(MIRAtoLR_Zg[:, i])
        R = df_mira.corrwith(df_LR)
        cor_coef = np.append(cor_coef, R)

    correlation_MIRA_LR = np.ma.masked_invalid(cor_coef)

    #print('correlation stuff: ')
    #print('   max,min,size             ',cor_coef.max(),cor_coef.min())#,cor_coef.shape())
    #print('   LRtoMIRA_Zemax,min,size  ',LRtoMIRA_Ze.max(),LRtoMIRA_Ze.min(),LRtoMIRA_Ze)#.shape())
    #print('   mira_Zg max,min,size     ',mira_Zg.max(),mira_Zg.min())#,mira_Zg.shape())



    fig, ((p1,p2), (p3,p4), (p5, p6)) = plt.subplots(nrows=3, ncols=2, figsize=(12,8))
    ########################################################################################################
    # LR_Zelectivity plot
    if pts: print('')
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)


    p1.set_title(r'Radar Reflectivity Factor $Z_{e}^{94}$')
    plot_data_set(p1, '',
                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
                  x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')



    p2.set_title(r'Radar Reflectivity Factor $Z_{g}^{35}$')
    plot_data_set(p2, '',
                  UTC_time_mira, mira_height, mira_Zg, vmi=-50, vma=20,
                  x_min=UTC_time_mira[0], x_max=UTC_time_mira[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen


    if pts: print('       -   synchronized Radar Reflectivity Factor   ', end='', flush=True)
    # MIRA reflectivit
    p3.set_title(r'Interpolation of LIMRad $Z_{e}^{94}$'
                 'onto MIRA grid resolution\n'
                 r'$\delta h = 30$m, $\delta t = 2$sec', multialignment='center')

    place_statistics(p3, stat_pos, [mean_diff_LRtoM, np.mean(correlation_LR_MIRA)], 'Ze')

    plot_data_set(p3, '',
                  UTC_time_mira, mira_height, LRtoMIRA_Ze, vmi=-50, vma=20,
                  x_min=UTC_time_mira[1], x_max=UTC_time_mira[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='dBZ')

    p4.set_title(r'Interpolation of MIRA $Z_{g}^{35}$'
                 'onto LIMRad grid resolution\n'
                 r'$\delta h = 30$m, $\delta t = 5$sec', multialignment='center')

    place_statistics(p4, stat_pos, [mean_diff_MtoLR, np.mean(correlation_MIRA_LR)], 'Ze')

    plot_data_set(p4, '',
                  UTC_time_LR, LR_height, MIRAtoLR_Zg, vmi=-50, vma=20,
                  x_min=UTC_time_LR[1], x_max=UTC_time_LR[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='dBZ')


    if pts: print('\u2713')  # #print checkmark (✓) on screen


    if pts: print('       -   correnlation of Radar Reflectivity Factors   ', end='', flush=True)
    # MIRA reflectivit
    p5.set_title(r'Correlation of LIMRad interp$(Z_{e}^{94})$ and $Z_e^{35}$', multialignment='center')


    plot_correlation(p5, '', UTC_time_mira, correlation_LR_MIRA, '', '.',
                     x_min=UTC_time_mira[1], x_max=UTC_time_mira[-1],
                     y_min=-1.1, y_max=1.1, x_lab='Time (UTC)', y_lab='correlation')



    p6.set_title(r'Correlation of LIMRad interp$(Z_{g}^{35})$ and $Z_e^{94}$', multialignment='center')


    plot_correlation(p6, '', UTC_time_LR, correlation_MIRA_LR, '', '.',
                     x_min=UTC_time_mira[1], x_max=UTC_time_mira[-1],
                     y_min=-1.1, y_max=1.1, x_lab='Time (UTC)', y_lab='correlation')
    #plot_correlation_matrix(p5, '', cor_coef, vmi=1.e-8, vma=1, x_lab='', y_lab='', z_lab='')
    #p5.plot(UTC_time_mira,cor_coef)

    p5.set_title(r'Correlation of LIMRad interp$(Z_{e}^{35})$ and $Z_e^{94}$', multialignment='center')


    plt.tight_layout(rect=[0, 0.05, 1, 0.93])
    plt.subplots_adjust(hspace=0.9)

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line = 'Comparison of LIMRad 94 GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC), '
    third_line = 'using: *.LV1.NC and *.mmclx data (unprocessed datasets)'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.90])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.75)

    file = date_str + '_Interpolation_LR_MIRA.png'
    if pts: print('')
    if pts: print('    Save Figure to File :: ' + file + '\n')
    fig.savefig(file, dpi=dpi_val)

    plt.close()


if plot_RectBivariateSpline:
    # scipy interpolation test

    x = LR_time[:]
    y = LR_height[:]

    Z = LR_Ze[:,:]

    #grid the data
    interp_spline = RectBivariateSpline(x, y ,Z)

    # define grid
    x2 = np.arange( LR_time[0], LR_time[-1], 1 )
    y2 = np.arange(LR_height[0], LR_height[-1], 0.005 )
    #print('LR_height first/last = ',LR_height[0], LR_height[-1])

    time_x2 = []
    for i in range(len(x2)):
        time_x2.append(  datetime.datetime(2001,1,1,0,0,0)
                         + datetime.timedelta(seconds=int(x2[i])))

    X , Y  = np.meshgrid(x,y)
    X2, Y2 = np.meshgrid(x2,y2)

    Z2_unmasked = interp_spline(x2,y2)
    Z2 = np.ma.masked_less_equal(Z2_unmasked, -70.).T
    #Z2 = np.ma.masked_greater_equal(Z2, -999.).T

    #print('x=',x)
    #print('y=',y)
    #print('z=',Z,Z.min(),Z.max())
    #print('x2=',x2)
    #print('y2=',y2)


    #print('x2',x2.shape)
    #print('y2',y2.shape)
    #print('z2',Z2.shape,Z2.min(),Z2.max())

    fig, (LR_Ze_plot1,LR_Ze_plot2) = plt.subplots(nrows=2, ncols=1, figsize=(12,8))


    plot_data_set( LR_Ze_plot1 , 'Radar Reflectivity Factor' ,
                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
                  x_min=time_x2[0] , x_max=time_x2[-1] ,
                  y_min=hmin, y_max=hmax,
                  x_lab='', y_lab='height (km)', z_lab='dBZ')

    plot_data_set(LR_Ze_plot2, 'Radar Reflectivity Factor',
                  time_x2, y2 , Z2 , vmi=-50 , vma=20 ,
                  x_min=time_x2[0] , x_max=time_x2[-1],
                  y_min=hmin , y_max=hmax ,
                  x_lab='Time (UTC)'   , y_lab='height (km)' , z_lab='dBZ' )


    #for axes in ax:
    #    axes.set_zlim(-0.2,1)
    #    axes.set_axis_off()

    fig.tight_layout()

    plt.show()



####################################################################################################################
#
#           ######## #### ##     ## ########          ######  ######## ########  #### ########  ######
#              ##     ##  ###   ### ##               ##    ## ##       ##     ##  ##  ##       ##    ##
#              ##     ##  #### #### ##               ##       ##       ##     ##  ##  ##       ##
#              ##     ##  ## ### ## ######   #######  ######  ######   ########   ##  ######    ######
#              ##     ##  ##     ## ##                     ## ##       ##   ##    ##  ##             ##
#              ##     ##  ##     ## ##               ##    ## ##       ##    ##   ##  ##       ##    ##
#              ##    #### ##     ## ########          ######  ######## ##     ## #### ########  ######
#
####################################################################################################################

if plot_radar_results:
    ### plot ###
    print('    Generate subplots:\n')

    #create figure
    font = FontProperties()
    fig, ((LR_Ze_plot, mira_Ze_plot,) ,
        (LR_mdv_plot, mira_mdv_plot,)  ,
        (LR_sw_plot,  mira_sw_plot)) = plt.subplots(3,2,figsize=(16,12))


    x_lim_left_LR  = UTC_time_LR[0]
    x_lim_right_LR = UTC_time_LR[-1]

    x_lim_left_mira  = UTC_time_mira[0]
    x_lim_right_mira = UTC_time_mira[-1]


    ########################################################################################################
    ########################################################################################################
    #LR_Zelectivity plot
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    # LIMRad reflectivity
    LR_Ze_plot.set_title(r'\textbf{LIMRad94')

    plot_data_set( LR_Ze_plot , 'Radar Reflectivity Factor' ,
                   UTC_time_LR , LR_height , LR_Ze , vmi=-50 , vma=20 ,
                   x_min=x_lim_left_LR , x_max=x_lim_right_LR ,
                   y_min=hmin , y_max=hmax ,
                   x_lab=''   , y_lab='height (km)' , z_lab='dBZ'         )

    # MIRA reflectivit
    mira_Ze_plot.set_title(r'\textbf{MIRA35}')

    plot_data_set( mira_Ze_plot , 'Radar Reflectivity Factor' ,
                   UTC_time_mira , mira_height , mira_Ze , vmi=-50 , vma=20 ,
                   x_min=x_lim_left_mira , x_max=x_lim_right_mira ,
                   y_min=hmin , y_max=hmax ,
                   x_lab=''   , y_lab='Height (km)' , z_lab='dBZ'               )

    if pts: print('\u2713')  # #print checkmark (✓) on screen


    ########################################################################################################
    #mean doppler velocity plot
    if pts: print('       -   Mean Doppler velocity   ', end='', flush=True)

    # LIMRad mean Doppler velocity
    plot_data_set( LR_mdv_plot , 'Mean Doppler Velocity' ,
                   UTC_time_LR , LR_height , LR_mdv , vmi=-4 , vma=2 ,
                   x_min=x_lim_left_LR , x_max=x_lim_right_LR ,
                   y_min=hmin , y_max=hmax ,
                   x_lab=''   , y_lab='Height (km)' , z_lab='m/s'       )

    # MIRA mean Doppler velocity
    plot_data_set( mira_mdv_plot , 'Mean Doppler Velocity' ,
                   UTC_time_mira , mira_height , mira_mdv , vmi=-4 , vma=2 ,
                   x_min=x_lim_left_mira , x_max=x_lim_right_mira ,
                   y_min=hmin , y_max=hmax ,
                   x_lab=''   , y_lab='Height (km)' , z_lab='m/s'             )

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    #spectral width plot
    if pts: print('       -   Spectral Width   ', end='', flush=True)

    # LIMRad spectral width
    plot_data_set( LR_sw_plot , 'Spectral Width' ,
                   UTC_time_LR , LR_height , LR_sw , vmi=10**(-1.5) , vma=10**0.5 ,
                   x_min=x_lim_left_LR , x_max=x_lim_right_LR ,
                   y_min=hmin , y_max=hmax ,
                   x_lab='Time (UTC)' , y_lab='Height (km)' , z_lab='m/s'            )

    # MIRA spectral widthAAAAAA
    plot_data_set( mira_sw_plot , 'Spectral Width' ,
                   UTC_time_mira , mira_height , mira_sw , vmi=10**(-1.5) , vma=10**0.5 ,
                   x_min=x_lim_left_mira , x_max=x_lim_right_mira ,
                   y_min=hmin , y_max=hmax ,
                   x_lab='Time (UTC)' , y_lab='Height (km)' , z_lab='m/s'                  )

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen


    ########################################################################################################
    # Save figure to file
    date_str = str(plotyear)+str(plotmonth).zfill(2)+str(plotday).zfill(2)
    first_line  = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: '+str(time_int[0])+' (UTC)  to:  '+str(time_int[3])+' (UTC),'
    third_line  = r'using: ' + LIMRad_file_extension + ' and ' + mira_file_extension + ' data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.15)


    file = date_str + '_MIRA_LIMRad94_profiles_ts_comp.png'
    fig.savefig(file, dpi=dpi_val)
    plt.close()

    print('    Save Figure to File :: ' + file + '\n')


######################################################################################################################################
#
#  ######## #### ##     ## ########         ##     ## ######## ####  ######   ##     ## ########            ###    ##     ##  ######
#     ##     ##  ###   ### ##               ##     ## ##        ##  ##    ##  ##     ##    ##              ## ##   ##     ## ##    ##
#     ##     ##  #### #### ##               ##     ## ##        ##  ##        ##     ##    ##             ##   ##  ##     ## ##
#     ##     ##  ## ### ## ######   ####### ######### ######    ##  ##   #### #########    ##    ####### ##     ## ##     ## ##   ####
#     ##     ##  ##     ## ##               ##     ## ##        ##  ##    ##  ##     ##    ##            #########  ##   ##  ##    ##
#     ##     ##  ##     ## ##               ##     ## ##        ##  ##    ##  ##     ##    ##            ##     ##   ## ##   ##    ##
#     ##    #### ##     ## ########         ##     ## ######## ####  ######   ##     ##    ##            ##     ##    ###     ######
#
######################################################################################################################################

if plot_comparisons:

    ########################################################################
    ### plot comparison ###
    from matplotlib import colors
    print('    Generate subplots:\n')

    #create figure
    font = FontProperties()
    fig, ((LR_Ze_plot,        mira_Ze_plot,) ,
          (Comp_avgT_Ze_plot, Comp_avgH_Ze_plot,),
          (Comp_avgT_mdv_plot, Comp_avgH_mdv_plot,),
          (Comp_avgT_sw_plot,  Comp_avgH_sw_plot)) = plt.subplots(4,2,figsize=(16,12))


    # calculate x axis limits (same for both time series)

    x_lim_left_time  = min(UTC_time_LR[0],UTC_time_mira[0])
    x_lim_right_time = max(UTC_time_LR[-1],UTC_time_mira[-1])


    ################################################################################################################
    ################################################################################################################
    #
    # PLOT LIMRAD and MIRA Reflectivity (Ze) DATASET
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    # LIMRad reflectivity
    LR_Ze_plot.set_title(r'\textbf{LIMRad94}')

    plot_data_set( LR_Ze_plot , 'Radar Reflectivity Factor' ,
                   UTC_time_LR , LR_height , LR_Ze , vmi=-50 , vma=20 ,
                   x_min=x_lim_left_time , x_max=x_lim_right_time ,
                   y_min=hmin , y_max=hmax ,
                   x_lab='Time (UTC)', y_lab='height (km)' , z_lab='dBZ'  )

    # MIRA reflectivit
    mira_Ze_plot.set_title(r'\textbf{MIRA35}')

    plot_data_set( mira_Ze_plot , 'Radar Reflectivity Factor' ,
                   UTC_time_mira , mira_height , mira_Ze , vmi=-50 , vma=20 ,
                   x_min=x_lim_left_time , x_max=x_lim_right_time ,
                   y_min=hmin , y_max=hmax ,
                   x_lab='Time (UTC)' , y_lab='height (km)' , z_lab='dBZ'       )

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   comp. Radar Reflectivity Factor   ', end='', flush=True)

    x_lim_left_Ze  = min([LR_timeavg_Ze.min(),mira_timeavg_Ze.min()])
    x_lim_right_Ze = max([LR_timeavg_Ze.max(),mira_timeavg_Ze.max()])


    Comp_avgT_Ze_plot.set_title(r'\textbf{Time-Mean}')

    plot_avg_data_set( Comp_avgT_Ze_plot , 'Radar Reflectivity Factor' ,
                       LR_timeavg_Ze   , LR_height ,
                       mira_timeavg_Ze , mira_height ,
                       label1='LIMRad' , marker1='+' , label2='MIRA' ,   marker2='o' ,
                       x_min=x_lim_left_Ze , x_max=x_lim_right_Ze ,
                       y_min=hmin , y_max=hmax , x_lab='dBZ' , y_lab='height (km)', ax='y' )


    Comp_avgH_Ze_plot.set_title(r'\textbf{Height-Mean}')

    plot_avg_data_set( Comp_avgH_Ze_plot , 'Radar Reflectivity Factor' ,
                       UTC_time_LR   , LR_heightavg_Ze ,
                       UTC_time_mira , mira_heightavg_Ze ,
                       label1='LIMRad' , marker1='.' , label2='MIRA' ,   marker2='.' ,
                       x_min=x_lim_left_time , x_max=x_lim_right_time ,
                       y_min=[] , y_max=[] , x_lab='Time (UTC)' , y_lab='dBZ' , ax='y' )

    if pts: print('\u2713')  # #print checkmark (✓) on screen


    ################################################################################################################
    #
    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   comp. Mean Doppler Velocity   ', end='', flush=True)

    x_lim_left_mdv  = min([LR_timeavg_mdv.min(),mira_timeavg_mdv.min()])
    x_lim_right_mdv = max([LR_timeavg_mdv.max(),mira_timeavg_mdv.max()])


    plot_avg_data_set( Comp_avgT_mdv_plot , 'Doppler Velocity' ,
                       LR_timeavg_mdv   , LR_height ,
                       mira_timeavg_mdv , mira_height ,
                       label1='LIMRad' , marker1='+' , label2='MIRA' ,   marker2='o' ,
                       x_min=x_lim_left_mdv , x_max=x_lim_right_mdv ,
                       y_min=hmin , y_max=hmax , x_lab='m/s' , y_lab='height (km)', ax='y' )


    plot_avg_data_set( Comp_avgH_mdv_plot , 'Doppler Velocity' ,
                       UTC_time_LR   , LR_heightavg_mdv ,
                       UTC_time_mira , mira_heightavg_mdv ,
                       label1='LIMRad' , marker1='.' , label2='MIRA' ,   marker2='.' ,
                       x_min=x_lim_left_time , x_max=x_lim_right_time ,
                       y_min=[] , y_max=[] , x_lab='Time (UTC)' , y_lab='m/s' , ax='y' )

    if pts: print('\u2713')  # #print checkmark (✓) on screen


    ################################################################################################################
    #
    # comparsion of Spectral Width LIMRAD-MIRA
    if pts: print('       -   comp. Spectral Width   ', end='', flush=True)

    x_lim_left_sw  = min([LR_timeavg_sw.min(),mira_timeavg_sw.min()])
    x_lim_right_sw = max([LR_timeavg_sw.max(),mira_timeavg_sw.max()])


    plot_avg_data_set( Comp_avgT_sw_plot , 'Spectral Width' ,
                       LR_timeavg_sw   , LR_height ,
                       mira_timeavg_sw , mira_height ,
                       label1='LIMRad' , marker1='+' , label2='MIRA' ,   marker2='o' ,
                       x_min=x_lim_left_sw , x_max=x_lim_right_sw ,
                       y_min=hmin , y_max=hmax , x_lab='m/s' , y_lab='height (km)', ax='y' )

    plot_avg_data_set( Comp_avgH_sw_plot , 'Spectral Width' ,
                       UTC_time_LR   , LR_heightavg_sw ,
                       UTC_time_mira , mira_heightavg_sw ,
                       label1='LIMRad' , marker1='.' , label2='MIRA' ,   marker2='.' ,
                       x_min=x_lim_left_time , x_max=x_lim_right_time ,
                       y_min=[] , y_max=[] , x_lab='Time (UTC)' , y_lab='m/s', ax='y' )

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)


    # Save figure to file
    date_str = str(plotyear)+str(plotmonth).zfill(2)+str(plotday).zfill(2)
    first_line  = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: '+str(time_int[0])+' (UTC)  to:  '+str(time_int[3])+' (UTC),'
    third_line  = r'using: ' + LIMRad_file_extension + ' and ' + mira_file_extension + ' data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)


    plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.35)

    file = date_str + '_MIRA_LIMRad94_t-h_comp.png'
    print('    Save Figure to File :: ' + file + '\n')
    fig.savefig(file, dpi=dpi_val)
    plt.close()




if plot_interpolation_scatter:
    ########################################################################
    ### plot comparison ###


    interp_meth = 'nearest'
    res_interp  = 5 # in [sec]
    head_pos    = [0.05]
    stat_pos    = [0.25 , -0.3]


    # preparations for interpolation plots

    LR_unix_t   = [(ts.replace(tzinfo=timezone.utc).timestamp()) for ts in UTC_time_LR]
    mira_unix_t = [(ts.replace(tzinfo=timezone.utc).timestamp()) for ts in UTC_time_mira]


    # create an array with evenly spaced gridsize
    xnew = np.arange(max(LR_unix_t[0], mira_unix_t[0]) ,
                     min(LR_unix_t[-1],mira_unix_t[-1]),
                     res_interp)


    # convert the x-axis unix time to readable date time format
    plot_time_xnew = [datetime.datetime(1970, 1, 1, 0, 0, 0)
                   + datetime.timedelta(seconds=int(xnew[i])) for i in range(len(xnew))]


    # what do you want to interpolate
    LR_x   = LR_unix_t;        mira_x = mira_unix_t
    LR_y   = LR_heightavg_Ze;  mira_y = mira_heightavg_Ze


    LR_ynew   = interpolate_data(LR_x, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(mira_x, mira_y, xnew, interp_meth)


    # y width +-5
    y_min, y_max = get_plot_ybounds(LR_y, mira_y, 7.0)

    # calculate the mean difference and covariance matrix
    mean_diff = np.mean(np.absolute(LR_ynew-mira_ynew))
    cor_coef  = np.corrcoef(LR_ynew,mira_ynew)


    print('    Generate subplots:\n')

    font = FontProperties()
    #fig, ((h_Ze_plot,  interp_h_Ze_plot,  scatter_Ze),
    #      (h_mdv_plot, interp_h_mdv_plot, scatter_mdv),
    #      (h_sw_plot,  interp_h_sw_plot,  scatter_sw)) = plt.subplots(3,3,figsize=(16,12))

    fig = plt.figure(figsize=(16,12))

    h_Ze_plot  = plt.subplot2grid((3, 3), (0, 0))
    h_mdv_plot = plt.subplot2grid((3, 3), (0, 1))#, rowspan=2)
    h_sw_plot  = plt.subplot2grid((3, 3), (0, 2))#, rowspan=2)
    interp_h_Ze_plot  = plt.subplot2grid((3, 3), (1, 0))#, colspan=2)
    interp_h_mdv_plot = plt.subplot2grid((3, 3), (1, 1))#, rowspan=2)
    interp_h_sw_plot  = plt.subplot2grid((3, 3), (1, 2))#, rowspan=2)
    scatter_Ze  = plt.subplot2grid((3, 3), (2, 0))#, colspan=2, rowspan=2)
    scatter_mdv = plt.subplot2grid((3, 3), (2, 1))#, rowspan=2)
    scatter_sw  = plt.subplot2grid((3, 3), (2, 2))#, rowspan=2)

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   Average reflectivity over height domain  ', end='', flush=True)


    h_Ze_plot.set_title(r' \textbf{Mean-Height Reflectivity}')
    plot_avg_data_set( h_Ze_plot , '' ,
                       UTC_time_LR   , LR_heightavg_Ze ,
                       UTC_time_mira , mira_heightavg_Ze ,
                       label1='LIMRad' , marker1='.' , label2='MIRA' ,   marker2='.' ,
                       x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                       y_min=y_min , y_max=y_max , x_lab='Time (UTC)' , y_lab='dBZ'  , ax='n' )


    interp_h_Ze_plot.set_title(r' \textbf{Mean-Height Reflectivity}')
    plot_interpol_data_set(interp_h_Ze_plot, '',
                           plot_time_xnew, LR_ynew,
                           plot_time_xnew, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                           y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='dBZ')

    place_statistics(interp_h_Ze_plot, stat_pos, [mean_diff,cor_coef[0,1]],'Ze')

    scatter_Ze.set_title('Scatter Plot of Height-Mean\n Reflectivity ')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 7.0)

    plot_scatter( scatter_Ze, '', LR_ynew, mira_ynew, '*'
                , x_min=xy_min, x_max=xy_max,
                  y_min=xy_min, y_max=xy_max,
                  x_lab='LIMRad data in dBZ', y_lab='MIRA data in dBZ')



    if pts: print('\u2713')  # #print checkmark (✓) on screen)


    # same for mean doppler velocity

    LR_y   = LR_heightavg_mdv
    mira_y = mira_heightavg_mdv

    LR_ynew   = interpolate_data(LR_x, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(mira_x, mira_y, xnew, interp_meth)
    LR_ynew   = np.ma.masked_equal(LR_ynew,0.0)
    mira_ynew = np.ma.masked_equal(mira_ynew,0.0)

    # y width +-5
    y_min, y_max = get_plot_ybounds(LR_y, mira_y, 0.1)

    # calculate the mean difference and covariance matrix
    mean_diff = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef = np.corrcoef(LR_ynew, mira_ynew)

    ################################################################################################################
    #
    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   Average mean doppler velocity over height domain  ', end='', flush=True)

    place_text(h_mdv_plot, [-1.1, 1.2], r'\LARGE{\textbf{Actual Dataset}}')

    h_mdv_plot.set_title(r' \textbf{Height-Mean Mean Doppler Velocity}')
    plot_avg_data_set( h_mdv_plot , '' ,
                       UTC_time_LR   , LR_heightavg_mdv ,
                       UTC_time_mira , mira_heightavg_mdv ,
                       label1='LIMRad' , marker1='.' , label2='MIRA' ,   marker2='.' ,
                       x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                       y_min=y_min , y_max=y_max , x_lab='Time (UTC)' , y_lab='m/s'  , ax='n' )

    place_text(interp_h_mdv_plot, [-1.1, 1.2], r'\LARGE{\textbf{Interpolated Dataset}}')

    interp_h_mdv_plot.set_title(r' \textbf{Height-Mean Mean Doppler Velocity}')
    plot_interpol_data_set(interp_h_mdv_plot, '',
                           plot_time_xnew, LR_ynew,
                           plot_time_xnew, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                           y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='m/s')


    place_statistics(interp_h_mdv_plot,  stat_pos, [mean_diff,cor_coef[0,1]],'mdv')
    scatter_mdv.set_title('Scatter Plot of Height-Mean\n Mean Doppler Velocity')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.1)

    plot_scatter(scatter_mdv, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max,
                 y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in m/s', y_lab='MIRA data in m/s')


    if pts: print('\u2713')  # #print checkmark (✓) on screen)




    # same for mean doppler velocity

    LR_y   = LR_heightavg_sw
    mira_y = mira_heightavg_sw

    LR_ynew   = interpolate_data(LR_x, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(mira_x, mira_y, xnew, interp_meth)

    # y width +-5
    y_min, y_max = get_plot_ybounds(LR_y, mira_y, 0.03)

    # calculate the mean difference and covariance matrix
    mean_diff = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef = np.corrcoef(LR_ynew, mira_ynew)

    ################################################################################################################
    #
    # comparsion of spectral width LIMRAD-MIRA
    if pts: print('       -   Average spectral width over height domain  ', end='', flush=True)

    h_sw_plot.set_title(r'\textbf{Height-Mean Spectral Width}')
    plot_avg_data_set( h_sw_plot , '' ,
                       UTC_time_LR   , LR_heightavg_sw ,
                       UTC_time_mira , mira_heightavg_sw ,
                       label1='LIMRad' , marker1='.' , label2='MIRA' ,   marker2='.' ,
                       x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                       y_min=y_min , y_max=y_max , x_lab='Time (UTC)' , y_lab='m/s' , ax='n' )


    interp_h_sw_plot.set_title(r'\textbf{Height-Mean Spectral Width}')
    plot_interpol_data_set(interp_h_sw_plot, '',
                           plot_time_xnew, LR_ynew,
                           plot_time_xnew, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                           y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='m/s')

    place_statistics(interp_h_sw_plot, stat_pos, [mean_diff,cor_coef[0,1]],'sw')

    scatter_sw.set_title('Scatter Plot of Height-Mean\n Spectral Width ')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.03)

    plot_scatter(scatter_sw, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max,
                 y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line  = 'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:x  ' + str(time_int[3]) + ' (UTC), '
    third_line  = 'using: ' + LIMRad_file_extension + ' and ' + mira_file_extension + ' data,  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{'+third_line+'}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.91])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.55)

    file = date_str + '_MIRA_LIMRad94_interp-avgheight_comp.png'
    print('    Save Figure to File :: ' + file + '\n')
    fig.savefig(file, dpi=dpi_val)
    plt.close()

    save_log_data(file[:-5], interp_meth, res_interp)


if plot_compare_mira_mmclx:





    #print('')
    #print('Zh_mira-Ze_mmclx = ',np.sum(miramira_Z-mira_Ze[:,:-1]))

    print('')
    print('    Generate subplots:\n')

    font = FontProperties()

    fig = plt.figure(figsize=(16, 9))

    Zh_plot = plt.subplot2grid((3, 4), (0, 0))
    Ze_plot = plt.subplot2grid((3, 4), (0, 1))
    Zg_plot  = plt.subplot2grid((3, 4), (0, 2))  # , rowspan=2)
    LR_Ze_plot = plt.subplot2grid((3, 4), (0, 3))  # , rowspan=2)

    v_plot     = plt.subplot2grid((3, 4), (1, 0))  # , colspan=2)
    VEL_plot   = plt.subplot2grid((3, 4), (1, 1))  # , colspan=2)
    VELg_plot = plt.subplot2grid((3, 4), (1, 2))  # , rowspan=2)
    LR_mdv_plot  = plt.subplot2grid((3, 4), (1, 3))  # , rowspan=2)

    sw_plot    = plt.subplot2grid((3, 4), (2, 0))  # , colspan=2, rowspan=2)
    RMS_plot   = plt.subplot2grid((3, 4), (2, 1))  # , colspan=2, rowspan=2)
    RMSg_plot = plt.subplot2grid((3, 4), (2, 2))  # , rowspan=2)
    LR_sw_plot  = plt.subplot2grid((3, 4), (2, 3))  # , rowspan=2)

    height_describtion = [-0.05, 1.3]
    place_text(Zh_plot, height_describtion,  r'\textbf{Zh} ... Calibrated reflectivity. Calibration\\'
                                              ' convention: in the absence of attenuation, \n'
                                              'a cloud at 273 K containing one million 100-micron\n'
                                              ' droplets per cubic metre will have a reflectivity \n'
                                              'of 0 dBZ at all frequencies)\n'
                                             r'\textbf{v}  ... Radial component of the velocity,\\'
                                              ' with positive velocities are away from the radar.)\n '
                                             r'\textbf{width} ... Standard deviation of the reflectivity-\\'
                                             r'weighted velocities in the radar pulse volume.)' )
    height_describtion = [0.1, 1.3]

    place_text(Ze_plot, height_describtion,  r'\textbf{Ze} ... Equivalent Radar Reflectivity\\'
                                             ' Factor Ze of Hydrometeors  \n'
                                             r'\textbf{VEL} ... Doppler Velocity VEL\\'
                                             r'\textbf{RMS} ... Peak Width RMS' )
    place_text(Zg_plot, height_describtion,  r'\textbf{Zg} ... Equivalent Radar Reflectivity\\'
                                             ' Factor Ze of all Targets \n'
                                             r'\textbf{VELg} ... Doppler Velocity VELg\\'
                                             r'\textbf{RMSg} ... Peak Width RMSg')
    place_text(LR_Ze_plot, height_describtion,   r'\textbf{Ze} ... Equivalent radar reflectivity factor\\'
                                             r'\textbf{mdv} ... Mean Doppler Velocity\\'
                                             r'\textbf{sw} ... Spectrum width')
    ################################################################################################################
    #
    ref_min = -60
    ref_max = 30
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   MIRA35 Reflectivity (Z)  ', end='', flush=True)
    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]

    Zh_plot.set_title(r'\textbf{Radar Reflectivity Factor (Zh)}')
    plot_data_set( Zh_plot , '' ,
                   UTC_time_mira , mira_height , miramira_Z , vmi=ref_min , vma=ref_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='dBZ' )

    Ze_plot.set_title(r'\textbf{Radar Reflectivity Factor (Ze)}')
    plot_data_set( Ze_plot , '' ,
                   UTC_time_mira , mira_height , mira_Ze , vmi=ref_min , vma=ref_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='dBZ' )


    Zg_plot.set_title(r'\textbf{Radar Reflectivity Factor (Zg)}')
    plot_data_set( Zg_plot , '' ,
                   UTC_time_mira , mira_height , mira_Zg , vmi=ref_min , vma=ref_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='dBZ' )

    LR_Ze_plot.set_title(r'\textbf{Radar Reflectivity Factor (Ze)}')
    plot_data_set( LR_Ze_plot , '' ,
                   UTC_time_LR , LR_height , LR_Ze , vmi=ref_min , vma=ref_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='dBZ' )

    if pts: print('\u2713')  # #print checkmark (✓) on screen)


    ################################################################################################################
    #
    vel_min = -5
    vel_max = 3
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   MIRA35 Mean Doppler Velocity (VEL)  ', end='', flush=True)
    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]

    v_plot.set_title(r'\textbf{Mean Doppler Velocity (v)}')
    plot_data_set( v_plot , '' ,
                   UTC_time_mira , mira_height , miramira_VEL , vmi=vel_min , vma=vel_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    VEL_plot.set_title(r'\textbf{Mean Doppler Velocity (VEL)}')
    plot_data_set( VEL_plot , '' ,
                   UTC_time_mira , mira_height , mira_VEL , vmi=vel_min , vma=vel_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    VELg_plot.set_title(r'\textbf{Mean Doppler Velocity (VELg)}')
    plot_data_set( VELg_plot , '' ,
                   UTC_time_mira , mira_height , mira_VELg , vmi=vel_min , vma=vel_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    LR_mdv_plot.set_title(r'\textbf{Mean Doppler Velocity (mdv)}')
    plot_data_set( LR_mdv_plot , '' ,
                   UTC_time_LR , LR_height , LR_mdv , vmi=vel_min , vma=vel_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    if pts: print('\u2713')  # #print checkmark (✓) on screen)



    ################################################################################################################
    #
    sw_min = 10**(-1.5)
    sw_max = 10**(0.5)
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   MIRA35 Spectral Width (RMS)  ', end='', flush=True)
    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]

    sw_plot.set_title(r'\textbf{Spectral Width (width)}')
    plot_data_set( sw_plot , 'sw' ,
                   UTC_time_mira , mira_height , miramira_RMS , vmi=sw_min , vma=sw_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    RMS_plot.set_title(r'\textbf{Spectral Width (RMS)}')
    plot_data_set( RMS_plot , 'sw' ,
                   UTC_time_mira , mira_height , mira_RMS , vmi=sw_min , vma=sw_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    RMSg_plot.set_title(r'\textbf{Spectral Width (RMSg)}')
    plot_data_set( RMSg_plot , 'sw' ,
                   UTC_time_mira , mira_height , mira_RMSg , vmi=sw_min , vma=sw_max ,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    LR_sw_plot.set_title(r'\textbf{Spectral Width (sw)}')
    plot_data_set( LR_sw_plot , 'sw' ,
                   UTC_time_LR , LR_height , LR_sw, vmi=sw_min , vma=sw_max,
                   x_min=x_lim[0] , x_max=x_lim[-1]  ,
                   y_min=hmin ,     y_max=hmax ,
                   x_lab='Time (UTC)'   ,     y_lab='height (km)' , z_lab='m/s' )

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)


    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line  = 'Comparison of processed and unprocessed MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC), '
    third_line  = 'using: *mira.nc data (1st column) and *.mmclx data (2nd - 4th column; no attenuation correction)'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{'+third_line+'}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.80])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)

    file = date_str + '_MIRA_mmclx_comparison.png'
    print('    Save Figure to File :: ' + file + '\n')
    fig.savefig(file, dpi=dpi_val)

    plt.close()






print('    Elapsed Time = {0:0.3f}'.format(time.clock()-start_time),'[sec]\n')
