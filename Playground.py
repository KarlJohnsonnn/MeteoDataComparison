import sys

import warnings

import modules.NetCDF_Mod as nc
from modules.PlotLibrary_Mod import *
from modules.Utility_Mod import *

'''
##################################################################################################

       ##     ##  ######  ######## ########     #### ##    ## ########  ##     ## ########
       ##     ## ##    ## ##       ##     ##     ##  ###   ## ##     ## ##     ##    ##
       ##     ## ##       ##       ##     ##     ##  ####  ## ##     ## ##     ##    ##
       ##     ##  ######  ######   ########      ##  ## ## ## ########  ##     ##    ##
       ##     ##       ## ##       ##   ##       ##  ##  #### ##        ##     ##    ##
       ##     ## ##    ## ##       ##    ##      ##  ##   ### ##        ##     ##    ##
        #######   ######  ######## ##     ##    #### ##    ## ##         #######     ##

##################################################################################################
'''
# Logicals for different tasks
calc_doppler_spectra = False
plot_radar_results = True
plot_compare_noise = False
plot_for_poster = False
plot_comparisons     = False
plot_interp2d        = False
plot_interpolation_scatter = False

interpolate_cn = False
create_nc_file = False

'''
####################################################################################################################

   ##     ##    ###    #### ##    ##         ########  ########   #######   ######   ########     ###    ##     ##
   ###   ###   ## ##    ##  ###   ##         ##     ## ##     ## ##     ## ##    ##  ##     ##   ## ##   ###   ###
   #### ####  ##   ##   ##  ####  ##         ##     ## ##     ## ##     ## ##        ##     ##  ##   ##  #### ####
   ## ### ## ##     ##  ##  ## ## ## ####### ########  ########  ##     ## ##   #### ########  ##     ## ## ### ##
   ##     ## #########  ##  ##  ####         ##        ##   ##   ##     ## ##    ##  ##   ##   ######### ##     ##
   ##     ## ##     ##  ##  ##   ###         ##        ##    ##  ##     ## ##    ##  ##    ##  ##     ## ##     ##
   ##     ## ##     ## #### ##    ##         ##        ##     ##  #######   ######   ##     ## ##     ## ##     ##

####################################################################################################################
'''
start_time = time.clock()

# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '          LIMRAD94 - MIRA35  Comparison          ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n'*2)


# gather arguments
if len(sys.argv) == 6:
    date = str(sys.argv[1])
    time_intervall = str(sys.argv[2]) + '-' + str(sys.argv[3])
    hmin, hmax = float(sys.argv[4]), float(sys.argv[5])

else:

    ## cirrus
    #hmin = 8.50  # (km)  - lower y-axis limit
    #hmax = 10.0  # (km) - upper y-axis limit, highest range gate may be higher
    #date = '180728'  # in YYMMDD
    #time_intervall = '0740-0805'  # in HHMM-HHMM

    ##cummulis
    #hmin = 5.0  #(km)  - lower y-axis limit
    #hmax = 12.0  #(km) - upper y-axis limit, highest range gate may be higher
    #date     = '180802'     # in YYMMDD
    #time_intervall = '0330-1200'  # in HHMM-HHM

    #hmin = 0.0  # (km)  - lower y-axis limit
    #hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    #date = '20180729'  # in YYMMDD
    #time_intervall = '000000-235959'  # in HHMM-HHMM

    hmin = 0.0  # (km)  - lower y-axis limit
    hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    date = '20181209'  # in YYMMDD
    time_intervall = '000000-030000'  # in HHMM-HHMM

    # hmin = 0.0  # (km)  - lower y-axis limit
    # hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    # date = '180810'  # in YYMMDD
    # time_intervall = '0500-0600'  # in HHMM-HHMM

    ##nimbus
    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180805'     # in YYMMDD
    # time_intervall = '0510-0620'  # in HHMM-HHMM

    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180805'     # in YYMMDD
    # time_intervall = '1030-1200'  # in HHMM-HHMM

    # hmin = 1.0 #(km)  - lower y-axis limit
    # hmax = 2.5 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180805'     # in YYMMDD
    # time_intervall = '0700-1210'  # in HHMM-HHMM

    ##nimbus
    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 12.00 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180808'     # in YYMMDD
    # time_intervall = '1330-1700'  # in HHMM-HHMM

warnings.filterwarnings("ignore")

'''
######################################################################################################

   ########     ###    ########    ###                  #### ##    ## ########  ##     ## ########
   ##     ##   ## ##      ##      ## ##                  ##  ###   ## ##     ## ##     ##    ##
   ##     ##  ##   ##     ##     ##   ##                 ##  ####  ## ##     ## ##     ##    ##
   ##     ## ##     ##    ##    ##     ##    #######     ##  ## ## ## ########  ##     ##    ##
   ##     ## #########    ##    #########                ##  ##  #### ##        ##     ##    ##
   ##     ## ##     ##    ##    ##     ##                ##  ##   ### ##        ##     ##    ##
   ########  ##     ##    ##    ##     ##               #### ##    ## ##         #######     ##

######################################################################################################
'''


# ----- LIMRAD 94GHz Radar data extraction
# special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
#hmin = 0.0  # (km)  - lower y-axis limit
#hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
#date = '20180810'  # in YYMMDD
#time_intervall = '050000-060000'  # in HHMM-HHMM

#  special case NoiseFac0_file = 'NOISEFAC0_180820_142451_P01_ZEN.LV0.NC'
# hmin = 0.0  # (km)  - lower y-axis limit
# hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
# date = '180820'  # in YYMMDD
# time_intervall = '1400-1500'  # in HHMM-HHMM


print('    date: ', date, time_intervall, hmin, hmax)

#LR_lv0 = nc.LIMRAD94_LV0(date, time_intervall, [hmin, hmax])
LR_lv1 = nc.LIMRAD94_LV1(LIMRAD_path, date, time_intervall, [hmin, hmax])

# if interpolate_cn: LR_lv1.interpolate_cn(t_res=interp_time_res, r_res=interp_range_res, method='constant')
# if create_nc_file: LR_lv1.save('/Users/willi/Desktop/')



# ----- MIRA 35GHz Raar data extraction
#MIRA_data  = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax])
MMCLX_data = nc.MIRA35_LV1(MIRA_path, date, time_intervall, [hmin, hmax], '*.mmclx')

#if interpolate_cn: MMCLX_data.interpolate_cn(t_res=interp_time_res, r_res=interp_range_res, method='constant')

if pts: print('')

'''
####################################################################################################################

        ######  ########    ###    ######## ####  ######  ######## ####  ######   ######
       ##    ##    ##      ## ##      ##     ##  ##    ##    ##     ##  ##    ## ##    ##
       ##          ##     ##   ##     ##     ##  ##          ##     ##  ##       ##
        ######     ##    ##     ##    ##     ##   ######     ##     ##  ##        ######
             ##    ##    #########    ##     ##        ##    ##     ##  ##             ##
       ##    ##    ##    ##     ##    ##     ##  ##    ##    ##     ##  ##    ## ##    ##
        ######     ##    ##     ##    ##    ####  ######     ##    ####  ######   ######

####################################################################################################################
'''
# averaged values
    # calculate Ze, mdv, sw (time averaged)
#LR_lv1.avg_time()
# MIRA_data.avg_time()
#MMCLX_data.avg_time()

    # calculate Ze, mdv, sw (height averaged)
#LR_lv1.avg_height()
# MIRA_data.avg_height()
#MMCLX_data.avg_height()


'''
####################################################################################################################

       ######   ########     ###    ########  ##     ## ####  ######   ######
      ##    ##  ##     ##   ## ##   ##     ## ##     ##  ##  ##    ## ##    ##
      ##        ##     ##  ##   ##  ##     ## ##     ##  ##  ##       ##
      ##   #### ########  ##     ## ########  #########  ##  ##        ######
      ##    ##  ##   ##   ######### ##        ##     ##  ##  ##             ##
      ##    ##  ##    ##  ##     ## ##        ##     ##  ##  ##    ## ##    ##
       ######   ##     ## ##     ## ##        ##     ## ####  ######   ######

####################################################################################################################
'''

if plot_for_poster:

    fig, plt = Plot_for_poster(LR_lv1)

    file = meteo_path + date + '_' + time_intervall + '_radar_LIMRAD94.png'
    fig.savefig(file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('')
    if pts: print('    Save Figure to File :: ' + file + '\n')

if plot_interp2d:

    fig, plt = Plot_2D_Interpolation(LR_lv1, MMCLX_data)

    file = date + '_2D-Interpolation.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('')
    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


if plot_radar_results:

    fig, plt = Plot_Radar_Results(LR_lv1, MMCLX_data)

    file = date + '_profiles_timeseries.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


if plot_comparisons:

    fig, plt = Plot_Comparison(LR_lv1, MMCLX_data)

    file = date + '_avg_time_height.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

if plot_interpolation_scatter:

    fig, plt = Plot_Scatter(LR_lv1, MMCLX_data)

    file = date + '_interp-avgheight_comp.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


#save_log_data(file[:-5], interp_meth, hmin, hmax, date, time_intervall)

if pts: print(f'    Total Elapsed Time = {time.clock()-start_time:.3f} sec.\n')
