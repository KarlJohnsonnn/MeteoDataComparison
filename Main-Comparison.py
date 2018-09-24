import NetCDF_Tool  as nc
from Graphics_Mod import *
from IO_Mod import *
from Parameter_Mod import *

rc('font', family='serif')
rc('text', usetex=True)

import sys, warnings, time



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
plot_radar_results = False
plot_comparisons = False
plot_interpolation_scatter = False
plot_interp2d = False
#plot_compare_mira_mmclx = False


os.chdir(meteo_path)


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
print('  \u250F' + 49 * '\u2501' + '\u2513')
print('  \u2503' + '          LIMRAD94 - MIRA35  Comparison          ' + '\u2503')
print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
print('')
print('')

# gather arguments
if len(sys.argv) == 6:

    hmin, hmax = float(sys.argv[1]), float(sys.argv[2])
    date = str(sys.argv[3])
    time_intervall = str(sys.argv[4]) + '-' + str(sys.argv[5])

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

    hmin = 3.0  # (km)  - lower y-axis limit
    hmax = 8.00  # (km) - upper y-axis limit, highest range gate may be higher
    date     = '180729'     # in YYMMDD
    time_intervall = '2200-2359'  # in HHMM-HHMM

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



# ----- LIMRad 94GHz Radar data extraction
LR_data = nc.LIMRad94_LV1(date, time_intervall, [hmin, hmax])

if interpolate_cn: LR_data.interpolate_cn(t_res=interp_time_res, r_res=interp_range_res, method='constant')
if create_nc_file:
    LR_data.save(LIMRad_path)
    exit(0)

# ----- MIRA 35GHz Radar data extraction
MIRA_data  = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax])
MMCLX_data = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax], '*.mmclx')


print('')


####################################################################################################################
#
#        ######  ########    ###    ######## ####  ######  ######## ####  ######   ######
#       ##    ##    ##      ## ##      ##     ##  ##    ##    ##     ##  ##    ## ##    ##
#       ##          ##     ##   ##     ##     ##  ##          ##     ##  ##       ##
#        ######     ##    ##     ##    ##     ##   ######     ##     ##  ##        ######
#             ##    ##    #########    ##     ##        ##    ##     ##  ##             ##
#       ##    ##    ##    ##     ##    ##     ##  ##    ##    ##     ##  ##    ## ##    ##
#        ######     ##    ##     ##    ##    ####  ######     ##    ####  ######   ######
#
####################################################################################################################

# time averaged values
    # calculate Ze, mdv, sw (time averaged)
LR_data.avg_time()
MIRA_data.avg_time()
MMCLX_data.avg_time()

    # calculate Ze, mdv, sw (height averaged)
LR_data.avg_height()
MIRA_data.avg_height()
MMCLX_data.avg_height()


####################################################################################################################
#
#       ######   ########     ###    ########  ##     ## ####  ######   ######
#      ##    ##  ##     ##   ## ##   ##     ## ##     ##  ##  ##    ## ##    ##
#      ##        ##     ##  ##   ##  ##     ## ##     ##  ##  ##       ##
#      ##   #### ########  ##     ## ########  #########  ##  ##        ######
#      ##    ##  ##   ##   ######### ##        ##     ##  ##  ##             ##
#      ##    ##  ##    ##  ##     ## ##        ##     ##  ##  ##    ## ##    ##
#       ######   ##     ## ##     ## ##        ##     ## ####  ######   ######
#
####################################################################################################################

if plot_interp2d:

    fig, plt = Plot_2D_Interpolation(LR_data, MMCLX_data)

    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_2D-Interpolation.png'
    if pts: print('')
    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(meteo_path + file, dpi=dpi_val)
    plt.close()


if plot_radar_results:

    fig, plt = Plot_Radar_Results(LR_data, MMCLX_data)

    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_profiles_timeseries.png'
    fig.savefig(meteo_path + file, dpi=dpi_val)
    plt.close()

    print('    Save Figure to File :: ' + meteo_path + file + '\n')


if plot_comparisons:

    fig, plt = Plot_Comparison(LR_data, MMCLX_data)

    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_avg_time_height.png'
    print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(meteo_path + file, dpi=dpi_val)
    plt.close()

if plot_interpolation_scatter:

    fig, plt = Plot_Scatter(LR_data, MMCLX_data)

    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_interp-avgheight_comp.png'
    print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(meteo_path + file, dpi=dpi_val)
    plt.close()



#save_log_data(file[:-5], interp_meth, hmin, hmax, date, time_intervall)
print('    Elapsed Time = {0:0.3f}'.format(time.clock() - start_time), '[sec]\n')
