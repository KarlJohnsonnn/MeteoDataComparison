import NetCDF_Tool  as nc
from Functions_Mod import *
from Graphics_Mod import *

rc('font', family='serif')
rc('text', usetex=True)

import sys, warnings, time, os



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
plot_radar_results = True
plot_for_poster = False
plot_comparisons = False
plot_interpolation_scatter = False
plot_interp2d = False
plot_doppler_spectra = False
#plot_compare_mira_mmclx = False

interpolate_cn = False
create_nc_file = False

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
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '          LIMRAD94 - MIRA35  Comparison          ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)


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

    hmin = 0.0  # (km)  - lower y-axis limit
    hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    date = '180729'  # in YYMMDD
    time_intervall = '0000-0100'  # in HHMM-HHMM

#    hmin = 0.0  # (km)  - lower y-axis limit
#    hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
#    date = '180810'  # in YYMMDD
#    time_intervall = '0500-0600'  # in HHMM-HHMM

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


# ----- LIMRAD 94GHz Radar data extraction
# NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'

hmin_lv0 = 0.0  # (km)  - lower y-axis limit
hmax_lv0 = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
date_lv0 = '180810'  # in YYMMDD
time_intervall_lv0 = '0500-0600'  # in HHMM-HHMM

LR_lv0 = nc.LIMRAD94_LV0(date_lv0, time_intervall_lv0, [hmin_lv0, hmax_lv0])

LR_lv1 = nc.LIMRAD94_LV1(date, time_intervall, [hmin, hmax])

# if interpolate_cn: LR_lv1.interpolate_cn(t_res=interp_time_res, r_res=interp_range_res, method='constant')
# if create_nc_file: LR_lv1.save(LIMRAD_path)



# ----- MIRA 35GHz Radar data extraction
# MIRA_data  = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax])
MMCLX_data = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax], '*.mmclx')

if pts: print('')


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
# LR_lv1.avg_time()
#MIRA_data.avg_time()
#MMCLX_data.avg_time()

    # calculate Ze, mdv, sw (height averaged)
#LR_lv1.avg_height()
#MIRA_data.avg_height()
#MMCLX_data.avg_height()


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

date_str = str(LR_lv1.year) + str(LR_lv1.month).zfill(2) + str(LR_lv1.day).zfill(2)


if plot_for_poster:

    fig, plt = Plot_for_poster(LR_lv1)

    file = meteo_path + date_str + '_radar_LIMRAD94_vergleich.png'
    fig.savefig(file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('')
    if pts: print('    Save Figure to File :: ' + file + '\n')

if plot_interp2d:

    fig, plt = Plot_2D_Interpolation(LR_lv1, MMCLX_data)

    file = date_str + '_2D-Interpolation.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('')
    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


if plot_radar_results:

    fig, plt = Plot_Radar_Results(LR_lv1, MMCLX_data)

    file = date_str + '_profiles_timeseries.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


if plot_comparisons:

    fig, plt = Plot_Comparison(LR_lv1, MMCLX_data)

    file = date_str + '_avg_time_height.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

if plot_interpolation_scatter:

    fig, plt = Plot_Scatter(LR_lv1, MMCLX_data)

    file = date_str + '_interp-avgheight_comp.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


#save_log_data(file[:-5], interp_meth, hmin, hmax, date, time_intervall)

if plot_doppler_spectra:
    # fix parameters for testing stuff
    c = 2
    h0 = 55
    t0 = 50

    i_png = 1
    n_spectra = LR_lv0.n_time * sum(LR_lv0.n_height)

    tstart = time.time()

    doppler_dBZ = np.multiply(np.ma.log10(LR_lv0.VHSpec[c]), 10.0)
    zbound = [doppler_dBZ.min(), doppler_dBZ.max()]

    for c in range(LR_lv0.no_c):
        for h0 in range(LR_lv0.n_height[c]):
            for t0 in range(LR_lv0.n_time):
                mean, threshold, var, nnoise = estimate_noise_hs74(LR_lv0.VHSpec[c][t0, h0, :], navg=LR_lv0.no_av[c])

                #        fig, plt = Plot_Doppler_Spectra(LR_lv0, c, t0, h0, zbound, threshold, mean)
                #
                #
                #        date_str = str(datetime_from_seconds(LR_lv0.t_unix[t0]))
                #        date_str = date_str[11:].replace(':', '_')
                #
                #        file = '/Users/willi/data/MeteoData/LIMRad94/PNG/spectra_' \
                #               + str(round(LR_lv0.height[c][h0], 4)) + '_' + date_str + '.png'
                #
                #        fig.savefig(file, dpi=100, format='png')
                #        plt.close()

                print(f'    File {i_png} of {n_spectra} written.', end='\r')
                i_png += 1

    print('\n' * 2)

    print(f'    Elapsed time for noise floor estimation and plotting = {time.time()-tstart:.3f} sec.')

if pts: print(f'    Total Elapsed Time = {time.clock()-start_time:.3f} sec.\n')
