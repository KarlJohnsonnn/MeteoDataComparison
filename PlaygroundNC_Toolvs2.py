import warnings

import modules.NetCDF_Mod2 as nc2
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
plot_radar_results = False
plot_compare_noise = False
plot_for_poster = False
plot_comparisons = False
plot_interp2d = False
plot_interpolation_scatter = False

interpolate_cn = False
create_nc_file = True

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
if pts: Print_Head()

hmin = 0.0  # (km)  - lower y-axis limit
hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
date = '180729'  # in YYMMDD
time_intervall = '000000-235959'  # in HHMMSS-HHMMSS

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

# LR_lv0 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/180729/LV0/', date, time_intervall, [hmin, hmax])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/180729/LV1/', date, time_intervall, [hmin, hmax])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/noise/180810/LV1/', date, time_intervall, [hmin, hmax])
LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/VdResDiff/180810/LV1/', '180810', time_intervall,
                      [hmin, hmax])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/VdResDiff/180810/LV1/VdRes2cms_180810_055219_P10_ZEN.LV1.NC')


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
# LR_lv1.avg_time()
# MIRA_data.avg_time()
# MMCLX_data.avg_time()

# calculate Ze, mdv, sw (height averaged)
# LR_lv1.avg_height()
# MIRA_data.avg_height()
# MMCLX_data.avg_height()


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

# remove noise from raw spectra and calculate radar moments
if calc_doppler_spectra:

    tstart = time.time()

    output = remove_noise(LR_lv0)
    if pts: print('    - all noise removed ')

    if include_noise:
        for ic in range(LR_lv0.no_c):
            output[4][ic][:, :, 0] = 0
            output[4][ic][:, :, 1] = -1

    output = spectra_to_moments(LR_lv0.VHSpec, LR_lv0.DopplerBins, output[4])

    if pts: print('    - moments calculated ')

    LR_lv0.save(LIMRAD_path, Ze=output[0], mdv=output[1], sw=output[2], skew=output[3], kurt=output[4])

    if pts: print('\n' * 2)
    if pts: print(f'    Elapsed time for noise floor estimation and plotting = {time.time()-tstart:.3f} sec.')

    compare_datasets(LR_lv0, LR_lv1)

if plot_for_poster:

    fig, plt = Plot_for_poster(LR_lv1)

    file = meteo_path + date + '_radar_LIMRAD94_vergleich.png'
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

if plot_compare_noise:
    # print(LR_lv0.Ze.shape, LR_lv1.Ze.shape)
    fig, plt = Plot_Compare_NoiseFac0(LR_lv1, LR_lv0)

    file = date + '_NoiseFac0_Lv1_Lv0moments.png'
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

# save_log_data(file[:-5], interp_meth, hmin, hmax, date, time_intervall)

if pts: print(f'    Total Elapsed Time = {time.clock()-start_time:.3f} sec.\n')
