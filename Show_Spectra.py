import sys
import warnings

import modules.NetCDF_Mod as nc
from modules.PlotLibrary_Mod import *
from modules.Utility_Mod import *

'''
####################################################################################################################

   ##     ##    ###    #### ##    ##         ########  ########   #######   ######   ########     ###    ##     ##
   ###   ###   ## ##    ##  ###   ##         ##     ## ##     ## ##     ## ##    ##  ##     ##   ## ##   ###   ###
   #### ####  ##   ##   ##  ####  ##         ##     ## ##     ## ##     ## ##        ##     ##  ##   ##  #### ####
   ## ### ## ##     ##  ##  ## ## ## ####### ########  ########  ##     ## ##   #### ########  ##     ## ## ### ##
   ##     ## #########  ##  ##  ####         ##        ##   ##   ##     ## ##    ##  ##   ##   ######### ##     ##
   ##     ## ##     ##  ##  ##   ###         ##        ##    ##  ##     ## ##    ##  ##    ##  ##     ## ##     ##
   ##     ## ##     ## #### ##    ##         ##        ##     ##  #######   ######   ##     ## ##     ## ##     ##




    The example call to the routine:    $  python Show_Spectra.py 180810 0500 0600 0.0 12.0
                                                                    |      |    |    |   |
                                                                   date   from to  from  to
                                                                            (UTC)     (km)

    The path to the netcdf files must contain: [...]/YYMMDD/LVx/

####################################################################################################################
'''
start_time = time.clock()


# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '          LIMRAD94 - Spectra to Moments          ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)

# gather arguments

if len(sys.argv) >= 6:
    date = str(sys.argv[1])
    time_intervall = str(sys.argv[2]) + '-' + str(sys.argv[3])
    h_min, h_max = float(sys.argv[4]), float(sys.argv[5])

else:

    # special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
    h_min = 0.0  # (km)  - lower y-axis limit
    h_max = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    date = '181203'  # in YYMMDD
    time_intervall = '0000-0200'  # in HHMM-HHMM

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
print('     date: ', date, time_intervall, h_min, h_max)
print('     standard deviations for moment calc: ', n_std_diviations, '\n')
print('     is this the correct folder??')

LR_lv0 = nc.LIMRAD94_LV0(date, time_intervall, [h_min, h_max])
LR_lv1 = nc.LIMRAD94_LV1(date, time_intervall, [h_min, h_max])

if pts: print('')

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

# Logicals for different tasks
calc_doppler_spectra = False
save_spectra_to_png = True
save_noise_comparison = False
save_moment_differences = False
save_moments_without_noise = False

########################################################################################################################
########################################################################################################################
########################################################################################################################


if save_spectra_to_png:
    n_png = sum(LR_lv0.n_height) * LR_lv0.Time
    n_png = LR_lv0.Time
    i_png = 0

    ic = 2
    h0 = 57

    str_time = input('     Enter time in (HH:MM:SS) ::  ')
    str_height = input('     Enter height in (km)     ::  ')

    bsp_time = string_to_datetime(LR_lv0, str_time)
    bsp_height = float(str_height)

    bsp_height0 = min(LR_lv0.height_all, key=lambda x: abs(x - bsp_height))
    bsp_time0 = min(LR_lv0.t_plt, key=lambda x: abs(x - bsp_time))

    itime = LR_lv0.t_plt.index(bsp_time0)

    for ic in range(LR_lv0.no_c):
        try:
            idx_height = list(LR_lv0.height[ic]).index(bsp_height0)
            if idx_height > 0:
                ichirp = ic
                iheight = idx_height

                break
        except:
            dummy = 0

    # for ic in range(LR_lv0.no_c):
    # for t0 in range(LR_lv0.Time):
    # itime = t0

    # show mean noise, threshold, and integration lines + spectrum
    #            fig, plt = Plot_Doppler_Spectra(LR_lv0, ic, t0, h0, [-40, 10],
    #                                            threshold[ic][t0, h0],
    #                                            mean_noise[ic][t0, h0],
    #                                            integration_bounds[ic][t0, h0, :])

    # show only spectra
    fig, plt = Plot_Doppler_Spectra(LR_lv0, ichirp, itime, iheight, [-60, 20])

    datestring = str(LR_lv0.t_plt[itime])
    idxSpace = str(datestring).find(' ')
    file = '/Users/willi/data/MeteoData/LIMRad94/PNG/' + date + '_' \
           + str(datestring[idxSpace + 1:]) + '_' + '{:.5f}'.format(bsp_height0) \
           + '_spectra_' + str(i_png).zfill(3) + '.png'

    fig.savefig(file, dpi=100, format='png')
    plt.close()
    if pts: print('    Save Figure to File :: ' + file + '\n')
    # if pts: print("    Save spectra: {} of {} ".format(i_png, n_png), end="\r")
    i_png += 1

########################################################################################################################
########################################################################################################################
########################################################################################################################


if pts: print(f'    Total Elapsed Time = {time.clock() - start_time:.3f} sec.\n')
