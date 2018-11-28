import sys
import warnings

import modules.NetCDF_Mod as nc
from modules.PlotLibrary_Mod import *
from modules.Utility_Mod import *

# Logicals for different tasks
calc_doppler_spectra = True
plot_spectra = False
plot_compare_noise = True

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
    print('  \u2503' + '          LIMRAD94 - Spectra to Moments          ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)

# gather arguments
if len(sys.argv) > 6:
    date = str(sys.argv[1])
    time_intervall = str(sys.argv[2]) + '-' + str(sys.argv[3])
    hmin, hmax = float(sys.argv[4]), float(sys.argv[5])

else:

    # special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
    hmin = 0.0  # (km)  - lower y-axis limit
    hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    date = '180810'  # in YYMMDD
    time_intervall = '0500-0600'  # in HHMM-HHMM

    #  special case NoiseFac0_file = 'NOISEFAC0_180820_142451_P01_ZEN.LV0.NC'
    # hmin = 0.0  # (km)  - lower y-axis limit
    # hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    # date = '180820'  # in YYMMDD
    # time_intervall = '1400-1500'  # in HHMM-HHMM

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
print('     date: ', date, time_intervall, hmin, hmax, '\n')
print('     is this the correct folder??')

LR_lv0 = nc.LIMRAD94_LV0(date, time_intervall, [hmin, hmax])
LR_lv1 = nc.LIMRAD94_LV1(date, time_intervall, [hmin, hmax])

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

# remove noise from raw spectra and calculate radar moments
if calc_doppler_spectra:

    tstart = time.time()

    # Estimate Noise Floor using Hildebrand & Sekhon Algorithm
    mean_noise, threshold, variance, numnoise, integration_bounds = remove_noise(LR_lv0)
    if pts: print('    - all noise removed ')

    include_noise = True
    if include_noise:
        for ic in range(LR_lv0.no_c):
            integration_bounds[ic][:, :, 0] = 0
            integration_bounds[ic][:, :, 1] = -1

    if plot_spectra:
        n_png = sum(LR_lv0.n_height) * LR_lv0.Time
        n_png = LR_lv0.Time
        i_png = 0

        ic = 2
        h0 = 57

        # ######  TESTING #######
        # save to .mat file to compare with matlab routine
        # import scipy.io
        # scipy.io.savemat('/Users/willi/data/MeteoData/LIMRad94/test_spectrum.mat',
        #                 {'spectrum': LR_lv0.VHSpec[ic][0, h0, :]})

        # mean, threshold, var, nnoise, left_intersec, right_intersece = \
        #    estimate_noise_hs74(LR_lv0.VHSpec[ic][0, h0, :], navg=64)

        # ######  TESTING #######

        # for ic in range(LR_lv0.no_c):
        for t0 in range(LR_lv0.Time):
            #       for h0 in range(LR_lv0.n_height[ic]):
            # print(f'         Noise Threshold = {threshold[ic][0,h0]:.15f}')
            # print(f'         Noise mean_noise= {mean_noise[ic][0,h0]:.15f}')
            # print('         integration_bounds= {}    {}'.format(integration_bounds[ic][0, h0, 0],
            #                                                     integration_bounds[ic][0, h0, 1]))

            fig, plt = Plot_Doppler_Spectra(LR_lv0, ic, t0, h0, [-40, 10],
                                            threshold[ic][t0, h0],
                                            mean_noise[ic][t0, h0],
                                            integration_bounds[ic][t0, h0, :])

            file = '/Users/willi/data/MeteoData/LIMRad94/PNG3/' + date + '_spectra_' + str(i_png).zfill(3) + '.png'
            fig.savefig(file, dpi=100, format='png')
            plt.close()
            if pts: print("    Save spectra: {} of {} ".format(i_png, n_png), end="\r")
            i_png += 1

    output = spectra_to_moments(LR_lv0.VHSpec, LR_lv0.DopplerBins, integration_bounds, mean_noise)

    if pts: print('    - moments calculated \n')
    if pts: print(f'    Elapsed time for noise floor estimation and plotting = {time.time()-tstart:.3f} sec.')

    LR_lv0.Ze = np.ma.log10(output[0].T) * 10.0
    LR_lv0.mdv = output[1].T
    LR_lv0.sw = output[2].T
    LR_lv0.skew = output[3].T
    LR_lv0.kurt = output[4].T

    Lv1ZELin = np.power(LR_lv1.Ze / 10, 10)

    LR_lv0.diffZe = np.ma.subtract(output[0].T, Lv1ZELin)

    compare_datasets(LR_lv0, LR_lv1)

if plot_compare_noise:
    fig, plt = Plot_Compare_NoiseFac0(LR_lv1, LR_lv0)

    file = date + '_NoiseFac0_Lv1_Lv0moments.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

    fig, plt = Plot_CalcMoments_minus_GivenMoments(LR_lv0)

    file = date + '_NoiseFac0_Lv0moments-Lv1.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()

    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

if pts: print(f'    Total Elapsed Time = {time.clock()-start_time:.3f} sec.\n')
