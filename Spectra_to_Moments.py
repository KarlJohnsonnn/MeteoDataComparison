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
if len(sys.argv) == 2:
    n_std_diviations = float(sys.argv[1])

else:
    n_std_diviations = 1.0

# special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
hmin = 0.0  # (km)  - lower y-axis limit
hmax = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
date = '180810'  # in YYMMDD
time_intervall = '0500-0600'  # in HHMM-HHMM


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
print('     date: ', date, time_intervall, hmin, hmax)
print('     standard deviations for moment calc: ', n_std_diviations, '\n')
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

# Logicals for different tasks
calc_doppler_spectra = True
save_spectra_to_png = False
save_noise_comparison = False
save_moment_differences = False
save_moments_without_noise = True

########################################################################################################################
########################################################################################################################
########################################################################################################################


# remove noise from raw spectra and calculate radar moments
if calc_doppler_spectra:

    tstart = time.time()

    # Estimate Noise Floor using Hildebrand & Sekhon Algorithm
    mean_noise, threshold, variance, numnoise, integration_bounds = remove_noise(LR_lv0, n_std_diviations)
    if pts: print('    - all noise removed ')

    include_noise = False
    if include_noise:
        for ic in range(LR_lv0.no_c):
            integration_bounds[ic][:, :, 0] = 0
            integration_bounds[ic][:, :, 1] = -1

    if save_spectra_to_png:
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

            fig, plt = Plot_Doppler_Spectra(LR_lv0, ic, t0, h0, [-40, 10],
                                            threshold[ic][t0, h0],
                                            mean_noise[ic][t0, h0],
                                            integration_bounds[ic][t0, h0, :])

            file = '/Users/willi/data/MeteoData/LIMRad94/PNG3/' + date + '_spectra_' + str(i_png).zfill(3) + '.png'
            fig.savefig(file, dpi=100, format='png')
            plt.close()
            if pts: print("    Save spectra: {} of {} ".format(i_png, n_png), end="\r")
            i_png += 1

    output = spectra_to_moments(LR_lv0.VHSpec, LR_lv0.DopplerBins, integration_bounds, LR_lv0.DoppRes)
    # output = spectra_to_moments(LR_lv0.VHSpec, LR_lv0.DopplerBins, integration_bounds, LR_lv0.DoppRes)

    if pts: print('    - moments calculated \n')
    if pts: print(f'    Elapsed time for noise floor estimation and plotting = {time.time()-tstart:.3f} sec.')

    LR_lv0.ZeLin = output[0].T
    LR_lv0.Ze = np.ma.log10(LR_lv0.ZeLin) * 10.0
    LR_lv0.mdv = output[1].T
    LR_lv0.sw = output[2].T
    LR_lv0.skew = output[3].T
    LR_lv0.kurt = output[4].T

    LR_lv0.diffZe = np.ma.subtract(LR_lv0.ZeLin, LR_lv1.ZeLin)
    LR_lv0.diffmdv = np.ma.subtract(LR_lv0.mdv, LR_lv1.mdv)
    LR_lv0.diffsw = np.ma.subtract(LR_lv0.sw, LR_lv1.sw)

#    for iT in range(LR_lv0.n_time):
#       for iR in range(len(LR_lv0.height_all)):
#            print(' difference l0mom - l1mom = {}:{}:{}'.format(LR_lv0.t_plt[iT].hour,
#                                                                LR_lv0.t_plt[iT].minute,
#                                                                LR_lv0.t_plt[iT].second),
#                  '  height = {:.5f} (km)    diffmdv '.format(LR_lv0.height_all[iR]), LR_lv0.diffmdv[iR, iT])

#    compare_datasets(LR_lv0, LR_lv1)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if save_noise_comparison:
    fig, plt = Plot_Compare_NoiseFac0(LR_lv1, LR_lv0)
    file = date + '_NoiseFac0_Lv1_Lv0moments.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()
    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

########################################################################################################################
########################################################################################################################
########################################################################################################################

if save_moment_differences:
    for mom in ['Ze', 'mdv', 'sw']:
        fig, plt = Plot_CalcMoments_minus_GivenMoments(LR_lv0, mom)

        file = date + '_NoiseFac0_Lv0moments-Lv1_' + mom + '.png'
        fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
        plt.close()

        if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

########################################################################################################################
########################################################################################################################
########################################################################################################################

if save_moments_without_noise:

    LR_lv0.n_std_div = n_std_diviations

    for mom in ['Ze']:
        fig, plt = Plot_moment_from_spectra(LR_lv0, mom)
        file = date + '_NoiseFac0_Lv0_to_moments_' + mom + '_nstddiv_' + str(int(n_std_diviations)).zfill(2) + '.png'
        fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
        plt.close()

        if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')


if pts: print(f'    Total Elapsed Time = {time.clock()-start_time:.3f} sec.\n')
