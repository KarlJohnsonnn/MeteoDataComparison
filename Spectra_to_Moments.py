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
   
   
   
   
    The example call to the routine:    $  python Spectra_to_Moments.py 180810 0500 0600 0.0 12.0 2.0
                                                                          |      |    |    |   |   |
                                                                         date   from to  from  to  std div
                                                                                  (UTC)     (km)   (noise est)

    The path to the netcdf files must contain: [...]/YYMMDD/LVx/

####################################################################################################################
'''
start_time = time.clock()

n_std_diviations = -1.0

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

    if len(sys.argv) == 7:
        n_std_diviations = float(sys.argv[6])

else:

    # special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
    h_min = 0.0  # (km)  - lower y-axis limit
    h_max = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
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
calc_doppler_spectra = True
save_spectra_to_png = False
save_noise_comparison = True
save_moment_differences = False
save_moments_without_noise = False

########################################################################################################################
########################################################################################################################
########################################################################################################################


# remove noise from raw spectra and calculate radar moments
if calc_doppler_spectra:

    tstart = time.time()

    for ic in range(LR_lv0.no_c):
        LR_lv0.VHSpec[ic] = np.ma.masked_less_equal(LR_lv0.VHSpec[ic], -999.0)

    # Estimate Noise Floor using Hildebrand & Sekhon Algorithm
    mean_noise, threshold, variance, numnoise, integration_bounds = remove_noise(LR_lv0, n_std_diviations)
    if pts: print('    - all noise removed ')

    include_noise = True
    if include_noise:
        for ic in range(LR_lv0.no_c):
            integration_bounds[ic][:, :, 0] = 0
            integration_bounds[ic][:, :, 1] = -1


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
#        for iR in range(len(LR_lv0.height_all)):
#            print(' difference l0mom - l1mom = {}:{}:{}'.format(LR_lv0.t_plt[iT].hour,
#                                                                LR_lv0.t_plt[iT].minute,
#                                                                LR_lv0.t_plt[iT].second),
#                  '  height = {:.5f} (km)    diffmdv '.format(LR_lv0.height_all[iR]), LR_lv0.diffmdv[iR, iT])

    compare_datasets(LR_lv0, LR_lv1)

########################################################################################################################
########################################################################################################################
########################################################################################################################


if save_spectra_to_png:
    n_png = sum(LR_lv0.n_height) * LR_lv0.Time
    n_png = LR_lv0.Time
    i_png = 0

    ic = 2
    h0 = 57

    bsp_time = string_to_datetime(LR_lv0, '05:00:00')
    bsp_height = 6.0

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
    for t0 in range(LR_lv0.Time):
        itime = t0

    # show mean noise, threshold, and integration lines + spectrum
        fig, plt = Plot_Doppler_Spectra(LR_lv0, ichirp, t0, iheight, [-40, 10],
                                        threshold[ichirp][t0, iheight],
                                        mean_noise[ichirp][t0, iheight],
                                        integration_bounds[ichirp][t0, iheight, :])

    # show only spectra
        #    fig, plt = Plot_Doppler_Spectra(LR_lv0, ichirp, itime, iheight, [-60, 20])

        datestring = str(LR_lv0.t_plt[itime])
        idxSpace = str(datestring).find(' ')
        file = '/Users/willi/data/MeteoData/LIMRad94/PNG/' + date + '_' \
               + str(datestring[idxSpace + 1:]) + '_' + '{:.5f}'.format(LR_lv0.height_all[iheight]) \
               + '_spectra_' + str(i_png).zfill(3) + '.png'

        fig.savefig(file, dpi=100, format='png')
        plt.close()
        if pts: print("    Save spectra: {} of {} ".format(i_png, n_png), end="\r")
        i_png += 1

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
