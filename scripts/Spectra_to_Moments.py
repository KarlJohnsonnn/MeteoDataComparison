########################################################################################################################
# THE FOLLOWING 3 LINES ARE NECESSARY FOR INPUT OF modules/ FOLDER !!!
#
import sys, os

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, '..')))
########################################################################################################################


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

import warnings
import modules.NetCDF_Mod as nc
from modules.PlotLibrary_Mod import *
from modules.Utility_Mod import *

start_time = time.clock()

n_std_diviations = 6.0

# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '          LIMRAD94 - Spectra to Moments          ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)

# gather arguments

h_min = 0.0  # (km)  - lower y-axis limit
h_max = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
if len(sys.argv) >= 4:
    date = str(sys.argv[1])
    time_intervall = str(sys.argv[2]) + '-' + str(sys.argv[3])
    if len(sys.argv) >= 6:
        h_min, h_max = float(sys.argv[4]), float(sys.argv[5])
    if len(sys.argv) == 7:
        n_std_diviations = float(sys.argv[6])

else:
    # special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
    date = '20180810'  # in YYYYMMDD
    time_intervall = '052000-053000'  # in HHMMSS-HHMMSS

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
print('    date: ', date, time_intervall, h_min, h_max)
print('    standard deviations for moment calc: ', n_std_diviations, '\n')

LR_lv0 = nc.LIMRAD94_LV0(LIMRAD_path, date, time_intervall, [h_min, h_max])
LR_lv1 = nc.LIMRAD94_LV1(LIMRAD_path, date, time_intervall, [h_min, h_max])

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
user_input = False
include_noise = False
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

    # fill values needs to be masked for noise removal otherwise wrong results
    for ic in range(LR_lv0.no_c):
        LR_lv0.VHSpec[ic] = np.ma.masked_less_equal(LR_lv0.VHSpec[ic], -999.0)

    # Estimate Noise Floor for all chirps, timesteps and heighsteps aka. for all pixels
    # Algorithm used: Hildebrand & Sekhon
    mean_noise, threshold, variance, numnoise, integration_bounds = remove_noise(LR_lv0, n_std_diviations)
    if pts: print('    - all noise removed ')

    if include_noise:
        for ic in range(LR_lv0.no_c):
            integration_bounds[ic][:, :, 0] = 0
            integration_bounds[ic][:, :, 1] = -1

    #  dimensions:
    #       -   LR_lv0.VHSpec       [Nchirps][Ntime,Nheight]
    #       -   LR_lv0.DopplerBins  [Nchirps]
    #       -   integration_bounds  [2]
    #       -   DoppRes             [Nchirps]
    output = spectra_to_moments(LR_lv0.VHSpec, LR_lv0.DopplerBins, integration_bounds, LR_lv0.DoppRes)

    if pts:
        print('    - moments calculated \n')
        print(f'    Elapsed time for noise floor estimation and plotting = {time.time() - tstart:.3f} sec.')

    LR_lv0.ZeLin = output[0].T
    LR_lv0.Ze = np.ma.log10(LR_lv0.ZeLin) * 10.0
    LR_lv0.mdv = output[1].T
    LR_lv0.sw = output[2].T
    LR_lv0.skew = output[3].T
    LR_lv0.kurt = output[4].T

    LR_lv0.diffZe = np.ma.subtract(LR_lv0.ZeLin, LR_lv1.ZeLin)
    LR_lv0.diffmdv = np.ma.subtract(LR_lv0.mdv, LR_lv1.mdv)
    LR_lv0.diffsw = np.ma.subtract(LR_lv0.sw, LR_lv1.sw)

    # print mean differences of Ze, mdv, sw
    compare_datasets(LR_lv0, LR_lv1)

    # this is just for debugging

#    for iT in range(LR_lv0.n_time):
#        for iR in range(len(LR_lv0.height_all)):
#            print(' difference l0mom - l1mom = {}:{}:{}'.format(LR_lv0.t_plt[iT].hour,
#                                                                LR_lv0.t_plt[iT].minute,
#                                                                LR_lv0.t_plt[iT].second),
#                  '  height = {:.5f} (km)    diffmdv '.format(LR_lv0.height_all[iR]), LR_lv0.diffmdv[iR, iT])


########################################################################################################################
########################################################################################################################
########################################################################################################################


if save_spectra_to_png:
    n_png = LR_lv0.Time
    i_png = 0

    if user_input:
        ichirp, itime, iheight = gather_user_input(LR_lv0)
    else:
        ichirp = 2
        iheight = 58
        itime = 50

    # for ic in range(LR_lv0.no_c):
    for t0 in range(LR_lv0.Time):
        itime = t0

        # show mean noise, threshold, and integration lines + spectrum
        fig, plt = Plot_Doppler_Spectra(LR_lv0, ichirp, itime, iheight, [-40, 10],
                                        thresh=threshold[ichirp][itime, iheight],
                                        mean=mean_noise[ichirp][itime, iheight],
                                        int_a=integration_bounds[ichirp][itime, iheight, 0],
                                        int_b=integration_bounds[ichirp][itime, iheight, 1])

        # show only spectra
        #    fig, plt = Plot_Doppler_Spectra(LR_lv0, ichirp, itime, iheight, [-60, 20])

        datestring = str(LR_lv0.t_plt[itime])
        idxSpace = str(datestring).find(' ')
        file = LIMRAD_path + 'PNG2/' + date + '_' \
               + str(datestring[idxSpace + 1:]) + '_' + '{:.5f}'.format(LR_lv0.height[ichirp][iheight]) \
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
    file = date + '_NoiseFac0_Lv1_Lv0moments' + '__nstddiv_' + str(int(n_std_diviations)).zfill(2) + '.png'
    fig.savefig(meteo_path + file, dpi=dpi_val, format='png')
    plt.close()
    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')

########################################################################################################################
########################################################################################################################
########################################################################################################################

if save_moment_differences:
    for mom in ['Ze', 'mdv', 'sw']:
        fig, plt = Plot_CalcMoments_minus_GivenMoments(LR_lv0, mom)

        file = date + '_NoiseFac0_Lv0moments-Lv1_' + mom + '_nstddiv_' + str(int(n_std_diviations)).zfill(2) + '.png'
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

if pts: print(f'    Total Elapsed Time = {time.clock() - start_time:.3f} sec.\n')
