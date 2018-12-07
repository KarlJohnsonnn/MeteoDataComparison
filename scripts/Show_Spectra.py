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



    The example call to the routine: user$  python Show_Spectra.py 20180810 050000 060000 0.0 12.0
                                                                      |       |      |     |    |
                                                                     date    from    to   from  to
                                                                                (UTC)        (km)

    The user has to specify the path to the LIMRAD94 LV0 and LV1 files under modules/Parameter_Mod.py !

####################################################################################################################
'''

# other imports
import warnings
from modules.NetCDF_Mod import *
from modules.PlotLibrary_Mod import *
from modules.Utility_Mod import *
from scipy import signal


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
    date = '20180810'  # in YYMMDD
    time_intervall = '050000-060000'  # in HHMM-HHMM

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

start_time = time.clock()

# ----- LIMRAD 94GHz Radar data extraction
print('    date: ', date, time_intervall, h_min, h_max)

LR_lv0 = LIMRAD94_LV0(LIMRAD_path, date, time_intervall, [h_min, h_max])
LR_lv1 = LIMRAD94_LV1(LIMRAD_path, date, time_intervall, [h_min, h_max])

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

calc_doppler_spectra = False
save_spectra_to_png = True
save_noise_comparison = False
save_moment_differences = False
save_moments_without_noise = False

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

    LR_lv0.VHSpec_nrm = []

    # for ic in range(LR_lv0.no_c):
    for itime in range(LR_lv0.Time):

        VHmin = LR_lv0.VHSpec[ichirp][itime, iheight].min()
        VHmax = LR_lv0.VHSpec[ichirp][itime, iheight].max()
        VHSpec_nrm = (LR_lv0.VHSpec[ichirp][itime, iheight]-VHmin)/(VHmax - VHmin)

        n_scales = 30

        widths = np.linspace(1, 7, n_scales)
        cwtmatr = signal.cwt(VHSpec_nrm, signal.ricker, widths)

        # show spectra, normalized spectra and wavlet transformation
        fig, plt = Plot_Doppler_Spectra_Wavelet_Transform(LR_lv0, VHSpec_nrm, ichirp, itime, iheight,
                                                          [-0.05, 1.05], cwtmatr, widths)

        # show mean noise, threshold, and integration lines + spectrum
        #fig, plt = Plot_Doppler_Spectra(LR_lv0, ic, t0, h0, [-40, 10],
        #                               threshold[ic][t0, h0],
        #                               mean_noise[ic][t0, h0],
        #                               integration_bounds[ic][t0, h0, :])

        # show only spectra
        #fig, plt = Plot_Doppler_Spectra(LR_lv0, ichirp, itime, iheight, [-60, 20])

        datestring = str(LR_lv0.t_plt[itime])
        idxSpace = str(datestring).find(' ')
        file = meteo_path + date + '_' \
               + str(datestring[idxSpace + 1:]) + '_' + '{:.5f}'.format(LR_lv0.height[ichirp][iheight]) \
               + '_spectra_' + str(i_png).zfill(3) + '.png'

        fig.savefig(file, dpi=100, format='png')
        if pts: print('    Save Figure to File :: ' + file + '   {} of {} ".format(i_png, n_png), end="\r")')
        i_png += 1

########################################################################################################################
########################################################################################################################
########################################################################################################################


if pts: print(f'    Total Elapsed Time = {time.clock() - start_time:.3f} sec.\n')
