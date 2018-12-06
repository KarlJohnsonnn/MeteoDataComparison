########################################################################################################################
# THE FOLLOWING 3 LINES ARE NECESSARY FOR INPUT OF modules/ FOLDER !!!
#
import sys, os
#SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
#sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, '..')))
########################################################################################################################

import warnings

import modules.NetCDF_Mod as nc
import modules.NetCDF_Mod2 as nc2

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
n_std_diviations = 2.0

# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '       compare   LIMRAD94 - MIRA    spectra      ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)

# gather arguments

if len(sys.argv) >= 6:
    date = str(sys.argv[1])
    time_intervall = str(sys.argv[2]) + '-' + str(sys.argv[3])
    height, time = float(sys.argv[4]), float(sys.argv[5])

    if len(sys.argv) == 7:
        n_std_diviations = float(sys.argv[6])

else:

    # special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
    height = 2  # (km)  - height of the spectrum to compare
    date = '181203'  # in YYMMDD
    time_intervall = '0100-0200'  # in HHMM-HHMM
    time ='0020' # time of the spectrum to compare
    spectra_height = [2.273, 1.648, 1.416]


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
print('     date: ', date, time_intervall, height)
print('     standard deviations for moment calc: ', n_std_diviations, '\n')
print('     is this the correct folder??')

LR_lv0 = nc.LIMRAD94_LV0(date, time_intervall, [height-0.1, height+0.1])
LR_lv1 = nc.LIMRAD94_LV1(date, time_intervall, [height-0.1, height+0.1])

spectra_time = [string_to_datetime(LR_lv0, '01:58:09'), string_to_datetime(LR_lv0, '01:30:56'),
                string_to_datetime(LR_lv0, '01:33:30')]

if pts: print('importing MIRA file...\n')

# ----- MIRA 35GHz Radar data extraction

MIRA_lv0 = nc2.MIRA35_spectra('/home/tvogl/PhD/comparison_limrad_mira/MIRA/spectra/D20181203_T0000_0230_Pun_zspc2nc_v1_02_standard.nc4')

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
save_spectra_to_png = True
save_noise_comparison = False
save_moment_differences = False
save_moments_without_noise = False


if save_spectra_to_png:
    n_png = LR_lv0.Time
    i_png = 0

    h0 = 57
    for bsp_height in spectra_height:
        bsp_time = spectra_time[i_png]
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
        mira_height0 = min(MIRA_lv0.variables['range'], key=lambda x: abs(x - (bsp_height * 1000)))
        mira_time0 = min(MIRA_lv0.variables['t_plt'], key=lambda  x: abs(x - bsp_time))
        mitime = MIRA_lv0.variables['t_plt'].index(mira_time0)
        miheight = list(MIRA_lv0.variables['range']).index(mira_height0)

        fig, plt, ax = Plot_Doppler_Spectra_LIMRad_MIRA(LR_lv0, ichirp, itime, iheight, [-60, 20], MIRA_lv0, mitime, miheight)
        datestring = str(LR_lv0.t_plt[itime])
        idxSpace = str(datestring).find(' ')
        file = '/home/tvogl/PhD/comparison_limrad_mira/' + date + '_' \
              + str(datestring[idxSpace + 1:]) + '_' + '{:.5f}'.format(LR_lv0.height_all[iheight]) \
              + 'LIMRad_MIRA_spectra_' + '.png'

        fig.savefig(file, dpi=100, format='png')
        plt.close()
        if pts:
           print("    Save spectra: {} ".format(i_png+1), end="\r")
        i_png += 1
