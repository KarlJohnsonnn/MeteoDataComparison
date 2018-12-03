import sys
import warnings

import modules.NetCDF_Mod2 as nc2
from modules.PlotLibrary_Mod2 import *
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
create_nc_file = False

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

# Print Head
if pts: Print_Head()

# gather arguments
if len(sys.argv) == 6:
    date = str(sys.argv[1])
    time_intervall = str(sys.argv[2]) + '-' + str(sys.argv[3])
    h_min, h_max = float(sys.argv[4]), float(sys.argv[5])

else:

    # h_min = 0.0  # (km)  - lower y-axis limit
    # h_max = 12.0  # (km) - upper y-axis limit, highest range gate may be higher
    # date = '180729'  # in YYMMDD
    # time_intervall = '000000-240000'  # in HHMMSS-HHMMSS

    # special case NoiseFac0_file = 'NoiseFac0/NoiseFac0_180810_052012_P01_ZEN.LV0.NC'
    h_min = 0.0  # (km)  - lower y-axis limit
    h_max = 12.00  # (km) - upper y-axis limit, highest range gate may be higher
    date = '181202'  # in YYMMDD
    time_intervall = '000000-240000'  # in HHMM-HHMM


warnings.filterwarnings("ignore")

start_time = time.time()

print('    Input values:')
print('         - date =       ', date)
print('         - time from:   ', time_intervall, ' (UTC)')
print('         - height from: ', h_min, '(km)  to: ', h_max, ' (km) \n')


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

# LR_lv0 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/noise/180810/LV0/', date, time_intervall, [h_min, h_max])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/180729/LV1/', date, time_intervall, [h_min, h_max])

#LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1/', date, time_intervall, [h_min, h_max])

# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1/', date, time_intervall, [h_min, h_max])
LR_lv1 = nc2.LIMRAD94('/Volumes/Data_Storag/MeteoData/LIMRad94/LV1/', date, time_intervall, [h_min, h_max])

# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/noise/180810/LV1/', date, time_intervall, [h_min, h_max])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/VdResDiff/180810/LV1/', '180810', time_intervall, [h_min, h_max])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/VdResDiff/180810/LV1/VdRes2cms_180810_055219_P10_ZEN.LV1.NC')


#LR_lv1.save('/Users/willi/Desktop/tmp/limrad_to_cloudnet/')

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

# fig, plt = Plot_Time_Series(LR_lv1, ['ZE'])
fig, plt = Plot_Time_Series(LR_lv1, ['ZE', 'MeanVel', 'SpecWidth', 'SLDR'])
#
file = meteo_path + date + '_' + time_intervall + '_time_series_LIMRAD94.png'
fig.savefig(file, bbox_inches='tight', dpi=dpi_val, format='png')
#plt.close()
#
if pts: print('')
if pts: print('    Save Figure to File :: ' + file + '\n')
if pts: print(f'    Total Elapsed Time = {time.time()-start_time:.3f} sec.\n')
