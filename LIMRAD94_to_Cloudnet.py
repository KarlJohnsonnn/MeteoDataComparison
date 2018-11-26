import sys
import warnings

import modules.NetCDF_Mod2 as nc2
from modules.PlotLibrary_Mod2 import *
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
   
   
   
   
    The example call to the routine:    $  python LIMRAD94_to_Cloudnet.py 180729 000000 240000 0.0 12.0 



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

    date = '180729'  # in YYMMDD
    time_intervall = '000000-240000'  # in HHMMSS-HHMMSS
    h_min = 0.0  # in X.XX     (float, unit: km - lower y-axis limit)
    h_max = 12.0  # in X.XX     (float, unit: km - upper y-axis limit)

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

LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1/', date, time_intervall, [h_min, h_max])

# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/noise/180810/LV1/', date, time_intervall, [h_min, h_max])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/VdResDiff/180810/LV1/', '180810', time_intervall, [h_min, h_max])
# LR_lv1 = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/VdResDiff/180810/LV1/VdRes2cms_180810_055219_P10_ZEN.LV1.NC')


LR_lv1.save('/Users/willi/Desktop/tmp/limrad_to_cloudnet/')

#
#
if pts: print('')
if pts: print(f'    Total Elapsed Time = {time.time()-start_time:.3f} sec.\n')
