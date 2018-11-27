import sys
import time
import warnings

import modules.NetCDF_Mod2 as nc2
from modules.Parameter_Mod import *

# from modules.PlotLibrary_Mod2 import *
# from modules.Utility_Mod import *

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

path_to_data = '/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1/'
path_to_output = '/Users/willi/Desktop/tmp/limrad_to_cloudnet/'

# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '       LIMRAD94 - Covert to CloudNet format      ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)

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

if pts:
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

LR_lv1 = nc2.LIMRAD94(path_to_data, date, time_intervall, [h_min, h_max])

LR_lv1.save(path_to_output)

#
#
if pts: print('')
if pts: print(f'    Total Elapsed Time = {time.time()-start_time:.3f} sec.\n')
