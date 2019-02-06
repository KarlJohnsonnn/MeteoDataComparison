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
   
   
   
   
    The example call to the routine:    $  python LIMRAD94_to_Cloudnet.py 20180729


####################################################################################################################
'''

import time
import warnings

import modules.NetCDF_Mod2 as nc2
from modules.Parameter_Mod import *


# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '       LIMRAD94 - Covert to CloudNet format      ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)

# gather argument
h_min = 0.0  # in X.XX     (float, unit: km - lower y-axis limit)
h_max = 12.0  # in X.XX     (float, unit: km - upper y-axis limit)

if len(sys.argv) >= 2:
    date = str(sys.argv[1])
else:
    date = '20181201'  # in YYMMDD

time_intervall = '000000-240000'  # in HHMMSS-HHMMSS
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

LR_lv1 = nc2.LIMRAD94(LIMRAD_path+date[:4]+'/', date, time_intervall, [h_min, h_max], 'LV1')

LR_lv1.save(meteo_path)

#
#
if pts: print('')
if pts: print(f'    Total Elapsed Time = {time.time()-start_time:.3f} sec.\n')
