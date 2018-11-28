import modules.NetCDF_Mod2 as nc2
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

# Print Head
if pts:
    print(' ')
    print('  \u250F' + 49 * '\u2501' + '\u2513')
    print('  \u2503' + '             LIMRAD94 - Compare MDFs             ' + '\u2503')
    print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
    print('\n' * 2)


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

LR_lv1_a = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1/180808_050001_P01_ZEN.LV1.NC')
LR_lv1_b = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1/180820_050000_P01_ZEN.LV1.NC')

constants_to_compare = ['AvgNum', 'NoiseFilt', 'SampDur', 'MaxVel', 'DoppRes']
variables_to_compare = ['SeqIntTime', 'QualFlag', 'Status', 'TPow']

print('')
print('    COMPARE SAME CHIRP TABLES DIFFERENT DAY :: 08.08.18  vs. 20.08.18 \n')

for item in constants_to_compare:
    eq = LR_lv1_a.dimensions[0][item] == LR_lv1_b.dimensions[0][item]
    print('')
    print(' Constants :: ', item, '  is equal ?  ', eq)

    if not eq:
        print('     --> a :: ', LR_lv1_a.dimensions[0][item])
        print('     --> b :: ', LR_lv1_b.dimensions[0][item])

for item in variables_to_compare:
    eq = np.array_equal(LR_lv1_a.time_series_1D[0][item]['Val'],
                        LR_lv1_b.time_series_1D[0][item]['Val'])
    print('')
    print(' Variables :: ', item, '  is equal ?  ', eq)

    if not eq:
        dif = Diff(LR_lv1_a.time_series_1D[0][item]['Val'], LR_lv1_b.time_series_1D[0][item]['Val'])
        print('     --> a :: ', LR_lv1_a.time_series_1D[0][item]['Val'])
        print('     --> b :: ', LR_lv1_b.time_series_1D[0][item]['Val'])
        if dif:
            for ival1, ival2 in zip(LR_lv1_a.time_series_1D[0][item]['Val'], LR_lv1_b.time_series_1D[0][item]['Val']):
                print('     different values of ', item, ' :: ', ival1, ival2)

print('')
