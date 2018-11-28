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

LR_lv0_b = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV0/180808_050001_P01_ZEN.LV0.NC')
LR_lv0_a = nc2.LIMRAD94('/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV0/180820_050000_P01_ZEN.LV0.NC')

constants_to_compare = ['AvgNum', 'NoiseFilt', 'SampDur', 'SampRate', 'MaxVel', 'DoppRes']
variables_to_compare = ['SeqIntTime', 'QualFlag', 'Status', 'TPow']

print('')
print('    COMPARE SAME CHIRP TABLES DIFFERENT DAY\n')

for item in constants_to_compare:
    eq = LR_lv0_a.dimensions[0][item] == LR_lv0_b.dimensions[0][item]
    print('')
    print(' Constants :: ', item, '  is equal ?  ', eq)

    if not eq:
        print('     --> a :: ', LR_lv0_a.dimensions[0][item])
        print('     --> b :: ', LR_lv0_b.dimensions[0][item])

for item in variables_to_compare:
    eq = np.array_equal(LR_lv0_a.time_series_1D[0][item]['Val'],
                        LR_lv0_b.time_series_1D[0][item]['Val'])
    print('')
    print(' Variables :: ', item, '  is equal ?  ', eq)

    if not eq:
        dif = Diff(LR_lv0_a.time_series_1D[0][item]['Val'], LR_lv0_b.time_series_1D[0][item]['Val'])
        print('     --> a :: ', LR_lv0_a.time_series_1D[0][item]['Val'])
        print('     --> b :: ', LR_lv0_b.time_series_1D[0][item]['Val'])
        if dif:
            for ival1, ival2 in zip(LR_lv0_a.time_series_1D[0][item]['Val'], LR_lv0_b.time_series_1D[0][item]['Val']):
                print('  diff values ::  ', ival1, ival2)

print('')

if pts: print('')
# if pts: print('    Save Figure to File :: ' + file + '\n')
if pts: print(f'    Total Elapsed Time = {time.time()-start_time:.3f} sec.\n')
