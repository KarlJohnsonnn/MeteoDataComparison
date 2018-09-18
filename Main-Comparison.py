from Functions_Mod import *
from Graphics_Mod import *
from IO_Mod       import *
from Parameter_Mod import *
from Interpolation_Tool import *
import NetCDF_Tool  as nc

rc('font', family='serif')
rc('text', usetex=True)

import datetime

import sys, warnings, time



##################################################################################################
#
#       ##     ##  ######  ######## ########     #### ##    ## ########  ##     ## ########
#       ##     ## ##    ## ##       ##     ##     ##  ###   ## ##     ## ##     ##    ##
#       ##     ## ##       ##       ##     ##     ##  ####  ## ##     ## ##     ##    ##
#       ##     ##  ######  ######   ########      ##  ## ## ## ########  ##     ##    ##
#       ##     ##       ## ##       ##   ##       ##  ##  #### ##        ##     ##    ##
#       ##     ## ##    ## ##       ##    ##      ##  ##   ### ##        ##     ##    ##
#        #######   ######  ######## ##     ##    #### ##    ## ##         #######     ##
#
##################################################################################################

# Logicals for different tasks
plot_interp2d = False
plot_RectBivariateSpline = False  # interp2 testen
plot_radar_results = True
plot_comparisons = True
plot_interpolation_scatter = True
plot_compare_mira_mmclx = False


os.chdir(meteo_path)


####################################################################################################################
#
#   ##     ##    ###    #### ##    ##         ########  ########   #######   ######   ########     ###    ##     ##
#   ###   ###   ## ##    ##  ###   ##         ##     ## ##     ## ##     ## ##    ##  ##     ##   ## ##   ###   ###
#   #### ####  ##   ##   ##  ####  ##         ##     ## ##     ## ##     ## ##        ##     ##  ##   ##  #### ####
#   ## ### ## ##     ##  ##  ## ## ## ####### ########  ########  ##     ## ##   #### ########  ##     ## ## ### ##
#   ##     ## #########  ##  ##  ####         ##        ##   ##   ##     ## ##    ##  ##   ##   ######### ##     ##
#   ##     ## ##     ##  ##  ##   ###         ##        ##    ##  ##     ## ##    ##  ##    ##  ##     ## ##     ##
#   ##     ## ##     ## #### ##    ##         ##        ##     ##  #######   ######   ##     ## ##     ## ##     ##
#
####################################################################################################################

start_time = time.clock()

# Print Head
print(' ')
print('  \u250F' + 49 * '\u2501' + '\u2513')
print('  \u2503' + '      LIMRAD94 - MIRA35  Comparison       ' + '\u2503')
print('  \u2517' + 49 * '\u2501' + '\u251B' + '\n')
print('')
print('')

# gather arguments
if len(sys.argv) == 6:

    hmin, hmax = float(sys.argv[1]), float(sys.argv[2])
    date = str(sys.argv[3])
    time_intervall = str(sys.argv[4]) + '-' + str(sys.argv[5])

else:

    ## cirrus
    #hmin = 8.50  # (km)  - lower y-axis limit
    #hmax = 10.0  # (km) - upper y-axis limit, highest range gate may be higher
    #date = '180728'  # in YYMMDD
    #time_intervall = '0740-0805'  # in HHMM-HHMM

    ##cummulis
    #hmin = 5.0  #(km)  - lower y-axis limit
    #hmax = 12.0  #(km) - upper y-axis limit, highest range gate may be higher
    #date     = '180802'     # in YYMMDD
    #time_intervall = '0330-1200'  # in HHMM-HHM

    hmin = 0.0 #(km)  - lower y-axis limit
    hmax = 12.00 #(km) - upper y-axis limit, highest range gate may be higher
    date     = '180729'     # in YYMMDD
    time_intervall = '2200-2359'  # in HHMM-HHMMM

    ##nimbus
    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180805'     # in YYMMDD
    # time_intervall = '0510-0620'  # in HHMM-HHMM

    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180805'     # in YYMMDD
    # time_intervall = '1030-1200'  # in HHMM-HHMM

    # hmin = 1.0 #(km)  - lower y-axis limit
    # hmax = 2.5 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180805'     # in YYMMDD
    # time_intervall = '0700-1210'  # in HHMM-HHMM

    ##nimbus
    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 12.00 #(km) - upper y-axis limit, highest range gate may be higher
    # date     = '180808'     # in YYMMDD
    # time_intervall = '1330-1700'  # in HHMM-HHMM

warnings.filterwarnings("ignore")

# calculate the time in decimal hours
comp_hours = [int(time_intervall[0:2]), int(time_intervall[5:7])]
comp_minutes = [int(time_intervall[2:4]), int(time_intervall[7:9])]

clock_time = np.array(comp_hours) + np.divide(comp_minutes, 60.)  # [hours] + [minutes]/60#

# -- gathering year, month, day for convertion to UTC time
plotyear = int('20' + date[:2])
plotmonth = int(date[2:4])
plotday = int(date[4:6])

time_int = [0, 0, 0, 0]
time_int[0] = datetime.datetime(plotyear, plotmonth, plotday,
                                hour=int(comp_hours[0]), minute=int(comp_minutes[0]))
time_int[1] = time_int[0] + datetime.timedelta(seconds=15)
time_int[3] = datetime.datetime(plotyear, plotmonth, plotday,
                                hour=int(comp_hours[1]), minute=int(comp_minutes[1]))
time_int[2] = time_int[3] - datetime.timedelta(seconds=15)

height = [hmin, hmax]

######################################################################################################
#
#   ########     ###    ########    ###                  #### ##    ## ########  ##     ## ########
#   ##     ##   ## ##      ##      ## ##                  ##  ###   ## ##     ## ##     ##    ##
#   ##     ##  ##   ##     ##     ##   ##                 ##  ####  ## ##     ## ##     ##    ##
#   ##     ## ##     ##    ##    ##     ##    #######     ##  ## ## ## ########  ##     ##    ##
#   ##     ## #########    ##    #########                ##  ##  #### ##        ##     ##    ##
#   ##     ## ##     ##    ##    ##     ##                ##  ##   ### ##        ##     ##    ##
#   ########  ##     ##    ##    ##     ##               #### ##    ## ##         #######     ##
#
######################################################################################################


# ----- LIMRad 94GHz Radar data extraction
#time1 = time.clock()
#LR_time, UTC_time_LR, LR_height, \
#LR_Ze, LR_mdv, LR_sw = extract_dataset(date, time_int, clock_time, height, '*.LV1.NC', '')
#
#
## ----- MIRA 35GHz Radar data extraction
#
## mira.nc data file (processed data)
#mira_timeNC, UTC_time_miraNC, mira_heightNC, \
#miraNC_Z, miraNC_VEL, miraNC_RMS = extract_dataset(date, time_int, clock_time, height, '*mira.nc', '')
#
#
## .mmclx data file (hydrometeors only)
#_, _, _, mira_Ze, mira_VEL, mira_RMS = extract_dataset(date, time_int, clock_time, height, '*.mmclx', '')
#
## .mmclx data file (all targets)
#mira_time, UTC_time_mira, mira_height, \
#mira_Zg, mira_VELg, mira_RMSg = extract_dataset(date, time_int, clock_time, height, '*.mmclx', 'g')
#
#print(' time with funcitons = ', time.clock()-time1)


time1 = time.clock()

LR_data  = nc.LIMRad94_LV1(date, time_intervall, [hmin, hmax])
MIRA_data  = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax])
MMCLX_data = nc.MIRA35_LV1(date, time_intervall, [hmin, hmax], '*.mmclx')

print(' time with classes = ', time.clock()-time1)

print('')

interp_meth = 'nearest'

####################################################################################################################
#
#        ######  ########    ###    ######## ####  ######  ######## ####  ######   ######
#       ##    ##    ##      ## ##      ##     ##  ##    ##    ##     ##  ##    ## ##    ##
#       ##          ##     ##   ##     ##     ##  ##          ##     ##  ##       ##
#        ######     ##    ##     ##    ##     ##   ######     ##     ##  ##        ######
#             ##    ##    #########    ##     ##        ##    ##     ##  ##             ##
#       ##    ##    ##    ##     ##    ##     ##  ##    ##    ##     ##  ##    ## ##    ##
#        ######     ##    ##     ##    ##    ####  ######     ##    ####  ######   ######
#
####################################################################################################################

# time averaged values
    # calculate Ze, mdv, sw (time averaged)
LR_data.avg_time()
MIRA_data.avg_time()

    # calculate Ze, mdv, sw (height averaged)
LR_data.avg_height()
MIRA_data.avg_height()


####################################################################################################################
#
#       ######   ########     ###    ########  ##     ## ####  ######   ######
#      ##    ##  ##     ##   ## ##   ##     ## ##     ##  ##  ##    ## ##    ##
#      ##        ##     ##  ##   ##  ##     ## ##     ##  ##  ##       ##
#      ##   #### ########  ##     ## ########  #########  ##  ##        ######
#      ##    ##  ##   ##   ######### ##        ##     ##  ##  ##             ##
#      ##    ##  ##    ##  ##     ## ##        ##     ##  ##  ##    ## ##    ##
#       ######   ##     ## ##     ## ##        ##     ## ####  ######   ######
#
####################################################################################################################

if plot_interp2d:

    npts = 500
    stat_pos = [0.2, -0.4]


    LRtoMIRA_Ze = Interpolate_2D_neu(LR_data, LR_data.Ze, MIRA_data, interp_meth)
    MIRAtoLR_Ze = Interpolate_2D_neu(MIRA_data, MIRA_data.Ze, LR_data, interp_meth)

    # converting back to linear untis for mean difference calculation
    LR_Ze_mm6m3_I = np.power(np.divide(LRtoMIRA_Ze, 10.0), 10.0)
    LR_Ze_mm6m3 = np.power(np.divide(LR_data.Ze, 10.0), 10.0)
    mira_Zg_mm6m3_I = np.power(np.divide(MIRAtoLR_Ze, 10.0), 10.0)
    mira_Zg_mm6m3 = np.power(np.divide(MIRA_data.Ze, 10.0), 10.0)

    # mean_diff_LRtoM = np.mean(LRtoMIRA_Ze - mira_Zg)
    # mean_diff_MtoLR = np.mean(LR_Ze - MIRAtoLR_Ze)

    # display mean difference in dBZ again
    mean_diff_LRtoM = 10 * np.log10(np.mean(mira_Zg_mm6m3 - LR_Ze_mm6m3_I))
    mean_diff_MtoLR = 10 * np.log10(np.mean(mira_Zg_mm6m3_I - LR_Ze_mm6m3))

    # correlation is calculated with logarithmic units
    correlation_LR_MIRA = correlation(LRtoMIRA_Ze, MIRA_data.Ze)
    correlation_MIRA_LR = correlation(LR_data.Ze, MIRAtoLR_Ze)

    fig, ((p1, p2), (p3, p4), (p5, p6)) = plt.subplots(nrows=3, ncols=2, figsize=(12, 8))

    ########################################################################################################
    # plot table
    cl = ['', r'time resolution $\delta t$', r'height resolution $\delta h$',
          r'height in [km]', r'mean difference in [dBZ]', r'correlation ']

    str_meanLRM = '{:5.2f}'.format(np.mean(mean_diff_LRtoM))
    str_meanMLR = '{:5.2f}'.format(np.mean(mean_diff_MtoLR))
    str_corrLRM = '{:5.2f}'.format(np.mean(correlation_LR_MIRA))
    str_corrMLR = '{:5.2f}'.format(np.mean(correlation_MIRA_LR))

    table = r'\renewcommand{\arraystretch}{1.25}' \
            r'\begin{tabular}{ | c | c | c | c | c | c |} \hline' \
            r'       & ' + str(cl[1]) + r' & ' + str(cl[2]) + r' & ' + str(cl[3]) + r' & ' + str(cl[4]) + r'  & ' + str(
        cl[5]) + r'  \\ \hline ' \
                 r'MIRA   & $\sim 3$ sec & $30$ m & $' + str(hmin) + r'-' + str(hmax) + \
            r' $ & $\overline{\mathcal{I}(Z_e^{35})-Z_e^{94}}=' + str_meanMLR + \
            r' $ & $\rho(\mathcal{I}(Z_e^{35}),Z_e^{94})=$' + \
            r' $' + str_corrMLR + r'$\\ \hline ' \
                                  r'LIMRad & $\sim 5$ sec & $30$ m & $' + str(hmin) + r'-' + str(hmax) + \
            r' $ & $\overline{Z_e^{35}-\mathcal{I}(Z_e^{94})}=' + str_meanLRM + \
            r' $ & $\rho(Z_e^{35},\mathcal{I}(Z_e^{94}))= $ ' \
            r' $ ' + str_corrLRM + r'$\\ \hline ' \
                                   r'\end{tabular} '

    fig.text(0.1, 0.8, table, size=12)

    ########################################################################################################
    # LR_Zelectivity plot
    if pts: print('')
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    p1.set_title(r'Radar Reflectivity Factor $Z_{e}^{94}$')
    plot_data_set(fig, p1, '',
                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
                  x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    p2.set_title(r'Radar Reflectivity Factor $Z_{g}^{35}$')
    plot_data_set(fig, p2, '',
                  UTC_time_mira, mira_height, mira_Zg, vmi=-50, vma=20,
                  x_min=UTC_time_mira[0], x_max=UTC_time_mira[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    if pts: print('       -   synchronized Radar Reflectivity Factor   ', end='', flush=True)
    # MIRA reflectivit

    p3.set_title(r'Interpolation of MIRA $Z_{g}^{35}$'
                 ' onto LIMRad grid resolution',
                 multialignment='center')

    plot_data_set(fig, p3, '',
                  UTC_time_LR, LR_height, MIRAtoLR_Ze, vmi=-50, vma=20,
                  x_min=UTC_time_LR[1], x_max=UTC_time_LR[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='dBZ')

    p4.set_title(r'Interpolation of LIMRad $Z_{e}^{94}$'
                 ' onto MIRA grid resolution', multialignment='center')

    plot_data_set(fig, p4, '',
                  UTC_time_mira, mira_height, LRtoMIRA_Ze, vmi=-50, vma=20,
                  x_min=UTC_time_mira[1], x_max=UTC_time_mira[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    if pts: print('       -   correnlation of Radar Reflectivity Factors   ', end='', flush=True)
    # MIRA reflectivit

    p5.set_title(r'Correlation of interp$(Z_{g}^{35})$ and $Z_e^{94}$', multialignment='center')
    plot_correlation(p5, '', UTC_time_LR, correlation_MIRA_LR, '', '.',
                     x_min=UTC_time_mira[1], x_max=UTC_time_mira[-1],
                     y_min=-1.1, y_max=1.1, x_lab='Time (UTC)', y_lab='correlation')

    p6.set_title(r'Correlation of interp$(Z_{e}^{94})$ and $Z_e^{35}$', multialignment='center')
    plot_correlation(p6, '', UTC_time_mira, correlation_LR_MIRA, '', '.',
                     x_min=UTC_time_mira[1], x_max=UTC_time_mira[-1],
                     y_min=-1.1, y_max=1.1, x_lab='Time (UTC)', y_lab='correlation')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line = 'Comparison of LIMRad 94 GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC), '
    third_line = 'using: *.LV1.NC and *.mmclx data (unprocessed datasets)'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.8])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.8)

    file = date_str + '_Interpolation_LR_MIRA.png'
    if pts: print('')
    if pts: print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(file, dpi=dpi_val)
    # plt.show()
    plt.close()

#if plot_RectBivariateSpline:
#    # scipy interpolation test
#
#    x = LR_time[:]
#    y = LR_height[:]
#
#    Z = LR_Ze[:, :]
#
#    # grid the data
#    interp_spline = interpolate.RectBivariateSpline(x, y, Z)
#
#    # define grid
#    x2 = np.arange(LR_time[0], LR_time[-1], 1)
#    y2 = np.arange(LR_height[0], LR_height[-1], 0.005)
#    # print('LR_height first/last = ',LR_height[0], LR_height[-1])
#
#    time_x2 = []
#    for i in range(len(x2)):
#        time_x2.append(datetime.datetime(2001, 1, 1, 0, 0, 0)
#                       + datetime.timedelta(seconds=int(x2[i])))
#
#    X, Y = np.meshgrid(x, y)
#    X2, Y2 = np.meshgrid(x2, y2)
#
#    Z2_unmasked = interp_spline(x2, y2)
#    Z2 = np.ma.masked_less_equal(Z2_unmasked, -70.).T
#    # Z2 = np.ma.masked_greater_equal(Z2, -999.).T
#
#    fig, (LR_Ze_plot1, LR_Ze_plot2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8))
#
#    plot_data_set(fig, LR_Ze_plot1, 'Radar Reflectivity Factor',
#                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
#                  x_min=time_x2[0], x_max=time_x2[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='', y_lab='height (km)', z_lab='dBZ')
#
#    plot_data_set(fig, LR_Ze_plot2, 'Radar Reflectivity Factor',
#                  time_x2, y2, Z2, vmi=-50, vma=20,
#                  x_min=time_x2[0], x_max=time_x2[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')
#
#    fig.tight_layout()
#
#    #plt.show()

####################################################################################################################
#
#           ######## #### ##     ## ########          ######  ######## ########  #### ########  ######
#              ##     ##  ###   ### ##               ##    ## ##       ##     ##  ##  ##       ##    ##
#              ##     ##  #### #### ##               ##       ##       ##     ##  ##  ##       ##
#              ##     ##  ## ### ## ######   #######  ######  ######   ########   ##  ######    ######
#              ##     ##  ##     ## ##                     ## ##       ##   ##    ##  ##             ##
#              ##     ##  ##     ## ##               ##    ## ##       ##    ##   ##  ##       ##    ##
#              ##    #### ##     ## ########          ######  ######## ##     ## #### ########  ######
#
####################################################################################################################

if plot_radar_results:

    fig, plt = Plot_Radar_Results(LR_data, MIRA_data)

    ########################################################################################################
    # Save figure to file
    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_MIRA_LIMRad_profiles_timeseries.png'
    fig.savefig('C:/'+meteo_path+file, dpi=dpi_val)
    plt.close()

    print('    Save Figure to File :: ' + meteo_path + file + '\n')

######################################################################################################################################
#
#  ######## #### ##     ## ########         ##     ## ######## ####  ######   ##     ## ########            ###    ##     ##  ######
#     ##     ##  ###   ### ##               ##     ## ##        ##  ##    ##  ##     ##    ##              ## ##   ##     ## ##    ##
#     ##     ##  #### #### ##               ##     ## ##        ##  ##        ##     ##    ##             ##   ##  ##     ## ##
#     ##     ##  ## ### ## ######   ####### ######### ######    ##  ##   #### #########    ##    ####### ##     ## ##     ## ##   ####
#     ##     ##  ##     ## ##               ##     ## ##        ##  ##    ##  ##     ##    ##            #########  ##   ##  ##    ##
#     ##     ##  ##     ## ##               ##     ## ##        ##  ##    ##  ##     ##    ##            ##     ##   ## ##   ##    ##
#     ##    #### ##     ## ########         ##     ## ######## ####  ######   ##     ##    ##            ##     ##    ###     ######
#
######################################################################################################################################

if plot_comparisons:

    fig, plt = Plot_Comparison(LR_data, MIRA_data)

    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_MIRA_LIMRad94_avg_time_height.png'
    print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(file, dpi=dpi_val)
    plt.close()

if plot_interpolation_scatter:
    ########################################################################
    ### plot comparison ###



    fig, plt = Plot_Scatter(LR_data, MIRA_data)

    date_str = str(LR_data.year) + str(LR_data.month).zfill(2) + str(LR_data.day).zfill(2)
    file = date_str + '_MIRA_LIMRad94_interp-avgheight_comp.png'
    print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(file, dpi=dpi_val)
    plt.close()

if plot_compare_mira_mmclx:

    # print('')
    # print('Zh_mira-Ze_mmclx = ',np.sum(miramira_Z-mira_Ze[:,:-1]))

    print('')
    print('    Generate subplots:\n')

    font = FontProperties()

    fig = plt.figure(figsize=(16, 9))

    Zh_plot = plt.subplot2grid((3, 4), (0, 0))
    Ze_plot = plt.subplot2grid((3, 4), (0, 1))
    Zg_plot = plt.subplot2grid((3, 4), (0, 2))  # , rowspan=2)
    LR_Ze_plot = plt.subplot2grid((3, 4), (0, 3))  # , rowspan=2)

    v_plot = plt.subplot2grid((3, 4), (1, 0))  # , colspan=2)
    VEL_plot = plt.subplot2grid((3, 4), (1, 1))  # , colspan=2)
    VELg_plot = plt.subplot2grid((3, 4), (1, 2))  # , rowspan=2)
    LR_mdv_plot = plt.subplot2grid((3, 4), (1, 3))  # , rowspan=2)

    sw_plot = plt.subplot2grid((3, 4), (2, 0))  # , colspan=2, rowspan=2)
    RMS_plot = plt.subplot2grid((3, 4), (2, 1))  # , colspan=2, rowspan=2)
    RMSg_plot = plt.subplot2grid((3, 4), (2, 2))  # , rowspan=2)
    LR_sw_plot = plt.subplot2grid((3, 4), (2, 3))  # , rowspan=2)

    height_describtion = [-0.05, 1.3]
    place_text(Zh_plot, height_describtion, r'\textbf{Zh} ... Calibrated reflectivity. Calibration\\'
                                            ' convention: in the absence of attenuation, \n'
                                            'a cloud at 273 K containing one million 100-micron\n'
                                            ' droplets per cubic metre will have a reflectivity \n'
                                            'of 0 dBZ at all frequencies)\n'
                                            r'\textbf{v}  ... Radial component of the velocity,\\'
                                            ' with positive velocities are away from the radar.)\n '
                                            r'\textbf{width} ... Standard deviation of the reflectivity-\\'
                                            r'weighted velocities in the radar pulse volume.)')
    height_describtion = [0.1, 1.3]

    place_text(Ze_plot, height_describtion, r'\textbf{Ze} ... Equivalent Radar Reflectivity\\'
                                            ' Factor Ze of Hydrometeors  \n'
                                            r'\textbf{VEL} ... Doppler Velocity VEL\\'
                                            r'\textbf{RMS} ... Peak Width RMS')
    place_text(Zg_plot, height_describtion, r'\textbf{Zg} ... Equivalent Radar Reflectivity\\'
                                            ' Factor Ze of all Targets \n'
                                            r'\textbf{VELg} ... Doppler Velocity VELg\\'
                                            r'\textbf{RMSg} ... Peak Width RMSg')
    place_text(LR_Ze_plot, height_describtion, r'\textbf{Ze} ... Equivalent radar reflectivity factor\\'
                                               r'\textbf{mdv} ... Mean Doppler Velocity\\'
                                               r'\textbf{sw} ... Spectrum width')
    ################################################################################################################
    #
    ref_min = -60
    ref_max = 30
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   MIRA35 Reflectivity (Z)  ', end='', flush=True)
    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]

    Zh_plot.set_title(r'\textbf{Radar Reflectivity Factor (Zh)}')
    plot_data_set(fig, Zh_plot, '',
                  UTC_time_mira, mira_height, mira_Ze, vmi=ref_min, vma=ref_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    Ze_plot.set_title(r'\textbf{Radar Reflectivity Factor (Ze)}')
    plot_data_set(fig, Ze_plot, '',
                  UTC_time_mira, mira_height, mira_Ze, vmi=ref_min, vma=ref_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    Zg_plot.set_title(r'\textbf{Radar Reflectivity Factor (Zg)}')
    plot_data_set(fig, Zg_plot, '',
                  UTC_time_mira, mira_height, mira_Zg, vmi=ref_min, vma=ref_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    LR_Ze_plot.set_title(r'\textbf{Radar Reflectivity Factor (Ze)}')
    plot_data_set(fig, LR_Ze_plot, '',
                  UTC_time_LR, LR_height, LR_Ze, vmi=ref_min, vma=ref_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    ################################################################################################################
    #
    vel_min = -5
    vel_max = 3
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   MIRA35 Mean Doppler Velocity (VEL)  ', end='', flush=True)
    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]

    v_plot.set_title(r'\textbf{Mean Doppler Velocity (v)}')
    plot_data_set(fig, v_plot, '',
                  UTC_time_mira, mira_height, miraNC_VEL, vmi=vel_min, vma=vel_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    VEL_plot.set_title(r'\textbf{Mean Doppler Velocity (VEL)}')
    plot_data_set(fig, VEL_plot, '',
                  UTC_time_mira, mira_height, mira_VEL, vmi=vel_min, vma=vel_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    VELg_plot.set_title(r'\textbf{Mean Doppler Velocity (VELg)}')
    plot_data_set(fig, VELg_plot, '',
                  UTC_time_mira, mira_height, mira_VELg, vmi=vel_min, vma=vel_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    LR_mdv_plot.set_title(r'\textbf{Mean Doppler Velocity (mdv)}')
    plot_data_set(fig, LR_mdv_plot, '',
                  UTC_time_LR, LR_height, LR_mdv, vmi=vel_min, vma=vel_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    ################################################################################################################
    #
    sw_min = 10 ** (-1.5)
    sw_max = 10 ** (0.5)
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   MIRA35 Spectral Width (RMS)  ', end='', flush=True)
    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]

    sw_plot.set_title(r'\textbf{Spectral Width (width)}')
    plot_data_set(fig, sw_plot, 'sw',
                  UTC_time_mira, mira_height, miraNC_RMS, vmi=sw_min, vma=sw_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    RMS_plot.set_title(r'\textbf{Spectral Width (RMS)}')
    plot_data_set(fig, RMS_plot, 'sw',
                  UTC_time_mira, mira_height, mira_RMS, vmi=sw_min, vma=sw_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    RMSg_plot.set_title(r'\textbf{Spectral Width (RMSg)}')
    plot_data_set(fig, RMSg_plot, 'sw',
                  UTC_time_mira, mira_height, mira_RMSg, vmi=sw_min, vma=sw_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    LR_sw_plot.set_title(r'\textbf{Spectral Width (sw)}')
    plot_data_set(fig, LR_sw_plot, 'sw',
                  UTC_time_LR, LR_height, LR_sw, vmi=sw_min, vma=sw_max,
                  x_min=x_lim[0], x_max=x_lim[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line = 'Comparison of processed and unprocessed MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC), '
    third_line = 'using: *mira.nc data (1st column) and *.mmclx data (2nd - 4th column; no attenuation correction)'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.80])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)

    file = date_str + '_MIRA_mmclx_comparison.png'
    print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(file, dpi=dpi_val)

    plt.close()

save_log_data(file[:-5], interp_meth, hmin, hmax, date, time_intervall)
print('    Elapsed Time = {0:0.3f}'.format(time.clock() - start_time), '[sec]\n')
