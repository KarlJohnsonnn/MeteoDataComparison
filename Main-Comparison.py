from Graphics_Mod import *
from IO_Mod import *
from Parameter_Mod import *

rc('font', family='serif')
rc('text', usetex=True)

import datetime

import sys, warnings, time

import pandas as pd
from scipy import interpolate

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
plot_comparisons = False
plot_interpolation_scatter = False
plot_compare_mira_mmclx = False


os.chdir(meteo_path)

#################################################################################################
#
#       ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
#       ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
#       ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
#       ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
#       ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
#       ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
#       ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######
#
#################################################################################################

# Subroutine definition and start of the program
def dim(a):
    if not type(a) == list:  return []
    return [len(a)] + dim(a[0])


def correlation(v1, v2):
    rho = np.array([[]])
    for i in range(len(v1[0, :])):
        df_v1 = pd.DataFrame(v1[:, i])
        df_v2 = pd.DataFrame(v2[:, i])
        rho = np.append(rho, df_v1.corrwith(df_v2))

    return np.ma.masked_invalid(rho)


def interpolate_data(x, y, xnew, method):
    # create a callable function from the actual data
    fnew = interpolate.interp1d(x, y, kind=method)

    # calculate the interpolation dataset
    ynew = fnew(xnew)

    # mask values ( minInterpol=minAcutalData, etc. )
    ynew = np.ma.masked_greater_equal(ynew, y.max())
    ynew = np.ma.masked_less_equal(ynew, y.min())

    return ynew


def Interpolate_2D(x1, y1, z1, x2, y2, method):
    len_x1 = len(x1)
    len_x2 = len(x2)
    len_y1 = len(y1)
    len_y2 = len(y2)

    if method in ['linear', 'quadratic', 'cubic', 'quintic']:
        fcn = interpolate.interp2d(x1, y1, z1, kind=method)

        interp_z = fcn(x2, y2)
        interp_z = np.ma.masked_equal(interp_z, 0.0)
        interp_z = np.ma.masked_less_equal(interp_z, -100.0)
        interp_z = np.ma.masked_invalid(interp_z)


    elif method in ['NearestNeighbour']:

        coord = np.empty((len_x1 * len_y1, 2))
        values = np.empty((len_x1 * len_y1, 1))
        cnt = 0
        for i in range(len_x1):
            for j in range(len_y1):
                coord[cnt, 0] = x1[i]
                coord[cnt, 1] = y1[j]
                values[cnt] = z1[j, i]
                cnt += 1
                print('in copy loop  ', i, j, x1[i], y1[j], z1[j, i])

        values = np.ma.masked_invalid(values)
        fcn = interpolate.NearestNDInterpolator(coord, values)

        interp_z = np.empty((len_x2 * len_y2))
        cnt = 0
        for xi in x2:
            for yi in y2:
                interp_z[cnt] = fcn(xi, yi)
                cnt += 1
                print('in fcn loop  ', cnt, xi, yi)

        interp_z = np.ma.masked_equal(interp_z, 0.0)
        interp_z = np.ma.masked_less_equal(interp_z, -100.0)
        interp_z = np.ma.masked_invalid(interp_z)
        interp_z = np.reshape(interp_z, (len_x2, len_y2)).T

    return interp_z


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
    comp_date = str(sys.argv[3])
    comp_time_int = str(sys.argv[4]) + '-' + str(sys.argv[5])

else:

    ## cirrus
    #hmin = 8.50  # (km)  - lower y-axis limit
    #hmax = 10.0  # (km) - upper y-axis limit, highest range gate may be higher
    #comp_date = '180728'  # in YYMMDD
    #comp_time_int = '0740-0805'  # in HHMM-HHMM

    ##cummulis
    #hmin = 5.0  #(km)  - lower y-axis limit
    #hmax = 12.0  #(km) - upper y-axis limit, highest range gate may be higher
    #comp_date     = '180802'     # in YYMMDD
    #comp_time_int = '0330-1200'  # in HHMM-HHM

    hmin = 0.0 #(km)  - lower y-axis limit
    hmax = 12.00 #(km) - upper y-axis limit, highest range gate may be higher
    comp_date     = '180729'     # in YYMMDD
    comp_time_int = '2200-2359'  # in HHMM-HHMMM

    ##nimbus
    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    # comp_date     = '180805'     # in YYMMDD
    # comp_time_int = '0510-0620'  # in HHMM-HHMM

    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 3.00 #(km) - upper y-axis limit, highest range gate may be higher
    # comp_date     = '180805'     # in YYMMDD
    # comp_time_int = '1030-1200'  # in HHMM-HHMM

    # hmin = 1.0 #(km)  - lower y-axis limit
    # hmax = 2.5 #(km) - upper y-axis limit, highest range gate may be higher
    # comp_date     = '180805'     # in YYMMDD
    # comp_time_int = '0700-1210'  # in HHMM-HHMM

    ##nimbus
    # hmin = 0.0 #(km)  - lower y-axis limit
    # hmax = 12.00 #(km) - upper y-axis limit, highest range gate may be higher
    # comp_date     = '180808'     # in YYMMDD
    # comp_time_int = '1330-1700'  # in HHMM-HHMM

warnings.filterwarnings("ignore")

# calculate the time in decimal hours
comp_hours = [int(comp_time_int[0:2]), int(comp_time_int[5:7])]
comp_minutes = [int(comp_time_int[2:4]), int(comp_time_int[7:9])]

clock_time = np.array(comp_hours) + np.divide(comp_minutes, 60.)  # [hours] + [minutes]/60#

# -- gathering year, month, day for convertion to UTC time
plotyear = int('20' + comp_date[:2])
plotmonth = int(comp_date[2:4])
plotday = int(comp_date[4:6])

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

LR_time, UTC_time_LR, LR_height, \
LR_Ze, LR_mdv, LR_sw = extract_dataset(comp_date, time_int, clock_time, height, '*.LV1.NC', '')

# ----- MIRA 35GHz Radar data extraction

# mira.nc data file (processed data)
mira_timeNC, UTC_time_miraNC, mira_heightNC, \
miraNC_Z, miraNC_VEL, miraNC_RMS = extract_dataset(comp_date, time_int, clock_time, height, '*mira.nc', '')

# .mmclx data file (hydrometeors only)
_, _, _, mira_Ze, mira_VEL, mira_RMS = extract_dataset(comp_date, time_int, clock_time, height, '*.mmclx', '')

# .mmclx data file (all targets)
mira_time, UTC_time_mira, mira_height, \
mira_Zg, mira_VELg, mira_RMSg = extract_dataset(comp_date, time_int, clock_time, height, '*.mmclx', 'g')

print('')

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
LR_timeavg_Ze = np.average(LR_Ze, axis=1)
LR_timeavg_mdv = np.average(LR_mdv, axis=1)
LR_timeavg_sw = np.average(LR_sw, axis=1)
mira_timeavg_Ze = np.average(mira_Zg, axis=1)
mira_timeavg_mdv = np.average(mira_VELg, axis=1)
mira_timeavg_sw = np.average(mira_RMSg, axis=1)

# height averaged values
LR_heightavg_Ze = np.average(LR_Ze, axis=0)
LR_heightavg_mdv = np.average(LR_mdv, axis=0)
LR_heightavg_sw = np.average(LR_sw, axis=0)
mira_heightavg_Ze = np.average(mira_Zg, axis=0)
mira_heightavg_mdv = np.average(mira_VELg, axis=0)
mira_heightavg_sw = np.average(mira_RMSg, axis=0)

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

    LRtoMIRA_Ze = Interpolate_2D(LR_time, LR_height, LR_Ze, mira_time, mira_height, interp_meth)
    MIRAtoLR_Zg = Interpolate_2D(mira_time, mira_height, mira_Ze, LR_time, LR_height, interp_meth)

    # converting back to linear untis for mean difference calculation
    LR_Ze_mm6m3_I = np.power(np.divide(LRtoMIRA_Ze, 10.0), 10.0)
    LR_Ze_mm6m3 = np.power(np.divide(LR_Ze, 10.0), 10.0)
    mira_Zg_mm6m3_I = np.power(np.divide(MIRAtoLR_Zg, 10.0), 10.0)
    mira_Zg_mm6m3 = np.power(np.divide(mira_Zg, 10.0), 10.0)

    # mean_diff_LRtoM = np.mean(LRtoMIRA_Ze - mira_Zg)
    # mean_diff_MtoLR = np.mean(LR_Ze - MIRAtoLR_Zg)

    # display mean difference in dBZ again
    mean_diff_LRtoM = 10 * np.log10(np.mean(mira_Zg_mm6m3 - LR_Ze_mm6m3_I))
    mean_diff_MtoLR = 10 * np.log10(np.mean(mira_Zg_mm6m3_I - LR_Ze_mm6m3))

    # correlation is calculated with logarithmic units
    correlation_LR_MIRA = correlation(LRtoMIRA_Ze, mira_Zg)
    correlation_MIRA_LR = correlation(LR_Ze, MIRAtoLR_Zg)

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
                  UTC_time_LR, LR_height, MIRAtoLR_Zg, vmi=-50, vma=20,
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

if plot_RectBivariateSpline:
    # scipy interpolation test

    x = LR_time[:]
    y = LR_height[:]

    Z = LR_Ze[:, :]

    # grid the data
    interp_spline = interpolate.RectBivariateSpline(x, y, Z)

    # define grid
    x2 = np.arange(LR_time[0], LR_time[-1], 1)
    y2 = np.arange(LR_height[0], LR_height[-1], 0.005)
    # print('LR_height first/last = ',LR_height[0], LR_height[-1])

    time_x2 = []
    for i in range(len(x2)):
        time_x2.append(datetime.datetime(2001, 1, 1, 0, 0, 0)
                       + datetime.timedelta(seconds=int(x2[i])))

    X, Y = np.meshgrid(x, y)
    X2, Y2 = np.meshgrid(x2, y2)

    Z2_unmasked = interp_spline(x2, y2)
    Z2 = np.ma.masked_less_equal(Z2_unmasked, -70.).T
    # Z2 = np.ma.masked_greater_equal(Z2, -999.).T

    fig, (LR_Ze_plot1, LR_Ze_plot2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8))

    plot_data_set(fig, LR_Ze_plot1, 'Radar Reflectivity Factor',
                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
                  x_min=time_x2[0], x_max=time_x2[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='', y_lab='height (km)', z_lab='dBZ')

    plot_data_set(fig, LR_Ze_plot2, 'Radar Reflectivity Factor',
                  time_x2, y2, Z2, vmi=-50, vma=20,
                  x_min=time_x2[0], x_max=time_x2[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    fig.tight_layout()

    #plt.show()

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
    ### plot ###
    print('    Generate subplots:\n')

    # create figure
    font = FontProperties()

    fig = plt.figure(figsize=(16, 10))

    LR_Ze_plot = plt.subplot2grid((3, 2), (0, 0))
    LR_mdv_plot = plt.subplot2grid((3, 2), (1, 0))  # , rowspan=2)
    LR_sw_plot = plt.subplot2grid((3, 2), (2, 0))  # , rowspan=2)
    mira_Zg_plot = plt.subplot2grid((3, 2), (0, 1))  # , colspan=2)
    mira_VELg_plot = plt.subplot2grid((3, 2), (1, 1))  # , rowspan=2)
    mira_RMSg_plot = plt.subplot2grid((3, 2), (2, 1))  # , rowspan=2)

    x_lim_left_LR = UTC_time_LR[0]
    x_lim_right_LR = UTC_time_LR[-1]

    x_lim_left_mira = UTC_time_mira[0]
    x_lim_right_mira = UTC_time_mira[-1]

    ########################################################################################################
    ########################################################################################################
    # LR_Zelectivity plot
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    # LIMRad reflectivity
    LR_Ze_plot.set_title(r'\textbf{LIMRad94')

    plot_data_set(fig, LR_Ze_plot, 'Radar Reflectivity Factor',
                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
                  x_min=x_lim_left_LR, x_max=x_lim_right_LR,
                  y_min=hmin, y_max=hmax,
                  x_lab='', y_lab='height (km)', z_lab='dBZ')

    # MIRA reflectivit
    mira_Zg_plot.set_title(r'\textbf{MIRA35}')

    plot_data_set(fig, mira_Zg_plot, 'Radar Reflectivity Factor',
                  UTC_time_mira, mira_height, mira_Zg, vmi=-50, vma=20,
                  x_min=x_lim_left_mira, x_max=x_lim_right_mira,
                  y_min=hmin, y_max=hmax,
                  x_lab='', y_lab='Height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # mean doppler velocity plot
    if pts: print('       -   Mean Doppler velocity   ', end='', flush=True)

    # LIMRad mean Doppler velocity
    plot_data_set(fig, LR_mdv_plot, 'Mean Doppler Velocity',
                  UTC_time_LR, LR_height, LR_mdv, vmi=-4, vma=2,
                  x_min=x_lim_left_LR, x_max=x_lim_right_LR,
                  y_min=hmin, y_max=hmax,
                  x_lab='', y_lab='Height (km)', z_lab='m/s')

    # MIRA mean Doppler velocity
    plot_data_set(fig, mira_VELg_plot, 'Mean Doppler Velocity',
                  UTC_time_mira, mira_height, mira_VELg, vmi=-4, vma=2,
                  x_min=x_lim_left_mira, x_max=x_lim_right_mira,
                  y_min=hmin, y_max=hmax,
                  x_lab='', y_lab='Height (km)', z_lab='m/s')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # spectral width plot
    if pts: print('       -   Spectral Width   ', end='', flush=True)

    # LIMRad spectral width
    plot_data_set(fig, LR_sw_plot, 'Spectral Width',
                  UTC_time_LR, LR_height, LR_sw, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=x_lim_left_LR, x_max=x_lim_right_LR,
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='m/s')

    # MIRA spectral widthAAAAAA
    plot_data_set(fig, mira_RMSg_plot, 'Spectral Width',
                  UTC_time_mira, mira_height, mira_RMSg, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=x_lim_left_mira, x_max=x_lim_right_mira,
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen

    ########################################################################################################
    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC),'
    third_line = r'using: ' + LIMRad_file_extension + ' and ' + mmclx_file_extension + ' data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)

    #plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.15)


    #plt.show()

    file = date_str + '_MIRA_LIMRad94_profiles_ts_comp.png'
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

    ########################################################################
    ### plot comparison ###

    print('    Generate subplots:\n')

    # create figure
    font = FontProperties()
    fig, ((LR_Ze_plot, mira_Ze_plot,),
          (Comp_avgT_Ze_plot, Comp_avgH_Ze_plot,),
          (Comp_avgT_mdv_plot, Comp_avgH_mdv_plot,),
          (Comp_avgT_sw_plot, Comp_avgH_sw_plot)) = plt.subplots(4, 2, figsize=(16, 12))

    # calculate x axis limits (same for both time series)

    x_lim_left_time = min(UTC_time_LR[0], UTC_time_mira[0])
    x_lim_right_time = max(UTC_time_LR[-1], UTC_time_mira[-1])

    ################################################################################################################
    ################################################################################################################
    #
    # PLOT LIMRAD and MIRA Reflectivity (Ze) DATASET
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    # LIMRad reflectivity
    LR_Ze_plot.set_title(r'\textbf{LIMRad94}')

    plot_data_set(fig, LR_Ze_plot, 'Radar Reflectivity Factor',
                  UTC_time_LR, LR_height, LR_Ze, vmi=-50, vma=20,
                  x_min=x_lim_left_time, x_max=x_lim_right_time,
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    # MIRA reflectivit
    mira_Ze_plot.set_title(r'\textbf{MIRA35}')

    plot_data_set(fig, mira_Ze_plot, 'Radar Reflectivity Factor',
                  UTC_time_mira, mira_height, mira_Ze, vmi=-50, vma=20,
                  x_min=x_lim_left_time, x_max=x_lim_right_time,
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   comp. Radar Reflectivity Factor   ', end='', flush=True)

    x_lim_left_Ze = min([LR_timeavg_Ze.min(), mira_timeavg_Ze.min()])
    x_lim_right_Ze = max([LR_timeavg_Ze.max(), mira_timeavg_Ze.max()])

    Comp_avgT_Ze_plot.set_title(r'\textbf{Time-Mean}')

    plot_avg_data_set(Comp_avgT_Ze_plot, 'Radar Reflectivity Factor',
                      LR_timeavg_Ze, LR_height,
                      mira_timeavg_Ze, mira_height,
                      label1='LIMRad', marker1='+', label2='MIRA', marker2='o',
                      x_min=x_lim_left_Ze, x_max=x_lim_right_Ze,
                      y_min=hmin, y_max=hmax, x_lab='dBZ', y_lab='height (km)', ax='y')

    Comp_avgH_Ze_plot.set_title(r'\textbf{Height-Mean}')

    plot_avg_data_set(Comp_avgH_Ze_plot, 'Radar Reflectivity Factor',
                      UTC_time_LR, LR_heightavg_Ze,
                      UTC_time_mira, mira_heightavg_Ze,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=x_lim_left_time, x_max=x_lim_right_time,
                      y_min=[], y_max=[], x_lab='Time (UTC)', y_lab='dBZ', ax='y')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   comp. Mean Doppler Velocity   ', end='', flush=True)

    x_lim_left_mdv = min([LR_timeavg_mdv.min(), mira_timeavg_mdv.min()])
    x_lim_right_mdv = max([LR_timeavg_mdv.max(), mira_timeavg_mdv.max()])

    plot_avg_data_set(Comp_avgT_mdv_plot, 'Doppler Velocity',
                      LR_timeavg_mdv, LR_height,
                      mira_timeavg_mdv, mira_height,
                      label1='LIMRad', marker1='+', label2='MIRA', marker2='o',
                      x_min=x_lim_left_mdv, x_max=x_lim_right_mdv,
                      y_min=hmin, y_max=hmax, x_lab='m/s', y_lab='height (km)', ax='y')

    plot_avg_data_set(Comp_avgH_mdv_plot, 'Doppler Velocity',
                      UTC_time_LR, LR_heightavg_mdv,
                      UTC_time_mira, mira_heightavg_mdv,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=x_lim_left_time, x_max=x_lim_right_time,
                      y_min=[], y_max=[], x_lab='Time (UTC)', y_lab='m/s', ax='y')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Spectral Width LIMRAD-MIRA
    if pts: print('       -   comp. Spectral Width   ', end='', flush=True)

    x_lim_left_sw = min([LR_timeavg_sw.min(), mira_timeavg_sw.min()])
    x_lim_right_sw = max([LR_timeavg_sw.max(), mira_timeavg_sw.max()])

    plot_avg_data_set(Comp_avgT_sw_plot, 'Spectral Width',
                      LR_timeavg_sw, LR_height,
                      mira_timeavg_sw, mira_height,
                      label1='LIMRad', marker1='+', label2='MIRA', marker2='o',
                      x_min=x_lim_left_sw, x_max=x_lim_right_sw,
                      y_min=hmin, y_max=hmax, x_lab='m/s', y_lab='height (km)', ax='y')

    plot_avg_data_set(Comp_avgH_sw_plot, 'Spectral Width',
                      UTC_time_LR, LR_heightavg_sw,
                      UTC_time_mira, mira_heightavg_sw,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=x_lim_left_time, x_max=x_lim_right_time,
                      y_min=[], y_max=[], x_lab='Time (UTC)', y_lab='m/s', ax='y')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC),'
    third_line = r'using: ' + LIMRad_file_extension + ' and ' + mira_file_extension + ' data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.35)

    file = date_str + '_MIRA_LIMRad94_t-h_comp.png'
    print('    Save Figure to File :: ' + meteo_path + file + '\n')
    fig.savefig(file, dpi=dpi_val)
    plt.close()

if plot_interpolation_scatter:
    ########################################################################
    ### plot comparison ###

    interp_meth = 'nearest'
    res_interp = 5  # in [sec]
    head_pos = [0.05]
    stat_pos = [0.25, -0.3]

    # create an array with evenly spaced gridsize
    xnew = np.arange(max(LR_time[0], mira_time[0]),
                     min(LR_time[-1], mira_time[-1]),
                     res_interp)

    # convert the x-axis unix time to readable date time format
    plot_time_xnew = [datetime.datetime(1970, 1, 1, 0, 0, 0)
                      + datetime.timedelta(seconds=int(xnew[i])) for i in range(len(xnew))]

    # what do you want to interpolate
    LR_x = LR_time;
    mira_x = mira_time
    LR_y = LR_heightavg_Ze;
    mira_y = mira_heightavg_Ze

    LR_ynew = interpolate_data(LR_x, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(mira_x, mira_y, xnew, interp_meth)

    # y width +-5
    y_min, y_max = get_plot_ybounds(LR_y, mira_y, 7.0)

    # calculate the mean difference and covariance matrix
    mean_diff = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef = np.corrcoef(LR_ynew, mira_ynew)

    print('    Generate subplots:\n')

    font = FontProperties()

    fig = plt.figure(figsize=(16, 12))

    h_Ze_plot = plt.subplot2grid((3, 3), (0, 0))
    h_mdv_plot = plt.subplot2grid((3, 3), (0, 1))  # , rowspan=2)
    h_sw_plot = plt.subplot2grid((3, 3), (0, 2))  # , rowspan=2)
    interp_h_Ze_plot = plt.subplot2grid((3, 3), (1, 0))  # , colspan=2)
    interp_h_mdv_plot = plt.subplot2grid((3, 3), (1, 1))  # , rowspan=2)
    interp_h_sw_plot = plt.subplot2grid((3, 3), (1, 2))  # , rowspan=2)
    scatter_Ze = plt.subplot2grid((3, 3), (2, 0))  # , colspan=2, rowspan=2)
    scatter_mdv = plt.subplot2grid((3, 3), (2, 1))  # , rowspan=2)
    scatter_sw = plt.subplot2grid((3, 3), (2, 2))  # , rowspan=2)

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   Average reflectivity over height domain  ', end='', flush=True)

    h_Ze_plot.set_title(r' \textbf{Mean-Height Reflectivity}')
    plot_avg_data_set(h_Ze_plot, '',
                      UTC_time_LR, LR_heightavg_Ze,
                      UTC_time_mira, mira_heightavg_Ze,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                      y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='dBZ', ax='n')

    interp_h_Ze_plot.set_title(r' \textbf{Mean-Height Reflectivity}')
    plot_interpol_data_set(interp_h_Ze_plot, '',
                           plot_time_xnew, LR_ynew,
                           plot_time_xnew, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                           y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='dBZ')

    place_statistics(interp_h_Ze_plot, stat_pos, [mean_diff, cor_coef[0, 1]], 'Ze')

    scatter_Ze.set_title('Scatter Plot of Height-Mean\n Reflectivity ')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 7.0)

    plot_scatter(scatter_Ze, '', LR_ynew, mira_ynew, '*'
                 , x_min=xy_min, x_max=xy_max,
                 y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in dBZ', y_lab='MIRA data in dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    # same for mean doppler velocity

    LR_y = LR_heightavg_mdv
    mira_y = mira_heightavg_mdv

    LR_ynew = interpolate_data(LR_x, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(mira_x, mira_y, xnew, interp_meth)
    LR_ynew = np.ma.masked_equal(LR_ynew, 0.0)
    mira_ynew = np.ma.masked_equal(mira_ynew, 0.0)

    # y width +-5
    y_min, y_max = get_plot_ybounds(LR_y, mira_y, 0.1)

    # calculate the mean difference and covariance matrix
    mean_diff = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef = np.corrcoef(LR_ynew, mira_ynew)

    ################################################################################################################
    #
    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   Average mean doppler velocity over height domain  ', end='', flush=True)

    place_text(h_mdv_plot, [-1.1, 1.2], r'\LARGE{\textbf{Actual Dataset}}')

    h_mdv_plot.set_title(r' \textbf{Height-Mean Mean Doppler Velocity}')
    plot_avg_data_set(h_mdv_plot, '',
                      UTC_time_LR, LR_heightavg_mdv,
                      UTC_time_mira, mira_heightavg_mdv,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                      y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='m/s', ax='n')

    place_text(interp_h_mdv_plot, [-1.1, 1.2], r'\LARGE{\textbf{Interpolated Dataset}}')

    interp_h_mdv_plot.set_title(r' \textbf{Height-Mean Mean Doppler Velocity}')
    plot_interpol_data_set(interp_h_mdv_plot, '',
                           plot_time_xnew, LR_ynew,
                           plot_time_xnew, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                           y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='m/s')

    place_statistics(interp_h_mdv_plot, stat_pos, [mean_diff, cor_coef[0, 1]], 'mdv')
    scatter_mdv.set_title('Scatter Plot of Height-Mean\n Mean Doppler Velocity')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.1)

    plot_scatter(scatter_mdv, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max,
                 y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    # same for mean doppler velocity

    LR_y = LR_heightavg_sw
    mira_y = mira_heightavg_sw

    LR_ynew = interpolate_data(LR_x, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(mira_x, mira_y, xnew, interp_meth)

    # y width +-5
    y_min, y_max = get_plot_ybounds(LR_y, mira_y, 0.03)

    # calculate the mean difference and covariance matrix
    mean_diff = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef = np.corrcoef(LR_ynew, mira_ynew)

    ################################################################################################################
    #
    # comparsion of spectral width LIMRAD-MIRA
    if pts: print('       -   Average spectral width over height domain  ', end='', flush=True)

    h_sw_plot.set_title(r'\textbf{Height-Mean Spectral Width}')
    plot_avg_data_set(h_sw_plot, '',
                      UTC_time_LR, LR_heightavg_sw,
                      UTC_time_mira, mira_heightavg_sw,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                      y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='m/s', ax='n')

    interp_h_sw_plot.set_title(r'\textbf{Height-Mean Spectral Width}')
    plot_interpol_data_set(interp_h_sw_plot, '',
                           plot_time_xnew, LR_ynew,
                           plot_time_xnew, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=UTC_time_LR[0], x_max=UTC_time_LR[-1],
                           y_min=y_min, y_max=y_max, x_lab='Time (UTC)', y_lab='m/s')

    place_statistics(interp_h_sw_plot, stat_pos, [mean_diff, cor_coef[0, 1]], 'sw')

    scatter_sw.set_title('Scatter Plot of Height-Mean\n Spectral Width ')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.03)

    plot_scatter(scatter_sw, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max,
                 y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file
    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
    first_line = 'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:x  ' + str(time_int[3]) + ' (UTC), '
    third_line = 'using: ' + LIMRad_file_extension + ' and ' + mira_file_extension + ' data,  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.91])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.55)

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

save_log_data(file[:-5], interp_meth, hmin, hmax, comp_date, comp_time_int)
print('    Elapsed Time = {0:0.3f}'.format(time.clock() - start_time), '[sec]\n')
