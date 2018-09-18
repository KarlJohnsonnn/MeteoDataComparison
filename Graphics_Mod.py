import datetime
import matplotlib.pyplot as plt
import matplotlib        as mpl

from matplotlib import rc
from matplotlib import colors

rc('font', family='serif')
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from matplotlib import dates
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Interpolation_Tool import *

from Parameter_Mod import pts


def plot_data_set(fig, axh, text, x, y, z, vmi, vma, x_min, x_max, y_min, y_max, x_lab, y_lab, z_lab):
    text = r'\textbf{' + text + '}'
    if text.find('sw') > 0 or text.find('Spectral Width') > 0:
        cp = axh.pcolormesh(x, y, z, norm=colors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
    else:
        place_text(axh, [.02, 1.05], text)
        cp = axh.pcolormesh(x, y, z, vmin=vmi, vmax=vma, cmap='jet')
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax0 = divider1.append_axes("right", size="3%", pad=0.1)
    cbar = fig.colorbar(cp, cax=cax0, ax=axh)
    cbar.set_label(z_lab)
    axh.axes.tick_params(axis='both', direction='inout', length=10, width=1.5)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.set_ylim(bottom=y_min, top=y_max, )

    # exceptions
    if x_lab == '':
        axh.axes.xaxis.set_ticklabels([])
    else:
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
        axh.set_xlabel(x_lab)


def plot_correlation_matrix(axh, text, z, vmi, vma, x_lab, y_lab, z_lab):
    from pylab import pcolor
    # text = r'\textbf{'+text+'}'
    # place_text(axh, [.02, 1.05], text )
    # z = z[1500:2000,1500:2000]
    cp = axh.pcolor(z, norm=colors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax0 = divider1.append_axes("right", size="3%", pad=0.1)
    # cbar= fig.colorbar(cp, cax=cax0, ax=axh)
    # cbar.set_label(z_lab)
    # Axh.axes.tick_params(axis='both', direction='inout', length=10, width=1.5)

    # exceptions
    # axh.axes.xaxis.set_ticklabels([]) if x_lab == '' else axh.set_xlabel(x_lab)
    # axh.axes.yaxis.set_ticklabels([]) if y_lab == '' else axh.set_ylabel(y_lab)

    # axh.set_aspect('equal', 'box')


def plot_correlation(axh, text, x, y, label, marker, x_min, x_max, y_min, y_max, x_lab, y_lab):
    place_text(axh, [.15, 1.1], text)
    axh.plot(x, y, marker, label=label)
    axh.set_xlabel(x_lab)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    divider1 = make_axes_locatable(axh)
    cax = divider1.append_axes("right", size="3%", pad=0.1)
    axh.legend(loc="upper right")

    cax.set_facecolor('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        cax.spines[axis].set_linewidth(0)
    cax.set_xticks([])
    cax.set_yticks([])

    # exceptions
    if not (y_min == y_max):
        axh.set_ylim(bottom=y_min, top=y_max)
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))


def plot_avg_data_set(axh, text, x1, y1, x2, y2, label1, marker1, label2, marker2, x_min, x_max, y_min, y_max, x_lab,
                      y_lab, ax):
    text = r'\textbf{' + text + '}'
    place_text(axh, [.02, 1.05], text)
    axh.scatter(x1, y1, marker=marker1, label=label1)
    axh.scatter(x2, y2, marker=marker2, label=label2)
    axh.set_xlabel(x_lab)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    if ax == 'y':
        divider1 = make_axes_locatable(axh)
        cax = divider1.append_axes("right", size="3%", pad=0.1)
        axh.legend(loc="upper right")

        cax.set_facecolor('none')
        for axis in ['top', 'bottom', 'left', 'right']:
            cax.spines[axis].set_linewidth(0)
        cax.set_xticks([])
        cax.set_yticks([])

    # exceptions
    if not (y_min == y_max):
        axh.set_ylim(bottom=y_min, top=y_max)
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))


def plot_phase_data_set(axh, text, x, y, marker, x_min, x_max, y_min, y_max, x_lab, y_lab):
    text = r'\textbf{' + text + '}'
    place_text(axh, [.02, 1.05], text)
    axh.scatter(x, y, marker=marker)
    axh.set_xlabel(x_lab)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax = divider1.append_axes("right", size="3%", pad=0.25)
    axh.legend(loc="upper right")

    cax.set_facecolor('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        cax.spines[axis].set_linewidth(0)
    cax.set_xticks([])
    cax.set_yticks([])

    # exceptions
    if not (y_min == y_max):
        axh.set_ylim(bottom=y_min, top=y_max)
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))


def plot_interpol_data_set(axh, text, x1, y1, x2, y2, label1, marker1, label2, marker2, x_min, x_max, y_min, y_max,
                           x_lab, y_lab):
    place_text(axh, [.15, 1.1], text)
    axh.plot(x1, y1, marker1, label=label1)
    axh.plot(x2, y2, marker2, label=label2)
    axh.set_xlabel(x_lab)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    # exceptions
    if not (y_max == y_min):
        axh.set_ylim(bottom=y_min, top=y_max)
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))


def plot_scatter(axh, text, x, y, marker, x_min, x_max, y_min, y_max, x_lab, y_lab):
    import matplotlib.ticker as ticker
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    # place_text(axh, [.05, 1.1], text )
    axh.plot(x, y, marker)
    axh.set_xlabel(x_lab)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.set_ylim(bottom=y_min, top=y_max)
    axh.xaxis.set_ticks((np.arange(x_min, x_max, (x_max - x_min) / 4.0)))
    axh.yaxis.set_ticks((np.arange(y_min, y_max, (y_max - y_min) / 4.0)))
    axh.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    axh.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    axh.grid(linestyle=':')
    # axh.legend(loc="upper right")
    axh.set_aspect('equal', 'box')

    # plot 1:1 line
    N = 100
    X = np.linspace(x_min, x_max, N)
    axh.plot(X, X, 'k--')

    # add patches
    colors = [np.divide([31, 119, 180], 255.), np.divide([255, 127, 14], 255.)]
    patches = [Polygon([[x_min, y_min], [x_max, y_min], [x_max, y_max]], facecolor='C0', fill=True),
               Polygon([[x_min, y_min], [x_min, y_max], [x_max, y_max]], facecolor='C1', fill=True)]

    p = PatchCollection(patches, alpha=0.25)
    p.set_color(colors)
    axh.add_collection(p)


def get_plot_ybounds(y1, y2, pm):
    return max(y1.min(), y2.min()) - pm, min(y1.max(), y2.max()) + pm


def place_text(plot, pos, text):
    plot.text(pos[0], pos[1], text,
              fontweight='bold',
              horizontalalignment='left',
              transform=plot.transAxes,
              bbox=dict(facecolor='white',
                        edgecolor='black',
                        pad=5.)
              )


def place_statistics(plot, pos, stat, vn):
    if vn == 'Ze':
        text = r'$\overline{  Z_g^{35} - Z_e^{94} }=$' + \
               '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'$\rho(Z_e^{94}, Z_g^{35}) = $' + '{:6.2f}'.format(stat[1])
    if vn == 'mdv':
        text = r'$\mathrm{mean}_h(\mathrm{mdv}_{\mathrm{lr}} - \mathrm{mdv}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)' + r'$=$' + '{:6.2f}'.format(stat[1])
    if vn == 'sw':
        text = r'$\mathrm{mean}_h(\mathrm{sw}_{\mathrm{lr}} - \mathrm{sw}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)' + r'$=$' + '{:6.2f}'.format(stat[1])

    plot.text(pos[0], pos[1], text, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))
    plot.text(pos[0] + 0.6, pos[1], text2, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))


def place_statistics_mean_corcoef(plot, pos, stat, vn):
    if vn == 'Ze':
        text = r'$\overline{Z_e^{94} - Z_e^{35} } =$' + \
               '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'corr$(Z_e^{94}, Z_e^{35}) = $' + r'$=$' + '{:6.2f}'.format(stat[1])
    if vn == 'mdv':
        text = r'$\mathrm{mean}_h(\mathrm{mdv}_{\mathrm{lr}} - \mathrm{mdv}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)' + r'$=$' + '{:6.2f}'.format(stat[1])
    if vn == 'sw':
        text = r'$\mathrm{mean}_h(\mathrm{sw}_{\mathrm{lr}} - \mathrm{sw}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)' + r'$=$' + '{:6.2f}'.format(stat[1])

    plot.text(pos[0], pos[1], text, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))
    plot.text(pos[0] + 0.6, pos[1], text2, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))


def Plot_Radar_Results(ds1, ds2):
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

    xb1 = [ds1.t_plt[0], ds1.t_plt[-1]]
    xb2 = [ds2.t_plt[0], ds2.t_plt[-1]]

    yb1 = [ds1.height[0], ds1.height[-1]]
    yb2 = [ds2.height[0], ds2.height[-1]]

    ########################################################################################################
    ########################################################################################################
    # LR_Zelectivity plot
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    LR_Ze_plot.set_title(r'\textbf{LIMRad94')

    plot_data_set(fig, LR_Ze_plot, 'Radar Reflectivity Factor',
                  ds1.t_plt, ds1.height, ds1.Ze, vmi=-50, vma=20,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab='height (km)', z_lab='dBZ')

    mira_Zg_plot.set_title(r'\textbf{MIRA35}')

    plot_data_set(fig, mira_Zg_plot, 'Radar Reflectivity Factor',
                  ds2.t_plt, ds2.height, ds2.Ze, vmi=-50, vma=20,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab='Height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # mean doppler velocity plot
    if pts: print('       -   Mean Doppler velocity   ', end='', flush=True)

    # LIMRad mean Doppler velocity
    plot_data_set(fig, LR_mdv_plot, 'Mean Doppler Velocity',
                  ds1.t_plt, ds1.height, ds1.mdv, vmi=-4, vma=2,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab='Height (km)', z_lab='m/s')

    # MIRA mean Doppler velocity
    plot_data_set(fig, mira_VELg_plot, 'Mean Doppler Velocity',
                  ds2.t_plt, ds2.height, ds2.mdv, vmi=-4, vma=2,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab='Height (km)', z_lab='m/s')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # spectral width plot
    if pts: print('       -   Spectral Width   ', end='', flush=True)

    # LIMRad spectral width
    plot_data_set(fig, LR_sw_plot, 'Spectral Width',
                  ds1.t_plt, ds1.height, ds1.mdv, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='m/s')

    # MIRA spectral widthAAAAAA
    plot_data_set(fig, mira_RMSg_plot, 'Spectral Width',
                  ds2.t_plt, ds2.height, ds2.mdv, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen

    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(xb1[0]) + ' (UTC)  to:  ' + str(xb1[1]) + ' (UTC),'
    third_line = r'using: LIMRad94 and MIRA35 data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)

    # plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.15)

    return fig, plt


def Plot_Comparison(ds1, ds2):
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

    xb1 = [ds1.t_plt[0], ds1.t_plt[-1]]
    xb2 = [ds2.t_plt[0], ds2.t_plt[-1]]

    yb1 = [ds1.height[0], ds1.height[-1]]
    yb2 = [ds2.height[0], ds2.height[-1]]

    ################################################################################################################
    ################################################################################################################
    #
    # PLOT LIMRAD and MIRA Reflectivity (Ze) DATASET
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    # LIMRad reflectivity
    LR_Ze_plot.set_title(r'\textbf{LIMRad94}')

    plot_data_set(fig, LR_Ze_plot, 'Radar Reflectivity Factor',
                  ds1.t_plt, ds1.height, ds1.Ze, vmi=-50, vma=20,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    # MIRA reflectivit
    mira_Ze_plot.set_title(r'\textbf{MIRA35}')

    plot_data_set(fig, mira_Ze_plot, 'Radar Reflectivity Factor',
                  ds2.t_plt, ds2.height, ds2.Ze, vmi=-50, vma=20,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   comp. Radar Reflectivity Factor   ', end='', flush=True)

    x_lim_left_Ze = min([ds1.timeavg_Ze.min(), ds2.timeavg_Ze.min()])
    x_lim_right_Ze = max([ds1.timeavg_Ze.max(), ds2.timeavg_Ze.max()])

    Comp_avgT_Ze_plot.set_title(r'\textbf{Time-Mean}')

    plot_avg_data_set(Comp_avgT_Ze_plot, 'Radar Reflectivity Factor',
                      ds1.timeavg_Ze, ds1.height, ds2.timeavg_Ze, ds2.height,
                      label1='LIMRad', marker1='+', label2='MIRA', marker2='o',
                      x_min=x_lim_left_Ze, x_max=x_lim_right_Ze,
                      y_min=yb1[0], y_max=yb1[1], x_lab='dBZ', y_lab='height (km)', ax='y')

    Comp_avgH_Ze_plot.set_title(r'\textbf{Height-Mean}')

    plot_avg_data_set(Comp_avgH_Ze_plot, 'Radar Reflectivity Factor',
                      ds1.t_plt, ds1.heightavg_Ze, ds2.t_plt, ds2.heightavg_Ze,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=[], y_max=[], x_lab='Time (UTC)', y_lab='dBZ', ax='y')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   comp. Mean Doppler Velocity   ', end='', flush=True)

    x_lim_left_mdv = min([ds1.timeavg_mdv.min(), ds2.timeavg_mdv.min()])
    x_lim_right_mdv = max([ds1.timeavg_mdv.max(), ds2.timeavg_mdv.max()])

    plot_avg_data_set(Comp_avgT_mdv_plot, 'Doppler Velocity',
                      ds1.timeavg_mdv, ds1.height, ds2.timeavg_mdv, ds2.height,
                      label1='LIMRad', marker1='+', label2='MIRA', marker2='o',
                      x_min=x_lim_left_mdv, x_max=x_lim_right_mdv,
                      y_min=yb1[0], y_max=yb1[1], x_lab='m/s', y_lab='height (km)', ax='y')

    plot_avg_data_set(Comp_avgH_mdv_plot, 'Doppler Velocity',
                      ds1.t_plt, ds1.heightavg_mdv, ds2.t_plt, ds2.heightavg_mdv,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=[], y_max=[], x_lab='Time (UTC)', y_lab='m/s', ax='y')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Spectral Width LIMRAD-MIRA
    if pts: print('       -   comp. Spectral Width   ', end='', flush=True)

    x_lim_left_sw = min([ds1.timeavg_sw.min(), ds2.timeavg_sw.min()])
    x_lim_right_sw = max([ds1.timeavg_sw.max(), ds2.timeavg_sw.max()])

    plot_avg_data_set(Comp_avgT_sw_plot, 'Spectral Width',
                      ds1.timeavg_sw, ds1.height,
                      ds2.timeavg_sw, ds2.height,
                      label1='LIMRad', marker1='+', label2='MIRA', marker2='o',
                      x_min=x_lim_left_sw, x_max=x_lim_right_sw,
                      y_min=yb1[0], y_max=yb1[1], x_lab='m/s', y_lab='height (km)', ax='y')

    plot_avg_data_set(Comp_avgH_sw_plot, 'Spectral Width',
                      ds1.t_plt, ds1.heightavg_sw,
                      ds2.t_plt, ds2.heightavg_sw,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=[], y_max=[], x_lab='Time (UTC)', y_lab='m/s', ax='y')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file

    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(xb1[0]) + ' (UTC)  to:  ' + str(xb1[1]) + ' (UTC), no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.35)

    return fig, plt


def Plot_Scatter(ds1, ds2):
    interp_meth = 'linear'
    res_interp = 5  # in [sec]
    head_pos = [0.05]
    stat_pos = [0.25, -0.3]


    # create an array with evenly spaced gridsize
    xnew = np.arange(max(ds1.t_unix[0], ds2.t_unix[0]),
                     min(ds1.t_unix[-1], ds2.t_unix[-1]),
                     res_interp)

    # convert the x-axis unix time to readable date time format
    t_plt_new = [datetime.datetime(1970, 1, 1, 0, 0, 0)
                 + datetime.timedelta(seconds=int(xnew[i])) for i in range(len(xnew))]

    # what do you want to interpolate
    LR_y = ds1.heightavg_Ze
    mira_y = ds2.heightavg_Ze

    LR_ynew = interpolate_data(ds1.t_unix, ds1.Ze, xnew, interp_meth)
    mira_ynew = interpolate_data(ds2.t_unix, ds2.Ze, xnew, interp_meth)

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

    xb1 = [ds1.t_plt[0], ds1.t_plt[-1]]
    xb2 = [ds2.t_plt[0], ds2.t_plt[-1]]

    yb1 = [ds1.height[0], ds1.height[-1]]
    yb2 = [ds2.height[0], ds2.height[-1]]

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   Average reflectivity over height domain  ', end='', flush=True)

    h_Ze_plot.set_title(r' \textbf{Mean-Height Reflectivity}')
    plot_avg_data_set(h_Ze_plot, '',
                      ds1.t_plt, ds1.heightavg_Ze, ds2.t_plt, ds2.heightavg_Ze,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1], ymin=yb1[0], y_max=yb1[1],
                      x_lab='Time (UTC)', y_lab='dBZ', ax='n')

    interp_h_Ze_plot.set_title(r' \textbf{Mean-Height Reflectivity}')
    plot_interpol_data_set(interp_h_Ze_plot, '',
                           t_plt_new, LR_ynew, t_plt_new, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=xb1[0], x_max=xb1[1],
                           y_min=yb1[0], y_max=yb1[1], x_lab='Time (UTC)', y_lab='dBZ')

    place_statistics(interp_h_Ze_plot, stat_pos, [mean_diff, cor_coef[0, 1]], 'Ze')

    scatter_Ze.set_title('Scatter Plot of Height-Mean\n Reflectivity ')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 7.0)

    plot_scatter(scatter_Ze, '', LR_ynew, mira_ynew, '*'
                 , x_min=xy_min, x_max=xy_max, y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in dBZ', y_lab='MIRA data in dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    # same for mean doppler velocity

    LR_y = ds1.heightavg_mdv
    mira_y = ds2.heightavg_mdv

    LR_ynew = interpolate_data(ds1.t_unix, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(ds2.t_unix, mira_y, xnew, interp_meth)
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
                      ds1.t_plt, ds1.heightavg_mdv, ds2.t_plt, ds2.heightavg_mdv,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=yb1[0], y_max=yb1[1], x_lab='Time (UTC)', y_lab='m/s', ax='n')

    place_text(interp_h_mdv_plot, [-1.1, 1.2], r'\LARGE{\textbf{Interpolated Dataset}}')

    interp_h_mdv_plot.set_title(r' \textbf{Height-Mean Mean Doppler Velocity}')
    plot_interpol_data_set(interp_h_mdv_plot, '',
                           t_plt_new, LR_ynew, t_plt_new, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=xb1[0], x_max=xb1[1],
                           y_min=yb1[0], y_max=yb1[1], x_lab='Time (UTC)', y_lab='m/s')

    place_statistics(interp_h_mdv_plot, stat_pos, [mean_diff, cor_coef[0, 1]], 'mdv')
    scatter_mdv.set_title('Scatter Plot of Height-Mean\n Mean Doppler Velocity')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.1)

    plot_scatter(scatter_mdv, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max, y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    # same for mean doppler velocity

    LR_y = ds1.heightavg_sw
    mira_y = ds2.heightavg_sw

    LR_ynew   = interpolate_data(ds1.t_unix, LR_y, xnew, interp_meth)
    mira_ynew = interpolate_data(ds2.t_unix, mira_y, xnew, interp_meth)

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
                      ds1.t_plt, ds1.heightavg_sw, ds2.t_plt, ds2.heightavg_sw,
                      label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=yb1[0], y_max=yb1[1], x_lab='Time (UTC)', y_lab='m/s', ax='n')

    interp_h_sw_plot.set_title(r'\textbf{Height-Mean Spectral Width}')
    plot_interpol_data_set(interp_h_sw_plot, '',
                           t_plt_new, LR_ynew, t_plt_new, mira_ynew,
                           label1='LIMRad', marker1='.', label2='MIRA', marker2='.',
                           x_min=xb1[0], x_max=xb1[1],
                           y_min=yb1[0], y_max=yb1[1], x_lab='Time (UTC)', y_lab='m/s')

    place_statistics(interp_h_sw_plot, stat_pos, [mean_diff, cor_coef[0, 1]], 'sw')

    scatter_sw.set_title('Scatter Plot of Height-Mean\n Spectral Width ')

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.03)

    plot_scatter(scatter_sw, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max, y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRad data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file

    first_line = 'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(xb1[0]) + ' (UTC)  to:x  ' + str(xb1[3]) + ' (UTC), no attenuation correction '

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.91])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.55)

    return fig, plt