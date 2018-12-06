import datetime
import math

import matplotlib        as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import LogFormatter

## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})


rc('text', usetex=True)

from matplotlib import dates
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
from modules.Interpolation_Mod import *

from modules.Parameter_Mod import pts, interp_meth
from modules.Utility_Mod import correlation
import numpy as np


def plot_data_set(fig, axh, text, x, y, z, vmi, vma, x_min, x_max, y_min, y_max, x_lab, y_lab, z_lab, p='l'):
    if text: text = r'\huge{\textbf{' + text + '}}'
    if text.find('sw') > 0 or text.find('Spectral Width') > 0:
        cp = axh.pcolormesh(x, y, z, norm=mcolors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
        if p in ['lr', 'r']:
            divider1 = make_axes_locatable(axh)
            cax3 = divider1.append_axes("right", size="3%", pad=0.1)
            formatter = LogFormatter(10, labelOnlyBase=False)
            cbar = fig.colorbar(cp, cax=cax3, ax=axh, format=formatter, ticks=[0.1, 0.2, 0.5, 1, 2])
            cbar.set_ticklabels([0.1, 0.2, 0.5, 1, 2])
            cbar.set_label(z_lab)
    elif text.find('ldr') > 0 or text.find('Linear Depolarisation Ratio') > 0:
        colors1 = plt.cm.binary(np.linspace(0.5, 0.5, 1))
        colors2 = plt.cm.jet(np.linspace(0, 0, 178))
        colors3 = plt.cm.jet(np.linspace(0, 1, 77))
        colors = np.vstack((colors1, colors2, colors3))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

        place_text(axh, [.02, 1.075], text)
        cp = axh.pcolormesh(x, y, z, vmin=vmi, vmax=vma, cmap=mymap)
        if p in ['lr', 'r']:
            divider1 = make_axes_locatable(axh)
            cax4 = divider1.append_axes("right", size="3%", pad=0.1)
            bounds = np.linspace(-30, 0, 500)
            cbar = fig.colorbar(cp, cax=cax4, ax=axh, boundaries=bounds, ticks=[-30, -25, -20, -15, -10, -5, 0])
            cbar.set_ticklabels([-30, -25, -20, -15, -10, -5, 0])
            cbar.set_label(z_lab)
    else:
        place_text(axh, [.02, 1.075], text)
        cp = axh.pcolormesh(x, y, z, vmin=vmi, vmax=vma, cmap='jet')
        if p in ['lr', 'r']:
            divider1 = make_axes_locatable(axh)
            cax0 = divider1.append_axes("right", size="3%", pad=0.1)
            cbar = fig.colorbar(cp, cax=cax0, ax=axh)
            cbar.set_label(z_lab)
    axh.grid(linestyle=':')
    axh.axes.tick_params(axis='x', direction='inout', length=10, width=1.5)
    if p in ['lr', 'r']: axh.set_yticklabels([])
    axh.set_xlim(left=x_min, right=x_max)
    if p == 'l':
        axh.set_ylabel(y_lab)
        axh.set_ylim(bottom=y_min, top=y_max)
        axh.axes.tick_params(axis='Y', direction='inout', length=10, width=1.5)

    # exceptions
    if x_lab == '':
        axh.axes.xaxis.set_ticklabels([])
    else:
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
        axh.set_xlabel(x_lab)


def plot_correlation_matrix(axh, text, z, vmi, vma, x_lab, y_lab, z_lab):
    # text = r'\textbf{'+text+'}'
    # place_text(axh, [.02, 1.05], text )
    # z = z[1500:2000,1500:2000]
    cp = axh.pcolor(z, norm=mcolors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
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
    axh.plot(x, y, marker)  # , label=label)
    axh.set_xlabel(x_lab)
    axh.set_ylabel(y_lab)
    axh.set_xlim(left=x_min, right=x_max)
    axh.grid(linestyle=':')
    # axh.legend(loc="upper right")

    divider1 = make_axes_locatable(axh)
    cax = divider1.append_axes("right", size="3%", pad=0.1)

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
    # text = r'\textbf{' + text + '}'
    # place_text(axh, [.02, 1.05], text)
    # axh.set_title(text)
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
    color1 = [np.divide([31, 119, 180], 255.), np.divide([255, 127, 14], 255.)]
    patches = [Polygon([[x_min, y_min], [x_max, y_min], [x_max, y_max]], facecolor='C0', fill=True),
               Polygon([[x_min, y_min], [x_min, y_max], [x_max, y_max]], facecolor='C1', fill=True)]

    p = PatchCollection(patches, alpha=0.25)
    p.set_color(color1)
    axh.add_collection(p)


def get_plot_ybounds(y1, y2, pm, n):
    bnd = round(min(y1.min(), y2.min()) - pm, n), math.ceil((max(y1.max(), y2.max()) + pm) * 10 ** n) / 10 ** n

    return bnd


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
    text = ''
    text2 = ''
    if vn == 'Ze':
        text = r'$\overline{  Z_g^{35} - Z_e^{94} }=$' + \
               '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'$\rho(Z_e^{94}, Z_g^{35}) = $' + '{:6.2f}'.format(stat[1])
    if vn == 'mdv':
        text = r'$\overline{v_m^{94} - v_m^{35} } =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = r'$\rho(v_m^{94},v_m^{35})=$' + '{:6.2f}'.format(stat[1])
    if vn == 'sw':
        text = r'$\overline{\mathrm{sw}^{94}  - \mathrm{sw}^{35}  } =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = r'$\rho(sw_m^{94},sw_m^{35})=$' + '{:6.2f}'.format(stat[1])

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
    text = ''
    text2 = ''
    if vn == 'Ze':
        text = r'$\overline{Z_e^{94} - Z_e^{35} } =$' + \
               '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'corr$(Z_e^{94}, Z_e^{35}) = $' + r'$=$' + '{:6.2f}'.format(stat[1])
    if vn == 'mdv':
        text = r'$\overline{\mathrm{mdv}^{94} - \mathrm{mdv}^{35} } =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)' + r'$=$' + '{:6.2f}'.format(stat[1])
    if vn == 'sw':
        text = r'$\overline{\mathrm{sw}^{94}  - \mathrm{sw}^{35}  } =$' \
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


def Plot_for_poster(ds):
    import matplotlib.style
    mpl.style.use('classic')

    rc('font', size=16)

    fig = plt.figure(figsize=(10, 10))

    plt_Ze = plt.subplot2grid((2, 1), (0, 0))
    plt_vm = plt.subplot2grid((2, 1), (1, 0))  # , rowspan=2)

    x_label = 'Time (UTC)'
    y_label = 'Height (km)'
    z_label = 'Reflectivity (dBZ)'

    xb = [ds.t_plt[0], ds.t_plt[-1]]

    yb = [ds.height[0], ds.height[-1]]

    plt_Ze.set_title('LIMRAD94, Leipzig, Germany', size=20)
    plot_data_set(fig, plt_Ze, '',
                  ds.t_plt, ds.height, ds.Ze, vmi=-50, vma=20,
                  x_min=xb[0], x_max=xb[1], y_min=yb[0], y_max=yb[1],
                  x_lab='', y_lab=y_label, z_lab=z_label)
    #    plt_Ze.tick_params(
    #        axis='x',  # changes apply to the x-axis
    #        which='both',  # both major and minor ticks are affected
    #        bottom=False,  # ticks along the bottom edge are off
    #        top=False,  # ticks along the top edge are off
    #        labelbottom=False)
    #    plt_Ze.tick_params(
    #        axis='y',  # changes apply to the x-axis
    #        which='both',  # both major and minor ticks are affected
    #        left=True,  # ticks along the bottom edge are off
    #        right=False,  # ticks along the top edge are off
    #        labelbottom=False)

    z_label = 'mean Doppler velocity (m/s)'

    plot_data_set(fig, plt_vm, '',
                  ds.t_plt, ds.height, ds.mdv, vmi=-4, vma=2,
                  x_min=xb[0], x_max=xb[1], y_min=yb[0], y_max=yb[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label)
    #    plt_vm.tick_params(
    #        axis='x',  # changes apply to the x-axis
    #        which='both',  # both major and minor ticks are affected
    #        bottom=True,  # ticks along the bottom edge are off
    #        top=True,  # ticks along the top edge are off
    #        labelbottom=True)
    #    plt_vm.tick_params(
    #        axis='y',  # changes apply to the x-axis
    #        which='both',  # both major and minor ticks are affected
    #        left=True,  # ticks along the bottom edge are off
    #        right=False,  # ticks along the top edge are off
    #        labelbottom=False)

    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.subplots_adjust(hspace=0.01)

    return fig, plt


def Plot_Radar_Results(ds1, ds2):
    ### plot ###
    print('    Generate subplots:\n')

    # create figure
    font = FontProperties()

    fig = plt.figure(figsize=(16, 10))

    LR_Ze_plot = plt.subplot2grid((4, 2), (0, 0))
    LR_mdv_plot = plt.subplot2grid((4, 2), (1, 0))  # , rowspan=2)
    LR_sw_plot = plt.subplot2grid((4, 2), (2, 0))  # , rowspan=2)
    LR_ldr_plot = plt.subplot2grid((4, 2), (3, 0))  # , rowspan=2)
    mira_Zg_plot = plt.subplot2grid((4, 2), (0, 1))  # , colspan=2)
    mira_VELg_plot = plt.subplot2grid((4, 2), (1, 1))  # , rowspan=2)
    mira_RMSg_plot = plt.subplot2grid((4, 2), (2, 1))  # , rowspan=2)
    mira_ldr_plot = plt.subplot2grid((4, 2), (3, 1))  # , rowspan=2)

    xb1 = [ds1.t_plt[0], ds1.t_plt[-1]]
    xb2 = [ds2.t_plt[0], ds2.t_plt[-1]]

    yb1 = [ds1.height[0], ds1.height[-1]]
    yb2 = [ds2.height[0], ds2.height[-1]]

    ########################################################################################################
    ########################################################################################################
    # LR_Zelectivity plot
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    x_label = r'\textbf{Time [UTC]}'
    y_label = r'\textbf{Height [km]}'
    z_label = r'\textbf{Reflectivity [dBZ]}'

    LR_Ze_plot.set_title(r'\large{\textbf{LIMRAD 94GHz Radar}}')
    plot_data_set(fig, LR_Ze_plot, '',
                  ds1.t_plt, ds1.height, ds1.Ze, vmi=-50, vma=20,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='l')

    mira_Zg_plot.set_title(r'\large{\textbf{MIRA 35GHz Radar}}')
    plot_data_set(fig, mira_Zg_plot, '',
                  ds2.t_plt, ds2.height, ds2.Ze, vmi=-50, vma=20,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab=y_label, z_lab=z_label)

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # mean doppler velocity plot
    if pts: print('       -   Mean Doppler velocity   ', end='', flush=True)

    z_label = r'\textbf{Mean Doppler}' + '\n' + r'\textbf{Velocity [m/s]}'

    # LIMRAD mean Doppler velocity
    plot_data_set(fig, LR_mdv_plot, '',
                  ds1.t_plt, ds1.height, ds1.mdv, vmi=-4, vma=2,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='l')

    # MIRA mean Doppler velocity
    plot_data_set(fig, mira_VELg_plot, '',
                  ds2.t_plt, ds2.height, ds2.mdv, vmi=-4, vma=2,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab=y_label, z_lab=z_label)

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # spectral width plot
    if pts: print('       -   Spectral Width   ', end='', flush=True)

    z_label = r'\textbf{Spectral Width [m/s]}'
    # LIMRAD spectral width
    plot_data_set(fig, LR_sw_plot, 'sw',
                  ds1.t_plt, ds1.height, ds1.sw, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='l')

    # MIRA spectral width
    plot_data_set(fig, mira_RMSg_plot, 'sw',
                  ds2.t_plt, ds2.height, ds2.sw, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab=y_label, z_lab=z_label)

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # linear depolarization ratio
    if pts: print('       -   Linear Depolarization Ratio  ', end='', flush=True)

    z_label = r'\textbf{Linear Depolarization}' + '\n' + r'\textbf{Ratio [dB]}'

    # LIMRAD linear depolarization ratio
    plot_data_set(fig, LR_ldr_plot, '',
                  ds1.t_plt, ds1.height, ds1.ldr, vmi=-30, vma=0,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label, p='l')

    # MIRA linear depolarization ratio
    plot_data_set(fig, mira_ldr_plot, '',
                  ds2.t_plt, ds2.height, ds2.ldr, vmi=-30, vma=0,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label)

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen

    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(xb1[0]) + ' (UTC)  to:  ' + str(xb1[1]) + ' (UTC),'
    third_line = r'using: LIMRAD94 and MIRA35 data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.90])
    plt.subplots_adjust(hspace=0.025, wspace=0.0075)

    return fig, plt


def Plot_CalcMoments_minus_GivenMoments(ds1, mom='Ze'):
    ### plot ###
    if pts: print('    Generate subplots:\n')

    # create figure

    fig = plt.figure(figsize=(16, 10))

    diff = plt.subplot2grid((1, 1), (0, 0))

    xb1 = [ds1.t_plt[0], ds1.t_plt[-1]]

    yb1 = [ds1.height_all[0], ds1.height_all[-1]]
    if mom == 'Ze':
        differ = 10 * np.log10(ds1.diffZe);
        vmin = np.min(differ);
        vmax = np.max(differ);
        z_label = r'\textbf{Reflectivity [dBZ]}'
        # differ = ds1.diffZe; vmin=0.0; vmax=np.max(ds1.diffZe); z_label=r'\textbf{[mm$^6$/m$^3$]}'
    elif mom == 'mdv':
        differ = ds1.diffmdv;
        vmin = np.min(ds1.diffmdv);
        vmax = np.max(ds1.diffmdv);
        z_label = r'\textbf{[m/s]}'
    elif mom == 'sw':
        differ = ds1.diffsw;
        vmin = 0.0;
        vmax = np.max(ds1.diffsw);
        z_label = r'\textbf{[m/s]}'

    ########################################################################################################
    ########################################################################################################
    # LR_Zelectivity plot
    if pts: print('       -  ', mom, '  ', end='', flush=True)

    x_label = r'\textbf{Time [UTC]}'
    y_label = r'\textbf{Height [km]}'

    diff.set_title(r'\large{\textbf{LIMRAD 94GHz Radar NoiseFac0 Lv1 (with noise)}}')
    plot_data_set(fig, diff, '',
                  ds1.t_plt, ds1.height_all, differ, vmi=vmin, vma=vmax,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label, p='lr')

    diff.set_ylabel(y_label)
    diff.set_ylim(bottom=yb1[0], top=yb1[1])
    # diff.axes.tick_params(axis='Y', direction='inout', length=10, width=1.5)

    first_line = r'Difference calculated moments from LV0 and given LV1 moments, Leipzig, Germany: ' + mom
    second_line = r'from: ' + str(xb1[0]) + ' (UTC)  to:  ' + str(xb1[1]) + ' (UTC),'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.90])
    plt.subplots_adjust(hspace=0.025, wspace=0.0075)

    return fig, plt


def Plot_Compare_NoiseFac0(ds1, ds2):
    ### plot ###
    if pts: print('    Generate subplots:\n')

    # create figure

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

    x_label = r'\textbf{Time [UTC]}'
    y_label = r'\textbf{Height [km]}'
    z_label = r'\textbf{Reflectivity [dBZ]}'

    LR_Ze_plot.set_title(r'\large{\textbf{LIMRAD 94GHz Radar NoiseFac0 Lv1 (with noise)}}')
    plot_data_set(fig, LR_Ze_plot, '',
                  ds1.t_plt, ds1.height, ds1.Ze, vmi=-50, vma=20,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='l')

    mira_Zg_plot.set_title(r'\large{\textbf{LIMRAD 94GHz Radar moments from spectra Lv0 (with noise)}}')
    plot_data_set(fig, mira_Zg_plot, '',
                  ds2.t_plt, ds2.height_all, ds2.Ze, vmi=-50, vma=20,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='r')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # mean doppler velocity plot
    if pts: print('       -   Mean Doppler velocity   ', end='', flush=True)

    z_label = r'\textbf{Mean Doppler}' + '\n' + r'\textbf{Velocity [m/s]}'

    # LIMRAD mean Doppler velocity
    plot_data_set(fig, LR_mdv_plot, '',
                  ds1.t_plt, ds1.height, ds1.mdv, vmi=-4, vma=2,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='l')

    # MIRA mean Doppler velocity
    plot_data_set(fig, mira_VELg_plot, '',
                  ds2.t_plt, ds2.height_all, ds2.mdv, vmi=-4, vma=2,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab='', y_lab=y_label, z_lab=z_label, p='r')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ########################################################################################################
    # spectral width plot
    if pts: print('       -   Spectral Width   ', end='', flush=True)

    z_label = r'\textbf{Spectral Width [m/s]}'
    # LIMRAD spectral width
    plot_data_set(fig, LR_sw_plot, 'sw',
                  ds1.t_plt, ds1.height, ds1.sw, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=xb1[0], x_max=xb1[1], y_min=yb1[0], y_max=yb1[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label, p='l')

    # MIRA spectral width
    plot_data_set(fig, mira_RMSg_plot, 'sw',
                  ds2.t_plt, ds2.height_all, ds2.sw, vmi=10 ** (-1.5), vma=10 ** 0.5,
                  x_min=xb2[0], x_max=xb2[1], y_min=yb2[0], y_max=yb2[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label, p='r')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(xb1[0]) + ' (UTC)  to:  ' + str(xb1[1]) + ' (UTC),'
    third_line = r'using: LIMRAD94 and MIRA35 data;  no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.90])
    plt.subplots_adjust(hspace=0.025, wspace=0.0075)

    return fig, plt


def Plot_Comparison(ds1, ds2):
    ########################################################################
    ### plot comparison ###

    print('    Generate subplots:\n')

    # create figure
    fig, ((LR_Ze_plot, mira_Ze_plot,),
          (Comp_avgT_Ze_plot, Comp_avgH_Ze_plot,),
          (Comp_avgT_mdv_plot, Comp_avgH_mdv_plot,),
          (Comp_avgT_sw_plot, Comp_avgH_sw_plot)) = plt.subplots(4, 2, figsize=(16, 10))

    # calculate x axis limits (same for both time series)
    tb = [ds1.t_plt[0], ds1.t_plt[-1]]
    hb = [ds1.height[0], ds1.height[-1]]

    ################################################################################################################
    ################################################################################################################
    #
    # PLOT LIMRAD and MIRA Reflectivity (Ze) DATASET
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    # LIMRAD reflectivity
    LR_Ze_plot.set_title(r'\textbf{LIMRAD94 - Reflectivity}')

    plot_data_set(fig, LR_Ze_plot, '',
                  ds1.t_plt, ds1.height, ds1.Ze, vmi=-50, vma=20,
                  x_min=tb[0], x_max=tb[1], y_min=hb[0], y_max=hb[1],
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    # MIRA reflectivity
    mira_Ze_plot.set_title(r'\textbf{MIRA35 - Reflectivity}')

    plot_data_set(fig, mira_Ze_plot, '',
                  ds2.t_plt, ds2.height, ds2.Ze, vmi=-50, vma=20,
                  x_min=tb[0], x_max=tb[1], y_min=hb[0], y_max=hb[1],
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Radar Reflectivities LIMRAD-MIRA
    if pts: print('       -   comp. Radar Reflectivity Factor   ', end='', flush=True)

    y1, y2 = get_plot_ybounds(ds1.timeavg_Ze, ds2.timeavg_Ze, 5, 0)
    y3, y4 = get_plot_ybounds(ds1.heightavg_Ze, ds2.heightavg_Ze, 5, 0)
    yb = [min(y1, y3), max(y2, y4)]

    Comp_avgT_Ze_plot.set_title(r'\textbf{Averaged over time}' + '\n' + r'\textbf{Reflectivity}')
    plot_avg_data_set(Comp_avgT_Ze_plot, '',
                      ds1.height, ds1.timeavg_Ze, ds2.height, ds2.timeavg_Ze,
                      label1='LIMRAD', marker1='+', label2='MIRA', marker2='o',
                      x_min=hb[0], x_max=hb[1],
                      y_min=yb[0], y_max=yb[1],
                      x_lab='Height [km]', y_lab='Reflectivity [dBZ]', ax='y')

    #    plot_avg_data_set(Comp_avgT_Ze_plot, '',
    #                      ds1.timeavg_Ze, ds1.height, ds2.timeavg_Ze, ds2.height,
    #                      label1='LIMRAD', marker1='+', label2='MIRA', marker2='o',
    #                      x_min=x_lim_left_Ze, x_max=x_lim_right_Ze,
    #                      y_min=yb1[0], y_max=yb1[1], x_lab='dBZ', y_lab='height (km)', ax='y')

    Comp_avgH_Ze_plot.set_title(r'\textbf{Averaged over range}' + '\n' + r'\textbf{Reflectivity}')
    plot_avg_data_set(Comp_avgH_Ze_plot, '',
                      ds1.t_plt, ds1.heightavg_Ze, ds2.t_plt, ds2.heightavg_Ze,
                      label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                      x_min=tb[0], x_max=tb[1],
                      y_min=yb[0], y_max=yb[1], x_lab='Time (UTC)', y_lab='Reflectivity [dBZ]', ax='y')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   comp. Mean Doppler Velocity   ', end='', flush=True)

    y1, y2 = get_plot_ybounds(ds1.timeavg_mdv, ds2.timeavg_mdv, 0.01, 10)
    y3, y4 = get_plot_ybounds(ds1.heightavg_mdv, ds2.heightavg_mdv, 0.01, 10)
    yb = [min(y1, y3), max(y2, y4)]

    Comp_avgT_mdv_plot.set_title(r'\textbf{ Mean Doppler Velocity}')

    plot_avg_data_set(Comp_avgT_mdv_plot, '',
                      ds1.height, ds1.timeavg_mdv, ds2.height, ds2.timeavg_mdv,
                      label1='LIMRAD', marker1='+', label2='MIRA', marker2='o',
                      x_min=hb[0], x_max=hb[1],
                      y_min=yb[0], y_max=yb[1],
                      x_lab='Height [km]', y_lab='Dopplar Velocity [m/s]', ax='y')

    #    plot_avg_data_set(Comp_avgT_mdv_plot, '',
    #                      ds1.timeavg_mdv, ds1.height, ds2.timeavg_mdv, ds2.height,
    #                      label1='LIMRAD', marker1='+', label2='MIRA', marker2='o',
    #                      x_min=x_lim_left_mdv, x_max=x_lim_right_mdv,
    #                      y_min=yb1[0], y_max=yb1[1], x_lab='m/s', y_lab='height (km)', ax='y')

    Comp_avgH_mdv_plot.set_title(r'\textbf{ Mean Doppler Velocity}')
    plot_avg_data_set(Comp_avgH_mdv_plot, '',
                      ds1.t_plt, ds1.heightavg_mdv, ds2.t_plt, ds2.heightavg_mdv,
                      label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                      x_min=tb[0], x_max=tb[1],
                      y_min=yb[0], y_max=yb[1], x_lab='Time (UTC)', y_lab='Doppler Velocity [m/s]', ax='y')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    ################################################################################################################
    #
    # comparsion of Spectral Width LIMRAD-MIRA
    if pts: print('       -   comp. Spectral Width   ', end='', flush=True)

    y1, y2 = get_plot_ybounds(ds1.timeavg_sw, ds2.timeavg_sw, 0.1, 10)
    y3, y4 = get_plot_ybounds(ds1.heightavg_sw, ds2.heightavg_sw, 0.1, 10)
    yb = [min(y1, y3), max(y2, y4)]

    Comp_avgT_sw_plot.set_title(r'\textbf{Spectral Width}')
    plot_avg_data_set(Comp_avgT_sw_plot, '',
                      ds1.height, ds1.timeavg_sw,
                      ds2.height, ds2.timeavg_sw,
                      label1='LIMRAD', marker1='+', label2='MIRA', marker2='o',
                      x_min=hb[0], x_max=hb[1],
                      y_min=yb[0], y_max=yb[1],
                      x_lab='Height [km]', y_lab='Spectral Width [m/s]', ax='y')

    #   plot_avg_data_set(Comp_avgT_sw_plot, '',
    #                     ds1.timeavg_sw, ds1.height,
    #                     ds2.timeavg_sw, ds2.height,
    #                     label1='LIMRAD', marker1='+', label2='MIRA', marker2='o',
    #                     x_min=x_lim_left_sw, x_max=x_lim_right_sw,
    #                     y_min=yb1[0], y_max=yb1[1], x_lab='m/s', y_lab='height (km)', ax='y')

    Comp_avgH_sw_plot.set_title(r'\textbf{Spectral Width}')
    plot_avg_data_set(Comp_avgH_sw_plot, '',
                      ds1.t_plt, ds1.heightavg_sw,
                      ds2.t_plt, ds2.heightavg_sw,
                      label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                      x_min=tb[0], x_max=tb[1],
                      y_min=yb[0], y_max=yb[1], x_lab='Time (UTC)', y_lab='Spectral Width [m/s]', ax='y')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # Save figure to file

    first_line = r'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = r'from: ' + str(tb[0]) + ' (UTC)  to:  ' + str(tb[1]) + ' (UTC), no attenuation correction'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}'
    plt.suptitle(file_name)

    plt.tight_layout(rect=[0, 0.01, 1, 0.93])
    plt.subplots_adjust(hspace=0.6)

    return fig, plt


def Plot_Scatter(ds1, ds2):
    interp_meth = 'linear'
    res_interp = 10  # in [sec]
    stat_pos = [0.2, -0.35]

    # create an array with evenly spaced gridsize
    xnew = np.arange(max(ds1.t_unix[0], ds2.t_unix[0]),
                     min(ds1.t_unix[-1], ds2.t_unix[-1]),
                     res_interp)

    # convert the x-axis unix time to readable date time format
    t_plt_new = [datetime.datetime(1970, 1, 1, 0, 0, 0)
                 + datetime.timedelta(seconds=int(xnew[i])) for i in range(len(xnew))]

    print('    Generate subplots:\n')

    fig = plt.figure(figsize=(16, 10))

    # fig, axis = plt.subplots(3, 3, figsize=(16, 10))

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

    ################################################################################################################
    #

    LR_ynew = interpolate_data(ds1.t_unix, ds1.heightavg_Ze, xnew, interp_meth)
    mira_ynew = interpolate_data(ds2.t_unix, ds2.heightavg_Ze, xnew, interp_meth)

    # calculate the mean difference and covariance matrix
    mean_diff_Ze = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef_Ze = np.corrcoef(LR_ynew, mira_ynew)

    # comparsion of Radar reflectivity LIMRAD-MIRA
    if pts: print('       -   Average reflectivity over height domain  ', end='', flush=True)

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 7.0, 0)

    h_Ze_plot.set_title(r' \textbf{Reflectivity}')
    plot_avg_data_set(h_Ze_plot, '',
                      ds1.t_plt, ds1.heightavg_Ze, ds2.t_plt, ds2.heightavg_Ze,
                      label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1], y_min=xy_min, y_max=xy_max,
                      x_lab='Time (UTC)', y_lab='dBZ', ax='n')

    interp_h_Ze_plot.set_title(r' \textbf{Reflectivity}')
    plot_interpol_data_set(interp_h_Ze_plot, '',
                           t_plt_new, LR_ynew, t_plt_new, mira_ynew,
                           label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                           x_min=xb1[0], x_max=xb1[1],
                           y_min=xy_min, y_max=xy_max, x_lab='Time (UTC)', y_lab='dBZ')

    scatter_Ze.set_title(r'\textbf{Scatter Plot of}' + '\n ' + r'\textbf{Reflectivity}')
    plot_scatter(scatter_Ze, '', LR_ynew, mira_ynew, '*'
                 , x_min=xy_min, x_max=xy_max, y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRAD data in dBZ', y_lab='MIRA data in dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    ################################################################################################################
    #

    # same for mean doppler velocity
    LR_ynew = np.ma.masked_equal(interpolate_data(ds1.t_unix, ds1.heightavg_mdv, xnew, interp_meth), 0.0)
    mira_ynew = np.ma.masked_equal(interpolate_data(ds2.t_unix, ds2.heightavg_mdv, xnew, interp_meth), 0.0)

    # calculate the mean difference and covariance matrix
    mean_diff_mdv = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef_mdv = np.corrcoef(LR_ynew, mira_ynew)

    # comparsion of Mean Doppler velocities LIMRAD-MIRA
    if pts: print('       -   Average mean doppler velocity over height domain  ', end='', flush=True)

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.1, 1)

    place_text(h_mdv_plot, [-1.1, 1.2], r'\LARGE{\textbf{original:}}' + ' 5 [sec] resolution')

    h_mdv_plot.set_title(r' \textbf{ Mean Doppler Velocity}')
    plot_avg_data_set(h_mdv_plot, '',
                      ds1.t_plt, ds1.heightavg_mdv, ds2.t_plt, ds2.heightavg_mdv,
                      label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=xy_min, y_max=xy_max, x_lab='Time (UTC)', y_lab='m/s', ax='n')

    place_text(interp_h_mdv_plot, [-1.1, 1.2],
               r'\LARGE{\textbf{interpolated:}}' + ' {0:.1f} [sec] resolution'.format(res_interp))

    interp_h_mdv_plot.set_title(r' \textbf{ Mean Doppler Velocity}')
    plot_interpol_data_set(interp_h_mdv_plot, '',
                           t_plt_new, LR_ynew, t_plt_new, mira_ynew,
                           label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                           x_min=xb1[0], x_max=xb1[1],
                           y_min=xy_min, y_max=xy_max, x_lab='Time (UTC)', y_lab='m/s')

    scatter_mdv.set_title(r'\textbf{Scatter Plot of}' + '\n ' + r'\textbf{Mean Doppler Velocity}')

    plot_scatter(scatter_mdv, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max, y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRAD data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713')  # #print checkmark (✓) on screen)

    ################################################################################################################
    #
    # comparsion of spectral width LIMRAD-MIRA

    LR_ynew = interpolate_data(ds1.t_unix, ds1.heightavg_sw, xnew, interp_meth)
    mira_ynew = interpolate_data(ds2.t_unix, ds2.heightavg_sw, xnew, interp_meth)

    # calculate the mean difference and covariance matrix
    mean_diff_sw = np.mean(np.absolute(LR_ynew - mira_ynew))
    cor_coef_sw = np.corrcoef(LR_ynew, mira_ynew)

    if pts: print('       -   Average spectral width over height domain  ', end='', flush=True)

    xy_min, xy_max = get_plot_ybounds(LR_ynew, mira_ynew, 0.03, 2)

    h_sw_plot.set_title(r'\textbf{ Spectral Width}')
    plot_avg_data_set(h_sw_plot, '',
                      ds1.t_plt, ds1.heightavg_sw, ds2.t_plt, ds2.heightavg_sw,
                      label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                      x_min=xb1[0], x_max=xb1[1],
                      y_min=xy_min, y_max=xy_max, x_lab='Time (UTC)', y_lab='m/s', ax='n')

    interp_h_sw_plot.set_title(r'\textbf{ Spectral Width}')
    plot_interpol_data_set(interp_h_sw_plot, '',
                           t_plt_new, LR_ynew, t_plt_new, mira_ynew,
                           label1='LIMRAD', marker1='.', label2='MIRA', marker2='.',
                           x_min=xb1[0], x_max=xb1[1],
                           y_min=xy_min, y_max=xy_max, x_lab='Time (UTC)', y_lab='m/s')

    scatter_sw.set_title(r'\textbf{Scatter Plot of}' + '\n ' + r'\textbf{Spectral Width} ')

    plot_scatter(scatter_sw, '', LR_ynew, mira_ynew, '*',
                 x_min=xy_min, x_max=xy_max, y_min=xy_min, y_max=xy_max,
                 x_lab='LIMRAD data in m/s', y_lab='MIRA data in m/s')

    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)

    # plot differences, correlation
    place_statistics(interp_h_Ze_plot, stat_pos, [mean_diff_Ze, cor_coef_Ze[0, 1]], 'Ze')
    place_statistics(interp_h_mdv_plot, stat_pos, [mean_diff_mdv, cor_coef_mdv[0, 1]], 'mdv')
    place_statistics(interp_h_sw_plot, stat_pos, [mean_diff_sw, cor_coef_sw[0, 1]], 'sw')

    # Save figure to file

    first_line = 'Comparison of LIMRAD 94GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(xb1[0]) + ' (UTC)  to:x  ' + str(xb1[1]) + ' (UTC), no attenuation correction '

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.91])
    plt.subplots_adjust(hspace=0.65)

    return fig, plt


def Plot_2D_Interpolation(ds1, ds2):
    LRtoMIRA_Ze = Interpolate_2D_neu(ds1, ds1.Ze, ds2, interp_meth)
    MIRAtoLR_Ze = Interpolate_2D_neu(ds2, ds2.Ze, ds1, interp_meth)

    # converting back to linear untis for mean difference calculation
    LR_Ze_mm6m3_I = np.power(np.divide(LRtoMIRA_Ze, 10.0), 10.0)
    LR_Ze_mm6m3 = np.power(np.divide(ds1.Ze, 10.0), 10.0)
    mira_Zg_mm6m3_I = np.power(np.divide(MIRAtoLR_Ze, 10.0), 10.0)
    mira_Zg_mm6m3 = np.power(np.divide(ds2.Ze, 10.0), 10.0)

    # display mean difference in dBZ again
    mean_diff_LRtoM = 10 * np.log10(np.mean(np.abs(mira_Zg_mm6m3 - LR_Ze_mm6m3_I)))
    mean_diff_MtoLR = 10 * np.log10(np.mean(np.abs(mira_Zg_mm6m3_I - LR_Ze_mm6m3)))

    # correlation is calculated with logarithmic units
    correlation_LR_MIRA = correlation(LRtoMIRA_Ze, ds2.Ze)
    correlation_MIRA_LR = correlation(ds1.Ze, MIRAtoLR_Ze)

    hmin, hmax = get_plot_ybounds(ds1.height, ds2.height, 0, 0)

    fig, ((p1, p2), (p3, p4), (p5, p6)) = plt.subplots(nrows=3, ncols=2, figsize=(16, 10))

    ########################################################################################################
    # plot table
    cl = ['', r'time resolution $\delta t$', r'height resolution $\delta h$',
          r'height in [km]', r'mean difference in [dBZ]', r'correlation ']

    str_meanLRM = '{:5.2f}'.format(mean_diff_LRtoM)
    str_meanMLR = '{:5.2f}'.format(mean_diff_MtoLR)
    str_corrLRM = '{:5.2f}'.format(np.mean(correlation_LR_MIRA))
    str_corrMLR = '{:5.2f}'.format(np.mean(correlation_MIRA_LR))

    str_rangeres_MIRA = '{:.2f}'.format(ds2.drg)
    str_rangeres_LR = ''
    for ic in range(len(ds1.range_res)):
        if ic + 1 == len(ds1.range_res):
            str_rangeres_LR += '{:.2f}'.format(ds1.range_res[ic])
        else:
            str_rangeres_LR += '{:.2f}/'.format(ds1.range_res[ic])

    table = r'\renewcommand{\arraystretch}{1.25}' \
            r'\begin{tabular}{ | c | c | c | c | c | c |} \hline' \
            r'       & ' + str(cl[1]) + r' & ' + str(cl[2]) + r' & ' + str(cl[3]) + r' & ' + str(cl[4]) + r'  & ' + str(
        cl[5]) + r'  \\ \hline ' \
                 r'MIRA   & $\sim 3$ sec & $' + str_rangeres_MIRA + r'$ m & $' + str(hmin) + r'-' + str(hmax) + \
            r' $ & $\overline{\mathcal{I}(Z_e^{35})-Z_e^{94}}=' + str_meanMLR + \
            r' $ & $\rho(\mathcal{I}(Z_e^{35}),Z_e^{94})=$' + \
            r' $' + str_corrMLR + r'$\\ \hline ' \
                                  r'LIMRAD & $\sim 5$ sec & $' + str_rangeres_LR + r'$ m & $' + str(hmin) + r'-' + str(
        hmax) + \
            r' $ & $\overline{Z_e^{35}-\mathcal{I}(Z_e^{94})}=' + str_meanLRM + \
            r' $ & $\rho(Z_e^{35},\mathcal{I}(Z_e^{94}))= $ ' \
            r' $ ' + str_corrLRM + r'$\\ \hline ' \
                                   r'\end{tabular} '

    fig.text(0.21, 0.8, table, size=12)

    ########################################################################################################
    # LIMRAD Reflectivity plot
    if pts: print('')
    if pts: print('       -   Radar Reflectivity Factor   ', end='', flush=True)

    p1.set_title(r'\textbf{Radar Reflectivity Factor} $Z_{e}^{94}$')
    plot_data_set(fig, p1, '',
                  ds1.t_plt, ds1.height, ds1.Ze, vmi=-50, vma=20,
                  x_min=ds1.t_plt[0], x_max=ds1.t_plt[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    p2.set_title(r'\textbf{Radar Reflectivity Factor} $Z_{g}^{35}$')
    plot_data_set(fig, p2, '',
                  ds2.t_plt, ds2.height, ds2.Ze, vmi=-50, vma=20,
                  x_min=ds2.t_plt[0], x_max=ds2.t_plt[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    if pts: print('       -   synchronized Radar Reflectivity Factor   ', end='', flush=True)
    # MIRA reflectivit

    p3.set_title(r'\textbf{Interpolation of MIRA} $Z_{g}^{35}$ '
                 r'\textbf{onto LIMRAD grid resolution}',
                 multialignment='center')

    plot_data_set(fig, p3, '',
                  ds1.t_plt, ds1.height, MIRAtoLR_Ze, vmi=-50, vma=20,
                  x_min=ds1.t_plt[1], x_max=ds1.t_plt[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='dBZ')

    p4.set_title(r'\textbf{Interpolation of LIMRAD} $Z_{e}^{94}$ '
                 r'\textbf{onto MIRA grid resolution}', multialignment='center')

    plot_data_set(fig, p4, '',
                  ds2.t_plt, ds2.height, LRtoMIRA_Ze, vmi=-50, vma=20,
                  x_min=ds2.t_plt[1], x_max=ds2.t_plt[-1],
                  y_min=hmin, y_max=hmax,
                  x_lab='Time (UTC)', y_lab='Height (km)', z_lab='dBZ')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    if pts: print('       -   correnlation of Radar Reflectivity Factors   ', end='', flush=True)
    # MIRA reflectivit

    p5.set_title(r'\textbf{Correlation of} $\mathcal{I}(Z_{g}^{35})$ \textbf{and} $Z_e^{94}$', multialignment='center')
    plot_correlation(p5, '', ds1.t_plt, correlation_MIRA_LR, 'LIMRAD', '.',
                     x_min=ds1.t_plt[1], x_max=ds1.t_plt[-1],
                     y_min=-1.1, y_max=1.1, x_lab='Time (UTC)', y_lab='correlation')

    p6.set_title(r'\textbf{Correlation of} $\mathcal{I}(Z_{e}^{94})$ \textbf{and} $Z_e^{35}$', multialignment='center')
    plot_correlation(p6, '', ds2.t_plt, correlation_LR_MIRA, 'MIRA', '.',
                     x_min=ds2.t_plt[1], x_max=ds2.t_plt[-1],
                     y_min=-1.1, y_max=1.1, x_lab='Time (UTC)', y_lab='correlation')

    if pts: print('\u2713')  # #print checkmark (✓) on screen

    # Save figure to file
    first_line = 'Comparison of LIMRAD 94 GHz and MIRA 35GHz Radar Data, Leipzig, Germany,'
    second_line = ' from: ' + str(ds1.t_plt[0]) + ' (UTC)  to:  ' + str(ds1.t_plt[3]) + ' (UTC), '
    third_line = 'using: *.LV1.NC and *.mmclx data (unprocessed datasets)'

    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
    plt.suptitle(file_name)  # place in title needs to be adjusted

    plt.tight_layout(rect=[0, 0.01, 1, 0.8])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.65)

    return fig, plt


def Plot_Doppler_Spectra(ds, c, t0, h0, zbound, thresh=0.0, mean=0.0, int_a=0.0, int_b=0.0):
    if thresh == 0.0 and mean == 0.0 and int_b == 0:
        plot_boundaries = False
    else:
        plot_boundaries = True

    # convert from linear units to logarithmic units
    doppler_spec = np.multiply(np.ma.log10(ds.VHSpec[c][t0, h0, :]), 10.0)

    x1, x2 = [ds.DopplerBins[c][0], ds.DopplerBins[c][-1]]

    # plot spectra
    fig, ax = plt.subplots(1, figsize=(10, 4))
    ax.plot(ds.DopplerBins[c], doppler_spec, color='blue', linestyle=':', label='Doppler Spec')

    if plot_boundaries:
        mean = np.multiply(np.ma.log10(mean), 10.0)
        thresh = np.multiply(np.ma.log10(thresh), 10.0)

        # plot mean noise line and threshold
        ax.plot([x1, x2], [thresh, thresh], color='k', linestyle='-', linewidth=2)
        ax.plot([x1, x2], [mean, mean], color='k', linestyle='--', linewidth=2)

        # plot integration boundaries
        if int_a > -1 and int_b > -1:
            x_0 = ds.DopplerBins[c][int(int_a)]
            x_1 = ds.DopplerBins[c][int(int_b)]
            ax.axvline(x_0, color='k', linestyle='--', linewidth=1)
            ax.axvline(x_1, color='k', linestyle='--', linewidth=1)

    ax.set_xlim(left=x1, right=x2)
    ax.set_ylim(bottom=zbound[0], top=zbound[1])
    ax.set_xlabel('Doppler Velocity (m/s)', fontweight='semibold', fontsize=13)
    ax.set_ylabel('Reflectivity (dBZ)', fontweight='semibold', fontsize=13)
    ax.grid(linestyle=':')
    plt.title("Height: " + str(round(ds.height[c][h0], 2)) + " (km);  Time: "
              + str(ds.t_plt[t0]) + ' (UTC)', fontweight='semibold', fontsize=13)
    ax.legend(fontsize=13)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

    return fig, plt, ax


def Plot_Doppler_Spectra_Wavelet_Transform(ds, vhspec_norm, c, t0, h0, zbound, cwtmatr, widths):
    fontsize = 12

    # convert from linear units to logarithic units
    doppler_spec = np.multiply(np.ma.log10(ds.VHSpec[c][t0, h0, :]), 10.0)
    cwtmatr_spec = cwtmatr

    # cwtmatr_spec = np.multiply(np.ma.log10(cwtmtr), 10.0)

    x1, x2 = [ds.DopplerBins[c][0], ds.DopplerBins[c][-1]]

    # plot spectra
    fig, ax = plt.subplots(3, figsize=(10, 10))

    ax[0].set_title('Doppler spectra, normalized and wavlet transformation, height: '
                    + str(round(ds.height[c][h0], 2)) + ' (km);  time: '
                    + str(ds.t_plt[t0]) + ' (UTC)', fontweight='bold', fontsize=fontsize)

    ax[0].plot(ds.DopplerBins[c], doppler_spec, marker='.', linestyle='-', color='blue', label='Doppler Spec')
    ax[0].set_xlim(left=x1, right=x2)
    ax[0].set_ylim(bottom=-55, top=20)
    ax[0].set_ylabel('Doppler spectrum (dBZ)', fontweight='bold', fontsize=fontsize)
    ax[0].grid(linestyle=':')

    ax[1].plot(ds.DopplerBins[c], vhspec_norm, marker='.', linestyle='-', color='blue', label='normalized Spec')
    ax[1].set_xlim(left=x1, right=x2)
    ax[1].set_ylim(bottom=zbound[0], top=zbound[1])
    ax[1].set_xlabel('Doppler Velocity (m/s)', fontweight='bold', fontsize=fontsize)
    ax[1].set_ylabel('normalized spectrum (-)', fontweight='bold', fontsize=fontsize)
    ax[1].grid(linestyle=':')

    img = ax[2].imshow(cwtmatr_spec, extent=[x1, x2, widths[-1], widths[0]],
                       cmap='gist_stern', aspect='auto', vmin=0.0, vmax=2.0)
    ax[2].set_ylabel('wavelet scale parameter', fontweight='bold', fontsize=fontsize)
    divider = make_axes_locatable(ax[2])
    cax = divider.new_vertical(size="5%", pad=0.5, pack_start=True)
    fig.add_axes(cax)
    cbar = fig.colorbar(img, cax=cax, orientation="horizontal")
    cbar.set_label('Magnitude', fontsize=fontsize)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95], h_pad=0.1)
    #plt.show()

    return fig, plt


def Plot_moment_from_spectra(ds, mom):
    import matplotlib.style
    mpl.style.use('classic')

    rc('font', size=16)

    fig = plt.figure(figsize=(16, 10))

    plt_Ze = plt.subplot2grid((1, 1), (0, 0))

    x_label = 'Time (UTC)'
    y_label = 'Height (km)'
    z_label = 'Reflectivity (dBZ)'

    xb = [ds.t_plt[0], ds.t_plt[-1]]

    yb = [ds.height_all[0], ds.height_all[-1]]

    if mom == 'Ze':
        moment = ds.Ze;
        vmin = -50;
        vmax = 20;
        z_label = '[dBZ]'
    elif mom == 'mdv':
        moment = ds.mdv;
        vmin = -5;
        vmax = 3;
        z_label = '[m/s]'
    elif mom == 'sw':
        moment = ds.sw;
        vmin = 0.0;
        vmax = 4;
        z_label = '[m/s]'

    plt_Ze.set_title('LIMRAD94, Leipzig, Germany', size=20)
    plot_data_set(fig, plt_Ze, '',
                  ds.t_plt, ds.height_all, moment, vmi=vmin, vma=vmax,
                  x_min=xb[0], x_max=xb[1], y_min=yb[0], y_max=yb[1],
                  x_lab=x_label, y_lab=y_label, z_lab=z_label, p='lr')

    plt_Ze.set_ylim(bottom=yb[0], top=yb[1])
    plt_Ze.set_xlabel('Doppler Velocity (m/s)', fontweight='semibold', fontsize=13)
    plt_Ze.set_ylabel('Reflectivity (dBZ)', fontweight='semibold', fontsize=13)
    plt_Ze.grid(linestyle=':')

    place_text(plt_Ze, [0.6, 0.9], 'Number std diviations:' + str(ds.n_std_div))

    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.subplots_adjust(hspace=0.01)

    return fig, plt

# if plot_compare_mira_mmclx:
#
#    # print('')
#    # print('Zh_mira-Ze_mmclx = ',np.sum(miramira_Z-mira_Ze[:,:-1]))
#
#    print('')
#    print('    Generate subplots:\n')
#
#    font = FontProperties()
#
#    fig = plt.figure(figsize=(16, 10))
#
#    Zh_plot = plt.subplot2grid((3, 4), (0, 0))
#    Ze_plot = plt.subplot2grid((3, 4), (0, 1))
#    Zg_plot = plt.subplot2grid((3, 4), (0, 2))  # , rowspan=2)
#    LR_Ze_plot = plt.subplot2grid((3, 4), (0, 3))  # , rowspan=2)
#
#    v_plot = plt.subplot2grid((3, 4), (1, 0))  # , colspan=2)
#    VEL_plot = plt.subplot2grid((3, 4), (1, 1))  # , colspan=2)
#    VELg_plot = plt.subplot2grid((3, 4), (1, 2))  # , rowspan=2)
#    LR_mdv_plot = plt.subplot2grid((3, 4), (1, 3))  # , rowspan=2)
#
#    sw_plot = plt.subplot2grid((3, 4), (2, 0))  # , colspan=2, rowspan=2)
#    RMS_plot = plt.subplot2grid((3, 4), (2, 1))  # , colspan=2, rowspan=2)
#    RMSg_plot = plt.subplot2grid((3, 4), (2, 2))  # , rowspan=2)
#    LR_sw_plot = plt.subplot2grid((3, 4), (2, 3))  # , rowspan=2)
#
#    height_describtion = [-0.05, 1.3]
#    place_text(Zh_plot, height_describtion, r'\textbf{Zh} ... Calibrated reflectivity. Calibration\\'
#                                            ' convention: in the absence of attenuation, \n'
#                                            'a cloud at 273 K containing one million 100-micron\n'
#                                            ' droplets per cubic metre will have a reflectivity \n'
#                                            'of 0 dBZ at all frequencies)\n'
#                                            r'\textbf{v}  ... Radial component of the velocity,\\'
#                                            ' with positive velocities are away from the radar.)\n '
#                                            r'\textbf{width} ... Standard deviation of the reflectivity-\\'
#                                            r'weighted velocities in the radar pulse volume.)')
#    height_describtion = [0.1, 1.3]
#
#    place_text(Ze_plot, height_describtion, r'\textbf{Ze} ... Equivalent Radar Reflectivity\\'
#                                            ' Factor Ze of Hydrometeors  \n'
#                                            r'\textbf{VEL} ... Doppler Velocity VEL\\'
#                                            r'\textbf{RMS} ... Peak Width RMS')
#    place_text(Zg_plot, height_describtion, r'\textbf{Zg} ... Equivalent Radar Reflectivity\\'
#                                            ' Factor Ze of all Targets \n'
#                                            r'\textbf{VELg} ... Doppler Velocity VELg\\'
#                                            r'\textbf{RMSg} ... Peak Width RMSg')
#    place_text(LR_Ze_plot, height_describtion, r'\textbf{Ze} ... Equivalent radar reflectivity factor\\'
#                                               r'\textbf{mdv} ... Mean Doppler Velocity\\'
#                                               r'\textbf{sw} ... Spectrum width')
#    ################################################################################################################
#    #
#    ref_min = -60
#    ref_max = 30
#    # comparsion of Radar Reflectivities LIMRAD-MIRA
#    if pts: print('       -   MIRA35 Reflectivity (Z)  ', end='', flush=True)
#    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]
#
#    Zh_plot.set_title(r'\textbf{Radar Reflectivity Factor (Zh)}')
#    plot_data_set(fig, Zh_plot, '',
#                  UTC_time_mira, mira_height, mira_Ze, vmi=ref_min, vma=ref_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')
#
#    Ze_plot.set_title(r'\textbf{Radar Reflectivity Factor (Ze)}')
#    plot_data_set(fig, Ze_plot, '',
#                  UTC_time_mira, mira_height, mira_Ze, vmi=ref_min, vma=ref_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')
#
#    Zg_plot.set_title(r'\textbf{Radar Reflectivity Factor (Zg)}')
#    plot_data_set(fig, Zg_plot, '',
#                  UTC_time_mira, mira_height, mira_Zg, vmi=ref_min, vma=ref_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')
#
#    LR_Ze_plot.set_title(r'\textbf{Radar Reflectivity Factor (Ze)}')
#    plot_data_set(fig, LR_Ze_plot, '',
#                  UTC_time_LR, LR_height, LR_Ze, vmi=ref_min, vma=ref_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='dBZ')
#
#    if pts: print('\u2713')  # #print checkmark (✓) on screen)
#
#    ################################################################################################################
#    #
#    vel_min = -5
#    vel_max = 3
#    # comparsion of Radar Reflectivities LIMRAD-MIRA
#    if pts: print('       -   MIRA35 Mean Doppler Velocity (VEL)  ', end='', flush=True)
#    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]
#
#    v_plot.set_title(r'\textbf{Mean Doppler Velocity (v)}')
#    plot_data_set(fig, v_plot, '',
#                  UTC_time_mira, mira_height, miraNC_VEL, vmi=vel_min, vma=vel_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    VEL_plot.set_title(r'\textbf{Mean Doppler Velocity (VEL)}')
#    plot_data_set(fig, VEL_plot, '',
#                  UTC_time_mira, mira_height, mira_VEL, vmi=vel_min, vma=vel_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    VELg_plot.set_title(r'\textbf{Mean Doppler Velocity (VELg)}')
#    plot_data_set(fig, VELg_plot, '',
#                  UTC_time_mira, mira_height, mira_VELg, vmi=vel_min, vma=vel_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    LR_mdv_plot.set_title(r'\textbf{Mean Doppler Velocity (mdv)}')
#    plot_data_set(fig, LR_mdv_plot, '',
#                  UTC_time_LR, LR_height, LR_mdv, vmi=vel_min, vma=vel_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    if pts: print('\u2713')  # #print checkmark (✓) on screen)
#
#    ################################################################################################################
#    #
#    sw_min = 10 ** (-1.5)
#    sw_max = 10 ** (0.5)
#    # comparsion of Radar Reflectivities LIMRAD-MIRA
#    if pts: print('       -   MIRA35 Spectral Width (RMS)  ', end='', flush=True)
#    x_lim = [UTC_time_mira[0], UTC_time_mira[-1]]
#
#    sw_plot.set_title(r'\textbf{Spectral Width (width)}')
#    plot_data_set(fig, sw_plot, 'sw',
#                  UTC_time_mira, mira_height, miraNC_RMS, vmi=sw_min, vma=sw_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    RMS_plot.set_title(r'\textbf{Spectral Width (RMS)}')
#    plot_data_set(fig, RMS_plot, 'sw',
#                  UTC_time_mira, mira_height, mira_RMS, vmi=sw_min, vma=sw_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    RMSg_plot.set_title(r'\textbf{Spectral Width (RMSg)}')
#    plot_data_set(fig, RMSg_plot, 'sw',
#                  UTC_time_mira, mira_height, mira_RMSg, vmi=sw_min, vma=sw_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    LR_sw_plot.set_title(r'\textbf{Spectral Width (sw)}')
#    plot_data_set(fig, LR_sw_plot, 'sw',
#                  UTC_time_LR, LR_height, LR_sw, vmi=sw_min, vma=sw_max,
#                  x_min=x_lim[0], x_max=x_lim[-1],
#                  y_min=hmin, y_max=hmax,
#                  x_lab='Time (UTC)', y_lab='height (km)', z_lab='m/s')
#
#    if pts: print('\u2713\n')  # #print checkmark (✓) on screen)
#
#    # Save figure to file
#    date_str = str(plotyear) + str(plotmonth).zfill(2) + str(plotday).zfill(2)
#    first_line = 'Comparison of processed and unprocessed MIRA 35GHz Radar Data, Leipzig, Germany,'
#    second_line = ' from: ' + str(time_int[0]) + ' (UTC)  to:  ' + str(time_int[3]) + ' (UTC), '
#    third_line = 'using: *mira.nc data (1st column) and *.mmclx data (2nd - 4th column; no attenuation correction)'
#
#    file_name = r'\textbf{' + first_line + '}\n' + r'\textbf{' + second_line + '}\n' + r'\textbf{' + third_line + '}'
#    plt.suptitle(file_name)  # place in title needs to be adjusted
#
#    plt.tight_layout(rect=[0, 0.01, 1, 0.80])
#    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
#
#    file = date_str + '_MIRA_mmclx_comparison.png'
#    print('    Save Figure to File :: ' + meteo_path + file + '\n')
#    fig.savefig(meteo_path + file, dpi=dpi_val)
#
#    plt.close()
