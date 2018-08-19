import numpy as np

import matplotlib.pyplot as plt
import matplotlib        as mpl

from matplotlib import rc
from matplotlib import colors

rc('font',family='serif')
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from matplotlib              import dates
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable



def plot_data_set(fig,axh,text,x,y,z,vmi,vma,x_min,x_max,y_min,y_max,x_lab,y_lab,z_lab):
    text = r'\textbf{'+text+'}'
    if text.find('sw') > 0 or text.find('Spectral Width')>0:
        cp = axh.pcolormesh(x, y, z,  norm=colors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
    else:
        place_text(axh, [.02, 1.05], text )
        cp = axh.pcolormesh(x, y, z, vmin=vmi, vmax=vma, cmap='jet')
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax0 = divider1.append_axes("right", size="3%", pad=0.1)
    cbar= fig.colorbar(cp, cax=cax0, ax=axh)
    cbar.set_label(z_lab)
    axh.axes.tick_params(axis='both', direction='inout', length=10, width=1.5)
    axh.set_ylabel(y_lab)
    axh.set_xlim( left = x_min, right = x_max )
    axh.set_ylim( bottom = y_min, top = y_max, )

    # exceptions
    if x_lab == '':
        axh.axes.xaxis.set_ticklabels([])
    else:
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
        axh.set_xlabel(x_lab)

def plot_correlation_matrix(axh,text,z,vmi,vma,x_lab,y_lab,z_lab):
    from pylab import pcolor
    #text = r'\textbf{'+text+'}'
    #place_text(axh, [.02, 1.05], text )
    #z = z[1500:2000,1500:2000]
    cp = axh.pcolor(z, norm=colors.LogNorm(vmin=vmi, vmax=vma), cmap='jet')
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax0 = divider1.append_axes("right", size="3%", pad=0.1)
    cbar= fig.colorbar(cp, cax=cax0, ax=axh)
    cbar.set_label(z_lab)
    axh.axes.tick_params(axis='both', direction='inout', length=10, width=1.5)

    # exceptions
    #axh.axes.xaxis.set_ticklabels([]) if x_lab == '' else axh.set_xlabel(x_lab)
    #axh.axes.yaxis.set_ticklabels([]) if y_lab == '' else axh.set_ylabel(y_lab)


    #axh.set_aspect('equal', 'box')

def plot_correlation(axh,text,x,y,label,marker,x_min,x_max,y_min,y_max,x_lab,y_lab):
    place_text(axh, [.15, 1.1], text )
    axh.plot(x,y,marker,label=label)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
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
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))

def plot_avg_data_set(axh,text,x1,y1,x2,y2,label1,marker1,label2,marker2,x_min,x_max,y_min,y_max,x_lab,y_lab,ax):
    text = r'\textbf{' + text + '}'
    place_text(axh, [.02, 1.05], text )
    axh.scatter(x1,y1,marker=marker1,label=label1)
    axh.scatter(x2,y2,marker=marker2,label=label2)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    if ax == 'y':
        divider1 = make_axes_locatable(axh)
        cax = divider1.append_axes("right", size="3%", pad=0.1)
        axh.legend(loc="upper right")

        cax.set_facecolor('none')
        for axis in ['top','bottom','left','right']:
            cax.spines[axis].set_linewidth(0)
        cax.set_xticks([])
        cax.set_yticks([])

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))

def plot_phase_data_set(axh,text,x,y,marker,x_min,x_max,y_min,y_max,x_lab,y_lab):
    text = r'\textbf{' + text + '}'
    place_text(axh, [.02, 1.05], text )
    axh.scatter(x,y,marker=marker)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    divider1 = make_axes_locatable(axh)
    cax = divider1.append_axes("right", size="3%", pad=0.25)
    axh.legend(loc="upper right")

    cax.set_facecolor('none')
    for axis in ['top','bottom','left','right']:
        cax.spines[axis].set_linewidth(0)
    cax.set_xticks([])
    cax.set_yticks([])

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))

def plot_interpol_data_set(axh,text,x1,y1,x2,y2,label1,marker1,label2,marker2,x_min,x_max,y_min,y_max,x_lab,y_lab):
    place_text(axh, [.15, 1.1], text )
    axh.plot(x1,y1,marker1,label=label1)
    axh.plot(x2,y2,marker2,label=label2)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min, right = x_max )
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")

    # exceptions
    if  not (y_min==y_max):
        axh.set_ylim( bottom = y_min, top = y_max )
    if x_lab == 'Time (UTC)':
        axh.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))



def plot_scatter(axh,text,x,y,marker,x_min,x_max,y_min,y_max,x_lab,y_lab):
    import matplotlib.ticker as ticker
    from matplotlib.patches     import Polygon
    from matplotlib.collections import PatchCollection

    #place_text(axh, [.05, 1.1], text )
    axh.plot(x,y,marker)
    axh.set_xlabel( x_lab )
    axh.set_ylabel( y_lab )
    axh.set_xlim( left = x_min,   right = x_max )
    axh.set_ylim( bottom = y_min, top = y_max)
    axh.xaxis.set_ticks((np.arange(x_min, x_max, (x_max-x_min)/4.0)))
    axh.yaxis.set_ticks((np.arange(y_min, y_max, (y_max-y_min)/4.0)))
    axh.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    axh.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    axh.grid(linestyle=':')
    axh.legend(loc="upper right")
    axh.set_aspect('equal', 'box')

    # plot 1:1 line
    N = 100
    X = np.linspace(x_min, x_max, N)
    axh.plot(X, X, 'k--')

    # add patches
    colors = [np.divide([31,119,180],255.) , np.divide([255,127,14],255.)]
    patches = [ Polygon( [[x_min,y_min],[x_max,y_min],[x_max,y_max]], facecolor='C0', fill=True),
                Polygon( [[x_min,y_min],[x_min,y_max],[x_max,y_max]], facecolor='C1', fill=True)]

    p = PatchCollection(patches, alpha=0.25)
    p.set_color(colors)
    axh.add_collection(p)

def get_plot_ybounds(y1,y2,pm):
    import math
    return max(y1.min(),y2.min())-pm, min(y1.max(),y2.max())+pm



def place_text(plot,pos,text):
    plot.text( pos[0], pos[1], text,
               fontweight = 'bold',
               horizontalalignment = 'left',
               transform = plot.transAxes,
               bbox = dict( facecolor = 'white',
                            edgecolor ='black',
                            pad = 5.)
               )

def place_statistics(plot,pos,stat,vn):
    if vn=='Ze':
        text  = r'$\overline{  Z_g^{35} - Z_e^{94} }=$' + \
                 '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'$\rho(Z_e^{94}, Z_g^{35}) = $' + '{:6.2f}'.format(stat[1])
    if vn=='mdv':
        text = r'$\mathrm{mean}_h(\mathrm{mdv}_{\mathrm{lr}} - \mathrm{mdv}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])
    if vn=='sw':
        text = r'$\mathrm{mean}_h(\mathrm{sw}_{\mathrm{lr}} - \mathrm{sw}_{\mathrm{mi}} ) =$' \
               +'{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])

    plot.text( pos[0], pos[1], text, fontweight = 'bold',
               horizontalalignment = 'center',
               transform = plot.transAxes,
               bbox=dict(facecolor='none',
                         edgecolor='black',
                         pad=5.))
    plot.text(pos[0]+0.6, pos[1], text2, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))

def place_statistics_mean_corcoef(plot,pos,stat,vn):
    if vn=='Ze':
        text  = r'$\overline{Z_e^{94} - Z_e^{35} } =$' + \
                 '{:6.2f}'.format(stat[0]) + ' dBZ'
        text2 = r'corr$(Z_e^{94}, Z_e^{35}) = $'+ r'$=$' + '{:6.2f}'.format(stat[1])
    if vn=='mdv':
        text = r'$\mathrm{mean}_h(\mathrm{mdv}_{\mathrm{lr}} - \mathrm{mdv}_{\mathrm{mi}} ) =$' \
               + '{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])
    if vn=='sw':
        text = r'$\mathrm{mean}_h(\mathrm{sw}_{\mathrm{lr}} - \mathrm{sw}_{\mathrm{mi}} ) =$' \
               +'{:6.2f}'.format(stat[0]) + ' m/s'
        text2 = 'corr(lr, mi)'+ r'$=$'  + '{:6.2f}'.format(stat[1])

    plot.text( pos[0], pos[1], text, fontweight = 'bold',
               horizontalalignment = 'center',
               transform = plot.transAxes,
               bbox=dict(facecolor='none',
                         edgecolor='black',
                         pad=5.))
    plot.text(pos[0]+0.6, pos[1], text2, fontweight='bold',
              horizontalalignment='center',
              transform=plot.transAxes,
              bbox=dict(facecolor='none',
                        edgecolor='black',
                        pad=5.))


