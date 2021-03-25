"""fuctions for plotting cfd run results against ref data"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline


def read_cfd_data(cfd_data_file):
    """read cfd results force coefficients data"""
    cf_array = []
    with open(cfd_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_count = 0

        for row in csv_reader:
            if line_count <= 14:
                line_count += 1
            else:
                cf_array.append([
                    float(row[0]),
                    float(row[1]),
                    float(row[2]),
                    float(row[3]),
                    float(row[4]),
                    float(row[5]),
                    float(row[6])
                ])
                line_count += 1

        print(f'Processed {line_count} lines in {cfd_data_file}')

    cf_array = np.array(cf_array)
    return cf_array


def cf_plotter(data_array, marks, time_to_plot, show_range, image_out_path,
               cycle_time, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (14, 8),
        'lines.linewidth': 2.0,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.9,
        'figure.subplot.bottom': 0.1,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    legendx = -0.6
    legendy = 1.27

    cf_array = np.array(data_array)
    markre = marks[0]
    markpa = marks[1]
    legends = marks[2]
    range_cl = show_range[0]
    range_cd = show_range[1]

    fig, axs = plt.subplots(3, 3)
    if plot_mode == 'against_t':
        for rei in range(len(markre)):
            ax = axs[:, rei]
            datano = rei * len(markpa) * len(legends)
            for i in range(len(markpa) * len(legends)):
                ax_mk = int(i / 4)
                legendi = np.mod(i, 4)
                ax[ax_mk].plot(cf_array[datano + i][:, 0] / cycle_time,
                               cf_array[datano + i][:, 3],
                               label=legends[legendi])

                if time_to_plot != 'all':
                    ax[ax_mk].set_xlim(time_to_plot)
                if range_cl != 'all':
                    ax[ax_mk].set_ylim(range_cl)

            for ai in range(len(ax)):
                if ai == 0:
                    marky_loc = ax[0].get_ylim()[1] + 0.27 * (
                        ax[0].get_ylim()[1] - ax[0].get_ylim()[0])
                    ax[0].annotate(s=markre[rei],
                                   xy=(1.1, marky_loc),
                                   ha='center',
                                   va='center',
                                   annotation_clip=False)

                ax[ai].axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
                ax[ai].axvline(x=1, color='k', linestyle='-', linewidth=0.5)

                if ax[0] == axs[:, 0][0]:
                    ax[ai].set_ylabel(r'$C_L$')
                else:
                    ax[ai].set_yticklabels([])

                if ai < len(ax) - 1:
                    ax[ai].set_xticklabels([])

                if ax[0] == axs[:, -1][0]:
                    markx_loc = ax[ai].get_xlim()[1] + 0.15 * (
                        ax[ai].get_xlim()[1] - ax[ai].get_xlim()[0])
                    marky_loc = ax[ai].get_ylim()[0] + 0.5 * (
                        ax[ai].get_ylim()[1] - ax[ai].get_ylim()[0])

                    ax[ai].annotate(s=markpa[ai],
                                    xy=(markx_loc, marky_loc),
                                    ha='center',
                                    va='center',
                                    annotation_clip=False)

                    ax[0].legend(loc='upper center',
                                 bbox_to_anchor=(legendx, legendy),
                                 ncol=4,
                                 fontsize='small',
                                 frameon=False)

            texty_loc = ax[-1].get_ylim()[0] + 0.05 * (ax[-1].get_ylim()[1] -
                                                       ax[-1].get_ylim()[0])
            ax[-1].annotate(s='stroke motion',
                            xy=(0.55, texty_loc),
                            ha='center',
                            va='center',
                            annotation_clip=False)
            ax[-1].annotate(s='wake effect',
                            xy=(1.55, texty_loc),
                            ha='center',
                            va='center',
                            annotation_clip=False)

            ax[-1].set_xlabel(r'$\^t$')

    title = 'force coefficients plot pa cl'
    out_image_file = os.path.join(image_out_path, title + '.png')
    fig.savefig(out_image_file)

    fig2, axs2 = plt.subplots(3, 3)
    if plot_mode == 'against_t':
        for rei in range(len(markre)):
            ax2 = axs2[:, rei]
            datano = rei * len(markpa) * len(legends)
            for i in range(len(markpa) * len(legends)):
                ax2_mk = int(i / 4)
                legendi = np.mod(i, 4)
                ax2[ax2_mk].plot(cf_array[datano + i][:, 0] / cycle_time,
                                 cf_array[datano + i][:, 1],
                                 label=legends[legendi])

                if time_to_plot != 'all':
                    ax2[ax2_mk].set_xlim(time_to_plot)
                if range_cd != 'all':
                    ax2[ax2_mk].set_ylim(range_cd)

            for ai in range(len(ax2)):
                if ai == 0:
                    marky_loc = ax2[0].get_ylim()[1] + 0.27 * (
                        ax2[0].get_ylim()[1] - ax2[0].get_ylim()[0])
                    ax2[0].annotate(s=markre[rei],
                                    xy=(1.1, marky_loc),
                                    ha='center',
                                    va='center',
                                    annotation_clip=False)

                ax2[ai].axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
                ax2[ai].axvline(x=1, color='k', linestyle='-', linewidth=0.5)

                if ax2[0] == axs2[:, 0][0]:
                    ax2[ai].set_ylabel(r'$C_D$')
                else:
                    ax2[ai].set_yticklabels([])

                if ai < len(ax2) - 1:
                    ax2[ai].set_xticklabels([])

                if ax2[0] == axs2[:, -1][0]:
                    markx_loc = ax2[ai].get_xlim()[1] + 0.15 * (
                        ax2[ai].get_xlim()[1] - ax2[ai].get_xlim()[0])
                    marky_loc = ax2[ai].get_ylim()[0] + 0.5 * (
                        ax2[ai].get_ylim()[1] - ax2[ai].get_ylim()[0])

                    ax2[ai].annotate(s=markpa[ai],
                                     xy=(markx_loc, marky_loc),
                                     ha='center',
                                     va='center',
                                     annotation_clip=False)

                    ax2[0].legend(loc='upper center',
                                  bbox_to_anchor=(legendx, legendy),
                                  ncol=4,
                                  fontsize='small',
                                  frameon=False)

            texty_loc = ax2[-1].get_ylim()[0] + 0.05 * (ax2[-1].get_ylim()[1] -
                                                        ax2[-1].get_ylim()[0])
            ax2[-1].annotate(s='stroke motion',
                             xy=(0.55, texty_loc),
                             ha='center',
                             va='center',
                             annotation_clip=False)
            ax2[-1].annotate(s='wake effect',
                             xy=(1.55, texty_loc),
                             ha='center',
                             va='center',
                             annotation_clip=False)

            ax2[-1].set_xlabel(r'$\^t$')

    title2 = 'force coefficients plot pa cd'
    out_image_file = os.path.join(image_out_path, title2 + '.png')
    fig2.savefig(out_image_file)
    plt.show()

    return fig, fig2
