"""fuctions for plotting cfd run results against ref data"""

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
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


def cf_plotter(data_array, legends, time_to_plot, coeffs_show_range,
               oimage_file, cycle_time, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (12, 10),
        'lines.linewidth': 1.8,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.85,
        'figure.subplot.top': 0.9,
        'figure.subplot.bottom': 0.125,
        'figure.subplot.wspace': 0.15,
        'figure.subplot.hspace': 0.1,
    })
    cf_array = np.array(data_array)
    legendre = legends[0]
    legendsc = legends[1]

    fig, ax = plt.subplots(3, 2)
    if plot_mode == 'against_t':
        rowno = 0
        for axrow in ax:
            datano = rowno * len(legendsc)
            for i in range(len(legendsc)):
                axrow[0].plot(cf_array[datano + i][:, 0] / cycle_time,
                              cf_array[datano + i][:, 3],
                              label=legendsc[i])
                axrow[1].plot(cf_array[datano + i][:, 0] / cycle_time,
                              cf_array[datano + i][:, 1],
                              label=legendsc[i])

            if time_to_plot != 'all':
                axrow[0].set_xlim(time_to_plot)
                axrow[1].set_xlim(time_to_plot)
            if coeffs_show_range != 'all':
                axrow[0].set_ylim(coeffs_show_range)
                axrow[1].set_ylim(coeffs_show_range)

            texty_loc1 = axrow[0].get_ylim()[0] + 0.05 * (
                axrow[0].get_ylim()[1] - axrow[0].get_ylim()[0])
            texty_loc2 = axrow[1].get_ylim()[0] + 0.05 * (
                axrow[1].get_ylim()[1] - axrow[1].get_ylim()[0])
            text_loc = [texty_loc1, texty_loc2]

            if axrow[0] == ax[-1][0]:
                for axc, text_loci in zip(axrow, text_loc):
                    if 'pa0' in oimage_file:
                        axc.annotate(s='stroke motion',
                                     xy=(0.55, text_loci),
                                     ha='center',
                                     va='center',
                                     annotation_clip=False)
                        axc.annotate(s='wake effect',
                                     xy=(1.55, text_loci),
                                     ha='center',
                                     va='center',
                                     annotation_clip=False)

                    axc.set_xlabel(r'$\^t$')
            else:
                axrow[0].set_xticklabels([])
                axrow[1].set_xticklabels([])

            axrow[1].set_yticklabels([])
            axrow[0].set_ylabel(r'$C_L$')
            axrow[1].set_ylabel(r'$C_D$')
            if axrow[0] == ax[0][0]:
                axrow[1].legend(loc='upper center',
                                bbox_to_anchor=(-0.05, 1.2),
                                ncol=3,
                                fontsize='small',
                                frameon=False)

            markx_loc = axrow[1].get_xlim()[1] + 0.18 * (
                axrow[1].get_xlim()[1] - axrow[1].get_xlim()[0])
            markymid_loc = axrow[1].get_ylim()[0] + 0.5 * (
                axrow[1].get_ylim()[1] - axrow[1].get_ylim()[0])

            axrow[1].annotate(s=legendre[rowno],
                              xy=(markx_loc, markymid_loc),
                              ha='center',
                              va='center',
                              annotation_clip=False)

            axrow[0].axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
            axrow[1].axhline(y=0, color='k', linestyle='-.', linewidth=0.5)

            axrow[0].axvline(x=1, color='k', linestyle='-', linewidth=0.5)
            axrow[1].axvline(x=1, color='k', linestyle='-', linewidth=0.5)

            v_lines = [1.2, 1.4, 1.6, 1.8]
            for line in v_lines:
                axrow[0].axvline(x=line,
                                 color='k',
                                 linestyle='-.',
                                 linewidth=0.5,
                                 alpha=0.25)
                axrow[1].axvline(x=line,
                                 color='k',
                                 linestyle='-.',
                                 linewidth=0.5,
                                 alpha=0.25)

            rowno += 1

    plt.savefig(oimage_file)
    plt.show()

    return fig
