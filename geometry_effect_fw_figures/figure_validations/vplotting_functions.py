"""fuctions for plotting cfd run results against ref data"""

import csv
import os

import matplotlib.patches as patches
import matplotlib.path as path
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
                cf_array.append([float(row[0]), float(row[3])])
                line_count += 1

        print(f'Processed {line_count} lines in {cfd_data_file}')

    cf_array = np.array(cf_array)

    return cf_array


def read_ref_data(ref_data_file):
    """read wing geometry data"""
    ref_array = []
    with open(ref_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            ref_array.append([float(row[0]), float(row[1])])
            line_count += 1

        print(f'Processed {line_count} lines in {ref_data_file}')

    ref_array = np.array(ref_array)

    return ref_array


def cf_plotter(data_array, time_scale, legends, time_to_plot, show_range,
               image_out_path):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 15,
        'figure.figsize': (14, 5),
        'lines.linewidth': 2.0,
        'lines.markersize': 4,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.8,
        'figure.subplot.bottom': 0.2,
        'figure.subplot.wspace': 0.25,
        'figure.subplot.hspace': 0.1,
    })

    mcf_mark = ['mcf_ad', 'mcf_sym', 'mcf_dl']
    fig_mark = ['(a)', '(b)', '(c)']
    fig, ax = plt.subplots(1, 3)

    for datai, legendi, axi, mcfname, fmk in zip(data_array, legends, ax,
                                                 mcf_mark, fig_mark):
        # if time_to_plot != 'all':
        # axi.set_xlim(time_to_plot)
        axi.set_xlim([0.0, 1.0])
        if show_range != 'all':
            axi.set_ylim(show_range)

        mcf_arr = []
        for cfdata, lgd in zip(datai, legendi):
            if lgd != 'Current':
                cfdata = np.array([
                    cfdata[:, 0] - np.rint(cfdata[0, 0]) + time_to_plot[0],
                    cfdata[:, 1]
                ])
                cfdata = np.transpose(cfdata)

            #----------------------

            cf_t = cfdata[:, 0] / time_scale - time_to_plot[0]
            axi.plot(cf_t, cfdata[:, 1], label=lgd)

            cf_spl = UnivariateSpline(cf_t, cfdata[:, 1], s=0)
            mcf = cf_spl.integral(0.0, 1.0)
            mcf_arr.append(mcf)

        axi.set_xlabel(r'$\^t$')
        # axi.axvline(x=0.5, color='k', linestyle='-', linewidth=0.5)
        # axi.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
        if axi == ax[0]:
            axi.set_ylabel(r'$C_L$')
        if axi == ax[1]:
            axi.legend(loc='upper center',
                       bbox_to_anchor=(0.5, 1.25),
                       ncol=4,
                       fontsize='small',
                       frameon=False)
        # axi.label_outer()

        marky_loc = axi.get_ylim()[1] + 0.07 * (axi.get_ylim()[1] -
                                                axi.get_ylim()[0])
        markx_loc = axi.get_xlim()[0] - 0.08 * (axi.get_xlim()[1] -
                                                axi.get_xlim()[0])
        axi.annotate(s=fmk,
                     xy=(markx_loc, marky_loc),
                     ha='center',
                     va='center',
                     annotation_clip=False)

        with open(mcfname + '.dat', 'w') as f:
            for item, lgd in zip(mcf_arr, legendi):
                f.write("%s:\n" % lgd)
                f.write("mcf = %s\n" % '{0:.8g}'.format(item))

    out_image_file = os.path.join(image_out_path, 'validation plot.png')
    fig.savefig(out_image_file)
    # plt.show()

    return fig
