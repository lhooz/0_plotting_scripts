"""fuctions for plotting cfd run results against ref data"""

import csv
import os

import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline


def read_profile_data(profile_data_file):
    """read kinematics from file"""
    profile_arr = []
    with open(profile_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count < 1:
                line_count += 1
            else:
                x_datai = row[1]
                # print(row)
                y_datai = row[0]
                profile_arr.append([float(x_datai), float(y_datai)])
                line_count += 1

        print(f'Processed {line_count} lines in {profile_data_file}')

    profile_arr = np.array(profile_arr)
    # print(profile_arr)

    return profile_arr


def single_plotter(ax_to_plot, data, marker):
    """
    single planform plotter
    """
    ax_to_plot.plot(data[:, 0], data[:, 1], linestyle='-', color='k')

    marky_loc = ax_to_plot.get_ylim()[1] + 0.08 * (ax_to_plot.get_ylim()[1] -
                                                   ax_to_plot.get_ylim()[0])
    markx_loc = ax_to_plot.get_xlim()[0] + 0.5 * (ax_to_plot.get_xlim()[1] -
                                                  ax_to_plot.get_xlim()[0])

    ax_to_plot.annotate(s=marker,
                        xy=(markx_loc, marky_loc),
                        ha='center',
                        va='center',
                        annotation_clip=False)


def g_plotter_Joukowsky(data_array, marks, x_range, y_range, image_out_path):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (13, 14),
        'lines.linewidth': 2,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.9,
        'figure.subplot.bottom': 0.15,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.15,
    })
    data_array = np.array(data_array)

    fig, axs = plt.subplots(nrows=4, ncols=2)
    axall = []
    for axr in axs:
        for axc in axr:
            axall.append(axc)

    for axi, data, marki in zip(axall, data_array, marks):
        axi.set_aspect('equal')
        if x_range != 'all':
            axi.set_xlim(x_range)
        if y_range != 'all':
            axi.set_ylim(y_range)

        # axi.axis('off')
        axi.set_xlabel(r'$r/R$')
        axi.set_ylabel(r'$c/R$')
        axi.label_outer()
        axi.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)

        single_plotter(axi, data, marki)

    title = 'planform plot'
    out_image_file = os.path.join(image_out_path, title + '.svg')

    plt.savefig(out_image_file)
    # plt.show()

    return fig
