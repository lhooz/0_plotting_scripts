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


def single_plotter(ax_to_plot, data, legends, marker):
    """
    single planform plotter
    """
    for datai, legendi in zip(data, legends):
        ax_to_plot.plot(datai[:, 0], datai[:, 1], label=legendi, linestyle='-')

    marky_loc = ax_to_plot.get_ylim()[1] + 0.15 * (ax_to_plot.get_ylim()[1] -
                                                   ax_to_plot.get_ylim()[0])
    markx_loc = ax_to_plot.get_xlim()[0] + 0.5 * (ax_to_plot.get_xlim()[1] -
                                                  ax_to_plot.get_xlim()[0])

    ax_to_plot.annotate(text=marker,
                        xy=(markx_loc, marky_loc),
                        ha='center',
                        va='center',
                        annotation_clip=False)


def g_plotter(data_array, legends, marks, x_range, y_range, image_out_path):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 30,
        'figure.figsize': (14, 8),
        'lines.linewidth': 5,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
    })
    data_array = np.array(data_array)

    gs_kw = dict(left=0.125,
                 right=0.8,
                 top=0.9,
                 bottom=0.1,
                 wspace=0.05,
                 hspace=0.0)

    fig, axs = plt.subplots(nrows=3, ncols=1, gridspec_kw=gs_kw)
    ax_all = axs
    for axi, datai, marker in zip(ax_all, data_array, marks):
        axi.set_aspect('equal')
        if x_range != 'all':
            axi.set_xlim(x_range)
        if y_range != 'all':
            axi.set_ylim(y_range)
        axi.axis('off')
        axi.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)

        single_plotter(axi, datai, legends, marker)

        if axi == ax_all[0]:
            axi.legend(loc='upper center',
                       bbox_to_anchor=(1.15, 1.0),
                       ncol=1,
                       fontsize='small',
                       frameon=False)

        if axi == ax_all[0]:
            axi.axvline(x=0, color='k', linestyle='-.', linewidth=0.5)

            axi.annotate(text=r'$r_R$',
                         xy=(0.03, 0.03),
                         ha='center',
                         va='center',
                         annotation_clip=False)

            axi.annotate(text=r'$o$',
                         xy=(-0.012, 0.0),
                         ha='center',
                         va='center',
                         annotation_clip=False)

            axi.annotate(text='',
                         xy=(0.0, 0.015),
                         xytext=(0.06, 0.015),
                         arrowprops=dict(arrowstyle='<-',
                                         linestyle='-.',
                                         facecolor='k',
                                         lw=0.5),
                         annotation_clip=False)

    title = 'geometry plot'
    out_image_file = os.path.join(image_out_path, title + '.svg')

    plt.savefig(out_image_file)
    # plt.show()

    return fig


def g_plotter_ofs(data_array, legends, x_range, y_range, image_out_path):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 14,
        'figure.figsize': (10, 6),
        'lines.linewidth': 2,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 200,
    })
    data_array = np.array(data_array)

    gs_kw = dict(left=0.125,
                 right=0.8,
                 top=0.8,
                 bottom=0.1,
                 wspace=0.05,
                 hspace=0.0)

    fig, axs = plt.subplots(nrows=1, ncols=1, gridspec_kw=gs_kw)
    axi = axs
    axi.set_aspect('equal')
    if x_range != 'all':
        axi.set_xlim(x_range)
    if y_range != 'all':
        axi.set_ylim(y_range)
    axi.axis('off')
    axi.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)

    single_plotter(axi, data_array, legends, 'Baseline wing')

    axi.legend(loc='upper center',
               bbox_to_anchor=(0.5, 1.15),
               ncol=4,
               fontsize='small',
               frameon=False)

    # if axi == ax_all[0]:
    # axi.axvline(x=0, color='k', linestyle='-.', linewidth=0.5)

    # axi.annotate(s=r'$r_R$',
    # xy=(0.03, 0.03),
    # ha='center',
    # va='center',
    # annotation_clip=False)

    # axi.annotate(s=r'$o$',
    # xy=(-0.012, 0.0),
    # ha='center',
    # va='center',
    # annotation_clip=False)

    # axi.annotate(s='',
    # xy=(0.0, 0.015),
    # xytext=(0.06, 0.015),
    # arrowprops=dict(arrowstyle='<-',
    # linestyle='-.',
    # facecolor='k',
    # lw=0.5),
    # annotation_clip=False)

    title = 'offset plot'
    out_image_file = os.path.join(image_out_path, title + '.svg')

    plt.savefig(out_image_file)
    # plt.show()

    return fig
