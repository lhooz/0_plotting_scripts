"""plotting functions for mesh and convergence figures"""

import os
import csv

import matplotlib.image as mpimg
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


def cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
               cycle_time, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 19,
        'figure.figsize': (10, 8),
        'lines.linewidth': 3,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.1,
        'figure.subplot.right': 0.95,
        'figure.subplot.top': 0.95,
        'figure.subplot.bottom': 0.1,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    cf_array = np.array(data_array)
    range_cl = show_range[0]
    range_cd = show_range[1]
    linstyles = ['-', '--', '--', '--']

    fig, axs = plt.subplots(2, 1)
    mcl_arr = []
    mcd_arr = []

    if plot_mode == 'against_t':
        for i in range(len(legends)):
            axs[0].plot(cf_array[i][:, 0] / cycle_time,
                        cf_array[i][:, 3],
                        label=legends[i],
                        linestyle=linstyles[i])
            axs[1].plot(cf_array[i][:, 0] / cycle_time,
                        cf_array[i][:, 1],
                        label=legends[i],
                        linestyle=linstyles[i])

            cl_spl = UnivariateSpline(cf_array[i][:, 0],
                                      cf_array[i][:, 3],
                                      s=0)
            mcl_s = cl_spl.integral(0.0, 1.0)
            mcl_w = cl_spl.integral(1.0, 2.0)
            mcl = cl_spl.integral(0.0, 2.0)

            cd_spl = UnivariateSpline(cf_array[i][:, 0],
                                      cf_array[i][:, 1],
                                      s=0)
            mcd_s = cd_spl.integral(0.0, 1.0)
            mcd_w = cd_spl.integral(1.0, 2.0)
            mcd = cd_spl.integral(0.0, 2.0)

            mcl_arr.append([mcl_s, mcl_w, mcl])
            mcd_arr.append([mcd_s, mcd_w, mcd])

            with open('meancf_convergence.dat', 'w') as f:
                for item, cf_lgd in zip(mcl_arr, legends):
                    f.write("%s:\n" % cf_lgd)
                    f.write("mcl_s = %s, mcl_w = %s, mcl = %s\n" %
                            ('{0:.8g}'.format(item[0]), '{0:.8g}'.format(
                                item[1]), '{0:.8g}'.format(item[2])))
                for item, ref_lgd in zip(mcd_arr, legends):
                    f.write("%s:\n" % ref_lgd)
                    f.write("mcd_s = %s, mcd_w = %s, mcd = %s\n" %
                            ('{0:.8g}'.format(item[0]), '{0:.8g}'.format(
                                item[1]), '{0:.8g}'.format(item[2])))

        if time_to_plot != 'all':
            axs[0].set_xlim(time_to_plot)
            axs[1].set_xlim(time_to_plot)
        if range_cl != 'all':
            axs[0].set_ylim(range_cl)
            axs[1].set_ylim(range_cd)

        for ax in axs:
            ax.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
            ax.axvline(x=1, color='k', linestyle='-', linewidth=0.5)
            ax.set_xlabel(r'$\^t$')
            ax.label_outer()

        axs[0].set_ylabel(r'$C_l$')
        axs[1].set_ylabel(r'$C_d$')

        axs[0].legend(loc='upper right',
                      ncol=1,
                      fontsize='small',
                      frameon=False)

    title = 'convergence_plot_new'
    out_image_file = os.path.join(image_out_path, title + '.svg')
    fig.savefig(out_image_file)
    # plt.show()

    return fig


def mesh_plotter(image_out_path):
    """plotting and annotating domain and mesh images"""
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 19,
        'figure.figsize': (10, 10),
        'lines.linewidth': 4,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'k',
        'figure.dpi': 300,
    })

    circle_bg = plt.Circle((0, 0), 20, color='k', linewidth=4, fill=False)
    circle_overset = plt.Circle((-3, 0),
                                10,
                                color='k',
                                linewidth=4,
                                fill=False)

    LE = [-3 + 0.5 * np.cos(45 * np.pi / 180), 0.5 * np.sin(45 * np.pi / 180)]
    TE = [-3 - 0.5 * np.cos(45 * np.pi / 180), -0.5 * np.sin(45 * np.pi / 180)]
    plate_line = [LE, TE]

    nverts = len(plate_line)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0::2] = path.Path.MOVETO
    plinepath = path.Path(plate_line, codes)
    platepatch = patches.PathPatch(plinepath,
                                   edgecolor='k',
                                   linewidth=2,
                                   linestyle='-')

    #-----------domain plot------------
    gs_d = dict(left=0.125,
                right=0.9,
                top=0.9,
                bottom=0.2,
                wspace=0.1,
                hspace=0.0)

    fig, axd = plt.subplots(nrows=1, ncols=1, gridspec_kw=gs_d)
    axd.spines['right'].set_visible(False)
    axd.spines['top'].set_visible(False)
    # make arrows
    axd.plot((1), (0),
             ls="",
             marker=">",
             ms=8,
             color="k",
             transform=axd.transAxes,
             clip_on=False)
    axd.plot((0), (1),
             ls="",
             marker="^",
             ms=8,
             color="k",
             transform=axd.transAxes,
             clip_on=False)

    axd.add_artist(circle_bg)
    axd.add_artist(circle_overset)
    axd.add_patch(platepatch)

    axd.annotate(s='Flat plate',
                 xy=(-5.5, -1.5),
                 ha='center',
                 va='center',
                 annotation_clip=False)
    # axd.annotate(s='Domain boundary',
                 # xy=(8.2, -8.0),
                 # ha='center',
                 # va='center',
                 # annotation_clip=False)
    axd.annotate(s='Overset domain',
                 xy=(-5.0, 4.5),
                 ha='center',
                 va='center',
                 annotation_clip=False)
    axd.annotate(s='Background domain',
                 xy=(0.0, 14.0),
                 ha='center',
                 va='center',
                 annotation_clip=False)
    axd.annotate(s='Far field',
                 xy=(-16.0, 16.0),
                 ha='center',
                 va='center',
                 annotation_clip=False)

    axd.set_xlim([-22, 22])
    axd.set_ylim([-22, 22])
    axd.set_aspect('equal')
    axd.set_xlabel(r'$x/c$')
    axd.set_ylabel(r'$y/c$')

    oimage_file = os.path.join(image_out_path, 'domain.svg')
    plt.savefig(oimage_file)
    # plt.show()
