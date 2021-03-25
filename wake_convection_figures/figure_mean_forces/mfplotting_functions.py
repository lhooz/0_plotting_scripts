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
    cl_spl = UnivariateSpline(cf_array[:, 0], cf_array[:, 3], s=0)
    cd_spl = UnivariateSpline(cf_array[:, 0], cf_array[:, 1], s=0)

    mcl_s = cl_spl.integral(0.0, 1.0)
    mcl_w = cl_spl.integral(1.0, 2.0)
    # mcl_w = cl_spl(1.2)
    mcd_s = cd_spl.integral(0.0, 1.0)
    mcd_w = cd_spl.integral(1.0, 2.0)

    ratio_l = mcl_w / mcl_s
    ratio_d = mcd_w / mcd_s
    ratio_ld = mcl_w / mcd_w

    mcf_array = [mcl_s, mcl_w, ratio_l, mcd_s, mcd_w, ratio_d, ratio_ld]

    return mcf_array


def cf_plotter(x_data, data_array, marks, legends, x_range, y_range, y_label,
               image_out_path):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (14, 8),
        'lines.linewidth': 0.5,
        'lines.markersize': 12,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 200,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.8,
        'figure.subplot.bottom': 0.2,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    markers = ['o', 'v', 's']
    x_array = np.array(x_data)
    cf_array = np.array(data_array)
    cf_legends = np.array(legends)
    markc = marks

    no_r = 2
    no_c = len(markc)
    no_legend = len(cf_legends)
    no_x = len(x_array)
    fig, ax = plt.subplots(no_r, no_c)
    fig2, ax2 = plt.subplots(no_r, no_c)
    fig3, ax3 = plt.subplots(1, no_c)
    for c in range(no_c):
        axr1i = ax[0][c]
        axr2i = ax[1][c]
        ax2r1i = ax2[0][c]
        ax2r2i = ax2[1][c]
        ax3i = ax3[c]
        axs = [axr1i, axr2i, ax2r1i, ax2r2i, ax3i]
        for lgd in range(no_legend):
            data_no = no_x * (no_legend * c + lgd)
            # print(data_no)
            mcl = []
            ratiol = []
            mcd = []
            ratiod = []
            ratio_ld = []
            for xi in range(no_x):
                mcl.append(cf_array[data_no + xi][1])
                ratiol.append(cf_array[data_no + xi][2])
                mcd.append(cf_array[data_no + xi][4])
                ratiod.append(cf_array[data_no + xi][5])
                ratio_ld.append(cf_array[data_no + xi][6])
                datatoplot = [mcl, mcd, ratiol, ratiod, ratio_ld]

            for i in range(len(datatoplot)):
                axs[i].plot(x_array,
                            datatoplot[i],
                            label=cf_legends[lgd],
                            marker=markers[lgd],
                            linestyle='-.')

                if x_range != 'all':
                    axs[i].set_xlim(x_range)
                if y_range != 'all':
                    axs[i].set_ylim(y_range[i])

                if c == 1 and i in [0, 2, 4]:
                    axs[i].legend(loc='upper center',
                                  bbox_to_anchor=(0.5, 1.2),
                                  ncol=3,
                                  fontsize='small',
                                  frameon=False)

                if lgd == 0:
                    marky_loc = axs[i].get_ylim()[1] + 0.2 * (
                        axs[i].get_ylim()[1] - axs[i].get_ylim()[0])
                    markx_loc = axs[i].get_xlim()[1] + 0.2 * (
                        axs[i].get_xlim()[1] - axs[i].get_xlim()[0])
                    markymid_loc = axs[i].get_ylim()[0] + 0.5 * (
                        axs[i].get_ylim()[1] - axs[i].get_ylim()[0])

                    if i in [0, 2, 4]:
                        axs[i].annotate(s=markc[c],
                                        xy=(6, marky_loc),
                                        ha='center',
                                        va='center',
                                        annotation_clip=False)

                    axs[i].set_ylabel(y_label[i])
                    axs[i].set_xlabel(r'$s/c$')
                    axs[i].label_outer()

                    axs[i].axhline(y=0,
                                   color='k',
                                   linestyle='-.',
                                   linewidth=0.5)

    fig3.set_size_inches(14, 4)
    title = 'mean force coefficients plot'
    title2 = 'mean force ratio plot'
    title3 = 'mean l_to_d plot '
    out_image_file = os.path.join(image_out_path, title + '.png')
    out_image_file2 = os.path.join(image_out_path, title2 + '.png')
    out_image_file3 = os.path.join(image_out_path, title3 + '.png')
    out_files = [out_image_file, out_image_file2, out_image_file3]
    figs = [fig, fig2, fig3]

    for i in range(len(figs)):
        figs[i].savefig(out_files[i])
        figs[i].show()

    return figs
