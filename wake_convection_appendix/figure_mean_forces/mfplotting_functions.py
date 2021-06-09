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
    mcl_w = cl_spl.integral(1.0, 1.2) / 0.2
    # mcl_w = cl_spl(1.2)
    mcd_s = cd_spl.integral(0.0, 1.0)
    mcd_w = cd_spl.integral(1.0, 1.2) / 0.2

    ratio_l = mcl_w / mcl_s
    ratio_d = mcd_w / mcd_s
    ratio_ld = mcl_w / mcd_w

    mcf_array = [mcl_s, mcl_w, ratio_l, mcd_s, mcd_w, ratio_d, ratio_ld]

    return mcf_array


def cf_plotter90(x_data, data_array, marks, legends, x_range, y_range, y_label,
                 image_out_path, figname):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 19,
        'figure.figsize': (10, 4),
        'lines.linewidth': 0.5,
        'lines.markersize': 8,
        # 'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.1,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.85,
        'figure.subplot.bottom': 0.15,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    markers = ['o', 'v', 's']
    colorlgd = ['tab:blue', 'tab:orange']
    x_array = np.array(x_data)
    cf_array = np.array(data_array)
    cf_legends = np.array(legends)
    markr = marks[0]
    markc = marks[1]

    # no_r = len(markr)
    no_r = 1
    no_c = len(markc)
    no_legend = len(cf_legends)
    no_x = len(x_array)
    fig, ax = plt.subplots(no_r, no_c)
    fig2, ax2 = plt.subplots(no_r, no_c)
    fig3, ax3 = plt.subplots(no_r, no_c)
    fig4, ax4 = plt.subplots(no_r, no_c)
    fig5, ax5 = plt.subplots(no_r, no_c)
    # for r in range(no_r):
    for r in range(3):
        for c in range(no_c):
            if no_r > 1:
                axi = ax[r][c]
                ax2i = ax2[r][c]
                ax3i = ax3[r][c]
                ax4i = ax4[r][c]
                ax5i = ax5[r][c]
            else:
                axi = ax[c]
                ax2i = ax2[c]
                ax3i = ax3[c]
                ax4i = ax4[c]
                ax5i = ax5[c]
            axs = [axi, ax2i, ax3i, ax4i, ax5i]
            for lgd in range(no_legend):
                data_no = no_x * (no_c * no_legend * r + no_legend * c + lgd)
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
                    datatoplot = [mcl, ratiol, mcd, ratiod, ratio_ld]

                for i in range(len(datatoplot)):
                    axs[i].plot(x_array,
                                datatoplot[i],
                                label=cf_legends[lgd],
                                color=colorlgd[lgd],
                                marker=markers[r],
                                linestyle='-.')

                    axs[i].set_xticks(np.arange(2, 12, step=2))
                    if x_range != 'all':
                        axs[i].set_xlim(x_range)
                    if y_range != 'all':
                        axs[i].set_ylim(y_range[i])

                    # if r == 0:
                    # if c == 1:
                    # axs[i].legend(loc='upper center',
                    # bbox_to_anchor=(0.5, 1.17),
                    # ncol=3,
                    # fontsize='small',
                    # frameon=False)

                    if lgd == 0:
                        marky_loc = axs[i].get_ylim()[1] + 0.1 * (
                            axs[i].get_ylim()[1] - axs[i].get_ylim()[0])
                        markx_loc = axs[i].get_xlim()[1] + 0.2 * (
                            axs[i].get_xlim()[1] - axs[i].get_xlim()[0])
                        markymid_loc = axs[i].get_ylim()[0] + 0.5 * (
                            axs[i].get_ylim()[1] - axs[i].get_ylim()[0])

                        if r == 0:
                            axs[i].annotate(s=markc[c],
                                            xy=(6, marky_loc),
                                            ha='center',
                                            va='center',
                                            annotation_clip=False)
                        if c == 0:
                            axs[i].set_ylabel(y_label[i])
                        else:
                            axs[i].set_yticklabels([])
                        # if c == no_c - 1:
                        # axs[i].annotate(s=markr[r],
                        # xy=(markx_loc, markymid_loc),
                        # ha='center',
                        # va='center',
                        # annotation_clip=False)
                        if r == no_r - 1:
                            axs[i].set_xlabel(r'$s/c$')
                        # else:
                        # axs[i].set_xticklabels([])

                        axs[i].axhline(y=0,
                                       color='k',
                                       linestyle='-.',
                                       linewidth=0.5)

    title = 'mean lift coefficients plot_ta ' + figname
    title2 = 'mean lift ratio plot ' + figname
    title3 = 'mean drag coefficients plot_ta ' + figname
    title4 = 'mean drag ratio plot ' + figname
    title5 = 'mean wake ld ratio plot ' + figname
    out_image_file = os.path.join(image_out_path, title + '.svg')
    out_image_file2 = os.path.join(image_out_path, title2 + '.svg')
    out_image_file3 = os.path.join(image_out_path, title3 + '.svg')
    out_image_file4 = os.path.join(image_out_path, title4 + '.svg')
    out_image_file5 = os.path.join(image_out_path, title5 + '.svg')
    out_files = [
        out_image_file, out_image_file2, out_image_file3, out_image_file4,
        out_image_file5
    ]
    figs = [fig, fig2, fig3, fig4, fig5]

    for i in [0, 2]:
        figs[i].savefig(out_files[i])
        # figs[i].show()

    return figs

def cf_plotter90V2(x_data, data_array, marks, legends, x_range, y_range, y_label,
                 image_out_path, figname):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 19,
        'figure.figsize': (10, 10),
        'lines.linewidth': 0.5,
        'lines.markersize': 8,
        # 'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.1,
        'figure.subplot.right': 0.85,
        'figure.subplot.top': 0.85,
        'figure.subplot.bottom': 0.15,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    markers = ['s', 's', 's']
    colorlgd = ['tab:blue', 'tab:orange']
    x_array = np.array(x_data)
    cf_array = np.array(data_array)
    cf_legends = np.array(legends)
    markr = marks[0]
    markc = marks[1]

    no_r = len(markr)
    no_c = len(markc)
    no_legend = len(cf_legends)
    no_x = len(x_array)
    fig, ax = plt.subplots(no_r, no_c)
    fig2, ax2 = plt.subplots(no_r, no_c)
    fig3, ax3 = plt.subplots(no_r, no_c)
    fig4, ax4 = plt.subplots(no_r, no_c)
    fig5, ax5 = plt.subplots(no_r, no_c)
    for r in range(no_r):
        for c in range(no_c):
            if no_r > 1:
                axi = ax[r][c]
                ax2i = ax2[r][c]
                ax3i = ax3[r][c]
                ax4i = ax4[r][c]
                ax5i = ax5[r][c]
            else:
                axi = ax[c]
                ax2i = ax2[c]
                ax3i = ax3[c]
                ax4i = ax4[c]
                ax5i = ax5[c]
            axs = [axi, ax2i, ax3i, ax4i, ax5i]
            for lgd in range(no_legend):
                data_no = no_x * (no_c * no_legend * r + no_legend * c + lgd)
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
                    datatoplot = [mcl, ratiol, mcd, ratiod, ratio_ld]

                for i in range(len(datatoplot)):
                    axs[i].plot(x_array,
                                datatoplot[i],
                                label=cf_legends[lgd],
                                color=colorlgd[lgd],
                                marker=markers[r],
                                linestyle='-.')

                    axs[i].set_xticks(np.arange(2, 12, step=2))
                    if x_range != 'all':
                        axs[i].set_xlim(x_range)
                    if y_range != 'all':
                        axs[i].set_ylim(y_range[i])

                    if r == 0:
                        if c == 1:
                            axs[i].legend(loc='upper center',
                            bbox_to_anchor=(0.5, 1.4),
                            ncol=3,
                            fontsize='small',
                            frameon=False)

                    if lgd == 0:
                        marky_loc = axs[i].get_ylim()[1] + 0.1 * (
                            axs[i].get_ylim()[1] - axs[i].get_ylim()[0])
                        markx_loc = axs[i].get_xlim()[1] + 0.25 * (
                            axs[i].get_xlim()[1] - axs[i].get_xlim()[0])
                        markymid_loc = axs[i].get_ylim()[0] + 0.5 * (
                            axs[i].get_ylim()[1] - axs[i].get_ylim()[0])

                        if r == 0:
                            axs[i].annotate(s=markc[c],
                                            xy=(6, marky_loc),
                                            ha='center',
                                            va='center',
                                            annotation_clip=False)
                        if c == 0:
                            axs[i].set_ylabel(y_label[i])
                        else:
                            axs[i].set_yticklabels([])
                        if c == no_c - 1:
                            axs[i].annotate(s=markr[r],
                            xy=(markx_loc, markymid_loc),
                            ha='center',
                            va='center',
                            annotation_clip=False)
                        if r == no_r - 1:
                            axs[i].set_xlabel(r'$s/c$')
                        else:
                            axs[i].set_xticklabels([])

                        axs[i].axhline(y=0,
                                       color='k',
                                       linestyle='-.',
                                       linewidth=0.5)

    title = 'mean lift coefficients plot_V2 ' + figname
    title2 = 'mean lift ratio plot ' + figname
    title3 = 'mean drag coefficients plot_V2 ' + figname
    title4 = 'mean drag ratio plot ' + figname
    title5 = 'mean wake ld ratio plot ' + figname
    out_image_file = os.path.join(image_out_path, title + '.svg')
    out_image_file2 = os.path.join(image_out_path, title2 + '.svg')
    out_image_file3 = os.path.join(image_out_path, title3 + '.svg')
    out_image_file4 = os.path.join(image_out_path, title4 + '.svg')
    out_image_file5 = os.path.join(image_out_path, title5 + '.svg')
    out_files = [
        out_image_file, out_image_file2, out_image_file3, out_image_file4,
        out_image_file5
    ]
    figs = [fig, fig2, fig3, fig4, fig5]

    for i in [0, 2]:
        figs[i].savefig(out_files[i])
        # figs[i].show()

    return figs
