"""fuctions for plotting cfd run results against ref data"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline


def read_kinematics_data(kinematics_data_file):
    """read kinematics from file"""
    u2data_file = kinematics_data_file + '.cf'
    kdata_file = kinematics_data_file + '.dat'
    with open(u2data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0

        for row in csv_reader:
            if line_count != 22:
                line_count += 1
            else:
                u2 = row[-1].split(';')[0]
                u2 = float(u2)
                line_count += 1

    kinematics_arr = []
    with open(kdata_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='(')
        line_count = 0

        for row in csv_reader:
            if line_count <= 1:
                line_count += 1
            elif row[0] == ')':
                line_count += 1
            else:
                t_datai = row[1]
                # print(row)
                trans_datai0 = row[3].split()[0]
                trans_datai1 = row[3].split()[1]
                trans_datai2 = row[3].split()[2].split(')')[0]

                rot_datai0 = row[4].split()[0]
                rot_datai1 = row[4].split()[1]
                rot_datai2 = row[4].split()[2].split(')')[0]
                kinematics_arr.append([
                    float(t_datai),
                    float(trans_datai0),
                    float(trans_datai1),
                    float(trans_datai2),
                    float(rot_datai0),
                    float(rot_datai1),
                    float(rot_datai2)
                ])
                line_count += 1

        print(f'Processed {line_count} lines in {kinematics_data_file}')

    kinematics_arr = np.array(kinematics_arr)

    dkinematics_arr = np.array([kinematics_arr[:, 0]])
    ddkinematics_arr = np.array([kinematics_arr[:, 0]])
    for i in range(6):
        spl = UnivariateSpline(kinematics_arr[:, 0],
                               kinematics_arr[:, i + 1],
                               s=0)
        dk = []
        ddk = []
        for t in kinematics_arr[:, 0]:
            dki = spl.derivatives(t)[1]
            ddki = spl.derivatives(t)[2]
            dk.append(dki)
            ddk.append(ddki)
        dk = np.array([dk])
        ddk = np.array([ddk])

        dkinematics_arr = np.append(dkinematics_arr, dk, axis=0)
        ddkinematics_arr = np.append(ddkinematics_arr, ddk, axis=0)
    dkinematics_arr = np.transpose(dkinematics_arr)
    ddkinematics_arr = np.transpose(ddkinematics_arr)

    return u2, kinematics_arr, dkinematics_arr, ddkinematics_arr


def read_cfd_data(cfd_data_file, u2, karr, dkarr):
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
    cmx_spl = UnivariateSpline(cf_array[:, 0], cf_array[:, 6], s=0)
    cmy_spl = UnivariateSpline(cf_array[:, 0], -cf_array[:, 5], s=0)
    cmz_spl = UnivariateSpline(cf_array[:, 0], cf_array[:, 4], s=0)

    omegax_spl = UnivariateSpline(dkarr[:, 0], dkarr[:, 4] * np.pi / 180, s=0)
    omegay_spl = UnivariateSpline(
        dkarr[:, 0],
        np.multiply(dkarr[:, 5], np.cos(karr[:, 4] * np.pi / 180)) * np.pi /
        180,
        s=0)
    omegaz_spl = UnivariateSpline(
        dkarr[:, 0],
        np.multiply(dkarr[:, 5], np.sin(karr[:, 4] * np.pi / 180)) * np.pi /
        180,
        s=0)

    t_arr = cf_array[:, 0]
    cpx_arr = [-1 * cmx_spl(t) * omegax_spl(t) / u2 for t in t_arr]
    cpy_arr = [-1 * cmy_spl(t) * omegay_spl(t) / u2 for t in t_arr]
    cpz_arr = [-1 * cmz_spl(t) * omegaz_spl(t) / u2 for t in t_arr]

    cpx_spl = UnivariateSpline(t_arr, cpx_arr, s=0)
    cpy_spl = UnivariateSpline(t_arr, cpy_arr, s=0)
    cpz_spl = UnivariateSpline(t_arr, cpz_arr, s=0)

    #--------plot axes-------
    # plt.plot(t_arr, cpx_spl(t_arr), t_arr, cpy_spl(t_arr) + cpz_spl(t_arr))
    # plt.plot(t_arr, cpx_spl(t_arr), t_arr, cpy_spl(t_arr), t_arr,
    # cpz_spl(t_arr))
    # plt.legend(['x', 'y', 'z'], loc='best')
    # plt.xlim([4.0, 5.0])
    # plt.show()
    #------------------------
    mcl = cl_spl.integral(4.0, 5.0)
    mcp = cpx_spl.integral(4.0, 5.0) + cpy_spl.integral(
        4.0, 5.0) + cpz_spl.integral(4.0, 5.0)

    mcf_array = [mcl, mcp]

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
        'font.size': 14,
        'figure.figsize': (14, 8),
        'lines.linewidth': 0.5,
        'lines.markersize': 12,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 200,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.85,
        'figure.subplot.top': 0.8,
        'figure.subplot.bottom': 0.2,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    markers = ['o', 'v', 's']
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
    for r in range(no_r):
        for c in range(no_c):
            axs = [ax[r][c], ax2[r][c]]
            for lgd in range(no_legend):
                data_no = no_x * (no_legend * no_c * r + no_legend * c + lgd)
                # print(data_no)
                mcl = []
                mcp = []
                for xi in range(no_x):
                    mcl.append(cf_array[data_no + xi][0])
                    mcp.append(cf_array[data_no + xi][0] /
                               (cf_array[data_no + xi][1]**(2 / 3)))

                datatoplot = [mcl, mcp]
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

                    if c == 1 and r == 0:
                        axs[i].legend(loc='upper center',
                                      bbox_to_anchor=(0.5, 1.2),
                                      ncol=3,
                                      fontsize='small',
                                      frameon=False)

                    if lgd == 0:
                        markx_loc = axs[i].get_xlim()[1] + 0.3 * (
                            axs[i].get_xlim()[1] - axs[i].get_xlim()[0])
                        markxmid_loc = axs[i].get_xlim()[0] + 0.5 * (
                            axs[i].get_xlim()[1] - axs[i].get_xlim()[0])
                        marky_loc = axs[i].get_ylim()[1] + 0.2 * (
                            axs[i].get_ylim()[1] - axs[i].get_ylim()[0])
                        markymid_loc = axs[i].get_ylim()[0] + 0.5 * (
                            axs[i].get_ylim()[1] - axs[i].get_ylim()[0])

                        if r == 0:
                            axs[i].annotate(s=markc[c],
                                            xy=(markxmid_loc, marky_loc),
                                            ha='center',
                                            va='center',
                                            annotation_clip=False)
                        if c == no_c - 1:
                            axs[i].annotate(s=markr[r],
                                            xy=(markx_loc, markymid_loc),
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

    title = 'mean lift coefficients plot'
    title2 = 'efficiency plot'
    out_image_file = os.path.join(image_out_path, title + '.png')
    out_image_file2 = os.path.join(image_out_path, title2 + '.png')
    out_files = [out_image_file, out_image_file2]
    figs = [fig, fig2]

    for i in range(len(figs)):
        figs[i].savefig(out_files[i])
        figs[i].show()

    return figs
