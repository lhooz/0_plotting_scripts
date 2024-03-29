"""plotting functions for mesh and convergence figures"""

import csv
import os

import matplotlib.image as mpimg
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline


def read_cfd_data(kinematics_file, cfd_data_file):
    """read cfd results force coefficients data"""
    kinematics_arr = []
    with open(kinematics_file) as csv_file:
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
                rot_datai0 = row[4].split()[0]
                kinematics_arr.append([
                    float(t_datai),
                    float(rot_datai0) * np.pi / 180,
                ])
                line_count += 1

        print(f'Processed {line_count} lines in {kinematics_file}')

    kinematics_arr = np.array(kinematics_arr)
    spl = UnivariateSpline(kinematics_arr[:, 0], kinematics_arr[:, 1], s=0)

    cf_array = []
    with open(cfd_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_count = 0

        for row in csv_reader:
            if line_count <= 14:
                line_count += 1
            else:
                ti = float(row[0])
                phii = spl(ti)
                dphii = -1.0 * spl.derivatives(ti)[1]
                cli = float(row[3])
                cdi = np.sign(dphii) * (np.sin(phii) * float(row[2]) +
                                        np.cos(phii) * float(row[1]))
                csi = np.cos(phii) * float(row[2]) - np.sin(phii) * float(
                    row[1])
                cf_array.append([ti, cli, cdi, csi])
                line_count += 1

        print(f'Processed {line_count} lines in {cfd_data_file}')

    cf_array = np.array(cf_array)
    return cf_array


def cf_plotter(AR, Re, data_array, legends, time_to_plot, show_range,
               image_out_path, cycle_time, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 24,
        'figure.figsize': (6, 10),
        'lines.linewidth': 4.0,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.2,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.9,
        'figure.subplot.bottom': 0.1,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })
    legendx = 0.5
    legendy = 1.15

    time_series = np.linspace(time_to_plot[0], time_to_plot[1], 1000)

    cf_array1S = np.array(data_array[0])
    cf_array5S = np.array(data_array[1])

    cl1s_spl = UnivariateSpline(cf_array1S[:, 0] / cycle_time,
                                cf_array1S[:, 1],
                                s=0)
    cd1s_spl = UnivariateSpline(cf_array1S[:, 0] / cycle_time,
                                cf_array1S[:, 2],
                                s=0)
    cl5s_spl = UnivariateSpline((cf_array5S[:, 0] - 4.0) / cycle_time,
                                cf_array5S[:, 1],
                                s=0)
    cd5s_spl = UnivariateSpline((cf_array5S[:, 0] - 4.0) / cycle_time,
                                cf_array5S[:, 2],
                                s=0)
    cf_1s = []
    cf_5s = []
    cf_wake = []
    for ti in time_series:
        cl_wake = cl5s_spl(ti) - cl1s_spl(ti)
        cd_wake = cd5s_spl(ti) - cd1s_spl(ti)

        cf_1s.append([ti, cl1s_spl(ti), cd1s_spl(ti)])
        cf_5s.append([ti, cl5s_spl(ti), cd5s_spl(ti)])
        cf_wake.append([ti, cl_wake, cd_wake])

    cf_1s = np.array(cf_1s)
    cf_5s = np.array(cf_5s)
    cf_wake = np.array(cf_wake)

    cf_array = [cf_1s, cf_5s, cf_wake]

    range_cl = show_range[0]
    range_cd = show_range[1]

    fig, axs = plt.subplots(2, 1)
    mcl_arr = []
    mcd_arr = []
    if plot_mode == 'against_t':
        for i in range(len(legends)):
            axs[0].plot(cf_array[i][:, 0], cf_array[i][:, 1], label=legends[i])
            axs[1].plot(cf_array[i][:, 0], cf_array[i][:, 2], label=legends[i])

            cl_spl = UnivariateSpline(cf_array[i][:, 0],
                                      cf_array[i][:, 1],
                                      s=0)
            mcl = cl_spl.integral(time_to_plot[0], time_to_plot[1])

            cd_spl = UnivariateSpline(cf_array[i][:, 0],
                                      cf_array[i][:, 2],
                                      s=0)
            mcd = cd_spl.integral(time_to_plot[0], time_to_plot[1])

            mcl_arr.append(mcl)
            mcd_arr.append(mcd)

        with open('meancf_AR.dat', 'w') as f:
            for item, cf_lgd in zip(mcl_arr, legends):
                f.write("%s:\n" % cf_lgd)
                f.write("mcl = %s\n" % '{0:.8g}'.format(item))
            f.write("\n")
            for item, ref_lgd in zip(mcd_arr, legends):
                f.write("%s:\n" % ref_lgd)
                f.write("mcd = %s\n" % '{0:.8g}'.format(item))
            f.write("\n")

        if time_to_plot != 'all':
            axs[0].set_xlim(time_to_plot)
            axs[1].set_xlim(time_to_plot)
        if range_cl != 'all':
            axs[0].set_ylim(range_cl)
            axs[1].set_ylim(range_cd)

        for ax in axs:
            ax.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
            ax.axvline(x=0.5, color='k', linestyle='-.', linewidth=0.5)
            ax.set_xlabel(r'$\^t$')

        axs[0].set_ylabel(r'$C_L$')
        axs[1].set_ylabel(r'$C_D$')

        axs[0].label_outer()
        axs[1].label_outer()

        axs[0].legend(
            loc='upper center',
            bbox_to_anchor=(legendx, legendy),
            # ncol=len(legends),
            ncol=2,
            fontsize='small',
            frameon=False)

    title = 'transient_force_AR' + '{0:.1f}'.format(
        AR) + '_Re' + '{0:.1f}'.format(Re)
    out_image_file = os.path.join(image_out_path, title + '.svg')
    fig.savefig(out_image_file)
    # plt.show()

    return fig
