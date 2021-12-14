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


def cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
               cycle_time, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 24,
        'figure.figsize': (15, 5),
        'lines.linewidth': 3.0,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.8,
        'figure.subplot.bottom': 0.2,
        'figure.subplot.wspace': 0.2,
        'figure.subplot.hspace': 0.1,
    })
    legendx = 1.1
    legendy = 1.2

    cf_array = np.array(data_array)
    range_cl = show_range[0]
    range_cd = show_range[1]

    fig, axs = plt.subplots(1, 2)
    mcl_arr = []
    mcd_arr = []
    if plot_mode == 'against_t':
        for i in range(len(legends)):
            axs[0].plot(cf_array[i][:, 0] / cycle_time,
                        cf_array[i][:, 1],
                        label=legends[i])
            axs[1].plot(cf_array[i][:, 0] / cycle_time,
                        cf_array[i][:, 2],
                        label=legends[i])

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

            with open('meancf_convergence.dat', 'w') as f:
                for item, cf_lgd in zip(mcl_arr, legends):
                    f.write("%s:\n" % cf_lgd)
                    f.write("mcl = %s\n" % '{0:.8g}'.format(item))
                for item, ref_lgd in zip(mcd_arr, legends):
                    f.write("%s:\n" % ref_lgd)
                    f.write("mcd = %s\n" % '{0:.8g}'.format(item))

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

        axs[0].set_ylabel(r'$C_L$')
        axs[1].set_ylabel(r'$C_D$')

        axs[0].legend(loc='upper center',
                      bbox_to_anchor=(legendx, legendy),
                      ncol=4,
                      fontsize='small',
                      frameon=False)

    title = 'convergence plot'
    out_image_file = os.path.join(image_out_path, title + '.svg')
    fig.savefig(out_image_file)
    # plt.show()

    return fig
