"""fuctions for plotting cfd run results against ref data"""

import csv
import os

import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline


def read_kinematics_data(kinematics_data_file):
    """read kinematics from file"""
    kinematics_arr = []
    with open(kinematics_data_file) as csv_file:
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

    return kinematics_arr, dkinematics_arr, ddkinematics_arr


def k_plotter(karray_all, time_scale, time_to_plot, show_range, image_out_path,
              plot_mode):
    """
    function to plot cfd force coefficients results
    """
    ini_t = 0
    acc_t = 0.125
    y1 = 1.2
    y2 = 1.3
    ymid = 0.5 * (y1 + y2)

    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (10, 10),
        'lines.linewidth': 2.0,
        'lines.markersize': 4,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.9,
        'figure.subplot.bottom': 0.15,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })

    fig1, ax = plt.subplots(2, 1)

    labels = [r'$\Omega/\Omega_M$', r'AoA']
    legendax2 = [
        r'$C_\theta = 0$', r'$C_\theta = 1$', r'$C_\theta = 5$',
        r'$C_\theta = 10$', r'$C_\theta = 100$'
    ]

    for kine_array, lgd2 in zip(karray_all, legendax2):
        if lgd2 == legendax2[0]:
            ax[0].plot(kine_array[:, 0] / time_scale, kine_array[:, 1])
        ax[1].plot(kine_array[:, 0] / time_scale, kine_array[:, 2], label=lgd2)

        vbar1 = [[ini_t, y1], [ini_t, y2]]
        vbar2 = [[ini_t + acc_t, y1], [ini_t + acc_t, y2]]
        vbars = vbar1 + vbar2
        arrow1 = [[ini_t, ymid], [ini_t + acc_t, ymid]]
        text_loc = [[ini_t + acc_t, ymid], [0.25, ymid]]
        arrows = [arrow1]
        annotate_text = r'$\frac{1}{2}\hat{t}_a = \frac{1}{2}\hat{t}_p = 0.125$'

    for axi, labeli in zip(ax, labels):
        axi.set_xticks(np.arange(0, 2.0, 0.25))
        axi.set_ylabel(labeli)
        if time_to_plot != 'all':
            axi.set_xlim(time_to_plot)
        if show_range != 'all':
            axi.set_ylim(show_range)

        axi.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
        axi.axvline(x=0.5, color='k', linestyle='-.', linewidth=0.5)
        axi.set_xlabel(r'$t/T$')

        axi.label_outer()

    ax[1].set_ylim([40, 100])
    ax[1].legend(
        # loc='upper center',
        # bbox_to_anchor=(0.5, 1.35),
        ncol=5,
        fontsize='small',
        frameon=False)

    out_image_file1 = os.path.join(image_out_path, 'kinematics.svg')

    fig1.savefig(out_image_file1)
    # plt.show()

    return fig1
