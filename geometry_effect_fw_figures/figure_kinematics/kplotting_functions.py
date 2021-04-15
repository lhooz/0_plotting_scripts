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


def k_plotter(kine_array, time_scale, time_to_plot, show_range, image_out_path,
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
        'font.size': 15,
        'figure.figsize': (10, 5),
        'lines.linewidth': 4.0,
        'lines.markersize': 4,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
        'figure.subplot.left': 0.125,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.9,
        'figure.subplot.bottom': 0.2,
        'figure.subplot.wspace': 0.1,
        'figure.subplot.hspace': 0.1,
    })

    fig1, ax1 = plt.subplots(1, 1)
    if time_to_plot != 'all':
        ax1.set_xlim(time_to_plot)
    if show_range != 'all':
        ax1.set_ylim(show_range)
    ax1.axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
    # ax1.axvline(x=0.125, color='k', linestyle='-.', linewidth=0.5)
    # ax1.axvline(x=0.375, color='k', linestyle='-.', linewidth=0.5)
    ax1.axvline(x=0.5, color='k', linestyle='-.', linewidth=0.5)
    # ax1.axvline(x=0.625, color='k', linestyle='-.', linewidth=0.5)
    # ax1.axvline(x=0.875, color='k', linestyle='-.', linewidth=0.5)

    labels = [r'$\phi/A$', r'$\alpha/\Theta$']
    ax1.plot(kine_array[:, 0] / time_scale, kine_array[:, 1], label=labels[0])
    ax1.plot(kine_array[:, 0] / time_scale, kine_array[:, 2], label=labels[1])
    ax1.legend(frameon=False)

    vbar1 = [[ini_t, y1], [ini_t, y2]]
    vbar2 = [[ini_t + acc_t, y1], [ini_t + acc_t, y2]]
    vbars = vbar1 + vbar2
    arrow1 = [[ini_t, ymid], [ini_t + acc_t, ymid]]
    text_loc = [[ini_t + acc_t, ymid], [0.25, ymid]]
    arrows = [arrow1]
    annotate_text = r'$\frac{1}{2}\hat{t}_a = \frac{1}{2}\hat{t}_p = 0.125$'

    nverts = len(vbars)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0::2] = path.Path.MOVETO
    vbarath = path.Path(vbars, codes)
    barpatch = patches.PathPatch(vbarath, edgecolor='k', linewidth=0.5)

    ax1.add_patch(barpatch)

    for arrow in arrows:
        ax1.annotate(s='',
                     xy=(arrow[0][0], arrow[0][1]),
                     xytext=(arrow[1][0], arrow[1][1]),
                     arrowprops=dict(arrowstyle='<->', facecolor='k', lw=0.5),
                     annotation_clip=False)
    ax1.annotate(s=annotate_text,
                 xy=(text_loc[0][0], text_loc[0][1]),
                 xytext=(text_loc[1][0], text_loc[1][1]),
                 arrowprops=dict(arrowstyle='->', facecolor='k', lw=0.5),
                 ha='center',
                 va='center',
                 annotation_clip=False)

    ax1.set_xlabel(r'Non-dimensional time $(\/\^t\/)$')

    ax1.set_ylabel(r'Normalized angle')
    ax1.axvline(x=1, color='k', linestyle='-', linewidth=0.5)

    out_image_file1 = os.path.join(image_out_path, 'kinematics.png')

    fig1.savefig(out_image_file1)
    plt.show()

    return fig1
