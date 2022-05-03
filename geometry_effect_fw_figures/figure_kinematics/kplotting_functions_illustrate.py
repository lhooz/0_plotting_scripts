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


def illustrationk_plotter(ax_toplot, kinematics_t, iarray, stroke):
    """
    function to plot cfd force coefficients results
    """
    LEsize = 100
    wing_thick = 3

    data_array = np.array(iarray)

    illustration_t_0_1 = kinematics_t

    t_spl = UnivariateSpline(iarray[:, 0], -iarray[:, 1], s=0)
    aoa_spl = UnivariateSpline(iarray[:, 0], -iarray[:, 2], s=0)

    #----------------wing patch-------------
    hinge = []
    LE = []
    wing = []
    for t in illustration_t_0_1:
        tt = t_spl(t)
        aoat = aoa_spl(t) * np.pi / 180
        LEt = np.array([tt - 0.25 * np.sin(aoat), 0.25 * np.cos(aoat)])
        TEt = np.array([tt + 0.75 * np.sin(aoat), -0.75 * np.cos(aoat)])
        hinge.append(tt)
        LE.append(LEt)
        wing.append(LEt)
        wing.append(TEt)

    nverts = len(wing)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0::2] = path.Path.MOVETO
    wingpath = path.Path(wing, codes)
    patch = patches.PathPatch(wingpath, edgecolor='k', linewidth=wing_thick)
    #--------------------------
    ax_toplot.add_patch(patch)
    LE = np.array(LE)
    ax_toplot.scatter(LE[:, 0], LE[:, 1], s=LEsize, c='k')
    #--------------------------------------
    ax_toplot.set_xlim([-0.6, stroke + 0.6])
    ax_toplot.set_ylim([-3.0, 3.0])
    ax_toplot.set_aspect('equal')
    ax_toplot.axis('off')

    return ax_toplot


def kf_plotter(iarray, kinematics_t, image_out_path, stroke):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 19,
        'figure.figsize': (6 * 2, 2.75 * 3),
        'lines.linewidth': 6,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
    })
    gs_i = dict(left=0.2,
                right=0.95,
                top=0.95,
                bottom=0.1,
                wspace=0.1,
                hspace=0.1)

    fig2, axi = plt.subplots(nrows=1, ncols=1, gridspec_kw=gs_i)
    #------at pt annotations--------------
    title = 'kinematics plot'
    out_image_file2 = os.path.join(image_out_path, title + '_illustration.svg')

    #----add illustration figure-----
    illustrationk_plotter(axi, kinematics_t, iarray, stroke)

    fig2.savefig(out_image_file2)
    # plt.show()

    return fig2
