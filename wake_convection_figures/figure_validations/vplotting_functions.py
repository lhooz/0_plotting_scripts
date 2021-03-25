"""fuctions for plotting cfd run results against ref data"""

import csv
import os
import matplotlib.patches as patches
import matplotlib.path as path
import numpy as np
import matplotlib.pyplot as plt
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


def read_ref_data(ref_data_file):
    """read wing geometry data"""
    ref_array = []
    with open(ref_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            ref_array.append([float(row[0]), float(row[1])])
            line_count += 1

        print(f'Processed {line_count} lines in {ref_data_file}')

    ref_array = np.array(ref_array)

    return ref_array


def cf_plotter(data_array, time_scale, legends, time_to_plot, show_range,
               image_out_path, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    datax_shift = 0.02 * time_scale - 0.05
    ref_shit_constant = 0.025
    ini_t = 0.02
    acc_t = 0.16 / time_scale
    y1 = 1.05
    y2 = 1.15
    ymid = 0.5 * (y1 + y2)

    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 15,
        'figure.figsize': (10, 6),
        'lines.linewidth': 2.0,
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
    kine_array = np.array(data_array[0])
    cf_array = np.array(data_array[1])
    ref_array = np.array(data_array[2])
    cf_legends = np.array(legends[0])
    ref_legends = np.array(legends[1])

    ref_array_shifted = []
    for ref_arrayi in ref_array:
        if time_to_plot == 'all':
            ref_arrayi = np.array([
                ref_arrayi[:, 0] - np.rint(ref_arrayi[0, 0]) +
                ref_shit_constant, ref_arrayi[:, 1]
            ])
        else:
            ref_arrayi = np.array([
                ref_arrayi[:, 0] - np.rint(ref_arrayi[0, 0]) +
                ref_shit_constant + time_to_plot[0], ref_arrayi[:, 1]
            ])
        ref_arrayi = np.transpose(ref_arrayi)

        ref_array_shifted.append(ref_arrayi)
        # print(ref_arrayi)

    fig1, ax1 = plt.subplots(1, 1)
    if time_to_plot != 'all':
        ax1.set_xlim(time_to_plot)
    if show_range != 'all':
        ax1.set_ylim(show_range)

    for datai in kine_array:
        ax1.plot((datai[:, 0] + datax_shift) / time_scale,
                 datai[:, 1],
                 linewidth=4)
        vbar1 = [[ini_t, y1], [ini_t, y2]]
        vbar2 = [[ini_t + acc_t, y1], [ini_t + acc_t, y2]]
        vbars = vbar1 + vbar2
        arrow1 = [[ini_t, ymid], [ini_t + acc_t, ymid]]
        text_loc = [[ini_t + acc_t, ymid], [0.5, ymid]]
        arrows = [arrow1]
        annotate_text = r'$\hat{a}_t = 0.16$'

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
                         arrowprops=dict(arrowstyle='<->',
                                         facecolor='k',
                                         lw=0.5),
                         annotation_clip=False)
        ax1.annotate(s=annotate_text,
                     xy=(text_loc[0][0], text_loc[0][1]),
                     xytext=(text_loc[1][0], text_loc[1][1]),
                     arrowprops=dict(arrowstyle='->', facecolor='k', lw=0.5),
                     ha='center',
                     va='center',
                     annotation_clip=False)

    ax1.set_xlabel(r'Non-dimensional time $(\/\^t\/)$')

    fig2, ax2 = plt.subplots(1, 2)

    if time_to_plot != 'all':
        ax2[0].set_xlim(time_to_plot)
        ax2[1].set_xlim(time_to_plot)

    mcf_arr = []
    mref_arr = []
    for i in range(len(cf_legends)):
        cf_t = (cf_array[i][:, 0] + datax_shift) / time_scale
        ax2[i].plot(cf_t, cf_array[i][:, 3], label=cf_legends[i])

        cf_spl = UnivariateSpline(cf_t, cf_array[i][:, 3], s=0)
        mcf_s = cf_spl.integral(0.0, 1.0)
        mcf_w = cf_spl.integral(1.0, 2.0)
        mcf = cf_spl.integral(0.0, 2.0)

        mcf_arr.append([mcf_s, mcf_w, mcf])

    for i in range(len(ref_legends)):
        ref_t = (ref_array_shifted[i][:, 0] + datax_shift) / time_scale
        ax2[i].plot(
            ref_t,
            ref_array_shifted[i][:, 1],
            label=ref_legends[i],
            # marker='s',
            linestyle='-.')

        ref_spl = UnivariateSpline(ref_t, ref_array_shifted[i][:, 1], s=0)
        mref_s = ref_spl.integral(0.0, 1.0)
        mref_w = ref_spl.integral(1.0, 2.0)
        mref = ref_spl.integral(0.0, 2.0)

        mref_arr.append([mref_s, mref_w, mref])

        ax2[i].set_xlabel(r'Non-dimensional time $(\/\^t\/)$')
        ax2[i].set_ylabel(r'$C_L$')
        ax2[i].axvline(x=1, color='k', linestyle='-', linewidth=0.5)
        ax2[i].axhline(y=0, color='k', linestyle='-.', linewidth=0.5)
        ax2[i].legend(frameon=False)

    ax1.set_ylabel(r'$u/U_T$')
    ax1.axvline(x=1, color='k', linestyle='-', linewidth=0.5)

    out_image_file1 = os.path.join(image_out_path, 'kinematics.png')
    out_image_file2 = os.path.join(image_out_path, 'coefficients.png')
    # ax.set_title(title)

    with open('meancf.dat', 'w') as f:
        for item, cf_lgd in zip(mcf_arr, cf_legends):
            f.write("%s:\n" % cf_lgd)
            f.write("mcf_s = %s, mcf_w = %s, mcf = %s\n" %
                    ('{0:.8g}'.format(item[0]), '{0:.8g}'.format(
                        item[1]), '{0:.8g}'.format(item[2])))
        for item, ref_lgd in zip(mref_arr, ref_legends):
            f.write("%s:\n" % ref_lgd)
            f.write("mcf_s = %s, mcf_w = %s, mref = %s\n" %
                    ('{0:.8g}'.format(item[0]), '{0:.8g}'.format(
                        item[1]), '{0:.8g}'.format(item[2])))

    fig1.set_size_inches(10, 4)
    fig2.set_size_inches(15, 4)
    fig1.savefig(out_image_file1)
    fig2.savefig(out_image_file2)
    plt.show()

    return fig1, fig2


def append_kinematics_array(cfd_arr, kinematics_data_file):
    """read stroke angle and append to cfd array"""
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
                phi_datai = row[-1].split()[0]
                kinematics_arr.append(
                    [float(t_datai), np.abs(float(phi_datai))])
                line_count += 1

        print(f'Processed {line_count} lines in {kinematics_data_file}')

    kinematics_arr = np.array(kinematics_arr)
    phi_spl = UnivariateSpline(kinematics_arr[:, 0], kinematics_arr[:, 1])

    phi = []
    for ti in cfd_arr[:, 0]:
        phii = phi_spl(ti)
        phi.append([phii])
    phi = np.array(phi)

    cfd_arr = np.append(cfd_arr, phi, axis=1)

    return cfd_arr
