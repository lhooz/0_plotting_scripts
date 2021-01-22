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


def kf_plotter(kinematic_data_list, data_array, legends, time_to_plot,
               show_range, image_out_path, time_scale, plot_mode):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (10, 6),
        'lines.linewidth': 4,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
    })
    data_array = np.array(data_array)
    legendt = legends[0]
    legendp = legends[1]

    fig, (ax1, ax2) = plt.subplots(2, 1)
    for kinematici, datai, legendti, legendpi in zip(kinematic_data_list,
                                                     data_array, legendt,
                                                     legendp):
        if plot_mode == 'against_t':
            if 'pf0.125' in kinematici:
                ax1.plot(datai[:, 0] / time_scale, datai[:, 1], label=legendti)

            if 'acf0.125' in kinematici:
                ax2.plot(datai[:, 0] / time_scale,
                         45 + datai[:, 2],
                         label=legendpi)

    # ax1.set_xlabel(r'$\^t$')
    ax1.set_xticklabels([])
    ax2.set_xlabel(r'Non-dimensional time $(\/\^t\/)$')

    ax1.set_ylabel(r'$u/U_T$')
    ax2.set_ylabel(r'$\alpha,\/\deg$')
    ax1.legend()
    ax2.legend()
    title = 'kinematics plot'
    out_image_file = os.path.join(image_out_path, title + '.png')
    # fig.suptitle(title)

    ax1.axvline(x=1.0, color='k', linestyle='-', linewidth=0.5)
    ax2.axvline(x=1.0, color='k', linestyle='-', linewidth=0.5)
    if time_to_plot != 'all':
        ax1.set_xlim(time_to_plot)
        ax2.set_xlim(time_to_plot)
    if show_range != 'all':
        ax1.set_ylim(show_range)
        ax2.set_ylim(show_range)

    plt.savefig(out_image_file)
    plt.show()

    return fig


def illustrationk_plotter(illustration_t, iarray, figure_parameters,
                          image_out_path):
    """
    function to plot cfd force coefficients results
    """
    width = 6
    whratio = 2.8
    vgap = 0.15
    vbarh = 0.3
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 25,
        'figure.figsize': (width * whratio, width),
        'lines.linewidth': 4,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
    })
    data_array = np.array(iarray)
    stroke = figure_parameters[0]
    at = figure_parameters[1] / 2 * stroke
    pt = figure_parameters[2] * stroke
    y1 = 0.25 + vgap
    y2 = y1 + vbarh
    yarrow1 = y1 + 0.5 * vbarh
    y3 = y2 + vgap
    y4 = y3 + vbarh
    yarrow2 = y3 + 0.5 * vbarh

    vbart1 = [[0, y1], [0, y2]]
    vbart2 = [[at, y1], [at, y2]]
    arrow1 = [[0, yarrow1], [at * 0.5, yarrow1 + 0.8 * vgap], [at, yarrow1]]
    vbart3 = [[stroke - at, y1], [stroke - at, y2]]
    vbart4 = [[stroke, y1], [stroke, y2]]
    arrow2 = [[stroke - at, yarrow1],
              [stroke - at * 0.5, yarrow1 + 0.8 * vgap], [stroke, yarrow1]]

    vbara1 = [[stroke - pt, y3], [stroke - pt, y4]]
    vbara2 = [[stroke, y3], [stroke, y4]]
    arrow3 = [[stroke - pt, yarrow2],
              [stroke - pt * 0.5, yarrow2 + 0.8 * vgap], [stroke, yarrow2]]
    vbars = vbart1 + vbart2 + vbart3 + vbart4 + vbara1 + vbara2
    arrows = [arrow1, arrow2, arrow3]
    arrowlegend = [r'$\frac{1}{2}\^a_t$', r'$\frac{1}{2}\^a_t$', r'$\^p_t$']

    nverts = len(vbars)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0::2] = path.Path.MOVETO
    vbarath = path.Path(vbars, codes)
    barpatch = patches.PathPatch(vbarath, edgecolor='k', linewidth=1)

    t_spl = UnivariateSpline(iarray[:, 0], -iarray[:, 1], s=0)
    aoa_spl = UnivariateSpline(iarray[:, 0], -iarray[:, 2], s=0)
    hinge = []
    LE = []
    wing = []
    for t in illustration_t:
        tt = t_spl(t)
        aoat = (aoa_spl(t) - 45) * np.pi / 180
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
    patch = patches.PathPatch(wingpath, edgecolor='k', linewidth=4)

    #--aoa line for annotation--
    line_length = 1.1
    annotation_plength = 0.89 * line_length
    aoa_line = []
    t_mid = illustration_t[int(len(illustration_t) / 2) - 2]
    aoat_mid = (aoa_spl(t_mid) - 45) * np.pi / 180
    tt_mid = t_spl(t_mid)

    #---points for line drawing and annotation---
    aoaline_origin = np.array([tt_mid, 0])
    annotation_p1 = np.array([tt_mid + annotation_plength, 0])
    aoalineLE = np.array([
        tt_mid - line_length * np.sin(aoat_mid), line_length * np.cos(aoat_mid)
    ])
    annotation_p2 = np.array([
        tt_mid - annotation_plength * np.sin(aoat_mid),
        annotation_plength * np.cos(aoat_mid)
    ])
    annotation_pmid = np.array([
        tt_mid - line_length * np.sin(aoat_mid * 1.5),
        line_length * np.cos(aoat_mid * 1.5)
    ])
    #-----------------------------

    aoa_line.append(aoalineLE)
    aoa_line.append(aoaline_origin)

    nverts = len(aoa_line)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0::2] = path.Path.MOVETO
    aoalinepath = path.Path(aoa_line, codes)
    aoapatch = patches.PathPatch(aoalinepath,
                                 edgecolor='k',
                                 linewidth=1,
                                 linestyle='-.')

    fig, ax1 = plt.subplots(1, 1)
    ax1.annotate(s='',
                 xy=(annotation_p1[0], annotation_p1[1]),
                 xytext=(annotation_p2[0], annotation_p2[1]),
                 arrowprops=dict(arrowstyle='<->',
                                 facecolor='k',
                                 lw=1,
                                 connectionstyle='arc3,rad=-0.4'),
                 annotation_clip=False)
    ax1.annotate(s=r'$\alpha$',
                 xy=(annotation_pmid[0], annotation_pmid[1]),
                 ha='center',
                 va='center',
                 annotation_clip=False)
    #--------------------------

    ax1.add_patch(patch)
    ax1.add_patch(aoapatch)
    ax1.add_patch(barpatch)
    LE = np.array(LE)
    ax1.scatter(LE[:, 0], LE[:, 1], s=250, c='k')

    for arrow, alegend in zip(arrows, arrowlegend):
        ax1.annotate(s='',
                     xy=(arrow[0][0], arrow[0][1]),
                     xytext=(arrow[2][0], arrow[2][1]),
                     arrowprops=dict(arrowstyle='<->', facecolor='k', lw=1),
                     annotation_clip=False)
        ax1.annotate(s=alegend,
                     xy=(arrow[1][0], arrow[1][1]),
                     ha='center',
                     va='center',
                     annotation_clip=False)

    #----time annotation---
    ax1.annotate(s=r'$\^t = 0$',
                 xy=(hinge[0], y2 + 0.5 * vgap),
                 xytext=(hinge[0], y4),
                 arrowprops=dict(arrowstyle='->', facecolor='k', lw=1),
                 ha='center',
                 va='center',
                 annotation_clip=False)
    ax1.annotate(s=r'$\^t\/=\/1$',
                 xy=(hinge[-1], 0.5 * vgap),
                 xytext=(hinge[-1] + 0.25, 0.75 * y1),
                 arrowprops=dict(
                     arrowstyle='->',
                     facecolor='k',
                     lw=1,
                     connectionstyle='angle,angleA=0,angleB=90,rad=0'),
                 ha='center',
                 va='center',
                 annotation_clip=False)
    #--------------------------------------
    ax1.set_xlim([-0.6, stroke + 0.6])
    ax1.set_ylim([-1, (stroke + 2) / whratio - 1])
    ax1.axhline(y=0, color='k', linestyle='-.', linewidth=1)

    out_image_file = os.path.join(image_out_path, 'wing_path.png')

    plt.axis('off')
    plt.savefig(out_image_file)
    plt.show()

    return fig
