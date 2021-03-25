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


def ktrace_plotter(ax_toplot, illustration_t, marker, iarray,
                   figure_parameters):
    """
    function to plot cfd force coefficients results
    """
    LEsize = 20
    wing_thick = 1.5

    vgap = 0.15

    data_array = np.array(iarray)
    stroke = figure_parameters[0]
    yt = 1.0 + vgap
    y_marker = 1.05 + vgap

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
    patch = patches.PathPatch(wingpath, edgecolor='k', linewidth=wing_thick)

    ax_toplot.add_patch(patch)
    LE = np.array(LE)
    ax_toplot.scatter(LE[:, 0], LE[:, 1], s=LEsize, c='k')

    #-----marks annotation--
    t_mid = illustration_t[int(len(illustration_t) / 2) - 1]
    tt_mid = t_spl(t_mid)

    ax_toplot.annotate(s=marker,
                       xy=(tt_mid, y_marker),
                       ha='center',
                       va='center',
                       annotation_clip=False)

    #----time annotation---
    ax_toplot.annotate(s=r'$\^t = 0$',
                       xy=(hinge[0], 1.7 * vgap),
                       xytext=(hinge[0], yt),
                       arrowprops=dict(arrowstyle='->', facecolor='k', lw=0.5),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    ax_toplot.annotate(s=r'$\^t\/=\/1$',
                       xy=(hinge[-1], 1.7 * vgap),
                       xytext=(hinge[-1], yt),
                       arrowprops=dict(arrowstyle='->', facecolor='k', lw=0.5),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    #--------------------------------------
    ax_toplot.set_xlim([-0.6, stroke + 0.6])
    ax_toplot.set_ylim([-0.8, 1.0])
    ax_toplot.axhline(y=0, color='k', linestyle='-.', linewidth=1)
    ax_toplot.axis('off')
    ax_toplot.set_aspect('equal')

    return ax_toplot


def illustrationk_plotter(ax_toplot, iarray, figure_parameters):
    """
    function to plot cfd force coefficients results
    """
    LEsize = 100
    wing_thick = 3.5

    vgap = 0.15
    vbarh = 0.2

    stroke = figure_parameters[0]
    data_array = np.array(iarray)
    at_time = 1.25 * figure_parameters[1] / 2
    pt_time = 1.0 - figure_parameters[2]

    illustration_t_0_1 = [0.0, pt_time, 1.0]
    illustration_t_mid = [pt_time]

    t_spl = UnivariateSpline(iarray[:, 0], -iarray[:, 1], s=0)
    aoa_spl = UnivariateSpline(iarray[:, 0], -iarray[:, 2], s=0)
    at = t_spl(at_time)
    pt = stroke - t_spl(pt_time)

    #----arrow placement---
    yp1 = 0.25 + vgap
    yp2 = yp1 + vbarh
    yparrow = yp1 + 0.5 * vbarh

    vbara1 = [[stroke - pt, yp1], [stroke - pt, yp2]]
    vbara2 = [[stroke, yp1], [stroke, yp2]]
    parrow = [[stroke - pt, yparrow],
              [stroke - pt * 0.5, yparrow + 0.8 * vgap], [stroke, yparrow]]
    #----------------------------
    vbart1 = [[0.0, -0.5 * vbarh], [0.0, 0.5 * vbarh]]
    vbart2 = [[stroke, -0.5 * vbarh], [stroke, 0.5 * vbarh]]
    sarrow1 = [[0.0, 0.0], [0.5 * stroke, 0.8 * vgap],
               [0.5 * stroke * 1.03, 0.0]]
    sarrow2 = [[0.5 * stroke, 0.0], [stroke, 0.8 * vgap], [stroke, 0.0]]

    #--side bars and arrows for motion illustration--
    vbars = vbara1 + vbara2 + vbart1 + vbart2
    arrows = [parrow, sarrow1, sarrow2]
    arrowlegend = ['plate pitch', 'plate translation', '']
    #----------------------------------

    nverts = len(vbars)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0::2] = path.Path.MOVETO
    vbarath = path.Path(vbars, codes)
    barpatch = patches.PathPatch(vbarath, edgecolor='k', linewidth=1)

    #----------------wing patch-------------
    hinge = []
    LE = []
    wing = []
    for t in illustration_t_0_1:
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
    patch = patches.PathPatch(wingpath, edgecolor='k', linewidth=wing_thick)

    #--aoa and pitch angle annotation--
    line_length = 1.1
    annotation_plength = 0.89 * line_length
    aoa_line = []
    t_aoa = illustration_t_0_1[0]
    aoat_aoa = (aoa_spl(t_aoa) - 45) * np.pi / 180
    tt_aoa = t_spl(t_aoa)

    #---points for line drawing and annotation---
    aoaline_origin = np.array([tt_aoa, 0])
    aoalineLE = np.array([
        tt_aoa - line_length * np.sin(aoat_aoa), line_length * np.cos(aoat_aoa)
    ])
    annotation_p1 = np.array([tt_aoa + annotation_plength, 0])
    annotation_p2 = np.array([
        tt_aoa - annotation_plength * np.sin(aoat_aoa),
        annotation_plength * np.cos(aoat_aoa)
    ])
    annotation_pmid = np.array([
        tt_aoa - 1.16 * annotation_plength * np.sin(aoat_aoa * 1.5),
        1.16 * annotation_plength * np.cos(aoat_aoa * 1.5)
    ])
    #----annotation points for pitch angle---
    annotation_palength = 0.2 * annotation_plength

    t_pa = [pt_time, 1.0]
    aoat_pa = np.array([
        aoa_spl(t_pa[0]),
        aoa_spl(t_pa[1]), 0.5 * (aoa_spl(t_pa[0]) + aoa_spl(t_pa[1]))
    ]) - 45 + 180
    aoat_pa = aoat_pa * np.pi / 180
    tt_pa = [
        t_spl(t_pa[0]),
        t_spl(t_pa[1]), 0.5 * (t_spl(t_pa[0]) + t_spl(t_pa[1]))
    ]

    annotation_pa1 = np.array([
        tt_pa[0] - annotation_palength * np.sin(aoat_pa[0]),
        annotation_palength * np.cos(aoat_pa[0])
    ])
    annotation_pa2 = np.array([
        tt_pa[1] - annotation_palength * np.sin(aoat_pa[1]),
        annotation_palength * np.cos(aoat_pa[1])
    ])
    annotation_pamid = np.array([
        tt_pa[2] - 2.6 * annotation_palength * np.sin(aoat_pa[2]),
        2.6 * annotation_palength * np.cos(aoat_pa[2])
    ])
    #----------------------------

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

    ax_toplot.annotate(s='',
                       xy=(annotation_p1[0], annotation_p1[1]),
                       xytext=(annotation_p2[0], annotation_p2[1]),
                       arrowprops=dict(arrowstyle='<-',
                                       facecolor='k',
                                       lw=1,
                                       connectionstyle='arc3,rad=-0.4'),
                       annotation_clip=False)
    ax_toplot.annotate(s=r'$\alpha$',
                       xy=(annotation_pmid[0], annotation_pmid[1]),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    ax_toplot.annotate(s='',
                       xy=(annotation_pa1[0], annotation_pa1[1]),
                       xytext=(annotation_pa2[0], annotation_pa2[1]),
                       arrowprops=dict(arrowstyle='<-',
                                       facecolor='k',
                                       lw=1,
                                       connectionstyle='arc3,rad=-0.4'),
                       annotation_clip=False)
    ax_toplot.annotate(s=r'$\theta$',
                       xy=(annotation_pamid[0], annotation_pamid[1]),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    #--------------------------

    ax_toplot.add_patch(aoapatch)
    ax_toplot.add_patch(patch)
    ax_toplot.add_patch(barpatch)
    LE = np.array(LE)
    ax_toplot.scatter(LE[:, 0], LE[:, 1], s=LEsize, c='k')

    for arrow, alegend in zip(arrows, arrowlegend):
        if arrow != arrows[-1]:
            arstyle = '<-'
        else:
            arstyle = '-'

        ax_toplot.annotate(s='',
                           xy=(arrow[0][0], arrow[0][1]),
                           xytext=(arrow[2][0], arrow[2][1]),
                           arrowprops=dict(arrowstyle=arstyle,
                                           linestyle='-.',
                                           facecolor='k',
                                           lw=1),
                           annotation_clip=False)
        ax_toplot.annotate(s=alegend,
                           xy=(arrow[1][0], arrow[1][1]),
                           ha='center',
                           va='center',
                           annotation_clip=False)
    #--------------------------------------
    ax_toplot.set_xlim([-0.6, stroke + 0.6])
    ax_toplot.set_ylim([-0.8, 1.0])
    ax_toplot.set_aspect('equal')
    ax_toplot.axis('off')

    return ax_toplot


def kf_plotter(kinematic_data_list, kinematics_t, data_array, idata_array,
               iarray, ifigure_parameters, marks, time_to_plot, show_range,
               image_out_path, time_scale):
    """
    function to plot cfd force coefficients results
    """
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 14,
        'figure.figsize': (14, 8),
        'lines.linewidth': 4,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 200,
    })
    data_array = np.array(data_array)

    gs_kw = dict(left=0.125,
                 right=0.9,
                 top=0.55,
                 bottom=0.2,
                 wspace=0.1,
                 hspace=0.0)

    fig, axs = plt.subplots(nrows=2, ncols=3, gridspec_kw=gs_kw)
    for axr1i, kinematic_listi, marker, iarrayi in zip(axs[0, :],
                                                       kinematic_data_list,
                                                       marks, idata_array):
        figure_parametersi = [
            float(kinematic_listi.split('_')[1].split('stroke')[1]),
            float(kinematic_listi.split('_')[2].split('acf')[1]),
            float(kinematic_listi.split('_')[3].split('pf')[1])
        ]
        ktrace_plotter(axr1i, kinematics_t, marker, iarrayi,
                       figure_parametersi)

    labels = [r'$u/U_T$', r'$\.\alpha/\.\alpha_M$']
    for axr2i, datai in zip(axs[1, :], data_array):
        axr2i.plot(datai[:, 0] / time_scale, datai[:, 1], label=labels[0])
        axr2i.plot(datai[:, 0] / time_scale, datai[:, 2], label=labels[1])
        axr2i.set_xlabel(r'Non-dimensional time $(\/\^t\/)$')
        axr2i.set_ylabel('Normalized velocity')
        axr2i.legend()
        axr2i.axvline(x=1.0, color='k', linestyle='-', linewidth=0.5)
        axr2i.label_outer()
        axr2i.set_aspect('auto')
        if time_to_plot != 'all':
            axr2i.set_xlim(time_to_plot)
        if show_range != 'all':
            axr2i.set_ylim(show_range)

    title = 'kinematics plot'
    out_image_file = os.path.join(image_out_path, title + '.png')

    bx_loc = axs[0, 0].get_xlim()[0] - 0.1 * (axs[0, 0].get_xlim()[1] -
                                              axs[0, 0].get_xlim()[0])
    by_loc = axs[0, 0].get_ylim()[1] + 0.2 * (axs[0, 0].get_ylim()[1] -
                                              axs[0, 0].get_ylim()[0])
    axs[0, 0].annotate(s='(b)',
                       xy=(bx_loc, by_loc),
                       ha='center',
                       va='center',
                       annotation_clip=False)

    #----add illustration figure-----
    gs_i = fig.add_gridspec(nrows=1,
                            ncols=1,
                            left=0.2,
                            right=0.825,
                            top=0.9,
                            bottom=0.6,
                            wspace=0.05,
                            hspace=0.0)

    axi = fig.add_subplot(gs_i[0, 0])
    illustrationk_plotter(axi, iarray, ifigure_parameters)

    ax_loc = axi.get_xlim()[0] - 0.0 * (axi.get_xlim()[1] - axi.get_xlim()[0])
    ay_loc = axi.get_ylim()[1] - 0.1 * (axi.get_ylim()[1] - axi.get_ylim()[0])
    axi.annotate(s='(a)',
                 xy=(ax_loc, ay_loc),
                 ha='center',
                 va='center',
                 annotation_clip=False)

    plt.savefig(out_image_file)
    # plt.show()

    return fig
