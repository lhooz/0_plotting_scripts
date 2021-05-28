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
    LEsize = 30
    wing_thick = 2.5

    vgap = 0.15

    data_array = np.array(iarray)
    stroke = figure_parameters[0]
    yt = 0.7 + vgap
    y_marker = 1.18 + vgap

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
    # t_mid = illustration_t[int(len(illustration_t) / 2) - 1]
    # tt_mid = t_spl(t_mid)

    ax_toplot.annotate(s=marker,
                       xy=(hinge[-1], y_marker),
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
    ax_toplot.set_ylim([-0.8, 1.5])
    ax_toplot.axhline(y=0, color='k', linestyle='-.', linewidth=1)
    ax_toplot.axis('off')
    ax_toplot.set_aspect('equal')

    return ax_toplot


def illustrationk_plotter(ax_toplot, iarray, figure_parameters):
    """
    function to plot cfd force coefficients results
    """
    LEsize = 250
    wing_thick = 6

    vgap = 0.65
    vbarh = 0.3

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
              [stroke - pt * 0.5, yparrow + 0.2 * vgap], [stroke, yparrow]]
    #----------------------------
    vbart1 = [[0.0, -0.5 * vbarh], [0.0, 0.5 * vbarh]]
    vbart2 = [[stroke, -0.5 * vbarh], [stroke, 0.5 * vbarh]]
    sarrow1 = [[0.0, 0.0], [0.5 * stroke, 0.2 * vgap],
               [0.5 * stroke * 1.03, 0.0]]
    sarrow2 = [[0.5 * stroke, 0.0], [stroke, 0.2 * vgap], [stroke, 0.0]]

    #--side bars and arrows for motion illustration--
    vbars = vbara1 + vbara2 + vbart1 + vbart2
    arrows = [parrow, sarrow1, sarrow2]
    arrowlegend = ['Plate rotation', 'Plate translation', '']
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

    #==========AOAs annotation===================
    #--aoa annotation--
    line_length = 1.22
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
        tt_aoa - 1.2 * annotation_plength * np.sin(aoat_aoa * 1.5),
        1.2 * annotation_plength * np.cos(aoat_aoa * 1.5)
    ])

    aoa_line.append(aoalineLE)
    aoa_line.append(aoaline_origin)

    #--aoaE annotation--
    aoaE_line = []
    t_aoa = illustration_t_0_1[-1]
    aoat_aoa = (aoa_spl(t_aoa) - 45) * np.pi / 180
    tt_aoa = t_spl(t_aoa)

    #---points for line drawing and annotation---
    aoaline_origin = np.array([tt_aoa, 0])
    aoalineLE = np.array([
        tt_aoa - line_length * np.sin(aoat_aoa), line_length * np.cos(aoat_aoa)
    ])
    annotation_p1E = np.array([tt_aoa - annotation_plength, 0])
    annotation_p2E = np.array([
        tt_aoa - annotation_plength * np.sin(aoat_aoa),
        annotation_plength * np.cos(aoat_aoa)
    ])
    annotation_pmidE = np.array([
        tt_aoa - 1.42 * annotation_plength * np.sin(aoat_aoa * 1.22),
        1.2 * annotation_plength * np.cos(aoat_aoa * 1.5)
    ])

    aoaE_line.append(aoalineLE)
    aoaE_line.append(aoaline_origin)

    #========================================================

    #----annotation points for pitch angle---
    annotation_palength = 0.25 * annotation_plength

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
        tt_pa[2] - 2.4 * annotation_palength * np.sin(aoat_pa[2]),
        2.4 * annotation_palength * np.cos(aoat_pa[2])
    ])
    #----------------------------

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
                                       connectionstyle='arc3,rad=-0.5'),
                       annotation_clip=False)
    ax_toplot.annotate(s=r'$\alpha$',
                       xy=(annotation_pmid[0], annotation_pmid[1]),
                       ha='center',
                       va='center',
                       annotation_clip=False)

    nvertsE = len(aoaE_line)
    codesE = np.ones(nverts, int) * path.Path.LINETO
    codesE[0::2] = path.Path.MOVETO
    aoalinepathE = path.Path(aoaE_line, codesE)
    aoapatchE = patches.PathPatch(aoalinepathE,
                                  edgecolor='k',
                                  linewidth=1,
                                  linestyle='-.')

    ax_toplot.annotate(s='',
                       xy=(annotation_p1E[0], annotation_p1E[1]),
                       xytext=(annotation_p2E[0], annotation_p2E[1]),
                       arrowprops=dict(arrowstyle='<-',
                                       facecolor='k',
                                       lw=1,
                                       connectionstyle='arc3,rad=0.5'),
                       annotation_clip=False)
    ax_toplot.annotate(s=r'$\alpha_{int}$',
                       xy=(annotation_pmidE[0], annotation_pmidE[1]),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    ax_toplot.annotate(s='',
                       xy=(annotation_pa1[0], annotation_pa1[1]),
                       xytext=(annotation_pa2[0], annotation_pa2[1]),
                       arrowprops=dict(arrowstyle='<-',
                                       linestyle='-',
                                       facecolor='k',
                                       lw=1,
                                       connectionstyle='arc3,rad=-0.5'),
                       annotation_clip=False)
    ax_toplot.annotate(s=r'$\theta$',
                       xy=(annotation_pamid[0], annotation_pamid[1]),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    #--------------------------
    #---------lift and drag frame annotation--------
    ax_toplot.annotate(s='Lift',
                       xy=(-1.0 + 0.5, -0.022),
                       xytext=(-1.0 + 0.5, 1.0),
                       arrowprops=dict(
                           arrowstyle='<-',
                           facecolor='k',
                           lw=1,
                       ),
                       ha='center',
                       va='center',
                       annotation_clip=False)

    ax_toplot.annotate(s='Drag',
                       xy=(-0.978 + 0.5, 0.0),
                       xytext=(-2.0 + 0.5, 0.0),
                       arrowprops=dict(
                           arrowstyle='<-',
                           facecolor='k',
                           lw=1,
                       ),
                       ha='center',
                       va='center',
                       annotation_clip=False)
    #------------------------------------------------

    ax_toplot.add_patch(aoapatch)
    ax_toplot.add_patch(aoapatchE)
    ax_toplot.add_patch(patch)
    ax_toplot.add_patch(barpatch)
    LE = np.array(LE)
    ax_toplot.scatter(LE[:, 0], LE[:, 1], s=LEsize, c='k')

    for arrow, alegend in zip(arrows, arrowlegend):
        if arrow != arrows[-1]:
            arstyle = '<-'
        else:
            arstyle = '-'
        if arrow == arrows[0]:
            lstyle = '-.'
        else:
            lstyle = '-.'

        ax_toplot.annotate(s='',
                           xy=(arrow[0][0], arrow[0][1]),
                           xytext=(arrow[2][0], arrow[2][1]),
                           arrowprops=dict(arrowstyle=arstyle,
                                           linestyle=lstyle,
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
    ax_toplot.set_ylim([-1.0, 2.0])
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
        'font.size': 19,
        'figure.figsize': (6 * 2, 2.75 * 3),
        'lines.linewidth': 6,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
    })
    data_array = np.array(data_array)

    gs_kw = dict(left=0.1,
                 right=0.95,
                 top=0.95,
                 bottom=0.1,
                 wspace=0.1,
                 hspace=0.1)

    gs_i = dict(left=0.2,
                right=0.95,
                top=0.95,
                bottom=0.1,
                wspace=0.1,
                hspace=0.1)

    fig1, axs = plt.subplots(nrows=3, ncols=2, gridspec_kw=gs_kw)
    fig2, axi = plt.subplots(nrows=1, ncols=1, gridspec_kw=gs_i)
    for axr1i, kinematic_listi, marker, iarrayi in zip(axs[:, 1],
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
    axr2 = axs[:, 0]
    for axr2i, datai in zip(axr2, data_array):
        axr2i.plot(datai[:, 0] / time_scale, datai[:, 1], label=labels[0])
        axr2i.plot(datai[:, 0] / time_scale, datai[:, 2], label=labels[1])
        axr2i.set_xlabel(r'$\^t$')
        if axr2i == axs[1, 0]:
            axr2i.set_ylabel('Normalized velocity')
        if axr2i == axs[0, 0]:
            axr2i.legend(fontsize='small', frameon=False)
        axr2i.axvline(x=1.0, color='k', linestyle='-', linewidth=0.5)
        axr2i.label_outer()
        if time_to_plot != 'all':
            axr2i.set_xlim(time_to_plot)
        if show_range != 'all':
            axr2i.set_ylim(show_range)

    #------at pt annotations--------------
    acc_t = 0.125
    ini_t_both = [0.02, 0.98 - acc_t]
    y1 = 1.05
    y2 = 1.15
    ymid = 0.5 * (y1 + y2)
    ytext = 1.25

    for ini_t in ini_t_both:
        vbar1 = [[ini_t, y1], [ini_t, y2]]
        vbar2 = [[ini_t + acc_t, y1], [ini_t + acc_t, y2]]
        vbars = vbar1 + vbar2
        arrow1 = [[ini_t, ymid], [ini_t + acc_t, ymid]]
        text_loc = [ini_t + 0.5 * acc_t, ytext]
        arrows = [arrow1]
        annotate_text = r'$\hat{t}_a$'

        nverts = len(vbars)
        codes = np.ones(nverts, int) * path.Path.LINETO
        codes[0::2] = path.Path.MOVETO
        vbarath = path.Path(vbars, codes)
        barpatch = patches.PathPatch(vbarath, edgecolor='k', linewidth=0.5)

        axr2[0].add_patch(barpatch)

        for arrow in arrows:
            axr2[0].annotate(s='',
                             xy=(arrow[0][0], arrow[0][1]),
                             xytext=(arrow[1][0], arrow[1][1]),
                             arrowprops=dict(arrowstyle='<-',
                                             facecolor='k',
                                             lw=0.5),
                             annotation_clip=False)
        axr2[0].annotate(s=annotate_text,
                         xy=(text_loc[0], text_loc[1]),
                         ha='center',
                         va='center',
                         annotation_clip=False)

    #--------------------pt2-------------
    pt2_both = [0.125, 0.24]
    for i in range(2):
        pt2 = pt2_both[i]
        ini_t_pt2 = 0.98 - pt2

        vbar1 = [[ini_t_pt2, y1], [ini_t_pt2, y2]]
        vbar2 = [[ini_t_pt2 + pt2, y1], [ini_t_pt2 + pt2, y2]]
        vbars = vbar1 + vbar2
        arrow1 = [[ini_t_pt2, ymid], [ini_t_pt2 + pt2, ymid]]
        text_loc = [ini_t_pt2 + 0.5 * pt2, ytext]
        arrows = [arrow1]
        annotate_text = r'$\hat{t}_p$'

        nverts = len(vbars)
        codes = np.ones(nverts, int) * path.Path.LINETO
        codes[0::2] = path.Path.MOVETO
        vbarath = path.Path(vbars, codes)
        barpatch = patches.PathPatch(vbarath, edgecolor='k', linewidth=0.5)

        axr2[i + 1].add_patch(barpatch)

        for arrow in arrows:
            axr2[i + 1].annotate(s='',
                                 xy=(arrow[0][0], arrow[0][1]),
                                 xytext=(arrow[1][0], arrow[1][1]),
                                 arrowprops=dict(arrowstyle='<-',
                                                 facecolor='k',
                                                 lw=0.5),
                                 annotation_clip=False)
        axr2[i + 1].annotate(s=annotate_text,
                             xy=(text_loc[0], text_loc[1]),
                             ha='center',
                             va='center',
                             annotation_clip=False)
    #---------------------------------------

    title = 'kinematics plot'
    out_image_file1 = os.path.join(image_out_path, title + '_case.svg')
    out_image_file2 = os.path.join(image_out_path, title + '_illustration.svg')

    #----add illustration figure-----
    illustrationk_plotter(axi, iarray, ifigure_parameters)

    fig1.savefig(out_image_file1)
    fig2.savefig(out_image_file2)
    # plt.show()

    return fig1, fig2
