"""circulation processing functions"""

import csv
import os
import shutil

import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from matplotlib import colors
from scipy.ndimage import zoom
from scipy.ndimage.filters import gaussian_filter


def read_sfield(field_data_file):
    """read field (vorticity or q) data"""
    vor_array = []
    with open(field_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                vor_array.append(
                    [-float(row[0]),
                     float(row[1]),
                     float(row[3])])
                line_count += 1

        print(f'Processed {line_count} lines in {field_data_file}')

    vor_array = np.array(vor_array)
    return vor_array


def read_wgeo(wgeo_data_file):
    """read wing geometry data"""
    wgeo_array = []
    with open(wgeo_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                wgeo_array.append([-float(row[0]), float(row[1])])
                line_count += 1

        print(f'Processed {line_count} lines in {wgeo_data_file}')

    wgeo_array = np.array(wgeo_array)

    # ----- sorting points in clockwise order ---
    x = wgeo_array[:, 0]
    y = wgeo_array[:, 1]
    cx = np.mean(x)
    cy = np.mean(y)
    a = np.arctan2(y - cy, x - cx)
    order = a.ravel().argsort()
    x = x[order]
    y = y[order]
    wgeo_array = np.vstack((x, y))

    wgeo_array = np.transpose(wgeo_array)
    w_centroid = np.array([cx, cy])

    return wgeo_array, w_centroid


def grid_vorz(window, resolution, vor_array):
    """grid interpolation for vorticity data"""
    grid_x, grid_y = np.mgrid[window[0]:window[1]:resolution[0] * 1j,
                              window[2]:window[3]:resolution[1] * 1j]

    grid_vz = scipy.interpolate.griddata(vor_array[:, 0:2],
                                         vor_array[:, 2], (grid_x, grid_y),
                                         method='nearest')

    return grid_x, grid_y, grid_vz


def single_plot_field(images, axto_plot, window, fdata, wdata, imnorm, levels):
    """plot one single field data"""
    images.append(
        axto_plot.imshow(fdata,
                         cmap='RdBu',
                         norm=imnorm,
                         aspect='equal',
                         extent=window,
                         origin='lower',
                         interpolation='bicubic'))
    axto_plot.contour(fdata,
                      levels,
                      linewidths=0.1,
                      colors='k',
                      extent=window,
                      origin='lower')

    nverts = len(wdata)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0] = path.Path.MOVETO
    codes[-1] = path.Path.CLOSEPOLY
    wgeopatch = path.Path(wdata, codes)
    patch = patches.PathPatch(wgeopatch,
                              linewidth=0.2,
                              facecolor='w',
                              edgecolor='k',
                              alpha=1.0)
    axto_plot.add_patch(patch)

    print('plotted image no = %s\n' % '{0:.0f}'.format(len(images)))

    return axto_plot


def field_plot(windows, field_data, wgeo_data, marks, oimage_file, mode):
    """plot field data"""
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (14, 12),
        'lines.linewidth': 0.5,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
    })
    marksc = marks[0]
    markpa = marks[1]
    no_r = len(markpa)
    no_c = len(marksc)
    imnorm = colors.Normalize(vmin=-100, vmax=100)
    levels = np.arange(-100.0, 100.0, 15)
    zoom_order = 5

    images = []
    ax_all = []
    #---plot t=0.75 data---
    gs_kw = dict(left=0.125,
                 right=0.9,
                 top=0.9,
                 bottom=0.72,
                 wspace=0.1,
                 hspace=0.0,
                 width_ratios=[6, 8, 10, 12])

    fig, axr1 = plt.subplots(nrows=1, ncols=no_c, gridspec_kw=gs_kw)
    for fdatai, wdatai, axr1i, window, marksci in zip(field_data, wgeo_data,
                                                      axr1, windows, marksc):
        fdatar1i = fdatai[0].T
        wdatar1i = wdatai[0]

        fdatar1i = zoom(fdatar1i, zoom_order)
        fdatar1i = gaussian_filter(fdatar1i, sigma=10.0)
        single_plot_field(images, axr1i, window, fdatar1i, wdatar1i, imnorm,
                          levels)

        axr1i.set_xticklabels([])
        axr1i.set_yticklabels([])

        markx_loc = axr1i.get_xlim()[0] + 0.5 * (axr1i.get_xlim()[1] -
                                                 axr1i.get_xlim()[0])
        marky_loc = axr1i.get_ylim()[1] + 0.1 * (axr1i.get_ylim()[1] -
                                                 axr1i.get_ylim()[0])
        axr1i.annotate(s=marksci,
                       xy=(markx_loc, marky_loc),
                       ha='center',
                       va='center',
                       annotation_clip=False)

        if marksci == marksc[0]:
            markx_loc = axr1i.get_xlim()[0] - 0.31 * (axr1i.get_xlim()[1] -
                                                      axr1i.get_xlim()[0])
            marky_loc = axr1i.get_ylim()[0] + 0.5 * (axr1i.get_ylim()[1] -
                                                     axr1i.get_ylim()[0])
            axr1i.annotate(s=r'$\^t$ = 0.75',
                           xy=(markx_loc, marky_loc),
                           ha='center',
                           va='center',
                           annotation_clip=False)

        ax_all.append(axr1i)

    #----plot different pa data----
    gs_i = fig.add_gridspec(nrows=no_r,
                            ncols=no_c,
                            left=0.125,
                            right=0.9,
                            top=0.7,
                            bottom=0.15,
                            wspace=0.1,
                            hspace=0.0,
                            width_ratios=[6, 8, 10, 12])

    for ci in range(no_c):
        for ri in range(no_r):
            axpa = fig.add_subplot(gs_i[ri, ci])
            fdatai = field_data[ci][ri + 1].T
            wdatai = wgeo_data[ci][ri + 1]

            fdatai = zoom(fdatai, zoom_order)
            fdatai = gaussian_filter(fdatai, sigma=10.0)
            single_plot_field(images, axpa, windows[ci], fdatai, wdatai,
                              imnorm, levels)

            axpa.set_xticklabels([])
            axpa.set_yticklabels([])
            if ci == no_c - 1:
                markx_loc = axpa.get_xlim()[1] + 0.12 * (axpa.get_xlim()[1] -
                                                         axpa.get_xlim()[0])
                marky_loc = axpa.get_ylim()[0] + 0.5 * (axpa.get_ylim()[1] -
                                                        axpa.get_ylim()[0])
                axpa.annotate(s=markpa[ri],
                              xy=(markx_loc, marky_loc),
                              ha='center',
                              va='center',
                              annotation_clip=False)

            if ci == 0:
                markx_loc = axpa.get_xlim()[0] - 0.31 * (axpa.get_xlim()[1] -
                                                         axpa.get_xlim()[0])
                marky_loc = axpa.get_ylim()[0] + 0.5 * (axpa.get_ylim()[1] -
                                                        axpa.get_ylim()[0])
                axpa.annotate(s=r'$\^t$ = 1.0',
                              xy=(markx_loc, marky_loc),
                              ha='center',
                              va='center',
                              annotation_clip=False)

            ax_all.append(axpa)

    cb = fig.colorbar(images[-1],
                      ax=ax_all,
                      orientation='horizontal',
                      fraction=0.1,
                      shrink=0.2,
                      pad=0.03)
    cb.ax.set_xlabel(r'$\omega$')
    cb.ax.xaxis.set_label_coords(-0.12, 1.5)

    if mode == 'save':
        plt.savefig(oimage_file)
    elif mode == 'show':
        plt.show()
