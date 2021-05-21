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


def read_vorfield(field_data_file, Uref):
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
                     float(row[1]), -float(row[3]) / Uref])
                line_count += 1

        print(f'Processed {line_count} lines in {field_data_file}')

    vor_array = np.array(vor_array)
    return vor_array


def read_qfield(field_data_file, Uref):
    """read field (vorticity or q) data"""
    refC = Uref * Uref
    q_array = []
    with open(field_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                q_array.append(
                    [-float(row[0]),
                     float(row[1]),
                     float(row[3]) / refC])
                line_count += 1

        print(f'Processed {line_count} lines in {field_data_file}')

    q_array = np.array(q_array)
    return q_array


def read_vfield(field_data_file, Uref):
    """read field (vorticity or q) data"""
    v_array = []
    with open(field_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                v_array.append([
                    -float(row[0]),
                    float(row[1]), -float(row[3]) / Uref,
                    float(row[4]) / Uref
                ])
                line_count += 1

        print(f'Processed {line_count} lines in {field_data_file}')

    v_array = np.array(v_array)
    return v_array


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


def grid_ufield(window, resolution, ufield_array):
    """grid interpolation for vorticity data"""
    grid_x, grid_y = np.mgrid[window[0]:window[1]:resolution[0] * 1j,
                              window[2]:window[3]:resolution[1] * 1j]

    grid_ux = scipy.interpolate.griddata(ufield_array[:, 0:2],
                                         ufield_array[:, 2], (grid_x, grid_y),
                                         method='nearest')

    grid_uy = scipy.interpolate.griddata(ufield_array[:, 0:2],
                                         ufield_array[:, 3], (grid_x, grid_y),
                                         method='nearest')

    return grid_x, grid_y, grid_ux, grid_uy


def single_plot_field(images, axto_plot, window, grid_x, grid_y, sImgdata,
                      sCtrdata, vdata, wdata, imnorm, levels):
    """plot one single field data"""
    images.append(
        axto_plot.imshow(
            sImgdata,
            # cmap='bwr',
            cmap='seismic',
            norm=imnorm,
            aspect='equal',
            extent=window,
            origin='lower',
            interpolation='bicubic'))
    axto_plot.contour(
        sCtrdata,
        levels,
        linewidths=1.5,
        # linestyles='--',
        colors='lime',
        extent=window,
        origin='lower')
    axto_plot.quiver(grid_x,
                     grid_y,
                     vdata[0],
                     vdata[1],
                     units='height',
                     scale=20,
                     width=0.0025,
                     headwidth=3.5)

    nverts = len(wdata)
    codes = np.ones(nverts, int) * path.Path.LINETO
    codes[0] = path.Path.MOVETO
    codes[-1] = path.Path.CLOSEPOLY
    wgeopatch = path.Path(wdata, codes)
    patch = patches.PathPatch(wgeopatch,
                              linewidth=1.5,
                              facecolor='k',
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
        'font.size': 14,
        'figure.figsize': (15, 11),
        'lines.linewidth': 4,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 300,
    })
    marksc = marks[0]
    markre = marks[1]
    no_r = len(markre)
    no_c = len(marksc)
    imnorm = colors.Normalize(vmin=-20, vmax=20)
    levels = [1.0]
    zoom_order = 6

    images = []
    ax_all = []
    #----plot different pa data----
    gs_kw = dict(left=0.2,
                 right=0.8,
                 top=0.8,
                 bottom=0.2,
                 wspace=0.05,
                 hspace=0.05,
                 width_ratios=[32, 40, 42])

    fig, ax = plt.subplots(nrows=no_r, ncols=no_c, gridspec_kw=gs_kw)
    for ci in range(no_c):
        axci = ax[:, ci]
        for ri in range(no_r):
            axre = axci[ri]
            grid_x = field_data[ci][ri][0]
            grid_y = field_data[ci][ri][1]
            sImgdatai = field_data[ci][ri][2]
            sCtrdatai = field_data[ci][ri][3]
            vdatai = field_data[ci][ri][4]

            wdatai = wgeo_data[ci][ri]

            sImgdatai = zoom(sImgdatai, zoom_order)
            sImgdatai = gaussian_filter(sImgdatai, sigma=5.0)
            sCtrdatai = zoom(sCtrdatai, zoom_order)
            sCtrdatai = gaussian_filter(sCtrdatai, sigma=5.0)
            single_plot_field(images, axre, windows[ci], grid_x, grid_y,
                              sImgdatai, sCtrdatai, vdatai, wdatai, imnorm,
                              levels)

            axre.set_xticklabels([])
            axre.set_yticklabels([])

            if ri == 0:
                markx_loc = axre.get_xlim()[0] + 0.5 * (axre.get_xlim()[1] -
                                                        axre.get_xlim()[0])
                marky_loc = axre.get_ylim()[1] + 0.1 * (axre.get_ylim()[1] -
                                                        axre.get_ylim()[0])
                axre.annotate(s=marksc[ci],
                              xy=(markx_loc, marky_loc),
                              ha='center',
                              va='center',
                              annotation_clip=False)

            if ci == no_c - 1:
                markx_loc = axre.get_xlim()[1] + 0.17 * (axre.get_xlim()[1] -
                                                         axre.get_xlim()[0])
                marky_loc = axre.get_ylim()[0] + 0.5 * (axre.get_ylim()[1] -
                                                        axre.get_ylim()[0])
                axre.annotate(s=markre[ri],
                              xy=(markx_loc, marky_loc),
                              ha='center',
                              va='center',
                              annotation_clip=False)

            ax_all.append(axre)

    cb = fig.colorbar(images[-1],
                      ax=ax_all,
                      orientation='horizontal',
                      fraction=0.1,
                      shrink=0.2,
                      pad=0.03)
    cb.ax.set_xlabel(r'$\omega$*')
    cb.ax.xaxis.set_label_coords(-0.12, 1.5)

    if mode == 'save':
        plt.savefig(oimage_file)
    elif mode == 'show':
        plt.show()
