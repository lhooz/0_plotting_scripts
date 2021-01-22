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
from scipy import ndimage
from scipy.interpolate import UnivariateSpline


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
                vor_array.append([float(row[0]), float(row[1]), float(row[3])])
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
                wgeo_array.append([float(row[0]), float(row[1])])
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


def field_plot(window, field_data, wgeo_data, oimage_file, mode):
    """plot field data"""
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (10, 14),
        'lines.linewidth': 1.0,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 200,
    })
    no_subfigs = len(field_data)
    imnorm = colors.Normalize(vmin=-100, vmax=100)

    images = []
    fig, axs = plt.subplots(no_subfigs, 1)
    for i in range(no_subfigs):
        images.append(axs[i].imshow(field_data[i].T,
                                    cmap='RdBu',
                                    norm=imnorm,
                                    aspect='equal',
                                    extent=(window[0], window[1], window[2],
                                            window[3]),
                                    origin='lower'))
        axs[i].label_outer()

        nverts = len(wgeo_data[i])
        codes = np.ones(nverts, int) * path.Path.LINETO
        codes[0] = path.Path.MOVETO
        codes[-1] = path.Path.CLOSEPOLY
        wgeopatch = path.Path(wgeo_data[i], codes)
        patch = patches.PathPatch(wgeopatch,
                                  linewidth=0.5,
                                  facecolor='w',
                                  edgecolor='k',
                                  alpha=1.0)
        axs[i].add_patch(patch)

    fig.colorbar(images[0],
                 ax=axs,
                 orientation='horizontal',
                 fraction=.2,
                 shrink=0.2,
                 pad=0.05)

    if mode == 'save':
        plt.savefig(oimage_file)
    elif mode == 'show':
        plt.show()
