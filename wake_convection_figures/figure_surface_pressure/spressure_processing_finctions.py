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


def read_geop(geop_data_file, Uref, pa, stroke):
    """read wing geometry data"""
    geop_array = []
    with open(geop_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                if stroke == 4.5:
                    pData = float(row[3]) / (0.5 * Uref**2)
                else:
                    pData = float(row[6]) / (0.5 * Uref**2)

                xData = -float(row[0])
                yData = float(row[1])

                angle = (-45 - pa) * np.pi / 180
                x_transform = xData * np.cos(angle) - yData * np.sin(angle)
                y_transform = xData * np.sin(angle) + yData * np.cos(angle)

                geop_array.append([x_transform, y_transform, pData])
                line_count += 1

        print(f'Processed {line_count} lines in {geop_data_file}')

    geop_array = np.array(geop_array)

    # ----- sorting points on each surface ---
    x = geop_array[:, 0]
    y = geop_array[:, 1]
    p = geop_array[:, 2]

    xMax = np.amax(x)
    xMin = np.amin(x)
    yMax = np.amax(y)
    yMin = np.amin(y)

    upper = []
    lower = []
    left = []
    right = []
    for xi, yi, pi in zip(x, y, p):
        if np.abs(xi - xMax) < 1e-5:
            right.append([xi, yi, pi])
        if np.abs(xi - xMin) < 1e-5:
            left.append([xi, yi, pi])
        if np.abs(yi - yMax) < 1e-5:
            upper.append([xi, yi, pi])
        if np.abs(yi - yMin) < 1e-5:
            lower.append([xi, yi, pi])

    upper = np.array(upper)
    lower = np.array(lower)
    left = np.array(left)
    right = np.array(right)
    # print(yMax, yMin)

    orderUp = upper[:, 0].ravel().argsort()
    orderedUp = upper[orderUp]

    orderLow = lower[:, 0].ravel().argsort()
    orderedLow = lower[orderLow]

    orderL = left[:, 1].ravel().argsort()
    orderedL = left[orderL]

    orderR = right[:, 1].ravel().argsort()
    orderedR = right[orderR]
    # print(orderedUp[:,0])
    #=========sorting in clockwise order============
    # cx = np.mean(x)
    # cy = np.mean(y)
    # a = np.arctan2(x - cx, y - cy)
    # order = a.ravel().argsort()
    # # print(a[order])
    # x = x[order]
    # y = y[order]
    # p = p[order]

    # w_centroid = np.array([cx, cy])
    # #------ sorting p in each surface----
    # x = np.append(x, x[0])
    # y = np.append(y, y[0])
    # p = np.append(p, p[0])
    # spressure_data = []
    # sip = [p[0]]
    # for i in range(len(p) - 2):
    # sip.append(p[i + 1])

    # tangend_diff = (y[i + 1] - y[i]) * (x[i + 2] - x[i + 1]) - (
    # y[i + 2] - y[i + 1]) * (x[i + 1] - x[i])
    # element_lsquare = (y[i + 1] - y[i])**2 + (x[i + 2] - x[i + 1])**2

    # if np.abs(tangend_diff / element_lsquare) > 0.5 or i == len(p) - 3:
    # spressure_data.append(sip)
    # sip = [p[i + 1]]

    # spressure_data = [spressure_data[-1] + spressure_data[0]
    # ] + spressure_data[1:-1]
    #==========================================================
    spressure_data = [
        orderedLow[:, 2], orderedUp[:, 2], orderedL[:, 2], orderedR[:, 2]
    ]

    return spressure_data


def spressure_plot(allSurfaceP, markt, oimage_file, show_pRange, mode):
    """plot field data"""
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (14, 12),
        'lines.linewidth': 3.0,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
    })
    allSurfaceP = np.array(allSurfaceP, dtype=object)

    nocol = len(allSurfaceP)
    norow = len(markt)

    fig, ax = plt.subplots(norow, nocol)
    for j in range(nocol):
        sorted_surfacep = allSurfaceP[j]
        for i in range(norow):
            x0 = np.linspace(0, 1, len(sorted_surfacep[i][0]))
            x1 = np.linspace(0, 1, len(sorted_surfacep[i][1]))
            ax[i][j].plot(x0, sorted_surfacep[i][0][:], label='Lower surface')
            ax[i][j].plot(x1, sorted_surfacep[i][1][:], label='Upper surface')
            ax[i][j].set_ylim(show_pRange)

            ax[i][j].set_xlabel(r'x/c')
            ax[i][j].set_ylabel(r'Cp')
            ax[i][j].legend()
            ax[i][j].axhline(y=0, linewidth=0.5, linestyle='-.', color='k')

    if mode == 'save':
        plt.savefig(oimage_file)
    elif mode == 'show':
        plt.show()
