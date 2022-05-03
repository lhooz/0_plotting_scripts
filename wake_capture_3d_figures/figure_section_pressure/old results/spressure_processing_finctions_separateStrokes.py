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


def read_geop(geop_data_file):
    """read wing geometry data"""
    geop_array = []
    with open(geop_data_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                pData = float(row[0])

                xData = float(row[2])
                yData = float(row[4])
                rData = float(row[3])

                geop_array.append([xData, yData, rData, pData])
                line_count += 1

        print(f'Processed {line_count} lines in {geop_data_file}')

    geop_array = np.array(geop_array)

    # ----- sorting points on each surface ---
    x = geop_array[:, 0]
    y = geop_array[:, 1]
    r = geop_array[:, 2]
    p = geop_array[:, 3]
    NoPoints = len(r)

    orderr = r.ravel().argsort()
    orderedx = x[orderr]
    orderedy = y[orderr]
    orderedr = r[orderr]
    orderedp = p[orderr]

    #----sorting different sections----
    section_all = []
    r_previous = r[0]
    seci = []
    tolerence = 1e-4
    for i in range(NoPoints):
        if np.abs(r[i] - r_previous) < tolerence:
            seci.append([x[i], y[i], p[i], r[i]])
        else:
            section_all.append(seci)
            seci = [[x[i], y[i], p[i], r[i]]]

        r_previous = r[i]
        # print(r_previous)

    section_all_sorted = []
    for sec in section_all:
        sec = np.array(sec)
        xMax = np.amax(sec[:, 0])
        xMin = np.amin(sec[:, 0])
        yMax = np.amax(sec[:, 1])
        yMin = np.amin(sec[:, 1])

        upper = []
        lower = []
        LE = []
        TE = []
        for xypi in sec:
            xi = xypi[0]
            yi = xypi[1]
            pi = xypi[2]
            ri = xypi[3]
            if np.abs(xi - xMax) < tolerence:
                TE.append([xi, yi, pi, ri])
            if np.abs(xi - xMin) < tolerence:
                LE.append([xi, yi, pi, ri])
            if np.abs(yi - yMax) < tolerence:
                upper.append([xi, yi, pi, ri])
            if np.abs(yi - yMin) < tolerence:
                lower.append([xi, yi, pi, ri])

        upper = np.array(upper)
        lower = np.array(lower)
        LE = np.array(LE)
        TE = np.array(TE)
        # print(yMax, yMin)

        orderUp = upper[:, 0].ravel().argsort()
        orderedUp = upper[orderUp]

        orderLow = lower[:, 0].ravel().argsort()
        orderedLow = lower[orderLow]

        orderL = LE[:, 1].ravel().argsort()
        orderedL = LE[orderL]

        orderR = TE[:, 1].ravel().argsort()
        orderedR = TE[orderR]
        #==========================================================
        spressure_data = [orderedLow, orderedUp, orderedL, orderedR]

        section_all_sorted.append(spressure_data)

    return section_all_sorted


def spressure_plot(allSurfaceP, markt, oimage_file, show_pRange, AR, mode):
    """plot field data"""
    plt.rcParams.update({
        # "text.usetex": True,
        'mathtext.fontset': 'stix',
        'font.family': 'STIXGeneral',
        'font.size': 18,
        'figure.figsize': (8, 15),
        'lines.linewidth': 5.0,
        'lines.markersize': 0.1,
        'lines.markerfacecolor': 'white',
        'figure.dpi': 100,
    })
    chord = 0.06  #---mean chord 0.06--
    R = AR * chord
    allSurfaceP = np.array(allSurfaceP, dtype=object)

    norow = len(markt)

    fig, ax = plt.subplots(norow)
    fig2, ax2 = plt.subplots(norow)
    for i in range(norow):
        sorted_surfacep = allSurfaceP[i]
        for j in range(len(sorted_surfacep)):
            # x0 = np.linspace(0, 1, len(seci[0]))
            # x1 = np.linspace(0, 1, len(seci[1]))
            # ax[i].plot(seci[0][:, 0], seci[0][:, 2], label='Lower surface')
            seci = sorted_surfacep[j]
            xMin = np.amin(seci[0][:, 0])
            xMax = np.amax(seci[0][:, 0])
            scale = xMax - xMin
            # print(seci[1][:, 3])
            x = (seci[0][:, 0] - seci[0][0, 0]) / scale
            r = seci[0][0, 3]

            if 0.1 * R < r < 3 * chord:
                ax[i].plot(x,
                           seci[0][:, 2],
                           label='r/c = ' + '{0:.2f}'.format(r / chord))
                ax[i].set_ylim(show_pRange)

                ax[i].set_xlabel(r'x/c')
                ax[i].set_ylabel(r'Cp')
                ax[i].label_outer()
                ax[i].axhline(y=0, linewidth=0.5, linestyle='-.', color='k')
            elif 3 * chord < r < 0.9 * R:
                ax2[i].plot(x,
                            seci[0][:, 2],
                            label='r/c = ' + '{0:.2f}'.format(r / chord))
                ax2[i].set_ylim(show_pRange)

                ax2[i].set_xlabel(r'x/c')
                ax2[i].set_ylabel(r'Cp')
                ax2[i].label_outer()
                ax2[i].axhline(y=0, linewidth=0.5, linestyle='-.', color='k')

    ax[0].legend(loc='upper center', ncol=3, fontsize='small', frameon=False)
    ax2[0].legend(loc='upper center', ncol=3, fontsize='small', frameon=False)
    if mode == 'save':
        fig.savefig(oimage_file + '_inboard.svg')
        fig2.savefig(oimage_file + '_outboard.svg')
    elif mode == 'show':
        plt.show()
