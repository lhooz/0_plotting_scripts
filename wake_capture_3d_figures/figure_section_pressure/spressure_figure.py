"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from spressure_processing_finctions import (spressure_plot, read_geop,
                                            Cn_image)

# --------data parameters------------
AR = 6
Re = 1000.0
# ---time instances to plot within local cycle--
time_to_plot = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

mag = 3
show_pRange = [-1.0 * mag, mag]
# -----------------------------------------------------
markt = [r'$\^t = ' + '{0:.2f}'.format(x) + '$' for x in time_to_plot]
# ---case dir setup---------------
cwd = os.getcwd()
data_dir0 = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'threeD_wakeCapture/FLOW_DATA_')
oimage_file = os.path.join(
    cwd, 'Cn' + '_AR' + '{0:.1f}'.format(AR) + '_Re' +
    '{0:.0f}'.format(Re))
odata_file = os.path.join(
    cwd,
    'sectionCn' + '_AR' + '{0:.1f}'.format(AR) + '_Re' + '{0:.0f}'.format(Re))
# -----------------------------------

case_name = 'phi160.0__ar' + '{0:.1f}'.format(
    AR) + '_ofs0.0_r1h0.5__Re' + '{0:.1f}'.format(Re) + '_pt0.25'

stroke = 1
data_dir = data_dir0 + '{0:.0f}'.format(stroke) + 'Stroke'
allSurfaceP1Stroke = []
allSectionCn1Stroke = []
for f in os.scandir(data_dir):
    if f.name.startswith(case_name):
        case_dir = os.path.join(data_dir, f.name)
        geop_folder = os.path.join(case_dir, 'processed_sec_data')
        # ----------------------------------------
        for timei in time_to_plot:
            # --------setting up file dirs-----------
            time_instance = '{0:.2f}'.format(stroke + timei - 1)
            geop_data_file = os.path.join(geop_folder,
                                          'geop_' + time_instance + '.csv')

            spressure_data, sectionCn, sec_r_c = read_geop(geop_data_file)
            # -----write wing centroid history-----
            # centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]
            # spressure_data = [spressure_data[0], spressure_data[3]]
            # print(len(spressure_data))

            allSurfaceP1Stroke.append(spressure_data)
            allSectionCn1Stroke.append(sectionCn)

stroke = 5
data_dir = data_dir0 + '{0:.0f}'.format(stroke) + 'Stroke'
allSurfaceP5Stroke = []
allSectionCn5Stroke = []
for f in os.scandir(data_dir):
    if f.name.startswith(case_name):
        case_dir = os.path.join(data_dir, f.name)
        geop_folder = os.path.join(case_dir, 'processed_sec_data')
        # ----------------------------------------
        for timei in time_to_plot:
            # --------setting up file dirs-----------
            time_instance = '{0:.2f}'.format(stroke + timei - 1)
            geop_data_file = os.path.join(geop_folder,
                                          'geop_' + time_instance + '.csv')

            spressure_data, sectionCn, sec_r_c = read_geop(geop_data_file)
            # -----write wing centroid history-----
            # centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]
            # spressure_data = [spressure_data[0], spressure_data[3]]
            # print(len(spressure_data))

            allSurfaceP5Stroke.append(spressure_data)
            allSectionCn5Stroke.append(sectionCn)

allSectionCn1Stroke = np.array(allSectionCn1Stroke)
allSectionCn5Stroke = np.array(allSectionCn5Stroke)

allSurfaceP = [allSurfaceP1Stroke, allSurfaceP5Stroke]
allSectionCn = [
    allSectionCn1Stroke, allSectionCn5Stroke,
    allSectionCn5Stroke - allSectionCn1Stroke
]

with open(odata_file + '.csv', 'w') as f:
    f.write("z - Cn\n")
    f.write("y - t_hat\n")
    f.write("x - r/c:," + ','.join([str(x)
                                    for x in sec_r_c]) + "\n")
    for st, datai in zip(['1', '5', 'wake'], allSectionCn):
        f.write("stroke = " + st + ":\n")
        for data_ti, ti in zip(datai, time_to_plot):
            f.write("%s," % '{0:.2f}'.format(ti))
            f.write(','.join([str(x) for x in data_ti]) + "\n")
# -------------plot all pressure-----------
# spressure_plot(allSurfaceP, markt, oimage_file, show_pRange, AR, 'save')
# plt.close()
Cn_image(time_to_plot, sec_r_c, allSectionCn, oimage_file)
# ----------------------------------------------
