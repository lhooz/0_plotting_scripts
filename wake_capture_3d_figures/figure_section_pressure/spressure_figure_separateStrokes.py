"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from spressure_processing_finctions_separateStrokes import (spressure_plot, read_geop)

# --------data parameters------------
AR = 6
Re = 1000.0
stroke = 5
time_to_plot = [0.05, 0.1,
                0.15]  # ---time instances to plot within local cycle--

mag = 8
show_pRange = [-1.0 * mag, mag]
# -----------------------------------------------------
markt = [r'$\^t = ' + '{0:.1f}'.format(x) + '$' for x in time_to_plot]
# ---case dir setup---------------
cwd = os.getcwd()
data_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'threeD_wakeCapture/FLOW_DATA_' + '{0:.0f}'.format(stroke) + 'Stroke')
oimage_file = os.path.join(
    cwd, 'surfacePressure' + '_AR' + '{0:.1f}'.format(AR) + '_Re' +
    '{0:.0f}'.format(Re) + '_Stroke' + '{0:.0f}'.format(stroke))
# -----------------------------------

case_name = 'phi160.0__ar' + '{0:.1f}'.format(
    AR) + '_ofs0.0_r1h0.5__Re' + '{0:.1f}'.format(Re) + '_pt0.25'

allSurfaceP = []
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

            spressure_data = read_geop(geop_data_file)
            # -----write wing centroid history-----
            # centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]
            # spressure_data = [spressure_data[0], spressure_data[3]]
            # print(len(spressure_data))

            allSurfaceP.append(spressure_data)

# -------------plot all vortices-----------
spressure_plot(allSurfaceP, markt, oimage_file, show_pRange, AR, 'save')
plt.close()
# ----------------------------------------------
