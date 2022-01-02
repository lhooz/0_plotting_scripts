"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from spressure_processing_finctions import (spressure_plot, read_geop)

#--------data parameters------------
pa = 90
Re = 1000.0
time_to_plot = [1.1, 1.2, 1.4]  #---time instances to plot--
stroke = [1.5, 3.0, 4.5, 6.0]
if pa == 0:
    if Re == 100.0:
        mag = 3
    else:
        mag = 3
if pa == 45:
    if Re == 100.0:
        mag = 6
    else:
        mag = 6
if pa == 90:
    if Re == 100.0:
        mag = 12
    else:
        mag = 12

show_pRange = [-1 * mag, mag]
#-----------------------------------
# Uref = [1.786, 3.571, 5.357, 7.143]  #--ref velocity--
Uref = [2.086, 3.571, 5.357, 6.143]  #--ref velocity--
data_time_increment = 0.1
#-----------------------------------------------------
markt = [r'$\^t = ' + '{0:.1f}'.format(x) + '$' for x in time_to_plot]
#---case dir setup---------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_convection_data/4_FIELD_RESULTS')
oimage_file = os.path.join(
    cwd, 'surfacePressure' + '_pa' + '{0:.0f}'.format(pa) + '_Re' +
    '{0:.0f}'.format(Re) + '.svg')
#-----------------------------------
if pa == 45.0:
    pf = 0.125
else:
    pf = 0.25

allSurfaceP = []
for st, Urefi in zip(stroke, Uref):
    case_name = 'Re' + '{0:.1f}'.format(Re) + '_stroke' + '{0:.1f}'.format(
        st) + '_acf0.25_pf' + '{0:.3g}'.format(pf) + '_pa' + '{0:.1f}'.format(
            pa)

    case_dir = os.path.join(data_dir, case_name)
    geop_folder = os.path.join(case_dir, 'geop_data')
    #----------------------------------------
    time_series_names = [
        f.name for f in os.scandir(geop_folder) if f.is_file()
    ]
    time_series = [x.split('_')[-1] for x in time_series_names]
    time_series = [int(x.split('.')[0]) for x in time_series]
    start_t = np.min(np.array(time_series))
    end_t = np.max(np.array(time_series))

    sorted_surfacep = []
    for ti in range(start_t, end_t + 1):
        timei = (ti + 1) * data_time_increment
        timei = round(timei, 1)
        if timei in time_to_plot:
            time_instance = str(ti).zfill(4)
            #--------setting up file dirs-----------
            geop_data_file = os.path.join(geop_folder,
                                          'geop_' + time_instance + '.csv')

            spressure_data = read_geop(geop_data_file, Urefi, pa, st)
            #-----write wing centroid history-----
            # centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]
            # spressure_data = [spressure_data[0], spressure_data[3]]
            print(len(spressure_data))

            sorted_surfacep.append(spressure_data)

    allSurfaceP.append(sorted_surfacep)
# -------------plot all vortices-----------
spressure_plot(allSurfaceP, markt, oimage_file, show_pRange, 'save')
plt.close()
# ----------------------------------------------
