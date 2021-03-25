"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from field_processing_finctions import (field_plot, grid_vorz, read_sfield,
                                        read_wgeo)

#--------------figure parameters------------
window_sc3 = [-1.0, 5.0, -3.0, 3.0]
window_sc5 = [-1.0, 7.0, -3.0, 3.0]
window_sc7 = [-1.0, 9.0, -3.0, 3.0]
window_sc9 = [-1.0, 11.0, -3.0, 3.0]
resolution = [400, 200]
data_time_increment = 0.25
#--------data parameters------------
Re = 1000.0
stroke = [3.0, 5.0, 7.0, 9.0]
pa = [0.0, 45.0, 90.0]
#-----------------------------------
windows = [window_sc3, window_sc5, window_sc7, window_sc9]
marksc = [r'$s/c = ' + '{0:.1f}'.format(x) + '$' for x in stroke]
markpa = [r'$\theta = ' + '{0:.0f}'.format(x) + '$' for x in pa]
marks = [marksc, markpa]
#-----------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_capture_2d_journal_data/4_FIELD_RESULTS')
oimage_file = os.path.join(cwd, 'flow_visulization.png')
#-----------------------------------
fdata_all = []
wdata_all = []
for si, window in zip(stroke, windows):
    #------read data for t=0.75 figures (1st row)-------
    case_name_t75 = 'Re' + '{0:.1f}'.format(Re) + '_stroke' + '{0:.1f}'.format(
        si) + '_acf0.25_pf0.25_pa90.0'

    case_dir = os.path.join(data_dir, case_name_t75)
    vorz_folder = os.path.join(case_dir, 'vorz_data')
    q_folder = os.path.join(case_dir, 'q_data')
    wgeo_folder = os.path.join(case_dir, 'wgeo_data')

    #----------------------------------------
    time_series_names = [
        f.name for f in os.scandir(vorz_folder) if f.is_file()
    ]
    time_series = [x.split('_')[-1] for x in time_series_names]
    time_series = [int(x.split('.')[0]) for x in time_series]
    start_t = np.min(np.array(time_series))
    end_t = np.max(np.array(time_series))

    fdata_si = []
    wdata_si = []
    for ti in range(start_t, end_t + 1):
        timei = (ti + 1) * data_time_increment
        if 0.7 <= timei <= 0.8:
            time_instance = str(ti).zfill(4)
            #--------setting up file dirs-----------
            vorz_data_file = os.path.join(vorz_folder,
                                          'vorz_' + time_instance + '.csv')
            q_data_file = os.path.join(q_folder, 'q_' + time_instance + '.csv')
            wgeo_data_file = os.path.join(wgeo_folder,
                                          'wgeo_' + time_instance + '.csv')

            vorz_data = read_sfield(vorz_data_file)
            # q_data = read_sfield(q_data_file)
            wgeo_data, w_centroid = read_wgeo(wgeo_data_file)
            #-------------------------------------
            grid_x, grid_y, grid_vz = grid_vorz(window, resolution, vorz_data)
            # grid_x, grid_y, grid_q = grid_vorz(window, resolution, q_data)
            #-----write wing centroid history-----
            centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]

            fdata_si.append(grid_vz)
            wdata_si.append(wgeo_data)
    #----------------------------------------
    #------read data for each pa-------------
    for pai in pa:
        if pai == 45.0:
            pf = 0.125
        else:
            pf = 0.25

        case_name_t10 = 'Re' + '{0:.1f}'.format(
            Re) + '_stroke' + '{0:.1f}'.format(
                si) + '_acf0.25_pf' + '{0:.3g}'.format(
                    pf) + '_pa' + '{0:.1f}'.format(pai)

        case_dir = os.path.join(data_dir, case_name_t10)
        vorz_folder = os.path.join(case_dir, 'vorz_data')
        q_folder = os.path.join(case_dir, 'q_data')
        wgeo_folder = os.path.join(case_dir, 'wgeo_data')

        #----------------------------------------
        time_series_names = [
            f.name for f in os.scandir(vorz_folder) if f.is_file()
        ]
        time_series = [x.split('_')[-1] for x in time_series_names]
        time_series = [int(x.split('.')[0]) for x in time_series]
        start_t = np.min(np.array(time_series))
        end_t = np.max(np.array(time_series))

        for ti in range(start_t, end_t + 1):
            timei = (ti + 1) * data_time_increment
            if 0.95 <= timei <= 1.05:
                time_instance = str(ti).zfill(4)
                #--------setting up file dirs-----------
                vorz_data_file = os.path.join(vorz_folder,
                                              'vorz_' + time_instance + '.csv')
                q_data_file = os.path.join(q_folder,
                                           'q_' + time_instance + '.csv')
                wgeo_data_file = os.path.join(wgeo_folder,
                                              'wgeo_' + time_instance + '.csv')

                vorz_data = read_sfield(vorz_data_file)
                # q_data = read_sfield(q_data_file)
                wgeo_data, w_centroid = read_wgeo(wgeo_data_file)
                #-------------------------------------
                grid_x, grid_y, grid_vz = grid_vorz(window, resolution,
                                                    vorz_data)
                # grid_x, grid_y, grid_q = grid_vorz(window, resolution, q_data)
                #-----write wing centroid history-----
                centroid_ti = [
                    str(timei),
                    str(w_centroid[0]),
                    str(w_centroid[1])
                ]

                fdata_si.append(grid_vz)
                wdata_si.append(wgeo_data)
        #----------------------------------------
    fdata_all.append(fdata_si)
    wdata_all.append(wdata_si)
# -------------plot all vortices-----------
field_plot(windows, fdata_all, wdata_all, marks, oimage_file, 'save')
plt.close()
# ----------------------------------------------
