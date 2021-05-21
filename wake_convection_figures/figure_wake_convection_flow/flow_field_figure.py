"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from field_processing_finctions import (field_plot, grid_vorz, grid_ufield,
                                        read_vorfield, read_qfield,
                                        read_vfield, read_wgeo)

#--------------figure parameters------------
#----pa = 0 window--
# window_sc15 = [-0.6, 2.6, -1.8, 1.8]
# window_sc30 = [-0.1, 3.9, -1.8, 1.8]
# window_sc60 = [2.7, 6.9, -1.8, 1.8]
#---pa = 45 window--
window_sc15 = [-0.1, 3.1, -2.0, 1.6]
window_sc30 = [0.4, 4.4, -2.0, 1.6]
window_sc60 = [3.2, 7.4, -2.0, 1.6]
#------------------------------------------
resolution = [250, 220]  #--interpolate res for vorticity field--
#--interpolate res for velocity field--
res_vector_sc15 = [20, 22]
res_vector_sc30 = [25, 22]
res_vector_sc60 = [26, 22]
resolution_vector = [res_vector_sc15, res_vector_sc30, res_vector_sc60]
time_to_plot = [1.2, 1.4, 1.6]  #---time instances to plot--
#--------------------------------------
Uref = [1.786, 3.571, 7.143]  #--ref velocity--
data_time_increment = 0.1
#--------data parameters------------
Re = 1000.0
stroke = [1.5, 3.0, 6.0]
pa = 45
#-----------------------------------
windows = [window_sc15, window_sc30, window_sc60]
marksc = [r'$s/c = ' + '{0:.1f}'.format(x) + '$' for x in stroke]
markt = [r'$\^t = ' + '{0:.1f}'.format(x) + '$' for x in time_to_plot]
marks = [marksc, markt]
#-----------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_convection_data/4_FIELD_RESULTS')
oimage_file = os.path.join(
    cwd, 'wake_flow_pa' + '{0:.0f}'.format(pa) + '_Re' + '{0:.0f}'.format(Re) +
    '.svg')
#-----------------------------------
fdata_all = []
wdata_all = []
for si, window, Urefi, res_vectori in zip(stroke, windows, Uref,
                                          resolution_vector):
    #------read data for each time-------------
    if pa == 45.0:
        pf = 0.125
    else:
        pf = 0.25

    case_name_t10 = 'Re' + '{0:.1f}'.format(Re) + '_stroke' + '{0:.1f}'.format(
        si) + '_acf0.25_pf' + '{0:.3g}'.format(pf) + '_pa' + '{0:.1f}'.format(
            pa)

    case_dir = os.path.join(data_dir, case_name_t10)
    vorz_folder = os.path.join(case_dir, 'vorz_data')
    q_folder = os.path.join(case_dir, 'q_data')
    ufield_folder = os.path.join(case_dir, 'ufield_data')
    wgeo_folder = os.path.join(case_dir, 'geop_data')

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
        timei = round(timei, 1)
        if timei in time_to_plot:
            time_instance = str(ti).zfill(4)
            #--------setting up file dirs-----------
            vorz_data_file = os.path.join(vorz_folder,
                                          'vorz_' + time_instance + '.csv')
            q_data_file = os.path.join(q_folder, 'q_' + time_instance + '.csv')
            ufield_data_file = os.path.join(ufield_folder,
                                            'ufield_' + time_instance + '.csv')
            wgeo_data_file = os.path.join(wgeo_folder,
                                          'geop_' + time_instance + '.csv')

            vorz_data = read_vorfield(vorz_data_file, Urefi)
            q_data = read_qfield(q_data_file, Urefi)
            ufield_data = read_vfield(ufield_data_file, Urefi)
            wgeo_data, w_centroid = read_wgeo(wgeo_data_file)
            #-------------------------------------
            grid_x, grid_y, grid_vz = grid_vorz(window, resolution, vorz_data)
            grid_x, grid_y, grid_q = grid_vorz(window, resolution, q_data)
            grid_x, grid_y, grid_ux, grid_uy = grid_ufield(
                window, res_vectori, ufield_data)
            #-----write wing centroid history-----
            centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]

            ufield_gridx = grid_x.T
            ufield_gridy = grid_y.T
            sImgdata = grid_vz.T
            sCtrdata = grid_q.T
            vdata = [grid_ux.T, grid_uy.T]

            fdata_si.append(
                [ufield_gridx, ufield_gridy, sImgdata, sCtrdata, vdata])
            wdata_si.append(wgeo_data)
        #----------------------------------------
    fdata_all.append(fdata_si)
    wdata_all.append(wdata_si)
# -------------plot all vortices-----------
field_plot(windows, fdata_all, wdata_all, marks, oimage_file, 'save')
plt.close()
# ----------------------------------------------
