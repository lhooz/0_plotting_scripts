"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from field_processing_finctions import (field_plot, grid_vorz, grid_ufield,
                                        read_sfield, read_vfield, read_wgeo)

#--------------figure parameters------------
# window_sc3 = [-1.0, 6.0, -3.0, 3.0]
# window_sc5 = [-1.0, 8.0, -3.0, 3.0]
# window_sc7 = [-1.0, 10.0, -3.0, 3.0]
# window_sc9 = [-1.0, 12.0, -3.0, 3.0]
#-------------------------------------------
window_sc = [-1.0, 12.0, -3.0, 3.0]
resolution = [400, 200]  #--interpolate res for vorticity field--
resolution_vector = [60, 30]  #--interpolate res for velocity field--
data_time_increment = 0.1
#--------data parameters------------
Re = 1400.0
stroke = 9.0
#-----------------------------------
#---case dir setup--
window = window_sc
#-----------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_capture_2d_reciprocal_data/FIELD_RESULTS')
oimage_file = os.path.join(cwd, 'flow_visulization.png')
#-----------------------------------
case_name = 'Re' + '{0:.1f}'.format(Re) + '_stroke' + '{0:.1f}'.format(
    stroke) + '_acf0.25_pf0.25_pa90.0'

case_dir = os.path.join(data_dir, case_name)
vorz_folder = os.path.join(case_dir, 'vorz_data')
q_folder = os.path.join(case_dir, 'q_data')
ufield_folder = os.path.join(case_dir, 'ufield_data')
wgeo_folder = os.path.join(case_dir, 'geop_data')

#----------------------------------------
time_series_names = [f.name for f in os.scandir(vorz_folder) if f.is_file()]
time_series = [x.split('_')[-1] for x in time_series_names]
time_series = [int(x.split('.')[0]) for x in time_series]
start_t = np.min(np.array(time_series))
end_t = np.max(np.array(time_series))

ufield_gridx = []
ufield_gridy = []
sImgdata = []
sCtrdata = []
vdata = []
wdata = []
markt = []
for ti in range(start_t, end_t + 1):
    timei = (ti + 1) * data_time_increment
    if 0.0 <= timei <= 1.05:
        time_instance = str(ti).zfill(4)
        #--------setting up file dirs-----------
        vorz_data_file = os.path.join(vorz_folder,
                                      'vorz_' + time_instance + '.csv')
        q_data_file = os.path.join(q_folder, 'q_' + time_instance + '.csv')
        ufield_data_file = os.path.join(ufield_folder,
                                        'ufield_' + time_instance + '.csv')
        geop_data_file = os.path.join(wgeo_folder,
                                      'geop_' + time_instance + '.csv')

        vorz_data = read_sfield(vorz_data_file)
        q_data = read_sfield(q_data_file)
        ufield_data = read_vfield(ufield_data_file)
        wgeo_data, w_centroid = read_wgeo(geop_data_file)
        #-------------------------------------
        grid_x, grid_y, grid_vz = grid_vorz(window, resolution, vorz_data)
        grid_x, grid_y, grid_q = grid_vorz(window, resolution, q_data)
        grid_x, grid_y, grid_ux, grid_uy = grid_ufield(window,
                                                       resolution_vector,
                                                       ufield_data)
        #-----write wing centroid history-----
        centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]

        ufield_gridx.append(grid_x.T)
        ufield_gridy.append(grid_y.T)
        sImgdata.append(grid_vz.T)
        sCtrdata.append(grid_q.T)
        vdata.append([grid_ux.T, grid_uy.T])
        wdata.append(wgeo_data)
        markt.append(timei)
# -------------plot all vortices-----------
field_plot(window, ufield_gridx, ufield_gridy, sImgdata, sCtrdata, vdata,
           wdata, markt, oimage_file, 'save')
plt.close()
# ----------------------------------------------
