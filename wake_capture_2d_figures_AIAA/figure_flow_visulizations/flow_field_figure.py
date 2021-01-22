"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from field_processing_finctions import field_plot, read_sfield, read_wgeo, grid_vorz

#--------------input parameters------------
window = [-9, -2, -2.5, 2.5]
resolution = [3500, 2500]
data_time_increment = 1e-1
#------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(cwd), '5_field_results')
case_dir = os.path.join(data_dir,
                        'FLOW_FIELDS_pa90/Re1000.0_stroke7.0_acf0.25_pf0.25')
vorz_folder = os.path.join(case_dir, 'vorz_data')
q_folder = os.path.join(case_dir, 'q_data')
wgeo_folder = os.path.join(case_dir, 'wgeo_data')

oimage_folder = os.path.join(cwd, 'oimages')
oimage_file = os.path.join(oimage_folder,
                           'Re1000.0_stroke7.0_acf0.25_pf0.25_pa90.png')
#----------------------------------------
if os.path.exists(oimage_folder):
    shutil.rmtree(oimage_folder)
os.mkdir(oimage_folder)
#----------------------------------------
time_series_names = [f.name for f in os.scandir(vorz_folder) if f.is_file()]
time_series = [x.split('_')[-1] for x in time_series_names]
time_series = [int(x.split('.')[0]) for x in time_series]
start_t = np.min(np.array(time_series))
end_t = np.max(np.array(time_series))

fdata_all = []
wdata_all = []
for ti in range(start_t, end_t + 1):
    timei = (ti + 1) * data_time_increment
    if 0.8 <= timei <= 1.1:
        time_instance = str(ti).zfill(4)
        #--------setting up file dirs-----------
        vorz_data_file = os.path.join(vorz_folder,
                                      'vorz_' + time_instance + '.csv')
        q_data_file = os.path.join(q_folder, 'q_' + time_instance + '.csv')
        wgeo_data_file = os.path.join(wgeo_folder,
                                      'wgeo_' + time_instance + '.csv')

        vorz_data = read_sfield(vorz_data_file)
        q_data = read_sfield(q_data_file)
        wgeo_data, w_centroid = read_wgeo(wgeo_data_file)
        #-------------------------------------
        grid_x, grid_y, grid_vz = grid_vorz(window, resolution, vorz_data)
        grid_x, grid_y, grid_q = grid_vorz(window, resolution, q_data)
        #-----write wing centroid history-----
        centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]

        fdata_all.append(grid_vz)
        wdata_all.append(wgeo_data)

# -------------plot all vortices-----------
field_plot(window, fdata_all, wdata_all, oimage_file, 'save')
plt.close()
# ----------------------------------------------
