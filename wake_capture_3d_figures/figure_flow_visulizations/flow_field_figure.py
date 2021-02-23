"""main script for circulation processing"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from field_processing_finctions import (field_plot, grid_vorz, grid_ufield,
                                        read_field, read_wgeo)

#--------------figure parameters------------
window_sc = [-0.25, 0.25, -0.25, 0.15]
resolution = [250, 200]  #--interpolate res for vorticity field--
resolution_vector = [25, 20]  #--interpolate res for velocity field--
#--------data parameters------------
ar = 5.0
ofs = 0.0
r1h = 0.5
re = 1000.0
p = 0.25

section = 0.8
#-----------------------------------
#---case dir setup--
window = window_sc
#-----------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_capture_3d_data/FIELD_RESULTS')
oimage_file = os.path.join(cwd, 'flow_visulization.png')
#-----------------------------------
case_name = 'ar' + '{0:.1f}'.format(ar) + '_ofs' + '{0:.1f}'.format(
    ofs) + '_r1h' + '{0:.1f}'.format(r1h) + '__Re' + '{0:.1f}'.format(
        re) + '_pt' + '{0:.3g}'.format(p)

case_dir = os.path.join(data_dir, case_name,
                        'section_' + '{0:.2f}'.format(section))
field_folder = os.path.join(case_dir, 'field_data')
geop_folder = os.path.join(case_dir, 'geop_data')
#----------------------------------------
time_series_names = [f.name for f in os.scandir(field_folder) if f.is_file()]
time_series = [x.split('_')[-1] for x in time_series_names]
time_series = np.sort(
    np.array([float(x.split('.csv')[0]) for x in time_series]))

ufield_gridx = []
ufield_gridy = []
sdata = []
vdata = []
wdata = []
markt = []
for timei in time_series:
    if 4.6 <= timei <= 5.0:
        time_instance = '{0:.2f}'.format(timei)
        #--------setting up file dirs-----------
        field_data_file = os.path.join(field_folder,
                                       'field_' + time_instance + '.csv')
        geop_data_file = os.path.join(geop_folder,
                                      'geop_' + time_instance + '.csv')

        vor_data, q_data, ufield_data = read_field(field_data_file)
        wgeo_data, w_centroid = read_wgeo(geop_data_file)
        #-------------------------------------
        vorz_data = vor_data[:, [0, 1, 4]]
        uxyfield_data = ufield_data[:, 0:4]
        grid_x, grid_y, grid_vz = grid_vorz(window, resolution, vorz_data)
        grid_x, grid_y, grid_ux, grid_uy = grid_ufield(window,
                                                       resolution_vector,
                                                       uxyfield_data)
        # grid_x, grid_y, grid_q = grid_vorz(window, resolution, q_data)
        #-----write wing centroid history-----
        centroid_ti = [str(timei), str(w_centroid[0]), str(w_centroid[1])]

        ufield_gridx.append(grid_x.T)
        ufield_gridy.append(grid_y.T)
        sdata.append(grid_vz.T)
        vdata.append([grid_ux.T, grid_uy.T])
        wdata.append(wgeo_data)
        markt.append(timei)
# -------------plot all vortices-----------
field_plot(window, ufield_gridx, ufield_gridy, sdata, vdata, wdata, markt,
           oimage_file, 'save')
plt.close()
# ----------------------------------------------
