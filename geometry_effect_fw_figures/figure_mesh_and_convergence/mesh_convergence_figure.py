"""script for mesh and convergence figure plots"""

import os
from mcplotting_functions import read_cfd_data, cf_plotter
#-------------input plot control----------
cfd_data_list = [
    'ar5.0_ofs0.0_r1h0.5__Re1000.0_pt0.25_basic',
    'ar5.0_ofs0.0_r1h0.5__Re1000.0_pt0.25_sc',
    'ar5.0_ofs0.0_r1h0.5__Re1000.0_pt0.25_tc'
]
#-----------------------------------------
legends = [
    r'mesh1, $\Delta t$ = 1e-3',
    r'mesh2, $\Delta t$ = 1e-3',
    r'mesh1, $\Delta t$ = 5e-4',
]
#-----------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [4.0, 5.0]
show_range_cl = [-1.0, 5.0]
show_range_cd = [-8.0, 8.0]
# show_range_cd = [-3, 8]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
#-----------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/2_mesh_and_convergence')
image_out_path = cwd
#-----------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(data_dir, 'SIM_RESULTS_convergence', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
           cycle_time, 'against_t')