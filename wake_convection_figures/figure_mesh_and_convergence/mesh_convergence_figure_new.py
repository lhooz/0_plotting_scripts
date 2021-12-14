"""script for mesh and convergence figure plots"""

import os
from mcplotting_functions_new import read_cfd_data, cf_plotter, mesh_plotter
#-------------input plot control----------
cfd_data_list = [
    'Re1000.0_stroke6.0_acf0.25_pf0.25_pa90.0_normal',
    'Re1000.0_stroke6.0_acf0.25_pf0.25_pa90.0_corase',
    'Re1000.0_stroke6.0_acf0.25_pf0.25_pa90.0_fine',
    'Re1000.0_stroke6.0_acf0.25_pf0.25_pa90.0_SmallDt'
]
#-----------------------------------------
legends = [
    r'140k cells, $\Delta t = 5$e-4', r'64k cells, $\Delta t = 5$e-4',
    r'310k cells, $\Delta t = 5$e-4', r'140k cells, $\Delta t = 2$e-4'
]
# legends = [
# r'mesh1, $\Delta t$ = 1e-3',
# r'mesh1, $\Delta t$ = 5e-4'
# ]
#-----------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [0, 2.0]
show_range_cl = [-2.0, 9.0]
show_range_cd = show_range_cl
# show_range_cd = [-2.0, 6.5]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
#-----------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_convection_data/1_mesh_and_convergence')
image_out_path = cwd
#-----------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(data_dir, 'SIM_RESULTS_new', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
           cycle_time, 'against_t')

mesh_plotter(image_out_path)
