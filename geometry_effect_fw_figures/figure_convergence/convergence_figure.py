"""script for mesh and convergence figure plots"""

import os
from cplotting_functions import read_cfd_data, cf_plotter
#-------------input plot control----------
kinematics_file = 'kinematics_convergence.dat'
cfd_data_list = [
    'ar3.0_ofs0.0_r1h0.5__Re100.0_pt0.25_Ro1.68_basic',
    'ar3.0_ofs0.0_r1h0.5__Re100.0_pt0.25_Ro1.68_sc',
    'ar3.0_ofs0.0_r1h0.5__Re100.0_pt0.25_Ro1.68_tc',
    'ar3.0_ofs0.0_r1h0.5__Re100.0_pt0.25_Ro1.68_sctc',
]
#-----------------------------------------
legends = [
    r'3 million, $\Delta t$ = 1e-3',
    r'6 million, $\Delta t$ = 1e-3',
    r'3 million, $\Delta t$ = 5e-4',
    r'6 million, $\Delta t$ = 5e-4',
]
#-----------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [4.0, 5.0]
show_range_cl = [-0.2, 2.2]
show_range_cd = [-1.2, 2.7]
# show_range_cd = [-1.0, 5.0]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
#-----------------------------------------
cwd = os.getcwd()
kinematics_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'geometry_effect_fw_data/1_kinematic_cases')
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/2_mesh_and_convergence')
image_out_path = cwd
#-----------------------------------------
kinematics_data = os.path.join(kinematics_dir, kinematics_file)
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(data_dir, 'SIM_RESULTS_convergence', cfi)
    cf_arrayi = read_cfd_data(kinematics_data, cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
           cycle_time, 'against_t')
