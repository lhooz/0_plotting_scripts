"""script for mesh and convergence figure plots"""

import os

from tplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
ar = 7.5
#-----------------------------------------
kinematics_file = 'kinematics.dat'
cfd_data_list = [
    '_ofs0.0_r1h0.4__Re100.0_pt0.25',
    '_ofs0.0_r1h0.5__Re100.0_pt0.25',
    '_ofs0.0_r1h0.6__Re100.0_pt0.25',
]
#-----------------------------------------
legends = [
    r'$\hat r_1$ = 0.4',
    r'$\hat r_1$ = 0.5',
    r'$\hat r_1$ = 0.6',
]
#-----------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [4.0, 5.0]
if ar == 1.5:
    show_range_cl = [-1.3, 6.0]
    show_range_cd = [-1.3, 6.0]
elif ar == 3.0:
    show_range_cl = [-1., 3.0]
    show_range_cd = [-1., 3.0]
elif ar == 7.5:
    show_range_cl = [-1., 3.0]
    show_range_cd = [-1., 3.0]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
#-----------------------------------------
cwd = os.getcwd()
kinematics_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'geometry_effect_fw_data/1_kinematic_cases')
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/5_SIM_RESULTS')
image_out_path = cwd
#-----------------------------------------
kinematics_data = os.path.join(kinematics_dir, kinematics_file)
CF_file_names = [f.name for f in os.scandir(data_dir) if f.is_file()]
cf_array = []
for cfi in cfd_data_list:
    cfi_name = 'ar' + '{0:.1f}'.format(ar) + cfi
    for CF_name in CF_file_names:
        if CF_name.startswith(cfi_name):
            cfd_datai = os.path.join(data_dir, CF_name)
            cf_arrayi = read_cfd_data(kinematics_data, cfd_datai)
            cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
           cycle_time, ar, 'against_t')
