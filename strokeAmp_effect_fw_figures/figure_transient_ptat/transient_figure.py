"""script for mesh and convergence figure plots"""

import os

from tplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
AR = 6
at_folder = 50
pt_all = [0.125, 0.25, 0.5]
#-----------------------------------------
cfd_data_list = []
legends = []
for pt in pt_all:
    cfd_data_list.append('phi160.0__ar' + '{:.1f}'.format(AR) + '_ofs0.0_r1h0.5__Re100.0_pt' + '{:.3g}'.format(pt))
    legends.append(r'$\^t_p$ = ' + '{:.3g}'.format(pt))
#-----------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [4.0, 5.0]
show_range_cl = [-0.3, 3.7]
show_range_cd = [-0.8, 4]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
#-----------------------------------------
cwd = os.getcwd()
kinematics_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'strokeAmpEffect_data/1_kinematic_cases/at' + '{:.0f}'.format(at_folder))
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'strokeAmpEffect_data/5_SIM_RESULTS/at' + '{:.0f}'.format(at_folder))
image_out_path = cwd
#-----------------------------------------
CF_file_names = [f.name for f in os.scandir(data_dir) if f.is_file()]
cf_array = []
for cfi in cfd_data_list:
    for CF_name in CF_file_names:
        if CF_name.startswith(cfi):
            kinematics_data = os.path.join(kinematics_dir, CF_name + '.dat')
            cfd_datai = os.path.join(data_dir, CF_name)
            cf_arrayi = read_cfd_data(kinematics_data, cfd_datai)
            cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, show_range, image_out_path,
           cycle_time, AR, pt_all, at_folder, 'against_t')
