"""plotting cfd run results"""
import os
import numpy as np

from vplotting_functions_flapping import cf_plotter, read_cfd_data, read_ref_data, read_kinematics_data

#-------------input plot control----------
kinematic_data_list = ['sinusoidal_2d']
cfd_data_list = ['CFD']
ref_data_lst = ['CL_LEE', 'CL_Wang', 'CD_LEE', 'CD_Wang']
#---------------------------------------
# time_to_plot = 'all'
# show_range = 'all'
time_to_plot = [3.0, 4.0]
show_range = [-0.05, 1.25]
#---------------------------------------
time_scale = 1.0
#---------------------------------------
validation_data_folder = '2_validation_results/flapping_validation'
cd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cd))),
                  'wake_convection_data')
image_out_path = cd
#---------------------------------------
# legendcfd = [
# r'simulation, $\alpha = 45 ^\circ$', r'simulation, $\alpha = 60 ^\circ$'
# ]
# legendexp = [
# r'experiment, $\alpha = 45 ^\circ$', r'experiment, $\alpha = 60 ^\circ$'
# ]
legendcfd = [r'Current CFD', r'Current CFD']
legendexp = [r'Lee et al.', r'Wang', r'Lee et al.', r'Wang']
legends = [legendcfd, legendexp]
#---------------------------------------
kine_array = []
for ki in kinematic_data_list:
    kinematics_datai = os.path.join(wd, validation_data_folder, ki + '.dat')
    k_arrayi, dk_arrayi, ddk_arrayi = read_kinematics_data(kinematics_datai)
    ti = np.array([dk_arrayi[:, 0]])
    dtransi = np.abs(np.array([dk_arrayi[:, 1]]))
    aoai = np.abs(np.array([k_arrayi[:, 6]]))

    #--scaling factors--
    vel_scalei = np.amax(dtransi)
    #-------------------
    dtransi = dtransi / vel_scalei

    arrayi = np.append(ti, dtransi, axis=0)
    arrayi = np.append(arrayi, aoai, axis=0)
    arrayi = np.transpose(arrayi)

    kine_array.append(arrayi)
#---------------------------------------

cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, validation_data_folder, cfi + '.dat')
    cf_arrayi = read_cfd_data(cfd_datai, time_to_plot)

    cf_array.append(cf_arrayi)

ref_array = []
for refi in ref_data_lst:
    ref_datai = os.path.join(wd, validation_data_folder, refi + '.csv')
    ref_arrayi = read_ref_data(ref_datai)

    ref_array.append(ref_arrayi)

data_array = [kine_array, cf_array, ref_array]
#---------------------------------------
cf_plotter(data_array, time_scale, legends, time_to_plot, show_range,
           image_out_path, 'against_t')
