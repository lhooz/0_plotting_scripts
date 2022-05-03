"""plotting cfd run results"""
import os
import numpy as np

from kplotting_functions_wave import k_plotter, read_kinematics_data

#-------------input plot control----------
kinematic_data_name = 'ar1.5_ofs0.0_r1h0.4__Re100.0_pt0.25_Ro0.71'
#---------------------------------------
time_to_plot = [0, 1.0]
show_range = [[-110, 110], [35, 100]]
#---------------------------------------
time_scale = 1.0
#---------------------------------------
cd = os.getcwd()
wd = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cd))),
    'geometry_effect_fw_data/1_kinematic_cases/3dbm_kinematic_cases')
image_out_path = cd
#---------------------------------------
kinematics_data = os.path.join(wd, kinematic_data_name + '.dat')
k_arrayi, dk_arrayi, ddk_arrayi = read_kinematics_data(kinematics_data)
ti = np.array([dk_arrayi[:, 0]])
phi = np.array([k_arrayi[:, 4]])
aoa = 90 - np.array([k_arrayi[:, 5]])

#--scaling factors--
phi_scale = np.amax(np.abs(phi)) / 2
phi += phi_scale
#-------------------
kine_array = np.append(ti, phi, axis=0)
kine_array = np.append(kine_array, aoa, axis=0)
kine_array = np.transpose(kine_array)

karray_all = [kine_array]
#---------------------------------------
k_plotter(karray_all, time_scale, time_to_plot, show_range, image_out_path,
          'against_t')
