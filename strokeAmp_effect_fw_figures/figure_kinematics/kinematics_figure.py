"""plotting cfd run results"""
import os
import numpy as np

from kplotting_functions import k_plotter, read_kinematics_data

#-------------input plot control----------
kinematic_data_list = [
    'phi124_at0.5_pt0.5_3d.dat',
    'phi160_at0.25_pt0.25_3d.dat',
    'phi178_at0.125_pt0.125_3d.dat',
]
#---------------------------------------
time_to_plot = [0, 1.0]
show_range = [-1.2, 1.2]
#---------------------------------------
time_scale = 1.0
#---------------------------------------
cd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cd))),
                  'strokeAmp_effect_fw_data/1_kinematic_cases')
image_out_path = cd
#---------------------------------------
karray_all = []
for kinematic_data_name in kinematic_data_list:
    kinematics_data = os.path.join(wd, kinematic_data_name)
    k_arrayi, dk_arrayi, ddk_arrayi = read_kinematics_data(kinematics_data)
    ti = np.array([dk_arrayi[:, 0]])
    dphi = np.array([k_arrayi[:, 4]])
    daoa = np.array([k_arrayi[:, 5]])

    #--scaling factors--
    vel_scale = np.amax(np.abs(dphi)) / 2
    avel_scale = np.amax(np.abs(daoa))
    dphi += vel_scale
    #-------------------
    dphi = dphi / vel_scale
    daoa = daoa / avel_scale

    kine_array = np.append(ti, dphi, axis=0)
    kine_array = np.append(kine_array, daoa, axis=0)
    kine_array = np.transpose(kine_array)

    karray_all.append(kine_array)
#---------------------------------------
k_plotter(karray_all, time_scale, time_to_plot, show_range, image_out_path,
          'against_t')
