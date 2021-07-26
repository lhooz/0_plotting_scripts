"""plotting cfd run results"""
import os
import numpy as np

from kplotting_functions import k_plotter, read_kinematics_data

#-------------input plot control----------
# kinematic_data_list = [
    # 'smooth_phi160_atf0_ptc0_3d.dat',
    # 'smooth_phi160_atf0_ptc1_3d.dat',
    # 'smooth_phi160_atf0_ptc5_3d.dat',
    # 'smooth_phi160_atf0_ptc10_3d.dat',
    # 'smooth_phi160_atf0_ptc100_3d.dat',
# ]
kinematic_data_list = [
    'smooth_phi160_atf0_ptc5_3d.dat',
]
#---------------------------------------
time_to_plot = [0, 1.0]
show_range = [-0.2, 1.2]
#---------------------------------------
time_scale = 1.0
#---------------------------------------
cd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cd))),
                  'Joukowsky_data/kinematics')
image_out_path = cd
#---------------------------------------
karray_all = []
for kinematic_data_name in kinematic_data_list:
    kinematics_data = os.path.join(wd, kinematic_data_name)
    k_arrayi, dk_arrayi, ddk_arrayi = read_kinematics_data(kinematics_data)
    ti = np.array([dk_arrayi[:, 0]])
    dphi = np.array([dk_arrayi[:, 4]]) * -1
    daoa = np.array([k_arrayi[:, 5]])

    #--scaling factors--
    vel_scale = np.amax(np.abs(dphi))
    avel_scale = np.amax(np.abs(daoa))
    # dphi += vel_scale
    #-------------------
    dphi = np.abs(dphi) / vel_scale
    daoa = 90 - np.abs(daoa)

    kine_array = np.append(ti, dphi, axis=0)
    kine_array = np.append(kine_array, daoa, axis=0)
    kine_array = np.transpose(kine_array)

    karray_all.append(kine_array)
#---------------------------------------
k_plotter(karray_all, time_scale, time_to_plot, show_range, image_out_path,
          'against_t')
