"""plotting cfd run results"""
import os
import numpy as np

from kplotting_functions_illustrate import kf_plotter, illustrationk_plotter, read_kinematics_data

#-------------input plot control----------
illustration_data = 'kinematics'
#-----------------------------------------
stroke = 4.5
kinematics_t = np.linspace(0.5, 1, 25)
#-----------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/1_kinematic_cases')
image_out_path = cwd
#---------------------------------------
illustration_file = os.path.join(data_dir, illustration_data + '.dat')
k_array, dk_array, ddk_array = read_kinematics_data(illustration_file)
ti = np.array([dk_array[:, 0]])
it = np.array([k_array[:, 0]])
itrans = np.array([k_array[:, 4]]) / 160 * stroke
iaoa = np.array([k_array[:, 5]])

iarray = np.append(ti, itrans, axis=0)
iarray = np.append(iarray, iaoa, axis=0)
iarray = np.transpose(iarray)
#---------------------------------------
kf_plotter(iarray, kinematics_t, image_out_path, stroke)
