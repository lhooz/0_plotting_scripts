"""plotting wing profiles"""
import os
import numpy as np

from gplotting_functions import g_plotter_Joukowsky, read_profile_data

#-------------input plot control----------
planform_list = [
    'ar2.8_ofs0.0_r1h0.46',
    'ar3.0_ofs0.0_r1h0.55',
    'ar3.2_ofs0.0_r1h0.5',
    'ar3.3_ofs0.0_r1h0.49',
    'ar3.5_ofs0.0_r1h0.47',
    'ar3.6_ofs0.0_r1h0.48',
    'ar4.2_ofs0.0_r1h0.52',
    'ar5.3_ofs0.0_r1h0.56',
]
#---------------------------------------
x_range = [-0.05, 1.05]
y_range = [-0.37, 0.15]
# y_range = 'all'
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'Joukowsky_data/wing_profiles')
image_out_path = cwd
#---------------------------------------
marks = ['FF', 'BB', 'HM', 'HB', 'CF', 'HF', 'DF', 'LB']

data_array = []
for data_name in planform_list:
    profile_datai = os.path.join(data_dir, data_name + '.csv')
    profile_arrayi = read_profile_data(profile_datai)

    maxR = np.amax(np.abs(profile_arrayi[:, 0]))

    data_array.append(profile_arrayi / maxR)
#---------------------------------------
g_plotter_Joukowsky(data_array, marks, x_range, y_range, image_out_path)
