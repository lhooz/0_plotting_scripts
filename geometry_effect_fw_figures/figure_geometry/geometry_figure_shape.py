"""plotting wing profiles"""
import os
import numpy as np

from gplotting_functions import g_plotter, read_profile_data

#-------------input plot control----------
r1h = [0.4, 0.5, 0.6]
ar = [1.5, 3.0, 4.5, 6.0, 7.5]
ofs = 1.0
#-----------------------------------------
x_range = [0, 0.54]
y_range = [-0.12, 0.05]
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/0_wing_profiles')
image_out_path = cwd
#---------------------------------------
marks = [r'$\^r_1$ = ' + '{0:.1f}'.format(x) for x in r1h]
legends = [r'$AR$ = ' + '{0:.1f}'.format(x) for x in ar]

data_array = []
for r1hi in r1h:
    r1h_arrayi = []
    for ari in ar:
        datai = 'ar' + '{0:.1f}'.format(ari) + '_ofs' + '{0:.1f}'.format(
            ofs) + '_r1h' + '{0:.1f}'.format(r1hi)
        profile_datai = os.path.join(data_dir, datai + '.csv')
        profile_arrayi = read_profile_data(profile_datai)

        r1h_arrayi.append(profile_arrayi)

    data_array.append(r1h_arrayi)
#---------------------------------------
g_plotter(data_array, legends, marks, x_range, y_range, image_out_path)
