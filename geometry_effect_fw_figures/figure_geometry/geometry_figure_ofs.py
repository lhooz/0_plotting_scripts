"""plotting wing profiles"""
import os
import numpy as np

from gplotting_functions import g_plotter_ofs, read_profile_data

#-------------input plot control----------
r1h = [0.5]
ar = [3.0]
offset = [0.0, 1.0, 2.0, 3.0]
#-----------------------------------------
x_range = [0, 0.54]
y_range = [-0.12, 0.05]
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/0_wing_profiles')
image_out_path = cwd
#---------------------------------------
legends = [r'$\^r_R$ = ' + '{0:.1f}'.format(x) for x in offset]

data_array = []
for r1hi in r1h:
    for ari in ar:
        for ofs in offset:
            datai = 'ar' + '{0:.1f}'.format(ari) + '_ofs' + '{0:.1f}'.format(
                ofs) + '_r1h' + '{0:.1f}'.format(r1hi)
            profile_datai = os.path.join(data_dir, datai + '.csv')
            profile_arrayi = read_profile_data(profile_datai)

            data_array.append(profile_arrayi)
#---------------------------------------
g_plotter_ofs(data_array, legends, x_range, y_range, image_out_path)
