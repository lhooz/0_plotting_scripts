"""plotting cfd run results"""
import os

from mfplotting_functions import cf_plotter, read_cfd_data, read_kinematics_data

#-------------input plot control----------
Re = [100, 1000, 10000]
AR = [2, 3, 4, 5]
r1hat = [0.4, 0.5, 0.6]
offset = [0.0, 0.1, 0.2]
pt = [0.25]
#-----------------------------------------
cfd_data_list = []
for r1h in r1hat:
    for ofs in offset:
        for re in Re:
            for ar in AR:
                for p in pt:
                    cfd_data_name = 'ar' + '{0:.1f}'.format(
                        ar) + '_ofs' + '{0:.1f}'.format(
                            ofs) + '_r1h' + '{0:.1f}'.format(
                                r1h) + '__Re' + '{0:.1f}'.format(
                                    re) + '_pt' + '{0:.3g}'.format(p)
                    cfd_data_list.append(cfd_data_name)
#-----------------------------------------
x_data = AR
markr = ['Area centroid = ' + '{0:.3g}'.format(x) for x in r1hat]
markc = ['Offset = ' + '{0:.3g}'.format(x) for x in offset]
legends = ['Re = ' + '{0:.1f}'.format(x) for x in Re]
marks = [markr, markc]
#---------------------------------------
# x_range = 'all'
# y_range = 'all'
x_range = [1, 6]
cl_range = [1.2, 3.3]
cp_range = [0.5, 1.2]
y_range = [cl_range, cp_range]
y_label = [r'$\bar{C_L}$', r'$\frac{1}{P^\ast}$']
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data')
image_out_path = cwd
#---------------------------------------
mcf_array = []
for cfi in cfd_data_list:
    kinematic_datai = os.path.join(data_dir, '0_kinematic_cases', cfi)
    cfd_datai = os.path.join(data_dir, '3_SIM_RESULTS', cfi)
    u2, karr, dkarr, ddkarr = read_kinematics_data(kinematic_datai)
    mcf_arrayi = read_cfd_data(cfd_datai, u2, karr, dkarr)
    mcf_array.append(mcf_arrayi)
#---------------------------------------
cf_plotter(x_data, mcf_array, marks, legends, x_range, y_range, y_label,
           image_out_path)
