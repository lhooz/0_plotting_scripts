"""plotting cfd run results"""
import os

from mfplotting_functions import cf_plotter, read_cfd_data, read_kinematics_data

#-------------input plot control----------
Re = [100]
AR = [1.5, 3.0, 4.5, 6.0, 7.5]
r1hat = [0.4, 0.5, 0.6]
offset = [0.0]
pt = [0.25]
#-----------------------------------------
kinematics_file = 'kinematics.dat'
cfd_data_list = []
for ar in AR:
    for r1h in r1hat:
        for ofs in offset:
            for re in Re:
                for p in pt:
                    cfd_data_name = 'ar' + '{0:.1f}'.format(
                        ar) + '_ofs' + '{0:.1f}'.format(
                            ofs) + '_r1h' + '{0:.1f}'.format(
                                r1h) + '__Re' + '{0:.1f}'.format(
                                    re) + '_pt' + '{0:.3g}'.format(p)
                    cfd_data_list.append(cfd_data_name)
#-----------------------------------------
x_data = r1hat
legends = ['AR = ' + '{0:.1f}'.format(ar) for ar in AR]
#---------------------------------------
# x_range = 'all'
# y_range = 'all'
x_range = [0.35, 0.65]
cl_range = [1.0, 1.8]
cd_range = [1.5, 3.0]
pf_range = [0.3, 1.0]
y_range = [cl_range, cd_range, pf_range]
y_label = [r'$\bar{C_L}$', r'$\bar{C_D}$', r'$\frac{1}{P^\ast}$']
#---------------------------------------
cwd = os.getcwd()
kinematics_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'geometry_effect_fw_data/1_kinematic_cases/3dbm_kinematic_cases')
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'geometry_effect_fw_data/5_SIM_RESULTS')
image_out_path = cwd
#---------------------------------------
CF_file_names = [f.name for f in os.scandir(data_dir) if f.is_file()]
mcf_array = []
for cfi in cfd_data_list:
    for CF_name in CF_file_names:
        if CF_name.startswith(cfi):
            kinematics_datai = os.path.join(kinematics_dir, CF_name)
            cfd_datai = os.path.join(data_dir, CF_name)
            u2, karr, dkarr, ddkarr = read_kinematics_data(kinematics_datai)
            mcf_arrayi = read_cfd_data(cfd_datai, u2, karr, dkarr)
            mcf_array.append(mcf_arrayi)
#---------------------------------------
cf_plotter(x_data, mcf_array, legends, x_range, y_range, y_label,
           image_out_path)
