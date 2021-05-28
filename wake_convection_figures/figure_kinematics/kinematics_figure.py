"""plotting cfd run results"""
import os
import numpy as np

from kplotting_functions import kf_plotter, illustrationk_plotter, read_kinematics_data

#-------------input plot control----------
illustration_data = 'Re100.0_stroke5.0_acf0.25_pf0.25_pa90.0'
kinematic_data_list = [
    'Re100.0_stroke5.0_acf0.25_pf0.25_pa0.0',
    'Re100.0_stroke5.0_acf0.25_pf0.125_pa45.0',
    'Re100.0_stroke5.0_acf0.25_pf0.25_pa90.0',
]
#-----------------------------------------
time_scale = 1.0
kinematics_t = np.linspace(0.0, 1.0, 30)
#-----------------------------------------
# time_to_plot = 'all'
# show_range = 'all'
time_to_plot = [0, 2.0]
show_range = [-0.1, 1.4]
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_capture_2d_journal_data/0_kinematic_cases')
image_out_path = cwd
#---------------------------------------
marks = [
    r'$\alpha_{int}$ = 135$^\circ$', r'$\alpha_{int}$ = 90$^\circ$',
    r'$\alpha_{int}$ = 45$^\circ$'
]

data_array = []
idata_array = []
for ki in kinematic_data_list:
    kinematics_datai = os.path.join(data_dir, ki + '.dat')
    k_arrayi, dk_arrayi, ddk_arrayi = read_kinematics_data(kinematics_datai)

    ti = np.array([dk_arrayi[:, 0]])
    dtransi = np.abs(np.array([dk_arrayi[:, 1]]))
    daoai = np.abs(np.array([dk_arrayi[:, 6]]))

    transi = np.array([k_arrayi[:, 1]])
    aoai = np.array([k_arrayi[:, 6]])

    #--scaling factors--
    vel_scalei = np.amax(dtransi)
    avel_scalei = np.amax(daoai)
    #-------------------
    dtransi = dtransi / vel_scalei
    if avel_scalei >= 1e-3:
        daoai = daoai / avel_scalei

    arrayi = np.append(ti, dtransi, axis=0)
    arrayi = np.append(arrayi, daoai, axis=0)
    arrayi = np.transpose(arrayi)

    iarrayi = np.append(ti, transi, axis=0)
    iarrayi = np.append(iarrayi, aoai, axis=0)
    iarrayi = np.transpose(iarrayi)

    data_array.append(arrayi)
    idata_array.append(iarrayi)
#---------------------------------------
illustration_file = os.path.join(data_dir, illustration_data + '.dat')
k_array, dk_array, ddk_array = read_kinematics_data(illustration_file)
it = np.array([k_array[:, 0]])
itrans = np.array([k_array[:, 1]])
iaoa = np.array([k_array[:, 6]])

iarray = np.append(ti, itrans, axis=0)
iarray = np.append(iarray, iaoa, axis=0)
iarray = np.transpose(iarray)

ifigure_parameters = [
    float(illustration_data.split('_')[1].split('stroke')[1]),
    float(illustration_data.split('_')[2].split('acf')[1]),
    float(illustration_data.split('_')[3].split('pf')[1])
]
#---------------------------------------
kf_plotter(kinematic_data_list, kinematics_t, data_array, idata_array, iarray,
           ifigure_parameters, marks, time_to_plot, show_range, image_out_path,
           time_scale)
