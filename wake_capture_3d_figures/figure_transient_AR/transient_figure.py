"""script for mesh and convergence figure plots"""

import os
import numpy as np
from tplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
AR = 3.0
Re = 1000.0
#-----------------------------------------
case_name = 'phi160.0__ar' + '{0:.1f}'.format(
    AR) + '_ofs0.0_r1h0.5__Re' + '{0:.1f}'.format(Re) + '_pt0.25'
#-----------------------------------------
legends = [
    r'1_stroke',
    r'5_stroke',
    r'wake_effect',
]
#-----------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [0.0, 0.5]
show_range_cl = [-1.5, 3.5]
show_range_cd = [-1.5, 3.5]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
#-----------------------------------------
# ---case dir setup---------------
cwd = os.getcwd()
data_dir0 = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'threeD_wakeCapture/SIM_RESULTS_')

kinematics_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
    'threeD_wakeCapture/3dbm_kinematic_cases_wakeCapture')

image_out_path = cwd
#-----------------------------------------
stroke = 1
data_dir = data_dir0 + '{0:.0f}'.format(stroke) + 'Stroke'
for f in os.scandir(data_dir):
    if f.name.startswith(case_name):
        kinematics_data = os.path.join(kinematics_dir, f.name + '.dat')
        cfd_datai = os.path.join(data_dir, f.name)
        cf_arrayi = read_cfd_data(kinematics_data, cfd_datai)
        cf_array1S = cf_arrayi

#-----------------------------------------
stroke = 5
data_dir = data_dir0 + '{0:.0f}'.format(stroke) + 'Stroke'
for f in os.scandir(data_dir):
    if f.name.startswith(case_name):
        kinematics_data = os.path.join(kinematics_dir, f.name + '.dat')
        cfd_datai = os.path.join(data_dir, f.name)
        cf_arrayi = read_cfd_data(kinematics_data, cfd_datai)
        cf_array5S = cf_arrayi

data_array = [cf_array1S, cf_array5S]
#---------------------------------------

cf_plotter(AR, Re, data_array, legends, time_to_plot, show_range, image_out_path,
           cycle_time, 'against_t')
