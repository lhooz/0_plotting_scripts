"""plotting cfd run results"""
import os
import numpy as np

from vplotting_functions import cf_plotter, read_cfd_data, read_ref_data

#-------------input plot control----------
# ad_data_list = ['ad_0', 'ad_exp', 'ad_cfd']
# sym_data_list = ['sym_0', 'sym_exp', 'sym_cfd1', 'sym_cfd2']
# dl_data_list = ['dl_0', 'dl_exp', 'dl_cfd']
ad_data_list = ['FF_advanced_0', 'ad_exp', 'ad_cfd', 'SunAndTang_advanced']
sym_data_list = ['FF_symmetric_0', 'sym_exp', 'sym_cfd1', 'sym_cfd2']
dl_data_list = ['FF_delayed_0', 'dl_exp', 'dl_cfd']
#---------------------------------------
time_to_plot = [4.0, 5.0]
# show_range = [-0.8, 1.2]
show_range = 'all'
#---------------------------------------
time_scale = 1.0
#---------------------------------------
cd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cd))),
                  'geometry_effect_fw_data/3_validation_data')
image_out_path = cd
#---------------------------------------
legendad = [
    'Current CFD', 'Experiment, Dickinson et al.', 'CFD, Erzincanli and Sahin', 'CFD, Sun and Tang'
]
legendsym = [
    'Current CFD', 'Experiment, Dickinson et al.', 'CFD, Erzincanli and Sahin',
    'CFD, Kweon and Choi'
]
legenddl = [
    'Current CFD', 'Experiment, Dickinson et al.', 'CFD, Erzincanli and Sahin'
]
legends = [legendad, legendsym, legenddl]
#---------------------------------------
ad_array = []
for cfi in ad_data_list:
    datai = os.path.join(wd, cfi + '.dat')
    if cfi.endswith('0'):
        arrayi = read_cfd_data(datai)
    else:
        arrayi = read_ref_data(datai)

    ad_array.append(arrayi)

sym_array = []
for cfi in sym_data_list:
    datai = os.path.join(wd, cfi + '.dat')
    if cfi.endswith('0'):
        arrayi = read_cfd_data(datai)
    else:
        arrayi = read_ref_data(datai)

    sym_array.append(arrayi)

dl_array = []
for cfi in dl_data_list:
    datai = os.path.join(wd, cfi + '.dat')
    if cfi.endswith('0'):
        arrayi = read_cfd_data(datai)
    else:
        arrayi = read_ref_data(datai)

    dl_array.append(arrayi)

data_array = [ad_array, sym_array, dl_array]
#---------------------------------------
cf_plotter(data_array, time_scale, legends, time_to_plot, show_range,
           image_out_path)
