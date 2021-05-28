"""plotting cfd run results"""
import os

from mfplotting_functions import cf_plotter90, read_cfd_data

#-------------input plot control----------
Re = [100, 1000]
stroke = [3, 5, 7, 9]
at = [0.125, 0.25, 0.5]
pt90 = [0.125, 0.25, 0.5]
pt0 = [0.25]
#-----------------------------------------
cfd_data_list_pa90 = []
for a in at:
    for p in pt90:
        for re in Re:
            for s in stroke:
                pafolder_pa90 = 'SIM_RESULTS_pa90'
                cfd_data_name = 'Re' + '{0:.1f}'.format(
                    re) + '_stroke' + '{0:.1f}'.format(
                        s) + '_acf' + '{0:.3g}'.format(
                            a) + '_pf' + '{0:.3g}'.format(p)
                cfd_data_path_pa90 = os.path.join(pafolder_pa90, cfd_data_name)
                cfd_data_list_pa90.append(cfd_data_path_pa90)
#-----------------------------------------
x_data = stroke
markp = [r'$\^p_t = ' + '{0:.3g}'.format(x) + '$' for x in pt90]
marka = [r'$\^a_t = ' + '{0:.3g}'.format(x) + '$' for x in at]
marks90 = [marka, markp]
legends = [r'$Re = 10^2$', r'$Re = 10^3$']
#---------------------------------------
# x_range = 'all'
# y_range = 'all'
x_range = [2, 10]
# cl_range = [-0.75, 0.75]
# rl_range = [-0.5, 0.5]
# cd_range = [-0.75, 0.75]
# rd_range = [-0.5, 0.5]
cl_range = [-3.5, 3.5]
rl_range = [-1, 1]
cd_range = [-3.5, 3.5]
rd_range = [-1, 1]
rld_range = [-15, 15]
y_range = [cl_range, rl_range, cd_range, rd_range, rld_range]
y_label = [
    r'$\bar{C_{lw}}$', r'$\bar{C_{lw}}/\bar{C_{ls}}$', r'$\bar{C_{dw}}$',
    r'$\bar{C_{dw}}/\bar{C_{ds}}$', r'$\bar{C_{lw}}/\bar{C_{dw}}$'
]
#---------------------------------------
cd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cd))),
                  'wake_capture_2d_AIAA_data')
image_out_path = cd
#---------------------------------------
mcf_array_pa90 = []
for cfi in cfd_data_list_pa90:
    cfd_datai = os.path.join(wd, '1_SIM_RESULTS', cfi)
    mcf_arrayi = read_cfd_data(cfd_datai)
    mcf_array_pa90.append(mcf_arrayi)

data_array_pa90 = mcf_array_pa90
figname90 = 'pa90'
#---------------------------------------

cf_plotter90(x_data, data_array_pa90, marks90, legends, x_range, y_range,
             y_label, image_out_path, figname90)
# cf_plotter0(x_data, data_array_pa0, marks0, legends, x_range, y_range,
# y_label, image_out_path, figname0)
