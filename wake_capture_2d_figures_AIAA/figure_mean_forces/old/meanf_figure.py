"""plotting cfd run results"""
import os

from mfplotting_functions import cf_plotter90, cf_plotter0, read_cfd_data

#-------------input plot control----------
Re = [100, 1000, 10000]
stroke = [3, 5, 7, 9]
at = [0.125, 0.25, 0.5]
pt90 = [0.125, 0.25, 0.5]
pt0 = [0.25]
#-----------------------------------------
cfd_data_list_pa0 = []
cfd_data_list_pa90 = []
for p in pt90:
    for a in at:
        for re in Re:
            for s in stroke:
                pafolder_pa90 = 'SIM_RESULTS_pa90'
                cfd_data_name = 'Re' + '{0:.1f}'.format(
                    re) + '_stroke' + '{0:.1f}'.format(
                        s) + '_acf' + '{0:.3g}'.format(
                            a) + '_pf' + '{0:.3g}'.format(p)
                cfd_data_path_pa90 = os.path.join(pafolder_pa90, cfd_data_name)
                cfd_data_list_pa90.append(cfd_data_path_pa90)

for p in pt0:
    for a in at:
        for re in Re:
            for s in stroke:
                pafolder_pa0 = 'SIM_RESULTS_pa0'
                cfd_data_name = 'Re' + '{0:.1f}'.format(
                    re) + '_stroke' + '{0:.1f}'.format(
                        s) + '_acf' + '{0:.3g}'.format(
                            a) + '_pf' + '{0:.3g}'.format(p)
                cfd_data_path_pa0 = os.path.join(pafolder_pa0, cfd_data_name)
                cfd_data_list_pa0.append(cfd_data_path_pa0)
#-----------------------------------------
x_data = stroke
markr90 = [r'$\^p_t = ' + '{0:.3g}'.format(x) + '$' for x in pt90]
markr0 = [r'$\^p_t = ' + '{0:.3g}'.format(x) + '$' for x in pt0]
markc = [r'$\^a_t = ' + '{0:.3g}'.format(x) + '$' for x in at]
marks90 = [markr90, markc]
marks0 = [markr0, markc]
legends = [r'$Re = ' + '{0:.0f}'.format(x) + '$' for x in Re]
#---------------------------------------
# x_range = 'all'
# y_range = 'all'
x_range = [2, 10]
cl_range = [-0.5, 0.6]
rl_range = [-0.4, 0.5]
cd_range = [-0.5, 0.6]
rd_range = [-0.4, 0.5]
rld_range = [-15, 15]
y_range = [cl_range, rl_range, cd_range, rd_range, rld_range]
y_label = [
    r'$\bar{C_{LW}}$', r'$\bar{C_{LW}}/\bar{C_{LS}}$', r'$\bar{C_{DW}}$',
    r'$\bar{C_{DW}}/\bar{C_{DS}}$', r'$\bar{C_{LW}}/\bar{C_{DW}}$'
]
#---------------------------------------
cd = os.getcwd()
wd = os.path.dirname(cd)
image_out_path = cd
#---------------------------------------
mcf_array_pa90 = []
for cfi in cfd_data_list_pa90:
    cfd_datai = os.path.join(wd, '1_SIM_RESULTS', cfi)
    mcf_arrayi = read_cfd_data(cfd_datai)
    mcf_array_pa90.append(mcf_arrayi)

data_array_pa90 = mcf_array_pa90
figname90 = 'pa90'
mcf_array_pa0 = []
for cfi in cfd_data_list_pa0:
    cfd_datai = os.path.join(wd, '1_SIM_RESULTS', cfi)
    mcf_arrayi = read_cfd_data(cfd_datai)
    mcf_array_pa0.append(mcf_arrayi)

data_array_pa0 = mcf_array_pa0
figname0 = 'pa0'
#---------------------------------------

cf_plotter90(x_data, data_array_pa90, marks90, legends, x_range, y_range,
             y_label, image_out_path, figname90)
cf_plotter0(x_data, data_array_pa0, marks0, legends, x_range, y_range,
             y_label, image_out_path, figname0)
