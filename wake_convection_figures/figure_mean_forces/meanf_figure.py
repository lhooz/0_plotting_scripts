"""plotting cfd run results"""
import os

from mfplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
pa = [0, 45, 90]
Re = [100, 1000, 10000]
stroke = [3, 5, 7, 9]
at = [0.25]
pt = [0.25]
#-----------------------------------------
cfd_data_list = []
for pai in pa:
    for re in Re:
        for s in stroke:
            for a in at:
                for p in pt:
                    if pai == 45:
                        p = 0.125
                    else:
                        p = 0.25
                    cfd_data_name = 'Re' + '{0:.1f}'.format(
                        re) + '_stroke' + '{0:.1f}'.format(
                            s) + '_acf' + '{0:.3g}'.format(
                                a) + '_pf' + '{0:.3g}'.format(
                                    p) + '_pa' + '{0:.1f}'.format(pai)
                    cfd_data_list.append(cfd_data_name)
#-----------------------------------------
x_data = stroke
markpa = [r'$\theta = ' + '{0:.0f}'.format(x) + '$' for x in pa]
legends = [r'$Re = ' + '{0:.0f}'.format(x) + '$' for x in Re]
#---------------------------------------
# x_range = 'all'
# y_range = 'all'
x_range = [2, 10]
cl_range = [-0.65, 0.5]
cd_range = cl_range
# cd_range = [-0.35, 0.55]
rl_range = [-0.4, 0.3]
rd_range = rl_range
# rd_range = [-0.2, 0.3]
rld_range = [-5, 2]
y_range = [cl_range, cd_range, rl_range, rd_range, rld_range]
y_label = [
    r'$\bar{C_{LW}}$', r'$\bar{C_{DW}}$', r'$\bar{C_{LW}}/\bar{C_{LS}}$',
    r'$\bar{C_{DW}}/\bar{C_{DS}}$', r'$\bar{C_{LW}}/\bar{C_{DW}}$'
]
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_capture_2d_journal_data/2_SIM_RESULTS')
image_out_path = cwd
#---------------------------------------
mcf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(data_dir, cfi)
    mcf_arrayi = read_cfd_data(cfd_datai)
    mcf_array.append(mcf_arrayi)
#---------------------------------------
cf_plotter(x_data, mcf_array, markpa, legends, x_range, y_range, y_label,
           image_out_path)
