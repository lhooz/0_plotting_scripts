"""plotting cfd run results"""
import os

from tfplotting_functions_pt import cf_plotter, read_cfd_data

#-------------input plot control----------
Re = [100, 1000, 10000]
stroke = [3, 5, 7, 9]
at = [0.25]
pt = [0.125, 0.25, 0.5]
#-----------------------------------------
cfd_data_list = []
for re in Re:
    for s in stroke:
        for a in at:
            for p in pt:
                cfd_data_name = 'Re' + '{0:.1f}'.format(
                    re) + '_stroke' + '{0:.1f}'.format(
                        s) + '_acf' + '{0:.3g}'.format(
                            a) + '_pf' + '{0:.3g}'.format(p)
                cfd_data_list.append(cfd_data_name)
#-----------------------------------------
markre = [r'$Re = ' + '{0:.0f}'.format(x) + '$' for x in Re]
marksc = [r'$s/c = ' + '{0:.1f}'.format(x) + '$' for x in stroke]
markpt = [r'$\^p_t = ' + '{0:.3g}'.format(x) + '$' for x in pt]
marks = [markre, marksc, markpt]
#---------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [0, 2.5]
coeffs_show_range = [-2.5, 10]
cycle_time = 1.1
#---------------------------------------
cd = os.getcwd()
wd = os.path.dirname(cd)
wd = os.path.dirname(wd)
image_out_path = cd
#---------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, '1_SIM_RESULTS', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, marks, time_to_plot, coeffs_show_range, image_out_path,
           cycle_time, 'against_t')
