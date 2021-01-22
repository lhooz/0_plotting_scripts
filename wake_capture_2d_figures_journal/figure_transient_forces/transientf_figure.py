"""plotting cfd run results"""
import os

from tfplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
pa = [0, 45, 90]
Re = [100, 1000, 10000]
stroke = [3, 5, 7, 9]
at = [0.25]
pt = [0.25]
#-----------------------------------------
cfd_data_list = []
for re in Re:
    for pai in pa:
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
markre = [r'$Re = ' + '{0:.0f}'.format(x) + '$' for x in Re]
marksc = [r'$s/c = ' + '{0:.1f}'.format(x) + '$' for x in stroke]
markpa = [r'$\theta = ' + '{0:.0f}'.format(x) + '$' for x in pa]
marks = [markre, markpa, marksc]
#---------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [0, 2.0]
show_range_cl = [-2.2, 7.5]
show_range_cd = [-13, 13]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
cd = os.getcwd()
wd = os.path.dirname(cd)
wd = os.path.dirname(wd)
image_out_path = cd
#---------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, '2_SIM_RESULTS', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, marks, time_to_plot, show_range, image_out_path,
           cycle_time, 'against_t')
