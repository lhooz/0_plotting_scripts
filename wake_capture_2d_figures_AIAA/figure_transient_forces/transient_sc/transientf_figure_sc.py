"""plotting cfd run results"""
import os

from tfplotting_functions_sc import cf_plotter, read_cfd_data

#-------------input plot control----------
Re = [100, 1000, 10000]
stroke = [3, 5, 7, 9]
at = [0.25]
pt = [0.25]
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
legendre = [r'$Re = 10^2$', r'$Re = 10^3$', r'$Re = 10^4$']
legendsc = [r'$s/c = ' + str(x) + '$' for x in stroke]
legends = [legendre, legendsc]
#---------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [0, 2.5]
coeffs_show_range = [-8, 13]
cycle_time = 1.0
#---------------------------------------
cd = os.getcwd()
wd = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(cd)))),
    'wake_capture_2d_AIAA_data')
image_out_path = cd
#---------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, '2_SIM_RESULTS', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, coeffs_show_range,
           image_out_path, cycle_time, 'against_t')
