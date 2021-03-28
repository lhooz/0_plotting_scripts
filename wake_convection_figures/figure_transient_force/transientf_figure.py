"""plotting cfd run results"""
import os

from tfplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
pa = 45
Re = [100, 1000, 10000]
stroke = [1.5, 3.0, 6.0]
at = [0.25]
pt = [0.25]
#-----------------------------------------
cfd_data_list = []
for re in Re:
    for s in stroke:
        for a in at:
            for p in pt:
                if pa == 45:
                    p = 0.125
                else:
                    p = 0.25
                cfd_data_name = 'Re' + '{0:.1f}'.format(
                    re) + '_stroke' + '{0:.1f}'.format(
                        s) + '_acf' + '{0:.3g}'.format(
                            a) + '_pf' + '{0:.3g}'.format(
                                p) + '_pa' + '{0:.1f}'.format(pa)
                cfd_data_list.append(cfd_data_name)
#-----------------------------------------
legendre = [r'$Re = ' + str(x) + '$' for x in Re]
legendsc = [r'$s/c = ' + str(x) + '$' for x in stroke]
legends = [legendre, legendsc]
#---------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
if pa == 0:
    time_to_plot = [0.0, 2.0]
else:
    time_to_plot = [1.0, 2.0]
# coeffs_show_range = [-7.0, 7.0]
coeffs_show_range = [-10.0, 4.0]
cycle_time = 1.0
#---------------------------------------
cwd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                  'wake_convection_data')
oimage_file = os.path.join(
    cwd, 'transient_force_pa' + '{0:.0f}'.format(pa) + '.png')
#---------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, '3_SIM_RESULTS', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, legends, time_to_plot, coeffs_show_range, oimage_file,
           cycle_time, 'against_t')
