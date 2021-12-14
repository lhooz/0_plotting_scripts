"""plotting cfd run results"""
import os

from tfplotting_functions import cf_plotter_wake, read_cfd_data

#-------------input plot control----------
pa = 90
Re = [100, 1000]
# stroke = [6]
stroke = [1.5, 3.0, 4.5, 6.0]
at = [0.25]
pt = [0.25]
#-----------------------------------------
cfd_data_list = []
for s in stroke:
    for re in Re:
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
legendre = [r'$Re = 10^2$', r'$Re = 10^3$']
legendsc = [r'$s/c = ' + str(x) + '$' for x in stroke]
legends = [legendre, legendsc]
#---------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'

time_to_plot = [1.0, 2.0]
# coeffs_show_range = [-5.0, 7.0] #--pa0--
# coeffs_show_range = [-42.0, 34.0] #--pa45--
# coeffs_show_range = [-16.0, 37.0] #--pa90--
if pa == 0:
    coeffs_show_range = [-1.5, 1.5]  #--pa0 wake--
elif pa == 45:
    coeffs_show_range = [-9.0, 3.0]  #--pa45 wake--
elif pa == 90:
    coeffs_show_range = [-3.5, 3.5]  #--pa90 wake--
cycle_time = 1.0
#---------------------------------------
cwd = os.getcwd()
wd = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                  'wake_convection_data')
oimage_file = os.path.join(
    cwd, 'transient_force_pa' + '{0:.0f}'.format(pa) + '_wake.svg')
#---------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, '3_SIM_RESULTS', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter_wake(pa, data_array, legends, time_to_plot, coeffs_show_range,
                oimage_file, cycle_time, 'against_t')
