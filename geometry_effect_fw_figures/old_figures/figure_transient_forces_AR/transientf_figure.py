"""plotting cfd run results"""
import os

from tfplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
Re = [100, 1000, 10000]
AR = [2, 3, 4, 5]
r1hat = [0.4, 0.5, 0.6]
offset = [0.0]
pt = [0.25]
#-----------------------------------------
cfd_data_list = []
for re in Re:
    for r1h in r1hat:
        for ar in AR:
            for ofs in offset:
                for p in pt:
                    cfd_data_name = 'ar' + '{0:.1f}'.format(
                        ar) + '_ofs' + '{0:.1f}'.format(
                            ofs) + '_r1h' + '{0:.1f}'.format(
                                r1h) + '__Re' + '{0:.1f}'.format(
                                    re) + '_pt' + '{0:.3g}'.format(p)
                    cfd_data_list.append(cfd_data_name)
#-----------------------------------------
markre = ['Re = ' + '{0:.1f}'.format(x) for x in Re]
markr1h = ['Area Centroid = ' + '{0:.1f}'.format(x) for x in r1hat]
markar = ['AR = ' + '{0:.1f}'.format(x) for x in AR]
marks = [markre, markr1h, markar]
#---------------------------------------
time_to_plot = 'all'
coeffs_show_range = 'all'
time_to_plot = [4.0, 5.0]
show_range_cl = [-2.0, 5.0]
show_range_cd = [-8.0, 8.0]
cycle_time = 1.0
#---------------------------------------
show_range = [show_range_cl, show_range_cd]
cd = os.getcwd()
wd = os.path.dirname(cd)
wd = os.path.dirname(wd)
wd = os.path.join(os.path.dirname(wd), 'geometry_effect_fw_data')
image_out_path = cd
#---------------------------------------
cf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(wd, '3_SIM_RESULTS', cfi)
    cf_arrayi = read_cfd_data(cfd_datai)
    cf_array.append(cf_arrayi)

data_array = cf_array
#---------------------------------------

cf_plotter(data_array, marks, time_to_plot, show_range, image_out_path,
           cycle_time, 'against_t')
