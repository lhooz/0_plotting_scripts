"""plotting cfd run results"""
import os

from mfplotting_functions import cf_plotter, read_cfd_data

#-------------input plot control----------
pa = [0, 45, 90]
Re = [100, 1000]
stroke = [1.5, 3.0, 6.0]
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
#----effects for each case-------
markEr1 = [['i', 's', 'i'], ['i', 's', 's']]
markEr2 = [['i', 'i', 'i'], ['i', 'i', 's']]
markEr3 = [['s', 'i', 's'], ['s', 'i', 's']]
markEffects = [markEr1, markEr2, markEr3]
#--------------------------------
x_data = stroke
markc = [
    r'$\alpha_E$ = 135$^\circ$', r'$\alpha_E$ = 90$^\circ$',
    r'$\alpha_E$ = 45$^\circ$'
]
legends = [r'$Re = 10^2$', r'$Re = 10^3$']
#---------------------------------------
# x_range = 'all'
# y_range = 'all'
x_range = [1, 6.5]
cl_range = [-4.5, 4.5]
cd_range = cl_range
rl_range = [-2, 2]
rd_range = rl_range
y_range = [cl_range, cd_range, rl_range, rd_range]
y_label = [
    r'$\bar{C}_{lw}$', r'$\bar{C}_{dw}$', r'$\bar{C}_{lw}/\bar{C}_{ls}$',
    r'$\bar{C}_{dw}/\bar{C}_{ds}$'
]
#---------------------------------------
cwd = os.getcwd()
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(cwd))),
                        'wake_convection_data/3_SIM_RESULTS')
image_out_path = cwd
#---------------------------------------
mcf_array = []
for cfi in cfd_data_list:
    cfd_datai = os.path.join(data_dir, cfi)
    mcf_arrayi = read_cfd_data(cfd_datai)
    mcf_array.append(mcf_arrayi)
#---------------------------------------
with open('meancf_all.dat', 'w') as f:
    for data_name, mvalues in zip(cfd_data_list, mcf_array):
        f.write("%s:\n" % data_name)
        f.write(
            "mcl_s = %s, mcl_w = %s, ratio_l = %s, mcd_s = %s, mcd_w = %s, ratio_d = %s, ratio_ld = %s\n"
            % (
                '{0:.8g}'.format(mvalues[0]),
                '{0:.8g}'.format(mvalues[1]),
                '{0:.8g}'.format(mvalues[2]),
                '{0:.8g}'.format(mvalues[3]),
                '{0:.8g}'.format(mvalues[4]),
                '{0:.8g}'.format(mvalues[5]),
                '{0:.8g}'.format(mvalues[6]),
            ))

cf_plotter(x_data, mcf_array, markc, markEffects, legends, x_range, y_range,
           y_label, image_out_path)
