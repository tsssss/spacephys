;+
; Test to load MMS data
;-

tr = time_double(['2016-08-04/22:15','2016-08-04/22:55'])
;tr = time_double(['2016-08-04/00:00','2016-08-04/02:00'])
data_tr = time_double(['2016-08-04','2016-08-06'])
;data_tr = time_double(['2016-08-09','2016-08-10'])
probe = '2'
prefix = 'mms'+probe+'_'

plot_vars = list()


; thermal.
foreach species, ['e','p','o'] do plot_vars.add, mms_read_en_spec(data_tr, probe=probe, species=species)


; keV electron.
foreach species, ['e','p','o'] do plot_vars.add, mms_read_en_spec_kev(data_tr, probe=probe, species=species)


; E field.
e_var = mms_read_efield(data_tr, probe=probe, coord='sm')
plot_vars.add, e_var
get_data, e_var, times, e_vec, limits=lim
spin_period = 20d
width = spin_period/sdatarate(times)
for ii=0,2 do e_vec[*,ii] = smooth(e_vec[*,ii],width,nan=1,edge_zero=1)
store_data, e_var+'_smooth', times, e_vec, limits=lim
plot_vars.add, lets_calc_vec_mag(e_var+'_smooth')

; Position.
r_var = mms_read_orbit(data_tr, probe=probe, coord='sm')

; B field.
b_var = mms_read_bfield(data_tr, probe=probe, coord='sm')
b_mag_var = lets_calc_vec_mag(b_var)
plot_vars.add, b_var
plot_vars.add, b_mag_var

tmp = ['_kev','']
plot_vars = prefix+['e_en_spec'+tmp,'p_en_spec'+tmp,'o_en_spec'+tmp,['e','b']+'_sm']

sgopen, 0, size=[6,8]
tplot, plot_vars, var_label=r_var, trange=data_tr
;tplot, [e_var,b_var,ele_var], var_label=r_var, trange=tr

end