;+
; Check hope energy bins over several months.
;-

time_range = ['2015-03-15','2015-03-18']
probe = 'a'
prefix = 'rbsp'+probe+'_'

time_range = time_double(time_range)

rbsp_read_en_spec, time_range, probe=probe

files = rbsp_load_hope(time_range, probe=probe, id='l2%sector')
var_list = list()
; flags_ion are just all 0's.
energy_bin_var = 'HOPE_ENERGY_Ion'
var = prefix+'energy_ion'
var_list.add, dictionary($
    'in_vars', [energy_bin_var], $
    'out_vars', [var], $
    'time_var_name', 'Epoch_Ion', $
    'time_var_type', 'epoch' )
energy_bin_var = 'HOPE_ENERGY_Ele'
var = prefix+'energy_ele'
var_list.add, dictionary($
    'in_vars', [energy_bin_var], $
    'out_vars', [var], $
    'time_var_name', 'Epoch_Ele', $
    'time_var_type', 'epoch' )
read_vars, time_range, files=files, var_list=var_list
get_data, var, times, data
nbin = n_elements(data[0,*])
store_data, var, times, data, findgen(nbin), limits={$
    spec:1, zrange: [1,5e4], zlog:1, ztitle:'(eV)', ystyle:1, no_interp:1 }
tplot, var, trange=time_range
end