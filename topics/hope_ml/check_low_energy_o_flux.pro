;+
; Check the low energy o+ flux data.
;-

time_range = time_double(['2013-02','2013-03'])
probe = 'a'
species = 'o'
energy_range = [10,1e3]

ae_var = ml_omni_read_ae(time_range)
dst_var = omni_read_symh(time_range)

var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, species=species, get_name=1)
if check_if_update(var, time_range) then begin
    var = ml_rbsp_hope_read_en_spec(time_range, probe=probe, species=species)
endif
options, var, 'color_table', 49
options, var, 'zrange', [1e3,1e6]

flux_var = var+'_test'
get_data, var, times, fluxs, energys
denergys = energys[1:-1]-energys[0:-2]
de_e = mean(denergys/((energys[1:-1]+energys[0:-2])*0.5))
diff_energy_fluxs = fluxs
foreach energy, energys, energy_id do diff_energy_fluxs[*,energy_id] *= energy
energy_index = lazy_where(energys, '[]', energy_range, count=nenergy)
the_fluxs = total(diff_energy_fluxs[*,energy_index],2)*de_e


dis_var = ml_rbsp_read_dis(time_range, probe=probe)
dis = get_var_data(dis_var, at=times)
index = where(dis le 3, count)
the_fluxs[index] = !values.f_nan

mlt_var = ml_rbsp_read_mlt(time_range, probe=probe)
mlts = get_var_data(mlt_var, at=times)
index = where(dis le 5.7)
mlts[index] = !values.f_nan
store_data, mlt_var, times, mlts

store_data, flux_var, times, the_fluxs
add_setting, flux_var, smart=1, dictionary($
    'display_type', 'scalar', $
    'short_name', 'O+ flux', $
    'unit', '#/cm!U2!N-s', $
    'ylog', 1, $
    'yrange', [1e5,1e10] )
    
end