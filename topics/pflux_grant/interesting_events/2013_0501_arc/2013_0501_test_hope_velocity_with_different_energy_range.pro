time_range = time_double(['2013-05-01/07:20','2013-05-01/07:50'])
probe = 'b'
prefix = 'rbsp'+probe+'_'

rgb = constant('rgb')

ion_energy_range = [50,1000]
rbsp_calc_hope_moments, time_range, probe=probe, $
    ion_energy_range=ion_energy_range
vars = prefix+['o','p']+'_vbulk_gse'
options, vars, 'colors', rgb
foreach var, vars do tmp = rename_var(var, output=var+'_low_energy')

ion_energy_range = [1000,5e5]
rbsp_calc_hope_moments, time_range, probe=probe, $
    ion_energy_range=ion_energy_range
vars = prefix+['o','p']+'_vbulk_gse'
options, vars, 'colors', rgb
foreach var, vars do tmp = rename_var(var, output=var+'_high_energy')

ion_energy_range = [50,5e5]
rbsp_calc_hope_moments, time_range, probe=probe, $
    ion_energy_range=ion_energy_range
vars = prefix+['o','p']+'_vbulk_gse'
options, vars, 'colors', rgb


vars = prefix+[$
    'o_vbulk_gse'+['','_'+['low','high']+'_energy'], $
    'p_vbulk_gse'+['','_'+['low','high']+'_energy']]
tplot, vars, trange=time_range
end