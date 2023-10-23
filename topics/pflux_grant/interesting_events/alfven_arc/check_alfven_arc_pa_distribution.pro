
units = ['energy','velocity']
test_species = ['p','o','e','he']

probe = 'a'
prefix = 'rbsp'+probe+'_'
time_range = time_double(['2015-04-16/07:40','2015-04-16/08:30'])
plot_dir = join_path([homedir(),'rbsp_hope','rbsp'+probe])

foreach species, test_species do begin
    case species of
        'o': zrange = [5e3,5e5]
        'p': zrange = [1e4,1e6]
        'e': zrange = [5e5,1e10]
        'he': zrange = [1e3,1e5]
    endcase

    foreach unit, units do begin
        suffix = ''
        the_dir = join_path([plot_dir,species+'_'+unit+suffix])
        var = rbsp_plot_pa2d(time_range, probe=probe, species=species, unit=unit, $
            zrange=zrange, contour=0, plot_dir=the_dir, combine_energy_bin=combine_energy_bin)
    endforeach
endforeach
stop




units = ['energy','velocity']
test_species = ['p','o','e','he']

event_info = alfven_arc_load_data('2015_0302_1100')
rbsp_info = event_info.rbsp.rbspa
prefix = rbsp_info['prefix']
probe = rbsp_info['probe']
time_range = event_info['time_range']
plot_dir = event_info.plot_dir
combine_energy_bin = 0
suffix = (combine_energy_bin)? '_energy_combined': ''
foreach species, test_species do begin
    case species of
        'o': zrange = [5e3,5e5]
        'p': zrange = [1e4,1e6]
        'e': zrange = [5e5,1e10]
        'he': zrange = [1e3,1e5]
    endcase

    foreach unit, units do begin
        the_dir = join_path([plot_dir,species+'_'+unit+suffix])
        var = rbsp_plot_pa2d(time_range, probe=probe, species=species, unit=unit, $
            zrange=zrange, contour=0, plot_dir=the_dir, combine_energy_bin=combine_energy_bin)
    endforeach
endforeach

stop

units = ['energy','velocity']
test_species = ['p','o','e','he']

event_info = alfven_arc_load_data('2015_0105_0100')
rbsp_info = event_info.rbsp.rbspa
prefix = rbsp_info['prefix']
probe = rbsp_info['probe']
time_range = event_info['time_range']
plot_dir = event_info.plot_dir
combine_energy_bin = 0
suffix = (combine_energy_bin)? '_energy_combined': ''
foreach species, test_species do begin
    case species of
        'o': zrange = [5e3,5e5]
        'p': zrange = [1e4,1e6]
        'e': zrange = [5e5,1e10]
        'he': zrange = [1e3,1e5]
    endcase

    foreach unit, units do begin
        the_dir = join_path([plot_dir,species+'_'+unit+suffix])
        var = rbsp_plot_pa2d(time_range, probe=probe, species=species, unit=unit, $
            zrange=zrange, contour=0, plot_dir=the_dir, combine_energy_bin=combine_energy_bin)
    endforeach
endforeach

stop


units = ['energy','velocity']
test_species = ['p','o','e','he']

if n_elements(event_info) eq 0 then event_info = _2013_0501_0850_load_data()
prefix = event_info['prefix']
probe = event_info['probe']
time_range = event_info['time_range']
time_range = time_double(['2013-05-01/08:30','2013-05-01/09:30'])
plot_dir = join_path([event_info['plot_dir'],'rbsp_hope','rbsp'+probe])

foreach species, test_species do begin
    case species of
        'o': zrange = [5e3,5e5]
        'p': zrange = [1e4,1e6]
        'e': zrange = [5e5,1e10]
        'he': zrange = [1e3,1e5]
    endcase

    foreach unit, units do begin
        foreach combine_energy_bin, [0,1] do begin
            suffix = (combine_energy_bin)? '_energy_combined': ''
            the_dir = join_path([plot_dir,species+'_'+unit+suffix])
            var = rbsp_plot_pa2d(time_range, probe=probe, species=species, unit=unit, $
                zrange=zrange, contour=0, plot_dir=the_dir, combine_energy_bin=combine_energy_bin)
        endforeach
    endforeach
endforeach

end