

units = ['energy','velocity']
test_species = ['p','o','e','he']
id = '2015_0312_0900'

if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)
rbsp_info = event_info.rbsp.rbspb
prefix = rbsp_info['prefix']
probe = rbsp_info['probe']
time_range = event_info['time_range']+[1,-1]*1800
plot_dir = join_path([event_info['plot_dir'],'rbsp_hope_pa','rbsp'+probe])

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