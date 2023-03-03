


units = ['energy','velocity']
test_species = ['p','o','e','he']

if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
prefix = event_info['prefix']
probe = event_info['probe']
time_range = event_info['time_range']
time_range = time_double(['2013-05-01/07:38','2013-05-01/07:50'])
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




; ; Check Chaston+2015 event on 2013-06-01 RBSP-A.
; time_range = time_double(['2013-06-01/03:00','2013-06-01/08:00'])
; probe = 'a'
; 
; log = 1
; unit = 'energy'
; species = ['electron','proton','oxygen']
; zrngs = [[4,8],[3.5,6],[3.5,6]]
; 
; hopel3 = sread_rbsp_hope_l3(time_range, probes=probe)
; times = make_bins(time_range, 12, inner=1)
; 
; foreach species_name, species, type_id do begin
;     foreach time, times do begin
;         plot_hope_l3_pitch2d, time, species_name, unit=unit, $
;             log=log, hopel3=hopel3, probe=probe, zrange=zrngs[*,type_id]
;     endforeach
; endforeach
