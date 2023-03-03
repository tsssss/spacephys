;+
; Try to plot and save pitch angle 2D in different ways.
;-

units = ['energy','velocity']
scale_methods = ['default']
test_species = ['p','o','e']
contour_flags = [0,1]

if n_elements(event_info) eq 0 then event_info = _2013_0501_load_data()
prefix = event_info['prefix']
probe = event_info['probe']
time_range = event_info['time_range']
plot_dir = join_path([event_info['plot_dir'],'pa_distr2d'])

foreach species, test_species do begin
    case species of
        'o': zrange = [5e3,5e5]
        'p': zrange = [1e4,1e6]
        'e': zrange = [5e5,1e10]
    endcase

    foreach unit, units do begin
        foreach scale_method, scale_methods do begin
            foreach use_contour, contour_flags do begin
                plot_type_str = (use_contour eq 1)? 'contour': 'polygon'
                the_dir = join_path([plot_dir,species])
                file_suffix = '_'+unit+'_'+scale_method+'_'+plot_type_str
                var = rbsp_plot_pa2d(time_range, probe=probe, species=species, $
                    scale_method=scale_method, unit=unit, zrange=zrange, contour=use_contour, plot_dir=the_dir, file_suffix=file_suffix)
            endforeach   
        endforeach
    endforeach
endforeach


end