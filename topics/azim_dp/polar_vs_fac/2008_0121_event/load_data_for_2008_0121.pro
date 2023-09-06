;+
; Load data for the 2008_0121 event.
;-

pro load_data_for_2008_0121_footpoint, default_settings


;---Settings.
    time_range = default_settings.event_time_range+[-1,1]*1800d
    probes = themis_probes()
    probes = ['a']
    models = ['t89','t96','t01','t04s']

    foreach probe, probes do begin
        prefix = 'th'+probe+'_'
        r_sm_var = prefix+'r_sm'
        get_data, r_sm_var, times, r_sm, limits=lim
        r_gsm = cotran(r_sm, times, 'sm2gsm')
        r_gsm_var = prefix+'r_gsm'        
        store_data, r_gsm_var, times, r_gsm
        add_setting, r_gsm_var, smart=1, dictionary($
            'short_name', 'R', $
            'unit', 'Re', $
            'display_type', 'vector', $
            'coord', 'GSM' )
        foreach model, models do begin
            read_geopack_info, r_gsm_var, model=model, prefix=prefix, suffix='_'+model, direction=-1
        endforeach
        
        colors = sgcolor(['black','red','green','blue'])
        foreach type, ['fmlt','fmlat','fmlon'] do begin
            var = prefix+type
            stplot_merge, prefix+type+'_'+models, newname=var, labels=strupcase(models), colors=colors
            options, var, 'ynozero', 1
            options, var, 'labflag', -1
        endforeach
    endforeach
    

end


pro load_data_for_2008_0121_themis_data, default_settings

;---Settings.
    time_range = default_settings.event_time_range+[-1,1]*1800d
    probes = themis_probes()
    coord = 'sm'
    var = 'themis_data'
    energy_range = [0,400d]
    themis_settings = dictionary($
        'time_range', time_range, $
        'energy_range', energy_range, $
        'var', var )
    default_settings[var] = themis_settings


    time_var = 'time_themis_data'
    data_file = default_settings.data_file
    if ~cdf_has_var(time_var, filename=data_file) then begin
        time_step = 3d
        common_times = make_bins(time_range, time_step)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=common_times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
    endif
    common_times = cdf_read_var(time_var, filename=data_file)


;---Loop through each probe.
    foreach probe, probes do begin
        prefix = 'th'+probe+'_'

    ;---Orbit.
        r_var = themis_read_orbit(time_range, probe=probe, coord=coord, get_name=1)
        if ~cdf_has_var(r_var, filename=data_file) then begin
            r_var = themis_read_orbit(time_range, probe=probe, coord=coord)
            interp_time, r_var, common_times
            get_data, r_var, times, r_sm, limits=lim
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            var = r_var
            cdf_save_var, var, value=r_sm, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, r_var, filename=data_file, time_var=time_var, time_type='unix'

    ;---MLT.
        mlt_var = prefix+'mlt'
        if ~cdf_has_var(mlt_var, filename=data_file) then begin
            get_data, r_var, times, r_sm
            mlt = pseudo_mlt(r_sm)
            settings = dictionary($
                'display_type', 'scalar', $
                'depend_0', time_var, $
                'unit', 'hr', $
                'short_name', 'MLT' )
            var = mlt_var
            cdf_save_var, var, value=mlt, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, mlt_var, filename=data_file, time_var=time_var, time_type='unix'
        

    ;---B field.
        b_var = themis_read_bfield(time_range, probe=probe, coord=coord, get_name=1)
        if ~cdf_has_var(b_var, filename=data_file) then begin
            b_var = themis_read_bfield(time_range, probe=probe, coord=coord)
            interp_time, b_var, common_times
            get_data, b_var, times, b_sm, limits=lim
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            var = b_var
            cdf_save_var, var, value=b_sm, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, b_var, filename=data_file, time_var=time_var, time_type='unix'
        
    ;---keV flux.
        if probe eq 'b' then continue
        e_energy_var = 'e_energy_'+probe
        e_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='e', get_name=1)
;cdf_del_var, e_flux_var, filename=data_file
        if ~cdf_has_var(e_flux_var, filename=data_file) then begin 
            e_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='e')
            interp_time, e_flux_var, common_times
            get_data, e_flux_var, times, flux, energys, limits=lim

            var = e_energy_var
            settings = dictionary($
                'unit', 'keV', $
                'data_type', 'metadata' )
            cdf_save_var, var, value=energys, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
            
            var = e_flux_var
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            settings['depend_1'] = e_energy_var
            cdf_save_var, var, value=flux, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, e_flux_var, filename=data_file, time_var=time_var, time_type='unix'


        p_energy_var = 'p_energy_'+probe
        p_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='p', get_name=1)
;cdf_del_var, p_flux_var, filename=data_file
        if ~cdf_has_var(p_flux_var, filename=data_file) then begin 
            p_flux_var = themis_read_kev_flux(time_range, probe=probe, no_spec=1, id='p')
            interp_time, p_flux_var, common_times
            get_data, p_flux_var, times, flux, energys, limits=lim

            var = p_energy_var
            settings = dictionary($
                'unit', 'keV', $
                'data_type', 'metadata' )
            cdf_save_var, var, value=energys, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
            
            var = p_flux_var
            settings = dictionary(lim)
            settings['depend_0'] = time_var
            settings['depend_1'] = p_energy_var
            cdf_save_var, var, value=flux, filename=data_file
            cdf_save_setting, settings, varname=var, filename=data_file
        endif
        cdf_load_var, p_flux_var, filename=data_file, time_var=time_var, time_type='unix'


    ;---Make some adjustments.
        flux_vars = prefix+['e','p']+'_flux'
        foreach flux_var, flux_vars do begin
            get_data, flux_var, times, fluxes, energys, limits=lim
            index = where_pro(energys, '[]', energy_range, count=count)
            if count eq 0 then message, 'Inconsistency...'
            store_data, flux_var, times, fluxes[*,index], energys[index], limits=lim
            options, flux_var, 'labels', lim.labels[index]
            options, flux_var, 'colors', lim.colors[index]
        endforeach
    endforeach

end


pro load_data_for_2008_0121_thg_mlt_image, default_settings

;---Settings.
    time_range = default_settings.event_time_range
    sites = ['fsim','fsmi','whit','gill']
    min_elev = 2.5  ; deg
    mlat_range = [50d,90]   ; deg
    merge_method = 'merge_elev'
    var = 'thg_asf_mlt_image'
    time_var = 'time_'+var
    themis_asi_settings = dictionary($
        'time_range', time_range, $
        'sites', sites, $
        'min_elev', min_elev, $
        'mlat_range', mlat_range, $
        'merge_method', merge_method, $
        'var', var )
    default_settings[var] = themis_asi_settings
    if tnames(var) ne '' then return

    
;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        themis_asf_read_mlt_image, time_range, sites=sites, merge_method=merge_method, min_elev=min_elev
        var = 'thg_asf_mlt_image'
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim

end



pro load_data_for_2008_0121_j_ver_mlt_image, default_settings
    
;---Settings.
    time_range = default_settings.event_time_range
    mlat_range = [50d,90]   ; deg
    var = 'thg_j_ver_mlt_image'
    time_var = 'time_weygand_j'
    j_ver_settings = dictionary($
        'time_range', time_range, $
        'mlat_range', mlat_range, $
        'var', var )
    default_settings[var] = j_ver_settings
    if tnames(var) ne '' then return
    
;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        var = themis_read_j_ver_mlt_image(time_range)
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif
    
    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim
end


pro load_data_for_2008_0121_po_mlt_image, default_settings

;---Settings.
    time_range = default_settings.event_time_range
    mlat_range = [50d,90]   ; deg
    var = 'po_mlt_image'
    time_var = 'time_'+var
    themis_asi_settings = dictionary($
        'time_range', time_range, $
        'sites', sites, $
        'mlat_range', mlat_range, $
        'var', var )
    default_settings[var] = themis_asi_settings
    if tnames(var) ne '' then return

;---Load data.
    data_file = default_settings.data_file
    if ~cdf_has_var(var, filename=data_file) then begin
        polar_calc_mlt_image, time_range
        get_data, var, times, mlt_images, limits=lim
        settings = dictionary(lim)
        time_settings = dictionary('unit', 'sec', 'data_type', 'metadata')
        cdf_save_var, time_var, value=times, filename=data_file
        cdf_save_setting, time_settings, varname=time_var, filename=data_file
        settings['depend_0'] = time_var
        cdf_save_var, var, value=mlt_images, filename=data_file
        cdf_save_setting, settings, varname=var, filename=data_file
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    get_data, var, limits=lim
    settings = dictionary(lim)
    image_size = settings.image_size
    npixel = product(image_size)
    foreach key, settings.keys() do begin
        val = settings[key]
        if n_elements(val) eq npixel then val = reform(val,image_size)
        settings[key] = val
    endforeach
    lim = settings.tostruct()
    store_data, var, limits=lim

end


function load_data_for_2008_0121

    root_dir = join_path([googledir(),'works','polar_vs_fac'])
    data_dir = join_path([root_dir,'data'])
    data_file = join_path([data_dir,'2008_0121_data.cdf'])
    plot_dir = join_path([root_dir,'plot'])
    if file_test(data_file) eq 0 then cdf_touch, data_file

;---Global attributes.
    event_time_range = time_double(['2008-01-21/07:00','2008-01-21/10:00'])
    aurora_times = time_double('2008-01-21/'+['07:13','07:25','07:32','07:58','08:15','08:34','08:42','09:02'])
    default_settings = dictionary($
        'event_time_range', event_time_range, $
        'aurora_times', aurora_times, $
        'root_dir', root_dir, $
        'data_dir', data_dir, $
        'plot_dir', plot_dir, $
        'data_file', data_file )

;---Polar UVI.
    load_data_for_2008_0121_po_mlt_image, default_settings

;---THEMIS ASI.
    load_data_for_2008_0121_thg_mlt_image, default_settings

;---THEMIS current.
    load_data_for_2008_0121_j_ver_mlt_image, default_settings

;---THEMIS data.
    load_data_for_2008_0121_themis_data, default_settings

;---Load footpoints.
    load_data_for_2008_0121_footpoint, default_settings

    return, default_settings
    
end


tmp = load_data_for_2008_0121()
end