;+
; Load all OMNI data, check for invalid data.
;-

pro cusp_ml_preprocess_omni

    data_file = cusp_ml_preprocess_data_file()
    common_times = cusp_ml_preprocess_common_times(time_var=time_var)
    time_range = minmax(common_times)
    vars = ['sw_b_sm','sw_v_sm', $
        'sw_dens','sw_pdyn','sw_temp']
    old_vars = [['bx','by','bz']+'_gsm',$
        ['vx','vy','vz']+'_gse','n','p','t']

    load_data = 0
    foreach var, vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach
    if load_data then begin
    ;---Read data.
        if check_if_update(old_vars[0], time_range) then begin
            omni_read, time_range, id='sw'

            ; These are 1 min data, checked to have valid values.
;            foreach var, old_vars do begin
;                interp_time, var, common_times
;            endforeach
        endif

        b_gsm = [$
            [get_var_data('bx_gsm')], $
            [get_var_data('by_gsm')], $
            [get_var_data('bz_gsm')]]
        b_sm = cotran(b_gsm, common_times, 'gsm2sm')
        store_data, 'sw_b_sm', common_times, b_sm
        add_setting, 'sw_b_sm', /smart, dictionary($
            'display_type', 'vector', $
            'short_name', 'SW B', $
            'unit', 'nT', $
            'coord', 'SM', $
            'coord_labels', ['x','y','z'])

        v_gse = [$
            [get_var_data('vx_gse')], $
            [get_var_data('vy_gse')], $
            [get_var_data('vz_gse')]]
        v_sm = cotran(v_gse, common_times, 'gse2sm')
        store_data, 'sw_v_sm', common_times, v_sm
        add_setting, 'sw_v_sm', /smart, dictionary($
            'display_type', 'vector', $
            'short_name', 'SW V', $
            'unit', 'km/m', $
            'coord', 'SM', $
            'coord_labels', ['x','y','z'])

        copy_data, 'n', 'sw_dens'
        copy_data, 't', 'sw_temp'
        copy_data, 'p', 'sw_pdyn'
        get_data, 'sw_temp', times, temp
        store_data, 'sw_temp', times, temp/11604.45
        
        options, 'sw_dens', 'unit', 'cm!U-3!N'
        options, 'sw_temp', 'unit', 'eV'
        options, 'sw_pdyn', 'unit', 'nPa'

        foreach var, 'sw_'+['dens','temp','pdyn'] do begin
            add_setting, var, /smart, dictionary($
                'diaplay_type', 'scalar' )
        endforeach


        ntime = n_elements(common_times)
        foreach var, vars do begin
            data = get_var_data(var, limits=settings)
            if n_elements(data[*,0]) ne ntime then stop
            cdf_save_var, var, value=data, filename=data_file
            settings = dictionary(settings)
            settings['depend_0'] = time_var
            settings['VAR_TYPE'] = 'data'
            cdf_save_setting, settings, filename=data_file, varname=var
        endforeach
    endif

    foreach var, vars do begin
        cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    endforeach



end
