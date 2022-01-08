;+
; Load orbit data and then calc dis, lat, lon.
;-

pro cusp_ml_preprocess_r_sm_sph

    data_file = cusp_ml_preprocess_data_file()
    tplot_var = 'r_sm_sph'
    common_times = cusp_ml_preprocess_common_times(time_var=time_var)
    time_range = minmax(common_times)
    vars = ['dis','lat','lon']
    load_data = 0
    foreach var, vars do begin
        if ~cdf_has_var(var, filename=data_file) then begin
            load_data = 1
            break
        endif
    endforeach
    if load_data then begin
    ;---Read data.
        if check_if_update(tplot_var, time_range) then begin
            cusp_ml_preprocess_r_sm

            deg = constant('deg')
            get_data, 'r_sm', common_times, r_sm
            dis = snorm(r_sm)
            lat = asin(r_sm[*,2]/dis)*deg
            lon = atan(r_sm[*,1], r_sm[*,0])*deg

            r_sm_sph = [[dis],[lat],[lon]]
            store_data, tplot_var, common_times, r_sm_sph
        endif

        get_data, tplot_var, common_times, r_sm_sph
        dis = r_sm_sph[*,0]
        cdf_save_var, 'dis', value=dis, filename=data_file
        settings = dictionary($
            'VAR_TYPE', 'data', $
            'depend_0', time_var, $
            'short_name', 'dis', $
            'display_type', 'scalar', $
            'unit', 'Re' )
        cdf_save_setting, settings, filename=data_file, varname=vars[0]

        lat = r_sm_sph[*,1]
        cdf_save_var, 'lat', value=lat, filename=data_file
        settings = dictionary($
            'VAR_TYPE', 'data', $
            'depend_0', time_var, $
            'display_type', 'scalar', $
            'short_name', 'latitude', $
            'unit', 'deg' )
        cdf_save_setting, settings, filename=data_file, varname=vars[1]

        lon = r_sm_sph[*,2]
        cdf_save_var, 'lon', value=lon, filename=data_file
        settings = dictionary($
            'VAR_TYPE', 'data', $
            'depend_0', time_var, $
            'display_type', 'scalar', $
            'short_name', 'longitude', $
            'range', [-180,180], $
            'notes', '0 deg at noon, continuous around noon', $
            'unit', 'deg' )
        cdf_save_setting, settings, filename=data_file, varname=vars[2]
    endif

    foreach var, vars do begin
        cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'
    endforeach

end
