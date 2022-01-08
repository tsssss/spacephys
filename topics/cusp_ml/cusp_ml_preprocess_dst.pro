;+
; Load all Dst data, check for invalid data.
;-

pro cusp_ml_preprocess_dst

    data_file = cusp_ml_preprocess_data_file()
    var = 'dst'
    common_times = cusp_ml_preprocess_common_times(time_var=time_var)
    time_range = minmax(common_times)
    if ~cdf_has_var(var, filename=data_file) then begin
    ;---Read data.
        if check_if_update(var, time_range) then begin
            omni_read_index, time_range
            interp_time, var, common_times
        endif

        get_data, var, common_times, dst
        cdf_save_var, var, value=dst, filename=data_file
        settings = dictionary($
            'VAR_TYPE', 'data', $
            'depend_0', time_var, $
            'display_type', 'scalar', $
            'unit', 'nT' )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'

end

cusp_ml_preprocess_dst
end