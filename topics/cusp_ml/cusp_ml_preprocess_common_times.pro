;+
; Get common_times saved and returned.
;-
function cusp_ml_preprocess_common_times, time_var=time_var

    data_file = cusp_ml_preprocess_data_file()
    time_var = 'time'

    load_data = 0
    if file_test(data_file) eq 0 then begin
        load_data = 1
    endif else begin
        if ~cdf_has_var(time_var, filename=data_file) then load_data = 1
    endelse

    if load_data then begin
        time_range = time_double(['1996-03-20','2008-06-14:23:59'])
        time_step = 60.
        common_times = make_bins(time_range, time_step)
        cdf_save_var, time_var, filename=data_file, value=common_times, cdf_type='CDF_DOUBLE'
        settings = dictionary($
            'VAR_TYPE', 'support_data', $
            'UNITS', 'sec', $
            'time_step', time_step, $
            'time_range', time_range, $
            'time_var_type', 'unix')
        cdf_save_setting, settings, filename=data_file, varname=time_var
    endif

    return, cdf_read_var(time_var, filename=data_file)

end
