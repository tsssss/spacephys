;+
; Load all orbit data, check for invalid data.
;-

pro cusp_ml_preprocess_r_sm

    data_file = cusp_ml_preprocess_data_file()
    var = 'r_sm'
    common_times = cusp_ml_preprocess_common_times(time_var=time_var)
    time_range = minmax(common_times)
    if ~cdf_has_var(var, filename=data_file) then begin
    ;---Read data.
        if check_if_update(var, time_range) then begin
            root_dir = join_path([googledir(),'works','works','global_efield','data'])

            probe = 'polar'
            files = join_path([root_dir,'orbit',probe+'_orbit_*.tplot'])
            files = file_search(files[0])
            nfile = n_elements(files)

            old_var = 'po_r_gsm'
            foreach file, files do begin
                del_data, old_var
                tplot_restore, filename=file
                get_data, old_var, times, data
                append_data, var, times, data
            endforeach
            get_data, var, times, r_gsm
            r_sm = cotran(r_gsm, times, 'gsm2sm')
            store_data, var, times, r_sm
            add_setting, var, /smart, dictionary($
                'display_type', 'vector', $
                'short_name', 'R', $
                'unit', 'Re', $
                'coord', 'SM', $
                'coord_labels', ['x','y','z'] )
            interp_time, var, common_times
        endif

        get_data, var, common_times, r_sm
        cdf_save_var, var, value=r_sm, filename=data_file
        settings = dictionary($
            'VAR_TYPE', 'data', $
            'depend_0', time_var, $
            'display_type', 'vector', $
            'short_name', 'R', $
            'unit', 'Re', $
            'coord', 'SM', $
            'coord_labels', ['x','y','z'] )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'

end

cusp_ml_preprocess_r_sm
end
