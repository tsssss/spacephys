;+
; Load all density data, check for invalid data.
;-

pro cusp_ml_preprocess_density

    data_file = cusp_ml_preprocess_data_file()
    var = 'density'
    common_times = cusp_ml_preprocess_common_times(time_var=time_var)
    time_range = minmax(common_times)
    if ~cdf_has_var(var, filename=data_file) then begin
    ;---Read data.
        if check_if_update(var, time_range) then begin
            root_dir = join_path([googledir(),'works','works','global_efield','data'])

            probe = 'polar'
            files = join_path([root_dir,'density',probe+'_density_*.tplot'])
            files = file_search(files[0])
            nfile = n_elements(files)

            old_var = 'po_ele_n'
            foreach file, files do begin
                del_data, old_var
                tplot_restore, filename=file
                get_data, old_var, times, data
                append_data, var, times, data
            endforeach
            add_setting, var, /smart, dictionary($
                'display_type', 'scalar', $
                'short_name', 'N', $
                'unit', 'cm!U-3!N', $
                'ylog', 1 )
            interp_time, var, common_times

            dens = get_var_data(var)
            dat = alog10(dens)
            index = where(finite(dat), complement=index2)
            dens[index2] = !values.f_nan
            store_data, var, common_times, dens
        endif

        get_data, var, common_times, dens
        cdf_save_var, var, value=dens, filename=data_file
        settings = dictionary($
            'VAR_TYPE', 'data', $
            'depend_0', time_var, $
            'display_type', 'scalar', $
            'unit', 'cm!U-3!N', $
            'ylog', 1 )
        cdf_save_setting, settings, filename=data_file, varname=var
    endif

    cdf_load_var, var, filename=data_file, time_var=time_var, time_type='unix'

end

cusp_ml_preprocess_density
end
