;+
; Load all density data, check for invalid data.
;-

pro cusp_ml_polar_test_preprocess_orbit

    probe = 'polar'

;---Read data.
    root_dir = join_path([googledir(),'works','works','global_efield','data'])

    files = join_path([root_dir,'orbit',probe+'_orbit_*.tplot'])
    files = file_search(files[0])
    nfile = n_elements(files)
    new_var = probe+'_r_sm'

    common_times = cusp_ml_polar_test_common_time()
    time_range = minmax(common_times)
    if check_if_update(new_var, time_range) then begin
        old_var = 'po_r_gsm'
        foreach file, files do begin
            del_data, old_var
            tplot_restore, filename=file
            get_data, old_var, times, data
            append_data, new_var, times, data
        endforeach
        get_data, new_var, times, r_gsm
        r_sm = cotran(r_gsm, times, 'gsm2sm')
        store_data, new_var, times, r_sm
        add_setting, new_var, /smart, dictionary($
            'display_type', 'vector', $
            'short_name', 'R', $
            'unit', 'Re', $
            'coord', 'SM', $
            'coord_labels', ['x','y','z'] )
        interp_time, new_var, common_times
    endif

    return

end

cusp_ml_polar_test_preprocess_orbit
end
