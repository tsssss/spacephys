;+
; Load all density data, check for invalid data.
;-

pro cusp_ml_polar_test_preprocess_density

    probe = 'polar'

;---Read data.
    root_dir = join_path([googledir(),'works','works','global_efield','data'])

    files = join_path([root_dir,'density',probe+'_density_*.tplot'])
    files = file_search(files[0])
    nfile = n_elements(files)
    new_var = probe+'_density'

    common_times = cusp_ml_polar_test_common_time()
    time_range = minmax(common_times)
    if check_if_update(new_var, time_range) then begin
        old_var = 'po_ele_n'
        foreach file, files do begin
            del_data, old_var
            tplot_restore, filename=file
            get_data, old_var, times, data
            append_data, new_var, times, data
        endforeach
        add_setting, new_var, /smart, dictionary($
            'display_type', 'scalar', $
            'short_name', 'Ne', $
            'unit', 'cm!U-3!N', $
            'ylog', 1 )
        interp_time, new_var, common_times

        get_data, new_var, times, dens
        dat = alog10(dens)
        index = where(finite(dat), complement=index2)
        dens[index2] = !values.f_nan
        store_data, new_var, times, dens
    endif
    return


    tpos = panel_pos(fig_size=fig_size, pansize=[4,3])
    sgopen, 0, xsize=fig_size[0], ysize=fig_size[1]
    cghistoplot, dat, binsize=0.1, /fill, position=tpos, $
        xtitle='Log!D10!N Density (cm!U-3!N)', ytitle='Count (#)'

end

cusp_ml_polar_test_preprocess_density
end
