;+
; Load all Dst data, check for invalid data.
;-

pro cusp_ml_polar_test_preprocess_dst

    probe = 'polar'
    new_var = 'dst'

;---Read data.
    common_times = cusp_ml_polar_test_common_time()
    time_range = minmax(common_times)
    if check_if_update(new_var, time_range) then begin
        omni_read_index, time_range
        interp_time, new_var, common_times
    endif
    return

    get_data, new_var, times, dat
    tpos = panel_pos(fig_size=fig_size, pansize=[4,3])
    sgopen, 0, xsize=fig_size[0], ysize=fig_size[1]
    cghistoplot, dat, binsize=0.1, /fill, position=tpos, $
        xtitle='Dst (nT)', ytitle='Count (#)'

end
