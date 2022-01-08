;+
; Test spin phase with different padding.
;-


;---Input.
    day = time_double('2014-08-28')
    probe = 'b'

;---Settings.
    day_time_range = day+[0,86400d]
    pad_times = [0,300,3600d]
    prefix = 'rbsp'+probe+'_'
    ndim = 3
    uvw = constant('uvw')
    xyz = constant('xyz')
    fillval = !values.f_nan
    common_time_step = 1d


    comp_var = prefix+'b_uvw'
    foreach pad_time, pad_times do begin
    ;---Read data.
        time_range = day_time_range+[-1,1]*pad_time
        common_times = make_bins(time_range, common_time_step)
        rbsp_read_emfisis, time_range, probe=probe, id='l2%magnetometer'

        interp_time, comp_var, common_times
        add_setting, comp_var, /smart, dictionary($
            'display_type', 'vector', $
            'unit', 'nT', $
            'short_name', 'B', $
            'coord', 'UVW', $
            'coord_labels', uvw )
        
        the_var = comp_var+'_'+string(pad_time,format='(I0)')+'_sec'
        tplot_rename, comp_var, the_var
    endforeach


    vars = comp_var+'_'+string(pad_times,format='(I0)')+'_sec'
    the_time_range = day_time_range
    data0 = get_var_data(vars[0], in=the_time_range, times=times)
    for jj=1,n_elements(vars)-1 do begin
        ddata = get_var_data(vars[jj], at=times)-data0
        store_data, comp_var+'_diff'+string(jj,format='(I0)'), times, ddata[*,0]
    endfor

end
