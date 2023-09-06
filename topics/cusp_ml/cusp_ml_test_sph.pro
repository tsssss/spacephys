;+
; Preprocess features and target to ensure all invalid data are marked by NaN.
;-

pro cusp_ml_test_sph, filename=filename

    root_dir = join_path([googledir(),'works','cusp_ml','data'])
    filename = join_path([root_dir,'test_sph.cdf'])
    if n_elements(filename) eq 0 then stop

;---These data should have invalid data marked by NaN.
    routines = ['density','r_sm','r_sm_sph','dst','omni']
    routines = 'cusp_ml_preprocess_'+routines
    foreach routine, routines do begin
        call_procedure, routine
    endforeach


;---Process data to make features and target.
    get_data, 'density', times, data
    store_data, 'density', times, alog10(data)

    stplot_split, 'sw_b_sm', newnames='sw_b'+['x','y','z']
    stplot_split, 'sw_v_sm', newnames='sw_v'+['x','y','z']
    get_data, 'sw_v_sm', times, data
    store_data, 'sw_vmag', times, snorm(data)
    add_setting, 'sw_vmag', /smart, dictionary($
        'display_type', 'scalar', $
        'short_name', 'SW |V|', $
        'unit', 'km/s' )


;---Mark all invalid data and remove them.
    target_var = 'density'
    feature_vars = ['dst','dis','lat','lon',$
        'sw_bx','sw_by','sw_bz','sw_vx','sw_vy','sw_vz',$
        'sw_vmag','sw_dens','sw_temp','sw_pdyn']
    all_vars = [target_var,feature_vars]
    fillval = !values.f_nan
    get_data, target_var, common_times
    ntime = n_elements(common_times)
    foreach var, all_vars do begin
        get_data, var, tmp, data
        index2 = where(finite(data), complement=index)
        common_times[index] = fillval
    endforeach
    get_data, 'r_sm', tmp, data
    index = where(data[*,0] le 0)   ; only consider x>0.
    common_times[index] = fillval

;---Test.
;    get_data, 'density', tmp, dens
;    get_data, 'lat', tmp, lat
;    dens[*] = 1e-2
;    index = where_pro(lat, '[]', [50,60])
;    dens[index] = 1e2
;    index = where_pro(lat, '[]', [-70,-60])
;    dens[index] = 1e1
;    store_data, 'density', tmp, dens

    good_index = where(finite(common_times))
    foreach var, all_vars do begin
        get_data, var, times, data
        store_data, var, times[good_index], data[good_index]
    endforeach

    stplot2cdf, all_vars, filename=filename, time_var='time'

end
