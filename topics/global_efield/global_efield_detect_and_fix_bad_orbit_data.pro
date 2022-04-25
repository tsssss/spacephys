;+
; Detect bad orbit data in the combined_file.
;
; Check for:
;   1. Missing data point, should be 1 min cadence.
;   2. Check for discontinuities, through checking large velocity.
;-

pro global_efield_detect_and_fix_bad_orbit_data, mission_probe, project=project

    global_efield_load_data, 'r_gsm', probe=mission_probe, project=project

    prefix = project[mission_probe].prefix
    var = prefix+'r_gsm'
    
    lprmsg, 'Detecting bad data in '+var+' ...'


;---Check cadence.
    get_data, var, times, orb_data
    dtimes = times[1:*]-times[0:-2]
    time_step = project['time_step']
    min_dtime = min(dtimes)
    if min_dtime ne time_step then message, 'Min dtime is: '+min_dtime
    max_dtime = min(dtimes)
    if max_dtime ne time_step then message, 'Min dtime is: '+max_dtime
    lprmsg, '    Cadence is good ...'

;---Check discontinuities.
    ntime = n_elements(times)
    ndim = 3
    vel_data = fltarr(ntime,ndim)
    re = 6378.
    for ii=0, ndim-1 do vel_data[*,ii] = deriv(orb_data[*,ii])/time_step*re
    vel_var = prefix+'v_gsm'
    store_data, vel_var, times, vel_data
    add_setting, vel_var, /smart, {$
        display_type: 'vector', $
        unit: 'km/s', $
        short_name: 'V', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}

    full_time_range = minmax(times)
    ;tplot, [var,vel_var], trange=full_time_range

    max_vel = 15.
    vmag = snorm(vel_data)
    bad_points = where(vmag ge max_vel, count)

    if count ne 0 then begin
        lprmsg, '    Discontinuity found ...'
        
        ; C2 have sections of step up/down in Ry and Rz.
        if mission_probe eq 'c2' then begin
            section_length = 600
            for ii=0, count-1 do begin
                ; Locate a section.
                pos0 = bad_points[ii]
                i1 = ii+1
                while i1 lt count-1 do begin
                    if bad_points[i1+1]-bad_points[i1] ge section_length then break
                    i1 += 1
                endwhile
                ; Found a section.
                pos1 = bad_points[i1]
                orb_data[pos0:pos1,1:2] = !values.f_nan
                ; Next section.
                ii = i1
            endfor
        ;---Fix Ry and Rz.
            for ii=1,2 do begin
                index = where(finite(orb_data[*,ii]))
                orb_data[*,ii] = interpol(orb_data[index,ii], times[index], times, /quadratic)
            endfor
        endif else begin
            ; Works for Polar, TH-[BC]. There are singular points of bad data.
            if count lt ntime then begin
                bad_points = where(vmag ge max_vel, count)
                npad = 5
                for ii=0, count-1 do begin
                    i0 = (bad_points[ii]-npad)>(0)
                    i1 = (bad_points[ii]+npad)<(ntime-1)
                    vmag[i0:i1] = !values.f_nan
                endfor
                index = where(finite(vmag))
                orb_data = sinterpol(orb_data[index,*], times[index], times, /quadratic)
            endif
        endelse

    ;---Update data.
        fix_var = prefix+'r_fix_gsm'
        store_data, fix_var, times, orb_data
        add_setting, fix_var, /smart, {$
            display_type: 'vector', $
            unit: 'Re', $
            short_name: 'R', $
            coord: 'GSM', $
            coord_labels: ['x','y','z']}
        tplot, [var, fix_var, vel_var], trange=full_time_range

    ;---Save results.
        cdf_file = join_path([project.data_dir,project[mission_probe].combined_file_suffix])
        the_var = prefix+'r_gsm'
        settings = cdf_read_setting(the_var, filename=cdf_file)
        get_data, fix_var, times, orb_data
        cdf_save_var, the_var, value=orb_data, filename=cdf_file
        cdf_save_setting, settings, varname=the_var, filename=cdf_file
        lprmsg, '    Bad points fixed ...'
    endif else lprmsg, '    No discontinuity found ...'


end