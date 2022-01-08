;+
; Level 2 efield includes e_sm.
;
; time. A time in UT sec. Only time[0] is used to determine the day.
; probe=. A string of 'a' or 'b'.
; filename=. A string to set the output file name.
;-

pro pflux_grant_gen_level2_efield, time, probe=probe, $
    filename=file, errmsg=errmsg

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    if n_elements(file) eq 0 then begin
        errmsg = handle_error('No output file ...')
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif

    if n_elements(time) eq 0 then begin
        errmsg = handle_error('No input time ...')
        return
    endif


;---Settings.
    secofday = constant('secofday')
    date = time[0]-(time[0] mod secofday)
    time_range = date+[0,secofday]
    prefix = 'rbsp'+probe+'_'
    ndim = 3
    xyz = constant('xyz')
    spin_period = 11d   ; sec, only a approximation.


;---Read data: r_gsm, v_gsm, q_uvw2gsm, b_gsm, bmod_gsm, e_uv.
    rbsp_read_orbit, time_range, probe=probe
    rbsp_read_sc_vel, time_range, probe=probe
    rbsp_read_quaternion, time_range, probe=probe
    pflux_grant_read_level2_bfield, time_range, probe=probe
    pflux_grant_read_level1_data, time_range, probe=probe, id='efield'


    e_uv = get_var_data(prefix+'e_uv', times=common_times)
    ncommon_time = n_elements(common_times)
    common_time_step = total(common_times[0:1]*[-1,1])
    e_uvw = [[e_uv[*,0]],[e_uv[*,1]],[fltarr(ncommon_time)]]
    width = spin_period/common_time_step
    for ii=0,1 do begin
        offset1 = smooth(e_uvw[*,ii], width, /nan, /edge_zero)
        offset2 = smooth(offset1, width, /nan, /edge_zero)
        e_uvw[*,ii] -= offset2
    endfor
    store_data, prefix+'e_uvw', common_times, e_uvw
    add_setting, prefix+'e_uvw', /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E', $
        coord: 'UVW', $
        coord_labels: ['u','v','w'] }
    rbsp_uvw2gsm, prefix+'e_uvw', prefix+'e_gsm'
   


;---Convert to MGSE.
    foreach var, prefix+['e','v','r','b0'] do begin
        get_data, var+'_gsm', times, data, limits=limits
        data = cotran(data, times, 'gsm2mgse', probe=probe)
        store_data, var+'_mgse', times, data
        limits.coord = 'MGSE'
        add_setting, var+'_mgse', /smart, limits
    endforeach


;---E_vxb.
    de_var = prefix+'de_mgse'
    rbsp_remove_efield_bg, e_var=prefix+'e_mgse', b_var=prefix+'b0_mgse', $
        v_var=prefix+'v_mgse', r_var=prefix+'r_mgse', save_to=de_var, probe=probe
    get_data, de_var, times, e_mgse
    e_mgse[*,0] = 0
    store_data, de_var, times, e_mgse


;---Remove bad E field.
    rbsp_efw_read_boom_flag, time_range, probe=probe
    boom_flag = get_var_data(prefix+'boom_flag', times=uts)
    flags = total(boom_flag,2) ne 4
    index = where(flags eq 1, count)
    if count eq 0 then begin
        bad_time_ranges = !null
    endif else begin
        time_step = total(uts[0:1]*[-1,1])
        bad_time_ranges = time_to_range(uts[index], time_step=time_step)
        nbad_time_range = n_elements(bad_time_ranges)*0.5
        pad_time = 300. ; sec.
        for ii=0,nbad_time_range-1 do flags[lazy_where(uts,'[]',bad_time_ranges[ii,*]+[-1,1]*pad_time)] = 1
        index = where(flags eq 1)
        bad_time_ranges = time_to_range(uts[index], time_step=time_step)
    endelse
    nbad_time_range = n_elements(bad_time_ranges)*0.5
    get_data, de_var, common_times, de_mgse
    for ii=0,nbad_time_range-1 do de_mgse[lazy_where(common_times,'[]',bad_time_ranges[ii,*]),*] = !values.f_nan
    store_data, de_var, common_times, de_mgse
    


stop

end


time = time_double('2013-07-21')
time = time_double('2013-05-01')
;time = time_double('2013-06-07')
probe = 'a'
time = time_double('2013-07-16')
data_file = join_path([homedir(),'test.cdf'])
pflux_grant_gen_level2_efield, time, probe=probe, filename=data_file
end