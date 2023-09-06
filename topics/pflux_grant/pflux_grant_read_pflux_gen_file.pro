;+
; Preprocess pflux for a given time range.
; Load data by day by checking buffer.
; Process data by orbit.
; Save data by day.
;-

pro pflux_grant_read_pflux_gen_file, time, probe=probe, filename=file, local_root=local_root

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
    day_time_range = date+[0,secofday]
    pad_time = 1800.
    time_range = day_time_range+[-1,1]*pad_time
    project = pflux_grant_load_project()
    common_time_step = project.common_time_step
    common_times = make_bins(time_range, common_time_step)
    ncommon_time = n_elements(common_times)
    ndim = 3
    prefix = 'rbsp'+probe+'_'
    xyz = constant('xyz')
    fac = ['earthward','westward','outward']
    fillval = !values.f_nan


;---Read B0, B1 GSM, E, Edot0 MGSE.
    b0_gsm_var = prefix+'b0_gsm'
    pflux_grant_read_bfield, time_range, probe=probe
    pflux_grant_read_efield, time_range, probe=probe
    pflux_grant_read_efield, time_range, probe=probe, id='edot0'


;---Normalization.
    rbsp_read_orbit, time_range, probe=probe
    r_var = prefix+'r_gse'
    r_gse = get_var_data(r_var, times=orbit_times)
    r_gsm = cotran(r_gse, orbit_times, 'gse2gsm')
    store_data, prefix+'r_gsm', orbit_times, r_gsm
    add_setting, prefix+'r_gsm', /smart, dictionary($
        'diaplay_type', 'vector', $
        'short_name', 'R', $
        'unit', 'Re', $
        'coord', 'GSM', $
        'coord_labels', xyz )
    norbit_time = n_elements(orbit_times)

    ndim = 3
    par = 2
    model = 't89'
    re = constant('re') ; in km.
    r0 = 100d/re+1
    bf_gsm = fltarr(norbit_time,ndim)
    r_sm = cotran(r_gsm, orbit_times, 'gsm2sm')
    dirs = float(r_sm[*,2] lt 0)*2-1    ; n-hem: -1, s-hem: 1.
    foreach time, orbit_times, ii do begin
        tilt = geopack_recalc(time)
        rx = r_gsm[ii,0]
        ry = r_gsm[ii,1]
        rz = r_gsm[ii,2]
        dir = dirs[ii]
        geopack_trace, rx,ry,rz, dir, par, fx,fy,fz, r0=r0, $
            /refine, /ionosphere, /t89
        f_gsm = [fx,fy,fz]
        geopack_igrf_gsm, fx,fy,fz, bx,by,bz
        bf_gsm[ii,*] = [bx,by,bz]
    endforeach

    cnorm = snorm(bf_gsm)/interpol(snorm(get_var_data(b0_gsm_var)),common_times,orbit_times)
    cnorm_var = prefix+'cnorm'
    store_data, cnorm_var, orbit_times, cnorm



;---Calc dE, dB FAC.
    define_fac, prefix+'b0_gsm', prefix+'r_gsm', time_var=prefix+'r_gsm'
    q_var = prefix+'q_gsm2fac'
    get_data, q_var, times, q_gsm2fac
    store_data, q_var+'_lowres', times, q_gsm2fac
    q_gsm2fac = qslerp(q_gsm2fac, times, common_times)
    store_data, q_var, common_times, q_gsm2fac

    foreach var, prefix+['e','edot0'] do begin
        vec = get_var_data(var+'_mgse')
        vec = cotran(vec, common_times, 'mgse2gsm', probe=probe)
        store_data, var+'_gsm', common_times, vec
        add_setting, var+'_gsm', /smart, dictionary($
            'display_type', 'vector', $
            'short_name', 'E', $
            'unit', 'mV/m', $
            'coord', 'GSM', $
            'coord_labels', xyz )
    endforeach

    foreach var, prefix+['b1','e','edot0'] do begin
        to_fac, var+'_gsm', to=var+'_fac', q_var=q_var
    endforeach



;---pflux in the wanted frequency range.
    if n_elements(project) eq 0 then project = pflux_grant_load_project()
    pflux_calc_setting = project['pflux_calc_setting']
    scale_info = pflux_calc_setting.scale_info
    fac = ['||','west','outward']
    rgb = constant('rgb')
    de_types = ['','dot0']
    db_var = prefix+'b1_fac'
    foreach de_type, de_types do begin
        de_var = prefix+'e'+de_type+'_fac'
        pf_var = prefix+'pf'+de_type+'_fac'
        stplot_calc_pflux_mor, de_var, db_var, pf_var, scaleinfo=scale_info
        add_setting, pf_var, /smart, {$
            display_type: 'vector', $
            short_name: 'S', $
            unit: 'mW/m!U2!N', $
            coord: 'FAC', $
            coord_labels: fac, $
            colors: rgb}

        pf_norm = get_var_data(pf_var, times=common_times, limits=lim)
        cnorm = get_var_data(cnorm_var, at=common_times)
        for ii=0, ndim-1 do pf_norm[*,ii] *= cnorm
        index = where(interpol(r_sm[*,2],orbit_times,common_times) lt 0, count)
        if count ne 0 then pf_norm[index,0] *= -1
        pf_map_var = prefix+'pf'+de_type+'_fac_norm'
        store_data, pf_map_var, common_times, pf_norm, limits=lim
        options, pf_map_var, 'labels', 'FAC S!D'+fac
    endforeach



;---Save data to file.
    cdf_id = (file_test(file))? cdf_open(file): cdf_create(file)

    time_var = 'ut_sec'
    time_index = where_pro(common_times, '[)', day_time_range)
    times = common_times[time_index]
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var


    ; Pflux.
    foreach var, prefix+'pf'+['','dot0']+'_fac_norm' do begin
        short_name = (var eq prefix+'pf_fac_norm')? 'S': 'S!S!UDot0!N!R'
        data = float((get_var_data(var))[time_index,*])
        settings = dictionary($
            'depend_0', time_var, $
            'display_type', 'vector', $
            'short_name', short_name, $
            'unit', 'mW/m!U2!N', $
            'coord', 'FAC', $
            'coord_labels', fac)
        cdf_save_var, var, value=data, filename=cdf_id
        cdf_save_setting, settings, filename=cdf_id, varname=var
    endforeach


    ; Orbit time.
    time_var = 'ut_orbit'
    time_index = where_pro(orbit_times, '[]', day_time_range)
    times = orbit_times[time_index]
    time_var_settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix')
    cdf_save_var, time_var, value=times, filename=cdf_id
    cdf_save_setting, time_var_settings, filename=cdf_id, varname=time_var

    ; Quaternion.
    var = prefix+'q_gsm2fac'
    data = float((get_var_data(var+'_lowres'))[time_index,*])
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'vector', $
        'short_name', 'Q', $
        'unit', '#', $
        'coord', 'GSM2FAC', $
        'coord_labels', ['a','b','c','d'], $
        'colors', sgcolor(['red','green','blue','black']) )
    cdf_save_var, var, value=data, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=var
        
    ; cnorm.
    var = cnorm_var
    data = float((get_var_data(var))[time_index,*])
    settings = dictionary($
        'depend_0', time_var, $
        'display_type', 'scalar', $
        'short_name', 'C!S!Unorm!N!R!DT89!N', $
        'unit', '#' )
    cdf_save_var, var, value=data, filename=cdf_id
    cdf_save_setting, settings, filename=cdf_id, varname=var

    cdf_close, cdf_id


end




probe = 'b'
time = time_double('2013-06-07')
file = join_path([homedir(),'test.cdf'])
if file_test(file) eq 1 then file_delete, file
tic
pflux_grant_read_pflux_gen_file, time, probe=probe, filename=file
toc
end
