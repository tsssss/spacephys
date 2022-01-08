;+
; Read level 2 data: MLT, B model SM, theta.
;-

pro azim_dp_read_level2_data_gen_file, time, probe=probe, filename=file, test=test

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


;---Constants and settings.
    secofday = constant('secofday')
    errmsg = ''
    xyz = constant('xyz')

    ; Derived settings.
    date = time[0]
    date = date-(date mod secofday)
    time_range = date+[0,secofday]
    prefix = probe+'_'
    probe_info = resolve_probe(probe)


;---Init file.
    out_dir = fgetpath(file)
    if file_test(out_dir,/directory) eq 0 then file_mkdir, out_dir
    data_file = file
    if file_test(data_file) eq 0 then begin
        ginfo = dictionary($
            'TITLE', 'MLT, B model SM, tilt-tilt_model for dipolarization study', $
            'TEXT', 'Generated by Sheng Tian at the University of Minnesota' )
        cdf_id = cdf_create(data_file)
        cdf_save_setting, ginfo, filename=cdf_id
    endif else cdf_id = cdf_open(data_file)


;---Common time.
    orbit_time_var = 'orbit_ut'
    if ~cdf_has_var(orbit_time_var, filename=cdf_id) then begin
        orbit_time_step = 60.
        orbit_times = make_bins(time_range+[0,-1]*orbit_time_step, orbit_time_step)
        cdf_save_var, orbit_time_var, value=orbit_times, filename=cdf_id
        setting = dictionary($
            'description', 'ut sec of orbit data', $
            'unit', 'sec', $
            'time_step', orbit_time_step )
        cdf_save_setting, setting, varname=orbit_time_var, filename=cdf_id
    endif

    bfield_time_var = 'bfield_ut'
    if ~cdf_has_var(bfield_time_var, filename=cdf_id) then begin
        bfield_time_step = 10.
        bfield_times = make_bins(time_range+[0,-1]*bfield_time_step, bfield_time_step)
        cdf_save_var, bfield_time_var, value=bfield_times, filename=cdf_id
        setting = dictionary($
            'description', 'ut sec of bfield data', $
            'unit', 'sec', $
            'time_step', bfield_time_step )
        cdf_save_setting, setting, varname=bfield_time_var, filename=cdf_id
    endif

    dummy_time_var = 'dummy_ut'
    if ~cdf_has_var(dummy_time_var, filename=cdf_id) then begin
        dummy_time_step = 10.
        cdf_save_var, dummy_time_var, value=time_range+[0,-1]*dummy_time_step, filename=cdf_id
        setting = dictionary($
            'description', 'ut sec of the time range', $
            'unit', 'sec' )
        cdf_save_setting, setting, varname=dummy_time_var, filename=cdf_id
    endif


;---Preparation.
    have_data = 1
    ndummy_time = 2
    ndim = 3
    b_sm_var = prefix+'b_sm'
    r_sm_var = prefix+'r_sm'
    azim_dp_read_bfield, time_range, probe=probe, errmsg=errmsg
    azim_dp_read_orbit, time_range, probe=probe, errmsg=errmsg
    get_data, b_sm_var, times
    have_data = n_elements(times) gt ndummy_time


;---MLT.
    the_var = prefix+'mlt'
    if have_data then begin
        get_data, r_sm_var, times, r_sm
        r_mag = cotran(r_sm, times, 'sm2mag')
        mlons = atan(r_mag[*,1],r_mag[*,0])*constant('deg')
        value = mlon2mlt(mlons, times)
        time_var = orbit_time_var
    endif else begin
        value = fltarr(ndummy_time)+!values.f_nan
        time_var = dummy_time_var
    endelse
    setting = dictionary($
        'have_data', have_data, $
        'display_type', 'scalar', $
        'short_name', 'MLT', $
        'unit', 'hr', $
        'depend_0', time_var )

    cdf_save_var, the_var, value=float(value), filename=cdf_id
    cdf_save_setting, setting, varname=the_var, filename=cdf_id


;---B model.
    the_var = prefix+'bmod_sm'
    model = 't89'
    if have_data then begin
        get_data, r_sm_var, times, r_sm
        r_gsm = cotran(r_sm, times, 'sm2gsm')
        ntime = n_elements(times)
        par = 2.
        bmod_gsm = fltarr(ntime,ndim)
        for ii=0, ntime-1 do begin
            tilt = geopack_recalc(times[ii])
            ; in-situ position
            rx = r_gsm[ii,0]
            ry = r_gsm[ii,1]
            rz = r_gsm[ii,2]
            ; in-situ B field.
            geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            geopack_t89, par, rx,ry,rz, dbx,dby,dbz
            bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor
        value = cotran(bmod_gsm, times, 'gsm2sm')
        time_var = orbit_time_var
    endif else begin
        value = fltarr(ndummy_time,ndim)+!values.f_nan
        time_var = dummy_time_var
    endelse
    setting = dictionary($
        'have_data', have_data, $
        'display_type', 'vector', $
        'short_name', strupcase(model)+' B', $
        'model', strupcase(model), $
        'unit', 'nT', $
        'coord', 'SM', $
        'coord_labels', constant('xyz'), $
        'depend_0', time_var )

    cdf_save_var, the_var, value=float(value), filename=cdf_id
    cdf_save_setting, setting, varname=the_var, filename=cdf_id



;---dtilt.
    bmod_sm_var = prefix+'bmod_sm'
    if have_data then begin
        cdf_load_var, bmod_sm_var, filename=cdf_id
        times = cdf_read_var(bfield_time_var, filename=cdf_id)
        interp_time, bmod_sm_var, times
        bmod_sm = get_var_data(bmod_sm_var)
        bmod_tilt = azim_df_calc_tilt(bmod_sm)
        b_sm_var = prefix+'b_sm'
        b_sm = get_var_data(b_sm_var)
        b_tilt = azim_df_calc_tilt(b_sm)

        the_var = prefix+'dtilt'
        value = b_tilt-bmod_tilt
        time_var = bfield_time_var
    endif else begin
        value = fltarr(ndummy_time)+!values.f_nan
        time_var = dummy_time_var
    endelse
    setting = dictionary($
        'have_data', have_data, $
        'display_type', 'scalar', $
        'short_name', tex2str('alpha')+'-'+tex2str('alpha')+'!D'+strupcase(model), $
        'unit', 'deg', $
        'depend_0', time_var )

    cdf_save_var, the_var, value=float(value), filename=cdf_id
    cdf_save_setting, setting, varname=the_var, filename=cdf_id



;---Clean up.
    cdf_close, cdf_id

end


time = time_double('2019-08-01')
file = join_path([homedir(),'test.cdf'])
file_delete, file, /allow_nonexistent
foreach probe, ['g13'] do $
    azim_dp_read_level2_data_gen_file, time, probe=probe, filename=file, test=1
end
