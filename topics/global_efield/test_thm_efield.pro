;+
; Read 1 day of data of ion_vel, bfield, efield.
; Coerce them to 1 min cadence. Calculate E_vxB = -VxB, then remove the spin-axis component.
; Then compare to the E_measured.
;-

    date = '2014-08-28'
    probe = 'd'
    
    ; Duan+2016.
    date = '2008-02-03'
    probe = 'd'
    date = '2008-02-03'
    probe = 'e'
    date = '2008-02-07'
    probe = 'd'
    
    secofday = 86400d
    time = time_double(date)+[0,secofday]
    themis_read_bfield, time, probe=probe
    themis_read_ion_vel, time, probe=probe
    themis_read_efield, time, probe=probe
    themis_read_orbit, time, probe=probe
    thm_load_state, trange=time, probe=probe, /get_support

    dt = 60.
    ntime = round((time[1]-time[0])/dt)
    times = smkarthm(time[0],dt,ntime, 'x0')

    pre0 = 'th'+probe+'_'
    vars = pre0+['r_gsm','u_gsm','b_gsm','edot0_gsm']
    foreach var, vars do begin
        get_data, var, uts, dat
;        if var eq pre0+'edot0_gsm' then begin
;            for ii=0, 2 do dat[*,ii] -= smooth(dat[*,ii],dt*50,/edge_truncate,/nan)
;        endif
        dims = size(dat,/dimensions)
        dims[0] = ntime
        data = make_array(dims, /float)
        for ii=0, ntime-1 do begin
            index = lazy_where(uts, '[)', times[ii]+[0,dt], count=count)
            if count eq 0 then continue
            vals = dat[0,*]
            foreach val, vals, ll do vals[ll] = mean(dat[index,ll])
            data[ii,*] = vals
        endfor
        get_data, var, limits=lim
        store_data, var+'_better', times, data, limits=lim
    endforeach


;---Calculate E = -uxB, the electric field due to plasma motion.
    u_gsm_var = pre0+'u_gsm_better'
    b_gsm_var = pre0+'b_gsm_better'
    get_data, u_gsm_var, times, ugsm
    get_data, b_gsm_var, times, bgsm
    euxb_gsm = -vec_cross(ugsm,bgsm)*1e-3   ; u in km/s, B in nT, e in mV/m.
    euxb_gsm_var = pre0+'euxb_gsm_better'
    store_data, euxb_gsm_var, times, euxb_gsm
    add_setting, euxb_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E!S!UUxB!N!R', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}

;---Calculate E = vxB, the electric field due to s/c motion.
    nvecdim = 3
    re = 6378d
    r_gsm_var = pre0+'r_gsm_better'
    get_data, r_gsm_var, times, rgsm
    vgsm = fltarr(ntime,nvecdim)
    for ii=0,nvecdim-1 do vgsm[*,ii] = deriv(rgsm[*,ii])/dt*re
    v_gsm_var = pre0+'v_gsm_better'
    store_data, v_gsm_var, times, vgsm
    add_setting, v_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'km/s', $
        short_name: 'V', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}
    v_gsm_var = pre0+'v_gsm_better'
    b_gsm_var = pre0+'b_gsm_better'
    get_data, v_gsm_var, times, vgsm
    get_data, b_gsm_var, times, bgsm
    evxb_gsm = vec_cross(vgsm,bgsm)*1e-3   ; u in km/s, B in nT, e in mV/m.
    evxb_gsm_var = pre0+'evxb_gsm_better'
    store_data, evxb_gsm_var, times, evxb_gsm
    add_setting, evxb_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E!S!UVxB!N!R', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}



;---Remove spin-axis for Edot0, E_uxb and E_vxb.
    foreach var, ['euxb','evxb','edot0'] do begin
        gsm_var = pre0+var+'_gsm_better'
        dsl_var = pre0+var+'_dsl_better'
        thm_cotrans, gsm_var, dsl_var, in_coord='gsm', out_coord='dsl'
        get_data, dsl_var, times, data
        data[*,2] = 0
        dsl0_var = pre0+var+'0_dsl_better'
        store_data, dsl0_var, times, data
        gsm0_var = pre0+var+'0_gsm_better'
        thm_cotrans, dsl0_var, gsm0_var, in_coord='dsl', out_coord='gsm'
        get_data, gsm_var, limits=lims
        store_data, gsm0_var, limits=lims
    endforeach

    e1_gsm_var = pre0+'e1_gsm_better'
    sys_subtract, pre0+'edot00_gsm_better', pre0+'evxb0_gsm_better', to=e1_gsm_var
    add_setting, e1_gsm_var, /smart, {$
        display_type: 'vector', $
        unit: 'mV/m', $
        short_name: 'E!S!UVxB!N!R', $
        coord: 'GSM', $
        coord_labels: ['x','y','z']}

end