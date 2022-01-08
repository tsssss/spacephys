;+
; Test the diff b/w two versions of B model:
;   1. Low-res R GSM then inteprolate B to high-res
;   2. Interpolate low-res R GSM to high-res
;-

;---Input.
    time_range = time_double(['2013-01-01','2013-01-02'])
    probe = 'b'

;---Settings.
    time_step = 1/16.
    prefix = 'rbsp'+probe+'_'
    model = 't89'
    par = 2.
    common_times = make_bins(time_range, time_step)
    xyz = constant('xyz')

;---Load orbit.
    rbsp_read_orbit, time_range, probe=probe
    r_gse = get_var_data(prefix+'r_gse', times=orbit_times)
    r_gsm = cotran(r_gse, orbit_times, 'gse2gsm')

;---Version 1.
    bmod_gsm = float(r_gsm)
    foreach time, orbit_times, ii do begin
        tilt = geopack_recalc(time)
        rx = r_gsm[ii,0]
        ry = r_gsm[ii,1]
        rz = r_gsm[ii,2]
        ; in-situ B field.
        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        geopack_t89, par, rx,ry,rz, dbx,dby,dbz
        bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
    endforeach
    bmod_gsm = sinterpol(bmod_gsm, orbit_times, common_times, /quadratic)
    store_data, prefix+'bmod_gsm1', common_times, bmod_gsm
    add_setting, prefix+'bmod_gsm1', /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'T89 B', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )

;---Version 2.
    r_gsm = sinterpol(r_gsm, orbit_times, common_times, /quadratic)
    bmod_gsm = float(r_gsm)
    foreach time, common_times, ii do begin
        tilt = geopack_recalc(time)
        rx = r_gsm[ii,0]
        ry = r_gsm[ii,1]
        rz = r_gsm[ii,2]
        ; in-situ B field.
        geopack_igrf_gsm, rx,ry,rz, bx,by,bz
        geopack_t89, par, rx,ry,rz, dbx,dby,dbz
        bmod_gsm[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
    endforeach
    store_data, prefix+'bmod_gsm2', common_times, bmod_gsm
    add_setting, prefix+'bmod_gsm2', /smart, dictionary($
        'display_type', 'vector', $
        'short_name', 'T89 B', $
        'unit', 'nT', $
        'coord', 'GSM', $
        'coord_labels', xyz )

;---Combine.
    bmod_gsm1 = get_var_data(prefix+'bmod_gsm1')
    bmod_gsm2 = get_var_data(prefix+'bmod_gsm2')
    ndim = 3
    for ii=0,ndim-1 do begin
        the_var = prefix+'bmod_'+xyz[ii]
        store_data, the_var, common_times, [[bmod_gsm1[*,ii]],[bmod_gsm2[*,ii]]]
        add_setting, the_var, /smart, dictionary($
            'ytitle', 'nT', $
            'colors', sgcolor(['red','blue']), $
            'labels', ['1','2'] )
    endfor

    vars = prefix+'bmod_'+xyz
    tplot, vars, trange=time_range
end
