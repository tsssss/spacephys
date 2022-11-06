function trace_to_equator, r_gsm, current_time

    ps = geopack_recalc(current_time)
    par = 2

    rx = r_gsm[0]
    ry = r_gsm[1]
    rz = r_gsm[2]
    geopack_trace, rx,ry,rz, 1, par, xf,yf,zf, fline=fline1, igrf=1 ; to S-hem.
    geopack_trace, rx,ry,rz,-1, par, xf,yf,zf, fline=fline2, igrf=1 ; to N-hem.

    fline = [reverse(fline2,1),fline1[1:*,*]]
    tmp = max(snorm(fline), index)
    return, fline[index,*]

end


pro get_bfield_info, r_gsm, current_time, b_gsm, gradb_gsm, model=model

    ps = geopack_recalc(current_time)
    par = 2
    rx = r_gsm[0]
    ry = r_gsm[1]
    rz = r_gsm[2]
    geopack_igrf_gsm, rx,ry,rz, bx,by,bz
    geopack_t89, par, rx,ry,rz, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    b_gsm = [bx,by,bz]+[dbx,dby,dbz]

    dis = snorm(r_gsm)
    dr = 200/constant('re')*dis
    geopack_igrf_gsm, rx+dr,ry,rz, bx,by,bz
    geopack_t89, par, rx+dr,ry,rz, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    bp = [bx,by,bz]+[dbx,dby,dbz]
    geopack_igrf_gsm, rx-dr,ry,rz, bx,by,bz
    geopack_t89, par, rx-dr,ry,rz, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    bm = [bx,by,bz]+[dbx,dby,dbz]
    gradb_x = (snorm(bp)-snorm(bm))/(2*dr)

    geopack_igrf_gsm, rx,ry+dr,rz, bx,by,bz
    geopack_t89, par, rx,ry+dr,rz, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    bp = [bx,by,bz]+[dbx,dby,dbz]
    geopack_igrf_gsm, rx,ry-dr,rz, bx,by,bz
    geopack_t89, par, rx,ry-dr,rz, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    bm = [bx,by,bz]+[dbx,dby,dbz]
    gradb_y = (snorm(bp)-snorm(bm))/(2*dr)

    geopack_igrf_gsm, rx,ry,rz-dr, bx,by,bz
    geopack_t89, par, rx,ry,rz-dr, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    bp = [bx,by,bz]+[dbx,dby,dbz]
    geopack_igrf_gsm, rx,ry,rz-dr, bx,by,bz
    geopack_t89, par, rx,ry,rz-dr, dbx,dby,dbz
    dbx = 0 & dby = 0 & dbz = 0
    bm = [bx,by,bz]+[dbx,dby,dbz]
    gradb_z = (snorm(bp)-snorm(bm))/(2*dr)

    gradb_gsm = [gradb_x,gradb_y,gradb_z]

end


function trace_injection_backward, r_gsm0, backward_times, energy, ion=ion, time_step=time_step

    re = constant('re')
    re1 = 1/re
    v_drift_coef = 1e6*re1  ; convert v to km/s.
    charge_state = 1   ; electron.
    if keyword_set(ion) then charge_state = -1
    if n_elements(time_step) eq 0 then time_step = abs(total(backward_times[0:1]*[-1,1]))

    nbackward_time = n_elements(backward_times)
    ndim = 3
    r_gsm_back = fltarr(nbackward_time,ndim)
    r_gsm = r_gsm0
    foreach current_time, backward_times, time_id do begin
        ; Get MLT.
        r_sm = cotran(r_gsm, current_time, 'gsm2sm')
        r_gsm_back[time_id,*] = r_gsm

        ; Trace to the equator.
        r_eq = trace_to_equator(r_gsm, current_time)

        ; Get B and gradB.
        get_bfield_info, r_eq, current_time, b_gsm, gradb_gsm, model='t89'
        bmag = snorm(b_gsm)

        ; assume all energy is perp.
        rr = 1
        en_perp = energy*rr     ; in keV.
        en_para = energy*(1-rr)

        ; Get the drift velocity in gsm.
        v_drift = (en_perp+2*en_para)*vec_cross(b_gsm, gradb_gsm)/bmag^3*v_drift_coef/charge_state
        v_sm = cotran(v_drift, current_time, 'gsm2sm')

        r_hat = sunitvec(r_sm)
        z_hat = [0,0,1]
        phi_hat = vec_cross(r_hat,z_hat)
        v_phi = vec_dot(v_sm, phi_hat)
        omega_phi = v_phi/snorm(r_sm)/re*constant('deg')

        ; Update position.
        r_gsm += v_drift*re1*time_step
    endforeach

    return, r_gsm_back

end