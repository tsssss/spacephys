;+
; Test to trace injeciton backward in time.
;-

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


;---Settings.
    time_range = time_double(['2008-01-19/05:00','2008-01-19/10:00'])
    probes = ['a','d','e']
    energy_range = [50d,250]; keV.
    test_times = time_double('2008-01-19/'+['07:17:46','07:16:00','07:15:23'])  ; 204 keV.
    test_times = time_double('2008-01-19/'+['07:19:22','07:16:00','07:15:23'])  ; 93 keV.
    target_time = time_double('2008-01-19/07:10')
    
    
    test_times = time_double('2008-01-19/'+['06:29:43','06:27:57','06:27:22'])  ; 204 keV.
    target_time = time_double('2008-01-19/06:26')    


;---Load data.
    foreach probe, probes do begin
        r_var = themis_read_orbit(time_range, probe=probe, coord='gsm', get_name=1)
        if check_if_update(r_var, time_range) then r_var = themis_read_orbit(time_range, probe=probe, coord='gsm')
        flux_var = themis_read_kev_flux(time_range, probe=probe, id='e', energy_range=energy_range, no_spec=1, get_name=1)
        if check_if_update(flux_var, time_range) then flux_var = themis_read_kev_flux(time_range, probe=probe, id='e', energy_range=energy_range, no_spec=1)
    endforeach


;---Trace back per probe.
    time_step = 9d
    psyms = [6,1,7]

    re = constant('re')
    re1 = 1/re
    v_drift_coef = 1e6*re1  ; make v in km/s.
    charge_state = 1   ; electron.

    sgopen, 0, xsize=5, ysize=5
    xrange = [target_time,max(test_times)]-target_time+[-1,1]*60
    yrange = [-5,3]
    plot, xrange, yrange, xstyle=1, ystyle=1, nodata=1, noerase=1, xrange=xrange, yrange=yrange

    foreach probe, probes, probe_id do begin
        prefix = 'th'+probe+'_'
        flux_var = prefix+'e_flux'
        get_data, flux_var, times, fluxs, energys
        r_gsm_var = prefix+'r_gsm'
        psym = psyms[probe_id]

        test_time = test_times[probe_id]
energys = max(energys)
energys = 204d
        foreach energy, energys do begin
            max_backward_time = test_time-target_time
            backward_times = test_time-smkarthm(0,max_backward_time, time_step, 'dx')

            r_gsm = get_var_data(r_gsm_var, at=test_time)
            foreach current_time, backward_times do begin    
                r_sm = cotran(r_gsm, current_time, 'gsm2sm')
                mlt = pseudo_mlt(r_sm)
                plots, current_time-target_time, mlt, psym=psym, symsize=0.2          

                ; Trace to the equator.
                r_eq = trace_to_equator(r_gsm, current_time)

                ; Get B and gradB.
                get_bfield_info, r_eq, current_time, b_gsm, gradb_gsm, model=model
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
            stop
        endforeach
        

    endforeach
end