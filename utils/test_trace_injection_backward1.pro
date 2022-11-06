;+
; Trace injection backward in time.
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

function trace_injection_backward, in_var, backward_times=backward_times, r_var=r_var, $
    model=model, target_time_range=target_time_range
    
    if n_elements(model) eq 0 then model = 't89'
    model = 't89'

;---Make a copy of the input flux and boost cadence.
    get_data, in_var, times, fluxs, energys, limits=lim
    
    if n_elements(target_time_range) eq 0 then target_time_range = minmax(times)
    
    out_var = in_var+'_backward'
    store_data, out_var, times, fluxs, energys, limits=lim
    nenergy = n_elements(energys)
    
    if n_elements(backward_times) eq 0 then begin
        max_backward_time = 1200d ; sec.
        time_step = 3d
        backward_times = smkarthm(time_step,max_backward_time, time_step, 'dx')
    endif
    nbackward_time = n_elements(backward_times)
    ndim = 3
    re = constant('re')
    re1 = 1d/re
    v_drift_coef = 1e6*re1  ; make v in km/s.
    charge_state = 1   ; electron.
    foreach time, times, time_id do begin
        mlt_list = list()
        dis_list = list()
        mlat_list = list()
        r_sm_list = list()
        r_gsm0 = get_var_data(r_var, at=time)
        foreach energy, energys, energy_id do begin
            flux0 = fluxs[time_id,energy_id]

            r_backs = fltarr(nbackward_time,ndim)
            r_eqs = fltarr(nbackward_time,ndim)
            v_drifts = fltarr(nbackward_time,ndim)
            b_gsms = fltarr(nbackward_time,ndim)
            foreach backward_time, backward_times, backward_id do begin
                time_step = backward_time
                if backward_id ne 0 then time_step -= backward_times[backward_id-1]
                current_time = time-backward_time+time_step
                if backward_id eq 0 then r_gsm = r_gsm0 else r_gsm = r_backs[backward_id-1,*]
                
                ; Trace to the equator.
                r_eq = trace_to_equator(r_gsm, current_time)
                r_eqs[backward_id,*] = r_eq

                ; Get B and gradB.
                get_bfield_info, r_eq, current_time, b_gsm, gradb_gsm, model=model
                b_gsms[backward_id,*] = b_gsm
                bmag = snorm(b_gsm)

                ; assume all energy is perp.
                rr = 1
                en_perp = energy*rr     ; in keV.
                en_para = energy*(1-rr)
                
                ; Get the drift velocity in gsm.
                v_drift = (en_perp+2*en_para)*vec_cross(b_gsm, gradb_gsm)/bmag^3*v_drift_coef/charge_state
                v_drifts[backward_id,*] = v_drift
;
;                ; Drift velocity using 3D tracing.
;                v_drift2 = calc_drift_velocity(rgsm=r_gsm, time=current_time, energy=energy, pitch_angle=90, species='e', model='igrf', bounce_period=bounce_period)
;                stop

                ; Get the previous location.
                r_gsm_new = r_gsm+v_drift*time_step*re1
                r_backs[backward_id,*] = r_gsm_new
            endforeach
            
            r_sm_backs = cotran(r_backs, time-backward_times, 'gsm2sm')
            mlt_backs = pseudo_mlt(r_sm_backs)
            mlt_list.add, mlt_backs
            dis_backs = snorm(r_sm_backs)
            dis_list.add, dis_backs
            mlat_backs = atan(r_sm_backs[*,2]/dis_backs)*constant('deg')
            mlat_list.add, mlat_backs
            r_sm_list.add, r_sm_backs
        endforeach
        
        xrange = []
        yrange = []
        zrange = []
        foreach data, r_sm_list, data_id do begin
            xrange = [xrange, minmax(data[*,0])]
            yrange = [yrange, minmax(data[*,1])]
            zrange = [zrange, minmax(data[*,2])]
        endforeach
        xrange = minmax(xrange)
        yrange = minmax(yrange)
        zrange = minmax(zrange)
        
        stop
    endforeach

stop

end


;---Settings.
    time_range = time_double(['2008-01-19/05:00','2008-01-19/10:00'])
    probes = ['a','d','e']
    energy_range = [50d,250]; keV.
    test_times = time_double('2008-01-19/'+['07:14','07:10','07:10'])
    test_duration = 20*60d  ; sec.

;---Load data.
    foreach probe, probes do begin
        r_var = themis_read_orbit(time_range, probe=probe, coord='gsm', get_name=1)
        if check_if_update(r_var, time_range) then r_var = themis_read_orbit(time_range, probe=probe, coord='gsm')
        flux_var = themis_read_kev_flux(time_range, probe=probe, id='e', energy_range=energy_range, no_spec=1, get_name=1)
        if check_if_update(flux_var, time_range) then flux_var = themis_read_kev_flux(time_range, probe=probe, id='e', energy_range=energy_range, no_spec=1, get_name=1)
    endforeach


backward_times = make_bins([3,600],3)
r_sm_var = themis_read_orbit(time_range, probe=probe, coord='sm')
r_sm = get_var_data(r_sm_var, times=uts)
mlt = pseudo_mlt(r_sm)
get_data, in_var, times, flux, energy, limits=lim
mlt = interpol(mlt, uts, times)
store_data, in_var+'_mlt', mlt, flux, energy, limits=lim
nenergy = n_elements(energy)
for ii=0,nenergy-1 do begin
    tmp = flux[*,ii]
    index = where(finite(tmp))
    tmp = interpol(tmp[index],times[index], times)
    flux[*,ii] = deriv(tmp)
endfor
lim.ylog = 0
store_data, in_var+'_deriv', times, flux, energy, limits=lim
stop
out_var = trace_injection_backward(in_var, backward_times=backward_times, r_var=r_var, model=model)

end