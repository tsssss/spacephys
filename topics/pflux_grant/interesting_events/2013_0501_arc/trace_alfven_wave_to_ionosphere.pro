
function model_ne, z

    n0 = 6e4    ; cc.
    n1 = 1.34e7 ; cc.
    z0 = 318.   ; km.
    h = 383.    ; km.

    re = constant('re')
    zz = z*re
    t = (z0-zz)/h
    n = n0*exp(t)+n1*zz^(-1.55)
    return, n

end


function trace_alfven_wave_to_ionosphere, time, model=model, mod_time=model_time, ion_mass=ion_mass, $
    r_gsm=r_gsm, igrf=igrf, stop_dis=stop_dis, _extra=ex

    tmp = geopack_resolve_model(model)
    t89 = tmp.t89
    t96 = tmp.t96
    t01 = tmp.t01
    ts04 = tmp.ts04
    storm = tmp.storm
    
    if n_elements(model_time) eq -1 then model_time = time
    
    trace_result = dictionary()
    foreach trace_dir, [-1,1] do begin
        msg = (trace_dir eq -1)? 'north': 'south'
        
        if n_elements(par) eq 0 then begin
            time_range = model_time+[-1,1]*600
            if keyword_set(t89) then t89_par = 2
            par_var = geopack_read_par(time_range, model=model, t89_par=t89_par, _extra=ex)
            par = get_var_data(par_var, at=model_time)
        endif


        ps = geopack_recalc(time)
        xp = r_gsm[0]
        yp = r_gsm[1]
        zp = r_gsm[2]

        geopack_trace, xp,yp,zp, trace_dir, par, $
            xf,yf,zf, r0=r0, refine=1, ionosphere=1, $
            t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf, $
            fline=fline

        ntrace = n_elements(fline[*,0])
        ndim = 3
        bline = fltarr(ntrace,ndim)
        for ii=0,ntrace-1 do begin
            rx = fline[ii,0]
            ry = fline[ii,1]
            rz = fline[ii,2]

            if keyword_set(igrf) then begin
                geopack_igrf_gsm, rx,ry,rz, bx,by,bz
            endif else begin
                geopack_dip, rx,ry,rz, bx,by,bz
            endelse

            if model eq 'dip' or model eq 'dipole' or model eq 'igrf' then begin
                dbx = -1
                dby = -1
                dbz = -1
            endif else begin
                routine = 'geopack_'+model
                if model eq 't04s' then routine = 'geopack_ts04'
                call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
            endelse

            bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor


        bmag = snorm(bline)
        if n_elements(ion_mass) eq 0 then ion_mass = 1


        re = constant('re')
        fline_dis = snorm(fline)

        index = where(fline_dis gt stop_dis, nstep)
        v_para = fltarr(nstep)
        b_stop = sinterpol(bmag, fline_dis, stop_dis)
        bmag_step = [bmag[0:nstep-1],b_stop]
        bmag_step = (bmag_step[0:nstep-1]+bmag_step[1:nstep])*0.5
        r_stop = sinterpol(fline, fline_dis, stop_dis)
        r_step = [fline[0:nstep-1,*],r_stop]
        dr_step = snorm(r_step[0:nstep-1,*]-r_step[1:nstep,*])

        n_step = model_ne(snorm(r_step)-1)
        va0 = 1e-9/sqrt(1e6*!dpi*4e-7*1.67e-27)*1e-3    ; km/s, B in nT, n in cc, m in m_p.
        va_step = va0*bmag_step/sqrt(n_step*ion_mass)

        dt_wave = total(dr_step*re/va_step)

        trace_result[msg] = dictionary($
            'va_dtime', dt_wave, $
            'stop_dis', stop_dis )
    endforeach

    return, trace_result
end