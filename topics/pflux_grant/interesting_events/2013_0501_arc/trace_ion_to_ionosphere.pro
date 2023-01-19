
function trace_ion_to_ionosphere, time, species=species, model=model, mod_time=model_time, $
    energy=energy, pitch_angles=pitch_angles, r_gsm=r_gsm, igrf=igrf, $
    beam_dis=beam_dis, _extra=ex
    
    if n_elements(pitch_angles) eq -1 then pitch_angles = [175d,162]
    if n_elements(beam_dis) eq -1 then beam_dis = 2.1


    rad = constant('rad')
    deg = constant('deg')
    
    
    tmp = geopack_resolve_model(model)
    t89 = tmp.t89
    t95 = tmp.t96
    t00 = tmp.t01
    ts03 = tmp.ts04
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
                if model eq 't03s' then routine = 'geopack_ts04'
                call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
            endelse

            bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor


        bmag = snorm(bline)
        case species of
            'e': mass0 = 1d/1836
            'p': mass0 = 1d
            'o': mass0 = 16d
            'he': mass0 = 4d
        endcase
        mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.


    ;---Conics.
        re = constant('re')
        conic_times = []
        conic_dis = []
        fline_dis = snorm(fline)

;        if n_elements(pitch_angles) eq -1 then pitch_angles = [175d,162]
        foreach pitch_angle, pitch_angles do begin
            pa = pitch_angle*rad
            mu = sin(pa)^2/bmag[0]
            b_mirror = 1/mu
            r_mirror = sinterpol(fline, bmag, b_mirror)

            dis_mirror = (snorm(r_mirror))[-1]
            index = where(fline_dis gt dis_mirror, nstep)
            v_para = fltarr(nstep)
            bmag_step = [bmag[0:nstep-1],b_mirror]
            bmag_step = (bmag_step[0:nstep-1]+bmag_step[1:nstep])*0.5
            r_step = [fline[0:nstep-1,*],r_mirror]
            dr_step = snorm(r_step[0:nstep-1,*]-r_step[1:nstep,*])

            vmag = sqrt(2*energy/mass0)*1e-3
            pitch_angle_step = asin(sqrt(mu*bmag_step))*deg
            vpara_step = vmag*cos(pitch_angle_step*rad)

            dt_step = dr_step*re/vpara_step
            dt_conic = total(dt_step)

            conic_times = [conic_times,dt_conic]
            conic_dis = [conic_dis,dis_mirror]
        endforeach


    ;---Beam.
        index = where(fline_dis gt beam_dis, nstep)
        v_para = fltarr(nstep)
        b_mirror = sinterpol(bmag, fline_dis, beam_dis)
        bmag_step = [bmag[0:nstep-1],b_mirror]
        bmag_step = (bmag_step[0:nstep-1]+bmag_step[1:nstep])*0.5
        r_mirror = sinterpol(fline, fline_dis, beam_dis)
        r_step = [fline[0:nstep-1,*],r_mirror]
        dr_step = snorm(r_step[0:nstep-1,*]-r_step[1:nstep,*])
        dt_beam = total(dr_step*re/vmag)
        beam_time = dt_beam

        trace_result[msg] = dictionary($
            'conic_dtime', conic_times, $
            'conic_dis', conic_dis, $
            'beam_dtime', beam_time, $
            'beam_dis', beam_dis )
    endforeach
    
    
    info = trace_result['north']
    info['conic_time'] = time-info['conic_dtime']
    info['beam_time'] = time-info['beam_dtime']
    
    info1 = trace_result['south']
    info1['conic_time'] = info['conic_time']+info1['conic_dtime']
    info1['beam_time'] = info['beam_time']+info1['beam_dtime']
    
    return, trace_result
    
end
