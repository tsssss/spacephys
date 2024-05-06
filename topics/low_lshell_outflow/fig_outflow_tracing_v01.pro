function get_mass0, species

    case species of
        'e': mass0 = 1d/1836
        'p': mass0 = 1d
        'o': mass0 = 16d
        'he': mass0 = 4d
    endcase
    mass0 = mass0*(1.67e-27/1.6e-19)   ; E in eV, mass in kg.
    return, mass0
    
end


;+
; Use one characteristic model_time to calculate the mirror distance 
; and travel times from both hemispheres. The individual times for each 
; test particle will be applied separately. This is essentially neglecting 
; the changes in the B field and SC position.
; 
; Thus this tracing only depends on the model_time (and internal and 
; external models), the energy and species of the test particles. Actually,
; the energy can be applied seperately as well because the travel time is
; inversely proportional to the parallel velocity, and thus the sqrt of energy.
;-
function trace_ion_conics, input

;---Model settings.
    internal_model = input.internal_model
    external_model = input.external_model

    tmp = geopack_resolve_model(external_model)
    t89 = tmp.t89
    t96 = tmp.t96
    t01 = tmp.t01
    ts04 = tmp.ts04
    storm = tmp.storm
    igrf = (internal_model eq 'igrf')? 1: 0
    r0 = 100d/constant('re')+1

;---Get the start position, and initial and opposite hemispheres.
    model_time = input.model_time
    r_gsm_var = input.r_gsm_var
    r_gsm = transpose(get_var_data(r_gsm_var, at=model_time))
    r_mag = cotran_pro(r_gsm, model_time, coord_msg=['gsm','mag'])
    mlat = calc_mlat(r_mag)
    init_hem = (mlat ge 0)? 'north': 'south'
    opposite_hem = (init_hem eq 'north')? 'south': 'north'
    ps = geopack_recalc(model_time)
    xp = r_gsm[0]
    yp = r_gsm[1]
    zp = r_gsm[2]

;---Get model parameters, field line positions and B at these positions.
    model_par_var = external_model+'_par'
    par = get_var_data(model_par_var, at=model_time)

    fline_info = dictionary()
    foreach hem, [init_hem,opposite_hem] do begin
        trace_dir = (hem eq 'north')? -1: 1

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

            if external_model eq 'dip' or external_model eq 'dipole' or external_model eq 'igrf' then begin
                dbx = 0
                dby = 0
                dbz = 0
            endif else begin
                routine = 'geopack_'+external_model
                if external_model eq 't04s' then routine = 'geopack_ts04'
                call_procedure, routine, par, rx,ry,rz, dbx,dby,dbz
            endelse

            bline[ii,*] = [bx,by,bz]+[dbx,dby,dbz]
        endfor

        fline_info[hem] = dictionary($
            'fline', fline, $
            'bline', bline )
    endforeach


;---Trace conics to initial hemisphere.
    the_info = fline_info[init_hem]
    bline = the_info['bline']
    fline = the_info['fline']

    ; Locate the mirror point's altitude.
    all_bmags = snorm(bline)
    pitch_angle = input.pitch_angle*constant('rad')
    mu = sin(pitch_angle)^2/all_bmags[0]
    bmag_mirror = 1d/mu
    fline_mirror = sinterpol(fline, all_bmags, bmag_mirror)
    mirror_dis = snorm(fline_mirror)

;---Calculate the travel time for init_hem and opposite_hem.
    trace_info = dictionary($
        'species', species, $
        'model_time', model_time, $
        'internal_model', internal_model, $
        'external_model', external_model, $
        'mirror_dis', mirror_dis, $
        'mirror_bmag', bmag_mirror, $
        'init_hem', init_hem, $
        'opposite_hem', opposite_hem )

    travel_time_info = dictionary()
    foreach hem, [init_hem,opposite_hem] do begin
        the_info = fline_info[hem]
        bline = the_info['bline']
        fline = the_info['fline']

        all_bmags = snorm(bline)
        index = where(all_bmags lt bmag_mirror, count)
        if count eq 0 then message, 'Inconsistency ...'
        bmags = [all_bmags[index],bmag_mirror]
        bline = sinterpol(bline, all_bmags, bmags)
        fline = sinterpol(fline, all_bmags, bmags)

        ; Calculate the travel time.
        npoint = n_elements(bmags)
        bmag_steps = (bmags[0:npoint-2]+bmags[1:npoint-1])*0.5
        pa_steps = asin(sqrt(mu*bmag_steps))

        target_energy = input.energy
        species = input.species
        mass0 = get_mass0(species)
        vmags = sqrt(2*target_energy/mass0)*1e-3

        pitch_angle_steps = asin(sqrt(mu*bmag_steps))
        vpara_steps = vmags*cos(pitch_angle_steps)

        r_steps = fline[0:npoint-2,*]-fline[1:npoint-1,*]
        dr_steps = snorm(r_steps)*constant('re')
        dt_steps = dr_steps/vpara_steps
        travel_time = total(dt_steps)

        travel_time_info[hem] = travel_time
;        trace_info[hem] = the_info
    endforeach
    trace_info['travel_time'] = travel_time_info
    
    return, trace_info

end



function fig_outflow_tracing_v01, event_info=event_info, test=test

    version = 'v01'
    id = '2015_0317'
    if n_elements(event_info) eq 0 then event_info = low_lshell_outflow_load_data(id)
    label_size = 0.8
    symsize = 0.3
    psym = 1
;    time_range = time_double(['2015-03-17/06:00','2015-03-17/15:00'])
    time_range = time_double(['2015-03-17/12:20','2015-03-17/13:20'])

    probe = 'a'
    prefix = 'rbsp'+probe+'_'

    p_en_var = rbsp_read_en_spec(time_range, probe=probe, species='p')
    p_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='p', $
        energy_range=[50,5e4], update=1)
    o_en_var = rbsp_read_en_spec(time_range, probe=probe, species='o')
    o_pa_var = rbsp_read_pa_spec(time_range, probe=probe, species='o', $
        energy_range=[50,5e4], update=1)

    the_vars = [p_pa_var,o_pa_var]
    options, the_vars, 'ytitle', 'PA!C(deg)'
    
    suffix = ['','_para','_anti','_perp']
    the_vars = [p_en_var+suffix,o_en_var+suffix]
    options, the_vars, 'ytitle', 'Energy!C(eV)'
    options, the_vars, ytickname=['10!U2','10!U3','10!U4'], $
        ytickv=[100,1000,10000], yticks=2, yminor=9, $
        yrange=[20,3e4]

    the_vars = o_pa_var
    options, the_vars, 'ytitle', 'PA!C(deg)'

    zticklen = -0.5
    the_vars = [p_en_var,o_en_var,o_pa_var]
    options, the_vars, zticklen=zticklen, zminor=9
    
    
;---Determine test particles.
    if probe eq 'a' then begin
        target_time_range = time_double(['2015-03-17/07:00','2015-03-17/15:00'])
    endif else begin
        target_time_range = time_double(['2015-03-17/12:00','2015-03-17/20:00'])
        target_time_range = time_double(['2015-03-17/03:00','2015-03-17/11:00'])
    endelse
    target_energy_range = [20,3e4]
    spec_var = o_en_var
    get_data, spec_var, times, fluxs, energys
    time_index = where_pro(times, '[]', target_time_range)
    target_times = times[time_index]
    energy_index = where_pro(energys[0,*], '[]', target_energy_range)
    the_energys = energys[0,energy_index]

    ntarget = n_elements(target_times)
    target_energys = fltarr(ntarget)
    target_fluxs = fltarr(ntarget)
    foreach tid, time_index, target_id do begin
        the_fluxs = reform(fluxs[tid,energy_index])
        the_fluxs = smooth(the_fluxs,3, edge_zero=1)
        target_fluxs[target_id] = max(the_fluxs, the_energy_index)
        target_energys[target_id] = the_energys[the_energy_index]
    endforeach
    var = spec_var+'_test_particle'
    
    lshells = get_var_data(prefix+'lshell', at=target_times)
    index = where(target_fluxs le 1e6 or lshells le 3, count)
    if count ne 0 then begin
        target_times[index] = !values.f_nan
        target_energys[index] = !values.f_nan
    endif
    store_data, var, target_times, target_energys, target_fluxs
    
    ; manual adjustments.
    manual_adjustments = hash($
        '2015-03-17/12:39:18', 701.2, $
        '2015-03-17/12:39:41', 565.5, $
        '2015-03-17/12:40:04', 565.5, $
        '2015-03-17/12:40:49', 456.0, $
        '2015-03-17/12:53:41', 47.6 )
    foreach the_time, manual_adjustments.keys() do begin
        tmp = min(target_times-time_double(the_time), tid, abs=1)
        target_energys[tid] = manual_adjustments[the_time]
        tmp = min(energys[tid,*]-target_energys[tid], eid, abs=1)
        target_fluxs[tid] = fluxs[tid,eid]
    endforeach
    store_data, var, target_times, target_energys, target_fluxs
    
    
    vars = prefix+'p_en_spec_'+['para','anti']
    options, vars, 'zrange', [1e5,1e8]
    
    vars = prefix+'o_en_spec_'+['para','anti']
    options, vars, 'zrange', [1e5,1e8]
    
    
    plot_dir = event_info['plot_dir']
    id = event_info['id']
    plot_file = join_path([plot_dir,'fig_outflow_tracing_'+id+'_'+version+'.pdf'])
    if keyword_set(test) then plot_file = 0
    fig_size = [6,7]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz

    plot_vars = prefix+['o_en_spec_'+['para','anti'],'o_pa_spec',$
        'p_en_spec_'+['para','anti'],'p_pa_spec','e_spec']
    nplot_var = n_elements(plot_vars)
    fig_labels = letters(nplot_var)+') '+['O+ '+['para','anti','PA'],'H+ '+['para','anti','PA'],'E spec']
    margins = [10,5,7,1]
    ypans = fltarr(nplot_var)+1
    ypans[-1] = 1.5
    poss = sgcalcpos(nplot_var, margins=margins, ypans=ypans)
    tplot_options, 'tickinterval', 600
    options, prefix+'lshell','ytitle', 'L-shell (#)'
    options, prefix+'mlat', 'ytitle', 'MLat (deg)'
    options, plot_vars, 'zticklen', -0.5
    
    tplot, plot_vars, trange=time_range, position=poss, $
        var_label=prefix+['lshell','mlat'], vlab_margin=9
    for pid=0,nplot_var-1 do begin
        tpos = poss[*,pid]
        ;var = plot_vars[pid]
        ;set_axis, var, position=tpos, xrange=time_range
        tx = tpos[0]-xchsz*9
        ty = tpos[3]-ychsz*0.7
        msg = fig_labels[pid]
        xyouts, tx,ty,msg, normal=1, charsize=label_size
    endfor
    
    

;---Solve for initial condition.
    test_time_range = time_double(['2015-03-17/12:12:27','2015-03-17/12:34:30'])
    test_time_range = time_double(['2015-03-17/12:37:00','2015-03-17/12:55:34'])
    test_time_range = time_double(['2015-03-17/12:39:18','2015-03-17/12:55'])
    test_particle_var = prefix+'o_en_spec_test_particle'
    target_energys = get_var_data(test_particle_var, times=target_times, in=test_time_range)

    
    
    species = 'o'
    para_var = prefix+species+'_en_spec_para'
    anti_var = prefix+species+'_en_spec_anti'
    
    pid = where(plot_vars eq para_var, count)
    if count ne 0 then begin
        tpos = poss[*,pid]
        set_axis, para_var, position=tpos, xrange=time_range
        plots, target_times, target_energys, psym=1, symsize=symsize, color=sgcolor('purple')
    endif
    

    ; Solver input.
    internal_model = 'igrf'
    external_model = 't89'
    par_var = external_model+'_par'
    if check_if_update(par_var, time_range) then par_var = geopack_read_par(time_range, model=external_model)

    solver_input = dictionary($
        'internal_model', internal_model, $
        'external_model', external_model, $
        'model_time', mean(test_time_range), $
        'r_gsm_var', prefix+'r_gsm', $
        'species', species, $
        'energy', 1e3, $    ; eV.
        'pitch_angle', 15d )
    
    ; Trace for ion conics.
    trace_result = trace_ion_conics(solver_input)

    ; Trace settings.
    plot_trace_setting = dictionary($
        para_var, dictionary($
            'n_init_hem', 2, $
            'n_second_hem', 2), $
        anti_var, dictionary($
            'n_init_hem', 2, $
            'n_second_hem', 2) )
    
    test_particle_var = prefix+'o_en_spec_test_particle'
    target_energys = get_var_data(test_particle_var, times=target_times, in=test_time_range)
    init_hem = trace_result['init_hem']
    opposite_hem = trace_result['opposite_hem']
    
    show_initial_energization = 1
    if keyword_set(show_initial_energization) then begin
        pid = 0
        tpos = poss[*,pid]
        the_var = plot_vars[pid]
        set_axis, the_var, position=tpos, xrange=time_range

        the_color = sgcolor('red')
        init_times = []
        foreach target_time, target_times, tid do begin
            target_energy = target_energys[tid]
            init_time = target_time-trace_result.travel_time[init_hem]*sqrt(solver_input.energy/target_energy)
            init_times = [init_times, init_time]
        endforeach
        init_tr = minmax(init_times)
        timebar, init_tr, linestyle=2, color=the_color
        msg = 'Initial time (UT): '+strjoin(time_string(init_tr,tformat='hh:mm'),'-')
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*2
        xyouts, tx,ty,msg, normal=1, charsize=label_size, color=the_color
        msg = 'Initial altitude (Re): '+string(trace_result.mirror_dis-1,format='(F3.1)')
        tx = tpos[0]+xchsz*0.5
        ty = tpos[3]-ychsz*1
        xyouts, tx,ty,msg, normal=1, charsize=label_size, color=the_color
        
        tx = tpos[2]-xchsz*0.5
        ty = tpos[3]-ychsz*1
        msg = strupcase('rbsp-'+probe)
        xyouts, tx,ty,msg, normal=1, alignment=1, color=sgcolor('white')
    endif

    foreach target_time, target_times, tid do begin
        target_energy = target_energys[tid]
        traced_travel_time = dictionary()
        foreach hem, [init_hem,opposite_hem] do begin
            traced_travel_time[hem] = trace_result.travel_time[hem]*sqrt(solver_input.energy/target_energy)
        endforeach
        init_time = target_time-traced_travel_time[init_hem]
        
        if target_energy le 50 then continue
        
        show_initial_energization = 1

        the_var = anti_var
        pid = where(plot_vars eq the_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            set_axis, the_var, position=tpos, xrange=time_range

;            the_time = init_time
;            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('red')

            the_time = init_time+traced_travel_time[opposite_hem]
            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('orange')            
        endif
        
        
        the_var = para_var
        pid = where(plot_vars eq the_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            set_axis, the_var, position=tpos, xrange=time_range

            the_time = init_time+traced_travel_time[init_hem]
            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('red')
            
;            the_time = init_time+traced_travel_time[opposite_hem]+traced_travel_time[init_hem]*2
;            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('red')
        endif
        
        
        ;if target_energy ge 1e3/4 then continue
        
        ; Test for H+
        foreach hem, [init_hem,opposite_hem] do begin
            traced_travel_time[hem] /= 4
        endforeach
        
        the_var = prefix+'p_en_spec_anti'
        pid = where(plot_vars eq the_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            set_axis, the_var, position=tpos, xrange=time_range

;            the_time = init_time
;            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('red')

            the_time = init_time+traced_travel_time[opposite_hem]
            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('orange')
        endif
        
        the_var = prefix+'p_en_spec_para'
        pid = where(plot_vars eq the_var, count)
        if count ne 0 then begin
            tpos = poss[*,pid]
            set_axis, the_var, position=tpos, xrange=time_range

;            the_time = init_time
;            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('red')

            the_time = init_time+traced_travel_time[init_hem]
            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('orange')

;            the_time = init_time+traced_travel_time[opposite_hem]+traced_travel_time[init_hem]*2
;            plots, the_time, target_energy, psym=psym, symsize=symsize, color=sgcolor('red')
        endif

    endforeach
    
    
    if keyword_set(test) then stop
    sgclose
    
    return, plot_file

end


print, fig_outflow_tracing_v01(event_info=event_info, test=0)
end