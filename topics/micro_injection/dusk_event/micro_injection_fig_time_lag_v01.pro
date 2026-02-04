;+
; Show time lag among MMS.
;-


function micro_injection_fig_time_lag_v01, input_event_id, test=test, get_name=get_name, update=update, errmsg=errmsg


    errmsg = ''
    retval = !null
    version = 'v01'
    project= micro_injection_load_project()
    project_id = project.id

    if n_elements(input_event_id) eq 2 then begin
        time_range = time_double(input_event_id)
        event_id = time_string(time_range[0],tformat='YYYY_MMDD_hh')
    endif else begin
        event_id = input_event_id
    endelse
    event = project_get_event(project, id=event_id)
    time_range = event.time_range
    if n_elements(event) eq 0 then message, 'Inconsistency ...'

    if n_elements(plot_dir) eq 0 then plot_dir = event.plot_dir
    base = project_id+'_fig_time_lag_'+event_id+'_'+version+'.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(get_name) then return, plot_file
    print, plot_file
    if keyword_set(update) then file_delete, plot_file, allow_nonexist=1
    if keyword_set(test) then begin
        plot_file = 0
    endif else begin
        if file_test(plot_file) eq 1 then begin
            print, plot_file+' exists, skip ...'
            return, plot_file
        endif
    endelse


;---Load data.
    default_coord = 'gsm'
    probes = string(findgen(4)+1,format='(I0)')
    nprobe = n_elements(probes)
    colors = sgcolor(['red','green','blue','purple'])
    labels = strupcase('mms'+probes)
    comps = constant('xyz')
    ncomp = n_elements(comps)
    
    foreach phys_quant, ['bfield','efield','orbit','ion_vel'] do begin
        vars = list()
        foreach probe, probes do begin
            vars.add, lets_read(phys_quant, time_range, source=['mms',probe], coord=default_coord)
        endforeach
        
        ; recombine according to component.
        vars = vars.toarray()
        for ii=0,ncomp-1 do begin
            var = 'mms_'+phys_quant+'_'+comps[ii]
            times = get_var_time(vars[0])
            ntime = n_elements(times)
            data = fltarr(ntime,nprobe)
            for jj=0,nprobe-1 do data[*,jj] = (get_var_data(vars[jj],at=times, limits=lim))[*,ii]
            store_data, var, times, data, limits=lim
            options, var, labels=labels, colors=colors
        endfor

        if phys_quant eq 'orbit' then begin
            for ii=0,ncomp-1 do begin
                var = 'mms_'+phys_quant+'_'+comps[ii]
                get_data, var, times, data

                var = 'mms_d'+phys_quant+'_'+comps[ii]
                del_data, var
                data0 = data[*,0]
                for jj=0,nprobe-1 do data[*,jj] -= data0
                ;data[*,0] = !values.f_nan
                data *= constant('re')
                store_data, var, times, data
                options, var, labels=labels, colors=colors, ytitle='(km)', labflag=-1
            endfor
        endif
    endforeach
    
    ; convert to fac.
    external_model = 't89'
    internal_model = 'igrf'
    b0_window = 1200d
    fac_coord = 'fac'
    fac_labels = ['b','w','o']
    foreach probe, probes, pid do begin
        prefix = 'mms'+probe+'_'
        b_var = prefix+'b_'+default_coord
        r_var = prefix+'r_'+default_coord
        bmod_var = lets_read_geopack_bfield(orbit_var=r_var, external_model=external_model, internal_model=internal_model)
        b_vars = lets_decompose_bfield(b_var=b_var, b0_window=b0_window, bmod_var=bmod_var)
        b0_var = b_vars['b0']
        q_fac_var = lets_define_fac(r_var=r_var,b_var=b0_var)
        foreach var, prefix+['u','b1','e']+'_'+default_coord do begin
            vec_coord = get_var_data(var, times=times, setting=setting)
            q_coord2fac = qslerp(get_var_data(q_fac_var, times=uts), uts, times)
            m_coord2fac = qtom(q_coord2fac)
            vec_fac = rotate_vector(vec_coord, m_coord2fac)
            var_info = streplace(var, default_coord, fac_coord)
            store_data, var_info, times, vec_fac
            setting['coord'] = fac_coord
            setting['coord_labels'] = fac_labels
            add_setting, var_info, smart=1, setting
        endforeach
    endforeach
    
    ; combine b1.
    phys_quant = 'b1_'+fac_coord
    vars = 'mms'+probes+'_b1_'+fac_coord
    for ii=0,ncomp-1 do begin
        var = 'mms_'+phys_quant+'_'+comps[ii]
        times = get_var_time(vars[0])
        ntime = n_elements(times)
        data = fltarr(ntime,nprobe)
        for jj=0,nprobe-1 do data[*,jj] = (get_var_data(vars[jj],at=times, limits=lim))[*,ii]
        store_data, var, times, data, limits=lim
        options, var, labels=labels, colors=colors
    endfor
    
    phys_quant = 'ion_vel_'+fac_coord
    vars = 'mms'+probes+'_u_'+fac_coord
    for ii=0,ncomp-1 do begin
        var = 'mms_'+phys_quant+'_'+comps[ii]
        times = get_var_time(vars[0])
        ntime = n_elements(times)
        data = fltarr(ntime,nprobe)
        for jj=0,nprobe-1 do data[*,jj] = (get_var_data(vars[jj],at=times, limits=lim))[*,ii]
        store_data, var, times, data, limits=lim
        options, var, labels=labels, colors=colors
    endfor
    
    phys_quant = 'e_'+fac_coord
    vars = 'mms'+probes+'_e_'+fac_coord
    for ii=0,ncomp-1 do begin
        var = 'mms_'+phys_quant+'_'+comps[ii]
        times = get_var_time(vars[0])
        ntime = n_elements(times)
        data = fltarr(ntime,nprobe)
        for jj=0,nprobe-1 do data[*,jj] = (get_var_data(vars[jj],at=times, limits=lim))[*,ii]
        store_data, var, times, data, limits=lim
        options, var, labels=labels, colors=colors
    endfor
    
    ; calc v_exb
    foreach probe, probes do begin
        prefix = 'mms'+probe+'_'
        e_var = prefix+'e_'+fac_coord
        b_var = prefix+'b0_gsm'
        v_var = prefix+'vexb_'+fac_coord
        e_vec = get_var_data(e_var, times=times)
        bmag = snorm(get_var_data(b_var, at=times))
        b_vec = fltarr(n_elements(times),ncomp)
        b_vec[*,0] = bmag
        v_vec = vec_cross(e_vec,b_vec)
        cc = 1d3/bmag^2
        for ii=0,ncomp-1 do v_vec[*,ii] *= cc
        store_data, v_var, times, v_vec
        add_setting, v_var, smart=1, dictionary($
            'display_type', 'vector', $
            'short_name', 'V', $
            'unit', 'km/s', $
            'coord', fac_coord, $
            'coord_labels', fac_labels )
    endforeach
    phys_quant = 'vexb_'+fac_coord
    vars = 'mms'+probes+'_vexb_'+fac_coord
    for ii=0,ncomp-1 do begin
        var = 'mms_'+phys_quant+'_'+comps[ii]
        times = get_var_time(vars[0])
        ntime = n_elements(times)
        data = fltarr(ntime,nprobe)
        for jj=0,nprobe-1 do data[*,jj] = (get_var_data(vars[jj],at=times, limits=lim))[*,ii]
        store_data, var, times, data, limits=lim
        options, var, labels=labels, colors=colors
    endfor

;---Make the plot.
    fig_size = [12,12]
    
    
    prefix = 'mms1_'
    plot_vars = prefix+['b','u','r']+'_'+default_coord
    plot_tr = time_range
    tickinterval = 10*60d

    
    zoom_vars = list()
    zoom_tr = mean(time_range)+[0,60]
    tickinterval = 6d
    foreach var, ['b1','e']+'_'+fac_coord do begin
        zoom_vars.add, 'mms_'+var+'_'+comps, extract=1
    endforeach
    zoom_vars.add, 'mms_ion_vel_fac_x'
    foreach ii, ['y','z'] do zoom_vars.add, 'mms_'+['vexb','ion_vel']+'_fac_'+ii, extract=1
    zoom_vars = zoom_vars.toarray()
    
    ;stop
    ;tmp = sgplot(zoom_vars, xrange=zoom_tr)
    tmp = sgplot(plot_vars, xrange=time_range, filename=1)
    
    stop
    

    ; Init plot_vars.
    plot_info = orderedhash()
    
    coord = default_coord

    ; B field.
    foreach comp, comps do begin
        plot_info['mms_bfield_'+comp] = dictionary($
            'routine', 'plot_line', $
            'panel_label_text', strupcase(coord+' B!D'+comp), $
            'ypan', 0.8, $
            'setting', dictionary() )
    endforeach

    ; Ion vel.
    foreach comp, comps do begin
        plot_info['mms_ion_vel_'+comp] = dictionary($
            'routine', 'plot_line', $
            'panel_label_text', strupcase(coord+' U!D'+comp), $
            'ypan', 0.8, $
            'setting', dictionary() )
    endforeach

    ; dorbit.
    foreach comp, comps do begin
        plot_info['mms_dorbit_'+comp] = dictionary($
            'routine', 'plot_line', $
            'panel_label_text', strupcase(coord+' dR!D'+comp), $
            'ypan', 0.8, $
            'setting', dictionary() )
    endforeach


    plot_vars = plot_info.keys()
    nplot_var = n_elements(plot_vars)
    ; Default settings.
    panel_letters = letters(nplot_var)
    foreach plot_var, plot_vars, pid do begin
        my_info = plot_info[plot_var]
        if ~my_info.haskey('ypan') then my_info['ypan'] = 1d
        if ~my_info.haskey('panel_label_text') then my_info['panel_label_text'] = ' '
        if ~my_info.haskey('panel_letter') then my_info['panel_letter'] = panel_letters[pid]
        if ~my_info.haskey('panel_label_msg') then my_info['panel_label_msg'] = my_info['panel_letter']+') '+my_info['panel_label_text']
    endforeach
    
    prefix = 'mms1_'
    var_labels = prefix+['mlat','dis','mlt']
    nvar_label = n_elements(var_labels)
    options, prefix+'mlat', ytitle='MLat (deg)'
    options, prefix+'dis', ytitle='|R| (Re)'
    options, prefix+'mlt', ytitle='MLT (h)'
    var = prefix+'mlt'
    get_data, var, times, data
    index = where(data le 0, count)
    if count ne 0 then begin
        data[index] += 24
        store_data, var, times, data
    endif
    
    margins = [12,3.5+nvar_label,8,2]
    ypans = []
    foreach plot_var, plot_vars, pid do begin
        ypans = [ypans,(plot_info[plot_var])['ypan']]
    endforeach
    plot_poss = panel_pos(plot_file, nypan=nplot_var, fig_size=fig_size, ypans=ypans, pansize=[12,0.8], margins=margins)
    
    ; Use positions to determine [x,y]ticklen.
    abs_ticklen = 0.3
    foreach plot_var, plot_vars, pid do begin
        my_info = plot_info[plot_var]
        my_info['position'] = plot_poss[*,pid]
        my_info['abs_ticklen'] = abs_ticklen
        ;plot_info[plot_var] = my_info
    endforeach
    
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    
    tplot_options, 'tickinterval', tickinterval
    foreach plot_var, plot_vars, pid do begin
        my_info = plot_info[plot_var]
    
        my_pos = my_info['position']
        plot_routine = my_info['routine']
        plot_setting = my_info['setting']
        plot_setting['position'] = my_pos
        plot_setting['noerase'] = (pid eq 0)? 0: 1
        plot_setting['xtickformat'] = (pid eq nplot_var-1)? '': '(A1)'
        plot_setting['novtitle'] = (pid eq nplot_var-1)? 0: 1
        plot_setting['time_range'] = plot_tr
        plot_setting['tickinterval'] = tickinterval
        plot_setting['panel_label_pos'] = [xchsz*1,my_pos[3]-ychsz*0.7]
        plot_setting['panel_label_msg'] = my_info['panel_label_msg']
        plot_setting['var_labels'] = var_labels
        plot_setting['vlab_margin'] = margins[0]-1
        plot_setting['trange'] = plot_tr
        foreach key, plot_setting.keys() do begin
            index = strpos(key,'tick_setting')
            if index[0] ne -1 then begin
                tick_setting = plot_setting[key]
                foreach comp, constant('xyz') do begin
                    the_comp = comp+'tickv'
                    if tick_setting.haskey(the_comp) then begin
                        tick_setting[comp+'ticks'] = n_elements(tick_setting[the_comp])-1
                    endif
                endforeach
            endif
        endforeach

        plot_setting = plot_setting.tostruct()        
        tmp = call_function(plot_routine, plot_var, _extra=plot_setting)
    endforeach


    if keyword_set(test) then stop
    sgclose
    

    return, plot_file



end


event_id = '2015_0901_11'
event_id = '2015_0901_18'
print, micro_injection_fig_time_lag_v01(event_id, test=1)
end