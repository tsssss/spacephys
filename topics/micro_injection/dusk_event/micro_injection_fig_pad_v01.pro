

function micro_injection_fig_pad_v01, input_event_id, probe=probe, $
    plot_dir=plot_dir, test=test, get_name=get_name, update=update

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
    base = project_id+'_fig_pad_'+event_id+'_mms'+probe+'_'+version+'.pdf'
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
    prefix = 'mms'+probe+'_'
    default_coord = 'gsm'
    mission_probe = 'mms'+probe

    ; Particle related vars.
    ele_pad_var = lets_read_this(func='mms_read_pad_ele_all', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load e thermal or kev data ...'
        return, retval
    endif

    ion_pad_var = lets_read_this(func='mms_read_pad_ion_all', $
;    ion_pad_var = lets_read_this(func='mms_read_pad_ion_thermal_fpi', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion thermal data ...'
        return, retval
    endif
    ion_kev_pad_var = lets_read_this(func='mms_read_pad_ion_kev', $
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion kev data ...'
        return, retval
    endif

    ion_pad_vars = list()
    foreach species, ['o','p'] do begin
        ion_pad_vars.add, lets_read_this(func='mms_read_pad_ion_thermal_hpca', $
            time_range, probe=mission_probe, errmsg=errmsg, species=species)
    endforeach

    
    ; en spec.
    en_spec_vars = list()
    pa_spec_vars = list()
    pad_vars = prefix+['e_pad_thermal','e_pad_kev','o_pad_thermal_hpca',$
        'p_pad_thermal_hpca','ion_pad_thermal','p_pad_kev']
    foreach pad_var, pad_vars do begin
        en_spec_vars.add, pad_get_en_spec(pad_var=pad_var)
        pa_spec_vars.add, pad_get_pa_spec(pad_var=pad_var)
    endforeach
    en_spec_vars = en_spec_vars.toarray()
    pa_spec_vars = pa_spec_vars.toarray()
    
    ele_pa_en_low = 80
    ele_pa_en = 800
    pad_var = prefix+'e_pad_thermal'
    energy_range = [ele_pa_en_low,ele_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_low')
    
    energy_range = [ele_pa_en,50000]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_high')
    
    ion_pa_en_low = 50
    ion_pa_en = 800
    pad_var = prefix+'ion_pad_thermal'
    energy_range = [ion_pa_en_low,ion_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_low')

    energy_range = [ion_pa_en,50000]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_high')
    
    vars = prefix+'e_pa_spec_thermal_low'
    options, vars, zrange=[1e6,1e7]*1.5
    vars = prefix+'e_pa_spec_thermal_high'
    options, vars, zrange=[1e5,1e6]*3
    vars = prefix+'e_en_spec_thermal'
    options, vars, zrange=[1e3,1e8], constant=[ele_pa_en_low,ele_pa_en]

    vars = prefix+'ion_pa_spec_thermal_low'
    options, vars, zrange=[1e4,1e5]*2
    vars = prefix+'ion_pa_spec_thermal_high'
    options, vars, zrange=[1e4,1e5]*1
    vars = prefix+'ion_en_spec_thermal'
    options, vars, zrange=[1e4,1e6], constant=[ion_pa_en_low,ion_pa_en]
    
    
    pa_suffix = '!C    PA'
    plot_vars = prefix+['e_en_spec_kev','e_pa_spec_kev',$
        'e_en_spec_thermal','e_pa_spec_thermal_high','e_pa_spec_thermal_low', $
        'ion_en_spec_thermal','ion_pa_spec_thermal_high','ion_pa_spec_thermal_low']
    panel_labels = ['Ele keV EN','Ele keV'+pa_suffix,'Ele EN','Ele mid'+pa_suffix,'Ele low'+pa_suffix,$
        'Ion EN','Ion mid'+pa_suffix,'Ion low'+pa_suffix]
    
    
    tplot_options, 'tickinterval', 1200
    fig_info = sgplot(plot_vars, panel_labels=panel_labels, xrange=time_range, filename=plot_file)
    
    if event.haskey('pad_times') then begin
        bar_times = event['pad_times']
        panel_info = fig_info.panel_info
        tpos = (panel_info[plot_vars[0]])['position']
        tpos[1] = ((panel_info[plot_vars[-1]])['position'])[1]
        yrange = [0,1]
        set_axis, position=tpos, xrange=time_range, yrange=yrange
        foreach bar_time, bar_times do begin
            txs = time_double(bar_time)+[0,0]
            plots, txs, yrange, color=sgcolor('red'), linestyle=3, data=1
        endforeach
    endif
    
    if keyword_set(test) then stop
    sgclose

stop

;---Make the plot.
    plot_tr = time_range
    pad_times = event.pad_times
    npad_time = n_elements(pad_times)
    
    poss = panel_pos(0, nxpan=2, nypan=1, pansize=[1,1]*3, xpad=15, fig_size=fig_size, margins=[10,4,10,4])
    ncolor = 25
    ct2 = 64
    ct = 40
    options, [ele_pad_var,ion_pad_var], 'mission', 'mms'
    for tid=0,npad_time-1 do begin
        pad_time = pad_times[tid]
        base = project_id+'_fig_pad_'+time_string(pad_time,tformat='YYYY_MMDD_hhmm')+'_v01.pdf'
        pad_plot_file = join_path([plot_dir,base])
        if keyword_set(test) then pad_plot_file = 0
        ;if file_test(pad_plot_file) eq 1 then continue
        
        sgopen, pad_plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
        tmp = plot_pad_polygon(ele_pad_var, test=test, plot_times=pad_time, zrange=[5e1,5e7], color_table=ct, ncolor=ncolor, position=poss[*,0])
        tmp = plot_pad_polygon(ion_pad_var, test=test, plot_times=pad_time, zrange=[5e0,5e6]*0.2, color_table=ct, ncolor=ncolor, position=poss[*,1])
        ;tmp = plot_pad_polygon(ion_kev_pad_var, test=test, plot_times=pad_time)

        tpos = poss[*,0]
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*0.7
        msg = 'a) Ele'
        xyouts, tx,ty,msg, normal=1
                
        tpos = poss[*,1]
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*0.7
        msg = 'a) Ion'
        xyouts, tx,ty,msg, normal=1        
        
        color = sgcolor('red')
        foreach pid, [0,1] do begin
            tpos = poss[*,pid]
            xrange = [-1,1]
            yrange = [-1,1]
            set_axis, position=tpos, xrange=xrange, yrange=yrange
            
            tr = 0.62
            tx =-tr
            ty = tr
            msg = 'FEEPS'
            xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
            
            tr = 0.28
            tx =-tr
            ty = tr
            msg = 'FPI'
            xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
        endforeach
        if keyword_set(test) then stop
        sgclose
    endfor

    return, plot_file

stop


;---Make the movie.
    poss = panel_pos(plot_file, nxpan=2, nypan=1, pansize=[1,1]*3, xpad=15, fig_size=fig_size, margins=[10,4,10,4])
    times = get_var_time(ele_pad_var, in=time_range+[1,-1]*60)

    plot_files = list()
    foreach pad_time, times do begin
        base = project_id+'fig_pad_'+time_string(pad_time,tformat='YYYY_MMDD_hhmm_ss')+'_v01.png'
        plot_file = join_path([plot_dir,'pad_figures',base])
        plot_files.add, plot_file
        if file_test(plot_file) eq 1 then continue

        if keyword_set(test) then plot_file = 0
        sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz, magnify=2
        tmp = plot_pad_polygon(ele_pad_var, test=test, plot_times=pad_time, zrange=[5e1,5e8], color_table=ct, ncolor=ncolor, position=poss[*,0])
        tmp = plot_pad_polygon(ion_pad_var, test=test, plot_times=pad_time, zrange=[5e0,5e6], color_table=ct, ncolor=ncolor, position=poss[*,1])

        tpos = poss[*,0]
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*0.7
        msg = 'a) Ele'
        xyouts, tx,ty,msg, normal=1

        tpos = poss[*,1]
        tx = tpos[0]-xchsz*6
        ty = tpos[3]-ychsz*0.7
        msg = 'a) Ion'
        xyouts, tx,ty,msg, normal=1

        color = sgcolor('red')
        foreach pid, [0,1] do begin
            tpos = poss[*,pid]
            xrange = [-1,1]
            yrange = [-1,1]
            set_axis, position=tpos, xrange=xrange, yrange=yrange

            tr = 0.62
            tx =-tr
            ty = tr
            msg = 'FEEPS'
            xyouts, tx,ty,msg, color=color, data=1, alignment=0.5

            tr = 0.28
            tx =-tr
            ty = tr
            msg = 'FPI'
            xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
        endforeach
        if keyword_set(test) then stop
        sgclose
    endforeach
    
    movie_file = join_path([plot_dir,event_id+'_pad_'+prefix+'_'+version+'.mp4'])
    fig2movie, movie_file, fig_files=plot_files.toarray()


    return, plot_file

end

probe = '1'
event_id = '2015_0901_18'
print, micro_injection_fig_pad_v01(event_id, probe=probe, test=0, update=1)
stop

probe = '1'
event_id = '2015_0901_11'
print, micro_injection_fig_pad_v01(event_id, probe=probe, test=1)
end
