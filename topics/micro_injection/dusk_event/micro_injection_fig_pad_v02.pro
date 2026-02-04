

function micro_injection_fig_pad_v02, input_event_id, probe=probe, $
    plot_dir=plot_dir, test=test, get_name=get_name, update=update

    errmsg = ''
    retval = !null
    version = 'v02'
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

;    ion_pad_vars = list()
;    foreach species, ['o','p'] do begin
;        ion_pad_vars.add, lets_read_this(func='mms_read_pad_ion_thermal_hpca', $
;            time_range, probe=mission_probe, errmsg=errmsg, species=species)
;    endforeach

    
    ; en spec.
    en_spec_vars = list()
    pa_spec_vars = list()
    pad_vars = prefix+['e_pad_thermal','e_pad_kev','ion_pad_thermal','p_pad_kev']
    foreach pad_var, pad_vars do begin
        en_spec_vars.add, pad_get_en_spec(pad_var=pad_var)
        var = pad_get_pa_spec(pad_var=pad_var)
        pa_spec_vars.add, var
        options, var, 'energy_range', minmax(get_var_setting(pad_var,'en_centers'))
    endforeach
    en_spec_vars = en_spec_vars.toarray()
    pa_spec_vars = pa_spec_vars.toarray()
    
    ele_pa_en_low = 80
    ele_pa_en = 800
    pad_var = prefix+'e_pad_thermal'
    ele_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ele_pa_en_low,ele_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range
    
    energy_range = [ele_pa_en,ele_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'e_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    
    ion_pa_en_low = 80
    ion_pa_en = 4000
    pad_var = prefix+'ion_pad_thermal'
    ion_en_max = max(get_var_setting(pad_var, 'en_centers'))
    energy_range = [ion_pa_en_low,ion_pa_en]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_low')
    options, tmp, 'energy_range', energy_range

    energy_range = [ion_pa_en,ion_en_max]
    tmp = pad_get_pa_spec(pad_var=pad_var, energy_range=energy_range, var_info=prefix+'ion_pa_spec_thermal_high')
    options, tmp, 'energy_range', energy_range
    

    log_ytickv = [2,3,4]
    yticks = n_elements(log_ytickv)-1
    zrange = [1e3,1e8]
    log_ztickv = [3,4,5,6,7,8]
    zticks = n_elements(log_ztickv)-1
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    ztickn[0:*:2] = ' '
    vars = prefix+'e_en_spec_thermal'
    options, vars, zrange=zrange, constant=[ele_pa_en_low,ele_pa_en], $
        ytickv=10d^log_ytickv, yticks=yticks, ytickname='10!U'+string(log_ytickv,format='(I0)'), $
        ztickv=10d^log_ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9
        
    
    vars = prefix+'e_en_spec_kev'
    log_ytickv = [4,5]
    yfactor = 5
    yrange = yfactor*10d^log_ytickv
    ytickn = string(yfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ytickv,format='(I0)')
    ytickv = 1e5
    yticks = n_elements(ytickv)-1
    ytickn = '10!U5'
    log_zrange = [1,5]-1
    zfactor = 5
    zrange = zfactor*10d^log_zrange
    log_ztickv = make_bins(log_zrange,1, inner=1)
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_ztickv,format='(I0)')
    index = where(log_ztickv eq 1, count)
    if count ne 0 then ztickn[index] = '10'
    ztickn[1:*:2] = ' '
    options, vars, zrange=zrange, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, ytickname=ytickn, $
        ztickv=ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9

    vars = prefix+'ion_pa_spec_thermal_low'
    zfactor = 1.3
    log_zrange = [4,5]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, vars, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    
    vars = prefix+'ion_pa_spec_thermal_high'
    zfactor = 8
    log_zrange = [4,5]-1
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(I0)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, vars, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    vars = prefix+'ion_en_spec_thermal'
    log_zrange = [4,6]
    zfactor = 1
    zrange = zfactor*10d^log_zrange
    log_yrange = [1,4]
    yfactor = 2
    yrange = yfactor*10d^log_yrange
    log_ytickv = make_bins(minmax(alog10(yrange)),1,inner=1)
    ytickv = 10d^log_ytickv
    yticks = n_elements(ytickv)-1
    ytickn = '10!U'+string(log_ytickv,format='(I0)')
    options, vars, constant=[ion_pa_en_low,ion_pa_en], $
        yrange=yrange, ytickv=ytickv, yticks=yticks, ytickname=ytickn, yminor=9, $
        zrange=zrange
    
    var = prefix+'e_pa_spec_kev'
    zfactor = 3.5
    log_zrange = [2,3]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    var = prefix+'e_pa_spec_thermal_high'
    zfactor = 2.8
    log_zrange = [5,6]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    var = prefix+'e_pa_spec_thermal_low'
    zfactor = 1.5
    log_zrange = [6,7]
    zrange = 10d^log_zrange*zfactor
    ztickv = zrange
    ztickn = string(zfactor,format='(F3.1)')+tex2str('times')+'10!U'+string(log_zrange,format='(I0)')
    zticks = n_elements(ztickv)-1
    zminor = 9
    options, var, zrange=zrange, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=zminor
    
    
    pa_suffix = ''
    plot_vars = prefix+[$
        'e_en_spec_kev','e_en_spec_thermal', $
        'e_pa_spec_kev','e_pa_spec_thermal_high','e_pa_spec_thermal_low', $
        'ion_en_spec_thermal', $
        'ion_pa_spec_thermal_high','ion_pa_spec_thermal_low']
    panel_labels = [$
        'Ele EN high','Ele EN',$
        'Ele high'+pa_suffix,'Ele mid'+pa_suffix,'Ele low'+pa_suffix,$
        'Ion EN', $
        'Ion mid'+pa_suffix,'Ion low'+pa_suffix]
    
    
    tplot_options, 'tickinterval', 1200
    fig_info = sgplot(plot_vars, panel_labels=panel_labels, xrange=time_range, filename=plot_file)
    
    pa_vars = prefix+[['e','ion']+'_pa_spec_kev','e_pa_spec_thermal_'+['low','high'],'ion_pa_spec_thermal_'+['low','high']]
    panel_info = fig_info['panel_info']
    xchsz = fig_info.xchsz
    ychsz = fig_info.ychsz
    foreach var, pa_vars do begin
        if not panel_info.haskey(var) then continue
        the_info = panel_info[var]
        tpos = the_info['position']
        energy_range = get_var_setting(var,'energy_range')*1e-3
        msg = strjoin(strtrim(string(energy_range,format='(F5.1)'),2),'-')+' keV'
        index = strpos(var, 'kev')
        if index[0] ne -1 then msg = strjoin(strtrim(string(energy_range,format='(I0)'),2),'-')+' keV'
        tx = tpos[2]-xchsz*4.5
        ty = tpos[3]-ychsz*1
        polyfill, tpos[2]-[0.5,8.5,8.5,0.5,0.5]*xchsz, tpos[3]-[0.2,0.2,1.2,1.2,0.2]*ychsz, normal=1, color=sgcolor('white')
        xyouts, tx,ty,msg, normal=1, alignment=0.5, color=sgcolor('black')
    endforeach
    
    
    
    pad_times = [$
        '2015-09-01/18:12:40',$
        '2015-09-01/18:18:00',$
        '2015-09-01/18:59:20',$
        '2015-09-01/19:06:00',$
        ;'2015-09-01/19:17:40',$
        '2015-09-01/19:22:00',$
        '2015-09-01/19:27:00',$
        '2015-09-01/19:58:00',$
        '2015-09-01/20:02:40',$
        '2015-09-01/20:30:20',$
        '2015-09-01/20:38:20' ]
    event['pad_times'] = time_double(pad_times)
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

probe = '4'
event_id = '2015_0901_18'
print, micro_injection_fig_pad_v02(event_id, probe=probe, test=0)
end
