
function micro_injection_fig_pad2d_v01, input_event_id, probe=probe, $
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
    base = project_id+'_fig_pad2d_'+event_id+'_mms'+probe+'_'+version+'.pdf'
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
        time_range, probe=mission_probe, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion thermal data ...'
        return, retval
    endif
    
    unit_type = 'nflux'
    
    foreach var, prefix+['e','p']+'_pad_all' do begin
        get_data, var, times, data, limits=lim
        vals = lim.en_centers
        foreach val, vals, vid do begin
            if unit_type eq 'nflux' then begin
                
            endif else if unit_type eq 'eflux' then begin
                data[*,*,vid] *= val*1e-3
            endif else if unit_type eq 'xflux' then begin
                data[*,*,vid] *= sqrt(val*1e-3)
            endif
        endforeach
        
        if var eq prefix+'p_pad_all' then begin
            index = where(vals ge 2.5e4)
            data[*,*,index] *= 5
        endif
        
        ; remove some useless energy bins
        index = where(vals le 1.5e5)
        data = data[*,*,index]
        vals = vals[index]
        
        var1 = var+'_eflux'
        store_data, var1, times, data, vals, limits=lim
        if unit_type eq 'nflux' then begin
            unit = '#/cm!U2!N-s-sr-keV'
        endif else if unit_type eq 'eflux' then begin
            unit = 'keV/cm!U2!N-s-sr-keV'
        endif else if unit_type eq 'xflux' then begin
            unit = 'keV!U1/2!N/cm!U2!N-s-sr-keV'
        endif
        options, var1, unit=unit, en_centers=vals
    endforeach
    
    
;---Settings.
    ele_pa_en_low = 80
    ele_pa_en = 800
    ion_pa_en_low = 80
    ion_pa_en = 4000
    ion_pa_high = 2.5e4
    ele_pa_high = ion_pa_high

    ele_pad_var = prefix+'e_pad_all_eflux'
    ion_pad_var = prefix+'p_pad_all_eflux'
    ct = 40
    ;ct = 64
    
    if unit_type eq 'nflux' then begin
        ele_zrange = [1e1,1e8]
        ion_zrange = [2e1,1e6]
    endif else if unit_type eq 'eflux' then begin
        ; for keV/xxx.
        ele_zrange = [1e3,1e7]*2
        ion_zrange = [1e2,1e6]*40
    endif else if unit_type eq 'xflux' then begin
        ; for keV^0.5/xxx.
        ele_zrange = [1e2,1e7]*4.5
        ion_zrange = [1e3,1e5]*5
    endif
    
    
    
    pad_times = [$
        '2015-09-01/19:58:00',$
        '2015-09-01/20:00:20',$
        '2015-09-01/20:02:40' ]
;        '2015-09-01/20:30:20',$
;        '2015-09-01/20:38:20' ]
    event['pad_times'] = time_double(pad_times)
    

;---Make the plot.
    plot_tr = time_range
    pad_times = event.pad_times
    npad_time = n_elements(pad_times)
    nypan = 2
    nxpan = npad_time
    
    poss = panel_pos(0, nxpan=nxpan, nypan=nypan, pansize=[1,1]*1.2, xpad=1, ypad=1, $
        fig_size=fig_size, margins=[8,4,8,1.5])

    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    ncolor = 25
    options, [ele_pad_var,ion_pad_var], 'mission', 'mms'
    
    letter = 'f'
    for tid=0,npad_time-1 do begin
        pad_time = pad_times[tid]
        if keyword_set(test) then pad_plot_file = 0
        ;if file_test(pad_plot_file) eq 1 then continue
        
        no_cb = (tid eq npad_time-1)? 0: 1
        ytitle = (tid eq 0)? 'Perp E (eV)': ' '
        ytickformat = (tid eq 0)? '': '(A1)'
        xtitle = ' '
        xtickformat = '(A1)'
        title = strupcase('mms-'+probe)+' '+time_string(pad_time)+' UT'
        constants = [ele_pa_en_low,ele_pa_en,ele_pa_high]
        tmp = plot_pad_polygon(ele_pad_var, test=0, plot_times=pad_time, $
            ytitle=ytitle, ytickformat=ytickformat, $
            xtitle=xtitle, xtickformat=xtickformat, $
            zrange=ele_zrange, color_table=ct, ncolor=ncolor, $
            constants=constants, title=title, $
            position=poss[*,tid,0], no_colorbar=no_cb)
        
        
        xtitle = 'Para E (eV)'
        xtickformat = ''
        title = ''
        constants = [ion_pa_en_low,ion_pa_en,ion_pa_high]
        tmp = plot_pad_polygon(ion_pad_var, test=0, plot_times=pad_time, $
            ytitle=ytitle, ytickformat=ytickformat, $
            xtitle=xtitle, xtickformat=xtickformat, $
            zrange=ion_zrange, color_table=ct, ncolor=ncolor, $
            constants=constants, title=title, $
            position=poss[*,tid,1], no_colorbar=no_cb)
        ;tmp = plot_pad_polygon(ion_kev_pad_var, test=test, plot_times=pad_time)     
        
        
    ;---Add label for FEEPS and FPI.
        if tid eq 1 then begin
        color = sgcolor('red')
            foreach pid, [0,1] do begin
                tpos = poss[*,tid,pid]
                xrange = [-1,1]
                yrange = [-1,1]
                set_axis, position=tpos, xrange=xrange, yrange=yrange
                
                tr = 0.6
                tx =-tr
                ty = tr
                msg = 'FEEPS'
                xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
                
                tr = 0.25
                tx =-tr
                ty = tr
                msg = 'FPI'
                xyouts, tx,ty,msg, color=color, data=1, alignment=0.5
            endforeach
        endif
        
        

    ;---Add label.
        foreach pid, [0,1] do begin
            tpos = poss[*,tid,pid]

            num = tid*2+pid+1
            num_str = string(num,format='(I0)')
            msg = letter+'-'+num_str+')'
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty,msg, normal=1
        endforeach
    endfor
    if keyword_set(test) then stop
    sgclose

    return, plot_file

end

probe = '4'
event_id = '2015_0901_18'
print, micro_injection_fig_pad2d_v01(event_id, probe=probe, test=1, update=1)
end

