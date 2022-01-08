;+
; Plot the result of velocity from 2d time lag.
;-

pro azim_prop_print_table, project

    if n_elements(project) eq 0 then project = azim_prop_load_project()
    events = project.events
    table_file = join_path([project.root_dir,'azim_prop_table.txt'])
    if file_test(table_file) then file_delete, table_file
    ftouch, table_file

    test = 0

;---Plot settings.
    re = project.constant.re
    deg = project.constant.deg
    rad = project.constant.rad
    tab = ' & '
    str_pm = '!9'+string(177b)+'!X'
    str_theta = '!9'+string(113b)+'!X'
    str_delta = '!9'+string(68b)+'!X'
    str_omega = '!9'+string(119b)+'!X'
    r_kms = 6.6
    coef2kms = r_kms*re*rad/60


    sorted_events = (events.keys()).toarray()
    index = sort(time_double(sorted_events,tformat='YYYY_MM_DD'))
    sorted_events = sorted_events[index]
    nevent = n_elements(sorted_events)


    msg_list = list()
    sep = '&'



;---Add event id.
    msgs = ['\multirow{2}{*}{Event ID}',' ']

    ; Get values.
    vals = strarr(nevent)
    foreach event_id, sorted_events, ii do begin
        tval = strmid(event_id,0,4)+'\_'+strmid(event_id,5,4)
        vals[ii] = tval
    endforeach

    ; Convert values to strings.
    ; do nothing.
    str_vals = vals
    str_median = 'Median'
    str_mean = 'Mean'
    str_vals = [str_vals,str_median,str_mean]

    ; Combine to title.
    msgs = extend_string([msgs,str_vals])+' '+sep+' '
    msg_list.add, msgs



;---Add extent.
    msgs = ['\multicolumn{6}{c}{DF extent seen by spacecraft}']
    keys = ['ref_mlt','ref_rsm','ref_time']
    formats = ['(F5.1)','(F3.1)','(I0)']
    nkey = n_elements(keys)
    titles = '\multicolumn{2}{c}'+['{dMLT (hr)}','{d$R_{xy}$ (Re)}','{dT (min)}']

    str_vals = strarr(nevent+2,nkey)
    vals = dblarr(nevent,nkey)
    foreach event_id, sorted_events, kk do begin
        event = project.events[event_id]
        probes = event.sorted_probes
        nprobe = n_elements(probes)

        ; Get values.
        tvals = dblarr(nprobe,nkey)
        foreach key, keys, ii do begin
            foreach probe, probes, jj do begin
                tval = (event[probe])[key]
                case key of
                    'ref_rsm': tval = snorm(tval[0:1])  ; r_xy.
                    'ref_time': tval = tval/60d   ; convert sec to hour.
                    else: ; do nothing.
                endcase
                tvals[jj,ii] = tval
            endforeach
        endforeach

        ; Convert values to string.
        foreach key, keys, ii do begin
            tval = max(tvals[*,ii])-min(tvals[*,ii])
            vals[kk,ii] = tval
            str_vals[kk,ii] = strtrim(string(tval,format=formats[ii]),2)
            index = strpos(str_vals[kk,ii],'.')
            if index[0] ne -1 then str_vals[kk,ii] = strjoin(strsplit(str_vals[kk,ii],'.',/extract),sep) else str_vals[kk,ii] += sep
        endforeach        
    endforeach
    
    foreach key, keys, ii do begin
        str_median = strtrim(string(median(vals[*,ii]),format=formats[ii]),2)
        index = strpos(str_median,'.')
        if index[0] ne -1 then str_median = strmid(str_median,0,index)+sep+strmid(str_median,index+1) else str_median += sep
        str_mean = strtrim(string(mean(vals[*,ii]),format=formats[ii]),2)
        index = strpos(str_mean,'.')
        if index[0] ne -1 then str_mean = strmid(str_mean,0,index)+sep+strmid(str_mean,index+1) else str_mean += sep
        str_mean += '$^{\pm'+strtrim(string(stddev(vals[*,ii]),format=formats[ii]),2)+'}$'
        str_vals[nevent:nevent+1,ii] = [str_median,str_mean]
    endforeach
    

    ; Combine to title.
    str_vals = [reform(titles,1,nkey),str_vals]
    foreach key, keys, ii do str_vals[*,ii] = extend_string(str_vals[*,ii])
    for ii=0, n_elements(str_vals)/nkey-1 do begin
        msgs = [msgs, strjoin(reform(str_vals[ii,*]),' '+sep+' ')]
    endfor

    msgs = extend_string(msgs)+' '+sep+' '
    msg_list.add, msgs



;---Add omega_azim and v_azim.
    msgs = ['']
    key = 'omega_azim'

    vals = dblarr(nevent)
    foreach event_id, sorted_events, ii do begin
        event = project.events[event_id]
        tval = abs(event[key])
        vals[ii] = tval
    endforeach

    str_vals = strarr(nevent)
    ;titles = ['\multicolumn{2}{c}{$|\omega_{\text{azim}}|$}','\multicolumn{2}{c}{($^{\circ}$/min)}']
    titles = ['\multicolumn{2}{c}{$|\omega_{\text{azim}}|$}','\multicolumn{2}{c}{(deg/min)}']
    format = '(F5.1)'
    for ii=0, nevent-1 do begin
        str_vals[ii] = strtrim(string(vals[ii], format=format),2)
        index = strpos(str_vals[ii],'.')
        if index[0] ne -1 then str_vals[ii] = strjoin(strsplit(str_vals[ii],'.',/extract),sep) else str_vals[ii] += sep
    endfor
    
    str_median = strtrim(string(median(vals),format=format),2)
    index = strpos(str_median,'.')
    if index[0] ne -1 then str_median = strmid(str_median,0,index)+sep+strmid(str_median,index+1)
    str_mean = strtrim(string(mean(vals),format=format),2)
    index = strpos(str_mean,'.')
    if index[0] ne -1 then str_mean = strmid(str_mean,0,index)+sep+strmid(str_mean,index+1)
    str_mean += '$^{\pm'+strtrim(string(stddev(vals),format=format),2)+'}$'
    str_vals = [str_vals,str_median,str_mean]
    
    msgs = extend_string([titles,str_vals])

    vals = vals*coef2kms
    str_vals = strarr(nevent)
    titles = ['\multicolumn{2}{c}{$|v_{\text{azim}}^{\text{'+string(r_kms,format='(F3.1)')+' Re}}|$}','\multicolumn{2}{c}{(km/s)}']
    format = '(I2)'
    for ii=0, nevent-1 do begin
        str_vals[ii] = string(vals[ii], format=format)
        index = strpos(str_vals[ii],'.')
        if index[0] ne -1 then str_vals[ii] = strjoin(strsplit(str_vals[ii],'.',/extract),sep) else str_vals[ii] += sep
    endfor
    
    str_median = strtrim(string(median(vals),format=format),2)
    index = strpos(str_median,'.')
    if index[0] ne -1 then str_median = strmid(str_median,0,index)+sep+strmid(str_median,index+1) else str_median += sep
    str_mean = strtrim(string(mean(vals),format=format),2)
    index = strpos(str_mean,'.')
    if index[0] ne -1 then str_mean = strmid(str_mean,0,index)+sep+strmid(str_mean,index+1) else str_mean += sep
    str_mean += '$^{\pm'+strtrim(string(stddev(vals),format=format),2)+'}$'
    str_vals = [str_vals,str_median,str_mean]
    
    msgs += ' '+sep+' '+extend_string([titles,str_vals])+ ' '+sep+' '
    msg_list.add, msgs


;---Add omega_2d.
    ;titles = ['\multicolumn{2}{c}{$|\omega_{\text{2D}}|$}','\multicolumn{2}{c}{($^{\circ}$/min)}']
    titles = ['\multicolumn{2}{c}{$|\omega_{\text{2D}}|$}','\multicolumn{2}{c}{(deg/min)}']
    keys = ['omega_2d','domega_2d']
    nkey = n_elements(keys)
    format = '(F4.1)'

    vals = dblarr(nevent,nkey)
    foreach event_id, sorted_events, ii do begin
        event = project.events[event_id]
        foreach key, keys, jj do if event.haskey(key) eq 0 then continue else vals[ii,jj] = abs(event[key])
    endforeach
    
    str_vals = strarr(nevent)
    foreach event_id, sorted_events, ii do begin
        if vals[ii,0]*vals[ii,1] eq 0 then begin
            str_vals[ii] = '-&'
            continue
        endif
        str_vals[ii] = strtrim(string(vals[ii,0],format=format),2)
        str_vals[ii] += '&$^{\pm'+strtrim(string(vals[ii,1],format=format),2)+'}$'
    endforeach
    
    str_median = '-&'
    str_mean = '-&'
    str_vals = [str_vals,str_median,str_mean]
    
    msgs = extend_string([titles,str_vals])+' '+sep+' '
    msg_list.add, msgs
    
    


;---Add scale.
    titles = ['\multicolumn{2}{c}{DF width}','\multicolumn{2}{c}{(Re)}']
    key = 'tilt_scale'
    format = '(F4.1)'
    
    vals = dblarr(nevent)
    foreach event_id, sorted_events, ii do begin
        event = project.events[event_id]
        vals[ii] = abs(event[key])
    endforeach
    
    str_vals = strarr(nevent)
    foreach event_id, sorted_events, ii do begin
        str_vals[ii] = strtrim(string(vals[ii],format=format),2)
        index = strpos(str_vals[ii],'.')
        if index[0] ne -1 then str_vals[ii] = strjoin(strsplit(str_vals[ii],'.',/extract),sep) else str_vals[ii] += sep
    endforeach
    str_median = strtrim(string(median(vals),format=format),2)
    index = strpos(str_median,'.')
    if index[0] ne -1 then str_median = strmid(str_median,0,index)+sep+strmid(str_median,index+1)
    str_mean = strtrim(string(mean(vals),format=format),2)
    index = strpos(str_mean,'.')
    if index[0] ne -1 then str_mean = strmid(str_mean,0,index)+sep+strmid(str_mean,index+1)
    str_mean += '$^{\pm'+strtrim(string(stddev(vals),format=format),2)+'}$'
    str_vals = [str_vals,str_median,str_mean]
    
    msgs = extend_string([titles,str_vals])+' \\'
    msg_list.add, msgs
    
    
    
;---Combine all columns.
    nrow = n_elements(msg_list[0])
    all_msgs = strarr(nrow)
    for ii=0, nrow-1 do begin
        foreach column, msg_list do all_msgs[ii] += column[ii]
    endfor
    
    all_msgs[0] += ' \cmidrule{2-7}'
    all_msgs[1] += ' \midrule'
    all_msgs[-3] += ' \midrule'
    all_msgs[-1] += ' \bottomrule'
    stop



end
