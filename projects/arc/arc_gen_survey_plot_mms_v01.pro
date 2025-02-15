

function arc_gen_survey_plot_mms_v01, input_event_id, test=test, get_name=get_name, errmsg=errmsg

    errmsg = ''
    retval = !null
    version = 'v01'
    project = arc_load_project()
    project_id = project.id

    if n_elements(input_event_id) eq 2 then begin
        time_range = time_double(input_event_id)
        event_id = time_string(time_range[0],tformat='YYYY_MMDD_hh')
        event = project_add_event(project, time_range=time_range, update=1)
    endif else begin
        event_id = input_event_id
        event = project_get_event(project, id=event_id)
    endelse
    time_range = event.time_range
    if n_elements(event) eq 0 then message, 'Inconsistency ...'

    if n_elements(plot_dir) eq 0 then plot_dir = event.plot_dir
    base = project_id+'_survey_plot_'+event_id+'_'+version+'.pdf'
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
    
    foreach phys_quant, ['bfield','efield','orbit'] do begin
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
                data0 = data[*,0]
                for jj=0,nprobe-1 do data[*,jj] -= data0
                ;data[*,0] = !values.f_nan
                data *= constant('re')
                store_data, var, times, data
                options, var, labels=labels, colors=colors, ytitle='(km)', labflag=-1
            endfor
        endif
    endforeach
    
    
    
    stop

    if keyword_set(test) then plot_file = 0
    return, plot_file

end


event_id = '2016_0902_00'
event_id = ['2016-09-07','2016-09-07/05:00']
print, arc_gen_survey_plot_mms_v01(event_id, test=1)
end