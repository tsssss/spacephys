
function micro_injection_survey_round1

    ; Survey round 1.
    project_id = 'micro_injection'
    project_info = project_load(project_id)
    root_dir = project_info.root_dir
    survey_file = join_path([project_info.data_dir,project_id+'_survey_round1_survey_info.sav'])
    if file_test(survey_file) eq 0 then begin
        survey_info = orderedhash()
        save, survey_info, filename=survey_file
    endif
    restore, survey_file

    trs = micro_injection_load_survey_time_range()
    probes = ['1','2','3','4']
    test = 0
    ntr = n_elements(trs[*,0])
    foreach probe, probes do begin
        mission_probe = 'mms'+probe
        for ii=0,ntr-1 do begin
            tr = reform(trs[ii,*])
            
            msg = 'Processing '+strupcase('mms'+probe)+' from '+strjoin(time_string(tr),' to ')+' ...'
            lprmsg, msg
            
            ; if the current time_range has been surveyed, then there should be an id.
            id = time_string(tr[0],tformat='YYYY_MMDD_hh')
            if survey_info.haskey(id) then begin
                my_info = survey_info[id]
            endif else begin
                my_info = dictionary()
                survey_info[id] = my_info
            endelse
            
            
            ; if the current mission_probe has been surveyed, then it should be a key there.
            if my_info.haskey(mission_probe) then begin
                probe_info = my_info[mission_probe]
            endif else begin
                probe_info = dictionary()
                my_info[mission_probe] = probe_info
            endelse
            
            ; if there is an errmsg, then it means it has been surveyed but there is an error.
            if probe_info.haskey('errmsg') then begin
                errmsg = probe_info['errmsg']
            endif else begin
                errmsg = !null
            endelse
            ; Already surveyed, pass.
            if n_elements(errmsg) ne 0 then continue
            
            ; check if a survey plot is there.
            if probe_info.haskey('survey_plot') then begin
                survey_plot = probe_info['survey_plot']
            endif else begin
                survey_plot = micro_injection_gen_survey_plot_v01(tr, probe=probe, test=test, errmsg=errmsg, get_name=1)
            endelse
            
            ; We have no errsmg or survey plot, just do the survey and get them.
            if file_test(survey_plot) eq 0 then begin
                survey_plot = micro_injection_gen_survey_plot_v01(tr, probe=probe, test=test, errmsg=errmsg)
                
                if errmsg ne '' then begin
                    probe_info['errmsg'] = errmsg
                endif else begin
                    if file_test(survey_plot) eq 0 then begin
                        probe_info['errmsg'] = 'No survey plot generated ...'
                    endif else begin
                        probe_info['survey_plot'] = survey_plot
                    endelse
                endelse
            endif else begin
                ; pass. We have the survey plot, meaning everything is good.
                probe_info['survey_plot'] = survey_plot
            endelse 
            
            ; Update the results.
            save, survey_info, filename=survey_file
            ; Start to survey next time range.
        endfor
    endforeach

    return, 'Done'

end


print, micro_injection_survey_round1()
end