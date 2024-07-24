
function mms_gen_survey_plot_v01, input_time_range, probe=probe, $
    plot_dir=plot_dir, position=full_pos, errmsg=errmsg, test=test, xpansize=xpansize, $
    version=version, local_root=local_root, get_name=get_name

    errmsg = ''
    retval = !null
    version = 'v01'

    time_range = time_double(input_time_range)
    ; This is just to use the new disk for thg b/c /data is almost full.
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'themis','thg','survey_plot','mms'])

    if n_elements(plot_dir) eq 0 then plot_dir = join_path([local_root,'%Y'])
    path = apply_time_to_pattern(plot_dir,time_range[0])
    base = 'mms_survey_plot_'+strjoin(time_string(time_range,tformat='YYYY_MMDD'),'_')+'_mms'+probe+'_'+version+'.pdf'
    plot_file = join_path([path,base])
    if keyword_set(get_name) then return, plot_file
    print, plot_file
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

    ; Orbit related vars.
    orbit_time_range = time_range+[-1,1]*30*60
    r_gsm_var = lets_read_this(func='mms_read_orbit', $
        orbit_time_range, probe=mission_probe, coord=default_coord)
    print, 'Loading '+r_gsm_var+' ...'
    mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
    foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'

    ; Model related vars.
    external_models = ['t89','t96','t01','t04s']
    internal_models = ['dipole','igrf']
    hemispheres = ['north','south']

    foreach external_model, external_models do begin
        foreach internal_model, internal_models do begin
            ; The B at sc position.
            suffix = '_'+internal_model+'_'+external_model
            bmod_var = prefix+'bmod_gsm'+suffix
            if tnames(bmod_var) ne '' then del_data, bmod_var
            bmod_var = lets_read_geopack_bfield(var_info=bmod_var, $
                orbit_var=r_gsm_var, time_var=orbit_time_var, $
                internal_model=internal_model, external_model=external_model, save_to=data_file, update=update)
            print, 'Loading '+bmod_var+' ...'
        endforeach
    endforeach

    
    ; B field related vars.
    field_time_range = time_range+[-1,1]*30.*60
    b_gsm_var = lets_read_this(func='mms_read_bfield', $
        field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load B field data ...'
        return, retval
    endif
    e_gsm_var = lets_read_this(func='mms_read_efield', $
        field_time_range, probe=mission_probe, coord=default_coord, errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load E field data ...'
        return, retval
    endif

    b0_window = 15.*60
    bmod_var = prefix+'bmod_gsm_igrf_t89'
    b_vars = lets_decompose_bfield(b0_window=b0_window, b_var=b_gsm_var, bmod_var=bmod_var)
    b0_gsm_var = b_vars['b0']
    b_elev_var = lets_calc_vec_elev(b_gsm_var, coord='sm')
    bmod_elev_var = lets_calc_vec_elev(bmod_var, coord='sm', var_info=prefix+'bmod_elev')
    db_elev_var = lets_subtract_vars(b_elev_var, bmod_elev_var, save_to=prefix+'db_elev')
    options, db_elev_var, constant=0, yrange=[-1,1]*90


    ; Particle related vars.
    ele_en_spec_var = lets_read_this(func='mms_read_en_spec_ele', $
        time_range, probe=mission_probe, id='thermal', errmsg=errmsg)
    if errmsg eq '' then begin
        errmsg = 'Failed to load e thermal data ...'
        return, retval
    endif
    options, ele_en_spec_var, $
        zrange=[1e4,1e7], zstyle=1, zlog=1, ztickv=[1e4,1e5,1e6,1e7], ztickname='10!U'+['4','5','6','7'], zticks=3, zminor=9, $
        yrange=[1.1e1,2.6e4], ystyle=1, ylog=1, ytickv=[1e2,1e3,1e4], ytickname='10!U'+['2','3','4'], yticks=2, yminor=9
    ele_kev_en_spec_var = lets_read_this(func='mms_read_en_spec_ele', $
        time_range, probe=mission_probe, id='kev', errmsg=errmsg)
    if errmsg eq '' then begin
        errmsg = 'Failed to load e kev data ...'
        return, retval
    endif
    options, ele_kev_en_spec_var, $
        zrange=[1e1,1e5], zstyle=1, zlog=1, ztickv=[1e1,1e2,1e3,1e4,1e5], ztickname=['10!U1',' ','10!U3',' ','10!U5'], zticks=4, zminor=9, $
        yrange=[4.7e4,5.2e5], ystyle=1, ylog=1, ytickv=[5e4,5e5], ytickname='10!U'+['4','5'], yticks=1, yminor=9

    
    


;---Make the plot.
    index = where(get_var_data(prefix+'dis',times=times) ge 9, count)
    if count eq 0 then return, retval
    plot_tr = minmax(times[index])
    tickinterval = 30*60d
    

    ; Init plot_vars.
    plot_info = orderedhash()
    
    
    plot_info[ele_kev_en_spec_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'e- high', $
        'ypan', 0.8, $
        'setting', dictionary( ) )
    plot_info[ele_en_spec_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'e- low', $
        'setting', dictionary( ) )
    
    plot_info[b_gsm_var] = dictionary($
        'routine', 'plot_linlog', $
        'panel_label_text', 'B GSM', $
        'setting', dictionary($
            'linear_yrange', [-1,1]*20, $
            'log_yrange', [-1,1]*200, $
            'linear_tick_setting', dictionary($
                'yticks', 2, $
                'ytickv', [-1,0,1]*10, $
                'yminor', 4 ), $
            'log_tick_setting', dictionary($
                'yticks', 1, $
                'ytickv', [1,10]*20, $
                'yminor', 9 ) $
        ) $
    )
    
    plot_info[e_gsm_var] = dictionary($
        'routine', 'plot_linlog', $
        'panel_label_text', 'E GSM', $
        'setting', dictionary($
            'linear_yrange', [-1,1]*5, $
            'log_yrange', [-1,1]*500, $
            'linear_tick_setting', dictionary($
                'yticks', 2, $
                'ytickv', [-1,0,1]*5, $
                'yminor', 5 ), $
            'log_tick_setting', dictionary($
                'yticks', 1, $
                'ytickv', [10,100]*5, $
                'yminor', 9 ) $
            ) $
        )

    
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
    plot_poss = panel_pos(plot_file, nypan=nplot_var, fig_size=fig_size, ypans=ypans, pansize=[12,1.2], margins=margins)
    
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

        plot_setting = plot_setting.tostruct()        
        tmp = call_function(plot_routine, plot_var, _extra=plot_setting)
    endforeach


    if keyword_set(test) then stop
    sgclose
    

    return, plot_file

end

; Survey round 1
root_dir = join_path(homedir(),'micro_injection_survey_round1')
survey_file = 'micro_injection_survey_round1_survey_info.sav'
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
        
        ; if the current time_range has been surveyed, then there should be an id.
        id = time_string(tr[0],tformat='YYYY_MMDD_hh')
        if survey_info.haskey(id) then begin
            my_info = survey_info[id]
        endif else my_info = dictionary()
        
        ; if the current mission_probe has been surveyed, then it should be a key there.
        if my_info.haskey(mission_probe) then begin
            probe_info = my_info[mission_probe]
        endif else begin
            probe_info = dictionary()
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
            survey_plot = mms_gen_survey_plot_v01(tr, probe=probe, test=test, errmsg=errmsg, get_name=1)
        endelse
        
        ; We have no errsmg or survey plot, just do the survey and get them.
        if file_test(survey_plot) eq 0 then begin
            survey_plot = mms_gen_survey_plot_v01(tr, probe=probe, test=test, errmsg=errmsg)
            
            if errmsg ne '' then begin
                probe_info['errmsg'] = errmsg
            endif else begin
                if file_test(survey_plot) eq 0 then begin
                    probe_info['errmsg'] = 'No survey plot generated ...'
                endif else begin
                    probe_info['survey_plot'] = survey_plot
                endelse
            endelse
            ; Update the results.
            save, survey_info, filename=survey_file
        endif else begin
            ; pass. We have the survey plot, meaning everything is good.
        endelse
        
        ; Start to survey next time range.
    endfor
endforeach
end