;+
; Generate overview figure for each event.
;-


function micro_injection_fig_overview_v03, input_event_id, probe=probe, $
    plot_dir=plot_dir, test=test, get_name=get_name, update=update

    errmsg = ''
    retval = !null
    version = 'v03'
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
    base = project_id+'_fig_overview_'+event_id+'_mms'+probe+'_'+version+'.pdf'
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
    
    
    ; B model.
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
    if errmsg ne '' then begin
        errmsg = 'Failed to load e thermal data ...'
        return, retval
    endif
    options, ele_en_spec_var, $
        zrange=[1e3,1e8], zstyle=1, zlog=1, ztickv=[1e4,1e5,1e6,1e7], ztickname='10!U'+['4','5','6','7'], zticks=3, zminor=9, $
        yrange=[1.1e1,2.6e4], ystyle=1, ylog=1, ytickv=[1e2,1e3,1e4], ytickname='10!U'+['2','3','4'], yticks=2, yminor=9
    ele_kev_en_spec_var = lets_read_this(func='mms_read_en_spec_ele', $
        time_range, probe=mission_probe, id='kev', errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load e kev data ...'
        return, retval
    endif


    ion_en_spec_var = lets_read_this(func='mms_read_en_spec_ion', $
        time_range, probe=mission_probe, id='thermal', errmsg=errmsg, instrument='fpi')
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion thermal data ...'
        return, retval
    endif
    options, ion_en_spec_var, $
        zrange=[1e4,1e6], zstyle=1, zlog=1, ztickv=[1e4,1e5,1e6], ztickname='10!U'+['4','5','6'], zticks=2, zminor=9, $
        yrange=[11,2.65e4], ystyle=1, ylog=1, ytickv=[1e2,1e3,1e4], ytickname='10!U'+['2','3','4'], yticks=2, yminor=9
    ion_kev_en_spec_var = lets_read_this(func='mms_read_en_spec_ion', $
        time_range, probe=mission_probe, id='kev', errmsg=errmsg)
    if errmsg ne '' then begin
        errmsg = 'Failed to load ion kev data ...'
        return, retval
    endif
    
    
    
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
    vars = [ele_kev_en_spec_var,ion_kev_en_spec_var]
    options, vars, zrange=zrange, yrange=yrange, $
        ytickv=ytickv, yticks=yticks, ytickname=ytickn, $
        ztickv=ztickv, zticks=zticks, ztickname=ztickn, zminor=9, yminor=9
;    options, ion_kev_en_spec_var, $
;        zrange=[1e1,1e4], zstyle=1, zlog=1, ztickv=[1e1,1e2,1e3,1e4], ztickname='10!U'+['1','2','3','4'], zticks=3, zminor=9, $
;        yrange=[4.7e4,5.2e5], ystyle=1, ylog=1, ytickv=[5e4,5e5], ytickname='5x10!U'+['4','5'], yticks=1, yminor=9
        
    ion_vel_var = lets_read_this(func='mms_read_ion_vel', $
        time_range, probe=mission_probe, errmsg=errmsg)

    
    tmp = omni_read_sw_b(time_range)
    get_data, tmp, times, b_gsm
    imf_bz_var = 'omni_imf_bz'
    store_data, imf_bz_var, times, b_gsm[*,2], limits={labels:'IMF Bz', ytitle:'(nT)', constant:0}
    
    
    
;    b0_var = b_vars['b0']
;    b1_var = b_vars['b1']
;
;    ; Convert to FAC.
;    q_fac_var = lets_define_fac(r_var=r_gsm_var, b_var=b0_var, fac_coord=fac_coord)
;    coord_msgs = [default_coord,fac_coord]
;    fac_vars = list()
;    foreach var, [b1_var,e_gsm_var,ion_vel_var] do begin
;        fac_vars.add, lets_cotran(coord_msgs, input=var, q_var=q_fac_var)
;    endforeach
;    b1_fac_var = fac_vars[0]
;    e1_fac_var = fac_vars[1]
;    b1_fac = get_var_data(b1_fac_var, times=times)
;    e1_fac = get_var_data(e1_fac_var, at=times)
;    pf_fac = spoynt(e1_fac,b1_fac)
;    pf_fac_var = prefix+'pf_fac'
;    store_data, pf_fac_var, times, pf_fac
;    fac_labels = ['||',tex2str('perp')+','+['west','out']]
;    add_setting, pf_fac_var, smart=1, dictionary($
;        'display_type', 'vector', $
;        'short_name', 'S', $
;        'coord', 'mms_fac', $
;        'coord_labels', fac_labels )



;---Make the plot.
    index = where(get_var_data(prefix+'dis',times=times) ge 9, count)
    if count eq 0 then return, retval
    plot_tr = time_range
    tickinterval = 60*60d
    

    ; Init plot_vars.
    plot_info = orderedhash()

    plot_info[imf_bz_var] = dictionary($
        'routine', 'plot_line', $
        'panel_label_text', 'IMF Bz', $
        'ypan', 0.6, $
        'setting', dictionary( $
            'yrange', [-1,1]*7, $
            'tick_setting', dictionary($
                'ytickv', [-1,0,1]*5, $
                'yminor', 5 ) ) )
    
    plot_info[ele_kev_en_spec_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'Ele high', $
        'ypan', 0.9, $
        'setting', dictionary( ) )
    plot_info[ele_en_spec_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'Ele low', $
        'setting', dictionary( ) )
    plot_info[ion_kev_en_spec_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'Ion high', $
        'ypan', 0.9, $
        'setting', dictionary( ) )
    plot_info[ion_en_spec_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'Ion low', $
        'setting', dictionary( ) )
    
    plot_info[b_gsm_var] = dictionary($
        'routine', 'plot_line', $
        'panel_label_text', 'B GSM', $
        'setting', dictionary($
            'yrange', [-10,50], $
            'plot_magnitude', 1, $
            'tick_setting', dictionary($
                'ytickv', [0,1,2]*20, $
                'yminor', 4 ) $
        ) $
    )
    
    plot_info[e_gsm_var] = dictionary($
        'routine', 'plot_line', $
        'panel_label_text', 'E GSM', $
        'setting', dictionary($
            'yrange', [-1,1]*5, $
            'tick_setting', dictionary($
                'ytickv', [-1,0,1]*4, $
                'yminor', 4 ) $
            ) $
        )
    plot_info[ion_vel_var] = dictionary($
        'routine', 'plot_line', $
        'panel_label_text', 'Ion Vel', $
        'setting', dictionary($
            'yrange', [-1,1]*140, $
            'plot_magnitude', 0, $
            'tick_setting', dictionary($
                'ytickv', [-1,0,1]*100, $
                'yminor', 5 ) $
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
    plot_poss = panel_pos(plot_file, nypan=nplot_var, fig_size=fig_size, ypans=ypans, pansize=[12,0.8]*0.55, margins=margins)
    
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
        
        if pid eq 0 then begin
            msg = strupcase('mms'+probe)
            color = sgcolor('white')
            tpos = my_pos
            tx = tpos[0]+xchsz*0.5
            ty = tpos[3]-ychsz*1
            xyouts, tx,ty,msg, normal=1, color=color
        endif
    endforeach

;    if event.haskey('pad_times') then begin
;        bar_times = event.pad_times
;        tpos = plot_poss[*,0]
;        tpos[1] = plot_poss[1,-1]
;        xrange = plot_tr
;        yrange = [0,1]
;        set_axis, xrange=xrange, yrange=yrange, position=tpos
;        color = sgcolor('red')
;        foreach bar_time, bar_times do begin
;            plots, bar_time+[0,0], yrange, linestyle=1, color=color
;        endforeach
;    endif

    

    tplot_options, get_options=opts
    str_element, opts, 'tickinterval', delete=1

    if keyword_set(test) then stop
    sgclose
    

    return, plot_file

end


event_id = '2015_0901_11'
probe = '1'

event_id = '2015_0901_10'
probe = '4'
print, micro_injection_fig_overview_v03(event_id, probe=probe, test=1, update=1)
end
