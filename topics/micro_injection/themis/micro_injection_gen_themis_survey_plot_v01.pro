
function micro_injection_gen_themis_survey_plot_v01, time_range, probe=probe, $
    filename=filename, test=test, get_name=get_name

    project_info = micro_injection_stat_load_project()
    my_name = get_filename()
    version = get_file_version(my_name)
    mission = 'themis'
    plot_root = join_path([diskdir('data'),'survey_plot'])
    plot_root = join_path([plot_root,'microinjection','micro_injection_themis_survey_'+version])

    base = 'micro_injection_themis_survey_plot_'+time_string(time_range[0],tformat='YYYYMMDD')+'_'+mission+'_'+probe+'_'+version+'.pdf'
    year_str = time_string(time_range[0], tformat='YYYY')
    plot_file = join_path([plot_root,'themis'+'_'+probe,year_str,base])
    if keyword_set(filename) then plot_file = filename
    if keyword_set(get_name) then return, plot_file
    if keyword_set(test) then plot_file = 0
    retval = !null
    
;    time_range = time_double('2009-06-19/'+['00:00','24:00'])
;    probe = 'a'

    
;---Load data.
    prefix = 'th'+probe+'_'
    mission_probe = 'th'+probe

    ; B field and ion velocity.
    b_gsm_var = themis_read_bfield(time_range, probe=probe, errmsg=errmsg, id='fgs')
    if errmsg ne '' then return, retval
    u_gsm_var = prefix+'u_gsm'
    del_data, u_gsm_var
    u_gsm_var = themis_read_ion_vel(time_range, probe=probe, errmsg=errmsg, id='peif')
    if errmsg ne '' then return, retval
    u_gsm_var1 = rename_var(u_gsm_var,output=u_gsm_var+'_lowres')
    u_gsm_var = themis_read_ion_vel(time_range, probe=probe, errmsg=errmsg, id='peir')    
    u_gsm_var2 = rename_var(u_gsm_var,output=u_gsm_var+'_highres')
    u1_gsm = var_get_data(u_gsm_var1, times=times1, settings=settings)
    u2_gsm = var_get_data(u_gsm_var2, times=times2)
    index_highres = where(finite(snorm(u2_gsm)),count, complement=index_lowres)
    if count ne 0 then begin
        lowres_trs = times2[time_to_range(index_lowres,time_step=1)]        
        ntr = n_elements(lowres_trs[*,0])
        index_lowres = []
        for tid=0,ntr-1 do begin
            index = where_pro(times1,'[]',lowres_trs[tid,*], count=count)
            if count eq 0 then continue
            index_lowres = [index_lowres,index]
        endfor
        times = [times1[index_lowres],times2[index_highres]]
        u_gsm = [u1_gsm[index_lowres,*],u2_gsm[index_highres,*]]
        index = sort_uniq(times, index=1)
        times = times[index]
        u_gsm = u_gsm[index,*]
    endif else begin
        times = times1
        u_gsm = u1_gsm
    endelse
    u_gsm_var = var_store(u_gsm_var, u_gsm, times, settings=settings)
    

    ; Orbit.
    r_var = themis_read_orbit(time_range, probe=probe)
    options, r_var, mission_probe=mission_probe
    mlat_vars = lets_read_mlat_vars(orbit_var=r_var)
    mlt_var = mlat_vars['mlt']
    dis_var = mlat_vars['dis']
    mlts = var_get_data(mlt_var, times=times)
    index = where(mlts lt 0, count)
    if count ne 0 then begin
        mlts[index] += 24
        mlt_var = var_store(mlt_var, mlts, times)
    endif


    ; SST.
    datatype = 'psef'
    prefix2 = prefix+datatype+'_'
    en_high_var = prefix2+'en_eflux'
    pa_high_var = prefix2+'an_eflux_pa'
    del_data, pa_high_var
    if check_if_update(en_high_var, time_range) then begin
        thm_part_load, data_type=datatype, probe=probe, trange=time_range
        thm_part_getspec, data_type=datatype, probe=probe, trange=time_range, outputs='energy'
        thm_part_getspec, data_type=datatype, probe=probe, trange=time_range, outputs='pa'
        options, en_high_var, requested_time_range=time_range
        options, pa_high_var, requested_time_range=time_range
    endif
    if check_if_update(pa_high_var) then return, retval

    unit = 'eV/cm!E2!N-s-sr-eV'
    zrange = [1e2,1e7]
    ztickv = [1e2,1e3,1e4,1e5,1e6,1e7]
    ztickv_log = alog10(ztickv)
    zticks = n_elements(ztickv)-1
    ztickn = '10!U'+string(ztickv_log,format='(I0)')
    ztickn[0:*:2] = ' '
    options, [en_high_var], color_table=40, no_interp=1, $
        ytitle='Energy!C(eV)', ztitle=unit, $
        zrange=zrange, zstyle=1, zlog=1, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=9, $
        yrange=[3.1e4,7.2e5], ystyle=1, ylog=1, ytickv=[5e4,5e5], ytickname='10!U'+['4','5'], yticks=1, yminor=9
    options, [pa_high_var], color_table=40, no_interp=1, $
        ytitle='PA!C(deg)', ztitle=unit, $
        zrange=zrange, zstyle=1, zlog=1, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=9, $
        yrange=[0,180], ystyle=1, ylog=0, ytickv=[30,90,150], ytickname=['30','90','150'], yticks=2, yminor=6

    
    ; ESA.
    datatype = 'peef'
    prefix2 = prefix+datatype+'_'
    zrange = [1e5,1e8]
    ztickv = [1e5,1e6,1e7,1e8]
    ztickn = '10!U'+['5','6','7','8']
;    ztickn[0:2:*] = ' '
    zticks = n_elements(ztickv)-1
    en_low_var = themis_read_en_spec(time_range, probe=probe, species='e', id='esa_l2')
    ;en_low_var = prefix2+'en_eflux'
    ;u_gsm_var = prefix2+'velocity'
;    if check_if_update(en_low_var, time_range) then begin
;        thm_part_load, data_type=datatype, probe=probe, trange=time_range
;        thm_part_getspec, data_type=datatype, probe=probe, trange=time_range, outputs='energy'
;        thm_part_products, data_type=datatype, probe=probe, trange=time_range, outputs='moments'
;        options, en_high_var, requested_time_range=time_range
;        options, pa_high_var, requested_time_range=time_range
;    endif
    options, en_low_var, color_table=40, no_interp=1, $
        zrange=zrange, zstyle=1, zlog=1, ztickv=ztickv, ztickname=ztickn, zticks=zticks, zminor=9, $
        yrange=[1.1e1,2.6e4], ystyle=1, ylog=1, ytickv=[1e2,1e3,1e4], ytickname='10!U'+['2','3','4'], yticks=2, yminor=9


;---Generate plot.
    tickinterval = 30*60d

    ; Init plot_vars.
    plot_info = orderedhash()
    
    
    plot_info[en_high_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'e- high', $
        'ypan', 0.8, $
        'setting', dictionary( ) )
    plot_info[pa_high_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'e- PA', $
        'ypan', 0.8, $
        'setting', dictionary( ) )
    plot_info[en_low_var] = dictionary($
        'routine', 'plot_spec', $
        'panel_label_text', 'e- low', $
        'setting', dictionary( ) )
    
    plot_info[b_gsm_var] = dictionary($
        'routine', 'plot_linlog', $
        'panel_label_text', 'B GSM', $
        'setting', dictionary($
            'constant', 0, $
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
    
    plot_info[u_gsm_var] = dictionary($
        'routine', 'plot_linlog', $
        'panel_label_text', 'U GSM', $
        'setting', dictionary($
            'constant', 0, $
            'linear_yrange', [-1,1]*50, $
            'log_yrange', [-1,1]*500, $
            'linear_tick_setting', dictionary($
                'yticks', 2, $
                'ytickv', [-1,0,1]*30, $
                'yminor', 3 ), $
            'log_tick_setting', dictionary($
                'yticks', 1, $
                'ytickv', [1,10]*50, $
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
    plot_tr = time_range
    
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


;time_range = time_double(['2008-12-31','2009-01-01'])
;probe = 'a'
;datatype = 'peef'
;thm_part_load, datatype=datatype, probe=probe, trange=time_range
;thm_part_products, datatype=datatype, probe=probe, trange=time_range, outputs='moments'
;
;stop

test = 0
years = make_bins([2008,2020],1)

year_strs = string(years,format='(I4)')
probes = ['d','e','a']
min_dis = 9d
mlt_range = 12+[-1,1]*8d
secofday = constant('secofday')
foreach year_str, year_strs do begin
    year_time_range = time_double(year_str+'/'+['01-01/00:00','12-31/24:00'])
    foreach probe, probes do begin
        print, 'Generating survey plot for THEMIS-'+strupcase(probe)+' in '+year_str+' ...'
        prefix = 'th'+probe+'_'
        mission_probe = 'th'+probe
        
        ; Read orbit var.
        r_var = themis_read_orbit(year_time_range+[-1,1]*secofday, probe=probe)
        options, r_var, mission_probe=mission_probe

        ; Read orbit.
;        dis_var = prefix+'dis'
;        if check_if_update(dis_var, year_time_range) then begin
;            diss = snorm(var_get_data(r_var, times=times))
;            dis_var = var_store(dis_var, diss, times, id='dis')
;            options, dis_var, requested_time_range=year_time_range
;        endif
        mlat_vars = lets_read_mlat_vars(orbit_var=r_var, update=1)
        mlt_var = mlat_vars['mlt']
        dis_var = mlat_vars['dis']
        mlts = var_get_data(mlt_var, times=times)
        index = where(mlts lt 0, count)
        if count ne 0 then begin
            mlts[index] += 24
            mlt_var = var_store(mlt_var, mlts, times)
        endif
        
        ; Find orbit sections.
        diss = var_get_data(dis_var, times=times)
        index = where(diss ge min_dis, count)
        if count eq 0 then continue
        trs = times[time_to_range(index,time_step=1)]
        ntr = n_elements(trs[*,0])
        
        ; Filter in mlt.
        tr_list = list()
        for tid=0,ntr-1 do begin
            tr = trs[tid,*]
            diss = var_get_data(dis_var, in=tr, times=times)
            
            ; Remove partial orbits.
            del_dis = 0.1
            if abs(diss[0]-min_dis) ge del_dis then continue
            if abs(diss[-1]-min_dis) ge del_dis then continue
            
            apogee_dis = max(diss, index)
            apogee_time = times[index]
            apogee_mlt = var_get_data(mlt_var, at=apogee_time)
            ;print, string(apogee_mlt,format='(F5.1)')+' hr'
            if product(apogee_mlt-mlt_range) le 0 then begin
                ;print, 'Inside mlt_range (hr): ['+strjoin(string(mlt_range,format='(I0)'),',')+']'                
                tr_list.add, tr
            endif
        endfor
        
        ; Generate plot.
        foreach tr, tr_list do begin
            plot_file = micro_injection_gen_themis_survey_plot_v01(tr, probe=probe, test=test, get_name=1)
            print, plot_file
            if file_test(plot_file) eq 1 then continue
            plot_file = micro_injection_gen_themis_survey_plot_v01(tr, probe=probe, test=test)
        endforeach
    endforeach
endforeach

end