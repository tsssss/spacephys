;+
; Check to use wavelet to select micro-injection events.
;-



function micro_injection_use_wavelet_to_select_event, time_range, mission_probe=mission_probe, errmsg=errmsg, test=test, gen_psd_plot=gen_psd_plot

    errmsg = ''
    retval = !null
    version = 'v01'

    mission = mission_probe[0]
    probe = mission_probe[1]
    prefix = mission+probe+'_'
    project_info = micro_injection_stat_load_project()
    sample_energys = project_info['sample_energys']
    nsample_energy = n_elements(sample_energys)
    scale_info = project_info['scale_info']
    significant_flux_levels = [1.5,0.3]
    

;---Load data.
    common_times = make_bins(time_double(time_range), 1d)

    log_flux_vars = list()
    routine = 'micro_injection_read_cwt_kev_electron'
    flux_vars = call_function(routine, time_range, mission_probe=mission_probe, errmsg=errmsg, log=1)
    spec_vars = flux_vars+'_mor'
    
    
    
;---Find micro injections.
    freq_range = [1,10.]*1e-3   ; Hz
    freq_range = minmax(1d/([1d,10]*60))
    significant_psd_levels = 1e-1*[1,significant_flux_levels[1]/significant_flux_levels[0]]
    significant_psd_levels = 1e-1*[1,1]


;---Collect info.
    energy_strs = list()
    search_info = orderedhash()
    foreach flux_var, flux_vars, vid do begin
        energy_str = get_var_setting(flux_var, 'labels')
        energy_strs.add, energy_str
        spec_var = spec_vars[vid]
        search_info[energy_str] = dictionary($
            'flux_var', flux_var, $
            'spec_var', spec_var, $
            'significant_flux_level', significant_flux_levels[vid], $
            'significant_psd_level', significant_psd_levels[vid], $
            'mi_trs', !null )
        options, [flux_var,spec_var], energy_str=energy_str
    endforeach
    search_settings = dictionary($
        'min_mi_dur', 600d, $
        'sample_energys', sample_energys, $
        'energy_strs', energy_strs, $
        'freq_range', freq_range )
    
;---Prepare for survey plot.
    routine = mission+'_read_kev_electron_cdaweb'
    en_spec_var = call_function(routine, time_range, probe=probe)
    plot_vars = list(en_spec_var)
    fig_labels = list('e- EN')
    foreach energy_str, search_settings['energy_strs'] do begin
        my_info = search_info[energy_str]
        flux_var = my_info['flux_var']
        spec_var = my_info['spec_var']
        plot_vars.add, [flux_var,spec_var], extract=1
        fig_labels.add, ['e- flux','Wavelet']+'!C    '+energy_str, extract=1
    endforeach
    plot_vars = plot_vars.toarray()
    nplot_var = n_elements(plot_vars)
    fig_labels = letters(nplot_var)+'. '+fig_labels.toarray()

    ; Settings for en_spec_var.
    zrange = [1e-1,1e5]
    log_ztickv = smkarthm(-1,5,1,'dx')
    ztickv = 10d^log_ztickv
    zticks = n_elements(ztickv)-1
    ztickn = '10!U'+string(log_ztickv,format='(I0)')
    ztickn[1:*:2] = ' '
    zminor = 9
    zticklen = -0.5
    yrange = [10,100]*5
    ytickv = yrange
    yticks = n_elements(ytickv)-1
    yminor = 9
    add_setting, en_spec_var, smart=1, dictionary($
        'yrange', yrange, $
        'ytickv', ytickv, $
        'yticks', yticks, $
        'yminor', yminor )
    options, en_spec_var, $
        ytitle='Energy!C(keV)', ysubtitle=' ', $
        zrange=zrange, ztickv=ztickv, zticks=zticks, zminor=zminor, ztickname=ztickn, zticklen=zticklen, $
        constant=sample_energys

    ; Settings for flux_vars.
    yrange = alog10(zrange)
    ystep = 2
    ytickv = smkarthm(yrange[0],yrange[1],ystep,'dx')
    yticks = n_elements(ytickv)-1
    options, flux_vars, yrange=yrange, ytickv=ytickv, yticks=yticks, yminor=ystep

    ; Settings for psd spec_vars.
    ; Change yticks from Hz to mHz.
    spec_plot_vars = spec_vars+'_plot'
    options, spec_plot_vars, zticklen=zticklen
    foreach spec_var, spec_vars, vid do begin
        specs = get_var_data(spec_var, times=times, freqs, settings=settings)
        settings['ytitle'] = 'Freq!C(mHz)'
        settings['yrange'] = settings['yrange']*1e3
        flux_var = flux_vars[vid]
        energy_str = get_var_setting(flux_var, 'labels')
        settings['ztitle'] = 'PSD of Log flux!C'+energy_str
        settings['constant'] = freq_range*1e3
        spec_plot_var = var_store(spec_plot_vars[vid], specs, times, freqs*1e3, settings=settings)
        pid = where(plot_vars eq spec_var)
        plot_vars[pid] = spec_plot_var
    endforeach

    ; Load orbit and related vars for var_labels.
    routine = mission+'_read_orbit'
    r_gsm_var = call_function(routine, time_range, probe=probe, coord='gsm')
    options, r_gsm_var, requested_time_range=time_range
    mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
    foreach var, mlat_vars.values() do print, 'Loading '+var+' ...'

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
    
    
    year_str = time_string(time_range[0],tformat='YYYY')
    plot_dir = join_path([project_info['plot_dir'],'survey_plot','wavelet_to_select_event',year_str])   ; This will fill up Google Drive's space.
    plot_dir = join_path([sdiskdir('data'),'survey_plot','microinjection','wavelet_to_select_event',year_str])
    
    
;---Find microinjections.
    foreach energy_str, search_settings['energy_strs'] do begin
        my_info = search_info[energy_str]
        flux_var = my_info['flux_var']
        spec_var = my_info['spec_var']
        significant_flux_level = my_info['significant_flux_level']
        significant_psd_level = my_info['significant_psd_level']
        options, spec_var, constant=freq_range
        options, flux_var, constant=significant_flux_level
        mi_search = micro_injection_use_wavelet_to_find_event(spec_var, flux_var, $
            significant_flux_level, significant_psd_level, freq_range, time_range, $
            test=test, plot_dir=plot_dir, gen_plot=gen_psd_plot)
        

        mi_times = list()
        my_info['mi_search'] = mi_search
        if n_elements(mi_search) ne 0 then begin
            foreach time, mi_search.keys() do begin
                info = mi_search[time]
                if not info['is_mi'] then continue
                mi_times.add, time
            endforeach
            mi_times = mi_times.toarray()
            data_rate = sdatarate(times)
            mi_trs = time_to_range(mi_times, time_step=data_rate)
            my_info['mi_trs'] = mi_trs
        endif else begin
            my_info['mi_trs'] = !null
        endelse

    endforeach
    
    
    
    
;---Survey plot.
    base = 'micro_injection_use_wavelet_to_select_event_survey_plot_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_'+mission+'_'+probe+'_v01.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0
    fig_size = [12,8]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    margins = [12,4+nvar_label,10,1]
    poss = sgcalcpos(nplot_var, margins=margins)
    xrange = time_double(time_range)
    xstep = 3600d*1
    tplot_options, 'tickinterval', xstep
    vlab_margin = margins[0]-1
    thick = keyword_set(test)? 4: 8
    
    abs_ticklen = -0.3*ychsz*fig_size[1]
    zcharsize = 1.
    foreach plot_var, plot_vars, pid do begin
        tpos = poss[*,pid]
        xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
        yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
        options, plot_var, xticklen=xticklen, yticklen=yticklen
        
        spec = get_var_setting(plot_var, 'spec', exist)
        if exist and spec then begin
            cbpos = tpos
            cbpos[0] = tpos[2]+xchsz*0.8
            cbpos[2] = cbpos[0]+xchsz*0.8
            zticklen = abs_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
            options, plot_var, zticklen=zticklen, zcharsize=zcharsize, zposition=cbpos
        endif
    endforeach
    
    tplot, plot_vars, trange=xrange, noerase=1, position=poss, $
        var_label=var_labels, vlab_margin=vlab_margin
    foreach plot_var, plot_vars, pid do begin
        tpos = poss[*,pid]
        tx = tpos[0]-xchsz*(margins[0]-2)
        ty = tpos[3]-ychsz*0.7
        xyouts, tx,ty,normal=1, fig_labels[pid]
    endforeach
    
    
    ; Add significant flux level.
    cbpos = get_var_setting(en_spec_var,'zposition')
    zrange = get_var_setting(en_spec_var,'zrange')
    zlog = 1
    set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
    sig_flux_color = sgcolor('red')
    foreach ty, significant_flux_levels do begin
        plots, [0,1], 10d^ty+[0,0], color=sig_flux_color, thick=thick, data=1
    endforeach

    ; Add contours and microinjection time ranges.
    contour_color = sgcolor('brown')
    foreach spec_var, spec_plot_vars, vid do begin
        pid = where(plot_vars eq spec_var, count)
        if count eq 0 then continue
        tpos = poss[*,pid]
        specs = get_var_data(spec_var, freqs, times=times, settings=settings)
        levels = [significant_psd_level]
        yrange = settings['yrange']

        contour, specs, times, freqs, position=tpos, $
            xlog=0, xstyle=5, xrange=xrange, $
            ylog=1, ystyle=5, yrange=yrange, levels=levels, noerase=1, color=contour_color

        energy_str = settings['energy_str']
        my_info = search_info[energy_str]
        mi_trs = my_info['mi_trs']
        if n_elements(mi_trs) ne 0 then begin
            nmi_tr = n_elements(mi_trs[*,0])
            for tid=0,nmi_tr-1 do begin
                tr = reform(mi_trs[tid,*])
                if tr[1] le xrange[0] then continue
                if tr[1] ge xrange[1] then continue
                if tr[0] lt xrange[0] then tr[0] = xrange[0]
                if tr[1] gt xrange[1] then tr[1] = xrange[1]
                ;oplot, tr, yrange[1]+[0,0], thick=thick, color=sgcolor('red')
                plots, tr, yrange[1]+[0,0], data=1, thick=thick, color=sgcolor('red')
            endfor
        endif
        
        
        ; Add significant level.
        cbpos = get_var_setting(spec_var,'zposition')
        zrange = get_var_setting(spec_var,'zrange')
        zlog = 1
        set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
        ty = significant_psd_levels[vid]
        plots, [0,1], ty+[0,0], color=contour_color, thick=thick, data=1
    endforeach
    
    if keyword_set(test) then stop
    sgclose

    return, search_info

end


test = 0
missions = ['mms']

search_trs = micro_injection_load_search_time_range()
nsearch_tr = n_elements(search_trs[*,0])
dates = list()
for ii=0,nsearch_tr-1 do begin
    tr = reform(search_trs[ii,*])
    nday = total(tr*[-1,1])/constant('secofday')
    dates.add, findgen(nday)*constant('secofday')+tr[0], extract=1
endfor

foreach date, dates do begin
    if date lt time_double('2015-09-20') then continue
    if date lt time_double('2016-10-01') then continue
    time_range = date+[0,constant('secofday')]
    foreach mission, missions do begin
        probes = call_function(mission+'_probes')
probes = ['4']
        foreach probe, probes do begin
            mission_probe = [mission,probe]
            search_info = micro_injection_use_wavelet_to_select_event(time_range, mission_probe=mission_probe, test=test, gen_psd_plot=1)
        endforeach
    endforeach
endforeach

stop

mission_probe = ['mms','4']
time_range = ['2015-09-01','2015-09-02']

time_range = ['2016-03-05/12:50','2016-03-06/02:40']
time_range = ['2016-03-05','2016-03-06']
time_range = ['2017-01-12/00:00','2017-01-12/14:00']
;time_range = ['2017-01-01/00:00','2017-01-02/00:00']
;mission_probe = ['mms','4']
print, micro_injection_use_wavelet_to_select_event(time_range, mission_probe=mission_probe, test=test)
end