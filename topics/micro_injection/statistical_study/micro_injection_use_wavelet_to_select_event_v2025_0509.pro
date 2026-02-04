;+
; Check to use wavelet to select micro-injection events.
;-

function find_micro_injection, spec_var, flux_var, $
    significant_flux_level, significant_psd_level, freq_range, plot_tr, gen_plot=gen_plot, test=test

    retval = !null
    
    specs = get_var_data(spec_var, times=times, freqs, settings=settings)
    ntime = n_elements(times)
    nfreq = n_elements(freqs)
    freq_index = where_pro(freqs, '[]', freq_range)
    the_specs = specs[*,freq_index]
    
;---Obtain the significant psd pixels and their times.
    significant_indexs = where(the_specs ge significant_psd_level, nsignificant_index)
    if nsignificant_index eq 0 then return, retval
    index_times = dblarr(nsignificant_index)
    index_freqs = dblarr(nsignificant_index)
    shape = [ntime,nfreq]
    for ii=0,nsignificant_index-1 do begin
        index2d = array_indices(shape, significant_indexs[ii], dimensions=1)
        index_times[ii] = index2d[0]
        index_freqs[ii] = index2d[1]
    endfor
    uniq_index = sort_uniq(index_times)
    uniq_times = times[uniq_index]

    
;---For each time, check if the flux is significant.
    mi_times = list()
    window = 600d   ; sec.
    foreach time, uniq_times, tid do begin
        the_tr = time+[-1,1]*window*0.5
        the_fluxs = get_var_data(flux_var, in=the_tr, times=the_times)
        index = where(the_fluxs ge significant_flux_level, count)
        if count eq 0 then continue
        mi_times.add, time
    endforeach

    
    
;---For each time, check if there is a spectral peak within the freq_range.
    uniq_times = mi_times.toarray()
    mi_times = list()
    foreach time, uniq_times, tid do begin
        the_tr = time+[-1,1]*window*0.5

        the_psds = get_var_data(spec_var, in=the_tr, times=the_times, freqs)
        the_psd = mean(the_psds,dimension=1)
        
        
        peak_indexs = []
        for fid=1,nfreq-2 do begin
            if the_psd[fid] gt the_psd[fid-1] and the_psd[fid] gt the_psd[fid+1] then peak_indexs = [peak_indexs,fid]
        endfor
        npeak_index = n_elements(peak_indexs)
        if npeak_index eq 0 then continue
        
        ; Need to be in the wanted freq range.
        index = where_pro(freqs[peak_indexs], '()', freq_range, count=npeak_index)
        if npeak_index eq 0 then continue
        peak_indexs = peak_indexs[index]
        
        ; Need to be a significant peak.
        the_flags = fltarr(npeak_index)
        rec_width = 4
        significant_ratio = 1.5
        foreach peak_index, peak_indexs, pid do begin
            i0 = peak_index-rec_width
            i1 = peak_index+rec_width
            if i0 le 0 or i1 lt npeak_index then continue
            if the_psd[peak_index] lt significant_psd_level then continue
            if the_psd[peak_index]/the_psd[i0] lt significant_ratio then continue
            if the_psd[peak_index]/the_psd[i1] lt significant_ratio then continue
            the_flags[pid] = 1
        endforeach
        index = where(the_flags eq 1, npeak_index)
        if npeak_index eq 0 then continue
        peak_indexs = peak_indexs[index]
        
        mi_times.add, time
    endforeach
    
    if n_elements(mi_times) eq 0 then return, retval
    
    
    mi_times = mi_times.toarray()
    time_step = sdatarate(times)
    mi_trs = time_to_range(mi_times, time_step=time_step)
    nmi_tr = n_elements(mi_trs[*,0])

    return, mi_trs
    
    
end


function micro_injection_use_wavelet_to_select_event, time_range, mission_probe=mission_probe, errmsg=errmsg, test=test

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
    
    plot_file = 0
    sgopen, plot_file, size=[12,8]
    thick = 2
    
    nplot_var = 3
    margins = [12,4,12,1]
    all_poss = sgcalcpos(nplot_var*nsample_energy, margins=margins)


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
    
    foreach energy_str, search_settings['energy_strs'] do begin
        my_info = search_info[energy_str]
        flux_var = my_info['flux_var']
        spec_var = my_info['spec_var']
        significant_flux_level = my_info['significant_flux_level']
        significant_psd_level = my_info['significant_psd_level']
        options, spec_var, constant=freq_range
        options, flux_var, constant=significant_flux_level
        mi_trs = find_micro_injection(spec_var, flux_var, $
            significant_flux_level, significant_psd_level, freq_range, time_range, test=1)
        
        ; Filter by duration.
        min_mi_dur = search_settings['min_mi_dur']
        nmi_tr = n_elements(mi_trs[*,0])
        if nmi_tr ge 1 then begin
            durs = mi_trs[*,1]-mi_trs[*,0]
            index = where(durs ge min_mi_dur, nmi_tr)
            if nmi_tr gt 0 then begin
                mi_trs = mi_trs[index,*]
            endif else begin
                mi_trs = !null
            endelse
        endif else begin
            mi_trs = !null
        endelse

        my_info['mi_trs'] = mi_trs
    endforeach
    
    
    
    
;---Survey plot.
    routine = mission+'_read_kev_electron_cdaweb'
    en_spec_var = call_function(routine, time_range, probe=probe)
    plot_vars = list(en_spec_var)
    foreach energy_str, search_settings['energy_strs'] do begin
        my_info = search_info[energy_str]
        flux_var = my_info['flux_var']
        spec_var = my_info['spec_var']
        plot_vars.add, [flux_var,spec_var], extract=1
    endforeach
    plot_vars = plot_vars.toarray()
    nplot_var = n_elements(plot_vars)
    
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
    plot_dir = join_path([project_info['plot_dir'],'survey_plot','wavelet_to_select_event',year_str])
    base = 'micro_injection_use_wavelet_to_select_event_survey_plot_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_'+mission+'_'+probe+'_v01.pdf'
    plot_file = join_path([plot_dir,base])
    if keyword_set(test) then plot_file = 0
    fig_size = [12,8]
    sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
    margins = [12,4+nvar_label,10,1]
    poss = sgcalcpos(nplot_var, margins=margins)
    xrange = time_double(time_range)
    xstep = 3600d*2
    tplot_options, 'tickinterval', xstep
    vlab_margin = margins[0]-1
    
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

        thick = 4
        energy_str = settings['energy_str']
        my_info = search_info[energy_str]
        mi_trs = my_info['mi_trs']
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
        
        ; Add significant level.
        cbpos = get_var_setting(spec_var,'zposition')
        zrange = get_var_setting(spec_var,'zrange')
        zlog = 1
        set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
        ty = significant_psd_levels[vid]
        plots, [0,1], ty+[0,0], color=contour_color, thick=thick, data=1
    endforeach
    
    ;if keyword_set(test) then stop
    sgclose
    
;---For each interval.
    psd_fig_size = [6,6]
    psd_margins = [12,4,10,1]
    foreach energy_str, search_settings['energy_strs'] do begin
        my_info = search_info[energy_str]
        flux_var = my_info['flux_var']
        spec_var = my_info['spec_var']
        spec_plot_var = spec_var+'_plot'
        mi_trs = my_info['mi_trs']
        nmi_tr = n_elements(mi_trs[*,0])
        if nmi_tr eq 0 then continue
        significant_psd_level = my_info['significant_psd_level']
        
        specs = get_var_data(spec_var, times=times, freqs, settings=settings)
        for tid=0,nmi_tr-1 do begin
            mi_tr = reform(mi_trs[tid,*])
            base = 'micro_injection_use_wavelet_to_select_event_survey_plot_'+energy_str+'_'+time_string(mi_tr[0],tformat='YYYY_MMDD')+'_'+mission+'_'+probe+'_v01.pdf'
            date_str = time_string(mi_tr[0],tformat='YYYY_MMDD')
            plot_file = join_path([plot_dir,date_str,base])
            if keyword_set(test) then plot_file = 0
            sgopen, plot_file, size=psd_fig_size, xchsz=xchsz, ychsz=ychsz

            index = where_pro(times, '[]', mi_tr, count=count)
            if count eq 0 then message, 'Inconsistency ...'
            plot_vars = [flux_var,spec_plot_var]
            nplot_var = n_elements(plot_vars)
            ypans = [fltarr(nplot_var)+1,1.5]
            ypads = [fltarr(nplot_var-1)+0.4,5]
            poss = sgcalcpos(nplot_var+1, margins=psd_margins, ypans=ypans,ypad=ypads)
            fig_labels = letters(nplot_var+1)+'. '+['e- flux','Wavelet','PSD']
            
            pid = where(plot_vars eq spec_plot_var)
            tpos = poss[*,pid]
            cbpos = tpos
            cbpos[0] = tpos[2]+xchsz*0.8
            cbpos[2] = cbpos[0]+xchsz*0.8
            zticklen = abs_ticklen/(cbpos[2]-cbpos[0])/fig_size[0]
            options, spec_plot_var, zticklen=zticklen, zcharsize=zcharsize, zposition=cbpos
            foreach plot_var, plot_vars, pid do begin
                tpos = poss[*,pid]
                xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
                yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
                options, plot_var, xticklen=xticklen, yticklen=yticklen
            endforeach
            
            xrange = time_double(time_range)
            tplot, plot_vars, position=poss[*,0:nplot_var-1], noerase=1, trange=xrange
            
            pid = where(plot_vars eq spec_plot_var)
            tpos = poss[*,pid]
            specs = get_var_data(spec_plot_var, freqs, times=times, settings=settings)
            levels = [significant_psd_level]
            yrange = settings['yrange']

            contour, specs, times, freqs, position=tpos, $
                xlog=0, xstyle=5, xrange=xrange, $
                ylog=1, ystyle=5, yrange=yrange, levels=levels, noerase=1, color=contour_color

            thick = 4
            energy_str = settings['energy_str']
            plots, mi_tr, yrange[1]+[0,0], data=1, thick=thick, color=sgcolor('red')

            ; Add significant level.
            cbpos = get_var_setting(spec_plot_var,'zposition')
            zrange = get_var_setting(spec_plot_var,'zrange')
            zlog = 1
            set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
            ty = significant_psd_level
            plots, [0,1], ty+[0,0], color=contour_color, thick=thick, data=1
            
            
            tpos = poss[*,nplot_var]
            xxs = freqs
            xtitle = 'Freq (mHz)'
            xrange = get_var_setting(spec_plot_var, 'yrange')
            yys = mean(specs[index,*], dimension=1)
            ytitle = 'PSD of Log flux!C'+energy_str
            yrange = get_var_setting(spec_plot_var, 'zrange')
            xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            plot, xxs, yys, $
                xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
                ystyle=1, ylog=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, $
                position=tpos, noerase=1
            ; Add significant level.
            plots, xrange, significant_psd_level+[0,0], data=1, linestyle=1
            foreach tx, freq_range do begin
                plots, tx*1e3+[0,0], yrange, data=1, linestyle=1
            endforeach
            
            ; Find spectral peaks.
            nfreq = n_elements(xxs)
            peak_indexs = []
            the_psd = yys
            for fid=1,nfreq-2 do begin
                if the_psd[fid] gt the_psd[fid-1] and the_psd[fid] gt the_psd[fid+1] then peak_indexs = [peak_indexs,fid]
            endfor
            npeak_index = n_elements(peak_indexs)
            if npeak_index eq 0 then message, 'Inconsistency ...'

            ; Need to be in the wanted freq range.
            index = where_pro(freqs[peak_indexs], '()', freq_range*1e3, count=npeak_index)
            if npeak_index eq 0 then message, 'Inconsistency ...'
            peak_indexs = peak_indexs[index]

            ; Need to be a significant peak.
            the_flags = fltarr(npeak_index)
            rec_width = 4
            significant_ratio = 1.5
            foreach peak_index, peak_indexs, pid do begin
                i0 = peak_index-rec_width
                i1 = peak_index+rec_width
                if i0 le 0 or i1 lt npeak_index then continue
                if the_psd[peak_index] lt significant_psd_level then continue
                if the_psd[peak_index]/the_psd[i0] lt significant_ratio then continue
                if the_psd[peak_index]/the_psd[i1] lt significant_ratio then continue
                the_flags[pid] = 1
            endforeach
            index = where(the_flags eq 1, npeak_index)
            if npeak_index eq 0 then message, 'Inconsistency ...'
            peak_indexs = peak_indexs[index]

            tmp = max(yys[peak_indexs], index)
            peak_index = peak_indexs[index]
            plots, xxs[peak_index], yys[peak_index], color=sgcolor('red'), data=1, psym=1
            
            sector_color = sgcolor('orange')
            i0 = peak_index-rec_width
            i1 = peak_index+rec_width
            plots, xxs[[i0,i1]], yrange[0], data=1, color=sector_color, thick=thick
            plots, xxs[i0:i1], yys[i0:i1], data=1, color=sector_color
            foreach ii, [i0,i1] do begin
                ty = yys[ii]*significant_ratio
                plots, xxs[[i0,i1]], ty+[0,0], linestyle=2, color=sector_color
            endforeach
            
            freq_index = where_pro(xxs,'[]',freq_range*1e3)
            if i0 lt freq_index[0] then i0 = freq_index[0]
            if i1 gt freq_index[-1] then i1 = freq_index[-1]
            
            peak_freq = freqs[peak_index]
            peak_psd = yys[peak_index]
            peak_power_ratio = total(yys[i0:i1])/total(yys[freq_index])
            msgs = list()
            msgs.add, 'Peak freq (mHz): '+string(peak_freq,format='(F4.1)')
            msgs.add, 'Peak period (min): '+string(1d3/peak_freq/60,format='(F4.1)')
            msgs.add, 'Power ratio : '+string(peak_power_ratio,format='(F5.2)')
            tx = tpos[0]+xchsz*1
            ty0 = tpos[3]-ychsz*1
            foreach msg, msgs, mid do begin
                ty = ty0-ychsz*mid
                xyouts, tx,ty,normal=1, msg
            endforeach
            
            foreach msg, fig_labels, pid do begin
                tpos = poss[*,pid]
                tx = tpos[0]-xchsz*(margins[0]-2)
                ty = tpos[3]-ychsz*0.7
                xyouts, tx,ty,msg, normal=1
            endforeach
            
            if keyword_set(test) then stop
            sgclose
        endfor
    endforeach

stop

end


test = 1
mission_probe = ['mms','4']
time_range = ['2015-09-01','2015-09-02']

;time_range = ['2016-03-05/12:50','2016-03-06/02:40']
;time_range = ['2017-01-12/00:00','2017-01-12/14:00']
;time_range = ['2017-01-01/00:00','2017-01-02/00:00']
;mission_probe = ['mms','4']
print, micro_injection_use_wavelet_to_select_event(time_range, mission_probe=mission_probe, test=test)
end