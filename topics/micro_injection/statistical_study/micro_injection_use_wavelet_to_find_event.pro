;+
; spec_var.
; flux_var.
; significant_flux_level.
; significant_psd_level.
; freq_range.
; plot_tr.
; gen_plot=.
; test=.
; plot_dir=.
;-

function micro_injection_use_wavelet_to_find_event, spec_var, flux_var, $
    significant_flux_level, significant_psd_level, freq_range, plot_tr, gen_plot=gen_plot, test=test, plot_dir=plot_dir

    retval = !null

;---Settings.
    window = 600d   ; sec. Look window in time around a given time.
    rec_width = 3   ; #. Look width in freq for spectral peak.
    psd_significant_ratio = 1.5d ; Ratio of the peak to the boundary of the look width.


    ; The wanted freq range.
    specs = get_var_data(spec_var, times=times, freqs, settings=settings)
    freq_index = where_pro(freqs, '[]', freq_range, count=count)
    if count eq 0 then message, 'Inconsistency ...'
    dfreq = freqs[1]/freqs[0]
    if dfreq lt 1 then dfreq = 1d/dfreq         ; ensure >1.
    freq_window_ratio = dfreq^rec_width
    actual_fr = [freq_range[0]*freq_window_ratio,freq_range[1]/freq_window_ratio]
    
    
;---Obtain the significant psd pixels and their times.
    index = where_pro(freqs, '[]', freq_range, count=ntarget_freq)
    if ntarget_freq eq 0 then return, retval
    target_specs = specs[*,index]
    target_freqs = freqs[index]
    ntime = n_elements(times)

    significant_indexs = where(target_specs ge significant_psd_level, nsignificant_index)
    if nsignificant_index eq 0 then return, retval
    index_times = dblarr(nsignificant_index)
    index_freqs = dblarr(nsignificant_index)
    shape = [ntime,ntarget_freq]
    for ii=0,nsignificant_index-1 do begin
        index2d = array_indices(shape, significant_indexs[ii], dimensions=1)
        index_times[ii] = index2d[0]
        index_freqs[ii] = index2d[1]
    endfor
    uniq_index = sort_uniq(index_times)
    uniq_times = times[uniq_index]

;---For each time, check if the flux is significant.
    mi_times = list()
    foreach time, uniq_times, tid do begin
        the_tr = time+[-1,1]*window*0.5
        the_fluxs = get_var_data(flux_var, in=the_tr, times=the_times)
        index = where(the_fluxs ge significant_flux_level, count)
        if count eq 0 then continue
        mi_times.add, time
    endforeach
    uniq_times = mi_times.toarray()
    nuniq_time = n_elements(uniq_times)

    spec_plot_var = spec_var+'_plot'
    freqs_mhz = get_var_value(spec_plot_var)

;---For each time, check if there is a spectral peak within the freq_range.
    search_info = orderedhash()
    data_rate = sdatarate(times)
    symsize = 0.2
    foreach time, uniq_times, tid do begin
;if time le time_double('2015-09-01/14:30:40') then continue
;if time le time_double('2015-09-01/15:30:40') then continue
;if time lt time_double('2015-09-01/16:27:20') then continue
        percentage = 100d*tid/nuniq_time
        msg = 'Processing '+time_string(time)+'...'+string(percentage,format='(F5.1)')+'%'
        lprmsg, msg

        the_tr = time+[-1,1]*data_rate
        the_psds = get_var_data(spec_var, in=the_tr, times=the_times, freqs)
        psds = mean(the_psds,dimension=1)
        msgs = list()
        msgs.add, time_string(time,tformat='YYYY-MM-DD/hh:mm:ss')+' UT'

    ;---Prepare for plot.
        energy_str = get_var_setting(flux_var,'energy_str')
        mission_info = resolve_probe((strsplit(flux_var,'_',extract=1))[0])
        mission = mission_info['mission']
        probe = mission_info['probe']
        if keyword_set(gen_plot) then begin
            base = 'micro_injection_use_wavelet_to_select_event_survey_plot_'+energy_str+'_'+time_string(time,tformat='YYYY_MMDD_hhmm_ss')+'_'+mission+'_'+probe+'_v01.pdf'
            date_str = time_string(time,tformat='YYYY_MMDD')
            plot_file = join_path([plot_dir,date_str,base])
            fig_size = [6d,6]
            margins = [12,4,10,1]
            if keyword_set(test) then plot_file = 0
            sgopen, plot_file, size=fig_size, xchsz=xchsz, ychsz=ychsz
            abs_ticklen = -0.3*ychsz*fig_size[1]

            
            plot_vars = [flux_var,spec_plot_var]
            nplot_var = n_elements(plot_vars)
            ypans = [fltarr(nplot_var)+1,1.5]
            ypads = [fltarr(nplot_var-1)+0.4,5]
            poss = sgcalcpos(nplot_var+1, margins=margins, ypans=ypans,ypad=ypads)
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
            
            time_range = plot_tr
            xrange = time_double(time_range)
            tplot_options, 'tickinterval', 3600d*2
            tplot, plot_vars, position=poss[*,0:nplot_var-1], noerase=1, trange=xrange
            
            pid = where(plot_vars eq spec_plot_var)
            tpos = poss[*,pid]
            ;specs = get_var_data(spec_plot_var, freqs_mhz, times=times, settings=settings)
            settings = get_var_setting(spec_plot_var)
            levels = [significant_psd_level]
            yrange = settings['yrange']

            contour_color = sgcolor('brown')
            contour, specs, times, freqs_mhz, position=tpos, $
                xlog=0, xstyle=5, xrange=xrange, $
                ylog=1, ystyle=5, yrange=yrange, levels=levels, noerase=1, color=contour_color
            plots, time+[0,0], yrange, data=1, color=sgcolor('red')

            ; Add significant level.
            cbpos = get_var_setting(spec_plot_var,'zposition')
            zrange = get_var_setting(spec_plot_var,'zrange')
            zlog = 1
            set_axis, position=cbpos, yrange=zrange, ylog=zlog, xrange=[0,1]
            ty = significant_psd_level
            plots, [0,1], ty+[0,0], color=contour_color, thick=thick, data=1
            
            
            tpos = poss[*,nplot_var]
            xxs = freqs_mhz
            xtitle = 'Freq (mHz)'
            xrange = get_var_setting(spec_plot_var, 'yrange')
            yys = psds
            ytitle = 'PSD of Log flux!C'+energy_str
            yrange = get_var_setting(spec_plot_var, 'zrange')
            xticklen = abs_ticklen/(tpos[3]-tpos[1])/fig_size[1]
            yticklen = abs_ticklen/(tpos[2]-tpos[0])/fig_size[0]
            plot, xxs, yys, $
                xstyle=1, xlog=1, xrange=xrange, xtitle=xtitle, xticklen=xticklen, $
                ystyle=1, ylog=1, yrange=yrange, ytitle=ytitle, yticklen=yticklen, $
                position=tpos, noerase=1
            oplot, xxs, yys, psym=1, symsize=symsize
            ; Add significant level.
            plots, xrange, significant_psd_level+[0,0], data=1, linestyle=1
            foreach tx, freq_range do begin
                plots, tx*1e3+[0,0], yrange, data=1, linestyle=1
            endforeach
            foreach tx, actual_fr do begin
                plots, tx*1e3+[0,0], yrange, data=1, linestyle=2
            endforeach

            
            foreach msg, fig_labels, pid do begin
                tpos = poss[*,pid]
                tx = tpos[0]-xchsz*(margins[0]-2)
                ty = tpos[3]-ychsz*0.7
                xyouts, tx,ty,msg, normal=1
            endforeach
        endif


    ;---Do work. 
        ; Check if the psd is significant.
        sig_psd_index = where(psds ge significant_psd_level, count)
        if count eq 0 then begin
            msgs.add, 'No significant psd found'
            search_info[time] = dictionary($
                'is_mi', 0, $
                'msgs', msgs)
            
            if keyword_set(gen_plot) then begin
                tpos = poss[*,nplot_var]
                tx = tpos[0]+xchsz*1
                ty0 = tpos[3]-ychsz*1
                foreach msg, msgs, mid do begin
                    ty = ty0-ychsz*mid
                    xyouts, tx,ty,normal=1, msg
                endforeach
                if keyword_set(test) then stop
                sgclose
            endif
            continue
        endif

        ; Check if significant psds are within the wanted freq range.
        sig_freq_index = where_pro(freqs[sig_psd_index], '[]', actual_fr, count=count)
        if count eq 0 then begin
            msgs.add, 'No significant psd in the wanted freq range'
            search_info[time] = dictionary($
                'is_mi', 0, $
                'msgs', msgs)
            
            if keyword_set(gen_plot) then begin
                tpos = poss[*,nplot_var]
                tx = tpos[0]+xchsz*1
                ty0 = tpos[3]-ychsz*1
                foreach msg, msgs, mid do begin
                    ty = ty0-ychsz*mid
                    xyouts, tx,ty,normal=1, msg
                endforeach
                if keyword_set(test) then stop
                sgclose
            endif
            continue
        endif
        psd_index = sig_psd_index[sig_freq_index]


        ; Obtain the peak psd and peak freq.
        peak_psd = max(psds[psd_index], index)
        peak_index = psd_index[index]
        peak_freq = freqs[peak_index]
        tmp_index = where_pro(freqs,'[]',peak_freq*[1d/freq_window_ratio,freq_window_ratio])
        i0 = peak_index-rec_width
        i1 = peak_index+rec_width
        if keyword_set(gen_plot) then begin
            tmp_index = smkarthm(i0,i1,1,'dx')
            plots, freqs[tmp_index]*1e3, psds[tmp_index], color=sgcolor('orange'), psym=-1, symsize=symsize
            plots, peak_freq*1e3, peak_psd, psym=6, color=sgcolor('red'), symsize=symsize*3
            plots, peak_freq*1e3+[0,0], yrange, linestyle=1, color=sgcolor('red')
        endif
        msgs.add, 'Peak freq (mHz): '+string(peak_freq*1e3,format='(F4.1)')
        msgs.add, 'Peak period (min): '+string(1/peak_freq/60,format='(F4.1)')
        psd_ratio = total(psds[i0:i1])/total(psds[freq_index])
        msgs.add, 'PSD ratio: '+string(psd_ratio,format='(F6.2)')
        
        is_significant_peak = 1
        range_list = list()
        range_list.add, [i0,peak_index]
        range_list.add, [peak_index,i1]
        foreach range, range_list, rid do begin
            if rid eq 1 then begin
                ti0 = peak_index
                ti1 = i1
                tcolor = sgcolor('blue')
            endif else begin
                ti0 = i0
                ti1 = peak_index
                tcolor = sgcolor('green')
            endelse
            min_psd = min(psds[ti0:ti1], index)
            min_freq = (freqs[ti0:ti1])[index]
            if keyword_set(gen_plot) then begin
                plots, min_freq*1e3, min_psd, color=tcolor, psym=6, symsize=symsize*2
                plots, freqs[[ti0,ti1]]*1e3, min_psd*psd_significant_ratio+[0,0], color=tcolor, linestyle=0
            endif
            if peak_psd/min_psd lt psd_significant_ratio then begin
                is_significant_peak = 0
                break
            endif
        endforeach
        if not is_significant_peak then begin
            msgs.add, 'Not a significant spectral peak'
            search_info[time] = dictionary($
                'is_mi', 0, $
                'msgs', msgs)

            if keyword_set(gen_plot) then begin
                tpos = poss[*,nplot_var]
                tx = tpos[0]+xchsz*1
                ty0 = tpos[3]-ychsz*1
                foreach msg, msgs, mid do begin
                    ty = ty0-ychsz*mid
                    xyouts, tx,ty,normal=1, msg
                endforeach
                if keyword_set(test) then stop
                sgclose
            endif
            continue
        endif
        
        ; Found a significant peak.
        msgs.add, 'Is a significant spectral peak'
        search_info[time] = dictionary($
            'is_mi', 1, $
            'freq', peak_freq, $
            'psd', peak_psd, $
            'psd_ratio', psd_ratio, $
            'msgs', msgs)

        if keyword_set(gen_plot) then begin
            tpos = poss[*,nplot_var]
            tx = tpos[0]+xchsz*1
            ty0 = tpos[3]-ychsz*1
            foreach msg, msgs, mid do begin
                ty = ty0-ychsz*mid
                xyouts, tx,ty,normal=1, msg
            endforeach
            if keyword_set(test) then stop
            sgclose
        endif
    ;---Done.
    endforeach
    
    return, search_info
    
    
end