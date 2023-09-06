;+
; Use the "moving cross correlation" to calculate the time lag. Specify the time range of the first spacecraft,
; cross correlate it to the next one, determine the time range of the second
; spacecraft, and so on.
;-

pro azim_prop_calc_time_lag, project

;---Make sure project has the latest event info, and all data are loaded.
    if n_elements(project) eq 0 then project = azim_prop_load_project()
    azim_prop_update_event_info, project
    events = project.events
    foreach event_id, events.keys() do begin
        event = events[event_id]
        if ~event.haskey('file') then azim_prop_load_data, project=project, event_id=event_id
    endforeach

    ;test_id = '2016_1025_08'
    ;test_id = '2017_0301_15'
    test_id = '2014_0828_10'
    test = 1

    str_pm = '!9'+string(177b)+'!X'
    event_ids = (events.keys()).sort()
    foreach event_id, event_ids do begin
        lprmsg, 'Processing '+event_id+' ...'
        if keyword_set(test_id) then if event_id ne test_id then continue

        event_info = events[event_id]
        cc_time_range = events[event_id].cc_time_range
        cc_time_lag = events[event_id].cc_time_lag
        ref_time = events[event_id].ref_time


    ;---Load data for the current event.
        azim_prop_load_data, project=project, event_id=event_id
        events = project.events
        data_file = events[event_id].file
        if file_test(data_file) eq 0 then message, 'Something is wrong when loading the data ...'


    ;---Sort probes by MLT.
        probes = events[event_id].probes
        event_time_range = events[event_id].time_range
        mean_time = mean(event_time_range)
        mean_mlts = list()
        foreach probe, probes do mean_mlts.add, get_var_data(probe+'_mlt', at=ref_time)
        mean_mlts = mean_mlts.toarray()
        sorted_probes = probes[sort(abs(mean_mlts))]    ; closest to the midnight comes first.


    ;---Unify the time rate for the dB tilt angle.
        vars = sorted_probes+'_db_tilt'
        tilt_data_rate = 5.
        project.common_data_rate = tilt_data_rate
        get_data, vars[0], times
        times = make_bins(minmax(times), tilt_data_rate)
        foreach var, vars do interp_time, var, times

    ;---Running cross correlation.
        cc_info = hash()
        shift_limit = cc_time_lag*60.
        shifts = smkarthm(-shift_limit*0.2,shift_limit*0.8,tilt_data_rate,'dx')
        nshift = n_elements(shifts)
        lags = round(shifts/tilt_data_rate)
        foreach probe, sorted_probes, ii do begin
            cc_info[probe] = dictionary()
            if ii eq 0 then begin
                cc_info[probe].time_shift = 0.
                cc_info[probe].max_corr = 0.
                running_time_range = cc_time_range
                cc_info[probe].time_shift_error = 0.
            endif else begin
                running_time_range += cc_info[sorted_probes[ii-1]].time_shift

                get_data, sorted_probes[ii-1]+'_db_tilt_sector', uts, f0s
                x0s = round((uts-event_time_range[0])/tilt_data_rate)
                nx0 = n_elements(x0s)
                get_data, sorted_probes[ii]+'_db_tilt', uts, f1s
                uts = uts[index]
                f1s = f1s[index]
                x1s = round((uts-event_time_range[0])/tilt_data_rate)
                nx1 = n_elements(x1s)
                corr_2d = fltarr(nshift)
                    ; Calculate the new time.
                    for kk=0, nshift-1 do begin
                        index = x0s[0]+lags[kk]
                        if index lt 0 then continue
                        if index+nx0 ge nx1 then continue
                        tf1 = f1s[index:index+nx0-1]
                        corr_2d[kk] = c_correlate(tf1,f0s,0)
                    endfor
                max_corr = max(corr_2d)
                index = where(corr_2d eq max_corr)
                the_shift = shifts[index[0]]

                ; Apply shift.
                get_data, sorted_probes[ii-1]+'_db_tilt_sector', uts, f0s, limits=lim
                uts += the_shift
                cc_info[probe].time_shift = the_shift
                cc_info[probe].max_corr = max_corr
                store_data, sorted_probes[ii-1]+'_db_tilt_sector', uts, f0s, limits=lim

                ; Error in the time lag (The cluster book, Equation 1.7)
                del_x = x1s-mean(x1s)
                dx_dt = deriv(x1s)/tilt_data_rate
                dt = sqrt(1./(nx1-1)*(1-max_corr)/max_corr*2*mean(del_x^2)/mean(dx_dt^2))
                ;dt = 1./(nx0-1)*(1-max_corr)/max_corr*2*mean((x0s-mean(x0s))^2)/mean((deriv(x0s)/tilt_data_rate)^2)
                cc_info[probe].time_shift_error = dt

                file = join_path([project.plot_dir,'cross_corr','fig_'+event_id+'_calc_cross_corr_'+strjoin(sorted_probes[ii-1:ii],'_')+'.pdf'])
                if keyword_set(test) then file = test+ii-1
                sgopen, file, xsize=5, ysize=6, /inch
                poss = sgcalcpos(4, xchsz=xchsz, ychsz=ychsz, tmargin=4, bmargin=10, rmargin=8)

                options, sorted_probes[ii-1]+'_db_tilt', 'labels', strupcase(sorted_probes[ii-1])
                options, sorted_probes[ii-1]+'_db_tilt_sector', 'labels', strupcase(sorted_probes[ii-1])+'!C  shifted'
                options, sorted_probes[ii]+'_db_tilt', 'labels', strupcase(sorted_probes[ii])

                tpos = poss[*,0:2]
                vars = [sorted_probes[ii-1]+'_db_tilt'+['','_sector'], sorted_probes[ii]+'_db_tilt']
                options, vars, 'xticklen', -0.02
                tplot, vars, /novtitle, $
                    trange=event_time_range, position=tpos

                ; Add labels.
                tpos = poss[*,0]
                tx = tpos[0]
                ty = tpos[3]+ychsz*0.5
                xyouts, tx,ty, /normal, event_id+', cross.corr.'
                tpos = poss[*,1]
                tx = tpos[0]+xchsz*1
                ty = tpos[3]-ychsz*1
                xyouts, tx,ty, /normal, 'dt = '+string(the_shift,format='(I0)')+' sec'

                tpos = poss[*,0]
                plot, event_time_range, [0,1], /nodata, /noerase, xstyle=5, ystyle=5, position=tpos
                tmp = convert_coord(running_time_range[0],0, /data, /to_normal)
                tx1 = tmp[0]
                tmp = convert_coord(running_time_range[1],0, /data, /to_normal)
                tx2 = tmp[0]
                ty = tpos[1]+ychsz*0.5
                plots, [tx1,tx2], [0,0]+ty, /normal
                plots, tx1+[0,0], ty+[-1,1]*ychsz*0.2, /normal
                plots, tx2+[0,0], ty+[-1,1]*ychsz*0.2, /normal
                xyouts, (tx1+tx2)*0.5, ty+ychsz*0.2, /normal, 'data for c.c.', alignment=0.5


                ; Add cc plot.
                tpos = poss[*,3]-[0,1,0,1]*ychsz*5
                xtitle = 'dt (sec)'
                xrange = minmax(shifts)
                ytitle = 'c.c. (#)'
                yrange = [0,1]
                plot, shifts, corr_2d, position=tpos, $
                    xstyle=1, xtitle=xtitle, xrange=xrange, xticklen=-0.02, $
                    ystyle=1, ytitle=ytitle, yrange=yrange, yticklen=-0.01, $
                    /noerase
                plots, the_shift+[0,0], yrange, linestyle=1
                plots, the_shift+[-1,1]*dt, mean(yrange)+[0,0]
                foreach tx, [-1,1] do plots, the_shift+tx*dt, mean(yrange)+[-1,1]*0.03
                tx = tpos[0]
                ty = tpos[3]+ychsz*0.5
                xyouts, tx,ty,/normal, $
                    'max c.c. = '+string(max_corr,format='(F5.2)')+', dt = '+string(the_shift,format='(I0)')+$
                    str_pm+string(dt,format='(I0)')+' sec'

                if keyword_set(test) then stop
                sgclose
            endelse

            var = probe+'_db_tilt'
            get_data, var, uts, dat, limits=lim
            index = where_pro(times, running_time_range+cc_info[probe].time_shift)
            store_data, var+'_sector', uts[index], dat[index], limits=lim
            cc_info[probe].var = probe+'_db_tilt_sector'
        endforeach


    ;---Save results.
        time_lags = fltarr(n_elements(sorted_probes))
        foreach probe, sorted_probes, ii do begin
            time_lags[ii:*] -= cc_info[probe].time_shift
            get_data, probe+'_db_tilt', uts, dat, limits=lim
            store_data, probe+'_db_tilt_new', uts+time_lags[ii], dat, limits=lim
        endforeach

        ; Get the MLT, R_SM, and R_GSM at the reforence time.
        ref_time0 = project.events[event_id].ref_time
        foreach probe, sorted_probes, ii do begin
            if ~project.events[event_id].haskey(probe) then (project.events[event_id])[probe] = dictionary()
            (project.events[event_id])[probe].time_lag = time_lags[ii]
            ref_time = ref_time0-time_lags[ii]
            (project.events[event_id])[probe].ref_time = ref_time
            (project.events[event_id])[probe].max_corr = cc_info[probe].max_corr
            (project.events[event_id])[probe].time_lag_error = cc_info[probe].time_shift_error
            (project.events[event_id])[probe].ref_mlt = get_var_data(probe+'_mlt', at=ref_time)
            rgsm = get_var_data(probe+'_r_gsm', at=ref_time)
            (project.events[event_id])[probe].ref_rsm = cotran(rgsm, ref_time, 'gsm2sm')
        endforeach


        ; Calculate omega_azim.
        nprobe = n_elements(sorted_probes)
        ref_times = dblarr(nprobe)
        ref_mlts = dblarr(nprobe)
        foreach probe, sorted_probes, ii do begin
            ref_times[ii] = (events[event_id])[probe].ref_time
            ref_mlts[ii] = (events[event_id])[probe].ref_mlt
        endforeach
        fit_result = linfit(ref_times, ref_mlts, yfit=yfit, sigma=sigma)
        omega_azim = fit_result[1]*15*60    ; in deg/min.
        project.events[event_id].omega_azim = omega_azim
        rsqr = 1-total((ref_mlts-yfit)^2)/total((ref_mlts-mean(ref_mlts))^2)
        project.events[event_id].rsqr = rsqr
        project.events[event_id].domega_azim = sigma[1]*15*60   ; in deg/min, 1 sigma error.

        ; Sort probes by MLT at ref_time.
        index = sort(abs(ref_mlts))
        sorted_probes = sorted_probes[index]
        project.events[event_id].sorted_probes = sorted_probes

        file = join_path([project.plot_dir,'cross_corr','fig_'+event_id+'_after_time_lag_applied.pdf'])
        if keyword_set(test) then file = 0
        sgopen, file, xsize=5, ysize=6
        tplot, sorted_probes+'_db_tilt_new', trange=event_time_range
        timebar, ref_time0, linestyle=1
        if keyword_set(test) then stop
        sgclose

        file = join_path([project.plot_dir,'cross_corr','fig_'+event_id+'_before_time_lag_applied.pdf'])
        if keyword_set(test) then file = 0
        sgopen, file, xsize=5, ysize=6
        tplot, sorted_probes+'_db_tilt', trange=event_time_range
        timebar, ref_time0, linestyle=1
        if keyword_set(test) then stop
        sgclose

    endforeach

    project.done_calc_time_lag = 1
    if ~keyword_set(test) then azim_prop_update_project, project

end
