;+
; Combine candidates from 2012-2017 and 2007-2009.
;
; For each candidate, we have a time_range, event_id, available_probes, midn_flag
; We check for if dipolarizations exist.
;-

pro azim_df_analysis_candidate, project=project

    if n_elements(project) eq 0 then project = azim_df_load_project()
    if ~project.haskey('done_search_candidate') then azim_df_search_candidate, project=project
    if ~project.haskey('done_search_candidate_themis') then azim_df_search_candidate_themis, project=project

reset = 0

;---Combine candidates from the 2 searches.
    candidates = hash()
    file = join_path([project.data_dir,project.candidate_suffix])
    lines = read_all_lines(file)
    all_probes = ['rbsp'+letters('b'),'g'+['13','14','15'],'th'+letters('e'),'mms1']
    foreach line, lines do begin
        info = strsplit(line,' ',/extract)
        time_range = time_double(info[[0,2]])
        event_id = time_string(mean(time_range), tformat='YYYY_MMDD_hh')
        candidates[event_id] = dictionary($
            'event_id', event_id, $
            'time_range', time_range, $
            'all_probes', all_probes)
    endforeach
;    file = join_path([project.data_dir,project.candidate_suffix_themis])
;    lines = read_all_lines(file)
;    all_probes = 'th'+letters('e')
;    foreach line, lines do begin
;        info = strsplit(line,' ',/extract)
;        time_range = time_double(info[[0,2]])
;        event_id = time_string(mean(time_range), tformat='YYYY_MMDD_hh')
;        candidates[event_id] = dictionary($
;            'event_id', event_id, $
;            'time_range', time_range, $
;            'all_probes', all_probes)
;    endforeach


;---Make sure data exists.
    data_file_flags = list()
    foreach candidate, candidates do begin
        time_range = candidate.time_range
        event_id = candidate.event_id
        data_file_suffix = event_id+'_candidate_data.tplot'
        candidate['data_file_suffix'] = data_file_suffix
        data_file = join_path([project.data_dir,data_file_suffix])
        if keyword_set(reset) then file_delete, data_file, /allow_nonexistent
        if file_test(data_file) eq 0 then begin
            azim_df_load_basic_data, time_range, event_id=event_id, project=project, data_file=data_file
            the_flag = file_test(data_file)
        endif else the_flag = 1
        data_file_flags.add, the_flag
    endforeach
    index = data_file_flags.where(0, count=count)
    if count ne 0 then begin
        event_ids = candidates.keys()
        event_ids = event_ids[index]
        foreach event_id, event_ids do lprmsg, event_id+' does not have data ...'
    endif


;---Identify dipolarizations.
    foreach candidate, candidates do begin
        time_range = candidate.time_range
        event_id = candidate.event_id
;if event_id ne '2014_0828_10' then continue
;if event_id ne '2016_1013_12' then continue
        
        data_file_suffix = candidate['data_file_suffix']
        data_file = join_path([project.data_dir,data_file_suffix])
        del_data, '*'
        tplot_restore, filename=data_file

        ; All available probes.
        all_available_probes = tnames('*_theta')
        foreach probe, all_available_probes, ii do all_available_probes[ii] = strmid(probe,0,strpos(probe,'_'))

        ; Side of the midnight.
        post_midn_triad_flag = get_var_data('post_midn_triad_flag')
        pre_midn_triad_flag = get_var_data('pre_midn_triad_flag')
        midn_flag = (mean(pre_midn_triad_flag) gt mean(post_midn_triad_flag))? 'pre': 'post'

        ; Remove data out of region of interest.
        mlt_limit = project.apogee_filter_setting.mlt_around_midnight
        mlt_limit = 9.
        mlt_range = (midn_flag eq 'pre')? [-mlt_limit,0]: [0,mlt_limit]
        dis_range = [4.,20]
        foreach probe, all_available_probes do begin
            prefix = probe+'_'
            the_var = prefix+'r_sm'
            get_data, the_var, times, rsm
            ntime = n_elements(times)
            flags = bytarr(ntime)+1
; if probe eq 'g13' then stop
            ; In magnetopause.
            rgsm = cotran(rsm, times, 'sm2gsm')
            index = where(check_if_in_magn(rgsm) eq 0, count)
            if count ne 0 then flags[index] = 0

            ; Within MLT range.
            index = lazy_where(get_var_data(prefix+'mlt'),'][', mlt_range, count=count)
            if count ne 0 then flags[index] = 0

            ; Within dis range.
            index = lazy_where(snorm(rsm), ']', min(dis_range), count=count)
            if count ne 0 then flags[index] = 0

            ; Apply to data.
            theta = get_var_data(prefix+'theta')
            index = where(flags eq 0, count)
            if count ne 0 then theta[index] = !values.f_nan
            the_var = prefix+'theta_good'
            store_data, the_var, times, theta
            add_setting, the_var, /smart, {$
                labels: strupcase(probe), $
                constant: 0, $
                ytitle: '(deg)'}
            options, prefix+'theta', 'constant', 0
        endforeach

        available_probes = list()
        min_good_ntime = ntime*0.5
        foreach probe, all_available_probes, ii do begin
            prefix = probe+'_'
            the_var = prefix+'theta_good'
            get_data, the_var, times, data
            index = where(finite(data,/nan), count)
            if count lt min_good_ntime then available_probes.add, probe
        endforeach
        if n_elements(available_probes) eq 0 then begin
            lprmsg, 'No available probes ...'
            tplot, all_available_probes+'_theta', trange=time_range
            stop
            continue
        endif
        available_probes = available_probes.toarray()
        mlts = list()
        mid_time = mean(time_range)
        foreach probe, available_probes do mlts.add, get_var_data(probe+'_mlt', at=mid_time)
        mlts = mlts.toarray()
        sorted_probes = available_probes[sort(abs(mlts))]

        tplot, sorted_probes+'_theta_good', trange=time_range
        print, available_probes
        stop
        
;        time_step = project.time_step
;        smooth_width = 10*60./time_step
;        foreach probe, sorted_probes do begin
;            prefix = probe+'_'
;            get_data, prefix+'theta_good', times, data
;            int_data = data
;            index = where(finite(int_data,/nan), count)
;            if count ne 0 then int_data[index] = 0
;            for ii=0, ntime-1 do begin
;                if ii eq 0 then int_data[ii] = 0
;                int_data[ii] = int_data[ii-1]+data[ii]
;            endfor
;            int_data = smooth(int_data, smooth_width, /edge_truncate)
;            store_data, prefix+'theta_good_int', times, int_data
;            
;            ; Find local minimum.
;            min_locations = list()
;            for ii=1, ntime-2 do begin
;                if int_data[ii] lt int_data[ii-1] and int_data[ii] lt int_data[ii+1] then min_locations.add, times[ii]
;            endfor
;            count = n_elements(min_locations)
;            if count ne 0 then begin
;                min_locations = min_locations.toarray()
;                tplot, prefix+['theta_good_int','theta','theta_good'], trange=time_range
;                timebar, min_locations, color=sgcolor('red')
;            endif
;            
;            ; Need to be isolated: 10 min.
;            ; Need to be large amplitude > 20 deg.
;            stop
;            
;        endforeach


    endforeach

end
