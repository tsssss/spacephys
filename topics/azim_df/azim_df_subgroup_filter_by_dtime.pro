;+
; Remove subgroups that contain large separations.
;-

function azim_df_subgroup_filter_by_dtime, subgroups, project=project, excluded=excluded, test_time=test_time


;test_time = time_double('2014-09-03/11:30')

    tab = constant('4space')
    retval = !null
    max_probe_length = 5

    ; No time lag should be larger than this.
    min_dt1 = 30.*60    ;sec.
    ; No DF on edge should be larger than this.
    min_dt2 = 20.*60    ; sec.

    bad_flags = bytarr(subgroups.length)
    foreach subgroup, subgroups, subgroup_id do begin
        if keyword_set(test_time) then begin
            if product(subgroup.time_range-test_time) gt 0 then continue
            stop
        endif
        lprmsg, 'Check dtime for event: '+time_string(subgroup.time_range)+' ...'
        
    ;---Gather info of the current subgroup.
        probes = subgroup.probes
        nprobe = n_elements(probes)
        df_observed_times = dblarr(nprobe)
        df_list = subgroup.df_list
        foreach df, df_list, ii do df_observed_times[ii] = df.arrival_time

    ;---Filter by dt1.
        df_times = df_observed_times[sort(df_observed_times)]
        dtimes = df_times[1:nprobe-1]-df_times[0:nprobe-2]
        index = where(dtimes ge min_dt1, count)
        if count ne 0 then begin
            lprmsg, tab+'dT too large ...'
            bad_flags[subgroup_id] = 1
            continue
        endif

    ;---Filter by dt2.
        dtimes = dtimes[[0,nprobe-2]]
        index = where(dtimes ge min_dt2, count)
        if count ne 0 then begin
            lprmsg, tab+'dT2 too large ...'
            bad_flags[subgroup_id] = 1
            continue
        endif
    endforeach

    index = where(bad_flags eq 1, nsubgroup)
    excluded = (nsubgroup eq 0)? list(): subgroups[index]
    index = where(bad_flags eq 0, nsubgroup)
    included = (nsubgroup eq 0)? list(): subgroups[index]
    return, included

end

routines = 'azim_df_search_'+['pre','post']+'_midn_events'
search_step = 'subgroup'
project = azim_df_load_project()
dirname = '2020_04_subgroups_bad_dtime'

foreach routine, routines do begin
    lprmsg, routine
    events = call_function(routine, search_step=search_step)
    goods = azim_df_subgroup_filter_by_dtime(events, excluded=bads)
    lprmsg, 'Total: '+string(events.length)
    lprmsg, 'Excluded: '+string(bads.length)
    lprmsg, 'Remained: '+string(goods.length)
;    stop
    azim_df_subgroup_gen_diagnostic_plot, bads, dirname=dirname, project=project
endforeach

end