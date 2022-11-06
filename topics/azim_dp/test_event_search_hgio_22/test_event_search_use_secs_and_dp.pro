

;---Read event list of SECS.
    event_list = join_path([googledir(),'works','azim_dp','hgio_22_event_list_v02.txt'])
    lines = read_all_lines(event_list)
    nline = n_elements(lines)
    nevent = nline-1
    event_times = dblarr(nevent)
    event_dict = orderedhash()
    for ii=0,nevent-1 do begin
        infos = strsplit(lines[ii+1], extract=1)
        event_id = infos[0]
        if ~event_dict.haskey(event_id) then event_dict[event_id] = list()
        event_dict[event_id].add, time_double(infos[1])
        if n_elements(infos) le 3 then continue
        infos = strsplit(strjoin(infos[3:*],' '),', ',extract=1)
        foreach info, infos do begin
    ;        if event_id eq '2015_0815' then stop
            if strpos(info,'/') eq -1 then continue
            info = strmid(info,0,16)
            event_dict[event_id].add, time_double(info)
        endforeach
    endfor


;---Test event list.
    if keyword_set(test_event_list) then begin
        event_ids = event_dict.keys()
        nevent = 0
        dtime = list()  ; check typos, should be good if dtime >0 and < several days.
        foreach event_id, event_ids do begin
            event_list = event_dict[event_id]
            print, event_id
            event_id_time = time_double(event_id,tformat='YYYY_MMDD')
            foreach event_time, event_list do begin
                if event_time-event_id_time then $
                    print, time_string(event_time)
            endforeach
            dtime.add, minmax(event_list.toarray()-event_id_time)/86400
            ;print, n_elements(event_list)
            nevent += n_elements(event_list)
        endforeach
        print, nevent
        dtime = dtime.toarray()
    endif


;---Load DP data.
test = 0
    plot_dir = join_path([googledir(),'works','azim_dp','plots','secs_dp'])
    event_ids = event_dict.keys()
    mlt_range = [-1,1]*12d
    rxy_range = [4.,15]
    mlat_range = [60,80]
    pdyn = 20d
    foreach event_id, event_ids do begin
        event_list = event_dict[event_id]
        foreach event_time, event_list do begin
            event_time_range = event_time+[-60,90]*60d  ; sec.
            ;j_up_ewo_var = 'thg_j_up_ewo'
            ;themis_read_upward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=mlat_range
            ;j_down_ewo_var = 'thg_j_down_ewo'
            ;themis_read_downward_current_ewo, time_range, mlt_range=mlt_range, mlat_range=mlat_range
            
if event_time ge time_double('2014-02-27') then continue
            base = 'secs_dp_storm_'+event_id+'_'+time_string(event_time,tformat='YYYY_MMDD_hh')+'_v01.pdf'
            plot_file = join_path([plot_dir,base])
            fig_ewogram_of_dp_and_up_down_current2, event_time_range, $
                test=test, mlt_range=mlt_range, mlat_range=mlat_range, filename=plot_file, rxy_range=rxy_range, pdyn=pdyn
            if keyword_set(test) then stop
        endforeach
    endforeach




end