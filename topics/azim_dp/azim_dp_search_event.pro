;+
; Search events in given time_range, probes, mlt_range, rxy_range.
;-

function azim_dp_search_event, time_range, probes=probes, $
    mlt_range=mlt_range, rxy_range=rxy_range, pdyn=pdyn, $
    test=test

    event_list = list()
    vars = ['roi_list','roi_dp_list','subseq_list','coherent_candidate_list','coherent_event_list']
    if ~keyword_set(test) then del_data, vars

;----Step 1: Search for DP in ROI.
    the_var = 'roi_list'
    ; Identify ROI.
    if check_if_update(the_var) then begin
        store_data, the_var, 0, $
            azim_dp_search_roi(time_range, probes=probes, $
                mlt_range=mlt_range, rxy_range=rxy_range, pdyn=pdyn, min_dp_count=min_dp_count)
    endif
    roi_list = get_var_data(the_var)
    if roi_list.length eq 0 then return, event_list
    
    ; Search DP in ROI.
    the_var = 'roi_dp_list'
    if check_if_update(the_var) then begin
        store_data, the_var, 0, $
            azim_dp_search_dp_in_roi(roi_list, min_dp_count=min_dp_count)
    endif
    roi_dp_list = get_var_data(the_var)
    if roi_dp_list.length eq 0 then return, event_list

;---Step 2: search for possible subsequence.
    the_var = 'subseq_list'
    if check_if_update(the_var) then begin
        store_data, the_var, 0, $
            azim_dp_search_subseq(roi_dp_list, min_dp_count=min_dp_count)
    endif
    subseq_list = get_var_data(the_var)
    
;---Step 3: search for coherent DP candidate.
    the_var = 'coherent_candidate_list'
    if check_if_update(the_var) then begin
        store_data, the_var, 0, $
            azim_dp_search_coherent_candidate(subseq_list, min_dp_count=min_dp_count)
    endif
    coherent_candidate_list = get_var_data(the_var)
    
;---Step 4: search for coherent DP event.
    the_var = 'coherent_event_list'
    if check_if_update(the_var) then begin
        store_data, the_var, 0, $
            azim_dp_search_coherent_event(coherent_candidate_list, min_dp_count=min_dp_count)
    endif
    coherent_event_list = get_var_data(the_var)


;---Done.
    return, coherent_event_list

end

;time_range = time_double(['2014-08-28/09:30','2014-08-28/11:00'])
time_range = time_double(['2014-08-27','2014-08-29'])
mlt_range = [0,9]


time_range = time_double(['2019-08-05','2019-08-07'])
time_range = time_double(['2019-08-01','2019-08-10'])
;time_range = time_double(['2019-08-25','2019-09-10'])
time_range = time_double(['2019-08-01','2019-09-10'])
;time_range = time_double(['2019-08-30','2019-09-04'])

time_range = time_double(['2019-05-12','2019-05-15'])
time_range = time_double(['2019-05-09','2019-05-12'])
time_range = time_double(['2019-05-11','2019-05-12'])

time_range = time_double(['2017-03-20','2017-03-23'])


probes = ['rbsp'+letters('b'),'th'+['a','d','e'],'g'+['13','14','15','16','17'],'mms1']
test = 0

event_list = list()
mlt_range_list = list()
mlt_range_list.add, [-9,0]
mlt_range_list.add, [0,9]
mlt_range_list.add, [-1,1]*5
rxy_range = [5.,20]
foreach mlt_range, mlt_range_list do begin
    the_list = azim_dp_search_event(time_range, probes=probes, $
        mlt_range=mlt_range, rxy_range=rxy_range, pdyn=pdyn, test=test)
    if the_list.length eq 0 then continue
    event_list.add, the_list, /extract
endforeach

print, event_list.length
foreach event, event_list do begin
    print, ''
    print, time_string(event.time_range)
    print, event.probes
    azim_dp_event_simple_overview, event
endforeach

end
