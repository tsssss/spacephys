;+
; From a given candidate_list, return another list, which has
;   either isolated candidates or overlapping candidates.
;-

function azim_df_sort_overlap_candidate, candidates

test = 0

;---Sort by time_range[0].
    ncandidate = candidates.length
    times = dblarr(ncandidate)
    foreach candidate, candidates, ii do times[ii] = candidate.time_range[0]
    index = sort(times)
    candidates = candidates[index]

;---Convert candidate_list to event_list.
    overlap_list = list()
    current_list = list()
    isolated_list = list()
    event_list = list()
    foreach the_group, candidates, id do begin
        str_id = string(id+1,format='(I0)')
        if keyword_set(test) then lprmsg, 'current candidate: '+str_id
        if current_list.length ne 0 then begin
            pre_group = current_list[-1]
            pre_time = pre_group.time_range
            the_time = the_group.time_range
            if min(the_time) ge max(pre_time) then begin
                if keyword_set(test) then lprmsg, 'no overlap with previous candidate, pop current_list'
                if current_list.length gt 1 then begin
                    if keyword_set(test) then lprmsg, 'add current_list to overlap_list'
                    overlap_list.add, current_list
                    event_list.add, current_list
                endif else begin
                    if keyword_set(test) then lprmsg, 'add current_list to isolated_list'
                    isolated_list.add, current_list
                    event_list.add, current_list
                endelse
                if keyword_set(test) then lprmsg, 'reset current_list'
                current_list = list()
            endif
        endif
        if keyword_set(test) then lprmsg, 'add '+str_id+' to current_list'
        current_list.add, the_group
        if id eq ncandidate-1 then begin
            if current_list.length gt 1 then begin
                if keyword_set(test) then lprmsg, 'add current_list to overlap_list'
                overlap_list.add, current_list
                event_list.add, current_list
            endif else begin
                if keyword_set(test) then lprmsg, 'add current_list to isolated_list'
                isolated_list.add, current_list
                event_list.add, current_list
            endelse
        endif
    endforeach
    
    return, event_list
end