function cusp_calc_id, id1s, id2s, id3s, id4s, exclude = exids, intersect = inids
        
    if n_elements(id1s) eq 0 then message, 'no id ...'
    ids = id1s
    
    if n_elements(id2s) ne 0 then ids = [ids,id2s]
    if n_elements(id3s) ne 0 then ids = [ids,id3s]
    if n_elements(id4s) ne 0 then ids = [ids,id4s]

    ids = ids[uniq(ids,sort(ids))]
    
    tids = []
    for i = 0, n_elements(inids)-1 do begin
        idx = where(ids eq inids[i], cnt)
        if cnt ne 0 then tids = [tids,inids[i]]
    endfor
    if n_elements(inids) ne 0 then ids = tids
    
    for i = 0, n_elements(exids)-1 do begin
        tid = exids[i]
        idx = where(ids eq tid, cnt)
        if cnt eq 0 then continue
        ids = ids[where(ids ne tid)]
    endfor
    
    return, ids

end
