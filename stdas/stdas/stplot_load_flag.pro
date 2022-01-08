;+
; mode: smaller/bigger/overlap/equal means minmax(uts) should be 
; within/enclosed/overlap/equal to utr.
; 
;-
function stplot_load_flag, var, trange = utr, $
    mode = mode0

    load = 0
    tvar = var[0]
    if tnames(tvar) eq '' then return, 1    ; no such variable.
    if n_elements(utr) ne 0 then return, load
    
    if n_elements(mode0) eq 0 then mode = mode0 else mode = 'within'
    
    get_data, tvar, uts
    case mode of
        'smaller': flag = (min(uts) ge min(utr)) and (max(uts) le max(utr))
        'bigger': flag = (min(uts) le min(utr)) and (max(uts) ge max(utr))
        'equal': flag = (min(uts) eq min(utr)) and (max(uts) eq max(utr))
        'overlap': flag = ~((min(uts) ge max(utr)) or (max(uts) le max(utr)))
    endcase
    return, load
end