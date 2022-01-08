function psbl_de_id, type

    if n_elements(type) eq 0 then type = 'default'
    type = strlowcase(type[0])
    
    rootdir = shomedir()+'/Google Drive/works/works/rbsp_de/eventids'
    
    if type eq 'detect_step' then begin
        ifn = rootdir+'/psbl_de_ids_detect_step.dat'
        restore, filename = ifn
        return, ids_detect_step
    endif
    
    
    if type eq 'spin_res' then begin    ; the first round, spin res data events.
        ifn = rootdir+'/psbl_de_ids_spin_res.dat'
        restore, filename = ifn
        return, ids_spin_res
    endif
    
    if type eq 'all_e' then begin       ; the 0th round, all valid E events.
        ifn = rootdir+'/psbl_de_ids_all_e.dat'
        restore, filename = ifn
        return, ids_all_e
    endif

end
