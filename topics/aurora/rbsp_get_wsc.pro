
function rbsp_get_wsc, t0, no_spice_load = no_spice_load, $
    no_interpol = no_interp, probe = probe
    
    ; probe, default to 'a', can be 'a', or 'b'.
    if n_elements(probe) eq 0 then prob = 'a' else prob = probe[0]
    
    ; time.
    nt0 = n_elements(t0)
    if nt0 eq 1 then ts = t0[0] $        ; one time.
    else if keyword_set(no_interp) then ts = t0 $   ; slow, be careful!
    else ts = t0[[0,nt0-1]]
    nt = n_elements(ts)
    tstrs = time_string(ts, tformat = 'YYYY-MM-DDThh:mm:ss.ffffff')
    
    ; load spice kernel.
    if not keyword_set(no_spice_load) then begin
        timespan, t0[0], t0[nt0-1]-t0[0], /second
        rbsp_load_spice_kernels
    endif
    
    wsc = dblarr(nt,3)
    for i = 0, nt-1 do begin
        cspice_str2et, tstrs[i], et
        cspice_pxform, 'RBSP'+prob+'_SCIENCE', 'GSE', et, pxform
        wsc[i,*] = pxform[2,*]
    endfor
    
    if keyword_set(no_interp) then return, wsc
    if n_elements(wsc) eq 3 then return, wsc    ; 1-record, no interpol.
    wsc1 = dblarr(nt0,3)
    for i = 0, 2 do wsc1[*,i] = interpol(wsc[*,i], ts, t0, /nan)
    return, wsc1

end