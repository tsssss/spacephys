
; type can be 'asc' for calibration data, 
;   'asf' for full resolution data,
;   'ast' for thumbnail data.
;
; when type = 'asc', et is optional, supply a dummy variable.


function sptn2fn_thm_asi, et, sites, root, type = type

    sat = 'thg'     ; for all asi types.
    ins = 'asi'     ; instrument.
    sep = '/'
    
    if n_elements(root) eq 0 then $
        root = 'http://themis.ssl.berkeley.edu/data/themis'
    
    if n_elements(type) eq 0 then type = 'asf'
    
    if type eq 'asc' then begin
        dir = strjoin([root,sat,'l2',ins,'cal',''],sep)
        fns = dir+sat+'_l2_'+type+'_'+sites+'_19700101_v01.cdf'
    endif else begin
        cdf_epoch, et, yr, mo, dy, hr, /breakdown_epoch
        yr = string(yr, format = '(I04)')
        mo = string(mo, format = '(I02)')
        dy = string(dy, format = '(I02)')
        hr = string(hr, format = '(I02)')
        ymd = yr+mo+dy
        if type eq 'asf' then ymd += hr
        dir = strjoin([root,sat,'l1',ins],sep)
        fns = dir+sep+sites+sep+yr+sep+mo+sep+$
            sat+'_l1_'+type+'_'+sites+'_'+ymd+'_v01.cdf'    ; version!!!
    endelse
    
    return, fns
    
end