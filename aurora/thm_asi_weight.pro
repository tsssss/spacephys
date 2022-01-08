
; extract weight for each site and write to
; variable multi in calibration file.

; fn is file name for true_weights variable.

; sites have 0 multi: pgeo, ekat, yknf, nrsq
; sites have ruggy multi: kuuj, snap, talo

pro thm_asi_weight, fn

    rootdir = spreproot()
    fillval = !values.d_nan

    ; load true weight.
    if n_elements(fn) eq 0 then $
        fn = rootdir+'/thg/l2/asi/cal/true_weights_new.sav'
    if file_search(fn) eq '' then message, 'file does not exist ...'
    restore, filename = fn
    tw = temporary(true_weights)    ; tw for true weight.
    dims = size(tw, /dimensions)
    nsite = dims[0]
    
    ; load site names.
    thm_asi_stations, sites, locs
    if n_elements(sites) ne nsite then message, 'site number inconsistent ...'
    sites = strlowcase(sites)

    for i = 0, nsite-1 do begin
        ; for this site.
        site = sites[i]
        sw = reform(tw[i,*,i,1])    ; sw for self weight.
        npx = 256
        
        ; load calibration file.
        ascfn = sptn2fn_thm_asi(0, site, rootdir, type = 'asc')
        cdfid = cdf_open(ascfn)
        vroot = strjoin(['thg','asf',site],'_')
        
        ; read elevation to get the sort index.
        cdf_varget, cdfid, vroot+'_elev', elev, /zvariable
        idx = sort(elev[*])
        
        ; reform sw to get multi.
        multi = reform(sw[sort(idx)],[npx,npx])
        tv, bytscl(multi, /nan), 0
        ; remove edge.
        edge = where(~finite(elev))
        multi[edge] = fillval
        tv, bytscl(multi, /nan), 1
;        if max(multi, /nan) gt 0 then $     ; some sites multi = 0.
            cdf_varput, cdfid, vroot+'_multiply', multi, /zvariable
        cdf_close, cdfid
        
    endfor

end