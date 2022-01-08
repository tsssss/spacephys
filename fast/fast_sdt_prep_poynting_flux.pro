;+
; Type: procedure.
; Purpose: FAST, sdt interface, prepare E/B fields for poynting flux.
;   Remove nan, interpolate to uniform data.
;   Interface: {out: 'fa_b0_gei','fa_b','db_facv','de_facv', $
;   'fa_pos','fa_vel','ilat','mlt','dis'}
;   Optional: {out: 'dB_fac_v','B_model','E_NEAR_B','E_ALONG_V'}
; Parameters: fn, in, string, req. Input filename.
; Keywords: facv, in, boolean, opt. Set to name the variable with 'facv'.
; Notes: FAST SDT E/B field are well preprocessed. Only minor adjustment here.
; Dependence: slib, tdas.
; Author: Sheng Tian.
; History: 2013-11-21, Sheng Tian, create.
;-
pro fast_sdt_prep_poynting_flux, fn, facv = facv, trplot = tr, $
    remove_bg=remove_bg

    if file_test(fn) eq 0 then message, 'file does not exist ...'
    tplot_restore, filename = fn
    pre = 'fa_'
    suf = keyword_set(facv)? 'v':''
    
    ; some quantities.
    rgb = [6,4,2]   ; r,g,b in colortable 43,45.
    labfacv = ['v','p (bxv)','b']
    
    ; truncate to given trange.
    if n_elements(tr) eq 2 then begin
        tr = minmax(tr)
        vars = ['dB_fac_v','B_model','E_ALONG_V','E_NEAR_B','ilat','mlt','dis']
        for i = 0, n_elements(vars)-1 do begin
            get_data, vars[i], t0, dat
            idx = where(t0 ge tr[0] and t0 le tr[1], cnt)
            if cnt eq 0 then message, 'invalid time range ...'
            store_data, vars[i], t0[idx], dat[idx,*]
        endfor
    endif
    
    ; get uniform time.
    maxnrec = 100000ul   ; prevent memory overflow.
    get_data, 'dB_fac_v', t0, dbfac
    dr = 1d/8
    t0 = smkarthm(t0[0], t0[-1], dr, 'dx')
    nrec = n_elements(t0)
    print, 'data rate: ', dr
    
    ; db field.
    get_data, 'dB_fac_v', data = tmp
    db_facv = sinterpol(tmp.y, tmp.x, t0, /nan)
    ; remove background field as in polar.
    if keyword_set(remove_bg) then begin
        db_bg = scalcbg(db_facv)
        db_facv = db_facv-db_bg
    endif
    store_data, pre+'db_fac'+suf, data = {x:t0, y:db_facv}, $
        limits = {ytitle: 'dB FAC!C(nT)', labels:labfacv, colors:rgb}
            
    ; model b field in GEI.
    get_data, 'B_model', tmp, bmod_gei
    bmod_gei = sinterpol(bmod_gei, tmp, t0)
    store_data, pre+'b0_gei', t0, bmod_gei, $
        limits = {ytitle:'B model GEI!C(nT)', labels:['x','y','z'], colors:rgb}

    ; total b field.
    get_data, pre+'db_fac'+suf, t0, btotal
    btotal[*,2]+= sqrt(total(bmod_gei^2,2)) & btotal = sqrt(total(btotal^2,2))
    store_data, pre+'b', t0, btotal, limits = $
        {ytitle:'B mag!C(nT)', ynozero:1}
    
    
    
    ; de field.
    get_data, 'E_ALONG_V', data = tmp & dev = sinterpol(tmp.y, tmp.x, t0, /nan)
    get_data, 'E_NEAR_B', data = tmp & deb = sinterpol(tmp.y, tmp.x, t0, /nan)
    ; fix sign error in SDT when FAST in south.
    get_data, 'ilat', tmp, dat
    idx = where(dat le 0, cnt)
    if cnt ne 0 then begin
        idx = where(t0 ge min(tmp[idx]) and t0 le max(tmp[idx]))
        dev[idx] = -dev[idx]
        deb[idx] = -deb[idx]
    endif
    de_facv = [[dev],[dblarr(nrec)],[deb]]      ; [v,bxv,b]
    store_data, pre+'de_fac'+suf, data = {x:t0, y:de_facv}, $
        limits = {ytitle: 'dE (mV/m)', labels:labfacv, colors:rgb}
        
        

    ; remove background field as in polar.
    if keyword_set(remove_bg) then begin
        tvar = pre+'db_fac'
        get_data, tvar, uts, dat
        bg = scalcbg(dat)
        store_data, tvar, uts, dat-bg
        
        tvar = pre+'de_fac'
        get_data, tvar, uts, dat
        bg = scalcbg(dat)
        store_data, tvar, uts, dat-bg
    endif
        
    
    ; ilat, mlt, dis.
    stplot_renew, 'ilat', newname = pre+'ilat'
    stplot_renew, 'mlt', newname = pre+'mlt'
    stplot_renew, 'dis', newname = pre+'dis'
    
    vars = ['dB_fac_v','B_model','E_ALONG_V','E_NEAR_B','ilat','mlt','dis']
    store_data, vars, /delete
end

fn = sdiskdir('Research')+'/data/cusp/fa_sdt_fld_1998_0925_08278.tplot'
fast_sdt_prep_poynting_flux, fn
end