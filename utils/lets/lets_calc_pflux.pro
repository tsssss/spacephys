;+
; Type: procedure.
; Purpose: Calculate poynting flux S and related variables in tplot.
;   Interface: {in: e_var, b_var}{optional in: filter, scales}
;   {out:pf_var}{bands of de, db, pf}.
; Parameters:
;   e_var, in, string, req. 3-d dE in <coord>.
;   b_var, in, string, req. 3-d dB in <coord>.
;   pf_var, in, string, req. 3-d S in <coord>.
; Keywords:
;   power, in, boolean, opt. Set to return S_hat, the 'power'. Default is 
;       to return S(t).
;   filter, in, dblarr[2], opt. Filter in time for dE and dB.
;   scale_info, in, struct, opt. Default is 
;       {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}.
;   tscale, in/out, dblarr[m], opt. Scales in time. Info used are
;       minscale, maxscale, nscale.
; Notes: Supply 3-d dE and dB in <coord> to calculate S in <coord>. For example 
;   e_var = 'de_facv', b_var = 'db_facv', pf_var = 'pf_facv'.
; Dependence: tdas,slib.
; History:
;   2017-07-26, Sheng Tian, create.
;-
function lets_calc_pflux, e_var=e_var, b_var=b_var, pf_var=pf_var, $
    scale_info=scale_info, $
    update=update, get_name=get_name, suffix=suffix, errmsg=errmsg, $
    save_to=data_file, time_var=time_var, $
    power=power, filter=filter, dump_spec=dump_spec, $
    tscale=tscales, $
    labels=labels, ytitle=ytitle, colors=colors

;----decompose field into components, check dimension.
    ; input: [de,db,pf]name.
    ; output: [de,db]xyz, nrec, dr0, ndim, dur.
    if n_elements(e_var) eq 0 then e_var = 'de_fac'
    if n_elements(b_var) eq 0 then b_var = 'db_fac'
    if n_elements(pf_var) eq 0 then pf_var = 'pf_fac'

    ; Check if update in memory.
    var_info = pf_var
    if keyword_set(update) then tmp = delete_var_from_memory(var_info)
    if ~check_if_update_memory(var_info, time_range) then return, var_info

    ; Check if loading from file.
    if keyword_set(update) then tmp = delete_var_from_file(var_info, file=data_file, errmsg=errmsg)
    is_success = read_var_from_file(var_info, file=data_file, errmsg=errmsg)
    if is_success then return, var_info

    ; Read var from routine.
    if tnames(e_var) eq '' then message, 'no dE data ...'
    if tnames(b_var) eq '' then message, 'no dB data ...'

    get_data, e_var, uts, dexyz
    get_data, b_var, uts, dbxyz

    nrec = n_elements(uts)              ; # of records.
    if nrec le 0 then message, 'no data ...'
    dr0 = sdatarate(uts) & dr1 = 1d/dr0 ; data rate.
    tmp = size(dexyz,/dimensions)
    if tmp[0] ne nrec then message, 'dE & dB have different times ...'
    ndim = tmp[1]       ; # of dimensions.
    dur = dr0*nrec      ; duration, in sec.
    
    ; remove nan.
    idx = where(finite(dexyz,/nan), cnt)
    if cnt ne 0 then dexyz[idx] = 0
    idx = where(finite(dbxyz,/nan), cnt)
    if cnt ne 0 then dbxyz[idx] = 0
    
;----prepare scales.
    ; input: tscales, scale_info.
    ; output: s0,s1,dj,ns,scale_info.
    if n_elements(scale_info) eq 0 then begin
        scale_info = dictionary($
            's0', 4d*dr0, $
            's1', 0.5d*dur, $
            'dj', 1d/8, $
            'ns', 0d )
    endif
    if ~scale_info.haskey('dj') then scale_info['dj'] = 1d/8
    if ~scale_info.haskey('ns') then scale_info['ns'] = 0d
    
    if n_elements(tscales) eq 0 then begin          ; no scales.
        s0 = scale_info.s0
        s1 = scale_info.s1
        dj = scale_info.dj
        ns = scale_info.ns

        if s0 eq 0 then s0 = 4*dr0
        if s1 eq 0 then s1 = 0.5d*dur
        if dj eq 0 and ns eq 0 then dj = 1d/8
        if ns eq 0 then begin
            j1 = floor(alog(s1/s0)/alog(2)/dj)  ; # of powers-of-two with dj
            s1 = s0*2d^(dj*j1)            
            ns = j1+1
        endif
        if dj eq 0 then dj = alog(s1/s0)/alog(2)/(ns-1)
    endif else begin
        s0 = min(tscales)
        s1 = max(tscales)
        ns = n_elements(tscales)
        dj = alog(s1/s0)/alog(2)/(ns-1)
    endelse
    j1 = ns-1
    w0 = 6d

    scale_info.s0 = s0   ; min scale.
    scale_info.s1 = s1   ; max scale.
    scale_info.dj = dj   ; 2^dj scale spacing.
    scale_info.ns = ns   ; # of scales.
    

;----calculate Morlet wavelet transform.
    ; input: [de,db]xyz.
    ; output: [de,db]mor.
    cpoynt = 1d/(400d*!dpi) ; from mV/m x nT -> mW/m^2.
    cdelta = 0.776d         ; constant for w0=6, for normalization.
    
    if keyword_set(power) then begin    ; time-averaged.
        demor = complexarr(nrec,ns,ndim)
        for i = 0, ndim-1 do begin
            demor[*,*,i] = wavelet(dexyz[*,i], dr0, /pad, s0=s0, dj=dj, j=j1, $
                mother='Morlet', param=w0, $
                period = ps, scale=ss, coi=coi, signif=signif)
        endfor
        dbmor = complexarr(nrec,ns,ndim)
        for i = 0, ndim-1 do begin
            dbmor[*,*,i] = wavelet(dbxyz[*,i], dr0, /pad, s0=s0, dj=dj, j=j1, $
                mother='Morlet', param=w0, $
                period = ps, scale=ss, coi=coi, signif=signif)
        endfor
        dbmor = conj(dbmor)
        pfmor = complexarr(nrec,ns,ndim)
        if ndim eq 1 then begin
            pfmor = demor*dbmor
        endif else begin
            for i = 0, ndim-1 do begin
                j = (i+1) mod ndim
                k = (i+2) mod ndim
                pfmor[*,*,i] = demor[*,*,j]*dbmor[*,*,k]-demor[*,*,k]*dbmor[*,*,j]
            endfor
        endelse
        pfmor *= cpoynt & demor = 0 & dbmor = 0

        ; normalize.
        for i = 0, ns-1 do pfmor[*,i,*] *= (1d/ss[i])
        pfmor *= (dr0*dj/cdelta)
        pfmor = real_part(pfmor)
        
        if n_elements(filter) eq 0 then idx = findgen(ns) else begin
            idx = where(ps ge min(filter) and ps le max(filter), cnt)
            if cnt eq 0 then message, 'wrong filter ...'
        endelse
        pfxyz = total(pfmor[*,idx,*],2)
    endif else begin                    ; instantaneous.
        demor = dblarr(nrec,ns,ndim)
        for i = 0, ndim-1 do begin
            tmp = wavelet(dexyz[*,i], dr0, /pad, s0=s0, dj=dj, j=j1, $
                mother='Morlet', param=w0, $
                period = ps, scale=ss, coi=coi, signif=signif)
            demor[*,*,i] = real_part(tmp)
        endfor
        dbmor = dblarr(nrec,ns,ndim)
        for i = 0, ndim-1 do begin
            tmp = wavelet(dbxyz[*,i], dr0, /pad, s0=s0, dj=dj, j=j1, $
                mother='Morlet', param=w0, $
                period = ps, scale=ss, coi=coi, signif=signif)
            dbmor[*,*,i] = real_part(tmp)
        endfor
        pfmor = dblarr(nrec,ns,ndim)
        if ndim eq 1 then begin
            pfmor = demor*dbmor
        endif else begin
            for i = 0, ndim-1 do begin
                j = (i+1) mod ndim
                k = (i+2) mod ndim
                pfmor[*,*,i] = demor[*,*,j]*dbmor[*,*,k]-demor[*,*,k]*dbmor[*,*,j]
            endfor
        endelse
        pfmor *= cpoynt & demor = 0 & dbmor = 0

        ; normalize.
        for i = 0, ns-1 do pfmor[*,i,*] *= (1d/ss[i])
        pfmor *= (dr0*dj/cdelta*4)

        if n_elements(filter) eq 0 then idx = findgen(ns) else begin
            idx = where(ps ge min(filter) and ps le max(filter), cnt)
            if cnt eq 0 then message, 'wrong filter ...'
        endelse
        pfxyz = total(pfmor[*,idx,*],2)
    endelse
    



;----options for pf_var.
    if n_elements(ytitle) eq 0 then ytitle = '(mW/m!U2!N)'
    if n_elements(labels) ne ndim then labels = ['x','y','z']
    if n_elements(colors) ne ndim then $
        colors = (ndim eq 1)? 0: constant('rgb')
    lims = {ytitle:ytitle, labels:labels, colors:colors}

    store_data, pf_var, uts, pfxyz, limits = lims
    tscales = ss

    if ~keyword_set(dump_spec) then begin
        for i = 0, ndim-1 do begin
            tvar = pf_var+'_mor_spec_'+string(i+1,format='(I0)')
            store_data, tvar, uts, pfmor[*,*,i], ps, limits = $
                {spec:1, no_interp:1, $
                ytitle:'Period!C(sec)', ylog:1, yrange:[s0,s1], ystyle:1, $
                ztitle:'(mW/m!U2!N)', zlog:0, cdelta:cdelta}
        endfor
    endif
    
    ; Save to file.
    if n_elements(data_file) ne 0 then begin
        is_success = save_var_to_file(var_info, file=data_file, time_var=time_var)
    endif

    return, pf_var


end