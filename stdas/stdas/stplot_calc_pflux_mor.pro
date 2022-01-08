;+
; Type: procedure.
; Purpose: Calculate poynting flux S and related variables in tplot.
;   Interface: {in: dename, dbname}{optional in: filter, scales}
;   {out:pfname}{bands of de, db, pf}.
; Parameters:
;   dename, in, string, req. 3-d dE in <coord>.
;   dbname, in, string, req. 3-d dB in <coord>.
;   pfname, in, string, req. 3-d S in <coord>.
; Keywords:
;   power, in, boolean, opt. Set to return S_hat, the 'power'. Default is 
;       to return S(t).
;   filter, in, dblarr[2], opt. Filter in time for dE and dB.
;   scaleinfo, in, struct, opt. Default is 
;       {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}.
;   tscale, in/out, dblarr[m], opt. Scales in time. Info used are
;       minscale, maxscale, nscale.
; Notes: Supply 3-d dE and dB in <coord> to calculate S in <coord>. For example 
;   dename = 'de_facv', dbname = 'db_facv', pfname = 'pf_facv'.
; Dependence: tdas,slib.
; History:
;   2017-07-26, Sheng Tian, create.
;-
pro stplot_calc_pflux_mor, dename, dbname, pfname, $
    power=power, filter=filter, dump_spec=dump_spec, $
    tscale = tscales, scaleinfo = scaleinfo, $
    labels = labels, ytitle = ytitle, colors = colors

;----decompose field into components, check dimension.
    ; input: [de,db,pf]name.
    ; output: [de,db]xyz, nrec, dr0, ndim, dur.
    if n_elements(dename) eq 0 then dename = 'de_fac'
    if n_elements(dbname) eq 0 then dbname = 'db_fac'
    if n_elements(pfname) eq 0 then pfname = 'pf_fac'
    if tnames(dename) eq '' then message, 'no dE data ...'
    if tnames(dbname) eq '' then message, 'no dB data ...'

    get_data, dename, uts, dexyz
    get_data, dbname, uts, dbxyz

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
    ; input: tscales, scaleinfo.
    ; output: s0,s1,dj,ns,scaleinfo.
    if n_elements(scaleinfo) eq 0 then $
        scaleinfo = {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}

    if n_elements(tscales) eq 0 then begin          ; no scales.
        s0 = scaleinfo.s0
        s1 = scaleinfo.s1
        dj = scaleinfo.dj
        ns = scaleinfo.ns

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

    scaleinfo.s0 = s0   ; min scale.
    scaleinfo.s1 = s1   ; max scale.
    scaleinfo.dj = dj   ; 2^dj scale spacing.
    scaleinfo.ns = ns   ; # of scales.
    

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
    



;----options for pfname.
    if n_elements(ytitle) eq 0 then ytitle = '(mW/m!U2!N)'
    if n_elements(labels) ne ndim then labels = ['x','y','z']
    if n_elements(colors) ne ndim then $
        colors = (ndim eq 1)? 0: [6,4,2]
    lims = {ytitle:ytitle, labels:labels, colors:colors}

    store_data, pfname, uts, pfxyz, limits = lims
    tscales = ss

    if ~keyword_set(dump_spec) then begin
        for i = 0, ndim-1 do begin
            tvar = pfname+'_mor_spec_'+string(i+1,format='(I0)')
            store_data, tvar, uts, pfmor[*,*,i], ps, limits = $
                {spec:1, no_interp:1, $
                ytitle:'Period!C(sec)', ylog:1, yrange:[s0,s1], ystyle:1, $
                ztitle:'(mW/m!U2!N)', zlog:0, cdelta:cdelta}
        endfor
    endif
    



end



_2013_0607_load_data
tplot_options, 'yticklen', -0.01
pre0 = 'rbspa_'
stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac_mor'

utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
vars = pre0+['pf_fac_mat','pf_fac_mor']
tplot, vars, trange = utr


;_2013_0607_0456_load_burst_data
;tplot_options, 'yticklen', -0.01
;pre0 = 'rbspa_'
;stplot_calc_pflux_mor, pre0+'eb1_fac', pre0+'mb1_fac', pre0+'pf1_fac'
;vars = pre0+['pf_fac_mor_mor_spec_1']
;tplot, vars
    
end
