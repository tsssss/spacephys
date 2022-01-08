;+
; Fix Timas epoch. Default is on universal time.
;-

function polar_ion_eflux_fix_timas_time, ut0s, epoch=epoch, error=err

    uts = ut0s
    if keyword_set(epoch) then uts = sfmepoch(ut0s,'unix')

;---Settings.
    dr0 = 12d
    ddr = 3d
    maxerr = 60  ; sec.

    

;---Preparations.
    uts = ut0s
    nrec = n_elements(uts)
    
    ; fix nan.
    idx = where(~finite(uts,/nan), cnt)
    tmp = findgen(nrec)
    uts = interpol(uts[idx],tmp[idx],tmp)
    
    dt0s = (uts-shift(uts,1))[1:*]  ; n-1 # of rec, dt0s = dt1s+ddts+dt1s.
    ddts = round(dt0s/dr0)*dr0      ; the bg of multiple of dr0.
    dt2s = dt0s-ddts
    dt1s = round(dt2s)              ; the rounded fluctuation.
    tdts = dt2s-dt1s                ; the small deviation.
    dt1s = dt1s+ddts
    
    
    test = 0
    if test then begin
        store_data, 'uts', uts, uts, limits={psym:-1, ynozero:1}
        store_data, 'dt1', uts, [0,dt1s], limits={psym:-1, ynozero:1, constant:[6,9,12,15,18,24]}
        tplot, ['uts','dt1']
        stop
    endif
    
   
;---fix dt1s, but make a copy of it.
    dt2s = dt1s
    
    ; case 1: 3-3-6 spikes. conserve time separation.
    idx = where(dt2s eq 18, cnt)
    for i=0, cnt-1 do begin
        ti = idx[i]
        if ti-1 le 0 then continue
        if ti+1 ge nrec-2 then continue
        if dt2s[ti-1] ne 9 then continue
        if dt2s[ti+1] ne 9 then continue
        dt2s[ti-1:ti+1] = 12
    endfor
    idx = where(dt2s eq 6, cnt)
    for i=0, cnt-1 do begin
        ti = idx[i]
        if ti-1 le 0 then continue
        if ti+1 ge nrec-2 then continue
        if dt2s[ti-1] ne 15 then continue
        if dt2s[ti+1] ne 15 then continue
        dt2s[ti-1:ti+1] = 12
    endfor
    
    
    if test then begin
        store_data, 'dt2', uts, [0,dt2s], limits={psym:-1, ynozero:1, constant:[6,9,12,15,18,24]}
        tplot, ['uts','dt1','dt2']
        stop
    endif
    
    
    ; case 2: 2-point bipolar spikes. conserve time separation.
    for i=0, nrec-3 do begin
        if dt2s[i] eq 12 then continue
        if dt2s[i]+dt2s[i+1] ne 24 then continue
        dt2s[i:i+1] = 12
        i = i+1
    endfor
    
    if test then begin
        store_data, 'dt2', uts, [0,dt2s], limits={psym:-1, ynozero:1, constant:[6,9,12,15,18,24]}
        tplot, ['uts','dt1','dt2']
        stop
    endif
    
    ; case 3: single 9/15 or 6/18 points.
    idx = where(abs(dt2s-12) eq 3, cnt)
    for i=0, cnt-1 do begin
        ti = idx[i]
        if ti eq 0 then begin
            if dt2s[1] eq 12 then continue
            tti = 1
        endif else if ti eq nrec-2 then begin
            if dt2s[nrec-3] then continue
            tti = nrec-3
        endif else begin
            if dt2s[ti-1] eq 12 and dt2s[ti+1] eq 12 then continue
            tti = (dt2s[ti-1] eq 12)? ti+1: ti-1
        endelse
        dt2s[tti] = dt2s[tti]+dt2s[ti]-12
        dt2s[ti] = 12
    endfor

    
    if test then begin
        store_data, 'dt2', uts, [0,dt2s], limits={psym:-1, ynozero:1, constant:[6,9,12,15,18,24]}
        tplot, ['uts','dt1','dt2']
        stop
    endif
    
    ; if there are only a few weird points left, wipe them out.
    idx1 = where(dt2s eq 09, cnt1)
    idx2 = where(dt2s eq 15, cnt2)
    idx3 = where(dt2s eq 06, cnt3)
    idx4 = where(dt2s eq 18, cnt4)
    err = (cnt4-cnt3)*6+(cnt2-cnt1)*3
    if err le maxerr then begin
        if cnt1 ne 0 then dt2s[idx1] = 12
        if cnt2 ne 0 then dt2s[idx2] = 12
        if cnt3 ne 0 then dt2s[idx3] = 12
        if cnt4 ne 0 then dt2s[idx4] = 12
    endif
    
    if test then begin
        store_data, 'dt2', uts, [0,dt2s], limits={psym:-1, ynozero:1, constant:[6,9,12,15,18,24]}
        tplot, ['uts','dt1','dt2']
        stop
    endif
 
    
;---Recover the times.
    dt3s = dt2s+tdts
    ut1s = ut0s
    for i=0,nrec-2 do ut1s[i+1] = ut1s[i]+dt3s[i]

    if keyword_set(epoch) then ut1s = stoepoch(ut1s,'unix')
    return, ut1s
 
end



