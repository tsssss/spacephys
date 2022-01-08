;+
; Type: procedure.
; Purpose: Remove spikes in given m-dimension line plots. 
; Parameters:
;   f0, in/out, dblarr[n,m]/dblarr[n]/string, req. Data, string for tplot var.
;   t0, in/out, dblarr[n], opt. Time or x data. If f0 is tplot var, then t0 
;       will be loaded, otherwise t0 can be set or generated with indgen.
; Keywords:
;   width, in, int, opt. Width to enclose spike for sure, default is 20.
;   nsigma, in, float, opt. Set threshold to select large deriv due to spikes.
; Notes: The width is usually good enough so you do not need to set. Do not set
;   it too short. The nsigma is set to 5 by default, lowing it will fix more
;   subtle spikes at the risk of including large amplitude and high freq wave.
; Dependence: none.
; History:
;   2014-10-01, Sheng Tian, create.
;-

pro sdespike_deriv, t0, f0, width = slotwd, nsigma = nsigma

    ; prepare t0 and f0.
    if n_params() eq 1 then begin
        if size(t0,/type) eq 7 then begin   ; tplot var.
            vname = t0 & get_data, vname, t0, f0
        endif else f0 = t0                  ; set f0 not t0.
    endif

    tmp = size(f0) & if tmp[0] gt 2 then message, 'more than 2 dimensions ...'
    ndim = (tmp[0] eq 1)? 1: tmp[2]
    nrec = n_elements(f0)/ndim
    if n_elements(vname) eq 0 and n_params() eq 1 then t0 = findgen(nrec)

    if n_elements(slotwd) eq 0 then slotwd = 20     ; empirical value.
    if n_elements(nsigma) eq 0 then nsigma = 5      ; empirical value.
    
    ; work on each dimension.
    for i = 0, ndim-1 do begin
        tf0 = f0[*,i]                               ; original data.
        df0 = [0d,tf0[1:nrec-1]-tf0[0:nrec-2]]      ; derivative.
        df_std = stddev(df0)

        ; flag0 is the crucial info. 1st pass: find large deriv.
        flag0 = abs(df0) lt nsigma*df_std           ; 1 for good, 0 for bad.

        ; 2nd pass: running standard deviation.
        df0_all = [reverse(df0[1:slotwd/2+1]),df0, $
            reverse(df0[nrec-2-slotwd/2:nrec-2])]
        df0_std = df0
        for j = 0, nrec-1 do $
            df0_std[j] = stddev(df0_all[j:j+slotwd])
        ; deviate enough locally, exclude large deriv from wave.
        flag0 = flag0 or (df0_std lt nsigma*df_std) 

        ; no spike.
        idx = where(flag0 eq 0, cnt)
        if cnt eq 0 then begin
            message, 'no spike...', /continue & continue
        endif

        ; extend slot to given width.
        flag1 = flag0       ; extended slot, wide enough to contain spikes.
        for j=0,cnt-1 do flag1[(idx[j]-slotwd/2>0):(idx[j]+slotwd/2<nrec-1)] = 0
        slotedges = where([0, flag1[0:nrec-2] xor flag1[1:nrec-1]] eq 1)
        slots = [[slotedges[0:*:2]],[slotedges[1:*:2]]] ; slot start&end index.
        nslot = n_elements(slots)/2
        tidx = lonarr(nslot,4)      ; save for some useful indices.
        ratios = fltarr(nslot)      ; type ratio.
        for j = 0, nslot-1 do begin
            ; type ratio, must be possitive and >= 1.
            segment = df0[slots[j,0]:slots[j,1]]
            ratios[j] = abs(max(segment)/min(segment))
            if ratios[j] lt 1 then ratios[j] = 1d/ratios[j]
            ; the larger range for deriv interpolation.
            idx1 = slots[j,0]-slotwd > 0
            idx2 = slots[j,1]+slotwd < nrec-1
            ; the strict range for bad deriv and data point.
            tmp = where(flag0[slots[j,0]:slots[j,1]] eq 0, cnt)
            idx3 = slots[j,0]-0.25*slotwd+tmp[0] > 0
            idx4 = slots[j,0]+0.25*slotwd+tmp[cnt-1] < nrec-1
            ; fix bad deriv using interpolation.
            df0[idx3:idx4] = interpol([df0[idx1:idx3],df0[idx4:idx2]], $
                [t0[idx1:idx3],t0[idx4:idx2]], t0[idx3:idx4])
            tidx[j,*] = [idx1,idx2,idx3,idx4]
        endfor
        j = 0   ; treat every spike.
        while j lt nslot do begin
            if ratios[j] le 5 then begin    ; type 1, real spike.
                for k = tidx[j,2], tidx[j,3]-1 do tf0[k] = tf0[k-1]+df0[k] 
                j = j+1
            endif else begin    ; type 2, pseudo spike.
                if j ge nslot-1 then break
                if tidx[j+1,3]-tidx[j,2] lt 5*slotwd then $
                    for k = tidx[j,2], tidx[j+1,3]-1 do tf0[k] = tf0[k-1]+df0[k] 
                j = j+2
            endelse
        endwhile
        f0[*,i] = tf0
    endfor
    if n_elements(vname) eq 1 then store_data, vname, t0, f0 $
    else if n_params() eq 1 then t0 = f0
end
