;+
; Type: function.
; Purpose: Auto calc limit of given data.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   2014-04-05, Sheng Tian, create.
;-

function sg_autolim, x0, ntick = ntick, maxntick = maxntick
    compile_opt idl2

    sg_prep_data, x0, px = px
    xmax0 = max(double(*px), min = xmin0, /nan)
    sg_free_data, px, x0 = x0

    if xmax0 eq xmin0 then return, [xmin0-1,xmax0+1]
    if finite(xmax0) eq 0 then return, [0,1]
    if finite(xmin0) eq 0 then return, [0,1]
    
    ; check spikes.
    mean_val = mean(x0,/nan)
    stddev_val = stddev(x0,/nan)
    index = where_pro(x0, mean_val+[-1,1]*stddev_val*3)
    xmax0 = max(x0[index], min=xmin0, /nan)
    
    ; calc range.
    cmin = 1.1     ; min dx coef, dmin_calc is larger than dmin*cmid.
    cmax = 1.3     ; max dx coef, dmin_calc will be around dmin*cmax.
    cmid = sqrt(cmin*cmax)  ; mid dx coef, dmin_calc is best around dmin*cmid.
    dx1 = xmax0-xmin0       ; diff b/w data's min/max.
    dx2 = (cmax-cmin)*dx1   ; diff b/w calc min/max limit.
    c0 = 1d
    ; # of candidates (10,20).
    if dx2 lt 10 then while dx2*c0 lt 10 do c0*= 2d
    if dx2 gt 10 then while dx2*c0 gt 20 do c0*= .5d
    dmin = dx1*cmin*c0      ; scaled data min.
    dmax = dx1*cmax*c0      ; scaled data max.
    dmid = dx1*cmid*c0      ; scaled data target.
    
    ; d0s: candidates for calc min/max limit after scale.
    ; d1s: candidates, prime number excluded.
    d0s = indgen(ceil(dmax-dmin))+floor(dmin)
    nd0 = n_elements(d0s)
    if n_elements(maxntick) eq 0 then maxntick = 5
    nticks = [2,3,4,5]
    nticks = nticks[where(nticks le maxntick)]
    nntick = n_elements(nticks)
    d1s = [0]
    for i = 0, nntick-1 do begin
        idx = where((d0s mod nticks[i]) eq 0, cnt)
        if cnt ne 0 then d1s = [d1s,d0s[idx]]
    endfor
    d1s = d1s[1:*]
    d1s = d1s[uniq(d1s, sort(d1s))]
    nd1 = n_elements(d1s)
    
    xmins = (xmin0-(d1s/c0-dx1)*0.5)*c0
    xmaxs = (xmax0+(d1s/c0-dx1)*0.5)*c0
    err = dx2/dx1/c0
    err = 0.05
    for i = 0, nd1-1 do xmins[i] = sround(xmins[i]/c0, error=err)*c0
    for i = 0, nd1-1 do xmaxs[i] = sround(xmaxs[i]/c0, error=err)*c0

    d2s = d1s

    ; weight, smallest wins.
    ; shorter digits are better.
    ws = intarr(nd1)
    for i = 0, nd1-1 do ws[i] = ssgnfig(xmins[i]/c0)+ssgnfig(xmaxs[i]/c0)
        
    ; closer to target middle is good.
    tmp = sqrt(abs(d2s-dmid))
    tmp = tmp/max(tmp)  ; normalize.
    ws+= tmp
    
    ; equal pad on both min/max sides is good.
    tmp = abs((xmins-xmin0*c0)/(xmaxs-xmax0*c0))
    idx = where(tmp lt 1, cnt)
    if cnt ne 0 then tmp[idx] = 1d/tmp[idx]
    ws*= tmp
    
    ; want the smallest ntick.
    tks = intarr(nd1)
    for i = 0, nd1-1 do $
        for j = 0, nntick-1 do $
            if (d2s[i] mod nticks[j]) eq 0 then begin
                tks[i] = nticks[j]
                break
            endif
    ws*= tks
    
    idx = where(ws eq min(ws,/nan), cnt)
    d2s = d2s[idx]
    ntick = tks[idx]
    xmin1 = xmins[idx]/c0
    xmax1 = xmaxs[idx]/c0
    if cnt gt 1 then begin
        tmp = min(tks, idx)
        xmin1 = xmin1[idx]
        xmax1 = xmax1[idx]
        ntick = ntick[idx]
    endif else begin
        xmin1 = xmin1[0]
        xmax1 = xmax1[0]
        ntick = ntick[0]
    endelse

;    print, xmin0, xmax0
;    print, xmin1, xmax1, ntick
;    stop

    return, [xmin1,xmax1]
end

;print, sg_autolim([-25,10])
;print, sg_autolim([-0.5,0.2])
;print, sg_autolim([40,45])
print, sg_autolim([-63.645,52.568])
print, sg_autolim([-15.597,17.931])
print, sg_autolim([-13.643,12.682])
print, sg_autolim([-11.643, 9.867])
print, sg_autolim([-10.972, 6.321])
print, sg_autolim([ -2.210,34.230])
end
