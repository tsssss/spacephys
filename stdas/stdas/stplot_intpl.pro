;+
; Type: procedure.
; Purpose: Interpolate a tplot var to given times.
; Parameters:
;   vars, in, strarr[m], req. tplot vars.
;   t0, in/out, strarr[n], opt. Given times.
; Keywords:
;   tinfo, in, double/dblarr[3], opt. Time info, datarate only or [t1,t2,dr].
; Notes: If set tinfo, then the first var's time will be used.
; Dependence: tdas,slib.
; History:
;   2014-02-11, Sheng Tian, create.
;-

pro stplot_intpl, vars, t0, tinfo = tinfo
    compile_opt idl2

    if n_elements(t0) eq 0 then begin
        if n_elements(tinfo) eq 1 then begin
            dr = tinfo[0]
            get_data, tnames(vars[0]), t0
            t1 = t0[0] & t2 = t0[-1]
            t0 = smkarthm(t1,t2,dr,'dx')
        endif else begin
            t0 = smkarthm(tinfo[0],tinfo[1],tinfo[2],'dx')
        endelse
    endif

    nvar = n_elements(vars)
    for i = 0, nvar-1 do begin
        var0 = tnames(vars[i])
        get_data, var0, t1, f1, limits = lm
        f0 = sinterpol(f1, t1, t0, /nan)
        store_data, var0, t0, f0, limits = lm
    endfor
end
