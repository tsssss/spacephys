;+
; Type: procedure.
; Purpose: Blend components to one tplot var.
; Parameter:
;   vars, in, strarr[m], req. Individual components.
; Keyword:
;   delete, in, boolean, opt. Delete old vars.
;   output, in, string, req. Combined tplot var name.
;   ytitle, in, string, opt. Default is ''.
;   labels, in, strarr[m], opt. Default is combining components' labels.
;   colors, in, intarr[m], opt. Colors for components.
;   limits, in, struct, opt. Limit structure received by store_data.
; Notes: none.
; Dependence: tdas,slib.
; History:
;   2013-11-21, Sheng Tian, create.
;-

function stplot_merge, vars, output=newname, ytitle=ytitle, labels=labels, $
    colors=colors, limits=lm, delete=delete
    
    if n_elements(newname) eq 0 then message, 'no newname ...'

    ndim = n_elements(vars)    
    if ndim eq 0 then message, 'no variables ...'
    if ndim eq 1 then return, ''    ; do nothing.
    get_data, vars[0], t0, f0
    nrec = n_elements(t0)
    
    ; ytitle and colors.
    if n_elements(ytitle) eq 0 then ytitle = ''
    if n_elements(colors) eq 0 then begin   ; [r,g,b] for ct42.
        colors = (ndim eq 3)? [6,4,2]: smkarthm(0,255,ndim,'n')
    endif

    ; labels.
    if n_elements(labels) eq 0 then begin
        labels = strarr(ndim)
        for i = 0, ndim-1 do begin
            get_data, vars[i], limits = tmp
            labels[i] = stagexist('labels',tmp)? tmp.labels: ''
        endfor
    endif
    
    ; limits.
    if n_elements(lm) eq 0 then $
        lm = {ytitle:ytitle, colors:colors, labels:labels}

    ; combined signal.
    s0 = make_array(nrec, ndim, type = size(f0[0],/type))
    for i = 0, ndim-1 do begin
        s0[*,i] = get_var_data(vars[i], at=t0)
    endfor

    store_data, newname, t0, s0, limits = lm
    if keyword_set(delete) then begin
        idx = where(vars eq newname, cnt)       ; ensure newname is different.
        if cnt eq 0 then store_data, vars, /delete
    endif
    
    return, newname
end
