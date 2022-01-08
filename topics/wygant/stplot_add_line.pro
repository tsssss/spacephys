
pro stplot_add_line, vxs, vys

    @tplot_com.pro

    wndw = struct_value(tplot_vars,'settings.window',default=-1)
    if wndw eq -1  then begin
        dprint,dlevel=0,'Tplot has not yet been initialized!,  Returning.'
        return
    endif
    tplot_d = tplot_vars.settings.d
    tplot_x = tplot_vars.settings.x
    tplot_y = tplot_vars.settings.y
    time_scale  = tplot_vars.settings.time_scale
    time_offset = tplot_vars.settings.time_offset
    tplot_var   = tplot_vars.options.varnames
    if keyword_set(cut) then routine_name='tplot_cut'

    wset, tplot_vars.settings.window

    ; get panel id.
    print, 'left click to select panel ...'
    cursor, cx0,cy0, /normal, /down
    pan = where(cy0 ge tplot_y[*].window[0] and cy0 le tplot_y[*].window[1])
    pan = pan[0]
    print, 'panel '+string(pan,format='(I0)')+' was selected ...'
    
    ; start to draw line.
    leftclick = 1
    rightclick = 4
    cxs = []
    cys = []
    color = 0   ; black.
    while !mouse.button ne rightclick do begin
        cursor, tcx,tcy, /normal, /down
        if !mouse.button ne leftclick then continue
        cxs = [cxs,tcx]
        cys = [cys,tcy]
        idx = sort(cxs)
        cxs = cxs[idx]
        cys = cys[idx]
        if n_elements(cxs) gt 1 then begin
            tplot
            plot, cxs, cys, /normal, /noerase, position = [0,0,1,1], color = color, $
                xrange = [0,1], yrange = [0,1], xstyle = 5, ystyle = 5
        endif
        plots, cxs,cys, psym = 1, /normal, color = color
    endwhile

    time_scale  = tplot_vars.settings.time_scale
    time_offset = tplot_vars.settings.time_offset
    vxs = normal_to_data(cxs,tplot_x)*time_scale+time_offset
    vys = normal_to_data(cys,tplot_y[pan])
    
    tvar = tplot_var[pan]
    get_data, tvar, t0
    vys = interpol(vys, vxs, t0)
    vxs = t0

end
