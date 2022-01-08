;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro cusp_fig_event_distribution, type, ids = ids

    sgindexcolor, 43
    blue = 2
    red = 6
    
    ; assume type must be there if no id is set.
    if n_elements(ids) eq 0 then ids = cusp_id(type) else type = 'usrdef'

    ; load all wanted data using cusp_gen_excel_form.
    get_data, 'cusp_stats', tmp, dat
    if n_elements(dat) ne n_elements(ids) then begin
        cusp_gen_excel_form, ids, /load
        get_data, 'cusp_stats', tmp, dat
    endif
    
    po_sb_wav = dat.polar.sb.fh
    po_kei = dat.polar.kei
    po_kee = dat.polar.kee

    po_dis = (dat.polar.cusp.entry.dis+dat.polar.cusp.exit.dis)*0.5
    po_mlt = (dat.polar.cusp.entry.mlt+dat.polar.cusp.exit.mlt)*0.5

    po_dis0 = dat.polar.cusp.entry.dis
    po_dis1 = dat.polar.cusp.exit.dis

    po_mlt0 = dat.polar.cusp.entry.mlt
    po_mlt1 = dat.polar.cusp.exit.mlt

    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)

    ofn = shomedir()+'/cusp/fig'+'/cusp_fig_polar_distribution_'+type+'.eps'
    sgpsopen, ofn, xsize = 4, ysize = 6, /inch
    erase

    tpos = sgcalcpos()

    xx = po_mlt
    yy = po_dis
    x0 = po_mlt0
    x1 = po_mlt1
    y0 = po_dis0
    y1 = po_dis1
    xr = [9,15]
    yr = [1,7]
    log = 0
    
    zz = alog10(po_sb_wav)
    zr = [floor(min(zz)),ceil(max(zz))]
    zz = bytscl(zz, min = zr[0], max = zr[1])
    ct = 55

    sgtruecolor

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:log, ylog:log, ytitle:'R (Re)', xtitle:'MLT (hr)'}
    plot, xr, yr, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar Wave S'

    nid = n_elements(ids)
    for i = 0, nid-1 do begin
;        plots, [x0[i],x1[i]], [y0[i],y1[i]], linestyle = 0
        usersym, txs, tys, color = sgcolor(zz[i], ct = ct), /fill
        plots, xx[i], yy[i], psym = 8, symsize = 0.7
    endfor
    
    loadct, ct, /silent
    pad = double(!d.x_ch_size)/!d.x_size
    width = pad*1.5
    tpos = [tpos[2]+pad,tpos[1],tpos[2]+pad+width,tpos[3]]
    sgcolorbar, indgen(256), zrange = zr, position = tpos, ztitle = 'Log10(W/m)'

    sgpsclose, /pdf

end

cusp_fig_event_distribution, 'default'
end
