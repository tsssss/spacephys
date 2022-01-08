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

pro cusp_prelim_fig_wavepflux, type, ids = ids

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
    fa_sb_wav = dat.fast.sb.fh
    po_kei = dat.polar.kei
    fa_kei = dat.fast.kei

    po_dis = (dat.polar.cusp.entry.dis+dat.polar.cusp.exit.dis)*0.5
    po_mlt = (dat.polar.cusp.entry.mlt+dat.polar.cusp.exit.mlt)*0.5
    fa_dis = (dat.fast.cusp.entry.dis+dat.fast.cusp.exit.dis)*0.5
    dmlt = dat.dmlt
    dr = dat.dr

    ofn = shomedir()+'/cusp/fig/'+'cusp_prelim_fig_wavepflux_'+type+'.eps'
    sgpsopen, ofn, xsize = 8, ysize = 8, /inch
    sgtruecolor
    erase

    poss = sgcalcpos(2,2, margin = [6,6,3,6], ypad = 8, xpad = 11)

    tpos = poss[*,0,0]
    xx = fa_sb_wav
    yy = po_sb_wav
    zz = dr
    xr = [2,2e4]
    yr = [2,2e4]
    zr = [0,5]
    ct = 1

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:1, ylog:1, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar vs FAST'
    zz = bytscl(zz, min = zr[0], max = zr[1])
    sz = 70
    for i = 0, n_elements(xx)-1 do begin
        tcolor = sgcolor(zz[i], ct = ct)
        coord = convert_coord(xx[i],yy[i], /data,/to_device)
        tx = coord[0]
        ty = coord[1]
        tmp = findgen(101)*2*!dpi/100
        txs = tx+sz*cos(tmp)
        tys = ty+sz*sin(tmp)
        polyfill, txs, tys, /device, color = tcolor
    endfor
    plots, xr, yr, linestyle = 2
    width = double(!d.y_ch_size)/!d.x_size*0.6
    pad = double(!d.y_ch_size)/!d.x_size*0.4
    tpos = [tpos[2]+pad,tpos[1],tpos[2]+pad+width,tpos[3]]
    loadct, ct
    sgcolorbar, zrange = zr, position = tpos, ztitle = 'dR (Re)'


    tpos = poss[*,0,1]
    xx = fa_sb_wav
    yy = po_sb_wav
    zz = dmlt
    xr = [2,2e4]
    yr = [2,2e4]
    zr = [0,1]
    ct = 1

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:1, ylog:1, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar vs FAST'
    zz = bytscl(zz, min = zr[0], max = zr[1])
    sz = 70
    for i = 0, n_elements(xx)-1 do begin
        tcolor = sgcolor(zz[i], ct = ct)
        coord = convert_coord(xx[i],yy[i], /data,/to_device)
        tx = coord[0]
        ty = coord[1]
        tmp = findgen(101)*2*!dpi/100
        txs = tx+sz*cos(tmp)
        tys = ty+sz*sin(tmp)
        polyfill, txs, tys, /device, color = tcolor
    endfor
    plots, xr, yr, linestyle = 2
    width = double(!d.y_ch_size)/!d.x_size*0.6
    pad = double(!d.y_ch_size)/!d.x_size*0.4
    tpos = [tpos[2]+pad,tpos[1],tpos[2]+pad+width,tpos[3]]
    loadct, ct
    sgcolorbar, zrange = zr, position = tpos, ztitle = 'dMLT (hr)'

    
    tpos = poss[*,1,0]
    xx = po_dis
    yy = po_kei
    xr = [1,6]
    yr = [-1200,300]

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:0, ylog:0, ytitle:'Polar KEi', xtitle:'Polar R'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar KEi vs R'
    plots, xr, [0,0], linestyle = 1
    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7
    
    
;    tpos = poss[*,1,1]
;    xx = fa_dis
;    yy = fa_kei
;    xr = [1,2]
;    yr = [-100,300]
;    
;    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:0, ylog:0, ytitle:'Polar KEi', xtitle:'Polar R'}
;    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar KEi vs R'
;    plots, xr, [0,0], linestyle = 1
;    tmp = findgen(11)*2*!dpi/10
;    txs = cos(tmp)
;    tys = sin(tmp)
;    usersym, txs, tys, color = 0
;    plots, xx, yy, psym = 8, symsize = 0.7

    tpos = poss[*,1,1]
    xx = po_dis
    yy = po_sb_wav
    xr = [1,6]
    yr = [1,6000]

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:0, ylog:0, ytitle:'Polar S', xtitle:'Polar R'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar S vs R'
    plots, xr, [0,0], linestyle = 1
    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7

    
    sgpsclose, /pdf
end


cusp_prelim_fig_wavepflux, 'south_imf'
end
