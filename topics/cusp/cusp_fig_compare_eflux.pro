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

pro cusp_fig_compare_eflux, type, ids = ids

    sgindexcolor, 43
    
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
    po_kee = dat.polar.kee
    fa_kee = dat.fast.kee
    
    po_eflux = po_sb_wav+po_kei+po_kee
    fa_eflux = fa_sb_wav+fa_kei+fa_kee

    po_dis = (dat.polar.cusp.entry.dis+dat.polar.cusp.exit.dis)*0.5
    po_mlt = (dat.polar.cusp.entry.mlt+dat.polar.cusp.exit.mlt)*0.5
    fa_dis = (dat.fast.cusp.entry.dis+dat.fast.cusp.exit.dis)*0.5
    dmlt = dat.dmlt
    dr = dat.dr
    
    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)

    ofn = shomedir()+'/cusp/fig'+'/cusp_fig_compare_eflux_'+type+'.eps'
    sgpsopen, ofn, xsize = 8, ysize = 8, /inch
    sgtruecolor
    erase

    poss = sgcalcpos(2,2, margin = [6,6,3,6], ypad = 8, xpad = 11)

    
    ; total eflux.
    tpos = poss[*,0,0]
    range = [10,1e4]
    log = 1
    
    xx = fa_eflux
    yy = po_eflux
;    xx = fa_eflux-po_kei-po_kee
;    yy = po_sb_wav
    xr = range
    yr = range
    
    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:log, ylog:log, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Total Eflux'
    plots, xr, yr, linestyle = 2
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7

    ; poynting flux.
    tpos = poss[*,0,1]
    range = [10,1e4]
    log = 1
    
    xx = fa_sb_wav
    yy = po_sb_wav
    xr = range
    yr = range

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:log, ylog:log, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Wave Poynting flux'
    plots, xr, yr, linestyle = 2
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7

    ; kei.
    tpos = poss[*,1,0]
    range = -1*[-200,1000]

    xx = fa_kei
    yy = po_kei
    xr = range
    yr = range

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:0, ylog:0, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'KEi'
    plots, xr, yr, linestyle = 2
    plots, xr, [0,0], linestyle = 1
    plots, [0,0], yr, linestyle = 1
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7
    

    ; kee.
    tpos = poss[*,1,1]
    range = [-200,500]

    xx = fa_kee
    yy = po_kee
    xr = range
    yr = range

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:0, ylog:0, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'KEe'
    plots, xr, yr, linestyle = 2
    plots, xr, [0,0], linestyle = 1
    plots, [0,0], yr, linestyle = 1
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7

    xyouts, 0.5, 0.95, /normal, alignment = 0.5, 'Polar vs FAST energy fluxes, '+type, charsize = 1.25
    
    sgpsclose, /pdf
end


cusp_fig_compare_eflux, 'default'
end
