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

pro cusp_fig_polar_eflux, type, ids = ids

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
    po_kei = dat.polar.kei
    po_kee = dat.polar.kee
    po_eflux = po_sb_wav+po_kei+po_kee

    po_dis = (dat.polar.cusp.entry.dis+dat.polar.cusp.exit.dis)*0.5
    po_mlt = (dat.polar.cusp.entry.mlt+dat.polar.cusp.exit.mlt)*0.5
    
    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)

    ofn = shomedir()+'/cusp/fig'+'/cusp_fig_polar_eflux_'+type+'.eps'
    sgpsopen, ofn, xsize = 8, ysize = 8, /inch
    sgtruecolor
    erase

    poss = sgcalcpos(2,2, margin = [6,6,3,6], ypad = 8, xpad = 11)

    ; pflux vs kei.
    tpos = poss[*,0,0]
    range = [10,1e4]
    log = 1
    
    xx = po_kei
    yy = po_sb_wav
    xr = range
    yr = range
    
    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:log, ylog:log, ytitle:'S (W/m)', xtitle:'KEi (W/m)'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar S vs KEi'
    plots, xr, yr, linestyle = 2
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7

    sgpsclose, /pdf

end

cusp_fig_polar_eflux, 'default'
end
