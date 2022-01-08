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

pro cusp_test_fig_facpflux, type, ids = ids

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

    po_sb_fac = abs(dat.polar.sb.fl)
    fa_sb_fac = abs(dat.fast.sb.fl)

    ofn = shomedir()+'/cusp/fig/'+'cusp_test_fig_facpflux_'+type+'.eps'
    sgpsopen, ofn, xsize = 4, ysize = 4, /inch
    sgtruecolor
    erase


    tpos = [0.15,0.15,0.85,0.85]
    xx = fa_sb_fac
    yy = po_sb_fac
    xr = [1e-1,1e5]
    yr = [1e-1,1e5]
    ct = 1

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:1, ylog:1, ytitle:'Polar', xtitle:'FAST'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar vs FAST FAC S, '+type
    plots, xr, yr, linestyle = 1
    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7
    
    sgpsclose, /pdf
end


cusp_test_fig_facpflux, 'south_imf'
end
