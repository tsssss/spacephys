
pro cusp_prelim_fig_all_poynting_flux

;    fn = shomedir()+'/cusp/data/'+'cusp_all_scidat.tplot'
;    tplot_restore, filename = fn
    
    sgindexcolor, 43
    blue = 2
    red = 6
    
    ids = cusp_id(/default)
        
    idx = sort(ids)
    ids = ids[idx]

    cusp_gen_excel_form, ids, /load    
    
    get_data, 'cusp_stats', tmp, dat
;    hem = dat.polar.cusp.entry.ilat gt 0
;    idx = where(hem eq 0, cnt)
;    if cnt ne 0 then hem[idx] = -1
    posbfl = (dat.polar.sb.fh)[idx]
    fasbfl = (dat.fast.sb.fh)[idx]
    podis = (dat.polar.cusp.entry.dis)[idx]
    sbfldiff = posbfl-fasbfl
    sbflratio = posbfl/fasbfl
    xx = findgen(n_elements(posbfl))
    
    vars = ['posbfl','fasbfl','sbfldiff']
    
    ofn = shomedir()+'/cusp/fig/'+'cusp_prelim_fig_all_poynting_flux_wave.eps'
    sgpsopen, ofn, xsize = 8, ysize = 8, /inch
    sgtruecolor
    erase
    poss = sgcalcpos(2,2, margin = [6,6,3,6], ypad = 8, xpad = 11)
    tpos = poss[*,0,0]
    yy = posbfl
    xr = [0,max(xx)+1]
    yr = [1,1e5]
    ex = {noerase:1, position:tpos, normal:1, ylog:1, xrange:xr, yrange:yr, xtitle:'Event', ytitle:'S para!C(W/m)'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Integrated Poynting flux all freq bands'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, color = blue, psym = 1
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, color = blue
    yy = fasbfl
    oplot, xx, yy, color = red, psym = 2
    oplot, xx, yy, color = red
    xyouts, tpos[2], tpos[1]+(tpos[3]-tpos[1])*0.67, '  Polar', color = blue, /normal
    xyouts, tpos[2], tpos[1]+(tpos[3]-tpos[1])*0.33, '  FAST', color = red, /normal
    
    tpos = poss[*,0,1]
    yy = sbfldiff
    xr = [0,max(xx)+1]
    yr = [-3e4,3e4]
    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xtitle:'Event', ytitle:'S para!C(W/m)'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar minus FAST'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, psym = 4
    plots, xr, [0,0], linestyle = 2
    
    tpos = poss[*,1,0]
    yy = sbflratio
    xr = [0,max(xx)+1]
    yr = [1e-2,1e3]
    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xtitle:'Event', ytitle:'S para!C(W/m)', ylog:1}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar over FAST'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, psym = 4   
    plots, xr, [1,1], linestyle = 2

    tpos = poss[*,1,1]
    yy = abs(posbfl)
    xx = abs(fasbfl)
    zz = podis
    xr = [1,1e5]
    yr = [1,1e5]
    zr = [1,7]
    ct = 55
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
    pad = double(!d.y_ch_size)/!d.x_size
    width = pad*0.5
    tpos = [tpos[2]+pad,tpos[1],tpos[2]+pad+width,tpos[3]]
    loadct, ct
    sgcolorbar, zrange = zr, position = tpos
    
    sgpsclose, /pdf
    
    fastlargers = where(sbflratio lt 1, nfastlarger)
    polarlargers = where(sbflratio ge 1, npolarlarger)
    print, snum2str(nfastlarger)
    print, snum2str(npolarlarger)

end
