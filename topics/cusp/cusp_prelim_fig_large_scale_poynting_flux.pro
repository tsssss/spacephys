
pro cusp_prelim_fig_large_scale_poynting_flux

    fn = shomedir()+'/cusp/data/'+'cusp_south_imf_scidat.tplot'
    fn = shomedir()+'/cusp/data/'+'cusp_all_scidat.tplot'
    tplot_restore, filename = fn 
    
    sgindexcolor, 43
    blue = 2
    red = 6
    
    get_data, 'cusp_stats', t0, dat
    
    ids = dat.id
    idx = sort(ids)
    ids = ids[idx]
    
    
;    hem = dat.polar.cusp.entry.ilat gt 0
;    idx = where(hem eq 0, cnt)
;    if cnt ne 0 then hem[idx] = -1
    posbfl = (dat.polar.sb.fl)[idx]
    fasbfl = (dat.fast.sb.fl)[idx]
    sbfldiff = posbfl-fasbfl
    sbflratio = posbfl/fasbfl
    xx = findgen(n_elements(posbfl))
    
    vars = ['posbfl','fasbfl','sbfldiff']
    
    ofn = shomedir()+'/cusp/fig/'+'cusp_prelim_fig_large_scale_poynting_flux.eps'
    sgpsopen, ofn, xsize = 6, ysize = 8, /inch
    erase
    poss = sgcalcpos(3, region = [0.1,0.4,.95,1], ypad = 5)
    ex = {noerase:1, position:poss[*,0], normal:1, yrange:[1e-1,2e4], xrange:[-5,50], ylog:1, xtitle:'Event', ytitle:'S para!C(W/m)'}
    yy = posbfl
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Integrated Poynting flux due to FAC'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, color = blue, psym = 1
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, color = blue
    yy = fasbfl
    oplot, xx, yy, color = red, psym = 2
    oplot, xx, yy, color = red
    xyouts, poss[2,0], poss[1,0]+(poss[3,0]-poss[1,0])*0.67, '  Polar', color = blue, /normal
    xyouts, poss[2,0], poss[1,0]+(poss[3,0]-poss[1,0])*0.33, '  FAST', color = red, /normal
    
    ex = {noerase:1, position:poss[*,1], normal:1, xrange:[-5,50], xtitle:'Event', ytitle:'S para!C(W/m)'}
    yy = sbfldiff
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar minus FAST'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, psym = 4
    plots, minmax(xx), [0,0], linestyle = 2
    
    ex = {noerase:1, position:poss[*,2], normal:1, xrange:[-5,50], xtitle:'Event', ytitle:'S para!C(W/m)', ylog:1, yrange:[1e-4,1e4]}
    yy = sbflratio
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar over FAST'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, psym = 4   
    plots, minmax(xx), [1,1], linestyle = 2

    xr = [1,1e5]
    yr = [1,1e5]
    region = [.1d,.1,.95,.35]
    poss = dblarr(4)
    xsz = !d.x_size*(region[2]-region[0])
    ysz = !d.y_size*(region[3]-region[1])
    len = xsz < ysz
    poss[[0,2]] = (region[0]+region[2])*.5+[-1,1]*(len/!d.x_size)*.5
    poss[[1,3]] = (region[1]+region[3])*.5+[-1,1]*(len/!d.y_size)*.5
    ex = {noerase:1, position:poss, normal:1, xlog:1, ylog:1, xrange:xr, yrange:yr, xtitle:'Polar', ytitle:'FAST'}
    xx = posbfl
    yy = fasbfl
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'Polar vs FAST'
    plot, xx, yy, _extra = ex, xstyle = 5, ystyle = 5, psym = 4
    plots, xr, yr, linestyle = 2
    
    
    sgpsclose, /pdf

end