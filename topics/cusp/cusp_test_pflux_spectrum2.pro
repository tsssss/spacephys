;+
; plot polar and fast pflux spectrum in conjunction pair.
;-


rad = !dpi/180
deg = 180/!dpi
re = 6378d

if n_elements(infos) eq 0 then begin
    infos = []
           
    ids = cusp_id('south_imf')
    foreach id, ids do begin
        tinfo = cusp_test_pflux_spectrum_read_scidat(id, forpolar = 1)
        infos = [infos,tinfo]
    endforeach
    
    ids = cusp_id('south_imf')
    foreach id, ids do begin
        tinfo = cusp_test_pflux_spectrum_read_scidat(id, forpolar = 0)
        infos = [infos,tinfo]
    endforeach
endif

ids = cusp_id('south_imf')
device, decomposed = 1
foreach id, ids do begin
    idx = where(infos.id eq id and infos.name eq 'polar')
    poinfo = infos[idx]
    idx = where(infos.id eq id and infos.name eq 'fast')
    fainfo = infos[idx]
    
    ; convert from deg to sec.
    pocoef = abs(re*poinfo.dis/poinfo.vsc/deg)
    facoef = abs(re*fainfo.dis/fainfo.vsc/deg)
    
    pocoef = abs(re*poinfo.dis/(poinfo.vsc-poinfo.vcv)/deg)
    facoef = abs(re*fainfo.dis/(fainfo.vsc-fainfo.vcv)/deg)
    
    ofn = 0
    sgopen, ofn, xsize = 7, ysize = 6, /inch
    
    ; pure temporal.
    xr = [1,1e4]
    yr = [1,1e4]
    tpos = [0.15,0.1,0.85,0.40]
    plot, xr, yr, xstyle = 1, ystyle = 1, xlog = 1, ylog = 1, /nodata, $
        position = tpos, /noerase, $
        xtitle = 'Temporal Scale, Period (sec)', ytitle = 'Line-integrated S!D||!N(W/m)'
    tmp = poinfo.nfilter
    oplot, poinfo.filters[0:tmp-1], poinfo.pfluxes[0:tmp-1], psym = -1, color = sgcolor('blue')
    faclim = abs(poinfo.faclim*pocoef)
    oplot, faclim+[0,0], yr, linestyle = 1, color = sgcolor('blue')
    tmp = fainfo.nfilter
    oplot, fainfo.filters[0:tmp-1], fainfo.pfluxes[0:tmp-1], psym = -4, color = sgcolor('red')
    faclim = abs(fainfo.faclim*facoef)
    oplot, faclim+[0,0], yr, linestyle = 1, color = sgcolor('red')
    
    ; pure spatial.
    xr = [1e-2,10]
    yr = [1,1e4]
    tpos = [0.15,0.6,0.85,0.9]
    plot, xr, yr, xstyle = 1, ystyle = 1, xlog = 1, ylog = 1, /nodata, $
        position = tpos, /noerase, $
        xtitle = 'Spatial Scale, ILat (deg)', ytitle = 'Line-integrated S!D||!N(W/m)', $
        title = 'Line-integrated S!D||!N spectrum, '+id+'!CPolar-blue, FAST-red'
    tmp = poinfo.nfilter
    oplot, poinfo.filters[0:tmp-1]/pocoef, poinfo.pfluxes[0:tmp-1], psym = -1, color = sgcolor('blue')
    faclim = abs(poinfo.faclim)
    oplot, faclim+[0,0], yr, linestyle = 1, color = sgcolor('blue')
    tmp = fainfo.nfilter
    oplot, fainfo.filters[0:tmp-1]/facoef, fainfo.pfluxes[0:tmp-1], psym = -4, color = sgcolor('red')
    faclim = abs(fainfo.faclim)
    oplot, faclim+[0,0], yr, linestyle = 1, color = sgcolor('red')
    
    sgclose
    
    sgopen, ofn, xsize = 6, ysize = 6, /inch
    
    xr = [1e-5,10]    ; k.
    yr = [1e-4,10]    ; omega.
    
    pocoef = abs(re*poinfo.dis/poinfo.vsc/deg)
    facoef = abs(re*fainfo.dis/fainfo.vsc/deg)
    
    plot, xr, yr, xstyle = 1, ystyle = 1, xlog = 1, ylog = 1, /nodata, $
        ytitle = 'omega (Hz)', xtitle = 'k!Dperp!N (km!U-1!N)', $
        title = sgnum2str(poinfo.dis,msgn=2)+' Re'
    
    ys = 1d/[poinfo.filters[[0,poinfo.nfilter-1]],fainfo.filters[[0,fainfo.nfilter-1]]]
    ys = 1d/[poinfo.filters[0],poinfo.faclim*pocoef,fainfo.filters[0],fainfo.faclim*facoef]
    xs = [ys[0:1]/abs(poinfo.vcv-poinfo.vsc),ys[2:3]/abs(fainfo.vcv-fainfo.vsc)]
    
    for i = 0, n_elements(xs)-1 do begin
        tc = (i le 1)? sgcolor('blue'): sgcolor('red')
        oplot, [xs[i],xr[0]], [yr[0],ys[i]], color = tc
    endfor
    oplot, 1d/(poinfo.faclim*rad*poinfo.dis*re)+[0,0], yr, linestyle = 1, color = sgcolor('blue')
    oplot, 1d/(fainfo.faclim*rad*fainfo.dis*re)+[0,0], yr, linestyle = 1, color = sgcolor('red')
    
    print, poinfo.vcv
    print, poinfo.vsc
    print, abs(poinfo.vcv-poinfo.vsc)

    print, fainfo.vcv
    print, fainfo.vsc
    print, abs(fainfo.vcv-fainfo.vsc)

    
    stop    
    sgclose
endforeach


end
