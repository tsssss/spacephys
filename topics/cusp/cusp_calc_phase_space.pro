;+
; Given a range in period, and spacecraft velocity, calc the range
; of real omega-k covered by the measurement.
;-

pro cusp_calc_phase_space, vsc, periods, range = pr, log = log

    omega = '!9'+string(119b)+'!X'
    !p.font = 1
    red = sgcolor('red')
    perp = '!9'+string(94b)+'!X'
    
    ws = 2*!dpi/periods
    ks = ws/vsc

    xr = [1e-5,max(ks)*1.5]
    yr = [1e-5,max(ws)*1.5]
    
    if n_elements(pr) ne 0 then begin
        yr = pr
        xr = pr/vsc
    endif
    
    dx = 0.09
    posl = [0.1,0.2,0.4,0.8]+dx*[1,0,1,0]
    posr = [0.6,0.2,0.9,0.8]-dx*[1,0,1,0]
    
    ofn = shomedir()+'/wk1.pdf'
;    ofn = 0
    sgopen, ofn, xsize = 6, ysize = 3, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    plot, xr, yr, psym = 1, xrange = xr, yrange = yr, xstyle = 1, ystyle = 1, $
        xlog = log, ylog = log, /nodata, xtitle = 'k!D'+perp+'!N (km!U-1!N)', ytitle = '', $
        position = posr, /noerase, ytickformat = '(A1)'
;        title = omega+'-k space covered by a spacecraft'
    
    nv0 = 200
    minv0 = 1e-3
    maxv0 = 1d/minv0
    vs = smkgmtrc(minv0,maxv0,nv0,'n')
    for i = 0, nv0-1 do begin
        dv = abs(vs[i]-vsc)
        w0s = abs(ws*vs[i]/dv)
        k0s = w0s/vs[i]
        idx = where(w0s ge yr[0] and w0s le yr[1] and k0s ge xr[0] and k0s le xr[1])
        w0s = w0s[idx]
        k0s = k0s[idx]
        plots, k0s, w0s, psym = 3, symsize = 0.2
;        oplot, xr, xr*vs[i], linestyle = 1
;        wait, 0.5
    endfor
    
;    plots, ks, ws, psym = 1, color = sgcolor('red')
    tc = sgcolor('black')
    oplot, min(ks)+[0,0], yr, color = sgcolor('grey')
    oplot, xr, yr, linestyle = 2, color = tc
    oplot, [xr[0],min(ks)], min(ws)+[0,0], linestyle = 2, color = tc
    oplot, [xr[0],max(ks)], max(ws)+[0,0], linestyle = 2, color = tc
    oplot, min(ks)+[0,0], [yr[0],min(ws)], linestyle = 2, color = tc
    oplot, max(ks)+[0,0], [yr[0],max(ws)], linestyle = 2, color = tc
    
    plots, minmax(ks), minmax(ws), thick = 8, color = sgcolor('red')
    p1 = convert_coord(min(ks),max(yr), /data,/to_normal)

    
    xr = reverse(xr)
    plot, xr, yr, psym = 1, xrange = xr, yrange = yr, xstyle = 1, ystyle = 1, $
        xlog = log, ylog = log, /nodata, xtitle = '-k!D'+perp+'!N (km!U-1!N)', ytitle = omega+' (Hz)', $
        position = posl, /noerase
    ;        title = omega+'-k space covered by a spacecraft'
    
    vs = smkgmtrc(minv0,maxv0,nv0,'n')
    for i = 0, nv0-1 do begin
        dv = abs(vs[i]+vsc)
        w0s = abs(ws*vs[i]/dv)
        k0s = w0s/vs[i]
        idx = where(w0s ge yr[0] and w0s le yr[1] and k0s ge min(xr) and k0s le max(xr))
        w0s = w0s[idx]
        k0s = k0s[idx]
        plots, k0s, w0s, psym = 3, symsize = 0.2
        ;        oplot, xr, xr*vs[i], linestyle = 1
        ;        wait, 0.5
    endfor
    
    ;    plots, ks, ws, psym = 1, color = sgcolor('red')
    tc = sgcolor('black')
    oplot, min(ks)+[0,0], yr, color = sgcolor('grey')
    oplot, xr, reverse(yr), linestyle = 2, color = tc
    oplot, [min(xr),min(ks)], min(ws)+[0,0], linestyle = 2, color = tc
    oplot, [min(xr),max(ks)], max(ws)+[0,0], linestyle = 2, color = tc
    oplot, min(ks)+[0,0], [yr[0],min(ws)], linestyle = 2, color = tc
    oplot, max(ks)+[0,0], [yr[0],max(ws)], linestyle = 2, color = tc
    
    plots, minmax(ks), minmax(ws), thick = 8, color = red
    p2 = convert_coord(min(ks),max(yr), /data,/to_normal)
    
    
    dy = 0.2
    arrow, 0.5*(p1[0]+p2[0]), p1[1]+(1+dy)*ychsz, p1[0], p1[1]+dy*ychsz, /normal, /solid, hsize = 100
    arrow, 0.5*(p1[0]+p2[0]), p2[1]+(1+dy)*ychsz, p2[0], p2[1]+dy*ychsz, /normal, /solid, hsize = 100
    xyouts, 0.5*(p1[0]+p2[0]), p1[1]+(1.2+dy)*ychsz, /normal, $
        'cusp spatial scale', alignment = 0.5, charsize = 0.7

    sgclose

end

vsc = 1 ; km/s.
ps = [10,40,100,500,2000]  ; sec.
ps = smkgmtrc(10,2000,50,'n')
log = 1
pr = 5*[1e-4,1e1]
;log = 0
;pr = [0,1]

cusp_calc_phase_space, vsc, ps, log = log, range = pr
end