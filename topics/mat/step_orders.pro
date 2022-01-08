;-
; Plot the time-period spectrogram for step function at different orders using MAT.
;+

    ofn = shomedir()+'/step_orders.pdf'
    ofn = 0

    pistr = '!9'+string(112b)+'!X'
    labs = ['a','b','c','d']
    labx = 0.05

    nrec = 5000d
    w2pi = 2/nrec   ; nrec -> pi.
    orders = [1,2,4]
    orders = 4
    norder = n_elements(orders)  ; order of mat.
    
    np = 4      ; # of fourier series included.
    ps = nrec/smkarthm(1,2,np,'x0')         ; wave periods in # of record.
    
    w0 = min(ps)/2>5                        ; min width.
    w1 = max(ps)*1.2                        ; max width.
    nw = 60                                 ; # of widths.
    ws = smkgmtrc(w0,w1,nw,'n')             ; smoothing widths.
    ws = fix(ws)
    ws = ws[uniq(ws)]
    
    
    sgopen, ofn, xsize = 5, ysize = 6, /inch
    sgtruecolor
    
    ; plot settings.
    xr = !dpi*[-1,1]
    xticks = 4
    xminor = 5
    
    
    pos = [0.3,0.15,0.85,0.9]
    poss = sgcalcpos(norder+1, position = pos)
    
    gstyle = 1                                  ; grid style.
    cgrid = sgcolor('silver')                   ; grid color.
    black = sgcolor('black')
    
    chxsz = double(!d.x_ch_size)/!d.x_size
    chysz = double(!d.y_ch_size)/!d.y_size
    
    ct = 66
    
    
    ; **** original function.
    ; f0 = / 1, [0,pi)
    ;      \-1, [-pi,0).
    ;    = f1+f3+f5+..., where fi = 4/pi/i*sin(x*i).    (1)
    nx = nrec*1.5
    x = (dindgen(nx)/(nx-1)-0.5)*3*!dpi     ; [-pi,pi].
    f0 = dblarr(nx)
    for i = 0, np-1 do f0 += 1d/(2*i+1)*sin((2*i+1)*x)
    
    tpos = poss[*,0]
    yr = 2*[-1,1]
    yticks = 4
    yminor = 2
    
    ; grid.
    plot, xr, yr, /nodata, /noerase, position = tpos, /normal, color = cgrid, $
        xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, xminor = 0, xticklen = 1, xgridstyle = gstyle, xtickformat = '(A1)', $
        yrange = yr, ystyle = 1, ylog = 0, yticks = yticks, yminor = 0, yticklen = 1, ygridstyle = gstyle, ytickformat = '(A1)'

    ; plot.
    plot, x, f0, /noerase, position = tpos, /normal, color = black, $
        xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, xminor = xminor, xtickformat='(A1)', $
        yrange = yr, ystyle = 1, ylog = 0, yticks = yticks, yminor = yminor, ytitle = 'Signal'
    xyouts, labx, tpos[3,0]-0.4*chysz, labs[0]+'. Approx.!C    Step Func', /normal, color = black



    ; **** spectrogram.
    xyouts, (poss[0,0]+poss[2,0])*0.5, poss[1,norder]-1.5*chysz, 'Time', /normal, $
        alignment = 0.5, color = black
        
    
    cticks = 2                                  ; # of colorbar major ticks.
    cminor = 0                                  ; # of colorbar minor ticks.
    ctitle = 'Waveform'                         ; colorbar title.
    
    
    yr = [w0,w1]*w2pi
    ytickv = [0.2,0.4,1,2]
    yticks = n_elements(ytickv)-1
    yminor = 5
    ytickn = strarr(yticks+1) & for i = 0, yticks do ytickn[i] = sgnum2str(ytickv[i])
    idx = where(stregex(ytickn, '1') eq 0, cnt)
    if cnt gt 0 then ytickn[idx] = ''
    ytickn+= pistr
    
    yr = [0.125,4]
    yticks = 4
    yminor = 5
    ytickv = [0.125,0.25,0.5,1,2]
    ytickn = [' ',' ',' ',['','2']+pistr]
    
    for i = 0, norder-1 do begin
        zr = 0.1*[-1,1]                   ; z range.
        zticks = 20                                 ; # of color levels.
        ztickv = smkarthm(zr[0],zr[1],zticks,'n')   ; the color levels.
        ctickv = smkarthm(0,255,zticks,'n')         ; the colors.
        
        mat = swvmat(f0, orders[i], scale = ws)
        stop
        z = mat             ; mat.
        x = x               ; time in pi.
        y = ws*w2pi       ; period/scale/width in pi.
        
        tpos = poss[*,i+1]
        
        ; contour.
        sgindexcolor, ct
        contour, z, x, y, position = tpos, /noerase, /fill, $
            xlog = 0, xstyle = 5, xrange = xr, $
            ylog = 1, ystyle = 5, yrange = yr, $
            zlog = 0, levels = ztickv, c_colors = ctickv
        sgtruecolor
        
        ; grid.
        plot, xr, yr, position = tpos, /nodata, /noerase, color = cgrid, $
            xlog = 0, xstyle = 1, xrange = xr, xticks = xticks, xminor = 0, xticklen = 1, xgridstyle = gstyle, xtickformat = '(A1)', $
            ylog = 1, ystyle = 1, yrange = yr, yticks = yticks, yminor = 0, yticklen = 0, ygridstyle = gstyle, ytickformat = '(A1)', $
            ytickv = ytickv
        
        ; box.
        title = ''
        xtitle = ''
        ytitle = 'Period'
        plot, xr, yr, position = tpos, /nodata, /noerase, color = black, title = title, $
            xlog = 0, xstyle = 1, xrange = xr, xticks = xticks, xminor = xminor, xtitle = xtitle, xtickformat = '(A1)', $
            ylog = 1, ystyle = 1, yrange = yr, yticks = yticks, yminor = yminor, ytitle = ytitle, ytickv = ytickv, ytickname = ytickn
        xyouts, labx, tpos[3]-0.4*chysz, /normal, labs[i+1]+'. MAT!C    order '+sgnum2str(orders[i]), color = sgcolor('black')
        
        sgindexcolor, ct2 = 43
        for j = 0, np-1 do begin
            tp = ps[j]*w2pi
            plots, xr, tp+[0,0], color = j, linestyle = 2
            tmp = convert_coord(xr[1],tp, /data, /to_normal)
            xyouts, tmp[0]+chxsz, tmp[1]-0.4*chysz, /normal, 'P'+sgnum2str(2*j+1), color = j
        endfor
        sgtruecolor
    endfor
   
    sgclose
end