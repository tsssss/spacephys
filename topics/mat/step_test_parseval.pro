;-
; test Parseval's theorem.
;+


ofn = shomedir()+'/fig_step_mat.pdf'
sgopen, ofn, xsize=6, ysize=3, /inch


;---constants and settings.
    pistr = '!9'+string(112b)+'!X'

    nrec = 5000d    ; # of record.
    w2pi = 2/nrec   ; nrec -> pi.
    order = 2       ; order of mat.

    nf = 2          ; # of fourier series included.


;---calculate the bands.
    fs = nrec/smkarthm(1,2,nf,'x0')     ; wave periods in # of record.

    w0 = min(fs)/2<500                  ; min width.
    w1 = max(fs)*1.2                    ; max width.
    nw = 60                             ; # of widths.
    case order of
        4: nw = 60
        2: nw = 40
    endcase
    ws = smkgmtrc(w0,w1,nw,'n')         ; smoothing widths.
    ws = fix(ws)
    ws = ws[uniq(ws)]

    ; amplitude coef.
    amp0s = dblarr(nf)+100              ; relative amplitude.
    amp1s = 100/smkarthm(1,2,nf,'x0')   ; absolute amplitude.
    amps = dblarr(nw,nf)

    for i = 0, nw-1 do begin
        u = !dpi*ws[i]/fs
        tamp = (1-sin(u)/u)^(order+1)   ; amplitudes of each band after mat.
        amps[i,*] = tamp*amp0s
        amp0s = amp0s-tamp*amp0s
    endfor

    ints = amps     ; integrated amplitudes.
    for i = 1, nw-1 do begin
        ints[i,*] = ints[i,*]+ints[i-1,*]
    endfor

    


;---plot settings.
    ytitle = 'Period'
    yr = [0.125,4]
    ytickv = [0.125,0.25,0.5,1,2]
    yticks = n_elements(ytickv)-1
    yminor = 1
    ytickn = [' ',' ',pistr+'/2',' ','2'+pistr]

    xticklen = -0.01
    yticklen = -0.01

    zr = 0.1*[-1,1]                             ; z range.
    zticks = 20                                 ; # of color levels.
    ztickv = smkarthm(zr[0],zr[1],zticks,'n')   ; the color levels.
    ctickv = smkarthm(0,255,zticks,'n')         ; the colors.
    cticks = 2                                  ; # of colorbar major ticks.
    cminor = 0                                  ; # of colorbar minor ticks.
    ctitle = 'Waveform (X)'                     ; colorbar title.

    gstyle = 1                                  ; grid style.
    cgrid = sgcolor('silver')                   ; grid color.
    black = sgcolor('black')

    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    symsz = 0.3

    ct = 66



;---spectrogram.
    xr = !dpi*[-1,1]

    tpos = sgcalcpos()
    nx = nrec*1.5
    xs = (dindgen(nx)/(nx-1)-0.5)*3*!dpi     ; [-pi,pi].
    f0 = dblarr(nx)
    for i = 0, nf-1 do f0 += 1d/(2*i+1)*sin((2*i+1)*xs)


    mat = swvmat(f0, order, scale = ws)
    matps = fltarr(nw)
    for i = 0, nw-1 do matps[i] = max(mat[*,i])
    
    matps = fltarr(nx,nw)
    for i = 0, nw-1 do matps[*,i] = sqrt((mat[*,i])^2+shift(mat[*,i],ws[i]*0.25)^2)
    matps = total(matps,1)/nx
    
    ;plot, matps
    for i=0, nf-1 do oplot, amps[*,i]*0.01*1d/(2*i+1), color=sgcolor('red')
    

;---Parseval's theorem.
    p0 = mean(f0^2)
    print, 'f0 total power:', p0
    
    
    ; filter.
    case nf of
        1: filters = [0,nw-1]
        2: filters = [0,19,nw-1]
        else: message, 'not implemented yet.'
    endcase
    
    matfs = dblarr(nx,nf)
    matp0 = 0
    for i=0, nf-1 do begin
        matfs[*,i] = total(mat[*,filters[i]:filters[i+1]-1],2)
        matp0+= mean(matfs[*,i]^2)
    endfor
    
    print, 'mat total power:', matp0
    print, 'energy loss:', 1-matp0/p0
    
    
;---plot the spectrogram and amplitude curve.
    z = mat        ; mat.
    x = xs         ; time in pi.
    y = ws*w2pi    ; period/scale/width in pi.


    pos1 = [0.10,0.2,0.70,0.8]
    pos2 = [0.75,0.2,0.95,0.8]
    pos1[2] = pos2[0]-xchsz*1

    tpos = pos1

    ; contour.
    sgindexcolor, ct
    contour, z, x, y, position = tpos, /noerase, /fill, $
        xlog = 0, xstyle = 5, xrange = xr, $
        ylog = 1, ystyle = 5, yrange = yr, $
        zlog = 0, levels = ztickv, c_colors = ctickv
    sgtruecolor
    
    ; grid.
    plot, xr, yr, position = tpos, /nodata, /noerase, color = cgrid, $
        xlog = 0, xstyle = 1, xrange = xr, xtickformat = '(A1)', $
        xticks = xticks, xminor = 1, xticklen = 1, xgridstyle = gstyle, $
        ylog = 0, ystyle = 1, yrange = yr, ytickformat = '(A1)', $
        yticks = yticks, yminor = 0, yticklen = 0, ygridstyle = gstyle
    
    ; box.
    title = ''
    xtitle = 'Time'
    xminor = 5
    plot, xr, yr, position = tpos, /nodata, /noerase, color = black, title = title, $
        xlog = 0, xstyle = 1, xrange = xr, xtickformat='(A1)', $
        xticks = xticks, xminor = xminor, xticklen=xticklen, $
        ylog = 1, ystyle = 1, yrange = yr, yticklen=yticklen, $
        yticks = yticks, yminor = yminor, ytitle = ytitle, ytickv = ytickv, ytickname = ytickn
    xyouts, tpos[0], tpos[1]-ychsz, /normal, '-'+pistr, alignment=0.5, color=black
    xyouts, tpos[2], tpos[1]-ychsz, /normal, pistr, alignment=1, color=black
    xyouts, (tpos[0]+tpos[2])*0.5, tpos[1]-ychsz*1.5, /normal, xtitle, alignment=0.5, color=black
    ; fig label.
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz, /normal, alignment=0, 'a'


    for i = 0,nf-1 do begin
        tp = fs[i]*w2pi
        device, decomposed=0
        loadct2, 43
        oplot, xr, tp+[0,0], linestyle = 2, color = i+1
        sgtruecolor
        loadct, ct
    endfor

    ; colorbar.
    tpos = pos1 & tpos[1] = tpos[3]+ychsz*0.5 & tpos[3] = tpos[1]+ychsz*0.5
    sgcolorbar, ctickv, zrange = zr, position = tpos, $
        zticks = cticks, zminor = cminor, ztitle = ctitle, /horizontal


    ; the amplitude.
    tpos = pos2
    
    y = ws*w2pi    ; period/scale/width in pi.
    x = matps      ; amplitude profile for the calculated mat spec.

    yticklen = xticklen*!d.y_size*(tpos[3]-tpos[1])/(!d.x_size*(tpos[2]-tpos[0]))
    xr = [0,0.15]
    xtickv = xr
    xticks = n_elements(xtickv)
    xminor = 5
    
    plot, x, y, position=tpos, /noerase, /nodata, color=black, $
        ystyle=1, ytickformat='(A1)', ylog=1, yrange=yr, yticklen=yticklen, $
        yticks=yticks, ytickv=ytickv, yminor=yminor, $
        xstyle=1, xtickformat='(A1)', xlog=0, xrange=xr, xticklen=xticklen, $
        xticks=xticks, xtickv=xtickv, xminor=xminor
    xyouts, tpos[0], tpos[1]-ychsz, alignment=0, /normal, sgnum2str(xr[0]), color=black
    xyouts, tpos[2], tpos[1]-ychsz, alignment=0.5, /normal, sgnum2str(xr[1]), color=black
    xyouts, (tpos[0]+tpos[2])*0.5, tpos[1]-ychsz*1.5, alignment=0.5, /normal, 'Amplitude (X)', color=black

    device, decomposed=0
    loadct2, 43
    for i=0, nf-1 do oplot, amps[*,i]*0.01/(2*i+1), y, color=i+1
    device, decomposed=1
    oplot, x, y, psym=1, symsize=0.3, color=black
    
    xyouts, tpos[2], tpos[3]+ychsz*1.5, /normal, alignment=1, $
        'MAT order = '+sgnum2str(order)
    xyouts, tpos[2], tpos[3]+ychsz*0.5, /normal, alignment=1, $
        'Energy loss = '+sgnum2str((1-matp0/p0)*100,ndec=1)+'%'
    ; fig label.
    xyouts, tpos[0]+xchsz*0.5, tpos[3]-ychsz, /normal, alignment=0, 'b'
    
    ; draw box again on top of the curves.
    plot, x, y, position=tpos, /noerase, /nodata, color=black, $
        ystyle=1, ytickformat='(A1)', ylog=1, yrange=yr, yticklen=yticklen, $
        yticks=yticks, ytickv=ytickv, yminor=yminor, $
        xstyle=1, xtickformat='(A1)', xlog=0, xrange=xr, xticklen=xticklen, $
        xticks=xticks, xtickv=xtickv, xminor=xminor
    
    sgclose

end
