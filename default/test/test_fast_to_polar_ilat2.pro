; test on dummy var.
sgpsopen, shomedir()+'/test_fast2polar_ilat.eps', xsize = 600, ysize = 900
;sgzopen, shomedir()+'/test_fast2polar_ilat.png', xsize = 600, ysize = 900
sgwopen, 0, xsize = 600, ysize = 900
sgtruecolor
red = sgcolor('red')
green = sgcolor('green')
white = sgcolor('white')
black = sgcolor('black')
erase, color = white

; fast var.
fy = [1,1.2,3.5,2.3]    ; y.
fx = [55,64,70,76]      ; x.
ft = [5d,4,3,2]         ; t.

; polar var.
py = [6,5,7.5,2.3]      ; y.
px = [60,62,70,80]      ; x.
pt = [1d,2,3,4]         ; t.

yr = [0,9]
xr = [55,80]

poss = sgcalcpos(3, ypad = 4, lmargin = 15, bmargin = 10)
xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

; **** fy on ft, with tick label of ft and fx.

xx = ft
yy = fy
tt = [[xx],[fx]]
pos = poss[0,*]
label = 'Fast'
xr = minmax(xx)
yr = sgcalcrange(minmax(yy))
; default settings.
ex0 = {noerase:1, position:[0d,0,0,0], xrange:[0d,0], yrange:[0d,0], $
    nodata:1, xstyle:5, ystyle:5, xtitle:'', ytitle:'', $
    xgridstyle:-1, ygridstyle:-1, xgridv:-1, ygridv:-1, $
    xticks:0, yticks:0, xminor:0, yminor:0, xticklen:0, yticklen:0, $
    xtickvar:'', xtickformat:'(A1)', $
    xticksize:!x.charsize, yticksize:!y.charsize, $
    color:black, background:white}
; plot setting.
ex1 = {noerase:1, position:pos, yrange:yr, xrange:xr, xtickvar:['UT','ILat']}
ex1 = sgstruct(ex0,ex1)
; set coord.
plot, xr, yr, _extra = ex1, ytick_get = ytickv, xtick_get = xtickv
ex1 = sgstruct(ex1, {xtickv:xtickv, ytickv:ytickv})
; plot data.
ex2 = sgstruct(ex1, 'nodata',0)
plot, xx, yy, _extra = ex2
; draw grid.
; draw label.
varlabeloffset = 2
varlabeloffset*= xchsz
nlabel = n_elements(label)
nlines = fltarr(nlabel)
for i = 0, nlabel-1 do nlines[i] = n_elements(strsplit(label[i],'!C',/extract))
labelpos = dblarr(2,nlabel)
labelpos[0,*] = ex1.position[2]+varlabeloffset
labelpos[1,*] = interpol(ex1.position[[1,3]],[0,nlabel+1],indgen(nlabel)+1)+$
    ychsz*(nlines*0.5-1)
xyouts, labelpos[0,*], labelpos[1,*], label, /normal, alignment = 0, color = black
; draw coord.
ex3 = sgstruct(ex1,{xstyle:1, ystyle:1, ytitle:'KE flux!C!C(mW/m!U2!N)'})
plot, xr, yr, _extra = ex3
; prep tick.
xticks = ex1.xtickvar
if xticks[0] ne '' then begin
    xtickvs = ex1.xtickv
    xmajor = n_elements(xtickvs)
    xtickpos = convert_coord(xtickvs, $
        dblarr(xmajor)+ex1.yrange[0], /data, /to_normal)
    ticklabeloffset = 8 > max(strlen(ex1.xtickvar))
    ticklabeloffset *= xchsz
    for i = 0, n_elements(xticks)-1 do begin
        xtickpos[1,*] -= 1.2*ychsz
        txtickvs = interpol(tt[*,i], tt[*,0],xtickvs)
        xticknames = snum2str(txtickvs,/short)
        xyouts, ex1.position[0]-ticklabeloffset, xtickpos[1,0], ex1.xtickvar[i], /normal, alignment = 0, color = black
        xyouts, reform(xtickpos[0,*]), reform(xtickpos[1,*]), xticknames, /normal, alignment = 0.5, color = black
    endfor
endif

; **** py on pt, with tick label of pt and px.

xx = pt
yy = py
tt = [[xx],[px]]
pos = poss[1,*]
label = 'Polar'
xr = minmax(xx)
yr = sgcalcrange(minmax(yy))
; default settings.
ex0 = {noerase:1, position:[0d,0,0,0], xrange:[0d,0], yrange:[0d,0], $
    nodata:1, xstyle:5, ystyle:5, xtitle:'', ytitle:'', $
    xgridstyle:-1, ygridstyle:-1, xgridv:-1, ygridv:-1, $
    xticks:0, yticks:0, xminor:0, yminor:0, xticklen:0, yticklen:0, $
    xtickvar:'', xtickformat:'(A1)', $
    xticksize:!x.charsize, yticksize:!y.charsize, $
    color:black, background:white}
; plot setting.
ex1 = {noerase:1, position:pos, yrange:yr, xrange:xr, xtickvar:['UT','ILat']}
ex1 = sgstruct(ex0,ex1)
; set coord.
plot, xr, yr, _extra = ex1, ytick_get = ytickv, xtick_get = xtickv
ex1 = sgstruct(ex1, {xtickv:xtickv, ytickv:ytickv})
; plot data.
ex2 = sgstruct(ex1, 'nodata',0)
plot, xx, yy, _extra = ex2
; draw grid.
; draw label.
varlabeloffset = 2
varlabeloffset*= xchsz
nlabel = n_elements(label)
nlines = fltarr(nlabel)
for i = 0, nlabel-1 do nlines[i] = n_elements(strsplit(label[i],'!C',/extract))
labelpos = dblarr(2,nlabel)
labelpos[0,*] = ex1.position[2]+varlabeloffset
labelpos[1,*] = interpol(ex1.position[[1,3]],[0,nlabel+1],indgen(nlabel)+1)+$
    ychsz*(nlines*0.5-1)
xyouts, labelpos[0,*], labelpos[1,*], label, /normal, alignment = 0, color = black
; draw coord.
ex3 = sgstruct(ex1,{xstyle:1, ystyle:1, ytitle:'KE flux!C!C(mW/m!U2!N)'})
plot, xr, yr, _extra = ex3
; prep tick.
xticks = ex1.xtickvar
if xticks[0] ne '' then begin
    xtickvs = ex1.xtickv
    xmajor = n_elements(xtickvs)
    xtickpos = convert_coord(xtickvs, $
        dblarr(xmajor)+ex1.yrange[0], /data, /to_normal)
    ticklabeloffset = 8 > max(strlen(ex1.xtickvar))
    ticklabeloffset *= xchsz
    for i = 0, n_elements(xticks)-1 do begin
        xtickpos[1,*] -= 1.2*ychsz
        txtickvs = interpol(tt[*,i], tt[*,0],xtickvs)
        xticknames = snum2str(txtickvs,/short)
        xyouts, ex1.position[0]-ticklabeloffset, xtickpos[1,0], ex1.xtickvar[i], /normal, alignment = 0, color = black
        xyouts, reform(xtickpos[0,*]), reform(xtickpos[1,*]), xticknames, /normal, alignment = 0.5, color = black
    endfor
endif


; **** map ft to pt using x as bridge.
ftp = interpol(pt,px,fx, /spline)
ftpp = interpol(ft,fx,px, /spline)
;fxp = interpol(fx,px,fx, /spline)
xr = sgcalcrange(minmax([ftp,pt]))
yr = sgcalcrange(minmax([fy,py]))
pos = reform(poss[2,*])
label = ['Polar','Fast']
x1 = pt
y1 = py
x2 = ftp
y2 = fy
tt = [[pt],[px],[ftpp]]

; default settings.
ex0 = {noerase:1, position:[0d,0,0,0], xrange:[0d,0], yrange:[0d,0], $
    nodata:1, xstyle:5, ystyle:5, xtitle:'', ytitle:'', $
    xgridstyle:-1, ygridstyle:-1, xgridv:-1, ygridv:-1, $
    xticks:0, yticks:0, xminor:0, yminor:0, xticklen:0, yticklen:0, $
    xtickvar:'', xtickformat:'(A1)', $
    xticksize:!x.charsize, yticksize:!y.charsize, $
    color:black, background:white}
; plot setting.
ex1 = {noerase:1, position:pos, yrange:yr, xrange:xr, xtickvar:['Polar UT','ILat','Fast UT']}
ex1 = sgstruct(ex0,ex1)
; set coord.
plot, xr, yr, _extra = ex1, ytick_get = ytickv, xtick_get = xtickv
ex1 = sgstruct(ex1, {xtickv:xtickv, ytickv:ytickv})
; plot data.
ex2 = sgstruct(ex1, 'nodata',0)
ex2.color = red & plot, x1, y1, _extra = ex2
ex2.color = green & plot, x2, y2, _extra = ex2
; draw grid.
; draw label.
varlabeloffset = 2
varlabeloffset*= xchsz
nlabel = n_elements(label)
labelcolors = [red,green]
nlines = fltarr(nlabel)
for i = 0, nlabel-1 do nlines[i] = n_elements(strsplit(label[i],'!C',/extract))
labelpos = dblarr(2,nlabel)
labelpos[0,*] = ex1.position[2]+varlabeloffset
labelpos[1,*] = interpol(ex1.position[[1,3]],[0,nlabel+1],indgen(nlabel)+1)+$
    ychsz*(nlines*0.5-1)
xyouts, labelpos[0,*], labelpos[1,*], label, /normal, alignment = 0, color = labelcolors
; draw coord.
ex3 = sgstruct(ex1,{xstyle:1, ystyle:1, ytitle:'KE flux!C!C(mW/m!U2!N)'})
plot, xr, yr, _extra = ex3
; prep tick.
xticks = ex1.xtickvar
if xticks[0] ne '' then begin
    xtickvs = ex1.xtickv
    xmajor = n_elements(xtickvs)
    xtickpos = convert_coord(xtickvs, $
        dblarr(xmajor)+ex1.yrange[0], /data, /to_normal)
    ticklabeloffset = 8 > max(strlen(ex1.xtickvar))
    ticklabeloffset *= xchsz
    for i = 0, n_elements(xticks)-1 do begin
        xtickpos[1,*] -= 1.2*ychsz
        txtickvs = interpol(tt[*,i], tt[*,0],xtickvs)
        xticknames = snum2str(txtickvs,/short)
        xyouts, ex1.position[0]-ticklabeloffset, xtickpos[1,0], ex1.xtickvar[i], /normal, alignment = 0, color = black
        xyouts, reform(xtickpos[0,*]), reform(xtickpos[1,*]), xticknames, /normal, alignment = 0.5, color = black
    endfor
endif
;sgzclose
;sgpsclose, /pdf
sgwclose

end
