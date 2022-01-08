; test on real data.
fn = sdiskdir('Works')+'/confs/seminar/seminar_2014_1021/data/seminar_dat.tplot'
tplot_restore, filename = fn
ofn = shomedir()+'/agu_2014_fig_eflux.eps'
eventid = '1998_1001_02'
logfile = sdiskdir('Works')+'/works/cusp/cusp_list_of_conjun.log'
info = cusp_read_conjun_list(logfile, event = eventid)
potr = info.polar.plot_time
fatr = info.fast.plot_time
potrcusp = info.polar.cusp_time
fatrcusp = info.fast.cusp_time


var = 'pf_fac_mat_para_map'
get_data, 'po_'+var, pot0, pvar
get_data, 'fa_'+var, fat0, fvar

get_data, 'po_ilat', tmp, poilat & poilat = interpol(poilat,tmp,pot0)
get_data, 'fa_ilat', tmp, failat & failat = interpol(failat,tmp,fat0)

sgpsopen, shomedir()+'/test_fast2polar_ilat3.eps', xsize = 600, ysize = 900
;sgzopen, shomedir()+'/test_fast2polar_ilat3.png', xsize = 600, ysize = 900
;sgwopen, 0, xsize = 600, ysize = 900
sgtruecolor
red = sgcolor('red')
green = sgcolor('green')
white = sgcolor('white')
black = sgcolor('black')
erase, color = white

; fast var.
fy = fvar
fx = failat      ; x.
ft = fat0        ; t.

; ensure that fast ilat is monotonic.
idx = where(ft ge fatr[0] and ft le fatr[1])
ft = ft[idx]
fx = fx[idx]
fy = fy[idx]

; polar var.
py = pvar        ; y.
px = poilat      ; x.
pt = pot0        ; t.

yr = [0,9]
xr = [55,80]

poss = sgcalcpos(3, ypad = 4, lmargin = 15, bmargin = 8)
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
        if stregex(xticks[i],'UT') ne -1 then xticknames = time_string(txtickvs,tformat='hh:mm')
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
        if stregex(xticks[i],'UT') ne -1 then xticknames = time_string(txtickvs,tformat='hh:mm')
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
plots, xr, [0,0], linestyle = 1
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
        if stregex(xticks[i],'UT') ne -1 then xticknames = time_string(txtickvs,tformat='hh:mm')
        xyouts, ex1.position[0]-ticklabeloffset, xtickpos[1,0], ex1.xtickvar[i], /normal, alignment = 0, color = black
        xyouts, reform(xtickpos[0,*]), reform(xtickpos[1,*]), xticknames, /normal, alignment = 0.5, color = black
    endfor
endif
;sgzclose
sgpsclose, /pdf
;sgwclose

end
