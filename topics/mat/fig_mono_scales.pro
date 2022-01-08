
ofn = shomedir()+'/mono_scales.pdf'
;ofn = 1
pistr = '!9'+string(112b)+'!X'

nrec = 5000d
w2pi = 2/nrec   ; nrec -> pi.
order = 4       ; order of mat.

np = 1      ; # of fourier series included.
ps = nrec/smkarthm(1,2,np,'x0')     ; wave periods in # of record.

nw = 20     ; # of width.
w0 = min(ps)/5>5    ; min width.
w1 = max(ps)*1.2    ; max width.
ws = smkarthm(w0,w1,nw,'n') ; wdiths for mat.



amp0s = dblarr(np)+100              ; relative amplitude.
amp1s = 100/smkarthm(1,2,np,'x0')   ; absolute amplitude.
amps = dblarr(nw,np)

for i = 0, nw-1 do begin
    u = !dpi*ws[i]/ps
    tamp = (1-sin(u)/u)^norder      ; amplitudes of each band after mat.
    amps[i,*] = tamp*amp0s
    amp0s = amp0s-tamp*amp0s
endfor

ints = amps     ; integrated amplitudes.
for i = 1, nw-1 do begin
    ints[i,*] = ints[i,*]+ints[i-1,*]
endfor




sgopen, ofn, xsize = 6, ysize = 2.5, /inch
sgtruecolor

; plot settings.
yr = [w0,w1]*w2pi
ytitle = 'Period'
ytickv = [1,2]
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

zr = 0.1*[-1,1]                             ; z range.
zticks = 20                                 ; # of color levels.
ztickv = smkarthm(zr[0],zr[1],zticks,'n')   ; the color levels.
ctickv = smkarthm(0,255,zticks,'n')         ; the colors.
cticks = 2                                  ; # of colorbar major ticks.
cminor = 0                                  ; # of colorbar minor ticks.
ctitle = 'Waveform'                         ; colorbar title.

gstyle = 1                                  ; grid style.
cgrid = sgcolor('silver')                   ; grid color.
black = sgcolor('black')

chxsz = double(!d.x_ch_size)/!d.x_size
chysz = double(!d.y_ch_size)/!d.y_size

ct = 66

; **** spectrogram.
pos = [0.15,0.15,0.45,0.7]
xr = !dpi*[-1,1]

tpos = pos
nx = nrec*1.5
x = (dindgen(nx)/(nx-1)-0.5)*3*!dpi     ; [-pi,pi].
f = dblarr(nx)
for i = 0, np-1 do f+= 1/(2*i+1)*sin((2*i+1)*x)
scls = ws
mat = swvmat(f, norder, scale = scls)

z = mat             ; mat.
x = x               ; time in pi.
y = scls*w2pi       ; period/scale/width in pi.

; contour.
sgindexcolor, ct
contour, z, x, y, position = tpos, /noerase, /fill, $
    xlog = 0, xstyle = 5, xrange = xr, $
    ylog = 1, ystyle = 5, yrange = yr, $
    zlog = 0, levels = ztickv, c_colors = ctickv
sgtruecolor

; grid.
xticks = 4
plot, xr, yr, position = tpos, /nodata, /noerase, color = cgrid, $
    xlog = 0, xstyle = 1, xrange = xr, xticks = xticks, xminor = 0, xticklen = 1, xgridstyle = gstyle, xtickformat = '(A1)', $
    ylog = 1, ystyle = 1, yrange = yr, yticks = yticks, yminor = 0, yticklen = 0, ygridstyle = gstyle, ytickformat = '(A1)', $
    ytickv = ytickv

; box.
title = 'Time-Period Spectrogram!CMAT Order = '+sgnum2str(order)
xtitle = ''
xminor = 5
plot, xr, yr, position = tpos, /nodata, /noerase, color = black, title = title, $
    xlog = 0, xstyle = 1, xrange = xr, xticks = xticks, xminor = xminor, xtitle = xtitle, xtickformat = '(A1)', $
    ylog = 1, ystyle = 1, yrange = yr, yticks = yticks, yminor = yminor, ytitle = ytitle, ytickv = ytickv, ytickname = ytickn
xtitle = 'Time'
xyouts, (tpos[0]+tpos[2])*0.5, tpos[1]-1.5*chysz, /normal, xtitle, alignment = 0.5, color = black

for i = 0,np-1 do begin
    tp = ps[i]*w2pi
    loadct2, 43
    oplot, xr, tp+[0,0], linestyle = 2, color = i
    loadct, ct
endfor

; colorbar.
tpos = pos & tpos[0] = tpos[2]+chysz*0.1 & tpos[2] = tpos[0]+chysz*0.2
sgcolorbar, ctickv, zrange = zr, position = tpos, zticks = cticks, zminor = cminor, ztitle = ctitle


; **** amplitude.
pos = [0.6,0.15,0.9,0.7]
poss = sgcalcpos(1,2, position = pos, xpad = 0.6)

; **** amplitude.
; x is amplitude
xr = [-5,20]
xticks = 5
tpos = poss[*,0]

; grid, set up coord.
plot, xr, yr, /nodata, /noerase, position = tpos, color = cgrid, $
    xrange = xr, xstyle = 1, xlog = 0, xticklen = 1, xgridstyle = 1, xticks = xticks, xminor = 0, xtickformat = '(A1)', $
    yrange = yr, ystyle = 1, ylog = 1, yticklen = 1, ygridstyle = 1, yticks = yticks, yminor = 0, ytickformat = '(A1)', $
    ytickv = ytickv

oplot, 0+[0,0], yr, linestyle = 2, color = black
oplot, 100+[0,0], yr, linestyle = 2, color = black
for i = 0, np-1 do begin
    tx = amps[*,i];*amp1s[i]
    ty = ws*w2pi
    tp = ps[i]*w2pi
    loadct2, 43
    plot, tx, ty, /noerase, position = tpos, $
        xrange = xr, xstyle = 5, xlog = 0, $
        yrange = yr, ystyle = 5, ylog = 1, $
        psym = -4, symsize = 0.5, color = black
    oplot, xr, tp+[0,0], linestyle = 2, color = i
endfor

; box.
xtickv = [-5,0,15]
xticks = n_elements(xtickv)
xminor = 0

axis, 0, yr[1], /xaxis, $
    xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, xtickv = xtickv, xminor = xminor, xtickformat = '(A1)', $
    xtitle = 'Amplitude!CPercent(%)', color = black
plot, xr, yr, /nodata, /noerase, position = tpos, color = black, $
    xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, xtickv = xtickv, xminor = xminor, $
    yrange = yr, ystyle = 1, ylog = 1, yticks = yticks, ytickv = ytickv, yminor = yminor, ytickformat = '(A1)'



; **** accumulation.
xr = [-5,105]
xtickv = [0,25,50,75,100]
xticks = n_elements(xtickv)-1
tpos = poss[*,1]

; grid, set up coord.
plot, xr, yr, /nodata, /noerase, position = tpos, color = cgrid, $
    xrange = xr, xstyle = 1, xlog = 0, xticklen = 1, xgridstyle = 1, xticks = xticks, xminor = 0, xtickformat = '(A1)', $
    yrange = yr, ystyle = 1, ylog = 1, yticklen = 1, ygridstyle = 1, yticks = yticks, yminor = 0, ytickformat = '(A1)', $
    ytickv = ytickv, xtickv = xtickv
    
oplot, 0+[0,0], yr, linestyle = 2, color = black
oplot, 100+[0,0], yr, linestyle = 2, color = black
for i = 0, np-1 do begin
    tx = ints[*,i]
    ty = ws*w2pi
    tp = ps[i]*w2pi
    loadct2, 43
    plot, tx, ty, /noerase, position = tpos, $
        xrange = xr, xstyle = 5, xlog = 0, $
        yrange = yr, ystyle = 5, ylog = 1, $
        psym = -4, symsize = 0.5, color = i
    oplot, xr, tp+[0,0], linestyle = 2, color = i
    tmp = convert_coord(xr[1],tp, /data, /to_normal)
    xyouts, tmp[0], tmp[1]-0.4*chysz, /normal, 'P'+sgnum2str(2*i+1), $
        alignment = -0.5, color = i
endfor

; box.
xtickv = [0,50,100]
xticks = n_elements(xtickv)
xminor = 0

axis, 0, yr[1], /xaxis, xtickformat = '(A1)', xticklen = 0, xticks = 1, $
    xtitle = 'Integrated!CPercent(%)', color = black
plot, xr, yr, /nodata, /noerase, position = poss[*,1], color = black, $
    xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, xtickv = xtickv, $
    yrange = yr, ystyle = 1, ylog = 1, yticks = yticks, ytickv = ytickv, yminor = 10, ytickformat = '(A1)'

sgclose

end