;-
; generate 3xn plots for step function example.
; the 3 columns are the time-period spectrogram, the amplitude persentage, and the accumulative amplitude.
; the n rows are for differet orders of approximation.
;+

pro step_amp, np, pos, figlab, single = single

if n_elements(np) eq 0 then np = 4
if n_elements(pos) eq 0 then pos = [0.15,0.15,0.9,0.7]
posl = pos & posl[0] = 0.25 & posl[2] = 0.60
posr = pos & posr[0] = 0.70 & posr[2] = 0.95


pistr = '!9'+string(112b)+'!X'
np0 = 2

nrec = 5000d
w2pi = 2/nrec   ; nrec -> pi.
order = 4       ; order of mat.

;np = 4      ; # of fourier series included.
ps = nrec/smkarthm(1,2,np,'x0')     ; wave periods in # of record.

w0 = min(ps)/2<500  ; min width.
w1 = max(ps)*1.2                        ; max width.
nw = 60                                 ; # of widths.
ws = smkgmtrc(w0,w1,nw,'n')             ; smoothing widths.
ws = fix(ws)
ws = ws[uniq(ws)]


amp0s = dblarr(np)+100                  ; relative amplitude.
amp1s = 100/smkarthm(1,2,np,'x0')       ; absolute amplitude.
amps = dblarr(nw,np)

for i = 0, nw-1 do begin
    u = !dpi*ws[i]/ps
    tamp = (1-sin(u)/u)^(order+1)       ; amplitudes of each band after mat.
    amps[i,*] = tamp*amp0s
    amp0s = amp0s-tamp*amp0s
endfor

ints = amps     ; integrated amplitudes.
for i = 1, nw-1 do begin
    ints[i,*] = ints[i,*]+ints[i-1,*]
endfor


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
yminor = 1
ytickv = [0.125,0.25,0.5,1,2]
ytickn = [' ',' ',pistr+'/2',' ','2'+pistr]

xticklen = -0.03
yticklen = -0.01

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
symsz = 0.3

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

ct = 66



; **** spectrogram.
xr = !dpi*[-1,1]

tpos = posl
nx = nrec*1.5
x = (dindgen(nx)/(nx-1)-0.5)*3*!dpi     ; [-pi,pi].
f0 = dblarr(nx)
for i = 0, np-1 do f0 += 1d/(2*i+1)*sin((2*i+1)*x)


mat = swvmat(f0, order, scale = ws)

z = mat             ; mat.
x = x               ; time in pi.
y = ws*w2pi         ; period/scale/width in pi.

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
title = '' & if single eq 1 then title = 'Time-Period Spectrogram!COrder = '+sgnum2str(order)
if np eq np0 then title = 'Time-Period Spectrogram!COrder = '+sgnum2str(order)
xtitle = ''
xminor = 5
plot, xr, yr, position = tpos, /nodata, /noerase, color = black, title = title, $
    xlog = 0, xstyle = 1, xrange = xr, xtickformat = '(A1)', $
    xticks = xticks, xminor = xminor, xtitle = xtitle, xticklen=xticklen, $
    ylog = 1, ystyle = 1, yrange = yr, yticklen=yticklen, $
    yticks = yticks, yminor = yminor, ytitle = ytitle, ytickv = ytickv, ytickname = ytickn
xtitle = 'Time' & if single ne 1 then xtitle = '' & if np eq 4 then xtitle = 'Time'
if single eq 1 then xtitle = 'Time'
xyouts, (tpos[0]+tpos[2])*0.5, tpos[1]-1.5*chysz, /normal, xtitle, alignment = 0.5, color = black

for i = 0,np-1 do begin
    tp = ps[i]*w2pi
    device, decomposed=0
    loadct2, 43
    oplot, xr, tp+[0,0], linestyle = 2, color = i
    sgtruecolor
    loadct, ct
endfor

; colorbar.
tpos = posl & tpos[0] = tpos[2]+chxsz*0.5 & tpos[2] = tpos[0]+chxsz
sgcolorbar, ctickv, zrange = zr, position = tpos, zticks = cticks, zminor = cminor, ztitle = ctitle

; label.
title = figlab+'. Step Func = '
for i = 0, np-1 do title+= 'f'+sgnum2str(2*i+1)+'+'
title = strmid(title,0,strlen(title)-1)
xyouts, 0.05, tpos[3]-chysz*0.5, /normal, title, color = black

; **** amplitude.
poss = sgcalcpos(1,2, position = posr, xpad = 0.2)


; **** amplitude.
; x is amplitude
xr = [-2,17]
xtickv = [0,5,10,15]
xticks = n_elements(xtickv)-1
tpos = poss[*,0]

; grid, set up coord.
plot, xr, yr, /nodata, /noerase, position = tpos, color = cgrid, $
    xrange = xr, xstyle = 1, xlog = 0, xticklen = 1, xgridstyle = 1, xticks = xticks, xminor = 0, xtickformat = '(A1)', $
    yrange = yr, ystyle = 1, ylog = 1, yticklen = 1, ygridstyle = 1, yticks = yticks, yminor = 0, ytickformat = '(A1)', $
    ytickv = ytickv, xtickv = xtickv

oplot, 0+[0,0], yr, linestyle = 2, color = black
oplot, 100+[0,0], yr, linestyle = 2, color = black
for i = np-1, 0, -1 do begin        ; draw larger period bands later.
    tx = amps[*,i];*amp1s[i]
    ty = ws*w2pi
    tp = ps[i]*w2pi
    device, decomposed=0
    loadct2, 43
    plot, tx, ty, /noerase, position = tpos, $
        xrange = xr, xstyle = 5, xlog = 0, $
        yrange = yr, ystyle = 5, ylog = 1, $
        psym = -4, symsize = symsz, color = i
    oplot, xr, tp+[0,0], linestyle = 2, color = i
    sgtruecolor
endfor

; box.
xticks = n_elements(xtickv)-1
xminor = 5
xtitle = 'Amplitude!CPercent(%)' & if single ne 1 then xtitle = ''
if np eq np0 then xtitle = 'Amplitude!CPercent(%)'
xtickf = '' & if np ne 4 then xtickf = '(A1)'
if single eq 1 then xtickf = ''

plot, xr, yr, /nodata, /noerase, position = tpos, color = black, $
    xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, $
    xtickv = xtickv, xminor = xminor, xtickformat = xtickf, xticklen=xticklen, $
    yrange = yr, ystyle = 1, ylog = 1, yticks = yticks, $
    ytickv = ytickv, yminor = yminor, ytickformat = '(A1)', yticklen=yticklen*2
xyouts, (tpos[0]+tpos[2])*0.5, tpos[3]+ychsz*2, /normal, alignment=0.5, xtitle, color=black


; **** accumulation.
xr = [-10,110]
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
for i = np-1, 0, -1 do begin
    tx = ints[*,i]
    ty = ws*w2pi
    tp = ps[i]*w2pi
    device, decomposed=0
    loadct2, 43
    plot, tx, ty, /noerase, position = tpos, $
        xrange = xr, xstyle = 5, xlog = 0, $
        yrange = yr, ystyle = 5, ylog = 1, $
        psym = -4, symsize = symsz, color = i
    oplot, xr, tp+[0,0], linestyle = 2, color = i
    tmp = convert_coord(xr[1],tp, /data, /to_normal)
    xyouts, tmp[0], tmp[1]-0.4*chysz, /normal, 'P'+sgnum2str(2*i+1), $
        alignment = -0.5, color = i
    sgtruecolor
endfor

; box.
xtickv = [0,50,100]
xticks = n_elements(xtickv)
xminor = 5
xtitle = 'Integrated!CPercent(%)' & if single ne 1 then xtitle = ''
if np eq np0 then xtitle = 'Integrated!CPercent(%)'
xtickf = '' & if np ne 4 then xtickf = '(A1)'
if single eq 1 then xtickf = ''


plot, xr, yr, /nodata, /noerase, position = poss[*,1], color = black, $
    xrange = xr, xstyle = 1, xlog = 0, xticks = xticks, $
    xtickv = xtickv, xminor = xminor, xtickformat = xtickf, xticklen=xticklen, $
    yrange = yr, ystyle = 1, ylog = 1, yticks = yticks, $
    ytickv = ytickv, yminor = yminor, ytickformat = '(A1)', yticklen=yticklen*2
xyouts, (tpos[0]+tpos[2])*0.5, tpos[3]+ychsz*2, /normal, alignment=0.5, xtitle, color=black

end



nps = 1
ofn = shomedir()+'/step_p1_scales.pdf'

nps = [2,4]
ofn = shomedir()+'/fig_step_mat.pdf'
ofn = 0



labs = ['a','b','c']
xsz = 9
ysz = 1.2
tmarg = 3.5d
bmarg = 2d
ypad = 1d

sgopen, 0, xsize = xsz, ysize = ysz, /inch
chysz = double(!d.y_ch_size)/!d.y_size*ysz
sgclose, /wdelete

nnp = n_elements(nps)
if nnp eq 1 then single = 1 else single = 0
ysz = ysz*nnp+(ypad*(nnp-1)+bmarg+tmarg)*chysz

sgopen, ofn, xsize = xsz, ysize = ysz, /inch
poss = sgcalcpos(nnp, ypad = ypad, tmargin = tmarg, bmargin = bmarg, region = [0,0,1,1])
for i = 0, nnp-1 do step_amp, nps[i], poss[*,i], labs[i], single = single
sgclose
end
