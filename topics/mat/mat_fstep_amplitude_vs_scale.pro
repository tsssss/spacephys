
ofn = shomedir()+'/step_amp.pdf'
ofn = 0
pistr = '!9'+string(112b)+'!X'

nrec = 5000d
np = 4      ; # of fourier series included.
nw = 60     ; # of scales.
norder = 4  ; order of mat.

ps = nrec/smkarthm(1,2,np,'x0')     ; wave periods in # of record.
xr = [min(ps)/2>5,max(ps)*1.2]
xtickv = [1,0.5]*nrec
xtickn = (['2',''])+pistr
xticks = n_elements(xtickv)-1
xminor = 10

ws = smkgmtrc(xr[0],xr[1],nw,'n')   ; smoothing widths.

amp0s = dblarr(np)+100              ; relative amplitude.
amp1s = 100/smkarthm(1,2,np,'x0')   ; absolute amplitude.
amps = dblarr(nw,np)

for i = 0, nw-1 do begin
    u = !dpi*ws[i]/ps
    tamp = (1-sin(u)/u)^norder      ; amplitudes of each band after mat.
;    plot, tamp
    amps[i,*] = tamp*amp0s
    amp0s = amp0s-tamp*amp0s
endfor

ints = amps     ; integrated amplitudes.
for i = 1, nw-1 do begin
    ints[i,*] = ints[i,*]+ints[i-1,*]
endfor

tpos = [0.15,0.15,0.9,0.9]
poss = sgcalcpos(2, position = tpos)

sgopen, ofn, xsize = 5, ysize = 3, /inch
sgindexcolor, 43

black = sgcolor('black')

chysz = double(!d.y_ch_size)/!d.y_size

yr = [-5,15]
xr = reverse(xr)

plot, xr, yr, /nodata, /noerase, position = poss[*,0], color = black, $
    xrange = xr, yrange = yr, xtickformat='(A1)', $
    xstyle = 1, ystyle = 1, xlog = 1, $
    ytitle = 'Amplitude!CPercent(%)'
oplot, xr, 0+[0,0], linestyle = 2
oplot, xr, 100+[0,0], linestyle = 2
for i = np-1, 0, -1 do begin
    tx = ws
    ty = amps[*,i];*amp1s[i]
    plot, tx, ty, /noerase, position = poss[*,0], $
        xrange = xr, xstyle = 5, xlog = 1, $
        yrange = yr, ystyle = 5, $
        psym = -4, symsize = 0.5, color = i
    oplot, ps[i]+[0,0], yr, linestyle = 2, color = i
    tmp = convert_coord(ps[i], yr[1], /data, /to_normal)
    xyouts, tmp[0], tmp[1]+0.2*chysz, /normal, 'P'+sgnum2str(2*i+1), $
        alignment = 0.5, color = i
endfor

yr = [-5,105]
plot, xr, yr, /nodata, /noerase, position = poss[*,1], color = black, $
    xrange = xr, yrange = yr, $
    xtickv = xtickv, xtickname = xtickn, xminor = xminor, xticks = xticks, $
    xstyle = 1, ystyle = 1, xlog = 1, $
    ytitle = 'Integrated!CPercent(%)', xtitle = 'Period'
;xyouts, (poss[0,1]+poss[2,1])*0.5, poss[1,1]-chysz*2.2, /normal, 'MAT Scales', alignment = 0.5
oplot, xr, 0+[0,0], linestyle = 2
oplot, xr, 100+[0,0], linestyle = 2
for i = np-1, 0, -1 do begin
    tx = ws
    ty = ints[*,i]
    plot, tx, ty, /noerase, position = poss[*,1], $
        xrange = xr, xstyle = 5, xlog = 1, $
        yrange = yr, ystyle = 5, $
        psym = -4, symsize = 0.5, color = i
    oplot, ps[i]+[0,0], yr, linestyle = 2, color = i
endfor

sgclose

end