
; settings.
utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.


tprobe='a'
pre0 = 'rbsp'+tprobe+'_'
dr0 = 1d/4096
tpad = 200*dr0 ; sec.
p0 = 0.02   ; sec.
utr = utr0+[-1,1]*tpad
timespan, utr[0], utr[1]-utr[0], /second


re = 6378d & re1 = 1d/re
rgb = sgcolor(['red','green','blue'])
fac = ['para','west','north']


tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 0
tplot_options, 'version', 1
tplot_options, 'ygridstyle', 0
tplot_options, 'yticklen', 0



_2013_0607_0456_load_burst_data



; **** settings.
get_data, 'morlet_info', tmp, info
filters = info.filters
nband = n_elements(filters)-1
bandids = string(findgen(nband)+1,format='(I0)')
periods = info.periods

w0 = 6
s2p = 4*!dpi/(w0+sqrt(2+w0^2))
s0 = info.s0
dj = info.dj
nj = info.nj

vars = pre0+['eb1_fac','mb1_fac']
ztits = ['(mV/m)','10!U-3!N (nT)']

freqs = findgen(nrec)/(nrec*dr0)
freqs = freqs[0:nrec/2]
xr = [2,200]


; calculate the spectrum using Morlet and FFT.
foreach tvar, vars do begin
    for i = 0, nband-1 do begin
        get_data, tvar+'_f'+bandids[i], uts, dat
        
        ; morlet power spectrum.
        morpows = dblarr(n_elements(periods),3)
        for j = 0, 2 do begin
            tdat = dat[*,j]
            mor = wv_cwt(tdat, 'Morlet', w0, /pad, $
                start_scale = s0/dr0, dscale = dj, nscale = nj, scale = recscls)
            timescls = recscls*dr0
            periods = timescls*s2p
            morpow = abs(mor^2)*2*!dpi*(uts[nrec-1]-uts[0])/nrec
            morspec = total(morpow,1)/nrec
            morpows[*,j] = morspec
        endfor
        store_data, tvar+'_f'+bandids[i]+'_morpow', 1d/periods, morpows, $
            limits = {colors:rgb, labels:fac, xrange:xr}
        
        ; fft power spectrum.
        fftpows = dblarr(nrec/2+1,3)
        for j = 0, 2 do begin
            tdat = dat[*,j]
            fftpow = (abs(fft(tdat)))^2*2*!dpi*(uts[nrec-1]-uts[0])
            fftpows[*,j] = fftpow[0:nrec/2]
        endfor
        store_data, tvar+'_f'+bandids[i]+'_fftpow', freqs, fftpows, $
            limits = {colors:rgb, labels:fac, xrange:xr}
         
;        plot, 1d/periods, snorm(morpows), /nodata, $
;            xrange = xr, xstyle = 1, xlog = 1
;        for j = 0, 2 do oplot, 1d/periods, morpows[*,j], color = rgb[j]
    endfor
endforeach







ofn = shomedir()+'/es_wave/fig_de_decompose.pdf'
;ofn = 0
sgopen, ofn, xsize = 5, ysize = 4, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


x0 = 0.0+xchsz*10
x1 = x0+0.55
x2 = x1+xchsz*1
x3 = 1-xchsz*5

poss = sgcalcpos(3,2, lmargin=10,rmargin=5)

    tvar = pre0+'eb1_fac'

    ; ** power spectrum.
    tpos = reform(poss[*,0,*])
    tpos[0,*] = x2
    tpos[2,*] = x3
    xr = [5,200]
    yr = [1e-4,1e2]
    wfs = dblarr(3) ; wave periods.
    for i = 0, 2 do begin
        ttvar = tvar+'_f'+bandids[i]+'_morpow'
        get_data, ttvar, freq, dat
        xtkfmt = (i eq 2)? '': '(A1)'
        xtitl = (i eq 2)? 'Freq (Hz)': ''
        plot, freq, snorm(dat), /nodata, /noerase, position = tpos[*,i], $
            xstyle = 1, xlog = 1, xrange = xr, xtickformat = xtkfmt, $
            xtickv = [10,100], xticks = 1, xminor = 5, xtitle = xtitl, $
            ystyle = 9, ylog = 1, yrange = yr, ytickformat = '(A1)', $
            ytickv = [1e-4,1e-1,1e2], yticks = 2, yminor = 5, $
            ytickname = string([-4,-1,2],format='(I0)')
        axis, yaxis = 1, $
            ystyle = 9, ylog = 1, yrange = yr, $
            ytickv = [1e-4,1e-1,1e2], yticks = 2, yminor = 5, $
            ytickname = string([-4,-1,2],format='(I0)')
            
        tdat = snorm(dat)
        tmp = max(tdat,idx)
        plots, freq[idx]+[0,0], yr, linestyle = 1
        wfs[i] = freq[idx]
        
        for j = 0, 2 do oplot, freq, dat[*,j], color = rgb[j]
        xyouts, tpos[0,i]+xchsz*0.5, tpos[3,i]-ychsz*1, /normal, $
            alignment = 0, 'b'+bandids[i]+'.'
    endfor
    
    xyouts, (tpos[0,0]+tpos[2,0])*0.5, tpos[3,0]+ychsz*0.5, /normal, $
        alignment = 0.5, 'PS (mV/m)!U2!N/Hz'



    ; ** waveform.
    tpos = reform(poss[*,0,*])
    tpos[0,*] = x0
    tpos[2,*] = x1
    ttvar = tvar+'_f'+bandids
    options, ttvar, 'ytitle', ''
    options, ttvar, 'yticks', 2
    options, ttvar, 'yminor', 5
    options, ttvar, 'labels', ''
    tplot, ttvar, position = tpos, /noerase
    
    for i = 0, 2 do begin
        xyouts, tpos[0,i]+xchsz*1, tpos[3,i]-ychsz*1, /normal, alignment = 0, $
            'a'+bandids[i]+'.'
        xyouts, tpos[0,i]-xchsz*9, (tpos[1,i]+tpos[3,i])*0.5, /normal, $
            alignment = 0, 'Wave '+bandids[i]+'!Cf'+bandids[i]+'='+$
            sgnum2str(round(wfs[i]))+' Hz'
    endfor
    xyouts, (tpos[0,0]+tpos[2,0])*0.5, tpos[3,0]+ychsz*0.5, /normal, $
        alignment = 0.5, 'dE FAC waveform (mV/m)'


    ; ** labels.
    ty = tpos[3,1]-ychsz*0.85
    tx = tpos[2,1]-xchsz*3.5
    for i = 0, 2 do begin
        tty = ty-ychsz*i*0.7
        xyouts, tx, tty, /normal, alignment = 0, $
            fac[i], color = rgb[i], charsize = 0.8
        plots, tx+[-2,-0.5]*xchsz, tty+ychsz*0.2, /normal, color = rgb[i]
    endfor
    
sgclose








ofn = shomedir()+'/es_wave/fig_db_decompose.pdf'
;ofn = 0
sgopen, ofn, xsize = 5, ysize = 4, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


x0 = 0.0+xchsz*10
x1 = x0+0.55
x2 = x1+xchsz*1
x3 = 1-xchsz*5

poss = sgcalcpos(3,2, lmargin=10,rmargin=5)

tvar = pre0+'mb1_fac'

; ** power spectrum.
tpos = reform(poss[*,0,*])
tpos[0,*] = x2
tpos[2,*] = x3
xr = [5,200]
yr = [1e-4,1e2]*1e-6
;wfs = dblarr(3) ; wave periods.
for i = 0, 2 do begin
    ttvar = tvar+'_f'+bandids[i]+'_morpow'
    get_data, ttvar, freq, dat
    xtkfmt = (i eq 2)? '': '(A1)'
    xtitl = (i eq 2)? 'Freq (Hz)': ''
    plot, freq, snorm(dat), /nodata, /noerase, position = tpos[*,i], $
        xstyle = 1, xlog = 1, xrange = xr, xtickformat = xtkfmt, $
        xtickv = [10,100], xticks = 1, xminor = 5, xtitle = xtitl, $
        ystyle = 9, ylog = 1, yrange = yr, ytickformat = '(A1)', $
        ytickv = [1e-4,1e-1,1e2]*1e-6, yticks = 2, yminor = 5, $
        ytickname = string([-4,-1,2],format='(I0)')
    axis, yaxis = 1, $
        ystyle = 9, ylog = 1, yrange = yr, $
        ytickv = [1e-4,1e-1,1e2]*1e-6, yticks = 2, yminor = 5, $
        ytickname = string([-4,-1,2],format='(I0)')
        
    plots, wfs[i]+[0,0], yr, linestyle = 1
    
    for j = 0, 2 do oplot, freq, dat[*,j], color = rgb[j]
    xyouts, tpos[0,i]+xchsz*0.5, tpos[3,i]-ychsz*1, /normal, $
        alignment = 0, 'b'+bandids[i]+'.'
endfor

xyouts, (tpos[0,0]+tpos[2,0])*0.5, tpos[3,0]+ychsz*0.5, /normal, $
    alignment = 0.5, 'PS (nT)!U2!N/Hz'
    
    
    
; ** waveform.
tpos = reform(poss[*,0,*])
tpos[0,*] = x0
tpos[2,*] = x1
ttvar = tvar+'_f'+bandids
options, ttvar, 'ytitle', ''
options, ttvar, 'yticks', 2
options, ttvar, 'yminor', 5
options, ttvar, 'labels', ''
tplot, ttvar, position = tpos, /noerase

for i = 0, 2 do begin
    xyouts, tpos[0,i]+xchsz*1, tpos[3,i]-ychsz*1, /normal, alignment = 0, $
        'a'+bandids[i]+'.'
    xyouts, tpos[0,i]-xchsz*9, (tpos[1,i]+tpos[3,i])*0.5, /normal, $
        alignment = 0, 'Wave '+bandids[i]+'!Cf'+bandids[i]+'='+$
        sgnum2str(round(wfs[i]))+' Hz'
endfor
xyouts, (tpos[0,0]+tpos[2,0])*0.5, tpos[3,0]+ychsz*0.5, /normal, $
    alignment = 0.5, 'dB FAC waveform (nT)'
    
    
; ** labels.
ty = tpos[3,1]-ychsz*0.85
tx = tpos[2,1]-xchsz*3.5
for i = 0, 2 do begin
    tty = ty-ychsz*i*0.7
    xyouts, tx, tty, /normal, alignment = 0, $
        fac[i], color = rgb[i], charsize = 0.8
    plots, tx+[-2,-0.5]*xchsz, tty+ychsz*0.2, /normal, color = rgb[i]
endfor

sgclose

end
