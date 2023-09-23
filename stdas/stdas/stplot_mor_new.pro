;+
; Type: procedure.
; Do not perform "scale bias" correction, info contains more useful info.
; Purpose: Wrap wavelet for quick spectrogram within tplot.
; Parameters: vname, in, type = string, required. Signal varname in tplot.
; Keywords:
;   period, in, boolean. Set to make ytitle to be period. Default is frequency.
;   scale_info, in, struct, opt. Default is
;       {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}.
;   tscale, in/out, dblarr[m], opt. Scales in time. Info used are
;       minscale, maxscale, nscale.
;   newname, in, string, optional. Output varname in tplot.
;   overwrite, in, boolean, optional. Set to overwrite vname.
;   ytitle, in, stirng, optional. Title for var.
; Notes: The variable should be 1-d data. This code forces label to be
;   specified by label keyword if vname contains no label.
; Dependence: tplot,slib.
; Author: Sheng Tian.
; History: 2017-07-25, Sheng Tian, create.
;   2023-09-22, Sheng Tian, change it to function, set frequency as default.
;-
function stplot_mor_new, vname, $
    period=period, $
    newname=newname, overwrite=overwrite, $
    tscale=tscales, scale_info=scale_info, $
    ytitle=ytitle, zrange=zrange

    get_data, vname, uts, f0, limits = lim
    if size(f0,/n_dimensions) eq 2 then f0 = snorm(f0)
    dr0 = sdatarate(uts) & dr1 = 1d/dr0
    nrec = n_elements(uts)
    dur = dr0*nrec
    ts = uts-uts[0]
    
    ; remove nan.
    idx = where(finite(f0,/nan), cnt)
    if cnt ne 0 then f0[idx] = 0

    
;----prepare scales.
    ; input: tscales, scale_info.
    ; output: s0,s1,dj,ns,scale_info.
    if n_elements(scale_info) eq 0 then $
        scale_info = {s0:4d*dr0, s1:0.5d*dur, dj:1d/8, ns:0d}

    if n_elements(tscales) eq 0 then begin          ; no scales.
        s0 = scale_info.s0
        s1 = scale_info.s1
        dj = scale_info.dj
        ns = scale_info.ns

        if s0 eq 0 then s0 = 4*dr0
        if s1 eq 0 then s1 = 0.5d*dur
        if dj eq 0 and ns eq 0 then dj = 1d/8
        if ns eq 0 then begin
            j1 = floor(alog(s1/s0)/alog(2)/dj)  ; # of powers-of-two with dj
            s1 = s0*2d^(dj*j1)            
            ns = j1+1
        endif
        if dj eq 0 then dj = alog(s1/s0)/alog(2)/(ns-1)
    endif else begin
        s0 = min(tscales)
        s1 = max(tscales)
        ns = n_elements(tscales)
        dj = alog(s1/s0)/alog(2)/(ns-1)
    endelse
    j1 = ns-1
    w0 = 6d
    cdelta = 0.776d     ; constant for w0=6, for normalization.

    scale_info.s0 = s0   ; min scale.
    scale_info.s1 = s1   ; max scale.
    scale_info.dj = dj   ; 2^dj scale spacing.
    scale_info.ns = ns   ; # of scales.



    if keyword_set(overwrite) then newname = vname
    if ~keyword_set(newname) then newname = vname+'_mor'
    ztitle = stagexist(lim,'ytitle')? lim.ytitle: ''

;----wavelet analysis.
    ; wavelet transform, in X.
    mor = wavelet(f0, dr0, /pad, s0=s0, dj=dj, j=j1, $
        mother='Morlet', param=w0, $
        period = ps, scale=ss, coi=coi)
    fs = 1d/ps
    
    ; power spectrogram, in X^2.
    pow = abs(mor)^2
    
    ; global wavelet spectrum, in X^2.
    gws = total(pow,1)/nrec
    
    ; power spectral density, in X^2/Hz.
    c_tau = mean(ps/ss)
    psd = 2*c_tau*dr0/cdelta*gws

    ; rectified, normalized power spectrogram, in X^2.
    npow = pow
    ;for i = 0, ns-1 do npow[*,i]/= ss[i]
    ;npow*= (dr0*dj/cdelta)
    
    ; rectified, normalized global wavelet spectrum, in, X^2.
    ngws = total(npow,1)/nrec
    
    ; significant level for global wavelet spectrum.
    dof = nrec-ss   ; the -scale corrects for padding at edges
    gsignif = wave_signif(f0, dr0, ss, 1, $
        lag1=0d, dof=dof, mother='Morlet')
    ngsignif = gsignif/ss*(dr0*dj/cdelta)
    
;----FFT.    
    ; Fourier transform, in X.
    ; boot definition is b_fft = fft*T, in X-s.
    ; current definition gives FFT indpendent of duration.
    fft = fft(f0)
    fft = fft[1:nrec/2] ; only use positive freq.
    
    ; Fourier power spectrum, in X^2.
    fftfps = 2*abs(fft)^2
    
    ; power spectral density, in X^2/Hz.
    T = ts[-1]-ts[0]
    fftpsd = fftfps*T

    ; freq spacing, freqs, periods.
    df = 1d/dur
    fftfs = (findgen(nrec)/(nrec*dr0))[1:nrec/2]
    fftps = 1d/fftfs
    
    sigma2 = (moment(f0))[1]
    
    if ~keyword_set(period) then coi = 1d/coi
    
;    sgopen, 0, xsize=8, ysize=5, /inch
;    plot, fftps, fftpsd, xlog=1
;    oplot, ps, gws/nrec, color=sgcolor('red')
;    oplot, ps, ngws, color=sgcolor('blue')

    ; info contain FFT and Morlet spectrum, etc.
    fftinfo = {$
        n:nrec, $       ; # of record.
        dt:dr0, $       ; delta t, in sec.
        cdelta:cdelta, $; constant.
        c_tau:c_tau, $  ; converts scale to period.
        sigma2:sigma2, $; variance.
        fft:fft, $      ; FFT, in X.
        fftps:fftps, $  ; periods, FFT periods, in sec.
        fftfs:fftfs, $  ; frequency, in Hz.
        fps:fftfps, $   ; Fourier power spectrum, in X^2.
        fftpsd:fftpsd, $; FFT PSD, in X^2/Hz.
        ps:ps, $        ; Morlet periods.
        fs:fs, $        ; Morlet frequency.
        ss:ss, $        ; Morlet scalea.
        dj:dj, $        ; 2^dj is spacing b/w scales.
        coi:coi, $      ; Morlet COI in period, or in Hz if frequency is set.
        gsignif:gsignif, $  ; Morlet significant level.
        ngsignif:ngsignif, $  ; Morlet normalized significant level.
        gws:gws, $      ; Morlet global wavelet spectrum, in X^2.
        psd:psd, $      ; Morlet power spectral density, in X^2/Hz.
        ngws:ngws, $    ; Morlet normalized GWS, in X^2.
        dof:dof}        ; degree of freedom.
    
    unit = ''
    var_unit = get_setting(vname, 'unit', exist)
    if exist then unit = var_unit+'!U2!N'
    if ~keyword_set(period) then begin
        val = fs
        yrange = minmax(1d/[s0,s1])
        ytitle = 'Freq!C(Hz)'
    endif else begin
        val = ps
        yrange = [s0,s1]
        ytitle = 'Period!C(sec)'
    endelse
    store_data, newname, uts, npow, val
    add_setting, newname, smart=1, dictionary($
        'spec', 1, $
        'ytitle', ytitle, $
        'yrange', yrange, $
        'ylog', 1, $
        'unit', unit, $
        'zlog', 1, $
        'color_table', 60 )
    store_data, newname+'_fft_info', 0, fftinfo
    return, newname
end

_2013_0607_load_data
pre0 = 'rbspa_'

var = pre0+'de_fac'
tvar = pre0+'de_mor'

;var = 'nino'
;tvar = 'nino2'


stplot_mor, var, newname=tvar


erase
device, decompose=1

get_data, var, uts, dat

get_data, tvar, uts, dat, val
get_data, tvar+'_fft_info', tmp, fftinfo

zr = [0,max(dat)*0.5]
ct = 40
top = 254
dat = bytscl(dat,min=zr[0],max=zr[1],top=top)
tpos = sgcalcpos() & tpos[2] = 0.7
sgtv, dat, position=tpos, ct=ct, /resize
ts = uts-uts[0]
yr = minmax(val)
plot, ts, val, /nodata, /noerase, position=tpos, $
    xstyle=1, xlog=0, xticklen=-0.01, xtitle='Time (sec)', $
    ystyle=1, ylog=1, yticklen=-0.01, yrange=yr, ytitle='Period!C(sec)'
oplot, ts, fftinfo.coi, color=sgcolor('white')

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size
tpos1 = [tpos[0],tpos[3]+ychsz*0.5,tpos[2],tpos[3]+ychsz*1]
sgcolorbar, findgen(top), /horizontal, zrange=zr, position=tpos1, ct=ct, $
    ztitle='Morlet |E|!U2!N (mV/m)!U2!N'

tpos2 = [tpos[2]+xchsz*1,tpos[1],tpos[2]+xchsz*1+(tpos[2]-tpos[0])*0.25,tpos[3]]
plot, fftinfo.ngws, fftinfo.ps, /noerase, position=tpos2, $
    xstyle=0, xlog=0, xticklen=-0.01, xticks=1, xtitl='PS (mV/m)!U2!N', $
    ystyle=1, ylog=1, yticklen=-0.01, yrange=yr, ytickformat='(A1)'
oplot, fftinfo.ngsignif, fftinfo.ps, linestyle=1
;oplot, fftinfo.psd, fftinfo.fftps, color=sgcolor('red')
;oplot, fftinfo.gws/n_elements(ts), fftinfo.ps, color=sgcolor('blue')


end
