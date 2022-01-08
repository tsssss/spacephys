
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


tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'version', 1



_2013_0607_0456_load_burst_data







; Fourier transform for the components of dE and dB.
vars = pre0+['mb1_uvw','eb1_uvw']
foreach tvar, vars do begin
    get_data, tvar, uts, dat
    
    tdr = sdatarate(uts)
    nrec = n_elements(uts)
    
    fs = findgen(nrec)/((nrec-1)*tdr)   ; freqs.
    fp = dblarr(nrec,3)                 ; fft power.
    for i = 0, 2 do begin
        tdat = dat[*,i]
;        tdat = tdat-smooth(tdat, 0.1/tdr, /edge_mirror)    ; smooth a little bit
        fp[*,i] = abs(fft(tdat))
    endfor
    store_data, tvar+'_fft', fs[0:nrec/2], fp[0:nrec/2,*]
endforeach


; plot the FFT power spectrums. 1 band peaks at 50 Hz, the other at 70 Hz.
ofn = 0
ofn = shomedir()+'/fig_mb1_fft.pdf'
sgopen, ofn, xsize = 11, ysize = 8.5, /inch

get_data, pre0+'mb1_uvw_fft', fs, fp
get_data, pre0+'mb1_uvw', uts, dat
tits = 'dB '+['U','V','W']
poss = sgcalcpos(3,2, ypad = 7)
for i = 0, 2 do begin
    plot, uts-uts[0], dat[*,i], position = poss[*,0,i]+[0,0,1,0]*0.1, /noerase, $
        xstyle = 1, xtitle = 'Sec from'+time_string(uts[0]), ytitle = '(nT)', title = tits[i]
    plot, fs, fp[*,i], xrange = [10,100], position = poss[*,1,i]+[1,0,0,0]*0.1, /noerase, $
        title = 'FFT power '+tits[i], xtitle = 'Frequency (Hz)', ytitle = 'PS (nT)!U2!N/Hz'
endfor
sgclose


ofn = 1
ofn = shomedir()+'/fig_eb1_fft.pdf'
sgopen, ofn, xsize = 11, ysize = 8.5, /inch

get_data, pre0+'eb1_uvw_fft', fs, fp
get_data, pre0+'eb1_uvw', uts, dat
tits = 'dE '+['U','V','W']
poss = sgcalcpos(3,2, ypad = 7)
for i = 0, 2 do begin
    plot, uts-uts[0], dat[*,i], position = poss[*,0,i]+[0,0,1,0]*0.1, /noerase, $
        xstyle = 1, xtitle = 'Sec from'+time_string(uts[0]), ytitle = '(mV/m)', title = tits[i]
    plot, fs, fp[*,i], xrange = [10,100], position = poss[*,1,i]+[1,0,0,0]*0.1, /noerase, $
        title = 'FFT power '+tits[i], xtitle = 'Frequency (Hz)', ytitle = 'PS (mV/m)!U2!N/Hz'
endfor
sgclose

stop



vars = pre0+[['u','v','w']+'_gsm']
muvw2gsm = dblarr(nrec,3,3)
for i = 0, 2 do begin
    get_data, vars[i], tuts, dat
    muvw2gsm[*,i,*] = sinterpol(dat,tuts,uts)
endfor
store_data, pre0+'uvw2gsm', uts, muvw2gsm



get_data, pre0+'pos_gsm', tuts, rgsm
rgsm = sinterpol(rgsm, tuts, uts)
get_data, pre0+'b_gsm', tuts, bgsm
bgsm = sinterpol(bgsm, tuts, uts)
rhat = sunitvec(rgsm)
bhat = sunitvec(bgsm)
phat = sunitvec(scross(rhat,bhat))
vhat = scross(bhat,phat)


vars = pre0+['mb1','eb1']
foreach tvar, vars do begin
    get_data, tvar+'_uvw', uts, dat, limits = lim
    ndim = (size(dat,/dimensions))[1]
    for i = 0, ndim-1 do begin
;        tdat = smooth(dat[*,i], tpad/dr0, /edge_truncate)
        tdat = smooth(dat[*,i], p0/dr0, /edge_truncate)
        dat[*,i] = dat[*,i]-tdat
;        dat[*,i] = dat[*,i]-dat[0,i]
    endfor
    store_data, tvar+'_uvw', uts, dat, limits = lim
    
    ; rotate to gsm.
    tdat = dat
    for i = 0, 2 do tdat[*,i] = $
        dat[*,0]*muvw2gsm[*,0,i]+dat[*,1]*muvw2gsm[*,1,i]+dat[*,2]*muvw2gsm[*,2,i]
    store_data, tvar+'_gsm', uts, tdat, limits = $
        {ytitle:lim.ytitle, colors:[6,4,2], labels:['x','y','z']}
    
    
    ; rotate to fac.
    dat = [[sdot(tdat,bhat)],[sdot(tdat,phat)],[sdot(tdat,vhat)]]
    store_data, tvar+'_fac', uts, tdat, limits = $
        {ytitle:lim.ytitle, colors:[6,4,2], labels:['b','p','v']}
    
    
endforeach

stop


twavpol, pre0+'mb1_fac'
options, pre0+'mb1_fac_*', 'yrange', [0,200]

sgopen, shomedir()+'/fig_mb1_fac_wave_pol.pdf', xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43
tplot, pre0+'mb1_fac*'
sgclose



twavpol, pre0+'eb1_fac'
options, pre0+'eb1_fac_*', 'yrange', [0,200] 

sgopen, shomedir()+'/fig_eb1_fac_wave_pol.pdf', xsize = 8.5, ysize = 11, /inch
device, decomposed = 0
loadct2, 43
tplot, pre0+'eb1_fac*'
sgclose





vars = pre0+['mb1_fac','eb1_fac','pf_fac']

get_data, vars[1], uts, de
get_data, vars[0], uts, db
pf = spoynt(de,db)
store_data, vars[2], uts, pf, limits = $
    {ytitle:'S!C(mW/m!U2!N)', colors:[6,4,2], labels:['b','p','v']}


ofn = shomedir()+'/fig_nice_waveform.pdf'
sgopen, ofn
device, decomposed = 0
loadct2, 43
tplot, vars
sgclose


sgopen, 0
device, decomposed = 0
loadct2, 43
tplot, vars


tmp = linfit(snorm(db)*1e3,snorm(de))

ofn = shomedir()+'/fig_ebratio.pdf'
sgopen, ofn, xsize = 5, ysize = 5, /inch
plot, snorm(db)*1e3, snorm(de), psym = 1, /iso, $
    xlog=1, ylog=1, xrange = [0.1,1000], yrange = [0.1,1000], $
    xtitle = 'dB x 1000 (nT)', ytitle = 'dE (mW/m)', $
    title = 'R(E/B) = '+sgnum2str(tmp[1]*1e6)+' km/s'
oplot, [0.1,1e3], [0.1,1e3]*tmp[1], color = sgcolor('red')
sgclose

end
