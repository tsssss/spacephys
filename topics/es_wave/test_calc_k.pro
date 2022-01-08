
; settings.
utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.36',':02.41'])  ; the largest sub-wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.35',':02.40'])  ; the largest sub-wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.35',':02.40'])  ; the largest sub-wave packet.
utr3 = time_double('2013-06-07/04:56'+[':02.40',':02.47'])  ; the largest sub-wave packet.



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
lmn = ['max','med','min']
wfs = [24d,48,74]   ; wave frequency.


tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'version', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'yticklen', 1



_2013_0607_0456_load_burst_data



get_data, 'morlet_info', tmp, info
filters = info.filters
nband = n_elements(filters)-1
bandids = string(findgen(nband)+1,format='(I0)')
periods = info.periods


tinfo = {l:dblarr(3),m:dblarr(3),n:dblarr(3),eigenvals:dblarr(3)}
lmninfos = replicate(tinfo,nband)


kvecs = dblarr(3,nband)

for i = 0, nband-1 do begin
    
    tutr = (i eq 0)? utr3: utr2
    w0 = 2*!dpi*wfs[i]  ; angular freq, in rad/Hz.
    p0 = 1d/wfs[i]      ; period, in sec.

    
    ; calc rotation matrix from certain coord to LMN.
    tvar = pre0+'mb1_fac_f'+bandids[i]
    get_data, tvar, uts, dat
    utidx = where(uts ge tutr[0] and uts le tutr[1], tnrec)
    tuts = uts[utidx]
    tdat = dat[utidx,*]
    
    eigvs = smva(tdat, rotm)
    dblmn = tdat # rotm  ; [max,median,min].


    ; rotate E into LMN.
    tvar = pre0+'eb1_fac_f'+bandids[i]
    get_data, tvar, uts, dat
    tuts = uts[utidx]
    tdat = dat[utidx,*]
    
    delmn = tdat # rotm
    
    
    ; step 1: determine k hat.
    kfac = rotm[*,2]
    khat = [0d,0,1]
    
    
    ; step 2: calculate relative phase.
    ntp = 50
    tps = smkarthm(0,p0,ntp,'n')
    tns = tps/dr0
    
;    ephas = dblarr(3)      ; relative phase, in deg.
;    dat = delmn
;    for i = 1, 2 do begin
;        cor = dblarr(ntp)
;        for j = 0, ntp-1 do cor[j] = c_correlate(dat[*,0],dat[*,i],-tns[j])
;        maxcor = max(cor, pidx)
;        ephas[i] = tps[pidx]*w0*deg
;    endfor
;    
;    
;    bphas = dblarr(3)      ; relative phase, in deg.
;    dat = dblmn
;    for i = 1, 2 do begin
;        cor = dblarr(ntp)
;        for j = 0, ntp-1 do cor[j] = c_correlate(dat[*,0],dat[*,i],-tns[j])
;        maxcor = max(cor, pidx)
;        bphas[i] = tps[pidx]*w0*deg
;    endfor
    
    
    ; step 3: calculate wave amplitude.
    dat = delmn
    eamps = dblarr(tnrec,3)
    for j = 0, 2 do begin
        tdat = deriv(dat[*,j])/(dr0*w0)
        eamps[*,j] = sqrt(tdat^2+dat[*,j]^2)
        eamps[*,j] = smooth(eamps[*,j], p0/dr0, /edge_truncate)
    endfor
    
    
    dat = dblmn
    bamps = dblarr(tnrec,3)
    for j = 0, 2 do begin
        tdat = deriv(dat[*,j])/(dr0*w0)
        bamps[*,j] = sqrt(tdat^2+dat[*,j]^2)
        bamps[*,j] = smooth(bamps[*,j], p0/dr0, /edge_truncate)
    endfor
    bamps[*,2] = 0  ; assume perfect hodogram.

    
    
    
    ; calculate |k|.
    ks = dblarr(tnrec)
    as = dblarr(tnrec)
    for j = 0, tnrec-1 do begin
        as[j] = sang(khat,reform(eamps[j,*]))
        kxe = scross(khat,reform(eamps[j,*]))
        wdb = w0*reform(bamps[j,*])
        tmp = linfit(abs(kxe), abs(wdb))
        ks[j] = tmp[1]
;        plot, abs(kxe), abs(wdb), psym = -1
;        stop
    endfor
    as = as*deg
    
    kvecs[*,i] = kfac*ks[tnrec/2]
    
    
    plot, snorm(eamps[*,0:1]), snorm(bamps[*,0:1]), psym = 1
    tmp = linfit(snorm(eamps[*,0:1]), snorm(bamps[*,0:1])*w0)
    k0 = tmp[1]*1e-3    ; |k| in rad/km.
    print, kfac
    print, k0
    print, w0/k0        ; in km/s.
    
    stop
    
    erase
    poss = sgcalcpos(2)
    plot, as, position = poss[*,0], /noerase
    plot, ks, position = poss[*,1], /noerase
    
endfor

stop

end
