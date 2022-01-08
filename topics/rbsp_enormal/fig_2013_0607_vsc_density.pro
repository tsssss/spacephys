
; plot vsc correlation.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
probes = ['a','b']



device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 5
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'xgridstyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1


foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'vsc', uts, vsc
    get_data, pre0+'n', tuts, hopen
    idx = where(finite(hopen))
    tuts = tuts[idx]
    hopen = hopen[idx]
    hopen = sinterpol(hopen, tuts, uts, /nan)
    
    
    ; hopen = a*exp(vsc*b)
    ; alog(hopen) = alog(a)+vsc*b
    
    tx0 = vsc
    ty0 = alog(hopen)
        
    a = regress(tx0, ty0, const = b)
    ty1 = a[0]*vsc+b[0]
    
    ybar = mean(ty0)
    sstot = total((ty0-ybar)^2)
    ssres = total((ty0-ty1)^2)
    
    r2 = 1-ssres/sstot
    print, a, b, r2
    

    efwn = exp(ty1)
    store_data, pre0+'n_efw', uts, efwn, limits = $
        {ytitle:'Density!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
        labels:'n!De,EFW Vsc'}
    store_data, pre0+'n_combine', uts, [[efwn],[hopen]], limits = $
        {ytitle:'RBSP-'+strupcase(tprobe)+'!CDensity!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
        labels:['n!De,EFW Vsc!N', 'n!De,HOPE > 200 eV'], colors:[0,6]}
endforeach

ofn = shomedir()+'/fig_vsc_density.pdf'
sgopen, ofn, xsize = 6, ysize = 3, /inch

device, decomposed = 0
loadct2, 43

vars = ['rbspa','rbspb']+'_n_combine'
nvar = n_elements(vars)
poss = sgcalcpos(nvar)

tplot, vars, position = poss, vlab_margin = 10, trange = utr
sgclose


end