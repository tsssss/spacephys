
; plot vsc and bmag.

_2013_0607_load_data


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:05'])
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
probes = ['a','b']
spinrate = 11




device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 8
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1



; calc b mag.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'b_gse', uts, bgse
    bmag = snorm(bgse)
    dr = sdatarate(uts)
    bmag = smooth(bmag,11/dr,/edge_truncate)
    store_data, pre0+'bmag', uts, bmag, limits = $
        {ytitle:'|B|!C(nT)', ynozero:1, labels:'|B| smoothed!C  at spin rate'}
endforeach


; scale Vsc to density.
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


; break pos_gsm into pos_gsm_[xyz].
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'pos_gsm', uts, posgsm
    store_data, pre0+'pos_gsm_z', uts, posgsm[*,2], limits = $
        {ytitle:'Z GSM'}
endforeach



tvar = 'rbspa_n_combine'
options, tvar, 'ytitle', 'Density!C(cm!U-3!N)'

tvar = 'rbspb_n_combine'
options, tvar, 'ytitle', 'Density!C(cm!U-3!N)'

tvar = 'rbspa_bmag'
options, tvar, 'yrange', [234,244]
options, tvar, 'yrange', [226,250]
options, tvar, 'yticks', 2
options, tvar, 'yminor', 6
options, tvar, 'ystyle', 1
options, tvar, 'ytitle', '|B|!C(nT)'


tvar = 'rbspb_bmag'
options, tvar, 'yrange', [220,230]
options, tvar, 'yrange', [210,234]
options, tvar, 'yticks', 2
options, tvar, 'yminor', 6
options, tvar, 'ystyle', 1
options, tvar, 'ytitle', '|B|!C(nT)'

tvar = 'rbspa_db_fac'
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5
options, tvar, 'yrange', [-40,40]
options, tvar, 'ystyle', 1


tvar = 'rbspb_db_fac'
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5
options, tvar, 'yrange', [-40,40]
options, tvar, 'ystyle', 1


tvar = 'rbspa_de_fac'
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5
options, tvar, 'yrange', [-1,1]*150
options, tvar, 'ystyle', 1

tvar = 'rbspb_de_fac'
options, tvar, 'yticks', 2
options, tvar, 'yminor', 5
options, tvar, 'yrange', [-1,1]*150
options, tvar, 'ystyle', 1


tvar = ['de_fac','db_fac']
tvar = ['rbspa_'+tvar,'rbspb_'+tvar]
options, tvar, 'labels', ['parallel','west','north']

tvar = 'rbspa_e_en'
options, tvar, 'ytitle', 'Ele energy!C(eV)'

tvar = 'rbspb_e_en'
options, tvar, 'ytitle', 'Ele energy!C(eV)'

ofn = 0
ofn = shomedir()+'/fig_rbsp_overview.pdf'
sgopen, ofn, xsize = 8.5, ysize = 11, /inch
;sgopen, ofn, xsize = 11, ysize = 8.5, /inch
device, decomposed = 0
loadct2, 43

!x.gridstyle = 1
!x.ticklen = 1

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

labs = ['pos_gsm_z','lshell','mlt']
tvar = ['e_en','n_combine','bmag','db_fac','de_fac']

vars = ['rbspa_'+tvar]
nvar = n_elements(vars)
poss = sgcalcpos(nvar, region = [0,0.5,1,1])
;poss = sgcalcpos(nvar, region = [0,0,0.5,1])
tplot, vars, trange = utr, position = poss, vlab_marg = 12, /noerase, $
    title = 'RBSP-A', var_label = 'rbspa_'+labs

tx = xchsz*4
tys = reform(poss[3,*])-ychsz
figlabs = ['a','b','c','d','e']+'-1.'
for i = 0, nvar-1 do xyouts, tx, tys[i], /normal, alignment = 0.5, figlabs[i]

vars = ['rbspb_'+tvar]
nvar = n_elements(vars)
poss = sgcalcpos(nvar, region = [0,0,1,0.5]+[0,1,0,1]*ychsz)

;poss = sgcalcpos(nvar, region = [0.5,0,1,1])
tplot, vars, trange = utr, position = poss, vlab_marg = 12, /noerase, $
    title = 'RBSP-B', var_label = 'rbspb_'+labs

tx = xchsz*4
tys = reform(poss[3,*])-ychsz
figlabs = ['a','b','c','d','e']+'-2.'
for i = 0, nvar-1 do xyouts, tx, tys[i], /normal, alignment = 0.5, figlabs[i]


!x.gridstyle = 0
!x.ticklen = 0

sgclose

end
