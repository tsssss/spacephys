; Vsc, HOPE density, S, keogram, 
; electron and ion energy spectrogram and MAGEIS electron energy spectrogram.


; var1 and var2 have to be on the same uniform times.
pro tplot_cross_corr, var1, var2, time1=utr1, time2=utr2

    if n_elements(utr1) eq 0 then message, 'set var1 range...'

    get_data, var1, uts, dat1
    get_data, var2, uts, dat2
    dr0 = sdatarate(uts)
    
    ; trancate data2.
    if n_elements(utr2) eq 2 then begin
        idx = where(uts ge utr2[0] and uts le utr2[1])
        uts = uts[idx]
        dat1 = dat1[idx]
        dat2 = dat2[idx]
    endif
    
    ; calc derivative.
;    dat1 = deriv(dat1)
;    dat2 = deriv(dat2)

    ; idx1, idx2: start and end index.
    tmp = min(uts-utr1[0],/absolute, idx1)
    tmp = min(uts-utr1[1],/absolute, idx2)
    dat1 = dat1[idx1:idx2]
    
    len1 = idx2-idx1+1
    len2 = n_elements(uts)
    ncorr = len2-len1
    corrs = dblarr(ncorr)
    
    dat1 = dat1-mean(dat1)
    dat2 = dat2-mean(dat2)
    for i=0, ncorr-1 do begin
        corrs[i] = c_correlate(dat1, dat2[i:i+len1-1], 0)
    endfor
    
    tmp = max(corrs, idx)
    print, 'max corr = ', corrs[idx]
    print, 'dT       = ', (idx-idx1)*dr0+' (sec)'
    print, 'vel      = ', 0.5*6400/(idx-idx1)/dr0+' (km/s)'    

    plot, dat2[idx:idx+len1-1], color=6, /ynozero
    plot, dat1, /ynozero, /noerase
    stop
end


_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'yticklen', -0.01
tplot_options, 'xticklen', -0.02
tplot_options, 'zcharsize', 0.8
tplot_options, 'ycharsize', 1
tplot_options, 'ystyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1




models = ['t89','t04s','t01','t96']
modelidx = 2
probes = ['a','b']
nprobe = n_elements(probes)

diputa = time_double('2013-06-07/04:54:30')
diputb = time_double('2013-06-07/04:54:59')

; make gsm.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'b_gse', uts, bgse
    bgsm = sgse2gsm(bgse, stoepoch(uts,'unix'))
    bcomb = bgsm
    ndim = 3
    b0s = dblarr(ndim)
    b0s = [200,-100,0,180]
    b0s = round(bcomb[0,*]/5)*5
    labs = ['B!D'+['X','Y','Z']]
    for i = 0, ndim-1 do begin
        bcomb[*,i] = bcomb[*,i]-b0s[i]
        if b0s[i] eq 0 then continue
        if b0s[i] gt 0 then begin
            labs[i] = labs[i]+'!N-'+sgnum2str(b0s[i])+' nT'
        endif else begin
            labs[i] = labs[i]+'!N+'+sgnum2str(-b0s[i])+' nT'
        endelse
    endfor
    
    store_data, pre0+'b_combine', uts, bcomb, limits = $
        {ytitle:'(nT)', colors:[6,4,2,0], labels:labs, $
        yrange:[-20,60], ystyle:1, yticks:2, yminor:5}
endforeach


yr = [60,68]
mlat0s = smkarthm(yr[0],yr[1],0.1,'dx')
nmlat0 = n_elements(mlat0s)

get_data, 'asf_info', tmp, info
mlats = info.mlats
mlts = info.mlts
imgsz = info.imgsz
get_data, 'asf_mos', uts, mos, pxidx
nrec = n_elements(uts)

; make keogram.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'fpt_mlat', tuts, scmlats
    get_data, pre0+'fpt_mlt', tuts, scmlts
    
    scmlats = reform(scmlats[*,modelidx])
    scmlts = reform(scmlts[*,modelidx])
    
    scmlats = interpol(scmlats,tuts,uts)
    scmlts = interpol(scmlts,tuts,uts)*15-360
    
    keos = fltarr(nrec,nmlat0)
    cnts = fltarr(nrec)
    for i = 0, nrec-1 do begin
        idx = where(abs(mlts-scmlts[i]) le 0.1, cnt)
        if cnt lt 3 then continue
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        
        timg = timg[idx]
        tmlats = mlats[idx]
        idx = sort(tmlats)
        timg = timg[idx]
        tmlats = tmlats[idx]
        
        keos[i,*] = interpol(timg,tmlats,mlat0s)
        
        tmp = min(abs(mlat0s-scmlats[i]), idx)
        cnts[i] = keos[i,idx]
        keos[i,idx] = !values.f_nan
    endfor
    
    tvar = (pre0+'keogram_'+models[modelidx])[0]
    store_data, tvar, uts, keos, mlat0s, $
        limits = {spec:1, no_interp:1, $
        ytitle:'MLat!C(deg)', $
        yrange:yr, ztitle:'Photon count', yticks:2, zrange:[50,300], zticks:1, ystyle:1}
    tvar = (pre0+'count_'+models[modelidx])[0]
    store_data, tvar, uts, cnts, $
        limits = {ytitle:'Photon count',yrange:[100,500], yminor:5, labels:'Photon!C  @SC footpnt'}
endforeach


; calc mapped pflux.
pflabs = 'S!D'+['||','!9^!X,West','!9^!X,North']

foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'pf_fac', uts, dat, limits = lim
    get_data, pre0+'map_coef', tuts, tmp
    tmp = interpol(tmp[*,modelidx],tuts, uts)
    for i = 0, 2 do dat[*,i] *= tmp
    store_data, pre0+'pf_fac_mor_map', uts, dat, limits = lim
    options, pre0+'pf_fac_mor_map', 'labels', pflabs
    
    store_data, pre0+'pf_fac_mor_map', uts, dat[*,0], limits = $
        {ytitle:'(mW/m!U2)', labels:pflabs[0]+'!N @100km', colors:0}
endforeach




vars = tnames('rbsp?_e_en')
options, vars, 'ytitle', '(eV)'
options, vars, 'zrange', [1e6,1e9]
options, vars, 'ytickname', '1e'+['2','3','4']
options, vars, 'ztickname', '1e'+['6','7','8','9']

vars = tnames('rbsp?_b_combine')
options, vars, 'constant', !values.d_nan

vars = tnames('rbsp?_n_combine')
options, vars, 'ytitle', '(cm!U-3!N)'

vars = tnames('rbsp?_pf_fac_mor_map')
options, vars, 'ytitle', '(mW/m!U2!N)'
options, vars, 'labels', 'Map to!C100 km'
options, vars, 'constant', !values.d_nan

vars = tnames('rbsp?_e_en')
options, vars, 'ytitle', '(eV)'
options, vars, 'yrange', [200,2e4]
options, vars, 'ytickv', [1000,10000]
options, vars, 'ytickname', '1e'+['3','4']
options, vars, 'zrange', [1e6,1e9]
options, vars, 'zticks', 3
options, vars, 'ztickname', '1e'+['6','7','8','9']

vars = tnames('rbsp?_o_en')
options, vars, 'ytitle', '(eV)'
options, vars, 'yrange', [40,4e4]
options, vars, 'ytickv', [100,1000,10000]
options, vars, 'ystyle', 1
options, vars, 'yminor', 5
options, vars, 'yticks', 2
options, vars, 'ytickname', '1e'+['2','3','4']
options, vars, 'zrange', [1e4,1e7]
options, vars, 'zticks', 3
options, vars, 'ztickname', '1e'+['4','5','6','7']

vars = 'rbsp?_keogram_t01'
options, vars, 'ytitle', 'MLat (deg)'
options, vars, 'yminor', 4


foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'

    ; |B|.
    get_data, pre0+'b_gse', uts, dat
    tmp = sqrt(dat[*,0]^2+dat[*,1]^2+dat[*,2]^2)
    ; smooth to remove spin modulation.
    dr0 = sdatarate(uts)
    tmp = smooth(tmp, 11d/dr0)
    store_data, pre0+'bmag', uts, tmp

    ; angle out of x-y plane.
    tmp = atan(dat[*,2],sqrt(dat[*,0]^2+dat[*,1]^2))*deg
    ;tmp = smooth(tmp, 11d/dr0)
    store_data, pre0+'db_angle', uts, tmp
    
    ; bz.
    store_data, pre0+'bz', uts, tmp
endforeach
options, 'rbsp?_bmag', 'ynozero', 1
options, 'rbspa_bmag', 'yrange', [236.5,244.5]
options, 'rbspb_bmag', 'yrange', [222.5,230.5]
options, 'rbsp?_bmag', 'ytitle', '(nT)'
options, 'rbsp?_db_angle', 'ytitle', '(deg)'
options, 'rbspa_db_angle', 'yrange', [0,10]
options, 'rbspb_db_angle', 'yrange', [2,12]
options, 'rbsp?_n_combine', 'ytitle', '(cm!U-3!N)'
options, 'rbsp?_bmag', 'labels', 'smoothed!C  @11 sec'


utr0 = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr1 = time_double(['2013-06-07/04:52','2013-06-07/05:02'])


vars = ['e_en','n_combine','b_combine','keogram_t01','pf_fac_mor_map','o_en']
logidx = 4
nvar = n_elements(vars)

labs = ['mlt','lshell','mlat']

tplot_cross_corr, 'rbspa_bmag', 'rbspb_bmag', $
    time1=time_double(['2013-06-07/04:54','2013-06-07/04:56']), $
    time2=time_double(['2013-06-07/04:52','2013-06-07/05:02'])
tplot_cross_corr, 'rbspa_bz', 'rbspb_bz', $
    time1=time_double(['2013-06-07/04:54','2013-06-07/04:55:10']), $
    time2=time_double(['2013-06-07/04:52','2013-06-07/05:02'])


ofn = 0
ofn = shomedir()+'/fig_psbl_b_timelag.pdf'
sgopen, ofn, xsize = 6, ysize = 8, /inch
device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

lmarg = 15
cols = [6,4,2]
pos1 = [0.1,0.5,1,1]+[0,1,0,1]*ychsz+[0,-1,0,1]*ychsz*1
pos2 = [0.1,0,1,0.5]+[0,1,0,1]*ychsz+[0,-1,0,1]*ychsz*1

vars = ['n_combine','db_angle','bmag']
nvar = n_elements(vars)
labs = ['a. Density', 'b. Ang.Bz/Bxy', 'c. |B|']

poss = sgcalcpos(nvar, region=pos1)
tplot, 'rbspa_'+vars, position=poss, trange=utr1, /noerase, title='RBSP-A'
for i=0, nvar-1 do xyouts, poss[0,i]-xchsz*15, poss[3,i]-ychsz*1, labs[i], /normal


poss = sgcalcpos(nvar, region=pos2)
tplot, 'rbspb_'+vars, position=poss, trange=utr1, /noerase, title='RBSP-B'
for i=0, nvar-1 do xyouts, poss[0,i]-xchsz*15, poss[3,i]-ychsz*1, labs[i], /normal


sgclose

end
