


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.
;utr = time_double(['2013-06-07/05:00','2013-06-07/05:30'])
reload = 0



datfn = shomedir()+'/psbl_de_32hz/2013_0607_data.tplot'
mapfn = shomedir()+'/psbl_de_32hz/2013_0607_map_data.tplot'
asifn = shomedir()+'/psbl_de_32hz/thg_asf_mosaic_2013_0607_0445.cdf'



device, decomposed = 0
loadct2, 43
tplot_options, 'ynozero', 1
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1




; **** load ASI.
load = 0
if n_elements(asifn) eq 0 then begin
    load = 1
endif else begin
    if file_test(asifn) eq 0 then load = 1
endelse

;load = 1
if load then begin
    sites = ['pina']
    tfn = sread_thg_mosaic(utr, sites, type = 'asf', minlat = 55, dark = 0, notop = 1)
    file_copy, tfn, asifn, /allow_same, /overwrite
    file_delete, tfn
endif

cdfs = scdfread(asifn)

uts = (*cdfs[0].value)
mos = (*cdfs[1].value)
midn= (*cdfs[2].value)
mlt = (*cdfs[3].value)
mlat= (*cdfs[4].value)
imgsz = (*cdfs[5].value)
pxidx = (*cdfs[6].value)
minlat = (*cdfs[7].value)[0]

_2013_0607_load_data



img = intarr(imgsz)
nrec = n_elements(uts)

txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
txs = txs/imgsz[0]*2
tys = tys/imgsz[0]*2
mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
mlts = atan(tys, txs)*deg+90    ; in deg.


;mlts = (mlts/15+24) mod 24
utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
idx = where(uts ge utr[0] and uts le utr[1], nrec)
uts = uts[idx]
mos = mos[idx,*]





; keogram at different mlts.
models = ['t89','t04s','t01','t96']
nmodel = n_elements(models)
probes = ['a','b']
nprobe = n_elements(probes)


yr = [60,70]
mlat0s = smkarthm(yr[0],yr[1],0.1,'dx')
nmlat0 = n_elements(mlat0s)

for modelidx = 0, nmodel-1 do begin

    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        get_data, pre0+'fpt_mlat', tuts, scmlats
        get_data, pre0+'fpt_mlt', tuts, scmlts
        
        scmlats = reform(scmlats[*,modelidx])
        scmlts = reform(scmlts[*,modelidx])
        
        scmlats = interpol(scmlats,tuts,uts)
        scmlts = interpol(scmlts,tuts,uts)*15-360
        
        keos = fltarr(nrec,nmlat0)
        for i = 0, nrec-1 do begin
            idx = where(abs(mlts-scmlts[i]) le 0.1, cnt)
            if cnt lt 3 then continue
            timg = img
            timg[pxidx] = mos[i,*]
            
            timg = timg[idx]
            tmlats = mlats[idx]
            idx = sort(tmlats)
            timg = timg[idx]
            tmlats = tmlats[idx]
            
            keos[i,*] = interpol(timg,tmlats,mlat0s)
        endfor
        
        tvar = (pre0+'keogram_'+models[modelidx])[0]
        store_data, tvar, uts, keos, mlat0s, $
            limits = {spec:1, no_interp:1, $
            ytitle:'RBSP-'+strupcase(tprobe)+'!CFptMLat!C(deg)', $
            yrange:yr, ztitle:'Photon count'}
    endforeach
    
    
    vars = ['rbspa_','rbspb_']+'keogram_'+models[modelidx]
;    vars = [vars,['rbspa_','rbspb_']+'pf_fac_mat']
    nvar = n_elements(vars)
    tutr = time_double(['2013-06-07/04:53:30','2013-06-07/04:58:30'])
    tutr = utr
    
    ofn = shomedir()+'/fig_keogram_'+models[modelidx]+'.pdf'
    sgopen, ofn, xsize = 8, ysize = 5, /inch
    device, decomposed = 0
    loadct2, 43
    
    poss = sgcalcpos(nvar)
    tplot, vars, position = poss, trange = tutr, title = 'Themis/ASI Keogram', vlab_margin = 10
    
    for i = 0, nprobe-1 do begin
        get_data, 'rbsp'+probes[i]+'_fpt_mlat', tuts, dat
        dat = reform(dat[*,modelidx])
        plot, tutr, yr, /nodata, /noerase, position = poss[*,i], $
            xstyle = 5, ystyle = 5
        oplot, tuts, dat, color = 255
    endfor
    
    sgclose
endfor
stop






; search for points in this region.
yr = [62,67]    ; in mlat in deg.
xr = [21.6,22.6]    ; in mlt in hr.

yr = 66.0+[-1,1]*0.1    ; in mlat in deg.
xr = 21.7+[-1,1]*0.5    ; in mlt in hr.

pts = where(mlats ge yr[0] and mlats le yr[1] and mlts ge xr[0] and mlts le xr[1], npt)

photoncnts = dblarr(nrec, npt)
for i = 0, nrec-1 do begin
    timg = img
    timg[pxidx] = mos[i,*]
    for j = 0, npt-1 do photoncnts[i,j] = timg[pts[j]]
endfor


vars = ['photoncnt','rbspb_pf']
get_data, 'rbspb_pf_fac_mat3', tuts, pflux
v1 = sinterpol(pflux[*,0], tuts, uts)
store_data, vars[1], uts, v1

corrs = dblarr(npt)
for i = 0, npt-1 do begin
;    store_data, vars[0], uts, photoncnts[*,i]
    corrs[i] = c_correlate(v1, photoncnts[*,i], 0)
;    tplot, vars, trange = utr
;    xyouts, 0.5, 0.95, /normal, alignment = 0.5, $
;        'MLat: '+sgnum2str(mlats[pts[i]])+', MLT: '+sgnum2str(mlts[pts[i]])
;    stop
endfor

idx = reverse(sort(corrs))
corrs = corrs[idx]
pts = pts[idx]
mlt1s = mlts[pts]
mlat1s = mlats[pts]
photoncnts = photoncnts[*,idx]

maxcorr = max(corrs,/nan,idx)

for i = idx[0], npt-1 do begin
    v2 = photoncnts[*,i]
;    v2 = v2-smooth(v2,20,/edge_mirror)
    store_data, vars[0], uts, v2
    tplot, vars, trange = utr
    xyouts, 0.5, 0.95, /normal, alignment = 0.5, $
        'Corr: '+sgnum2str(corrs[i])+', MLat: '+sgnum2str(mlats[pts[i]])+', MLT: '+sgnum2str(mlts[pts[i]])
    stop
endfor




stop




end