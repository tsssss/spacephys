


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.
;utr = time_double(['2013-06-07/05:00','2013-06-07/05:30'])
reload = 0


; 't01s' gives same result as 't04s'.
; the B0 model ~ measured. tried t04s, t96, t89, B0 are different.
models = ['t89','t04s','t01','t96']
nmodel = n_elements(models)
probes = ['a','b']
nprobe = n_elements(probes)


dr0 = 1d/16
site = 'pina'            ; the site under conjunction.
sites = ['pina']         ; sites for asf.
scsyms = [1,6]           ; symbols for rbsp-a and -b.
scclrs = [6,255]         ; colors for rbsp-a and -b.


datfn = shomedir()+'/psbl_de_32hz/2013_0607_data.tplot'
mapfn = shomedir()+'/psbl_de_32hz/2013_0607_map_data.tplot'
asifn = shomedir()+'/psbl_de_32hz/thg_asf_mosaic_2013_0607_0445.cdf'


vars = ['b0_gsm','pos_gse','fpt_'+['mlat','mlon','mlt'],'bmod_gsm_'+models]
mapvars = ['rbspa_'+vars, 'rbspb_'+vars]


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

load = 0
foreach tvar, mapvars do if tnames(tvar) eq '' then load = 1
if load then begin
    if file_test(mapfn) then begin
        load = 0
        tplot_restore, filename = mapfn
    endif else reload = 1
endif



; load position, map using different models.
if reload eq 1 then begin
    
    store_data, '*', /delete
    if file_test(mapfn) then file_delete, mapfn
    
    defsysv,'!rbsp_spice', exists = flag
    if flag eq 0 then rbsp_load_spice_kernels, trange = utr0, probes = probes

    tplot_restore, filename = datfn
    

    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        uts = smkarthm(utr0[0],utr0[1],3,'dx')
        get_data, pre0+'b0gse', tuts, dat
        dat = interpol(dat, tuts, uts)
        dat = sgse2gsm(stoepoch(uts,'unix'),dat)
        store_data, pre0+'b0_gsm', uts, dat, limits = $
            {ytitle:'B0 GSM!C(nT)', colors:[6,4,2], labels:'B GSM'+['x','y','z']}

        tvar = pre0+'pos_gse'
        rbsp_load_spice_state, probe = tprobe, coord = 'gse', times = uts, /no_spice_load
        get_data, pre0+'state_pos_gse', uts, posgse
        posgse = posgse*re1
        tdat = posgse
        store_data, pre0+'state_*', /delete
        
        store_data, tvar, uts, tdat, limits = $
            {ytitle:'R GSE (Re)', colors:[6,4,2], labels:'GSE '+['x','y','z']}
            
            
        tvar = pre0+'mlt'
        tdat = atan(posgse[*,1],posgse[*,0])*deg
        tdat = (tdat+360) mod 360   ; convert to 0-360.
        tdat = (tdat/15 + 12) mod 24
        store_data, tvar, uts, tdat, limits = {ytitle:'MLT (hr)'}
        
        
        tvar = pre0+'lshell'
        possm = sgse2sm(posgse, stoepoch(tuts,'unix'))
        mlat = atan(possm[*,2],sqrt(possm[*,0]^2+possm[*,1]^2)) ; in rad.
        dis = sqrt(possm[*,0]^2+possm[*,1]^2+possm[*,2]^2)
        tdat = dis/(cos(mlat)^2)
        store_data, tvar, uts, tdat, limits = {ytitle:'L-shell'}
        
        store_data, pre0+'dis', uts, dis


        ; map for each model.
        foreach model, models do begin
            suf0 = '_'+model
            ; map to 100 km, northern hemisphere.
            scalc_map_coef, pre0+'pos_gse', pre0+'b0_gsm', model = model, $
                coord = 'gse', /igrf, prefix = pre0, suffix = suf0, dir = -1
        endforeach
        
        vars = pre0+'fpt_'+['mlat','mlon','mlt']
        foreach tvar, vars do begin
            get_data, tvar+'_'+models[0], limits = lim
            stplot_merge, tvar+'_'+models, newname = tvar, limits = $
                {ytitle:lim.ytitle, labels:models, colors:findgen(nmodel)}
        endforeach
    endforeach
    
    tplot_save, mapvars, filename = mapfn

endif

stop





; **** load ASI.
load = 0
if n_elements(asifn) eq 0 then begin
    load = 1
endif else begin
    if file_test(asifn) eq 0 then load = 1
endelse

if load then begin
    tfn = sread_thg_mosaic(utr, sites, type = 'asf', minlat = 55, dark = 0)
    file_copy, tfn, asifn
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


img = bytarr(imgsz)
nrec = n_elements(uts)

txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
txs = txs/imgsz[0]*2
tys = tys/imgsz[0]*2
mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
mlts = atan(tys, txs)*deg+90    ; in deg.


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
    nvar = n_elements(vars)
    tutr = time_double(['2013-06-07/04:53:30','2013-06-07/04:58:30'])
    tutr = utr
    
    ofn = shomedir()+'/fig_keogram_'+models[modelidx]+'.pdf'
    sgopen, ofn, xsize = 8, ysize = 5, /inch
    device, decomposed = 0
    loadct2, 43
    
    poss = sgcalcpos(nvar)
    tplot, vars, position = poss, trange = tutr, title = 'Themis/ASI Keogram', vlab_margin = 10
    
    for i = 0, nvar-1 do begin
        get_data, 'rbsp'+probes[i]+'_fpt_mlat', tuts, dat
        dat = reform(dat[*,modelidx])
        plot, tutr, yr, /nodata, /noerase, position = poss[*,i], $
            xstyle = 5, ystyle = 5
        oplot, tuts, dat, color = 255
    endfor
    
    sgclose
endfor
stop


tpos = [0,0,1,1]
xr = [-90,0]
yr = [minlat,90]
yticks = [55,65,75]

for i = 0, nrec-1 do begin
    tut = uts[i]
    timg = img
    timg[pxidx] = mos[i,*]
    timg = timg[0:imgsz[0]/2,*]
    
    ofn = shomedir()+'/thg_asi2/thg_rb_'+ $
        time_string(tut,tformat='YYYY_MMDD_hhmm_ss')+'.png'
;    ofn = 0
    sgopen, ofn, xsize = 6, ysize = 6, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    device, decomposed = 0
    loadct2, 1
    sgtv, timg, position = tpos


    loadct2, 43
    sgset_map, position = tpos, color = 255, xrange = xr, $
        yrange = yr, ytickv = yticks
    xyouts, tpos[0]+xchsz, tpos[1]+ychsz, /normal, $
        time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
        
    
    for j = 0, nprobe-1 do begin
        tsym = scsyms[j]
        
        xyouts, tpos[0]+xchsz, tpos[1]+ychsz*(j+2), /normal, $
            'RBSP-'+strupcase(probes[j]), color = 255
        plots, tpos[0]+xchsz*7, tpos[1]+ychsz*(j+2.3), /normal, $
            psym = tsym, color = 255, symsize = 0.8

        pre0 = 'rbsp'+probes[j]+'_'
        get_data, pre0+'fpt_mlt', tuts, mlts
        get_data, pre0+'fpt_mlat', tuts, mlats
        
        tmlats = sinterpol(mlats, tuts, tut, /quadratic)
        tmlts = sinterpol(mlts, tuts, tut, /quadratic)
        
        clrs = 6-indgen(nmodel)
        for k = 0, nmodel-1 do begin
            plots, tmlts[k]*15, tmlats[k], /data, psym = tsym, color = clrs[k]
            xyouts, tpos[0]+xchsz, tpos[1]+ychsz*(k+4), /normal, $
                strupcase(models[k]), color = clrs[k]
        endfor
    endfor
    sgclose
endfor


end
