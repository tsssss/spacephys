;+
; load ASI data.
;-

pro test_asi_donovan_load_data, utr, reload = reload

; settings.
sites = ['gill','snkq']
exclude = ['atha','tpas','fsmi','fsim','kapu','gbay','kuuj','chbg']
dr0 = 3d

model0 = 't01'
model0 = 't89'

if ~keyword_set(reload) then reload = 0


; prepare var names.
allvars = []

asivars = 'asf_'+['info','mos']
allvars = [allvars,asivars]

;mapvars = []
;vars = ['grid_vertical','grid_equatorial']
;mapvars = model0+'_'+vars
;allvars = [allvars,mapvars]




device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1
dir = -1        ; trace to northern hemisphere.



datdir = shomedir()+'/psbl_asi'
datfn = datdir+'/psbl_asi_'+time_string(utr[0],tformat='YYYY_MMDD_hh')+'.tplot'
load = 0
foreach tvar, allvars do if tnames(tvar) eq '' then load = 1
if load then if file_test(datfn) ne 0 then begin
    store_data, '*', /delete
    tplot_restore, filename = datfn
endif



; **** load ASI.
load = 0
foreach tvar, asivars do if tnames(tvar) eq '' then load = 1
if reload then load = 1
if load then begin
    asifn = sread_thg_mosaic(utr, sites, type = 'asf', minlat = 55, dark = 0, /notop)
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
    
    store_data, 'asf_mos', uts, mos, pxidx
    store_data, 'asf_info', 0, {imgsz:imgsz, minlat:minlat, mlats:mlats, mlts:mlts}
    
    ; update data file.
    tplot_save, allvars, filename = datfn
endif



;load = 0
;foreach tvar, mapvars do if tnames(tvar) eq '' then load = 1
;if reload then load = 1


;;load = 1
;if load then begin
;    
;    ; prepare model parameters.
;    foreach model0, models do sgeopack_par, utr, model0
;
;    
;    ; use the start, middle and end times.
;    uts = smkarthm(utr[0], utr[1], dr0, 'dx')
;    nrec = n_elements(uts)
;    utidx = fix([0,(nrec-1)/2,nrec-1])
;    uts = uts[utidx]
;    nrec = n_elements(uts)
;    
;    
;    ; ** map the vertical grid to ionosphere.
;    nsct = 60+1   ; # of grid in azimuthal direction.
;    nscz = 10+1   ; # of grid in z direction.
;    
;    scr0 = 5d     ; the L of vertical plane.
;    sctrng = [2.8,-1.7] ; 4.5 MLT wide.
;    sczrng = [1.5,4.5]
;    scts = smkarthm(sctrng[0],sctrng[1], nsct, 'n')*15*rad
;    sczs = smkarthm(sczrng[0],sczrng[1], nscz, 'n')
;    scxs = -scr0*cos(scts)
;    scys =  scr0*sin(scts)
;    
;    scfmlats = dblarr(nrec,nsct,nscz)
;    scfmlons = dblarr(nrec,nsct,nscz)
;    scfmlts = dblarr(nrec,nsct,nscz)
;
;    
;    foreach tmodel, models do begin
;        get_data, tmodel+'_par', tmp, dat
;        pars = sinterpol(dat, tmp, uts)
;        ets = stoepoch(uts,'unix')
;        
;        t89 = (tmodel eq 't89')? 1: 0
;        t96 = (tmodel eq 't96')? 1: 0
;        t01 = (tmodel eq 't01')? 1: 0
;        t04s = (tmodel eq 't04s')? 1: 0
;        storm = (tmodel eq 't04s')? 1: 0
;
;        ; map in situ grid to ionosphere.
;        for i = 0, nrec-1 do begin
;            tet = ets[i]
;            geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
;            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
;            
;            par = reform(pars[i,*])
;            for j = 0, nsct-1 do begin
;                for k = 0, nscz-1 do begin
;                    xp = scxs[j] & yp = scys[j] & zp = sczs[k]
;                    
;                    geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
;                        /refine, /ionosphere, $
;                        t89 = t89, t96 = t96, t01 = t01, ts04 = ts04, storm = storm
;                        
;                    ; convert from gsm to mag.
;                    geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
;                    
;                    scfmlats[i,j,k] = asin(tzf/r0)*deg
;                    scfmlons[i,j,k] = atan(tyf,txf)*deg
;                    
;                    scfmlts[i,j,k] = slon2lt(scfmlons[i,j,k], tet, /mag, /deg)/15
;                endfor
;            endfor
;        endfor
;        scfmlts = (scfmlts+24) mod 24
;        ; grid mlt and mlat in ionosphere; grid r0, azimuthal angle, z in tail.
;        info = {fmlt:scfmlts, fmlat:scfmlats, r0:scr0, ts:scts, zs:sczs}
;        
;        tvar = tmodel+'_grid_vertical'
;        store_data, tvar, uts, info
;    endforeach
;    
;    
;    
;    ; ** map the equatorial grid to ionosphere.
;    neqt = 60+1   ; # of grid in MLT.
;    neqr = 10+1   ; # of grid in R.
;
;    tilt0 = 13      ; deg, based both on model (12.99) and eyeballing (12).    
;    eqlrng = [4.5,6.5]      ; L-shell range.
;    eqtrng = 24-[2.8,-1.7]  ; MLT range.
;    eqrs = smkarthm(eqlrng[0],eqlrng[1],neqr, 'n')
;    eqts = smkarthm(eqtrng[0],eqtrng[1],neqt, 'n')*15*rad
;    
;    scfmlats = dblarr(nrec,neqt,neqr)
;    scfmlons = dblarr(nrec,neqt,neqr)
;    scfmlts = dblarr(nrec,neqt,neqr)
;    
;    
;    foreach tmodel, models do begin
;        get_data, tmodel+'_par', tmp, dat
;        pars = sinterpol(dat, tmp, uts)
;        ets = stoepoch(uts,'unix')
;        
;        t89 = (tmodel eq 't89')? 1: 0
;        t96 = (tmodel eq 't96')? 1: 0
;        t01 = (tmodel eq 't01')? 1: 0
;        t04s = (tmodel eq 't04s')? 1: 0
;        storm = (tmodel eq 't04s')? 1: 0
;        
;        ; map in situ grid to ionosphere.
;        for i = 0, nrec-1 do begin
;            tet = ets[i]
;            geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
;            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
;            
;            par = reform(pars[i,*])
;            for k = 0, neqr-1 do begin    
;                for j = 0, neqt-1 do begin
;                    xp =-eqrs[k]*cos(tilt0*rad)*cos(eqts[j])
;                    yp =-eqrs[k]*cos(tilt0*rad)*sin(eqts[j])
;                    zp = eqrs[k]*sin(tilt0*rad)
;                    
;                    geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
;                        /refine, /ionosphere, $
;                        t89 = t89, t96 = t96, t01 = t01, ts04 = ts04, storm = storm
;                        
;                    ; convert from gsm to mag.
;                    geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
;                    
;                    scfmlats[i,j,k] = asin(tzf/r0)*deg
;                    scfmlons[i,j,k] = atan(tyf,txf)*deg
;                    
;                    scfmlts[i,j,k] = slon2lt(scfmlons[i,j,k], tet, /mag, /deg)/15
;                endfor
;            endfor
;        endfor
;        scfmlts = (scfmlts+24) mod 24
;        ; grid mlt and mlat in ionosphere; grid x, y, z in tail.
;        info = {fmlt:scfmlts, fmlat:scfmlats, z0:tilt0, rs:eqrs, ts:eqts}
;        
;        tvar = tmodel+'_grid_equatorial'
;        store_data, tvar, uts, info
;    endforeach
;    
;    ; update data file.
;    tplot_save, allvars, filename = datfn
;endif

end


; settings.
utr = time_double(['2008-03-05/06:00','2008-03-05/06:15'])
sites = ['*']
dr0 = 3d
ext0 = 'png'
tstr = time_string(utr[0],tformat='YYYY_MMDD_hh')
rootdir = shomedir()+'/psbl_asi'


minlat = 55d
top = 254
deg = 180d/!dpi
rad = !dpi/180d

;store_data, '*', /delete
test_asi_donovan_load_data, utr
tplot_options, 'num_lab_min', 6


get_data, 'asf_mos', uts, mos, pxidx


; **** ewogram.
mincnt = 50d
maxcnt = 800d
mlat0 = 67.5d
dmlat = 0.1d
mltrg = [-1.5,1.5]
ewomlts = smkarthm(mltrg[0],mltrg[1],0.1,'dx')
newomlt = n_elements(ewomlts)

figfn = rootdir+'/ewo_'+tstr+'.pdf'
tvar = 'donovan_ewo'
;if tnames(tvar) eq '' then begin
if 1 then begin


    get_data, 'asf_mos', uts, mos, pxidx
    nrec = n_elements(uts)
    
    ewo = fltarr(nrec,newomlt)
    
    get_data, 'asf_info', tmp, info
    imgsz = info.imgsz
    
    txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
    tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
    txs = txs/imgsz[0]*2
    tys = tys/imgsz[0]*2
    mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
    mlts = atan(tys, txs)*deg+90    ; in deg.
    mlts = mlts/15          ; in hr.
    idx = where(mlts gt 12, cnt)    ; make mlt in [-12,12]
    if cnt ne 0 then mlts[idx] = mlts[idx]-12
    
    
    ewoidx = where(abs(mlats-mlat0) le dmlat, cnt)
    if cnt eq 0 then message, 'no pixel found ...'
    
    for i = 0, nrec-1 do begin
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        
        tewo = timg[ewoidx]
        tmlt = mlts[ewoidx]
        idx = sort(tmlt)
        tewo = tewo[idx]
        tmlt = tmlt[idx]
        ewo[i,*] = interpol(tewo,tmlt,ewomlts)
    endfor
    
    store_data, tvar, uts, ewo, ewomlts, limits = $
        {spec:1, no_interp:1, $
        ytitle:'EWO!CMLT (hr)', yrange:mltrg, ystyle:1, $
        ztitle:'Count (#)', zrange:[mincnt, maxcnt]}
        
    sgopen, figfn, xsize = 5, ysize = 2.5, /inch
    device, decomposed = 0
    loadct, 40
    tplot, tvar, trange = utr
    sgclose
endif

stop

; **** generate plots and make a movie.
mincnt = 50d
maxcnt = 400d
picdir = roodir+'/donovan'
movfn = rootdir+'/thg_'+tstr+'.mp4'
if file_test(movfn) eq 0 then begin
    get_data, 'asf_info', tmp, info
    mlats = info.mlats
    mlts = info.mlts
    imgsz = info.imgsz
    
    get_data, 'asf_mos', uts, mos, pxidx
    
    
    tpos = [0,0,1,1]
    xr = [-90,90]
    ;xr = [-90,0]
    yr = [minlat,90]
    yticks = [55,65,75]
    
    
    nrec = n_elements(uts)
    for i = 0, nrec-1 do begin
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        timg = bytscl(timg, min = mincnt, max = maxcnt, top = top)
        tut = uts[i]
        
        ofn = picdir+'/thg_rb_'+ $
            time_string(tut,tformat='YYYY_MMDD_hhmm_ss')+'.'+ext0
            
        sgopen, ofn, xsize = 10, ysize = 5, /inch
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        device, decomposed = 0
        loadct, 1
        sgtv, timg, position = tpos
        
        
        loadct2, 43
        sgset_map, position = tpos, color = 255, xrange = xr, $
            yrange = yr, ytickv = yticks
        xyouts, tpos[0]+xchsz, tpos[1]+ychsz, /normal, $
            time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color = 255
            
        sgclose
    endfor
    spic2movie, picdir, movfn, ext0
endif







end

