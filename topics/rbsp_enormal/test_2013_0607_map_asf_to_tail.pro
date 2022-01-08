

; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1
dr0 = 3

top = 254
maxcnt = 800
maxcnt = 500
mincnt = 50
asifn = shomedir()+'/psbl_de_32hz/thg_asf_mosaic_2013_0607_0445_3sites.cdf'


tutr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
model = 't01'
dir = -1
sgeopack_par, tutr, model



; **** load asf data.
_2013_0607_load_data

load = 0
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


idx = where(uts ge tutr[0] and uts le tutr[1], idx)
uts = uts[idx]
mos = mos[idx,*]
ets = stoepoch(uts,'unix')



; x-y mesh in the ionosphere plane at 100 km.
; and corresponding mlat and mlt.
ionoxs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & ionoxs = ionoxs-imgsz[0]/2
ionoys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & ionoys = ionoys-imgsz[0]/2
ionoxs = ionoxs/imgsz[0]*2
ionoys = ionoys/imgsz[0]*2
mlats = 90-sqrt(ionoxs^2+ionoys^2)*(90-minlat)
mlts = atan(ionoys, ionoxs)*deg+90    ; in deg.
; mlts = (mlts/15+24) mod 24






; **** make a meshed grid at the equatorial plane tilted at 13 deg.
nrec = n_elements(uts)
utidx = [0,nrec/2,nrec-1]
utidx = where(uts eq time_double('2013-06-07/04:55:06'))
nrec = n_elements(utidx)

neqx = 30/3             ; along L-shell.
neqy = 30*2             ; along MLT.
eqlrng = [3.5,5.5]      ; L-shell range.
eqtrng = 24-[2.8,-1.7]  ; MLT range.
tilt0 = 13  ; deg, tilt.
eqls = smkarthm(eqlrng[0],eqlrng[1],neqx, 'n')
eqts = smkarthm(eqtrng[0],eqtrng[1],neqy, 'n')*15*rad


scfmlats = dblarr(nrec,neqx,neqy)
scfmlons = dblarr(nrec,neqx,neqy)
scfmlts = dblarr(nrec,neqx,neqy)

get_data, model+'_par', tmp, dat
pars = sinterpol(dat, tmp, uts[utidx])
ets = stoepoch(uts[utidx],'unix')

for i = 0, n_elements(utidx)-1 do begin
    i0 = utidx[i]
    tet = ets[i]
    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
    
    par = reform(pars[i,*])
    for j = 0, neqx-1 do begin
        for k = 0, neqy-1 do begin
            xp =-eqls[j]*cos(tilt0*rad)*cos(eqts[k])
            yp =-eqls[j]*cos(tilt0*rad)*sin(eqts[k])
            zp = eqls[j]*sin(tilt0*rad)
            
            geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
                /refine, /ionosphere, /t01
                
            ; convert from gsm to mag.
            geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
            
            scfmlats[i,j,k] = asin(tzf/r0)*deg
            scfmlons[i,j,k] = atan(tyf,txf)*deg
            
            scfmlts[i,j,k] = slon2lt(scfmlons[i,j,k], tet, /mag, /deg)/15
        endfor
    endfor
endfor
scfmlts = (scfmlts+24) mod 24


tpos = [0.1,0.1,0.9,0.9]
ofn = shomedir()+'/test_map_asf_to_tail_equatorial_grid.pdf'
;ofn = 0
sgopen, ofn, xsize = 8, ysize = 4, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

device, decomposed = 0
loadct2, 43

i0 = utidx[0]
timg = fltarr(imgsz)
timg[pxidx] = mos[i0,*]
timg = bytscl(timg, max=maxcnt, min=mincnt, top=top, /nan)
sgtv, timg, position = tpos, ct = 1

sgset_map, position = tpos, color = 255, xrange = [-90,90], $
    yrange = [minlat,90], ytickv = [55,65,75]
xyouts, tpos[0]+xchsz, tpos[1]+ychsz*0.2, /normal, $
    time_string(uts[i0],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
    
plot, [-1,1], [-1,0], xstyle=5, ystyle=5, /iso, position = tpos, /nodata, /noerase

scfrs = (90-scfmlats)/(90-minlat)
scfts = (24-scfmlts)*15*rad
scfxs =-scfrs*sin(scfts)
scfys =-scfrs*cos(scfts)

plots, scfxs, scfys, psym = 3, color = 6

xyouts, 0.5, tpos[3]+ychsz*0.2, /normal, alignment = 0.5, $
    'Model: '+strupcase(model)+'. Red dot: equatorial grid at tilt = '+$
    sgnum2str(tilt0)+' deg, R in ('+$
    sgnum2str(eqlrng[0])+','+sgnum2str(eqlrng[1])+') Re , MLT in ('+$
    sgnum2str(eqtrng[0] mod 24)+','+sgnum2str(eqtrng[1] mod 24)+') Re'
    
sgclose

stop






; **** make a meshed grid at the common L = 4.96.
nsct = 30*2
nscz = 30/3
nrec = n_elements(uts)
utidx = [0,nrec/2,nrec-1]
utidx = where(uts eq time_double('2013-06-07/04:55:06'))
nrec = n_elements(uts)

scr0 = 4.96
sctrng = [2.5,-1.5]
sczrng = [1,4]
scts = smkarthm(sctrng[0],sctrng[1], nsct, 'n')*15*rad
sczs = smkarthm(sczrng[0],sczrng[1], nscz, 'n')
scxs = -scr0*cos(scts)
scys =  scr0*sin(scts)

scfmlats = dblarr(nrec,nsct,nscz)
scfmlons = dblarr(nrec,nsct,nscz)
scfmlts = dblarr(nrec,nsct,nscz)

get_data, model+'_par', tmp, dat
pars = sinterpol(dat, tmp, uts[utidx])
ets = stoepoch(uts[utidx],'unix')

for i = 0, n_elements(utidx)-1 do begin
    i0 = utidx[i]
    tet = ets[i]
    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
    
    par = reform(pars[i,*])
    for j = 0, nsct-1 do begin
        for k = 0, nscz-1 do begin
            xp = scxs[j] & yp = scys[j] & zp = sczs[k]

            geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
                /refine, /ionosphere, /t01

            ; convert from gsm to mag.
            geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag

            scfmlats[i,j,k] = asin(tzf/r0)*deg
            scfmlons[i,j,k] = atan(tyf,txf)*deg

            scfmlts[i,j,k] = slon2lt(scfmlons[i,j,k], tet, /mag, /deg)/15
        endfor
    endfor
endfor
scfmlts = (scfmlts+24) mod 24


tpos = [0.1,0.1,0.9,0.9]
ofn = shomedir()+'/test_map_asf_to_tail_vertical_grid.pdf'
sgopen, ofn, xsize = 8, ysize = 4, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

device, decomposed = 0
loadct2, 43

i0 = utidx[0]
timg = fltarr(imgsz)
timg[pxidx] = mos[i0,*]
timg = bytscl(timg, max=maxcnt, min=mincnt, top=top, /nan)
sgtv, timg, position = tpos, ct = 1

sgset_map, position = tpos, color = 255, xrange = [-90,90], $
    yrange = [minlat,90], ytickv = [55,65,75]
xyouts, tpos[0]+xchsz, tpos[1]+ychsz*0.2, /normal, $
    time_string(uts[i0],tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz

plot, [-1,1], [-1,0], xstyle=5, ystyle=5, /iso, position = tpos, /nodata, /noerase

scfrs = (90-scfmlats)/(90-minlat)
scfts = (24-scfmlts)*15*rad
scfxs =-scfrs*sin(scfts)
scfys =-scfrs*cos(scfts)

plots, scfxs, scfys, psym = 3, color = 6

xyouts, 0.5, tpos[3]+ychsz*0.2, /normal, alignment = 0.5, $
    'Model: '+strupcase(model)+'. Red dot: vertical grid at R = '+$
    sgnum2str(scr0)+' Re, MLT in ('+$
    sgnum2str(-sctrng[0])+','+sgnum2str(-sctrng[1])+') hr, Z GSM in ('+$
    sgnum2str(sczrng[0])+','+sgnum2str(sczrng[1])+') Re'

sgclose

stop






; calc where the pixels within the grid at ionosphere, save to file.
scewos = fltarr(nrec,nsct)
zgsm0 = (2.84+2.80)*0.5 ; A: 2.84, B:2.80.
tmp = sczs-zgsm0
scewoidx = where(tmp eq min(tmp,/absolute))
noplot = 1

for i = 0, nrec-1 do begin
    timg = fltarr(imgsz)
    timg[pxidx] = mos[i,*]
    
    ; plot the grid on the aurora field of view.
;    window, 1, xsize = imgsz[0], ysize = imgsz[1]
;    tv, timg
;    plot, [0,imgsz[0]],[0,imgsz[1]], /nodata, /noerase, position = [0,0,1,1]
;    plots, fptxs[i,*,*], fptys[i,*,*], psym = 1

    
    scasi = fltarr(nsct,nscz)
    for j = 0, nsct-1 do begin
        for k = 0, nscz-1 do begin
            tj = fptxs[i,j,k]
            tk = fptys[i,j,k]
            scasi[j,k] = timg[tj,tk]
        endfor
    endfor
    
    scewos[i,*] = scasi[*,scewoidx]
    
    if noplot then continue
    
    ofn = shomedir()+'/scasi_'+model+'/scasi_'+model+'_'+ $
        time_string(uts[i],tformat='YYYY_MMDD_hhmm_ss')+'.png'
;    ofn = 0
    sgopen, ofn, xsize = 600, ysize = 450
    erase, sgcolor('white')
    device, decomposed = 0
    loadct2, 43
    
    tpos = [0.15,0.15,0.9,0.9]
    scimg = congrid(scasi, 600, 450)
    sgtv, scimg, position = tpos, ct = 1
    xr = [max(scts),min(scts)]*scr0
    yr = minmax(sczs)
    plot, xr, yr, xrange = xr, yrange = yr, xstyle = 1, ystyle = 1, /iso, /nodata, /noerase, $
        position = tpos, xtitle = 'Azimuthal distance from midnight (Re)', ytitle = 'z!DGSM!N (Re)'
    
    
    loadct2, 43
    foreach tprobe, probes do begin
        tc = (tprobe eq 'a')? 6: 4
        get_data, 'rbsp'+tprobe+'_pos_gsm', tuts, dat
        tposgsm = sinterpol(dat, tuts, uts[i])
        
        tz = tposgsm[2]
        tt = scr0*atan(tposgsm[1],-tposgsm[0])
        plots, tt, tz, psym = 1, color = tc
    endforeach
    
    sgclose
endfor

store_data, 'sc_ewogram', uts, scewos, scts*scr0, limits = $
    {ytitle:'SC EWOgram!C@L = 4.96, z = 2.82 Re!CAzm.Dist from midnight!C(Re)', $
    spec:1, no_interp:1, yrange:[1,4]}

stop


sz = size(fptmlats,/dimensions)
nrec = sz[0]
nmlt = sz[1]
nmlat = sz[2]


;for i = 0, 0 do begin
for i = 200, nrec-1 do begin
    timg = fltarr(imgsz)
    timg[pxidx] = mos[i,*]
    
    scimgsz = [nmlt,nmlat]-1
    scimg = fltarr(scimgsz)
    for j = 0, scimgsz[0]-1 do begin
        for k = 0, scimgsz[1]-1 do begin
            txrg = minmax(fptmlts[i,j:j+1,k])
            tyrg = minmax(fptmlats[i,j,k:k+1])
            idx = where(mlts ge txrg[0] and mlts le txrg[1] and $
                mlats ge tyrg[0] and mlats le tyrg[1], cnt)
            if cnt eq 0 then begin
                txrg = mean(txrg)
                tyrg = mean(tyrg)
                idx = where(abs(mlts-txrg) le 0.1/15 and $
                    abs(mlats-tyrg) le 0.1, cnt)
            endif
            if cnt eq 0 then stop
            scimg[j,k] = mean(timg[idx])
        endfor
    endfor
    
    window, 1, xsize = 800, ysize = 400
    device, decomposed = 0
    loadct, 1
    tv, congrid(scimg, 400, 400), 0
    tv, congrid(timg[0:imgsz[0]/2,*], 400, 400), 1
    xyouts, 0.01,0.9, /normal, time_string(uts[i]), color = 255
    stop
endfor



stop





; fpt photon count.
foreach tprobe, probes do begin
    pre0 = 'rbsp'+tprobe+'_'
    get_data, pre0+'fpt_mlt', tuts, scmlts
    get_data, pre0+'fpt_mlat', tuts, scmlats
    scmlts = interpol(scmlts[*,0], tuts, uts)*15-360
    scmlats = interpol(scmlats[*,0], tuts, uts)
    
    phcnts = fltarr(nrec)
    for i = 0, nrec-1 do begin
        tmlat = scmlats[i]
        tmlt = scmlts[i]
        idx = where(abs(mlats-tmlat) le 0.1 and abs(mlts-tmlt) le 0.1, cnt)
        if cnt eq 0 then continue
        
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        phcnts[i] = max(timg[idx])
    endfor
    store_data, pre0+'asf_cnt', uts, phcnts, limits = $
        {ytitle:'Photon cnt!C@Footprint', yrange:[0,300], yticks:2, ystyle:1, yminor:5}
    get_data, pre0+'pf_fac_mat_map', tuts, dat
    store_data, pre0+'pfb_fac_mat_map', tuts, dat[*,0], limits = $
        {ytitle:'S@100 km!C0.1-812 sec!C(mV/m!U2!N)', yrange:[-15,15], ystyle:1, yticks:2, yminor:5, labels:'S!D||!N'}
endforeach


vars = ['asf_cnt','pfb_fac_mat_map']
vars = ['rbspa_'+vars, 'rbspb_'+vars]
nvar = n_elements(vars)

ofn = shomedir()+'/fig_asf_pf.pdf'
sgopen, ofn, xsize = 8.5, ysize = 5, /inch
device, decomposed = 0

poss = sgcalcpos(nvar)
loadct2, 43
tplot, vars, vlab_margin = 12, position = poss
sgclose

stop



; ewogram.
tmlat = 65.1
yr = [21,23]*15
mlt0s = smkarthm(yr[0],yr[1],0.1,'dx')-360
nmlt0 = n_elements(mlt0s)


ewos = fltarr(nrec,nmlt0)
for i = 0, nrec-1 do begin
    idx = where(abs(mlats-tmlat) le 0.1, cnt)
    if cnt lt 3 then continue
    timg = fltarr(imgsz)
    timg[pxidx] = mos[i,*]
    
    timg = timg[idx]
    tmlts = mlts[idx]
    idx = sort(tmlts)
    timg = timg[idx]
    tmlts = tmlts[idx]
    
    ewos[i,*] = interpol(timg,tmlts,mlt0s)
endfor

tvar = 'ewogram_'+sgnum2str(tmlat)
store_data, tvar, uts, ewos, mlt0s/15+24, $
    limits = {spec:1, no_interpo:1, $
    ytitle:'Ewogram!C@'+sgnum2str(tmlat)+' deg!CMLT (hr)', yrange:yr/15, $
    yticks:2, yminor:5, $
    ztitle:'Norm.count', zrange:[0,256], zticks:2}





ofn = 0
ofn = shomedir()+'/fig_ewogram.pdf'
sgopen, ofn, xsize = 8.5, ysize = 5, /inch
device, decomposed = 0
loadct2, 43

tutr = time_double(['2013-06-07/04:54','2013-06-07/04:58'])

poss = sgcalcpos(2, ypad = 5)

tvar = 'ewogram_'+sgnum2str(tmlat)
tplot, tvar, trange = utr, position = poss[*,0], /noerase, vlab_margin = 12
timebar, tutr, color = 6

tplot, tvar, trange = tutr, position = poss[*,1], /noerase, vlab_margin = 12

txs = time_double(['2013-06-07/04:54:45','2013-06-07/04:55:53'])
tys = [22.59,21.65]
plot, tutr, yr/15, xstyle = 5, ystyle = 5, /nodata, /noerase, position = poss[*,1]
plots, txs, tys, color = 255, thick = 4

sgclose















if keyword_set(test_single_insitu_point) then begin
    ; **** track a point on the vertical plane at L = 4.96.
    scr0 = 4.96
    scmlt0 = 1.9*15*rad     ; in rad, hr before mid-night.
    scx0 =-scr0*cos(scmlt0)
    scy0 = scr0*sin(scmlt0)
    scz0 = 2.8
    
    
    fptmlats = dblarr(nrec)
    fptmlons = dblarr(nrec)
    fptmlts = dblarr(nrec)
    
    sgeopack_par, utr, model, /delete
    if model ne '' then begin
        get_data, model+'_par', tmp, dat
        pars = sinterpol(dat, tmp, uts)
    endif
    
    for i = 0, nrec-1 do begin
        tet = ets[i]
        geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
        
        par = pars[i,*]
        xp = scx0 & yp = scy0 & zp = scz0
        
        geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
            /refine, /ionosphere, /t01
            
        ; convert from gsm to mag.
        geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
        
        fptmlats[i] = asin(tzf/r0)*deg
        fptmlons[i] = atan(tyf,txf)*deg
        
        fptmlts[i] = slon2lt(fptmlons[i], tet, /mag, /deg)/15
    endfor
    fptmlts = (fptmlts+24) mod 24
    
    idx = [0, nrec/2, nrec-1]
    
    ofn = shomedir()+'/test_map_asf_to_tail_use_interpol.pdf'
    ;ofn = 0
    sgopen, ofn, xsize = 4, ysize = 3, /inch
    
    device, decomposed = 0
    loadct2, 43
    tplot_options, 'num_lab_min', 4
    
    poss = sgcalcpos(2,lmargin=10,rmargin=2,tmargin=3,bmargin=3)
    
    tdat = fptmlats
    yr = minmax(tdat) & yr = (yr[1]+yr[0])*0.5+0.5*[-1,1]*0.14
    store_data, 'tmp', uts, tdat, limits = {ytitle:'Fpt Mlat!C(deg)', yrange:yr, ystyle:1}
    tplot, 'tmp', position = poss[*,0], /noerase, /nouttick
    plot, [0,nrec-1], yr, position = poss[*,0], /noerase, /nodata, $
        yrange=yr, ystyle=5, xstyle=5
    plots, idx, tdat[idx], psym = 1, color = 6
    tdat = interpol(tdat[idx],idx,findgen(nrec),/quadratic)
    plots, findgen(nrec), tdat, color = 6
    
    
    tdat = fptmlts
    yr = minmax(tdat) & yr = (yr[1]+yr[0])*0.5+0.5*[-1,1]*0.004
    store_data, 'tmp', uts, tdat, limits = {ytitle:'Fpt MLT!C(hr)', yrange:yr, ystyle:1}
    tplot, 'tmp', position = poss[*,1], /noerase, /novtitle
    plot, [0,nrec-1], yr, position = poss[*,1], /noerase, /nodata, $
        yrange=yr, ystyle=5, xstyle=5
    plots, idx, tdat[idx], psym = 1, color = 6
    tdat = interpol(tdat[idx],idx,findgen(nrec),/quadratic)
    plots, findgen(nrec), tdat, color = 6
    
    
    xyouts, poss[0,0], 0.92, /normal, alignment = 0, 'Test point at: ('+ $
        sgnum2str(scx0,ndec=1)+','+sgnum2str(scy0,ndec=1)+','+sgnum2str(scz0,ndec=1)+') Re GSM!CModel:'+model+', red: quadratic interpol'
        
    sgclose
endif
stop

end
