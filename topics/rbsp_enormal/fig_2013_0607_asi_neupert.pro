
;_2013_0607_load_data


device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'xgridstyle', 1
tplot_options, 'zcharsize', 0.8



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1


get_data, 'asf_info', tmp, info
mlats = info.mlats
mlts = info.mlts
imgsz = info.imgsz
get_data, 'asf_mos', uts, mos, pxidx


utr = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
idx = where(uts ge utr[0] and uts le utr[1], nrec)
uts = uts[idx]
mos = mos[idx,*]

nrec = n_elements(uts)
phcnts = dblarr(nrec)


txs = findgen(imgsz[0]) # (fltarr(imgsz[1])+1) & txs = txs-imgsz[0]/2
tys = (fltarr(imgsz[0])+1) # findgen(imgsz[1]) & tys = tys-imgsz[0]/2
txs = txs/imgsz[0]*2
tys = tys/imgsz[0]*2
mlats = 90-sqrt(txs^2+tys^2)*(90-minlat)
mlts = atan(tys, txs)*deg+90    ; in deg.
mlts = (mlts+360) mod 360

mltrng = [21.6,22.6]     ; in hr.

mltsecs = smkarthm(mltrng[0],mltrng[1],0.1,'dx')
mltsecs*= 15        ; in deg.


tmltrng = mltrng*15

for i = 0, nrec-1 do begin
    timg = fltarr(imgsz)
    timg[pxidx] = mos[i,*]
    idx = where(mlts lt tmltrng[0] or mlts gt tmltrng[1])
    
    ;        tmp = timg
    ;        tmp[idx] = 0
    ;        tv, tmp
    ;        wait, 0.05
    
    timg[idx] = 0
    ; remove bg.
    idx = where(timg le 50)
    timg[idx] = 0
    phcnts[i] = total(timg)
endfor

titl = 'Integrated photon count at PINA between '+sgnum2str(tmltrng[0]/15)+' to '+sgnum2str(tmltrng[1]/15)+' MLT'

tvar = 'int_photon_cnt_3sites'
store_data, tvar, uts, phcnts
options, tvar, 'ytitle', 'Photon count'

ofn = shomedir()+'/maybe_neupert'+sgnum2str(tmltrng[0]/15)+'-'+sgnum2str(tmltrng[1]/15)+'.pdf'
;ofn = 0
sgopen, ofn, xsize = 6, ysize = 4, /inch

tpos = sgcalcpos(1)
tplot, tvar, /novtitle, position = tpos
xyouts, 0.5*(tpos[0]+tpos[2]), tpos[3]+0.01, /normal, alignment = 0.5, titl

sgclose
stop

for j = 0, n_elements(mltsecs)-2 do begin
    
    tmltrng = mltsecs[j:j+1]

    for i = 0, nrec-1 do begin
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        idx = where(mlts lt tmltrng[0] or mlts gt tmltrng[1])
        
;        tmp = timg
;        tmp[idx] = 0
;        tv, tmp
;        wait, 0.05
        
        timg[idx] = 0
        ; remove bg.
        idx = where(timg le 50)
        timg[idx] = 0
        phcnts[i] = total(timg)
    endfor
    
    titl = 'Integrated photon count at PINA between '+sgnum2str(tmltrng[0]/15)+' to '+sgnum2str(tmltrng[1]/15)+' MLT'
    
    tvar = 'int_photon_cnt_3sites'
    store_data, tvar, uts, phcnts
    options, tvar, 'ytitle', 'Photon count'
    
    ofn = shomedir()+'/maybe_neupert'+sgnum2str(tmltrng[0]/15)+'-'+sgnum2str(tmltrng[1]/15)+'.pdf'
    ;ofn = 0
    sgopen, ofn, xsize = 6, ysize = 4, /inch
    
    tpos = sgcalcpos(1)
    tplot, tvar, /novtitle, position = tpos
    xyouts, 0.5*(tpos[0]+tpos[2]), tpos[3]+0.01, /normal, alignment = 0.5, titl
    
    sgclose

end





end
