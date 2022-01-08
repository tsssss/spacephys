


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.
;utr = time_double(['2013-06-07/05:00','2013-06-07/05:30'])
reload = 0

; adjust the range to control color.
maxcnt = 500
mincnt = 50


datfn = shomedir()+'/psbl_de_32hz/2013_0607_data.tplot'
mapfn = shomedir()+'/psbl_de_32hz/2013_0607_map_data.tplot'
asifn = shomedir()+'/psbl_de_32hz/thg_asf_mosaic_2013_0607_0445_3sites.cdf'



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
    sites = ['pina','kapu','chbg']
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





; **** plot aurora images.
scsyms = [1,6]           ; symbols for rbsp-a and -b.
scclrs = [6,255]         ; colors for rbsp-a and -b.

tpos = [0,0,1,1]
xr = [-90,90]
;xr = [-90,0]
yr = [minlat,90]
yticks = [55,65,75]

for i = 0, nrec-1 do begin
    tut = uts[i]
    timg = img
    timg[pxidx] = mos[i,*]
;    timg = timg[0:imgsz[0]/2,*]
    
    timg = bytscl(timg, max=maxcnt, min=mincnt, /nan, top=254)
    
    ofn = shomedir()+'/thg_asi_3sites/thg_rb_'+ $
        time_string(tut,tformat='YYYY_MMDD_hhmm_ss')+'.png'
    ;    ofn = 0
    sgopen, ofn, xsize = 6, ysize = 6, /inch
    sgopen, ofn, xsize = 12, ysize = 6, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    device, decomposed = 0
    loadct, 1
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