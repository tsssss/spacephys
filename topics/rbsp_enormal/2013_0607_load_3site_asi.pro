


; settings.
utr = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
utr0 = utr+[-1,1]*300   ; 5 min pad time.

; utr = time_double(['2013-06-07/04:50','2013-06-07/04:50:03'])
reload = 0


; 't01s' gives same result as 't04s'.
; the B0 model ~ measured. tried t04s, t96, t89, B0 are different.
models = ['t89','t04s','t01','t96']
nmodel = n_elements(models)
probes = ['a','b']
nprobe = n_elements(probes)


dr0 = 1d/16
site = 'pina'            ; the site under conjunction.
sites = ['pina','kapu','chbg']         ; sites for asf.
scsyms = [1,6]           ; symbols for rbsp-a and -b.
scclrs = [6,255]         ; colors for rbsp-a and -b.


datfn = shomedir()+'/psbl_de_32hz/2013_0607_data.tplot'
mapfn = shomedir()+'/psbl_de_32hz/2013_0607_map_data.tplot'
asifn = shomedir()+'/psbl_de_32hz/thg_asf_mosaic_2013_0607_0445_3sites.cdf'


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



tpos = [0,0,1,1]
xr = [-90,90]
yr = [minlat,90]
yticks = [55,65,75]

for i = 0, nrec-1 do begin
    tut = uts[i]
    timg = img
    timg[pxidx] = mos[i,*]
    
    ofn = shomedir()+'/thg_asi3/thg_rb_'+ $
        time_string(tut,tformat='YYYY_MMDD_hhmm_ss')+'.png'
;    ofn = 0
    sgopen, ofn, xsize = 12, ysize = 6, /inch
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    device, decomposed = 0
    loadct2, 40
    sgtv, timg, position = tpos


    loadct2, 43
    sgset_map, position = tpos, color = 255, xrange = xr, $
        yrange = yr, ytickv = yticks
    xyouts, tpos[0]+xchsz, tpos[1]+ychsz, /normal, $
        time_string(tut,tformat='YYYY-MM-DD/hh:mm:ss'), color = 255, charsize = chsz
        
    sgclose
    
endfor


end
