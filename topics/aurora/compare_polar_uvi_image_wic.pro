
pro compare_polar_uvi_image_wic, utr, rootdir = rootdir


; root directory, in YYYY_MMDD.
if n_elements(rootdir) eq 0 then rootdir = shomedir()+'/aurora_compare2'
figdir = rootdir+'/'+time_string(utr[0],tformat='YYYY_MMDD')
if file_test(figdir,/directory) eq 0 then file_mkdir, rootdir

logfn = rootdir+'/polar_uvi_image_wic.log'
if file_test(logfn) eq 0 then stouch, logfn
openw, loglun, logfn, /append, /get_lun
printf, loglun, ''
printf, loglun, time_string(utr)
free_lun, loglun

minlat = 50
height = 110
xr = [-90,90]       ; longitude range.  
yr = [minlat,90]    ; latitude range.
ytickpos = -45      ; latitude tick position.
white = sgcolor('white')
xtickname = ['18','00','06']
ct = 39             ; color table 39.
xsz = 400
ysz = 400
console = -1
top = 254

; polar uvi on the left.
uvipos = [0.0,0,1,0.5]
wicpos = [0.0,0.5,1,1]
uvitxtpos = [0.01,0.01]
wictxtpos = [0.01,0.51]

; read image data.
wic = sread_image_fuv(utr, minlat = minlat, /half, height = 130)
if size(wic,/type) ne 8 then begin
    openw, loglun, logfn, /append, /get_lun
    printf, loglun, 'Image/WIC no data ...'
    free_lun, loglun
    return
endif

uvi = sread_polar_uvi(utr, minlat = minlat, /half, height = 110)
if size(uvi,/type) ne 8 then begin
    openw, loglun, logfn, /append, /get_lun
    printf, loglun, 'Polar/UVI no data ...'
    free_lun, loglun
    return
endif

ets = wic.epoch
uts = sfmepoch(ets,'unix')
nrec = n_elements(ets)

for i = 1, nrec-1 do begin

    ; test aurora image.
    img1 = reform(wic.mltimg[i,*,*])
    idx = where(img1 eq 0, cnt)
    if cnt gt 0.8*n_elements(img1) then begin
        printf, console, 'image view is too small ...'
        continue   ; mostly 0.
    endif

    ; read polar data.
    tutr = uts[i-1:i]
    tetr = stoepoch(tutr,'unix')
    idx = where(uvi.epoch gt tetr[0] and uvi.epoch le tetr[1], cnt)
    if cnt eq 0 then continue
    img2 = total(uvi.mltimg[idx,*,*],1,/nan)/cnt
    idx = where(img2 eq 0, cnt)
    if cnt gt 0.8*n_elements(img2) then begin
        printf, console, 'polar view is too small ...'
        continue   ; mostly 0.
    endif
    
    tet = ets[i]
    fn = figdir+'/polar_image_aurora_'+sfmepoch(tet,'YYYY_MMDD_hhmm_ss')+'.png'
    sgopen, fn, xsize = xsz, ysize = ysz
    
    ; image on the bottom.
    timg = bytscl(img1, min = 100, max = 3000)<top
    sgtv, timg, position = wicpos, ct = ct
    sgset_map, xrange = xr, yrange = yr, pos = wicpos, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    xyouts, wictxtpos[0],wictxtpos[1],/normal, 'Image/WIC '+sfmepoch(tet), color = white
    
    ; polar on the top.
    timg = img2<top
    sgtv, timg, position = uvipos, ct = ct
    sgset_map, xrange = xr, yrange = yr, pos = uvipos, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    xyouts, uvitxtpos[0],uvitxtpos[1],/normal, 'Polar/UVI '+sfmepoch(tet), color = white
    
    sgclose
endfor

end

tr = ['2001-05-28','2002-01-01']
tr = ['2002-01-01','2003-01-01']
tr = ['2003-01-01','2004-01-01']
tr = ['2003-09-09','2004-01-01']
tr = ['2004-01-01','2005-01-01']
tr = ['2005-01-01','2005-07-09']

tr = ['2003-09-09','2005-07-09']


omni = sread_omni(tr)
ets = omni.epoch
aes = omni.ae_index

etr = stoepoch(tr)
idx = where(ets ge etr[0] and ets le etr[1])
ets = ets[idx]
aes = aes[idx]
omni = 0

ae0 = 1000      ; nT.
maxdt = 30*60   ; 30 min.
minnfig = 10    ; want more than 10 figures.
rootdir = shomedir()+'/aurora_compare2'

idx = where(aes ge ae0, cnt)
uts = sfmepoch(ets[idx],'unix')
idx = where(uts[1:*]-uts[0:cnt-2] ge maxdt)
utrs = [[uts[[0,idx+1]]],[uts[[idx,cnt-1]]]]

dt = maxdt/2    ; 15 min.
idx = where(utrs[*,1]-utrs[*,0] gt dt, cnt)
if cnt gt 0 then utrs = utrs[idx,*]
nutr = n_elements(utrs)/2

for i = 0, nutr-1 do begin
    ; clean up.
    if nutr eq 1 then checkcleanup = 1 else $
        if i eq 0 then checkcleanup = 0 else begin
        tmp = (utrs[i,0]-(utrs[i,0] mod 86400))- $
            (utrs[i-1,0]-(utrs[i-1,0] mod 86400))
        if tmp eq 0 then checkcleanup = 0 else checkcleanup = 1
    endelse
    if checkcleanup eq 0 then continue
    
    if i ne 0 then tutr = reform(utrs[i-1,*])   ; prevous time.
    figdir = rootdir+'/'+time_string(tutr[0],tformat='YYYY_MMDD')
    figfns = file_search(figdir+'/*', count = nfigfn)
    if nfigfn le minnfig then begin
        ; remove the folder.
        for j = 0, nfigfn-1 do file_delete, figfns[j], /allow_nonexistent
        file_delete, figdir, /allow_nonexistent
        ; remove image data and polar data.
        tfns = sread_image_fuv(tutr, /get_filename)
        for j = 0, n_elements(tfns)-1 do if tfns[j] ne '' then file_delete, tfns[j]
        tfns = sread_polar_uvi(tutr, /get_filename)
        for j = 0, n_elements(tfns)-1 do if tfns[j] ne '' then file_delete, tfns[j]
        ; write log.
        logfn = rootdir+'/polar_uvi_image_wic.log'
        openw, loglun, logfn, /append, /get_lun
        printf, loglun, 'Not enough plots ...'
        free_lun, loglun
    endif
    
    ; new plots.
    tutr = reform(utrs[i,*])
    tutr = tutr-(tutr mod dt) & tutr[1]+= dt
    print, time_string(tutr)
    compare_polar_uvi_image_wic, tutr, rootdir = rootdir
    
endfor

end