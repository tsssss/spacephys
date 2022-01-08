

;t0 = '2007-03-24/07:40:56'
;minlat = 50
;xsz = 500 & ysz = xsz
;window, 0, xsize = 2*xsz, ysize = ysz
;white = 255
;uvipos = [0.05,0.15,0.45,0.55]
;asipos = [0.55,0.15,0.95,0.55]
;xr = [-90.1,90.1]
;yr = [minlat,90]
;ytickpos = -45
;xtickname = ['18','00','06']
;
;uvi = sread_polar_uvi(t0, minlat = minlat, /half)
;sz = size(uvi.mltimg,/dimensions)
;img = congrid(uvi.mltimg, xsz, xsz*sz[1]/sz[0], /interp)
;sz = size(img,/dimensions)
;sgtv, img, position = uvipos
;sgset_map, xrange = xr, yrange = yr, pos = uvipos, color = white, $
;    ytickpos = ytickpos, xtickname = xtickname
;xyouts, 0.05,0.05,/normal, 'Polar/UVI '+sfmepoch(uvi.epoch), color = white
;
;;sites = ['gill','tpas','atha','gako','whit','inuv']
;asi = sread_thm_asi(t0, type = 'asf', minlat = minlat, sites, /half)
;sz = size(uvi.mltimg,/dimensions)
;img = congrid(asi.mltimg, xsz, xsz*sz[1]/sz[0], /interp)
;sz = size(img,/dimensions)
;sgtv, img, position = asipos
;sgset_map, xrange = xr, yrange = yr, pos = asipos, color = white, $
;    ytickpos = ytickpos, xtickname = xtickname
;xyouts, 0.55,0.05,/normal, 'THEMIS/ASI '+sfmepoch(asi.epoch), color = white

;tr = ['2006-04-14/00:00','2006-04-15/00:00']
;tr = ['2006-08-20/00:00','2006-08-21/00:00']
;tr = ['2006-04-05/14:00','2006-04-05/17:00']
;tr = ['2006-04-09/07:00','2006-04-09/11:00']
;tr = ['2006-11-30/00:00','2006-12-01/00:00']
;tr = ['2008-03-09/00:00','2008-03-10/00:00']
;tr = ['2008-03-28/00:00','2008-03-29/00:00']
tr = ['2006-10-13/23:00','2006-10-13/24:00']
tr = ['2006-10-21/15:00','2006-10-21/16:00']
tr = ['2006-10-29/07:00','2006-10-29/08:00']
tr = ['2006-11-10/13:00','2006-11-10/14:00']
tr = ['2006-11-24/08:00','2006-11-24/09:00']
tr = ['2007-03-13/09:00','2007-03-13/10:00']
tr = ['2007-03-24/09:00','2007-03-24/10:00']
tr = ['2007-04-01/12:00','2007-04-01/13:00']
tr = ['2007-04-12/11:00','2007-04-12/12:00']
tr = ['2007-05-23/08:00','2007-05-23/09:00']
tr = ['2007-09-29/08:00','2007-09-29/09:00']
tr = ['2007-09-02/03:00','2007-09-02/04:00']
tr = ['2007-09-02/06:00','2007-09-02/07:00']

;tr = ['2007-09-23/08:00','2007-09-23/09:00']
;sites = ['fsmi','tpas','pina','inuv','gill']
;exclude = ['atha','pgeo','fsim','whit','kian','mcgr','fykn']

;tr = ['2007-09-23/06:00','2007-09-23/07:00']
;sites = ['fsmi','tpas','pina','inuv','gill']

;tr = ['2007-10-29/09:00','2007-10-29/09:03']

tr = ['2007-11-14/12:13','2007-11-14/13:00']

tr = ['2006-11-24/08:00','2006-11-24/09:00']
;sites = ['kian','fykn','whit','fsim','fsmi']
;exclude = ['mcgr','ekat','atha','gill','tpas','atha']

tr = ['2006-10-29/07:00','2006-10-29/08:00']
;sites = ['gill','pina','whit','fsim','fsmi']

tr = ['2007-11-05/10:00','2007-11-05/11:00']
;sites = ['fsim','whit','fykn']

tr = ['2007-11-20/13:00','2007-11-20/14:00']
sites = ['*']

tr = ['2007-11-21/09:00','2007-11-21/10:00']
sites = ['*']

tr = ['2008-02-10/07:00','2008-02-10/08:00']
tr = ['2008-02-10/08:00','2008-02-10/09:00']
tr = ['2008-02-10/09:00','2008-02-10/10:00']
tr = ['2008-02-10/11:00','2008-02-10/12:00']
tr = ['2008-02-10/12:00','2008-02-10/13:00']
tr = ['2008-02-10/13:00','2008-02-10/14:00']
sites = ['*']

tr = ['2008-02-04/09:00','2008-02-04/10:00']
tr = ['2008-02-04/10:00','2008-02-04/11:00']
tr = ['2008-02-04/11:00','2008-02-04/12:00']
sites = ['*']

tr = ['2007-11-14/12:00','2007-11-14/13:00']
sites = ['mcgr','whit']

tr = ['2008-02-02/03:00','2008-02-02/04:00']
tr = ['2008-02-02/04:00','2008-02-02/05:00']
tr = ['2008-02-02/05:00','2008-02-02/06:00']


tr = ['2008-02-02/13:00','2008-02-02/14:00']

tr = ['2008-02-03/08:00','2008-02-03/11:00']
tr = ['2008-02-03/11:00','2008-02-03/13:00']

tr = ['2008-01-05/09:00','2008-01-05/10:00']

tr = ['2008-01-06/13:00','2008-01-06/15:00']

tr = ['2008-01-07/11:00','2008-01-07/12:00']

tr = ['2008-01-08/09:45','2008-01-08/10:00']

tr = ['2008-01-14/11:00','2008-01-14/12:00']

tr = ['2008-01-13/13:00','2008-01-13/14:00']

tr = ['2008-01-15/06:00','2008-01-15/07:00']

tr = ['2008-02-01/22:00','2008-02-02/03:00']

tr = ['2008-02-10/15:00','2008-02-10/18:00']

tr = ['2007-12-18/11:00','2007-12-18/12:00']
tr = ['2007-12-18/12:00','2007-12-18/18:00']

tr = ['2007-12-11/08:00','2007-12-11/10:00']
tr = ['2007-12-11/10:00','2007-12-11/13:00']

sites = ['gbay','kuuj']
;excludes = ['whit','fsim','fsmi']
tr = ['2008-01-05/21:00','2008-01-05/23:00']

; settings.
utr = time_double(tr)

; root directory, in YYYY_MMDD.
rootdir = shomedir()+'/aurora_compare'+'/'+time_string(utr[0],tformat='YYYY_MMDD')
if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir

minlat = 50
height = 110
xr = [-90,90]       ; longitude range.  
yr = [minlat,90]    ; latitude range.
ytickpos = -45      ; latgitude tick position.
white = sgcolor('white')
xtickname = ['18','00','06']
ct = 39             ; color table 39.
xsz = 400
ysz = 400
console = -1

; polar uvi on t he left.
uvipos = [0.0,0,1,0.5]
asfpos = [0.0,0.5,1,1]
uvitxtpos = [0.01,0.01]
asftxtpos = [0.01,0.51]

; read polar data.
uvi = sread_polar_uvi(tr, minlat = minlat, /half, height = height)

ets = uvi.epoch
uts = sfmepoch(ets,'unix')
nrec = n_elements(ets)

; plot polar uvi, read and plot themis asf at the nearest time.
for i = 1, nrec-1 do begin
;for i = 10,11 do begin

    ; test polar image.
    timg = uvi.mltimg[i,*,*]
    idx = where(timg eq 0, cnt)
    if cnt gt 0.8*n_elements(timg) then begin
        printf, console, 'polar view is too small ...'
        continue   ; mostly 0.
    endif
    
    ; read themis asf.
    printf, console, 'reading themis/asf '+time_string(uts[i])+' ...'
    asi = sread_thg_mlt(uts[i-1:i], sites, minlat = minlat, /half, /nomoon, type = 'asf', height = height)
    if size(asi,/type) ne 8 then continue
    
    tmp = asi.mltimg
    tmp = where(finite(tmp,/nan), cnt)
    if cnt eq n_elements(asi.mltimg) then continue  ; all nan.
    
    tet = ets[i]
    fn = rootdir+'/polar_themis_aurora_'+sfmepoch(tet,'YYYY_MMDD_hhmm_ss')+'.png'
    sgopen, fn, xsize = xsz, ysize = ysz
;    sgindexcolor, ct = ct, /silent
    
    ; themis on the bottom.
    timg = total(asi.mltimg,1,/nan)/n_elements(asi.epoch)
    sgtv, timg, position = asfpos, ct = ct
    sgset_map, xrange = xr, yrange = yr, pos = asfpos, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    xyouts, asftxtpos[0],asftxtpos[1],/normal, 'Themis/ASI '+sfmepoch(tet), color = white
    
    ; polar on the top.
    timg = reform(uvi.mltimg[i,*,*])
    sgtv, timg, position = uvipos, ct = ct
    sgset_map, xrange = xr, yrange = yr, pos = uvipos, $
        ytickpos = ytickpos, xtickname = xtickname, color = white
    xyouts, uvitxtpos[0],uvitxtpos[1],/normal, 'Polar/UVI '+sfmepoch(tet), color = white
    
    sgclose
endfor

end