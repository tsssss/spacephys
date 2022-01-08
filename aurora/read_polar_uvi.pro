

function read_polar_uvi, et0, l1fn, atfn, orfn, pafn, $
    height = height, minlat = minlat
    
    compile_opt idl2
    on_error, 0
    
    ; emission height and min geomagnetic latitude.
    if n_elements(height) eq 0 then height = 120D   ; km in altitude.
    if n_elements(minlat) eq 0 then minlat = 50D    ; degree.
    
    ; prepare root dir.
    if n_elements(locroot) eq 0 then $
        locroot = spreproot('polar')
    if n_elements(remroot) eq 0 then $
        remroot = 'ftp://cdaweb.gsfc.nasa.gov/pub/data/polar'
        
    ; find local file.
    if n_params() eq 2 then begin
        l1fn = et0; & rec = rec
    endif else if n_params() ne 5 then begin
        ptn = locroot+'/uvi/uvi_level1/YYYY/po_level1_uvi_YYYYMMDD_v[0-9]{2}.cdf'
        l1fn = sprepfile(sptn2fn(ptn,et0), locroot, remroot)
;        ptn = locroot+'/uvi_k0/YYYY/po_k0_uvi_YYYYMMDD_v[0-9]{2}.cdf'
        ptn = locroot+'/orbit/pre_at/YYYY/po_at_pre_YYYYMMDD_v[0-9]{2}.cdf'
        atfn = sprepfile(sptn2fn(ptn,et0), locroot, remroot)
        ptn = locroot+'/orbit/pre_or/YYYY/po_or_pre_YYYYMMDD_v[0-9]{2}.cdf'
        orfn = sprepfile(sptn2fn(ptn,et0), locroot, remroot)
        ptn = locroot+'/orbit/def_pa/YYYY/po_pa_def_YYYYMMDD_v[0-9]{2}.cdf'
        pafn = sprepfile(sptn2fn(ptn,et0), locroot, remroot)
    endif

    
    ; get cdf id.
    if ~file_test(l1fn) then $
        message, 'file '+l1fn+' does not exist ...'
    cdfid = cdf_open(l1fn)
    
    ; record.
    if n_elements(rec) eq 0 then begin
        vname = 'EPOCH'
        cdf_control, cdfid, variable = vname, get_var_info = tmp, /zvariable
        cdf_varget, cdfid, vname, et, rec_count = tmp.maxrec+1
        et = reform(et[0,0,*])
        tmp = min(et-et0, rec, /absolute)
        dr = sdatarate(et)      ; datarate in msec.
        if abs(et[rec]-et0) ge dr then begin
            message, 'no polar/uvi at the time ...', /continue
            cdf_close, cdfid
            return, -1
        endif
    endif
    
    ; epoch.
    cdf_varget1, cdfid, 'EPOCH', et, rec_start = rec
    print, sfmepoch(et, 'yyyy-MM-dd hh:mm:ss.fff')
    
    ; image.
    cdf_varget, cdfid, 'INT_IMAGE', img, rec_start = rec
    
    ; variables for housekeeping.
    ; door, filter, gain, integration period, system, uvi mode.
    cdf_varget1, cdfid, 'DOORPOS', door, rec_start = rec    ; 0:unk 1:fullex 2:open 3:closed
    cdf_varget1, cdfid, 'FILTER', filter, rec_start = rec   ; 0:home 1:1304 2:1356 3:lbhs 4:lbhl 5:solr 6:shutter 7:unk
    cdf_varget1, cdfid, 'GAIN', gain, rec_start = rec       ; ?.
    cdf_varget1, cdfid, 'FRAMERATE', rate, rec_start = rec  ; 1,2,4.
    cdf_varget1, cdfid, 'SYSTEM', system, rec_start = rec & system += 1   ; 0:primary, 1:secondary.
    cdf_varget1, cdfid, 'UVIMODE', mode, rec_start = rec    ; mode.
   
    cdf_close, cdfid
    
    fet = et - floor((rate+4)*9.2)*1000d        ; frame time.
    
    ; despun angle at frame time in deg.
    cdfid = cdf_open(pafn)
    cdfinfo = cdf_inquire(cdfid)
    cdf_varget, cdfid, 'Epoch', et1, rec_count = cdfinfo.maxrec+1
    et1 = reform(et1[0,*])
    dr = 60000d     ; 1 min max.
    tmp = min(et1-et, rec, /absolute)
    if abs(et1[rec]-et) ge dr then begin
        message, 'no platform attitude info at the time ...', /continue
        cdf_close, cdfid
        return, -1
    endif
    cdf_varget, cdfid, 'DSP_ANGLE', dsp, rec_count = cdfinfo.maxrec+1
    dsp = reform(dsp[0,*])/!dtor
    ; remove invalid angle > 180 deg.
    idx = where(abs(dsp) le 180, cnt)
    if cnt eq 0 then message, 'no valid despun angle ...'
    et1 = et1[idx]
    dsp = dsp[idx]
    ; locate the wanted section.
    idx = where(et1 ge fet and et1 le et, cnt)
    if cnt le 3 then idx = [rec-1,rec,rec+1]
    if idx[0] lt 0 then idx-= idx[0]
    if idx[-1] ge cdfinfo.maxrec then idx-= (idx[-1]-cdfinfo.maxrec)
    dsp_angle = spl_interp(et1[idx],dsp[idx],spl_init(et1[idx],dsp[idx]),fet)
    cdf_close, cdfid
    ; use psl file to correct? see get_dsp_angles.pro.

    ; read GCI pos at frame time in km.
    cdfid = cdf_open(orfn)
    cdfinfo = cdf_inquire(cdfid)
    cdf_varget, cdfid, 'Epoch', et1, rec_count = cdfinfo.maxrec+1
    et1 = reform(et1[0,*])
    dr = 60000d     ; 1 min max.
    tmp = min(et1-fet, rec, /absolute)
    if abs(et1[rec]-fet) ge dr then begin
        message, 'no orbit info at the time ...', /continue
        cdf_close, cdfid
        return, -1
    endif
    cdf_varget, cdfid, 'GCI_POS', orbit, rec_start = rec 
    cdf_close, cdfid

    ; read attitude at frame time. no interpolation?
    cdfid = cdf_open(atfn)
    cdfinfo = cdf_inquire(cdfid)
    cdf_varget, cdfid, 'Epoch', et1, rec_count = cdfinfo.maxrec+1
    et1 = reform(et1[0,*])
    dr = 600000d     ; 1 min max.
    tmp = min(et1-fet, rec, /absolute)
    if abs(et1[rec]-fet) ge dr then begin
        message, 'no attitude info at the time ...', /continue
        cdf_close, cdfid
        return, -1
    endif
    cdf_varget1, cdfid, 'GCI_R_ASCENSION', ra, rec_start = rec
    cdf_varget1, cdfid, 'GCI_DECLINATION', dec, rec_start = rec
    att = [cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec)]
    cdf_close, cdfid
    
    ; calculate glat/glon, from looking direction at frame time.
    polar_uvilook, orbit, att, dsp_angle, filter, lookdir, system = system
    polar_ptg, fet, height, att, orbit, system, lookdir, glat, glon, /geodetic 
    
    ; do line-of-sight and dayglow correction.
;    polar_uvi_corr, fet, orbit, system, glat, glon, img
        
    ; other info.
    sphere = orbit[2] gt 0
    if n_elements(imgsz) eq 0 then imgsz = 4*(90-minlat) else imgsz -= 1
        
    ; get mlat/mlon. method 1: geotoapex.
    rootdir = srootdir()+'/image'
    apexfile = rootdir+'/support/mlatlon.1997a.xdr'
    geotoapex, glat, glon, apexfile, mlat, mlon
    get_local_time, fet, glat, glon, apexfile, glt, mlt
    get_mlt_image, img, mlat, mlt, minlat, sphere, mltimg, ncell = imgsz
    
;    print, minmax(mltimg)
    if max(mltimg,/nan) eq 0 then return, img
    
    rt = shomedir()+'/uvi'
    if file_test(rt) eq 0 then file_mkdir, rt
    ofn = rt+'/po_uvi_'+sfmepoch(fet,'YYYY_MMDD_HHmm_ss')+'.pdf'
;    sgopen, 0, xsize = imgsz, ysize = imgsz
    sgopen, ofn, xsize = imgsz, ysize = imgsz
    device, decomposed = 0
    loadct, 39, /silent
    tvscl, mltimg, /nan
    !p.charsize = 0.5
    sgset_map, yrange = [minlat,90], position = [0d,0,1,1], $
        ytickpos = 135, color = 255, xtickpos = [54,52,49,52]
    !p.charsize = 0
    sgclose
    return, img
    
end

l1fn = '/Volumes/Research/data/polar/uvi/uvi_level1/1997/po_level1_uvi_19970501_v01.cdf'
atfn = '/Volumes/Research/data/polar/orbit/def_at/1997/po_at_def_19970501_v01.cdf'
orfn = '/Volumes/Research/data/polar/orbit/def_or/1997/po_or_def_19970501_v01.cdf'
pafn = '/Volumes/Research/data/polar/orbit/def_pa/1997/po_pa_def_19970501_v01.cdf'
et = stoepoch('1997-05-01/20:22:20')    ; wygant.
uvi = read_polar_uvi(et, l1fn, atfn, orfn, pafn, minlat=50)
end
;;et = stoepoch('1997-05-01/20:25:20')    ; wygant.
;;et = stoepoch('1997-05-01/20:32:00')    ; wygant.
;;et = stoepoch('1997-05-09/05:44:50')    ; wygant.
;;et = stoepoch('1997-05-09/05:42:50')    ; wygant.
;;et = stoepoch('1997-05-01/10:25:29')
;;et = stoepoch('1997-05-01/17:18:52')
;et = stoepoch('1997-05-01/19:31:60')
;;et = stoepoch('2008-01-09/00:01:44')
;et = stoepoch('1997-05-01/20:22:30')    ; wygant.
;uvi = read_polar_uvi(et)
;
;locroot = spreproot('polar')
;remroot = 'ftp://cdaweb.gsfc.nasa.gov/pub/data/polar'
;et0 = stoepoch('2007-03-24')
;ptn = locroot+'/uvi/uvi_level1/YYYY/po_level1_uvi_YYYYMMDD_v[0-9]{2}.cdf'


;l1fn = '/Volumes/Research/data/polar/uvi/uvi_level1/2006/po_level1_uvi_20061124_v01.cdf'
;atfn = '/Volumes/Research/data/polar/orbit/pre_at/2006/po_at_pre_20061124_v02.cdf'
;orfn = '/Volumes/Research/data/polar/orbit/pre_or/2006/po_or_pre_20061124_v03.cdf'
;pafn = '/Volumes/Research/data/polar/orbit/def_pa/2006/po_pa_def_20061124_v01.cdf'
;
;cdfid = cdf_open(l1fn)
;vname = 'EPOCH'
;cdf_control, cdfid, variable = vname, get_var_info = tmp
;cdf_varget, cdfid, vname, ets, rec_count = tmp.maxrec+1
;cdf_close, cdfid
;ets = reform(ets[0,0,*])
;
;for i = 0, n_elements(ets)-1 do $
;    uvi = read_polar_uvi(ets[i], minlat = 50, l1fn, atfn, orfn, pafn)
;end
