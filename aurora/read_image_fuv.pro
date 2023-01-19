;+
; Type: function.
; Purpose: Read aurora image on IMAGE/FUV, and calc the mapped image.
; Parameters:
;   et0, in, string/double, req. String for filename, double for epoch.
;   rec, in, int, opt. Supply record index w/ filename. Omit if epoch is given.
; Keywords:
;   type, in, string, opt. 'sie', 'sip', 'wic'. 'wic' by default.
;   spin, in, float, opt. Camera spin phase in degree. Adjust when cdf contains
;       wrong spin phase info.
;   height, in, float, opt. Assumed emission height, default is 130 km.
;   minlat, in, float, opt. Minimum mlat for mapped image, default is 50 deg.
;   ncell, in, int, opt. Size of the mapped image, default is 161 pixel.
;   locroot, in, string, opt. Local root directory.
;   remroot, in, string, opt. Remote root directory.
; Return: structure. Contain raw image, mapped image, glat/glon at each pixel.
; Notes: Set horizontal field of view equal to vertical, because sometimes hfov
;   is wrong. Not sure should multiply or divide fov scale.
;       azim, elev, roll are contained in the cdf file but sometimes are wrong,
;   better to get them using
; Dependence: slib.
; History:
;   2009-04, Sheng Tian, create.
;   2015-03-04, Sheng Tian, re-write.
;-

function read_image_fuv, et0, rec, spin = spinphase, $
    height = height, minlat = minlat, type = type, $
    locroot = locroot, remroot = remroot, $
    ncell = imgsz

    compile_opt idl2
    on_error, 0

    ; emission height and min geomagnetic latitude.
    if n_elements(height) eq 0 then height = 130D   ; km in altitude.
    if n_elements(minlat) eq 0 then minlat = 50D    ; degree.

    ; choose ele/ion aurora or wic.
    if n_elements(type) eq 0 then type = 'wic'
    
    ; prepare root dir.
    if n_elements(locroot) eq 0 then $
        locroot = spreproot('image/fuv')
    if n_elements(remroot) eq 0 then $
        remroot = 'ftp://cdaweb.gsfc.nasa.gov/pub/data/image/fuv'
    
    ; find local file.
    if n_params() eq 2 then begin
        fuvfn = et0; & rec = rec
    endif else begin
        ptn = locroot+'/'+type+'_k0/YYYY/im_k0_'+type+'_YYYYMMDD_v01.cdf'
        fuvfn = sptn2fn(ptn, et0)
;        fuvfn = sprepfile(fuvfn, locroot, remroot)
    endelse
    
    ; get cdf id.
    if ~file_test(fuvfn) then $
        message, 'file '+fuvfn+' does not exist ...'
    cdfid = cdf_open(fuvfn)
    cdfinq = cdf_inquire(cdfid)
    
    ; record.
    if n_elements(rec) eq 0 then begin
        vname = 'EPOCH'
        cdf_control, cdfid, variable = vname, get_var_info = tmp, /zvariable
        cdf_varget, cdfid, vname, et, rec_count = tmp.maxrec+1, /zvariable
        et = reform(et[0,0,*])
        tmp = min(et-et0, rec, /absolute)
        dr = sdatarate(et)      ; datarate in msec.
        if abs(et[rec]-et0) ge dr then begin
            message, 'no image/fuv at the time ...', /continue
            return, -1
        endif
    endif
    
    ; epoch.
    cdf_varget1, cdfid, 'EPOCH', et, rec_start = rec
    print, sfmepoch(et, 'yyyy-MM-dd hh:mm:ss.fff')
    
    ; instrument id.
    cdf_varget1, cdfid, 'INSTRUMENT_ID', instr, /string, rec_start = rec
    print, instr
    instr = strtrim(instr,2)
    vname = instr
    if vname ne 'WIC' then vname = 'SI'
    
    ; image.
    cdf_varget, cdfid, vname+'_PIXELS', img, rec_start = rec
    npx = (size(img, /dimensions))[0]
    case vname of
        'WIC': rawimg = float(reverse(rotate(img,1),2))
        'SI': rawimg = float(rotate(img,3))
    endcase        
        
    ; horizontal and vertical field of view.
    cdf_varget1, cdfid, 'VFOV', vfov, rec_start = rec
    hfov = vfov
    
    ; field of view scale.
    cdf_varget1, cdfid, 'FOVSCALE', fovscl, rec_start = rec
    
    ; horizontal and vertical angular resolution, degree per pixel.
    hfov = hfov*fovscl/npx
    vfov = vfov*fovscl/npx
    
    ; spin axis in GCI.
    cdf_varget1, cdfid, 'SV_X', x, rec_start = rec
    cdf_varget1, cdfid, 'SV_Y', y, rec_start = rec
    cdf_varget1, cdfid, 'SV_Z', z, rec_start = rec
    sv = [x,y,z]
    
    ; spin phase.
    if n_elements(spinphase) eq 0 then $
        cdf_varget1, cdfid, 'SPINPHASE', spinphase, rec_start = rec
    if abs(spinphase) gt 180 then begin
        message, 'invalid spin phase ...', /continue
        return, -1
    endif
    print, 'spin phase: ', spinphase

    ; spacecraft position in GCI.
    cdf_varget1, cdfid, 'ORB_X', x, rec_start = rec
    cdf_varget1, cdfid, 'ORB_Y', y, rec_start = rec
    cdf_varget1, cdfid, 'ORB_Z', z, rec_start = rec
    orbit = [x,y,z]
    
    ; spin attitude in spacecraft system.
    cdf_varget1, cdfid, 'SCSV_X', x, rec_start = rec
    cdf_varget1, cdfid, 'SCSV_Y', y, rec_start = rec
    cdf_varget1, cdfid, 'SCSV_Z', z, rec_start = rec
    scsv = [x,y,z]
    cdf_close, cdfid

    ; instrument azimuth, co elevation, and roll angle.
    cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
    rootdir = srootdir()+'/image'
    tmp = where(instr eq ['WIC','SI1356','SI1218'])
    image_get_inst_angles_p, tmp, yr, stodoy(yr,mo,dy), tmp, rootdir
    azim = tmp[0] & elev = tmp[1] & roll = tmp[2]

    ; get pointing info, glat/glon.
    time = lonarr(2)
    time[0] = yr*1000L+stodoy(yr,mo,dy)
    time[1] = 1000L*(sc+60*(mi+60*hr))+msc     ; mili second of day.
    image_ptg, height, sv, orbit, scsv, $
        spinphase, time, npx, npx, vfov, hfov, azim, elev, roll, $
        glat, glon  ; output
    glat[where(abs(glat) gt 1e20)] = !values.d_nan
    glon[where(abs(glon) gt 1e20)] = !values.d_nan
    
    ; other info.
    sphere = orbit[2] gt 0
    if n_elements(imgsz) eq 0 then imgsz = 4*(90-minlat) else imgsz -= 1
        
    ; get mlat/mlon. method 1: geo2apex.
    apexfile = rootdir+'/support/mlatlon.1997a.xdr'
    geo2apex, glat, glon, apexfile, mlat, mlon
    get_local_time, et, glat, glon, glt, mlt
    get_mlt_image, rawimg, mlat, mlt, minlat, sphere, mltimg, ncell = imgsz
    
    return, {rawimg: rawimg, mltimg: mltimg, glat:glat, glon:glon}

end

;fpn = spreproot('image')+'/fuv/wic_k0/2000/im_k0_wic_20001203_v01.cdf'
;wic = read_image_fuv(fpn, 320)
wic = read_image_fuv(stoepoch('2001-10-22/08:19:39'))
tvscl, wic.mltimg
end