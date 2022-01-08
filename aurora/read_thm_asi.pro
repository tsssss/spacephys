;+
; et0, in/out
; minlat, strongly suggest to ignore or set to 50 (deg).
; height, strongly suggest to ignore or set to 110 (km).
; minelev, no need to set in most cases.
; weight, set to do illumination correction, on trial.
; dark, set to 30-50 to exclude trees, etc. default is 50.
;-
function read_thm_asi, et0, site0, exclude = exclude, $
    height = height, minlat = minlat, type = type, $
    midn = midn, minelev = minelev, dark = dark, weight = weight, $
    locroot = locroot, remroot = remroot, $
    plot = plot, half = half, save = save

    compile_opt idl2
    on_error, 0

    top = 254

    ; emission height and min geomagnetic latitude.
    if n_elements(height) eq 0 then height = 110D   ; km in altitude.
    if n_elements(minlat) eq 0 then minlat = 60D    ; degree.

    ; choose thumbnail ast, or full res asf.
    if n_elements(type) eq 0 then type = 'asf'
    
    ; min elevation.
    if n_elements(minelev) eq 0 then begin
        if type eq 'ast' then minelev = -10
        if type eq 'asf' then minelev = 5
    endif
    
    ; dark threshold.
    if n_elements(dark) eq 0 then dark = keyword_set(weight)? 5d: 50d

    if n_elements(locroot) eq 0 then $
        locroot = spreproot('themis')
    if n_elements(remroot) eq 0 then $
        remroot = 'http://themis.ssl.berkeley.edu/data/themis'

    ; **** prepare file names.
    nsite = n_elements(site0)
    if nsite eq 0 then begin
        site0 = ['atha','chbg','ekat','fsmi','fsim','fykn',$
            'gako','gbay','gill','inuv','kapu','kian',$
            'kuuj','mcgr','pgeo','pina','rank','snkq',$
            'tpas','whit','yknf','nrsq','snap','talo']
        nsite = n_elements(site0)
    endif
    sites = site0
    
    ; exclude sites.
    for i = 0, n_elements(exclude)-1 do begin
        idx = where(sites eq exclude[i], tmp)
        if tmp gt 0 then sites[idx] = ''
    endfor
    idx = where(sites ne '', tmp)
    if tmp eq 0 then begin
        message, 'no site ...', /continue
        return, -1
    endif

    ; find local files, if not found then download remote files.
    sites = sites[idx]
    nsite = n_elements(sites)
    asifns = sptn2fn_thm_asi(et0, sites, locroot, type = type)
    ascfns = sptn2fn_thm_asi(et0, sites, locroot, type = 'asc')
;    for i = 0, nsite-1 do begin
;        asifns[i] = sprepfile(asifns[i], locroot, remroot)
;        ascfns[i] = sprepfile(ascfns[i], locroot, remroot)
;    endfor
    idx = where(asifns ne '' and ascfns ne '', cnt)
    if cnt eq 0 then begin
        message, 'no remote data available ...', /continue
        return, -1
    endif
    asifns = asifns[idx]
    ascfns = ascfns[idx]
    sites = sites[idx]
    nsite = n_elements(sites)
    
    ; find midnight mlon, in degree, use all available calib files.
    tascfns = file_search(count = cnt, $
        sptn2fn_thm_asi(et0, '????', locroot, type = 'asc'))
    cmlon = fltarr(cnt)
    cmidn = fltarr(cnt)
    for i = 0, cnt-1 do begin
        cdfid = cdf_open(tascfns[i])
        cdf_varget, cdfid, 25, tmp, /string, /zvariable ; 'thg_asc_site_midn'.
        if tmp eq '     ' then begin
            cmidn[i] = !values.d_nan
        endif else $
            cmidn[i] = float(strmid(tmp,0,2))+float(strmid(tmp,3,2))/60
        cdf_varget, cdfid, 22, tmp, /zvariable  ; 'thg_asc_site_mlon', 23:mlat.
        cmlon[i] = tmp
        cdf_close, cdfid
    endfor
    idx = sort(cmidn)
    midn = interpol(cmlon[idx], cmidn[idx], (et0/86400000D mod 1)*24, /nan)

    ; final image.
    npx = (type eq 'asf')? 256: 32
    isz = (type eq 'asf')? 3*npx: 8*npx
    del = 0.5*isz
    imf = dblarr(isz, isz)
    imc = dblarr(isz, isz)

    ; loop each site.
    for ii = 0, nsite-1 do begin
        ; read data.
        tfn = asifns[ii]
        if file_test(tfn) eq 0 then continue
        print, 'reading '+tfn+' ...'
        vroot = strjoin(['thg',type,sites[ii]],'_')
        cdfid = cdf_open(tfn)

        ; read epoch.
        vname = vroot+'_time'
        cdf_control, cdfid, variable = vname, get_var_info = tmp, /zvariable
        cdf_varget, cdfid, vname, ut, rec_count = tmp.maxrec+1, /zvariable
        ut = reform(ut, /overwrite)
        dr = sdatarate(ut)

        ; find record.
        ut0 = 0.001D*et0-62167219200D   ; convert et to ut.
        tmp = min(ut-ut0, rec, /absolute)
        if abs(ut[rec]-ut0) ge dr then begin
            message, 'no asi at site: '+sites[ii]+' ...', /continue
            continue
        endif
        
        ; raw image and epoch.
        et = 1000D*ut[rec]+62167219200000D
        cdf_varget, cdfid, vroot, img, rec_start = rec, /zvariable
        npx = (size(img, /dimensions))[0]   ; image size.
        et0 = et    ; overwrite original et.

        cdf_close, cdfid

        ; calibration data.
        tfn = ascfns[ii]
        cdfid = cdf_open(tfn)

        ; read elevation, mlat, mlon, glat, glon in degree.
        rr = 1
        cdf_varget, cdfid, vroot+'_elev', elev, rec_start = rr, /zvariable
        cdf_varget, cdfid, vroot+'_mlat', mlat, rec_start = rr, /zvariable
        cdf_varget, cdfid, vroot+'_mlon', mlon, rec_start = rr, /zvariable
        cdf_varget, cdfid, vroot+'_glat', glat, rec_start = rr, /zvariable
        cdf_varget, cdfid, vroot+'_glon', glon, rec_start = rr, /zvariable
        if type eq 'ast' then begin
            elev = reform(elev, [npx,npx])
        endif else begin
            mlon = reform(mlon[1,*,*],[npx+1,npx+1])
            mlat = reform(mlat[1,*,*],[npx+1,npx+1])
            glon = reform(glon[1,*,*],[npx+1,npx+1])
            glat = reform(glat[1,*,*],[npx+1,npx+1])
        endelse

        ; luminosity multiplier, site midnight.
        if type eq 'ast' then $
            multi = make_array(size(elev,/dimensions), value = 1d) $
        else cdf_varget, cdfid, vroot+'_multiply', multi, /zvariable
        cdf_close, cdfid
        if keyword_set(weight) then multi = sin(elev*!dtor)

        ; image processing.
        img = double(img)
        
        ; remove edge.
        if type eq 'asf' then begin
            edge = where(elev lt minelev or ~finite(elev))
            img[edge] = mean(img[0:10,0:10], /nan)
            ; crude corner average subtraction. use thg_asf_site_offset?
            img = img - img[0] > 0
        endif else begin
            ;img = rotate(img, 2)
        endelse
        ; scale luminosity to color index. adopted from thm_asi_merge_mosaic.
        img *= 64D/(median(img) > 1)
        img <= top
        ; apply true weight.
        img *= multi

        lat = mlat & lon = mlon
        ; convert lat/lon corner to x/y corner grid.
        r = (90-lat)/(90-minlat)
        t = (lon-(90+midn))*(!dpi/180)  ; need rot mid night to -y.
        xc = r*cos(t) & yc = r*sin(t)   ; xc,yc in [-1,1].
        ; bin x/y corner grid to uniform x/y center grid.
        xc = floor((xc+1)*del+0.5)
        yc = floor((yc+1)*del+0.5)
        ; get illuminated pixel index.
        idx = where(finite(elev) and elev gt minelev)
        for k = 0, n_elements(idx)-1 do begin
            tk = idx[k]
            if img[tk] le dark then continue
            ; extract pixel corner.
            if type eq 'ast' then begin
                xcor = xc[*,tk] & ycor = yc[*,tk]
            endif else begin
                tidx = array_indices([npx,npx],tk,/dimensions)
                xcor = xc[tidx[0]:tidx[0]+1,tidx[1]:tidx[1]+1]
                ycor = yc[tidx[0]:tidx[0]+1,tidx[1]:tidx[1]+1]
            endelse
            ib = max(xcor, min = ia)
            jb = max(ycor, min = ja)
            if max([ia,ib,ja,jb]) gt isz then continue  ; nan.
            ia = ia > 0
            ja = ja > 0
            ib = ib < (isz-1) > ia
            jb = jb < (isz-1) > ja
            imc[ia:ib,ja:jb] += multi[tk]
            imf[ia:ib,ja:jb] += img[tk]
        endfor
    endfor

    if n_elements(et) eq 0 then return, -1      ; no data.
    et0 = et
    imf /= imc
    xsz = isz
    ysz = xsz
    if keyword_set(half) then ysz = ceil(ysz*0.5)
    imf = imf[*,0:ysz-1]
    if ~keyword_set(plot) then return, imf

    ; plot and return wid.
    rt = shomedir()+'/'+type
    if file_test(rt) eq 0 then file_mkdir, rt
    fn = rt+'/thm_'+type+'_'+sfmepoch(et, 'YYYY_MMDD_hhmmss')+'.png'
    
    white = 255
    xsz = npx > 600
    ysz = xsz
    lim = [minlat-0.01,-180,90,180]
    if keyword_set(half) then begin
        ysz *= 0.5 & lim[1] = -90 & lim[3] = 90
    endif
    
    sgopen, fn, xsize=xsz, ysize=ysz

    device, decomposed = 0
    loadct, 1

    imf = congrid(imf, xsz, ysz)
    tv, imf
    map_set, name = 'AzimuthalEquidistant', 90, 0, 0, $
        /noborder, position = [0,0,1,1], limit = lim, $
        latdel = 10, londel = 45, glinestyle = 1, $
        /noerase, color = white, /isotropic
    xyouts, [0,88,180,272], minlat+[2,2,3,2], ['12','18','00','06'], $
        alignment = 0.5, /data, color = white
    xyouts, 45+intarr(5), [5,6,7,8,9]*10-2, ['50','60','70','80','90'], $
        /data, color = white
    xyouts, 10, 10, sfmepoch(et), /device, color = white

    sgclose
end

;asi = read_thm_asi(stoepoch('2008-01-09/21:10'), /plot, /half, type='ast')
;asi = read_thm_asi(stoepoch('2013-05-01/07:38'), /plot, /half, type='ast', /weight, exclude = ['fsmi','fsim','snkq'])

uts = smkarthm(time_double('2014-08-28/04:50'),time_double('2014-08-28/05:30'),3,'dx')
ets = stoepoch(uts, 'unix')
foreach et, ets do asi = read_thm_asi(et, /plot, /half, type='ast', /weight, /save)
end