
function read_thm_asi_site, et0, site, $
    height = height, minlat = minlat, type = type, $
    midn = midn, minelev = minelev, dark = dark, weight = weight, $
    locroot = locroot, remroot = remroot, $
    png = png, pdf = pdf, save = save
    
    compile_opt idl2
    on_error, 0
    
    ; settings.
    top = 254
    dymsc = 86400000d
    hrmsc = 3600000d
    dr = 3000d          ; data rate, 3 sec.
    loadct2, 40
    
    ; prepare etr, ets, nrec.
    case n_elements(et0) of
        0: message, 'no epoch info ...'
        1: etr = [et0,et0]
        2: etr = et0
        else: ets = et0
    endcase
    
    if n_elements(etr) eq 0 then etr = [min(ets),max(ets)] else $
        ets = smkarthm(etr[0]-(etr[0] mod dr),etr[1]-(etr[1] mod dr),dr, 'dx')
        
    nrec = n_elements(ets)
    
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
        
    tmsc = (type eq 'asf')? hrmsc: dymsc
    npx = (type eq 'asf')? 256: 32          ; single site image size.
    isz = (type eq 'asf')? 3*npx: 8*npx     ; mosaic image size.
    
    ; final site.
    site = site[0]
    
    ; collect asc data: mlon/mlat/elev.
    case type of
        'asf': ascinfo = {site:'',elev:fltarr(npx,npx), $
            mlon:fltarr(npx+1,npx+1),mlat:fltarr(npx+1,npx+1)}
        'ast': ascinfo = {site:'',elev:fltarr(npx,npx), $
            mlon:fltarr(4,npx*npx),mlat:fltarr(4,npx*npx)}
    endcase
    
    ascfn = sptn2fn_thm_asi(et0, site, locroot, type = 'asc')
    ascfn = sprepfile(ascfn, locroot, remroot)
    vroot = strjoin(['thg',type,site],'_')
    rr = 1
    cdfid = cdf_open(ascfn)
    cdf_varget, cdfid, vroot+'_mlon', mlon, rec_start = rr, /zvariable
    cdf_varget, cdfid, vroot+'_mlat', mlat, rec_start = rr, /zvariable
    cdf_varget, cdfid, vroot+'_elev', elev, rec_start = rr, /zvariable
    if type eq 'ast' then begin
        mlon = 0.5*(mlon[0,*]+mlon[3,*])    ; [ll,ul,ur,lr].
        mlat = 0.5*(mlat[0,*]+mlat[1,*])
        mlons = mlon[sort(mlon)] & mlons = mlons[uniq(mlons)]
        mlats = mlat[sort(mlat)] & mlats = mlats[uniq(mlats)]
        nmlon = n_elements(mlons)
        nmlat = n_elements(mlats)
        mlons = mlons # (bytarr(nmlat)+1)
        mlats = (bytarr(nmlon)+1) # mlats
        astidx = intarr(npx*npx)
        for i = 0,npx*npx-1 do $            ; mapping index.
            astidx[i] = where(mlon[i] eq mlons and mlat[i] eq mlats)
    endif else begin
        mlon = reform(mlon[1,*,*],[npx+1,npx+1])
        mlat = reform(mlat[1,*,*],[npx+1,npx+1])
    endelse
    ascinfo.site = site
    ascinfo.mlon = mlon
    ascinfo.mlat = mlat
    ascinfo.elev = elev
    
    ; determine how many loops are needed.
    etrs = [etr[0]-(etr[0] mod tmsc),etr[1]-(etr[1] mod tmsc)]
    if (etr[1] mod tmsc) eq 0 then etrs[1] -= tmsc
    etrs = smkarthm(etrs[0],etrs[1],tmsc,'dx')
    netr = n_elements(etrs)
    
    foreach et0, etrs do begin
        print, 'processing for '+sfmepoch(et0)+' ...'
        ; for each loop, prepare asi files for site.
        asifn = sptn2fn_thm_asi(et0, site, locroot, type = type)
        asifn = sprepfile(asifn, locroot, remroot)
        if asifn eq '' then continue   ; no file for this loop.
        
        cdfid = cdf_open(asifn)
        
        idx = where(ets ge et0 and ets lt et0+tmsc, cnt)
        if cnt eq 0 then continue   ; no epoch for this loop.
        tets = ets[idx]
        
        ; prepare the epoch for site.
        vname = strjoin(['thg',type,site],'_')+'_time'
        cdf_control, cdfid, variable = vname, get_var_info = tmp
        cdf_varget, cdfid, vname, tmp, rec_count = tmp.maxrec+1
        utptr = ptr_new(reform(tmp))
        
        ; loop each epoch.
        foreach ttet, tets do begin
            ttut = 0.001d*ttet-62167219200d
            rec = where(*utptr eq ttut, cnt)
            if cnt eq 0 then continue   ; no asi for this site at ttet.
            print, sfmepoch(ttet)+' ...'
            
            mlon = ascinfo.mlon
            mlat = ascinfo.mlat
            elev = ascinfo.elev
            
            ; read raw image.
            vname = strjoin(['thg',type,site],'_')
            cdf_varget, cdfid, vname, tmp, rec_start = rec[0]
            img = double(tmp)
            
            ; remove edge.
            if type eq 'asf' then begin
                edge = where(elev lt minelev or ~finite(elev))
                img[edge] = mean(img[0:10,0:10], /nan)
                ; crude corner average subtraction. use thg_asf_site_offset?
                img = img - img[0] > 0
            endif else begin
                timg = dblarr(nmlon,nmlat)
                for i = 0, npx*npx-1 do timg[astidx[i]] = img[i]
                img = timg
            endelse
            
            ; scale luminosity to color index. adopted from thm_asi_merge_mosaic.
            idx = where(img ne 0)
            img[idx] *= 64D/(median(img[idx]) > 1)
            img <= top

            ; apply true weight.
            if keyword_set(weight) then multi = sin(elev*!dtor) $
            else multi = bytarr(npx,npx)+1b
            if type eq 'asf' then img *= multi
            
            ; output.
            fn = shomedir()+'/'+site+'/thg_'+type+'_'+site+'_'+$
                time_string(ttut,tformat='YYYY_MMDD_hhmm_ss')
            if keyword_set(png) then begin
                fn = fn+'.png'
                sgopen, fn, xsize = npx, ysize = npx
                sgtv, img, position = [0d,0,1,1], /nan
                sgclose
            endif else if keyword_set(pdf) then begin
                fn = fn+'.pdf'
                sgopen, fn, xsize = npx, ysize = npx
                sgtv, img, position = [0d,0,1,1], /nan
                sgclose
            endif else begin
                sgopen, 0, xsize = npx, ysize = npx
                sgtv, img, position = [0d,0,1,1], /nan
            endelse
        endforeach
        
        ; release epoch ptrarr.
        ptr_free, utptr
        ; close cdfid.
        cdf_close, cdfid
    endforeach
end

loadct, 2
etr = stoepoch(['2013-04-14/05:10','2013-04-14/05:11'])
type = 'asf'
sites = ['fsim']
fn = read_thm_asi_site(etr, sites, type=type, /pdf)

etr = stoepoch(['2013-06-01/02:00','2013-06-01/03:00'])
type = 'asf'
sites = ['chbg']
fn = read_thm_asi_site(etr, sites, type=type, /pdf)
end