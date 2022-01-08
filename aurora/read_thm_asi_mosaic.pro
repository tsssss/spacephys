;+
; This version is for batch mosaic generation.
; 
; et0, in, epochs for mosaic, or epoch range.
; minlat, strongly suggest to ignore or set to 50 (deg).
; height, strongly suggest to ignore or set to 110 (km).
; minelev, no need to set in most cases.
; weight, set to do illumination correction, on trial.
; dark, set to 30-50 to exclude trees, etc. default is 50.
;-

function read_thm_asi_mosaic, et0, site0, exclude = exclude, $
    height = height, minlat = minlat, type = type, $
    midn = midn, minelev = minelev, dark = dark, weight = weight, $
    locroot = locroot, remroot = remroot, $
    plot = plot, half = half, save = save

    compile_opt idl2
    on_error, 0
    
    ; settings.
    top = 254
    dymsc = 86400000d
    hrmsc = 3600000d
    dr = 3000d          ; data rate, 3 sec.
    
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

    ; final sites.
    sites = sites[idx]
    nsite = n_elements(sites)
    
    ; calibration info for 24 curent sites: name, center mlon/mlat, midnigh hour.
    ascinfo = [$
        {site:'atha',cmlon:-55.0432d,cmlat:62.3675d,midn:8.18333d},$
        {site:'chbg',cmlon: 2.93617d,cmlat:60.5017d,midn:4.51667d},$
        {site:'ekat',cmlon:-57.1432d,cmlat:72.5488d,midn:8.31667d},$
        {site:'fsim',cmlon:-68.4154d,cmlat:67.6388d,midn:9.08333d},$
        {site:'fsmi',cmlon:-55.5899d,cmlat:67.8262d,midn:8.21667d},$
        {site:'fykn',cmlon:-96.1186d,cmlat:67.3232d,midn:10.9167d},$
        {site:'gako',cmlon:-92.9370d,cmlat:63.1690d,midn:10.7000d},$
        {site:'gbay',cmlon: 23.2090d,cmlat:61.7209d,midn:3.30000d},$
        {site:'gill',cmlon:-28.9852d,cmlat:66.9531d,midn:6.48333d},$
        {site:'inuv',cmlon:-87.3175d,cmlat:71.3305d,midn:10.3500d},$
        {site:'kapu',cmlon:-9.12731d,cmlat:60.6328d,midn:5.25000d},$
        {site:'kian',cmlon:-108.494d,cmlat:65.0508d,midn:11.6833d},$
        {site:'kuuj',cmlon: 13.1838d,cmlat:67.8662d,midn:3.90000d},$
        {site:'mcgr',cmlon:-101.031d,cmlat:61.9210d,midn:11.2167d},$
        {site:'nrsq',cmlon: 43.9138d,cmlat:66.4655d,midn:!values.d_nan},$
        {site:'pgeo',cmlon:-65.9696d,cmlat:59.4184d,midn:!values.d_nan},$
        {site:'pina',cmlon:-30.1183d,cmlat:60.7356d,midn:6.55000d},$
        {site:'rank',cmlon:-26.4093d,cmlat:73.2397d,midn:6.31667d},$
        {site:'snap',cmlon:-55.8303d,cmlat:70.4358d,midn:!values.d_nan},$
        {site:'snkq',cmlon:-3.88124d,cmlat:67.4142d,midn:4.93333d},$
        {site:'talo',cmlon:-32.7407d,cmlat:79.2442d,midn:6.71667d},$
        {site:'tpas',cmlon:-38.0022d,cmlat:63.8674d,midn:7.05000d},$
        {site:'whit',cmlon:-82.8524d,cmlat:63.8677d,midn:10.0500d},$
        {site:'yknf',cmlon:-60.4068d,cmlat:69.7826d,midn:8.53333d}]

;    ; this block prints the above structure.
;    tascfns = file_search(count = cnt, $
;        sptn2fn_thm_asi(et0, '????', locroot, type = 'asc'))
;    for i = 0, cnt-1 do begin
;        cdfid = cdf_open(tascfns[i])
;        cdf_varget, cdfid, 25, tmp, /string, /zvariable
;        if tmp eq '     ' then begin
;            cmidn = '!values.d_nan'
;        endif else $
;            cmidn = strtrim(string($
;                float(strmid(tmp,0,2))+float(strmid(tmp,3,2))/60),2)+'d'
;        cdf_varget, cdfid, 22, cmlon, /zvariable
;        cdf_varget, cdfid, 23, cmlat, /zvariable
;        cdf_close, cdfid
;        print, "{site:'"+strmid(file_basename(tascfns[i]),11,4)+"',mlon:"+$
;            strtrim(string(cmlon),2)+'d,mlat:'+$
;            strtrim(string(cmlat),2)+'d,midn:'+cmidn+'},$'
;    endfor

    ; collect asc data: cmlon, cmidn, mlon/mlat/elev.
    cmlons = ascinfo.cmlon
    cmidns = ascinfo.midn
    idx = sort(cmidns)
    cmlons = cmlons[idx]
    cmidns = cmidns[idx]
    case type of
        'asf': ascinfo = replicate({site:'',elev:fltarr(npx,npx), $
            mlon:fltarr(npx+1,npx+1),mlat:fltarr(npx+1,npx+1)},nsite)
        'ast': ascinfo = replicate({site:'',elev:fltarr(npx,npx), $
            mlon:fltarr(4,npx*npx),mlat:fltarr(4,npx*npx)},nsite)
    endcase
    
    for i = 0, nsite-1 do begin
        ascfn = sptn2fn_thm_asi(et0, sites[i], locroot, type = 'asc')
        ascfn = sprepfile(ascfn, locroot, remroot)
        vroot = strjoin(['thg',type,sites[i]],'_')
        rr = 1
        cdfid = cdf_open(ascfn)
        cdf_varget, cdfid, vroot+'_mlon', mlon, rec_start = rr, /zvariable
        cdf_varget, cdfid, vroot+'_mlat', mlat, rec_start = rr, /zvariable
        cdf_varget, cdfid, vroot+'_elev', elev, rec_start = rr, /zvariable
        if type eq 'ast' then begin
            elev = reform(elev, [npx,npx])
        endif else begin
            mlon = reform(mlon[1,*,*],[npx+1,npx+1])
            mlat = reform(mlat[1,*,*],[npx+1,npx+1])
        endelse
        ascinfo[i].site = sites[i]
        ascinfo[i].mlon = mlon
        ascinfo[i].mlat = mlat
        ascinfo[i].elev = elev
    endfor
    
    ; determine how many loops are needed.
    etrs = [etr[0]-(etr[0] mod tmsc),etr[1]-(etr[1] mod tmsc)]
    if (etr[1] mod tmsc) eq 0 then etrs[1] -= tmsc
    etrs = smkarthm(etrs[0],etrs[1],tmsc,'dx')
    netr = n_elements(etrs)

    ; mosaic image.
    cnt = dblarr(isz,isz/2)
    del = 0.5*isz
    uts = []        ; ut.
    mos = []        ; mosaic images.
    midns = []      ; midnight mlt.
    tmo = fltarr(isz,isz/2)
    
    foreach et0, etrs do begin
        print, 'processing for '+sfmepoch(et0)+' ...'
        ; for each loop, prepare asi files for sites.
        asifns = sptn2fn_thm_asi(et0, sites, locroot, type = type)
        for i = 0, nsite-1 do asifns[i] = sprepfile(asifns[i], locroot, remroot)
        idx = where(asifns ne '', cnt)
        if cnt eq 0 then continue   ; no file for this loop.
        tsites = sites[idx]
        ntsite = cnt
        asifns = asifns[idx]
        cdfids = lonarr(ntsite)
        for i = 0, ntsite-1 do cdfids[i] = cdf_open(asifns[i])
        
        idx = where(ets ge et0 and ets lt et0+tmsc, cnt)
        if cnt eq 0 then continue   ; no epoch for this loop.
        tets = ets[idx]
        
        ; prepare the epoch for all sites.
        utptrs = ptrarr(ntsite, /allocate_heap)
        for i = 0, ntsite-1 do begin
            vname = strjoin(['thg',type,tsites[i]],'_')+'_time'
            cdf_control, cdfids[i], variable = vname, get_var_info = tmp
            cdf_varget, cdfids[i], vname, tmp, rec_count = tmp.maxrec+1
            *utptrs[i] = reform(tmp)
        endfor
        
        ; loop each epoch.
        foreach ttet, tets do begin
            print, sfmepoch(ttet)+' ...'
            ttut = 0.001d*ttet-62167219200d
            ; refresh mosaic and count.
            imf = fltarr(isz,isz)
            imc = fltarr(isz,isz)
            ; loop each site.
            for i = 0, ntsite-1 do begin
                rec = where(*utptrs[i] eq ttut, cnt)
                if cnt eq 0 then continue   ; no asi for this site at ttet.
                
                ; load mlon/mlat/elev.
                idx = where(ascinfo.site eq tsites[i])
                mlon = ascinfo[idx].mlon
                mlat = ascinfo[idx].mlat
                elev = ascinfo[idx].elev
                
                ; calc midnight.
                midn = interpol(cmlons, cmidns, (ttet/dymsc mod 1)*24, /nan, /spline)
                midn = midn[0]
                
                ; read raw image.
                vname = strjoin(['thg',type,tsites[i]],'_')
                cdf_varget, cdfids[i], vname, tmp, rec_start = rec[0]
                img = double(tmp)
                
                ; remove edge.
                if type eq 'asf' then begin
                    edge = where(elev lt minelev or ~finite(elev))
                    img[edge] = mean(img[0:10,0:10], /nan)
                    ; crude corner average subtraction. use thg_asf_site_offset?
                    img = img - img[0] > 0
                endif
                ; scale luminosity to color index. adopted from thm_asi_merge_mosaic.
                img *= 64D/(median(img) > 1)
                img <= top
                ; apply true weight.
                if keyword_set(weight) then multi = sin(elev*!dtor) $
                else multi = bytarr(npx,npx)+1b
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
            
            ; normalize use count.
            imf /= imc
            idx = where(finite(imf,/nan))
            imf = byte(imf)
            imf[idx] = 0b
            
            ; store image.
            uts = [uts,ttut]
            midns = [midns,midn]
            tmo += imf[*,0:isz/2]
            mos = [[[temporary(mos)]],[[imf[*,0:isz/2]]]]
        endforeach
        
        ; release epoch ptrarr.
        for i = 0, ntsite-1 do ptr_free, utptrs[i]
        ; close cdfids.
        for i = 0, ntsite-1 do cdf_close, cdfids[i]
    endforeach

    ; shrink data.
    pxidx = where(tmo ne 0,npxid)
    nrec = n_elements(uts)
    mos1 = bytarr(nrec,npxid)
    for i = 0, nrec-1 do mos1[i,*] = (mos[*,*,i])[pxidx]
    mos = temporary(mos1)

    ; make mlat, mlt. don't save mlon b/c it is changing, mlt+midn -> mlon.
    xc = smkarthm(-1d,1,isz,'n') # (dblarr(isz)+1)
    yc = smkarthm(-1d,1,isz,'n') ##(dblarr(isz)+1)
    r = sqrt(xc^2+yc^2)
    t = atan(yc, xc)
    mlat = 90-r*(90-minlat)         ; in deg.
    mlt = (t/!dpi*12+24) mod 24     ; in hour.
    mlat = (rotate(mlat,3))[*,0:isz/2]
    mlt = (rotate(mlt,3))[*,0:isz/2]
    mlat = mlat[pxidx]
    mlt = mlt[pxidx]
    
    ; save data.
    ; uts, mos, cmidns, minlat,  mlt, mlat.
    locroot = shomedir()
    fn = locroot+'/thg_'+type+'_mosaic_'+sfmepoch(etr[0],'YYYY_MMDD_hhmm')+'.cdf'
    print, 'saving data to '+fn+' ...'
    cdfid = cdf_create(fn, /clobber)
    cdf_compression, cdfid, set_compression = 5, set_gzip_level = 9
    ; create vatt.
    attid = cdf_attcreate(cdfid, 'FIELDNAM', /variable_scope)
    attid = cdf_attcreate(cdfid, 'UNITS', /variable_scope)
    attid = cdf_attcreate(cdfid, 'DEPEND_0', /variable_scope)
    attid = cdf_attcreate(cdfid, 'DEPEND_1', /variable_scope)
    attid = cdf_attcreate(cdfid, 'DEPEND_2', /variable_scope)
    attid = cdf_attcreate(cdfid, 'DEPEND_3', /variable_scope)
    attid = cdf_attcreate(cdfid, 'DEPEND_4', /variable_scope)
    attid = cdf_attcreate(cdfid, 'CATDESC', /variable_scope)
    ; create var.
    vname = 'time' & dimvary = 0 & var = transpose(uts)
    extra = create_struct('cdf_epoch',1, 'recvary',1, 'zvariable',1, $
        'dimension',1)
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'ut'
    cdf_attput, cdfid, 'UNITS', vname, 'sec'
    cdf_attput, cdfid, 'CATDESC', vname, 'unit time in sec'
    
    vname = 'thg_mosaic' & dimvary = 1. & var = transpose(temporary(mos))
    extra = create_struct('cdf_uint1',1, 'recvary',1, 'zvariable',1, $
        'dimensions',n_elements(pxidx))
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_compression, cdfid, variable = vname, $
        set_compression = 5, set_gzip_level = 9
    cdf_varput, cdfid, vname, var & var = 0
    cdf_attput, cdfid, 'FIELDNAM', vname, 'thg_mosaic'
    cdf_attput, cdfid, 'UNITS', vname, 'count'
    cdf_attput, cdfid, 'DEPEND_0', vname, 'time'
    cdf_attput, cdfid, 'DEPEND_1', vname, 'image_size'
    cdf_attput, cdfid, 'DEPEND_2', vname, 'pixel_index'
    cdf_attput, cdfid, 'DEPEND_3', vname, 'mlt'
    cdf_attput, cdfid, 'DEPEND_4', vname, 'mlat'
    cdf_attput, cdfid, 'CATDESC', vname, 'thg mosaic in mlt/mlat coord'
    
    vname = 'midn' & dimvary = 0 & var = transpose(midns)
    extra = create_struct('cdf_float',1, 'recvary',1, 'zvariable',1, $
        'dimensions',1)
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'midn'
    cdf_attput, cdfid, 'UNITS', vname, 'hr'
    cdf_attput, cdfid, 'DEPEND_0', vname, 'time'
    cdf_attput, cdfid, 'CATDESC', vname, 'mlt of midnight'
    
    vname = 'mlt' & dimvary = 1 & var = mlt
    extra = create_struct('cdf_float',1, 'recvary',0, 'zvariable',1, $
        'dimensions',size(var,/dimensions))
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'mlt'
    cdf_attput, cdfid, 'UNITS', vname, 'hr'
    cdf_attput, cdfid, 'CATDESC', vname, 'pixel mlt'
    
    vname = 'mlat' & dimvary = 1 & var = mlat
    extra = create_struct('cdf_float',1, 'recvary',0, 'zvariable',1, $
        'dimensions',size(var,/dimensions))
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'mlat'
    cdf_attput, cdfid, 'UNITS', vname, 'deg'
    cdf_attput, cdfid, 'CATDESC', vname, 'pixel mlat'
    
    vname = 'image_size' & dimvary = 1 & var = [isz,isz/2+1]
    extra = create_struct('cdf_int4',1, 'recvary',0, 'zvariable',1, $
        'dimensions',size(var,/dimensions))
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'image size in pixel'
    cdf_attput, cdfid, 'CATDESC', vname, '[xsize,ysize]'
    
    vname = 'pixel_index' & dimvary = [1,1] & var = pxidx
    extra = create_struct('cdf_uint4',1, 'recvary',0, 'zvariable',1, $
        'dimensions',size(var,/dimensions))
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'pixel index'
    cdf_attput, cdfid, 'CATDESC', vname, 'index of non-zero pixels'
    
    vname = 'minlat' & dimvary = 0 & var = minlat
    extra = create_struct('cdf_float',1, 'recvary',0, 'zvariable',1, $
        'dimensions',1)
    varid = cdf_varcreate(cdfid, vname, dimvary, _extra = extra)
    cdf_varput, cdfid, vname, var
    cdf_attput, cdfid, 'FIELDNAM', vname, 'minlat'
    cdf_attput, cdfid, 'UNITS', vname, 'deg'
    cdf_attput, cdfid, 'CATDESC', vname, 'min latitude in deg'
    
    cdf_close, cdfid
    
    return, fn
end

type = 'asf'

;etr = stoepoch(['2013-05-01/04:00','2013-05-01/10:00'])
;etr = stoepoch(['2013-05-01/07:00','2013-05-01/07:00:06'])
;sites = ['tpas','atha']


etr = stoepoch(['2013-04-14/04:00','2013-04-14/10:00'])
;etr = stoepoch(['2013-04-14/08:30','2013-04-14/08:30:03'])
sites = ['tpas','gill']

fn = read_thm_asi_mosaic(etr, sites, type=type, /weight)

cdfs = scdfread(fn)

uts = (*cdfs[0].value)
mos = (*cdfs[1].value)
midn= (*cdfs[2].value)
mlt = (*cdfs[3].value)
mlat= (*cdfs[4].value)
imgsz = (*cdfs[5].value)
pxidx = (*cdfs[6].value)
minlat = (*cdfs[7].value)

img = bytarr(imgsz)
nrec = n_elements(uts)

sgwopen, 1, xsize = imgsz[0], ysize = imgsz[1]
loadct, 1
for i = 0, nrec-1 do begin
    timg = img
    timg[pxidx] = mos[i,*]
    sgtv, timg
endfor
sgwclose
end
