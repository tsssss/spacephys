; some date, polar orbit data at cdaweb has data gap, use sscweb data as save.
; fn should contain 1 day of polar data from sscweb.

pro polar_fill_orb_gap, fn, outdir, overwrite = overwrite
    compile_opt idl2

    re = 6378.137d
    if n_elements(outdir) eq 0 then outdir = shomedir()+'/sdata/polar/orbit/'
    if ~file_test(outdir,/directory) then file_mkdir, outdir
    
    vars = ['Epoch','GM_LAT','GM_LCT_T','XYZ_GSE','INVAR_LAT']
    data = scdfread(fn, vars)
    et0 = *data[0].value
    mlat = *data[1].value       ; in deg, has sign.
    mlt = *data[2].value        ; in hr.
    pos_gse = *data[3].value    ; in km.
    ilat = *data[4].value       ; in deg, no sign.
    
    ets = sepochfloor(et0[0])+dindgen(1440)*60000d
    mlat = sinterpol(mlat, et0, ets)
    mlt = sinterpol(mlt, et0, ets)
    pos_gse = sinterpol(pos_gse, et0, ets)
    ilat = sinterpol(ilat, et0, ets)
    dis = sqrt(total(pos_gse^2,2))*(1d/re)
    
    date = sfmepoch(et0[0],'yyyyMMdd')
    odir = outdir+strmid(date,0,4)+'/'
    if ~file_test(odir,/directory) then file_mkdir, odir
    ofn = odir+'po_or_def_'+date+'_v01.cdf'
    if keyword_set(overwrite) then if file_test(ofn) then file_delete, ofn
    printf, -1, 'write to file: '+ofn+' ...'
    cdfid = cdf_create(ofn)
    unit = 'UNIT' & dummy = cdf_attcreate(cdfid, unit, /variable)
    fnam = 'FIELDNAM' & dummy = cdf_attcreate(cdfid, fnam, /variable)
    ; write epoch.
    vname = 'Epoch'
    tval = ets
    varid = cdf_varcreate(cdfid, vname, /zvariable, /cdf_epoch)
    cdf_varput, cdfid, vname, tval
    ; write ilat.
    vname = 'ilat'
    tval = ilat
    idx = where(mlat lt 0, cnt)
    if cnt gt 0 then tval[idx] = -tval[idx]
    varid = cdf_varcreate(cdfid, vname, /zvariable)
    cdf_varput, cdfid, vname, tval
    cdf_attput, cdfid, unit, vname, 'deg'
    ; write mlt.
    vname = 'mlt'
    tval = mlt
    varid = cdf_varcreate(cdfid, vname, /zvariable)
    cdf_varput, cdfid, vname, tval
    cdf_attput, cdfid, unit, vname, 'hr'
    ; write re.
    vname = 'Re'
    tval = re
    varid = cdf_varcreate(cdfid, vname, /zvariable, /rec_novary, $
        /cdf_double)
    cdf_varput, cdfid, vname, tval
    cdf_attput, cdfid, unit, vname, 'km'
    cdf_attput, cdfid, fnam, vname, 'Polar Earth radius'
    ; write dis.
    vname = 'dis'
    tval = dis
    varid = cdf_varcreate(cdfid, vname, /zvariable)
    cdf_varput, cdfid, vname, tval
    cdf_attput, cdfid, unit, vname, 'Re'
    ; write pos.
    vname = 'pos_gse'
    tval = transpose(pos_gse)
    varid = cdf_varcreate(cdfid, vname, 1, dimensions = 3, /zvariable, $
        /cdf_double)
    cdf_varput, cdfid, vname, tval
    cdf_attput, cdfid, unit, vname, 'km'
    ; write mlat.
    vname = 'mlat'
    tval = mlat
    varid = cdf_varcreate(cdfid, vname, /zvariable)
    cdf_varput, cdfid, vname, tval
    cdf_attput, cdfid, unit, vname, 'deg'
    ; finish and close.
    cdf_close, cdfid

end