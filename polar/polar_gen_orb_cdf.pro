;+
; Type: procedure.
; Purpose: Generate date-based position file, with file name pattern
;   po_or_def_yyyyMMdd_v09.cdf The info are: mlt, ilat, distance, pos in GSE.
;   Supposed to run once.
; Parameters:
;   ets, in, dblarr[2], opt. Epoch range, default to 1996-02-27/00:00 to
;       2008-06-14/00:00. That's when polar has orbital info on cdaweb.
; Keywords: none.
; Notes:
;   * The output file has pattern 'po_or_def_yyyymmdd_v01.cdf'. The vars are:
;       'Epoch': epoch of data.
;       'ilat': invariant latitude in degree, has sign.
;       'mlt': magnetic local time in hour.
;       'dis': distance in Re.
;       'Re': 6378.137 km.
;       'pos_gse': position in GSE coord in km.
;       'mlat': magnetic latitude in degree, has sign.
; Dependence: slib.
; History:
;   2014-05-24, Sheng Tian, create.
;-
pro polar_gen_orb_cdf, ets, outdir, overwrite = overwrite
    compile_opt idl2
    
    re = 6378.137d
    if n_elements(ets) eq 0 then ets = stoepoch(['1996-02-27','2008-06-14'])
    ets = stoepoch(ets)
    et0 = sepochfloor(ets[0]) & et1 = sepochfloor(ets[1])
    det = 86400000d
    
    if n_elements(outdir) eq 0 then outdir = shomedir()+'/sdata/polar/orbit/'
    if ~file_test(outdir,/directory) then file_mkdir, outdir
    
    for i = et0, et1, det do begin
        data = sread_polar_pos(t=[i,i+det], 'cdaweb')
        ; write cdf file.
        date = sfmepoch(i,'yyyyMMdd')
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
        tval = data.epoch
        varid = cdf_varcreate(cdfid, vname, /zvariable, /cdf_epoch)
        cdf_varput, cdfid, vname, tval
        ; write ilat.
        vname = 'ilat'
        tval = acos(sqrt(1D/data.l_shell))*(180d/!dpi)   ; in degree.
        idx = where(data.mag_latitude lt 0, cnt)
        if cnt gt 0 then tval[idx] = -tval[idx]
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'deg'
        ; write mlt.
        vname = 'mlt'
        tval = data.edmlt_time
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
        tval = sqrt(total(data.gse_pos^2,2))*(1d/re)
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'Re'
        ; write pos.
        vname = 'pos_gse'
        tval = transpose(data.gse_pos)
        varid = cdf_varcreate(cdfid, vname, 1, dimensions = 3, /zvariable, $
            /cdf_double)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'km'
        ; write mlat.
        vname = 'mlat'
        tval = data.mag_latitude
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'deg'
        ; finish and close.
        cdf_close, cdfid
    endfor
end