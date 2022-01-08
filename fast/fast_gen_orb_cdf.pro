;+
; Type: procedure.
; Purpose: Generate date-based position file, with file name pattern 
;   fa_or_def_yyyyMMdd_v02.cdf. The info are: mlt, ilat, distance, 
;   orbit number, pos in GSE, calculated from get_fa_orbit.
;   Supposed to run once.
; Parameters:
;   ets, in, dblarr[2], opt. Epoch range, default 1996-09-20/00:00 to 
;       2009-05-01/00:00, corresponds to the date of [330,51315], 
;       that's when fast has esa data.
; Keywords: none.
; Notes:
;   * The output file has pattern 'fa_or_def_yyyymmdd_v02.cdf'. The vars are:
;       'Epoch': epoch of data.
;       'ilat': invariant latitude in degree, has sign.
;       'mlt': magnetic local time in hour.
;       'dis': distance in Re.
;       'Re': 6378.137 km.
;       'pos_gse': fast position in GSE coord in km.
;       'orbit': orbit number.
; Dependence: slib.
; History:
;   2013-04-09, Sheng Tian, create.
;-
pro fast_gen_orb_cdf, ets, outdir, overwrite = overwrite
    compile_opt idl2

    re = 6378.137d
    if n_elements(ets) eq 0 then ets = stoepoch(['1996-09-20','2009-05-01'])
    et0 = sepochfloor(ets[0]) & et1 = sepochfloor(ets[1])
    uts = sfmepoch([et0,et1],'unix')
    dut = 86400d & uts[1] -= dut
    if n_elements(outdir) eq 0 then outdir = '~/works/sdata/fast/orbit/'
    if ~file_test(outdir,/directory) then file_mkdir, outdir

    for i = uts[0], uts[1], dut do begin
        get_fa_orbit, i, i+dut-60, delta_t = 60, struc = data, $
            /no_store, /definitive
        ; write cdf file.
        date = sfmepoch(stoepoch(i,'unix'),'yyyyMMdd')
        odir = outdir+strmid(date,0,4)+'/'
        if ~file_test(odir,/directory) then file_mkdir, odir
        ofn = odir+'fa_or_def_'+date+'_v02.cdf'
        if keyword_set(overwrite) then if file_test(ofn) then file_delete, ofn
        printf, -1, 'write to file: '+ofn+' ...'
        cdfid = cdf_create(ofn)
        unit = 'UNIT' & dummy = cdf_attcreate(cdfid, unit, /variable)
        fnam = 'FIELDNAM' & dummy = cdf_attcreate(cdfid, fnam, /variable)
        ; write epoch.
        vname = 'Epoch'
        tval = stoepoch(i,'unix')+dindgen(n_elements(data.time))*60000d
        varid = cdf_varcreate(cdfid, vname, /zvariable, /cdf_epoch)
        cdf_varput, cdfid, vname, tval
        ; write ilat.
        vname = 'ilat'
        tval = data.ilat
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'deg'
        ; write mlt.
        vname = 'mlt'
        tval = data.mlt
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'hr'
        ; write re.
        vname = 'Re'
        tval = re
        varid = cdf_varcreate(cdfid, vname, /zvariable, /rec_novary)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'km'
        cdf_attput, cdfid, fnam, vname, 'Polar Earth radius'
        ; write dis.
        vname = 'dis'
        tval = data.alt*(1d/re)+1
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'Re'
        ; write pos.
        vname = 'pos_gse'
        tval = transpose(data.fa_pos)
        varid = cdf_varcreate(cdfid, vname, 1, dimensions = 3, /zvariable, $
            /cdf_double)
        cdf_varput, cdfid, vname, tval
        cdf_attput, cdfid, unit, vname, 'km'
        ; write orb.
        vname = 'orbit'
        tval = data.orbit
        varid = cdf_varcreate(cdfid, vname, /zvariable)
        cdf_varput, cdfid, vname, tval
        ; finish and close.
        cdf_close, cdfid
    endfor
end
