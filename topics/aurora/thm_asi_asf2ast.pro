
pro thm_asi_asf2ast, sites, ut0

    compile_opt idl2
    on_error, 0
    
    if n_elements(rootdir) eq 0 then $
        rootdir = spreproot('themis')

    et0 = stoepoch(ut0, 'unix')
    asffns = sptn2fn_thm_asi(et0, sites, rootdir, type = 'asf')
    ascfns = sptn2fn_thm_asi(et0, sites, rootdir, type = 'asc')
    astfns = sptn2fn_thm_asi(et0, sites, rootdir, type = 'ast')
    
    i = 0
    
    ; read calibration file.
    cdfid = cdf_open(ascfns[i])
    vroot = 'thg_asf_'+sites[i]
    cdf_varget, cdfid, vroot+'_bin', asfbin, /zvariable
    vroot = 'thg_ast_'+sites[i]
    cdf_varget, cdfid, vroot+'_binr', astbinr, /zvariable
    cdf_varget, cdfid, vroot+'_binc', astbinc, /zvariable
    cdf_varget, cdfid, vroot+'_elev', elev, rec_start = 1, /zvariable
    cdf_varget, cdfid, vroot+'_mlat', mlat, rec_start = 1, /zvariable
    cdf_varget, cdfid, vroot+'_mlon', mlon, rec_start = 1, /zvariable
    cdf_close, cdfid
    
    ; read asf file.
    cdfid = cdf_open(asffns[i])
    vroot = 'thg_asf_'+sites[i]
    cdf_control, cdfid, variable = vroot+'_time', get_var_info = tmp, /zvariable
    cdf_varget, cdfid, vroot+'_time', uts, rec_count = tmp.maxrec+1, /zvariable
    uts = reform(uts, /overwrite)
    tmp = min(ut0-uts, rec, /absolute)
    cdf_varget, cdfid, vroot, asfimg, rec_start = rec, /zvariable
    cdf_close, cdfid
    
    ; read ast file.
    cdfid = cdf_open(astfns[i])
    vroot = 'thg_ast_'+sites[i]
    cdf_control, cdfid, variable = vroot+'_time', get_var_info = tmp, /zvariable
    cdf_varget, cdfid, vroot+'_time', uts, rec_count = tmp.maxrec+1, /zvariable
    uts = reform(uts, /overwrite)
    tmp = min(ut0-uts, rec, /absolute)
    cdf_varget, cdfid, vroot, astimg0, rec_start = rec, /zvariable
    cdf_close, cdfid
    
    ; asf to ast.
    npx = 32
    astimg = fltarr(npx,npx)
    for j = 0, npx*npx-1 do begin
        idx = where(asfbin eq j)
        astimg[j] = mean(asfimg[idx])
    endfor
    
    ; recover ast.
    lats = mlat[sort(mlat)]
    lons = mlon[sort(mlon)]
    lats = lats[uniq(lats)]
    lons = lons[uniq(lons)]
    nlat = n_elements(lats)-1
    nlon = n_elements(lons)-1
    maps = intarr(nlon, nlat)
    for i = 0, npx*npx-1 do begin
        ; [lower left, upper left, upper right, lower right].
        plat = reform(mlat[*,i])
        plon = reform(mlon[*,i])
        pi = where(lats eq min(plat), cnt)
        pj = where(lons eq min(plon), cnt)
        maps[pj,pi] = i
    endfor
    astimg = astimg[maps]
    window, 0, xsize = 500*3, ysize = 500
    tvscl, congrid(rotate(astimg,2), 500, 500), 0
    tvscl, congrid(asfimg, 500, 500), 1
    tvscl, congrid(astimg0[maps], 500, 500), 2
    stop

end

ut = time_double('2008-03-09/06:30:00')
thm_asi_asf2ast, 'chbg', ut
end