
pro thm_asi_show_ast, sites, ut0

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
    vroot = 'thg_ast_'+sites[i]
    cdf_varget, cdfid, vroot+'_binr', astbinr, /zvariable
    cdf_varget, cdfid, vroot+'_binc', astbinc, /zvariable
    cdf_varget, cdfid, vroot+'_elev', elev, rec_start = 1, /zvariable
    cdf_varget, cdfid, vroot+'_mlat', mlat, rec_start = 1, /zvariable
    cdf_varget, cdfid, vroot+'_mlon', mlon, rec_start = 1, /zvariable
    cdf_close, cdfid
    
    ; read ast file.
    cdfid = cdf_open(astfns[i])
    vroot = 'thg_ast_'+sites[i]
    cdf_control, cdfid, variable = vroot+'_time', get_var_info = tmp, /zvariable
    cdf_varget, cdfid, vroot+'_time', uts, rec_count = tmp.maxrec+1, /zvariable
    uts = reform(uts, /overwrite)
    tmp = min(ut0-uts, rec, /absolute)
    cdf_varget, cdfid, vroot, astimg, rec_start = rec, /zvariable
    cdf_close, cdfid
    npx = 32

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
    mlat = mlat[maps]
    mlon = mlon[maps]
    
    stop
end

ut = time_double('2013-04-14/07:00')
thm_asi_show_ast, 'tpas', ut
end