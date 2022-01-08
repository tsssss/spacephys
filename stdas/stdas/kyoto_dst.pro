
pro kyoto_dst, tr, t0, dat

    if n_elements(tr) eq 0 then get_timespan, tr
    
    fn = file_dailynames(file_format = 'YYYYMM/dstyyMM', trange = tr, $
        times = times, /unique)+'.for.request'
    source = file_retrieve(/struct)
    source.local_data_dir = root_data_dir()+'geom_indices/kyoto/dst/'
    
    ; realtime dst.
    source.remote_data_dir = 'http://wdc.kugi.kyoto-u.ac.jp/dst_realtime/'
    lfns = file_retrieve(fn,_extra = source)
    
    ; provisional dst.
    source.remote_data_dir = 'http://wdc.kugi.kyoto-u.ac.jp/dst_provisional/'
    lfns = file_retrieve(fn,_extra = source)
    
    ; final dst.
    source.remote_data_dir = 'http://wdc.kugi.kyoto-u.ac.jp/dst_final/'
    lfns = file_retrieve(fn,_extra = source)
    
    nlfn = n_elements(lfns) & tmp = ''
    dat = !null & t0 = !null
    for i = 0, nlfn-1 do begin
        tfn = lfns[i]
        if ~file_test(tfn, /regular) then continue
        openr, lun, tfn, /get_lun
        nline = file_lines(tfn)-2
        tdat = fltarr(nline,24)
        tt0 = dblarr(nline,24)
        for j = 0, nline-1 do begin
            readf, lun, tmp
            tt0[j,*] = time_double(strmid(tmp,14,2)+strmid(tmp,3,2)+$
                '-'+strmid(tmp,5,2)+'-'+strmid(tmp,8,2))+dindgen(24)*3600D
            tdat[j,*] = float(strmid(tmp, indgen(24)*4+20,4))
        endfor
        free_lun, lun
        dat = [dat,reform(transpose(tdat),nline*24)]
        t0 = [t0,reform(transpose(tt0),nline*24)]
    endfor
    
    idx = where(t0 ge tr[0] and t0 le tr[1], cnt)
    dat = dat[idx]
    t0 = t0[idx]
    
    idx = where(dat eq 9999, cnt)
    if cnt ge 0 then dat[idx] = !values.f_nan
    
end

tr = time_double(['2014-05-08','2014-05-13'])
kyoto_dst, tr, t0, dat
store_data, 'dst', t0, dat
tplot, 'dst', trange = tr
end