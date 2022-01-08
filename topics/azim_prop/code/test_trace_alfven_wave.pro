;---Load basic data.
    id = '2014_0828_10'
    info_var = id+'_event_info'
    if tnames(info_var) eq '' then _2014_0828_10_load_data

    test = 1

;---Settings.
    einfo = get_var_data(info_var)
    pres = ['thd','the','tha']+'_'
    trace_times = time_double('2014-08-28/'+['10:10','10:10','10:20'])
    trace_dir = -1
    trace_height = 100. ; km.
    re = 6378.  ; km.
    r0 = trace_height/re+1
    model = 't89'
    va0 = 1000. ; km/s.

    foreach pre0, pres, ii do begin
        rvar = pre0+'r_gsm'
        time = trace_times[ii]
        rgsm = get_var_data(rvar, at=time)
        tilt = geopack_recalc(time)
        par = 2
        geopack_trace, rgsm[0],rgsm[1],rgsm[2], trace_dir, par, tilt=tilt, r0=r0, /igrf, /t89, $
            fx,fy,fz, fline=fline
        ; fline in [n,3], in Re, in GSM.
        nstep = n_elements(fline)/3
        bmag = fltarr(nstep)
        for jj=0, nstep-1 do begin
            ; Need to recalc as time elapse???
            geopack_igrf_gsm, fline[jj,0],fline[jj,1],fline[jj,2], bx,by,bz
            geopack_t89, par, fline[jj,0],fline[jj,1],fline[jj,2], dbx,dby,dbz
            bmag[jj] = snorm([bx,by,bz]+[dbx,dby,dbz])
        endfor

        duration = 0
        length = 0
        for jj=0, nstep-2 do begin
            tva = va0*(bmag[jj]+bmag[jj+1])*0.5/bmag[0]
            tlength = snorm(fline[jj,*]-fline[jj+1,*])
            tduration = tlength/tva*re
            duration += tduration
            length += tlength
;            print, tva
;            print, tlength
;            print, tduration
        endfor
        print, pre0
        print, duration
        print, length
    endforeach
end
