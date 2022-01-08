;+
; Type: procedure.
; Purpose: Find FAST and Polar conjunction candidate. Only on umn computers.
;   Default is find satellite in: mlt: [6,18]; * ilat: [50,90].
;   Other thresholds are: * AE, Dst, Polar distance.
;   For all thresholds, [0,0] means no check. The relations are 'and', NOT 'or'.
;       Results are written to, for example, 
;   '~/conj_fa_po_mlt_06-18_ilat_50-90_ae_300-600_dst_0-500_dist_0-3.dat',
;   '~/conj_fa_po_mlt_10-14_ilat_60-80_ae_300-900_dst_0-0_dist_0-0.dat'.
; Parameters:
;   tr, in, dblarr[2], req. Time range in unix time.
;       FAST has ESA data from 
;       orbit 330 (start from 1996-09-20 21:47 UT) to 
;       orbit 51315 (start from 2009-04-30 09:23 UT).
; Keywords: none.
; Notes: none.
; Dependence: none.
; History:
;   2012-07-28, Sheng Tian, create.
;   2013-03-23, Sheng Tian, revise.
;-
pro sdt_get_po_pos, fn, t0, ilat, mlt, dist
    compile_opt idl2

    ; preparation.
    re = 6356.7523      ; polar radius, 6356.7523 km.

    cdfid = cdf_open(fn)

    ; nrec.
    cdf_control, cdfid, variable = 'Epoch', get_var_info = varinfo
    nrec = varinfo.maxrec+1

    ; read epoch.
    vname = 'Epoch'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    t0 = sfmepoch(reform(tval[0,*]),'unix')

    ; read mlt.
    vname = 'EDMLT_TIME'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    mlt = reform(tval[0,*])

    ; get ilat: from L shell, get sign from mlat.
    vname = 'L_SHELL'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    ilat = acos(sqrt(1D/reform(tval[0,*])))*(1D/!dtor)      ; in deg.
    vname = 'MAG_LATITUDE'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    mlat = reform(tval[0,*])
    idx = where(mlat lt 0)
    if idx[0] ne -1 then ilat[idx] = -ilat[idx]

    ; get dist: from gse pos, convert to unit in Re.
    vname = 'GSE_POS'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    dist = sqrt(total(tval^2, 1))*(1D/re)

    cdf_close, cdfid

end


pro sdt_get_fa_pos, fn, t0, ilat, mlt, dist, orb, uniform = uniform
    compile_opt idl2

    ; preparation.
    re = 6356.7523      ; polar radius, 6356.7523 km.

    cdfid = cdf_open(fn)

    ; nrec.
    cdf_control, cdfid, variable = 'time', get_var_info = varinfo
    nrec = varinfo.maxrec+1

    ; read time.
    vname = 'time'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    times = reform(tval)

    ; read mlt.
    vname = 'mlt'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    mlt = reform(tval)

    ; read ilat.
    vname = 'ilat'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    ilat = reform(tval)

    ; get dist: from altitude.
    vname = 'alt'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    dist = reform(tval)*(1D/re)+1

    ; get orbit.
    vname = 'orbit'
    cdf_varget, cdfid, vname, tval, rec_count = nrec
    orb = reform(tval)

    cdf_close, cdfid

    if ~keyword_set(uniform) then return

    ; remove time overlap.
    idx0 = sort(times)
    idx1 = uniq(times[idx0])
    times = (times[idx0])[idx1]
    mlt = (mlt[idx0])[idx1]
    ilat = (ilat[idx0])[idx1]
    dist = (dist[idx0])[idx1]
    orb = (orb[idx0])[idx1]

    ; treat the gaps between orbit.
    nrec0 = n_elements(times)       ; original nrec.
    nrec = 1441
    flag = bytarr(nrec)
    tval = 86400D*floor(times[0]/86400)+dindgen(nrec)*60
    for i = 0, nrec0-1 do begin
        idx = where(times[i] eq tval)
        flag[idx] = 1B
    endfor
    idx1 = where(flag eq 1)
    idx0 = where(flag eq 0)
    ; time.
    t0 = tval[0:nrec-2]
    ; mlt.
    tval = findgen(nrec)
    tval[idx1] = mlt
    if idx0[0] ne -1 then tval[idx0] = !values.f_nan
    mlt = tval[0:nrec-2]
    ; ilat.
    tval = findgen(nrec)
    tval[idx1] = ilat
    if idx0[0] ne -1 then tval[idx0] = !values.f_nan
    ilat = tval[0:nrec-2]
    ; dist.
    tval = findgen(nrec)
    tval[idx1] = dist
    if idx0[0] ne -1 then tval[idx0] = !values.f_nan
    dist = tval[0:nrec-2]
    ; orbit.
    tval = findgen(nrec)
    tval[idx1] = orb
    if idx0[0] ne -1 then tval[idx0] = !values.f_nan
    orb = tval[0:nrec-2]

end


pro sdt_get_ae, fn, tr
    compile_opt idl2

    restore, filename = fn
    t0 = ind_data.time
    ae = ind_data.data
    if n_elements(tr) ne 0 then idx = where(t0 ge tr[0] and t0 le tr[1]) $
    else idx = indgen(n_elements(t0))
    if idx[0] ne -1 then store_data, 'ae', data = {x:t0[idx], y:ae[idx]}

end


pro sdt_get_dst, fn, tr
    compile_opt idl2

    restore, filename = fn
    t0 = ind_data.time
    ae = ind_data.data
    if n_elements(tr) ne 0 then idx = where(t0 ge tr[0] and t0 le tr[1]) $
    else idx = indgen(n_elements(t0))
    if idx[0] ne -1 then store_data, 'dst', data = {x:t0[idx], y:ae[idx]}

end


pro plot_polar_fast_line_summary, tr
    compile_opt idl2

    ; check time range.
    if n_elements(tr) eq 0 then message, 'no time range ...'
    eps = sbreaktr(stoepoch(tr, 'unix'))        ; in epoch.
    nep = n_elements(eps)*0.5

    ; data source.
    ; polar orbit:      '/data/polar/oa/orbit/po_or_def_yyyymmdd_v*.cdf'
    ; fast orbit time:  '/data1/fast_oa/orbit/definitive'
    ; ae index:         '/data1/indices/AE/ae_data_yyyy_mm.dat'
    ; dst index:        '/data1/indices/DST/dst_data_yyyy_mm.dat'
    poptn = '/data/polar/oa/orbit/po_or_*_yyyymmdd_v*.cdf'
    faptn = '~/works/data/fast/orbit/yyyy/fa_or_def_yyyymmdd_v01.cdf'
    aeptn = '/data1/indices/AE/ae_data_yyyy_mm.dat'
    dstptn= '/data1/indices/DST/dst_data_yyyy_mm.dat'

    ;plot settings.
    device, decomposed = 0
    loadct2, 42

    ; loop each day, assuming whole day.
    for i = 0, nep-1 do begin
        ep0 = eps[*,i]                  ; in epoch.
        tr0 = sfmepoch(ep0, 'unix')     ; in unix time.
        nrec = 1440
        t0 = tr0[0] + findgen(nrec)*60

        ; read polar position: ilat, mlt, dist.
        pofn = file_search(sptn2fn(poptn, ep0[0]))
        sdt_get_po_pos, pofn[0], pot0, poilat, pomlt, podist

        ; read fast position: ilat, mlt, dist.
        fafn = file_search(sptn2fn(faptn, ep0[0]))
        if fafn[0] eq '' then continue      ; 1997-05-18 no data.
        sdt_get_fa_pos, fafn[0], fat0, failat, famlt, fadist, faorb, /uniform

        ; read ae index.
        aefn = sptn2fn(aeptn, ep0[0])
        sdt_get_ae, aefn[0], tr0

        ; read dst index.
        dstfn = sptn2fn(dstptn, ep0[0])
        sdt_get_dst, dstfn[0], tr0

        ; ilat.
        ilat = abs([[poilat],[failat]])
        idx = where(poilat*failat lt 0)     ; different hemisphere.
        if idx[0] ne -1 then ilat[idx,*] = !values.d_nan
        store_data, 'ilat', data = {x:t0, y:ilat}
        ylim, 'ilat', 60, 90, 0

        ; hemisphere.
        hem = (failat ge 0)*2-1D
        if idx[0] ne -1 then hem[idx] = !values.d_nan
        store_data, 'hem', data = {x:t0, y:hem}
        ylim, 'hem', -2, 2, 0

        ; mlt.
        store_data, 'mlt', data = {x:t0, y:[[pomlt],[famlt]]}
        ylim, 'mlt', 6, 18, 0

        ; dist.
        store_data, 'dist', data = {x:t0, y:podist}

        ; orbit.
        store_data, 'orbit', data = {x:t0, y:faorb}
        ylim, 'orbit', min(faorb, /nan)-1, max(faorb, /nan)+1, 0

        vars = ['ilat','mlt']
        options, vars, 'labels', ['polar','fast']
        options, vars, 'colors', [2,4]

        ; plot data.
        vars = ['ilat','mlt','hem','orbit','dist','ae','dst']
        tplot, vars, trange = tr0
        stop

    endfor

end
