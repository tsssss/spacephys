
pro swfaorbtr, orbs, ofn, iptn

    if n_elements(orbs) eq 0 then orbs = indgen(50986, /ulon)+330   ;[330,51315]
    norb = n_elements(orbs)

    if n_elements(iptn) eq 0 then $
        iptn = ['/data1/fast_cdf/ies/fa_k0_ies_','_v*.cdf']

    if n_elements(ofn) eq 0 then ofn = '~/fa_orb_tr.dat'
    openw, lun, ofn, /get_lun, /append
    term = -1       ; terminal lun.

    for i = 0UL, norb-1 do begin
        torb = orbs[i]
        printf, term, 'orbit: ', torb
        ifn = file_search(iptn[0]+string(torb, format = '(I05)')+iptn[1])
        if not file_test(ifn[0]) then continue      ; skip empty orbit.
        printf, term, 'reading file: ' + ifn[0]
        cdfid = cdf_open(ifn[0])    ; use the lowest version, ok for time range.
        vname = 'Epoch'
        cdf_control, cdfid, variable = vname, get_var_info = vinfo
        cdf_varget1, cdfid, vname, ep1, rec_start = 0
        cdf_varget1, cdfid, vname, ep2, rec_start = vinfo.maxrec
        cdf_epoch, ep1, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        ep1 = string(yr, mo, dy, hr, mi, sc, msc, $
            format = '(I04,"-",I02,"-",I02,"/",I02,":",I02,":",I02,".",I03)')
        cdf_epoch, ep2, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        ep2 = string(yr, mo, dy, hr, mi, sc, msc, $
            format = '(I04,"-",I02,"-",I02,"/",I02,":",I02,":",I02,".",I03)')
        printf, lun, torb, ep1, ep2, format = '(I05, 4X, A23, 4X, A23)'
        cdf_close, cdfid
    endfor

    close, lun

end
