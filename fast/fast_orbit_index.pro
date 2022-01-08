;+
; Grab orbit, orbit start time, ees time range and ies time range.
; Print the info out to file and save the info to 'sav' file.
;-
pro fast_orbit_index, orbfn, datadir, outdir = outdir
    compile_opt idl2
    
    ; path setting module.
    spawn, 'hostname', usrhost
    case !version.os_family of
        'unix': usrhost = getenv('USER')+'@'+usrhost
        'Windows': usrhost = getenv('username')+'@'+usrhost
    endcase
    ; read definitive or predicted orbit file.
    if n_elements(orbfn) eq 0 then begin
        case usrhost of
            'sheng@XpsMintv': $
                orbfn = sdiskdir('Research')+'/data/fast_cdf/orbit/definitive'
            'Sheng@Xps': $
                orbfn = sdiskdir('Research')+'/data/fast_cdf/orbit/definitive'
            'sheng@dhcp-131.spa.umn.edu': $
                orbfn = sdiskdir('Tian')+'/data/fast_cdf/orbit/definitive'
            'shengt@strelka': $
                orbfn = '/data1/fast_oa/orbit/definitive'
            'shengt@mushka': $
                orbfn = '/data1/fast_oa/orbit/definitive'
        endcase
    endif
    ; data root directory. file pattern is
    ; <datadir>/ees/fa_k0_ees_orbit_v??.cdf.
    ; <datadir>/ies/fa_k0_ies_orbit_v??.cdf.
    if n_elements(datadir) eq 0 then begin
        case usrhost of
            'sheng@XpsMintv': $
                datadir = sdiskdir('Research')+'/data/fast_cdf/'
            'Sheng@Xps': $
                datadir = sdiskdir('Research')+'/data/fast_cdf/'
            'sheng@dhcp-131.spa.umn.edu': $
                datadir = sdiskdir('Tian')+'/data/fast_cdf/'
            'shengt@strelka': $
                datadir = '/data1/fast_cdf/'
            'shengt@mushka': $
                datadir = '/data1/fast_cdf/'
        endcase
    endif
    ; filename to save orbit index variable.
    if n_elements(outdir) eq 0 then outdir = '~/'
    
    ; [orb#,orbbegtime,iesbegtime,iesendtime,eesbegtime,eesendtime].
    norb = file_lines(orbfn)/7
    orbit_index = dblarr(norb,6)
    
    ; read definitive orbit file.
    ; 7 lines per orbit. Use the 3rd line, which is
    ; ORBIT: 37 EPOCH: 1996 237 19:08:59.476280
    if file_test(orbfn) eq 0 then message, 'orbit file dooes not exist ...'
    print, 'read orbit epoch ...'
    openr, lun, orbfn, /get_lun & line = ''
    for i = 0, norb-1 do begin
        skip_lun, lun, 2, /lines
        readf, lun, line
        tmp = strsplit(strcompress(line),' :.',/extract)
        tmp = float(tmp[[1,3,4,5,6,7,8]])
        orb = tmp[0] & yr = tmp[1] & doy = tmp[2]
        hr = tmp[3] & mi = tmp[4] & sc = tmp[5] & ms = tmp[6]*1e-3
        skip_lun, lun, 4, /lines
        ; convert doy to month/date.
        days = [0,31,59,90,120,151,181,212,243,273,304,334,365]
        if yr mod 400 eq 0 or yr mod 100 ne 0 and yr mod 4 eq 0 then $
            days[2:*] += 1
        mo = (where(doy le days))[0]
        dy = doy-days[mo-1]
        cdf_epoch, et, yr, mo, dy, hr, mi, sc, ms, /compute_epoch
        orbit_index[i,0:1] = [orb, et]
    endfor
    free_lun, lun
    
    ; read each orbit's ees data available time.
    if file_test(datadir) eq 0 then message, 'data directory does not exist ...'
    print, 'read ees and ies time range ...'
    !quiet = 1
    for i = 0, norb-1 do begin
        orbstr = string(i,format='(I05)') & vname = 'Epoch'
        print, 'processing orbit '+orbstr+' ...'
        eesfn = datadir+'ees/fa_k0_ees_'+orbstr+'_v??.cdf'
        esafn = file_search(eesfn, count=cnt)
        if cnt gt 0 then begin
            cdfid = cdf_open(esafn[cnt-1])
            cdf_control, cdfid, variable = vname, /zvariable, $
                get_var_info = vinfo
            maxrec = vinfo.maxrec
            if maxrec gt 1 then begin
                cdf_varget, cdfid, vname, tr, /zvariable, $
                    rec_start = 0, rec_interval = maxrec, rec_count = 2
                orbit_index[i,2:3] = tr
            endif
            cdf_close, cdfid
        endif
        iesfn = datadir+'ies/fa_k0_ies_'+orbstr+'_v??.cdf'
        esafn = file_search(iesfn, count=cnt)
        if cnt gt 0 then begin
            cdfid = cdf_open(esafn[cnt-1])
            cdf_control, cdfid, variable = vname, /zvariable, $
                get_var_info = vinfo
            maxrec = vinfo.maxrec
            if maxrec gt 1 then begin
                cdf_varget, cdfid, vname, tr, /zvariable, $
                    rec_start = 0, rec_interval = maxrec, rec_count = 2
                orbit_index[i,4:5] = tr
            endif
            cdf_close, cdfid
        endif
    endfor
    
    ; save orbit index.
    outfn = outdir+'fast_orbit_index.sav'
    save, orbit_index, filename = outfn
    
    ; print orbit index to file.
    norb = n_elements(orbit_index)/6
    outfn = outdir+'fast_orbit_index.dat'
    tfmt = 'yyyy-MM-dd/hh:mm:ss.fff'
    openw, lun, outfn, /get_lun
    for i = 0, norb-1 do begin
        orb = string(orbit_index[i,0], format='(I05)')
        orbt0 = sfmepoch(orbit_index[i,1],tfmt)
        if orbit_index[i,2] eq 0 then begin
            eest0 = 'nan' & eest1 = 'nan'
        endif else begin
            eest0 = sfmepoch(orbit_index[i,2],tfmt)
            eest1 = sfmepoch(orbit_index[i,3],tfmt)
        endelse
        if orbit_index[i,4] eq 0 then begin
            iest0 = 'nan' & iest1 = 'nan'
        endif else begin
            iest0 = sfmepoch(orbit_index[i,4],tfmt)
            iest1 = sfmepoch(orbit_index[i,5],tfmt)
        endelse
        printf, lun, orb, orbt0, eest0, eest1, iest0, iest1, $
            format='(A,T10,A,T37,A,T64,A,T91,A,T118,A)'
    endfor
    free_lun, lun

end
