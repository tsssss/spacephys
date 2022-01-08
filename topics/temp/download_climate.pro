
pro download_climate, yr1

    yr0 = 1870
    if n_elements(yr1) eq 0 then spawn, 'date +%Y', yr1 & yr1 = fix(yr1[0])
    
    remdir = 'http://climate.umn.edu/doc/twin_cities'
    locdir = shomedir()+'/temp/twin_cities'
    if file_test(locdir) eq 0 then file_mkdir, locdir
    
    nfile = (yr1-yr0)/10+1
    files = 'msp'+string(yr0+findgen(nfile)*10,format='(I4)')+"'s.htm"
    for i = 0, nfile-1 do $
        spawn, 'curl "'+remdir+'/'+files[i]+'" -o "'+locdir+'/'+files[i]+'"'
        
    ; read data.
    header = strarr(17)
    footer = strarr(5)
    ut = [] & tmn = [] & tmx = []
    for i = 0, nfile-1 do begin
        tfile = locdir+'/'+files[i]
        if file_test(tfile) eq 0 then continue
        nline = file_lines(tfile)-n_elements(header)-n_elements(footer)
        lines = strarr(nline)
        openr, lun, tfile, /get_lun
        readf, lun, header
        readf, lun, lines
        free_lun, lun
        
        for j = 0, nline-1 do begin
            if lines[j] eq '' then continue ; empty line.
            tmp = strsplit(lines[j],/extract)
            if tmp[0] eq 0 then continue    ; illegal start.
            if strupcase(tmp[3]) eq 'M' then continue   ; no high temp.
            if strupcase(tmp[4]) eq 'M' then continue   ; no low temp.
            tmp = double(tmp[0:4])
            cdf_epoch, ep, tmp[0], tmp[1], tmp[2], /compute_epoch
            ut = [ut,ep] & tmx = [tmx,tmp[3]] & tmn = [tmn,tmp[4]]
        endfor
    endfor
    ut = sfmepoch(ut,'unix')
    tavg = (tmx+tmn)*0.5    ; in F.
    tmx = (tmx-32)*5d/9     ; in C.
    tmn = (tmn-32)*5d/9     ; in C.
    tavg = (tavg-32)*5d/9   ; in C.
    
    tlimit, time_double('1870-01-01'), time_double('2020-01-01')
    store_data, 'tmn', ut, tmn
    store_data, 'tmx', ut, tmx
    store_data, 'tavg', ut, tavg

    ; read sunsplot number data.
    spfile = locdir+'/dayssn.dat'
    nline = file_lines(spfile)
    lines = strarr(nline)
    openr, lun, spfile, /get_lun
    readf, lun, lines
    free_lun, lun
    
    ut = [] & ssn = []
    for j = 0, nline-1 do begin
        if lines[j] eq '' then continue ; empty line.
        tmp = strsplit(lines[j],/extract)
        if tmp[2] eq '?' then continue  ; no ssn.
        ut = [ut,stoepoch(tmp[0],'YYYYMMDD')] & ssn = [ssn,double(tmp[2])]
    endfor
    ut = sfmepoch(ut,'unix')
    store_data, 'ssn', ut, ssn
    
    temp = {ut:ut, tmx:tmx, tmn:tmn, tavg:tavg, ssn:ssn}
    save, temp, filename = locdir+'/temp.sav'
end

download_climate
end