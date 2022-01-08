
type = 'detect_step'
type = 'all_e'

ifn = shomedir()+'/psbl_de_32hz/asi_availability_'+type+'.log'
if file_test(ifn) eq 0 then message, 'file not found ...'

ofn = shomedir()+'/psbl_de_32hz/asi_list_'+type+'.log'
if file_test(file_dirname(ofn),/directory) eq 0 then file_mkdir, ofn
if file_test(ofn) eq 1 then file_delete, ofn


nline = file_lines(ifn)
lines = strarr(nline)
openr, lun, ifn, /get_lun
readf, lun, lines
free_lun, lun


for i = 0, nline-1 do begin
    tline = lines[i]
    if tline eq '' then continue
    if strmid(tline,0,4) eq '****' then begin
        idline = tline
        
        i = i+3     ; jump 3 lines to read time.
        asiut = time_double(lines[i])
        i = i+1
        tline = lines[i]

        sites = []
        
        while tline ne '' do begin
            i = i+1
            if i ge nline then break
            tline = lines[i]
            sites = [sites,tline]
        endwhile
        
        idx = stregex(sites, 'file exists')
        idx = where(idx ne -1, cnt)
        if cnt eq 0 then continue   ; no site has file.
        sites = sites[idx]
        
        idx = stregex(sites, 'no moon')
        idx = where(idx ne -1, cnt)
        if cnt eq 0 then continue   ; all sites have moon.
        sites = sites[idx]
        
        nsite = n_elements(sites)
        for j = 0, nsite-1 do sites[j] = strmid(sites[j],0,4)
        if nsite gt 1 then sites = strjoin(sites,',') 
       

        openw, lun, ofn, /get_lun, /append
        printf, lun, ''
        printf, lun, idline
        printf, lun, time_string(asiut)
        printf, lun, sites
        free_lun, lun
    endif
endfor

end
