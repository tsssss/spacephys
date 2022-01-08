
type = 'detect_step'

ifn = shomedir()+'/psbl_de_32hz/rbsp_burst_availability_'+type+'.log'
if file_test(ifn) eq 0 then message, 'file not found ...'

ofn = shomedir()+'/psbl_de_32hz/rbsp_burst_list_'+type+'.log'
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
        vb1s = []
        vb2s = []
        i = i+1
        while tline ne '' do begin
            i = i+1
            if i ge nline then break
            tline = lines[i]
            if strpos(tline, 'vb1') ne -1 then vb1s = [vb1s,tline]
            if strpos(tline, 'vb2') ne -1 then vb2s = [vb2s,tline]
        endwhile

        idx = stregex(vb1s, 'no')
        idx = where(idx eq -1, cnt)
        if cnt eq 0 then hasvb1 = 0 else begin
            vb1s = vb1s[idx]
            hasvb1 = 1
        endelse

        idx = stregex(vb2s, 'no')
        idx = where(idx eq -1, cnt)
        if cnt eq 0 then hasvb2 = 0 else begin
            vb2s = vb2s[idx]
            hasvb2 = 1
        endelse
        
        maxdt0 = 60     ; 1 min.
        tut = time_double(strmid(idline,16,14), tformat='YYYY_MMDD_hhmm')

        if hasvb1 then begin
            nvb1 = n_elements(vb1s)
            flags = bytarr(nvb1)+1
            for j = 0, nvb1-1 do begin
                tmp = strmid(vb1s[j],5)
                tmp = strsplit(tmp,' ',/extract)
                tutr = time_double(tmp,tformat='YYYY-MM-DD/hh:mm:ss.ffffff')
                if tut-tutr[1] gt maxdt0 then flags[j] = 0
                if tutr[0]-tut gt maxdt0 then flags[j] = 0
            endfor
            idx = where(flags eq 1, cnt)
            if cnt eq 0 then hasv12 = 0 else vb1s = vb1s[idx]
        endif
        if hasvb1 then vb1s = vb1s[0]   ; only keep one, enough for vb1.
        
        if hasvb2 then begin
            nvb2 = n_elements(vb2s)
            flags = bytarr(nvb2)+1
            for j = 0, nvb2-1 do begin
                tmp = strmid(vb2s[j],5)
                tmp = strsplit(tmp,' ',/extract)
                tutr = time_double(tmp,tformat='YYYY-MM-DD/hh:mm:ss.ffffff')
                if tut-tutr[1] gt maxdt0 then flags[j] = 0
                if tutr[0]-tut gt maxdt0 then flags[j] = 0
            endfor
            idx = where(flags eq 1, cnt)
            if cnt eq 0 then hasvb2 = 0 else vb2s = vb2s[idx]
        endif


        if hasvb1 eq 0 and hasvb2 eq 0 then continue
        openw, lun, ofn, /get_lun, /append
        printf, lun, ''
        printf, lun, idline
        if hasvb1 then printf, lun, vb1s
        if hasvb2 then printf, lun, vb2s
        free_lun, lun
    endif
endfor

end
