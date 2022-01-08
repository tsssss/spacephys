;+
; generate a structure var by reading the event list log files.
;-

pro psbl_de_detect_psbl_gen_id_list

    rootdir = shomedir()+$
        '/Google Drive/works/works/rbsp_de'
    probes = ['a','b']
    nprobe = n_elements(probes)

    logfns = rootdir+'/psbl_de_detect_step/psbl_de_detect_psbl_info_'+ $
        probes+'.log'
    
    infos = []
    info0 = {id:'', name:'psbl_de_detect_step', probe:'', utr:[0d,0d], $
        dn_jump:-1d, dt_jump:-1d}   ; -1 means no jump, otherwise step value.

    for i = 0, nprobe-1 do begin
        nline = file_lines(logfns[i])
        lines = strarr(nline)
        openr, lun, logfns[i], /get_lun
        readf, lun, lines
        free_lun, lun
        
        for j = 0, nline-1 do begin
            if lines[j] eq '' then continue
            if strmid(lines[j], 0,1) eq '*' then begin
                thead = strmid(lines[j], 1)
                tmp = strpos(thead, '.')
                if tmp ne -1 then thead = strmid(thead, 0, tmp)
                heads = strsplit(thead, ',', /extract)
                nhead = n_elements(heads)
                for k = 0, nhead-1 do begin
                    thead = strtrim(strsplit(heads[k], ':', /extract), 2)
                    nentry = fix(strmid(thead[0],0,strpos(thead[0],' ')))
                    step_type = strlowcase(thead[1])    ; dn, dT, or dn & dT.
                    for l = 0, nentry-1 do begin
                        j = j+1
                        tmp = strsplit(lines[j],' ',/extract)
                        tprb = strmid(tmp[2],4,1)
                        tutr = time_double(tmp[0:1])
                        tval = float(tmp[3])

                        infos = [infos, info0]
                        
                        infos[-1].id = time_string(tutr[0], $
                            tformat='YYYY_MMDD_hhmm_'+tprb)
                        infos[-1].utr = tutr
                        infos[-1].probe = tprb
                        if strpos(step_type, 'dn') ne -1 then $
                            infos[-1].dn_jump = 0
                        if strpos(step_type, 'dt') ne -1 then $
                            infos[-1].dt_jump = 0
                        if strmid(tmp[2],6,1) eq 'n' then $
                            infos[-1].dn_jump = tval $
                        else infos[-1].dt_jump = tval
                    endfor
                endfor
            endif
        endfor
    endfor
    
    ofn = rootdir+'/eventids/psbl_de_ids_detect_step.dat'
    ids_detect_step = infos
    save, ids_detect_step, filename = ofn

end
