;+
; generate a structure for all large E field events, 
; include those without hope data, exclude bad E data, wake effect events.
;-

pro psbl_de_all_e_gen_id_list

    rootdir = shomedir()+$
        '/Google Drive/works/works/rbsp_de'
    probes = ['a','b']
    nprobe = n_elements(probes)

    logfns = [rootdir+'/psbl_de_round3/list_rbsp'+probes]+ $
        '_large_efield_32hz_quick.log'
    
    

    infos = []
    info0 = {id:'', name:'psbl_de_all_e', probe:'', utr:[0d,0d], $
        posgsm:[0d,0,0], emax:0d, mlt:0d, lshell:0d, mlat:0d}

    for i = 0, nprobe-1 do begin
        tprb = probes[i]
        
        nline = file_lines(logfns[i])
        lines = strarr(nline)
        openr, lun, logfns[i], /get_lun
        readf, lun, lines
        free_lun, lun
        lines = lines[3:*]

        ; remove wake events.
        idx = stregex(lines, 'wake')
        idx = where(idx eq -1, nline)
        lines = lines[idx]
        
        for j = 0, nline-1 do begin
            tline = lines[j]
            tid = strmid(tline,0,14)+'_'+tprb
            tdate = strmid(tline,0,9)
            tutr = time_double(tdate+[strmid(tline,18,8),strmid(tline,30,8)], $
                tformat='YYYY_MMDDhh:mm:ss')
            emax = double(strmid(tline,42,6))
            tpos = double([strmid(tline,52,4),strmid(tline,60,4),strmid(tline,68,4)])
            tnote = strmid(tline,76)
            
            infos = [infos, info0]
            infos[-1].id = tid
            infos[-1].utr = tutr
            infos[-1].probe = tprb
            infos[-1].posgsm = tpos
            infos[-1].emax = emax
            
            efw = sread_rbsp_efw_l3(tutr, probes = tprb, $
                vars = ['epoch','mlt_lshell_mlat'])
            if size(efw,/type) eq 8 then begin
                uts = sfmepoch(efw.epoch,'unix')
                tut = mean(tutr)
                pos = sinterpol(efw.mlt_lshell_mlat, uts, tut)
                infos[-1].mlt = pos[0]
                infos[-1].lshell = pos[1]
                infos[-1].mlat = pos[2]
            endif else begin
                infos[-1].mlt = !values.d_nan
                infos[-1].lshell = !values.d_nan
                infos[-1].mlat = !values.d_nan
            endelse
        endfor
    endfor
    
    ofn = rootdir+'/eventids/psbl_de_ids_all_e.dat'
    ids_all_e = infos
    save, ids_all_e, filename = ofn


end
