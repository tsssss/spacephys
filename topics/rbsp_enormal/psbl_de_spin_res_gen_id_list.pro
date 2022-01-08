;+
; generate a structure var by reading the event list log files.
;-

pro psbl_de_spin_res_gen_id_list

    rootdir = shomedir()+$
        '/Google Drive/works/works/rbsp_de'

    logfns = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
        'dipolarization/list_large_de_round3.log'
    
    infos = []
    info0 = {id:'', name:'psbl_de_spin_res', probe:'', utr:[0d,0d], $
        dn_jump:-1d, dt_jump:-1d}   ; -1 means no jump, otherwise step value.

    nline = file_lines(logfns)
    nhead = 3
    nline = nline-nhead
    lines = strarr(nline)
    heads = strarr(nhead)
    openr, lun, logfns, /get_lun
    readf, lun, heads
    readf, lun, lines
    free_lun, lun
    
    for j = 0, nline-1 do begin
        tline = lines[j]
        if tline eq '' then continue
        
        id = strmid(tline,0,16)
        tprobe = strmid(id,strlen(id)-1)
        
        t1 = time_double(strmid(id,0,9)+strmid(tline,20,5),tformat='YYYY_MMDDhh:mm')
        t2 = time_double(strmid(id,0,9)+strmid(tline,29,5),tformat='YYYY_MMDDhh:mm')
        if t2 lt t1 then t2+= 86400d
        tutr = [t1,t2]
        
        infos = [infos, info0]
        
        infos[-1].id = id
        infos[-1].utr = tutr
        infos[-1].probe = tprobe
    endfor
    
    ofn = rootdir+'/eventids/psbl_de_ids_spin_res.dat'
    ids_spin_res = infos
    save, ids_spin_res, filename = ofn

end
