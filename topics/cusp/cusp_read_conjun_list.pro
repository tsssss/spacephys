;+
; Purpose: Read the event info from formatted conjunction list.
; may contain empty values for certain column.
; set nofast for polar only events.
;-

function cusp_parse_conjun_list, str0, colmask, nofast = nofast
    str = str0
    if n_elements(colmask) eq 0 then begin
        info0 = strsplit(str,' ',/extract)
    endif else begin
        idx = strsplit(colmask,' ',length = lens)
        ncol = n_elements(lens)
        info0 = strarr(ncol)
        ; no need to worry about idx[i] is longer than strlen(str).
        for i = 0, ncol-1 do info0[i] = strmid(str,idx[i],lens[i])
;        foreach tmp, info0 do print, '"'+tmp+'"'
    endelse
    id = info0[0]
    date = strmid(id,0,10)
    ; polar.
    poidx = 1
    dedb = double(strsplit(info0[poidx+8],',',/extract))
    if n_elements(dedb) ne 2 then dedb = [0d,0d]
    keinfo = strsplit(info0[poidx+2],',',/extract)
    polar = {$
        plot_time: time_double(date+strsplit(info0[poidx+0],',',/extract), $
            tformat='YYYY_MMDD_hh:mm'), $
        cusp_time: time_double(date+strsplit(info0[poidx+1],',',/extract), $
            tformat='YYYY_MMDD_hh:mm'), $
        field_time: time_double(date+strsplit(info0[poidx+5],',',/extract), $
            tformat='YYYY_MMDD_hh:mm'), $
        keflux: double(keinfo[0:1]), $          ; [ele,ion] max value.
        ketype: keinfo[2], $    ; etp,p2p,ptp for epstopdf, ps2pdf, pstopdf.
        nflux: double(strsplit(info0[poidx+3],',',/extract)), $
        n: double(strsplit(info0[poidx+4],',',/extract)), $
        e56: info0[poidx+6] eq 'ok', $
        scaleinfo: double(strsplit(info0[poidx+7],',',/extract)), $
        de: [-1d,1]*dedb[0], $
        db: [-1d,1]*dedb[1], $
        filters: double(strsplit(info0[poidx+9],',',/extract)), $
        faclim: double(info0[poidx+10]), $
        noisedelim: double(strsplit(info0[poidx+11],',',/extract))}

    if keyword_set(nofast) then return, {id:id, polar:polar}

    ; fast.
    faidx = 13
    dedb = double(strsplit(info0[faidx+4],',',/extract))
    if n_elements(dedb) ne 2 then dedb = [0d,0d]
    fast = {$
        plot_time: time_double(date+strsplit(info0[faidx+0],',',/extract), $
            tformat='YYYY_MMDD_hh:mm'), $
        cusp_time: time_double(date+strsplit(info0[faidx+1],',',/extract), $
            tformat='YYYY_MMDD_hh:mm:ss'), $
        orbit: long(info0[faidx+2]), $
        scaleinfo: double(strsplit(info0[faidx+3],',',/extract)), $
        de: [-1d,1]*dedb[0], $
        db: [-1d,1]*dedb[1], $
        filters: double(strsplit(info0[faidx+5],',',/extract)), $
        faclim: double(info0[faidx+6]), $
        noisedelim: double(strsplit(info0[faidx+7],',',/extract)), $
        ionrange: double(strsplit(info0[faidx+8],',',/extract))}
    tmp = strmid(id,10,2)
    if fix(tmp) ge 24 then dt = 86400d else dt = 0d
    fast.cusp_time+= dt
    fast.plot_time+= dt
    ; combine.
    return, {id:id, polar:polar, fast:fast}
end

function cusp_read_conjun_list, fn, event = id, nofast = nofast
    compile_opt idl2
    if ~file_test(fn) then message, 'file does not exist ...'
    header = strarr(4)
    
    ; read all event info.
    nevent = (file_lines(fn)-4)[0]  ; 4 lines of header.
    if nevent le 0 then message, 'no event ...'
    events = strarr(nevent)
    openr, lun, fn, /get_lun
    readf, lun, header
    readf, lun, events
    free_lun, lun
    
    colmask = header[3]    ; column mask.
    eventids = strarr(nevent)
    for i = 0, nevent-1 do eventids[i] = strmid(events[i],0,12)
    
    idx = where(id[0] eq eventids, cnt)
    if cnt eq 0 then return, -1
    return, cusp_parse_conjun_list(events[idx],colmask, nofast = nofast)
end

fn = shomedir()+'/Google Drive/works/works/cusp/cusp_list_of_conjun_9_10_all.log'
eventid = '1998_0918_13'
eventid = '1998_0925_05'
eventid = '1998_0922_23'
infos = cusp_read_conjun_list(fn, event = eventid)
end