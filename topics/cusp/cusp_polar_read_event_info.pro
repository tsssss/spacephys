;+
; Purpose: Read the event info from formatted list.
; may contain empty values for certain column.
;-

function cusp_polar_parse_event_list_parse_time, date, str
    the_time = strsplit(str, ',', /extract)
    if n_elements(the_time) eq 1 then the_time = ['','']
    time_range = time_double(date+the_time, tformat='YYYY_MMDD_hh:mm')
    if time_range[1] lt time_range[0] then time_range[1] += constant('secofday')
    return, time_range
end

function cusp_polar_parse_event_list, str0, colmask

    secofday = constant('secofday')

    str = str0
    if n_elements(colmask) eq 0 then begin
        info0 = strsplit(str,' ',/extract)
    endif else begin
        idx = strsplit(colmask,' ',length=lens)
        ncol = n_elements(lens)
        info0 = strarr(ncol)
        ; no need to worry about idx[i] is longer than strlen(str).
        for i = 0, ncol-1 do info0[i] = strmid(str,idx[i],lens[i])
;        foreach tmp, info0 do print, '"'+tmp+'"'
    endelse
    id = info0[0]
    date = strmid(id,0,10)

;---Polar.
    poidx = 1

    plot_time = cusp_polar_parse_event_list_parse_time(date,info0[poidx+0])
    cusp_time = cusp_polar_parse_event_list_parse_time(date,info0[poidx+1])
    field_time = cusp_polar_parse_event_list_parse_time(date,info0[poidx+5])

    dedb = double(strsplit(info0[poidx+8],',',/extract))
    if n_elements(dedb) ne 2 then dedb = [0d,0d]
    keinfo = strsplit(info0[poidx+2],',',/extract)

    polar = {$
        plot_time: plot_time, $
        cusp_time: cusp_time, $
        field_time: field_time, $
        keflux: double(keinfo[0:1]), $          ; [ele,ion] max value.
        ketype: strtrim(keinfo[2],2), $    ; etp,p2p,ptp for epstopdf, ps2pdf, pstopdf.
        nflux: double(strsplit(info0[poidx+3],',',/extract)), $
        n: double(strsplit(info0[poidx+4],',',/extract)), $
        e56: info0[poidx+6] eq 'ok', $
        scaleinfo: double(strsplit(info0[poidx+7],',',/extract)), $
        de: [-1d,1]*dedb[0], $
        db: [-1d,1]*dedb[1], $
        filters: double(strsplit(info0[poidx+9],',',/extract)), $
        faclim: double(info0[poidx+10]), $
        noisedelim: double(strsplit(info0[poidx+11],',',/extract))}

    return, {id:id, polar:polar}

end

function cusp_polar_read_event_info, fn, event=id

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
    return, cusp_polar_parse_event_list(events[idx],colmask)
end

fn = googledir()+'/works/works/cusp/cusp_list_of_polar_2-4Re.log'
eventid = '2000_0914_08'
infos = cusp_polar_read_event_info(fn, event=eventid)
end