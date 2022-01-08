
    server = 'themis.ssl.berkeley.edu'
    port = 80
    url = 'http://themis.ssl.berkeley.edu/data/themis/thg/l2/asi/cal/thm_map_add.sav'
    purl = '/data/themis/thg/l2/asi/cal/thm_map_add.sav'
    socket, unit, server, port, /get_lun, /swap_if_little_endian, error = error, $
        read_timeout = read_timeout, connect_timeout = connect_timeout

    printf, unit, 'GET '+purl +  ' HTTP/1.0'
    printf, unit, 'Host: ' + server
    printf, unit, ''

    ; read header.
    header = strarr(256)
    text = 'xxx'
    linesread = 0
    WHILE text NE '' do begin
        readf, unit, text
        Header[LinesRead] = text
        LinesRead = LinesRead+1
        IF LinesRead MOD 256 EQ 0 THEN $
            Header=[Header, StrArr(256)]
    ENDWHILE
    
    if LinesRead eq 0 then begin
        free_lun, unit
        url_info.io_error = 1
        dprint,dlevel=0,verbose=verbose,!error_state.msg
    endif

    Header = Header[0:LinesRead-1]
    
    localname = 'M:/data/themis/thg/l2/asi/cal/thm_map_add.sav'
    openw, wunit, localname, /get_lun
    free_lun, wunit
    
end