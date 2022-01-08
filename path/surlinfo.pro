;+
; get info for given remote file, including server, path, port, 
; file size, last-modified time.
; 
; return -1 if the given remote file doesn't exist.
; 
;-
function surlinfo, remfn0

;    isunix = !version.os_family eq 'unix'
    remfn = strimpath(remfn0,/keep_trailing_slash)   ; for directory, include trainling /.
    
    ; server.
    idx = strpos(remfn,'://')+3
    server = strmid(remfn, idx)
    server = strmid(server, 0, strpos(server,'/'))
    
    ; path.
    purl = strmid(remfn,strlen(server)+idx)
    
    ; port.
    port = 80
    
    socket, remlun, server, port, error = err, /get_lun, connect_timeout = 10
    
    if err ne 0 then begin       ; cannot open server.
        if n_elements(remlun) ne 0 then free_lun, remlun
        return, -1
    endif
    
    ; send message.
    printf, remlun, 'GET '+purl+' HTTP/2.0'
    printf, remlun, 'Host: '+server
    printf, remlun, ''
    
    ; read header.
    header = ''
    tline = 'xxx'
    while tline ne '' do begin
        readf, remlun, tline
        header = [header,tline]
    endwhile
    header = header[1:*]
    
    idx = where(stregex(header[0],'404') ne -1, cnt)
    if cnt ne 0 then begin
        free_lun, remlun
        return, -1     ; file not found.
    endif
    
    ; filesize.
    idx = where(stregex(header,'Content-Length:') ne -1, cnt)
    fsize = cnt? ulong64(strmid(header[idx],strpos(header[idx],':')+1)): 0ull
    
    ; last modified time, in universal time.
    idx = where(stregex(header,'Last-Modified') ne -1, cnt)
    mtime0 = cnt? (strtrim(strmid(header[idx],strpos(header[idx],':')+1),2))[0]: ''
    
    free_lun, remlun
    
    ut0 = (mtime0 eq '')? 0d: sfmdate(mtime0,'%a, %d %b %Y %H:%M:%S %Z')
    info = {server:server, path:purl, port:port, mtime:ut0, size:fsize[0]}
    
    return, info
    
end

url = 'http://themis.ssl.berkeley.edu/data/rbsp/rbspb/l1/vb1/2015/'
;url = 'http://themis.ssl.berkeley.edu/data/rbsp/rbspb/l1/vb1/2015/rbspb_l1_vb1_20151219_v02.cdf'
info = surlinfo(url)
end
