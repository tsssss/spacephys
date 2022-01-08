idinfos = psbl_de_id('detect_step')
idinfos = psbl_de_id('all_e')
nidinfo = n_elements(idinfos)

rootdir = shomedir()+'/psbl_de_32hz/'
if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
logfn = rootdir+'/rbsp_burst_availability.log'

dt0 = -600 ; min overlap time b/w given time and burst data.

for j = 0, nidinfo-1 do begin
    id = idinfos[j].id
    utr = idinfos[j].utr
    tprobe = strmid(id,0,1,/reverse)
    print, 'processing '+id+' ...'
    openw, loglun, logfn, /get_lun, /append
    printf, loglun, ''
    printf, loglun, '**** processing '+id+' ...'
    printf, loglun, ''
    free_lun, loglun

    remroot = 'http://rbsp.space.umn.edu/data/rbsp/burst_playback'
    locroot = spreproot('rbsp')
    prb = tprobe
    vsn = 'v[0-9.]{2}'
    ext = 'txt'
    types = ['vb1','vb2']
    
    foreach type, types do begin
        baseptn = 'rbsp'+prb+'_efw_'+type+'_playback_YYYYMMDD_'+vsn+'.'+ext
        rempaths = [remroot,'rbsp'+prb,type+'_playback','YYYY',baseptn]
        locpaths = [locroot,'rbsp'+prb,'efw/burst_playback',type,'YYYY',baseptn]
        
        remfns = sprepfile(utr, paths = rempaths)
        locfns = sprepfile(utr, paths = locpaths)
        nfn = n_elements(locfns)
        for i = 0, nfn-1 do begin
            basefn = file_basename(locfns[i])
            locpath = file_dirname(locfns[i])
            rempath = file_dirname(remfns[i])
            locfns[i] = sgetfile(basefn, locpath, rempath)
        endfor
        idx = where(locfns ne '', nfn)
        if nfn eq 0 then begin
            openw, loglun, logfn, /get_lun, /append
            printf, loglun, 'no '+type+' data ...'
            free_lun, loglun
        endif else begin
            butrs = []
            locfns = locfns[idx]
            for i = 0, nfn-1 do begin
                nline = file_lines(locfns[i])
                lines = strarr(nline)
                openr, lun, locfns[i], /get_lun
                readf, lun, lines
                free_lun, lun
                butrs = [butrs,lines[4:*]]  ; exclude header.
            endfor
            
            nbutr = n_elements(butrs)
            lines = butrs
            flags = bytarr(nbutr)+1
            
            for i = 0, nbutr-1 do begin
                tbutr = time_double(strsplit(lines[i],' ',/extract))
                if max(tbutr)-min(utr) le dt0 then flags[i] = 0
                if min(tbutr)-max(utr) ge dt0 then flags[i] = 0
            endfor
            
            idx = where(flags eq 1, cnt)
            if cnt eq 0 then begin
                openw, loglun, logfn, /get_lun, /append
                printf, loglun, type+': no overlap within -/+'+ $
                    snum2str(abs(dt0))+' sec ...'
                free_lun, loglun
            endif else begin
                lines = lines[idx]
                openw, loglun, logfn, /get_lun, /append
                foreach tline, lines do printf, loglun, type+': '+tline
                free_lun, loglun
            endelse
        endelse
    endforeach
endfor


end
