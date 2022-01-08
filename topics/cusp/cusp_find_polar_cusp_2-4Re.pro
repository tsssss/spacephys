
; pro cusp_find_polar_cusp_2-3Re

    ; set months to read.
    yrs = [1996,1997,1998,1999,2000]
    mos = [9,10,11]
    
    ; range in distance and mlt.
    mindis = 2
    maxdis = 4
    minmlt = 10
    maxmlt = 14
    
    infos = []
    logs = []
 
    ; pattern for cusp_list_polar_yyyy_mm.log.
    rootdir = shomedir()+'/Google Drive/works/works/cusp/cusp list polar'
    baseptn = 'cusp_list_polar_YYYY_MM.log'

    ; read cusp_list_polar_yyyy_mm.log.
    foreach tyr, yrs do begin
        foreach tmo, mos do begin
            ; prepare filename.
            tut = string(tyr,format='(I4)')+string(tmo,format='(I02)')
            tut = sfmdate(tut,'%Y%m')
            tfn = sprepfile(tut, paths = [rootdir,baseptn])
            
            ; read content.
            if file_test(tfn) eq 0 then continue
            nline = file_lines(tfn)
            lines = strarr(nline)
            openr, lun, tfn, /get_lun
            readf, lun, lines
            free_lun, lun
            
            ; parse lines.
            header = lines[0:3]
            for i = 4, nline-1 do begin ; skip 4 lines of header.
                tline = lines[i]
                tmp = strmid(tline,0,13)
                ut = sfmdate(tmp,'%Y-%m%d  %H')
                id = stodate(ut, '%Y_%m%d_%H')
                
                print, 'processing '+id+'...'
                
                tmp = strmid(tline,24,7)
                tmp = strsplit(tmp,'-',/extract)
                diss = float(strtrim(tmp,2))
                if min(diss) gt maxdis then continue
                if max(diss) lt mindis then continue
                
                tmp = strmid(tline,45,10)
                tmp = strsplit(tmp,'-',/extract)
                mlts = float(strtrim(tmp,2))
                if min(mlts) gt maxmlt then continue
                if max(mlts) lt minmlt then continue
                
                tmp = strmid(tline,33,10)
                hem = (strmid(tmp,0,1) eq 'N')? 1: -1
                tmp = strsplit(strmid(tmp,1),'-',/extract)
                ilats = float(strtrim(tmp,2))*hem
                
                dst = float(strtrim(strmid(tline,56,4),2))
                ae = float(strtrim(strmid(tline,62,4),2))
                
                tinfo = {id:id, diss:diss, mlts:mlts, ilats:ilats, $
                    dst:dst, ae:ae}
                infos = [infos,tinfo]
                logs = [logs,tline]
            endfor
        endforeach
    endforeach
    
    ; save to tplot.
    store_data, 'polar_2-4Re', 0, infos
    
    ofn = shomedir()+'/cusp_list_polar_2-4Re.tplot'
    tplot_save, 'polar_2-4Re', filename = ofn
    
    ofn = shomedir()+'/cusp_list_polar_2-4Re.log'
    openw, lun, ofn, /get_lun
    printf, lun, header
    printf, lun, logs
    free_lun, lun
    
end