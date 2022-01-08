
pro cusp_gen_excel_form, id, filename = ofn, load = load, test = test

    rootdir = shomedir()+'/Google Drive/works/data/cusp'
    
    if keyword_set(test) then rootdir = shomedir()+'/cusp/data'
    
    lun = -1
    if n_elements(ofn) ne 0 then openw, lun, ofn, /get_lun
    
    fns = rootdir+'/*_all_data.tplot'
    if n_elements(id) ne 0 then fns = rootdir+'/'+id+'_all_data.tplot'
    fns = file_search(fns, count = nfn)
    if nfn eq 0 then message, 'no data ...'
    
    get_data, 'cusp_stats', tmp
    if n_elements(tmp) eq 0 then load = 1

    tinfo = cusp_info_struct()
    infos = replicate(tinfo, nfn)

    if n_elements(load) eq 0 then load = 0
    
    if load then begin
        for i = 0, nfn-1 do begin
            print, fns[i]
            tplot_restore, filename = fns[i]
            get_data, 'scidat', tmp, info
            infos[i] = info
        endfor
        store_data, '*', /delete
        store_data, 'cusp_stats', 0, infos
    endif

    get_data, 'cusp_stats', tmp, infos
    ninfo = n_elements(infos)
    idx = sort(infos.id)
    infos = infos[idx]

    ; if output to console, use '  ', else use tab.
    sep = (lun eq -1)? '  ': string(9b)
    sep2 = (lun eq -1)? ',': string(9b)
    sep3 = (lun eq -1)? ' ': string(9b)


    for i = 0, ninfo-1 do begin

        cmd = ''
        tit = ''

        tinfo = infos[i]

        ; event id.
        tit+= 'yyyy_mmdd_hh'+sep
        cmd+= tinfo.id+sep

        ; ae.
        tit+= ' AE '+sep
        cmd+= string(tinfo.ae, format='(I4)')+sep

        ; symh.
        tit+= 'Symh'+sep
        cmd+= string(tinfo.dst, format='(I4)')+sep

        ; dt.
        dt = (tinfo.fast.cusp.entry.ut-tinfo.polar.cusp.entry.ut)/3600d
        tit+= ' dt '+sep
        cmd+= string(dt, format='(F4.1)')+sep

        ; dmlt.
        tit+= 'dMLT'+sep
        cmd+= string(tinfo.dmlt, format='(F4.1)')+sep

        ; imf bz, by.
        tit+= ' Bz'+sep
        cmd+= string(round(max(tinfo.imfbz,/absolute)), format='(I3)')+sep
        tit+= ' By'+sep
        cmd+= string(round(max(tinfo.imfby,/absolute)), format='(I3)')+sep

        ; max upward kei.
        tit+= 'KEi'+sep
        tmp = sgnum2str(tinfo.max_kei, msgn=2)
        cmd+= string(tmp,format='(A3)')+sep
        
        ; ion ratio.
        tit+= ' ion'+sep
        cmd+= string(tinfo.ratio_ion, format='(F4.1)')+sep

        ; poynting flux ratio.
        tit+= ' Sb '+sep
        cmd+= string(tinfo.ratio_pflux, format='(F4.1)')+sep

        ; energy flux ratio.
        tit+= 'consv'+sep
        cmd+= string(tinfo.ratio_eflux, format='(F5.2)')+sep
        
        ; polar poynting flux ratio.
        tit+= 'poynt'+sep
        tmp = (tinfo.fast.eflux-tinfo.polar.kei-tinfo.polar.kee)/tinfo.polar.pflux
        cmd+= string(tmp, format='(F5.2)')+sep


    ; **** cusp info.
        tit+= ' R1'+sep
        cmd+= string(tinfo.polar.cusp.entry.dis, format='(F3.1)')+sep
        tit+= ' mlt'+sep
        cmd+= string(tinfo.polar.cusp.entry.mlt, format='(F4.1)')+sep
        tit+= ' dR'+sep
        cmd+= string(tinfo.dr, format='(F3.1)')+sep


    ; **** polar eflux.
        tmp = tinfo.polar
        
        tit+= ' KEi'+sep2+' KEe'+sep3
        cmd+= string(tmp.kei, format='(I4)')+sep2+string(tmp.kee, format='(I4)')+sep3
        
        tit+= '  S  '+sep2
        cmd+= string(round(tmp.pflux), format='(I5)')+sep2
        tit+= '  Sb '+sep2
        cmd+= string(round(tmp.sb.fh), format='(I5)')+sep2
        tit+= ' SvSp'+sep2
        cmd+= string(round(tmp.sv.fh+tmp.sp.fh), format='(I5)')+sep2
        tit+= 'SbFac'+sep2
        cmd+= string(round(tmp.sb.fl), format='(I5)')+sep2
        tit+= 'SbWav'+sep2
        cmd+= string(round(tmp.sb.fh), format='(I5)')+sep2


    ; **** fast eflux.
        tmp = tinfo.fast
        
        tit+= ' KEi'+sep2+' KEe'+sep3
        cmd+= string(tmp.kei, format='(I4)')+sep2+string(tmp.kee, format='(I4)')+sep3
        
        tit+= '  S  '+sep2
        cmd+= string(round(tmp.pflux), format='(I5)')+sep2
        tit+= '  Sb '+sep2
        cmd+= string(round(tmp.sb.fh+tmp.sb.fl), format='(I5)')+sep2
        tit+= ' SvSp'+sep2
        cmd+= string(round(tmp.sv.fh+tmp.sv.fl+tmp.sp.fh+tmp.sp.fl), format='(I5)')+sep2
        tit+= 'SbFac'+sep2
        cmd+= string(round(tmp.sb.fl), format='(I5)')+sep2
        tit+= 'SbWav'+sep2
        cmd+= string(round(tmp.sb.fh), format='(I5)')+sep2

        cmd = strmid(cmd,0,strlen(cmd)-strlen(sep2))
        printf, lun, cmd
    endfor

    tit = strmid(tit,0,strlen(tit)-strlen(sep2))
    printf, lun, tit
    
    if n_elements(ofn) ne 0 then free_lun, lun
end

ids = cusp_id('low_altitude')
cusp_gen_excel_form, ids, /load
end
