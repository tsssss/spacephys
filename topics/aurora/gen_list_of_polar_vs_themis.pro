
pro gen_list_of_polar_vs_themis

    spc = '  '
    console = -1
    
    rootdir = shomedir()+'/polar_uvi_vs_themis_asi'
    ofn = shomedir()+'/list_polar_vs_themis_aurora.log'
    sites = ['atha','chbg','ekat','fsmi','fsim','fykn',$
        'gako','gbay','gill','inuv','kapu','kian',$
        'kuuj','mcgr','nrsq','pgeo','pina','rank',$
        'snap','snkq','talo','tpas','whit','yknf',$
        'resu']
    nsite = n_elements(sites)
    thmdir = 'http://themis.ssl.berkeley.edu/data/themis/thg/l1/asi'
    
    if file_test(rootdir,/directory) eq 0 then message, 'no such directory ...'
    tmp = file_search(rootdir+'/*')
    dirs = ['']
    for i = 0, n_elements(tmp)-1 do begin
        tdirs = file_search(tmp[i]+'/po_uvi_*', count = cnt)
        if cnt eq 0 then continue   ; wrong directory.
        dirs = [dirs,tdirs]
    endfor
    dirs = dirs[where(dirs ne '', ndir)]
    
    if ndir eq 0 then return
    openw, lun, ofn, /get_lun
    printf, lun, 'YYYY_MMDD_HH'+spc+'Symh'+spc+' AE '+spc+'ASI sites'
    free_lun, lun
    
    for i = 0, ndir-1 do begin
        tdir = dirs[i]
        uviflags = bytarr(24)
        
        yr = strmid(tdir,51,4)
        mo = strmid(tdir,56,2)
        dy = strmid(tdir,58,2)
        
        ; test polar availability.
        tfns = file_search(tdir+'/po_uvi_mlt_*.png', count = cnt)
        for j = 0, cnt-1 do begin
            tmp = file_basename(tfns[j])
            hr = strmid(tmp,21,2)
            uviflags[fix(hr)] = 1b
        endfor
        
        idx = where(uviflags eq 1b, cnt)
        if cnt eq 0 then continue
        
        ; test themis asf.
        for j = 0, cnt-1 do begin
            hr = string(idx[j],format='(I02)')
            asfflags = bytarr(nsite)
            for k = 0, nsite-1 do begin
                tfn = thmdir+'/'+sites[k]+'/'+yr+'/'+mo+'/thg_l1_asf_'+$
                    sites[k]+'_'+yr+mo+dy+hr+'_v01.cdf'
;                file_http_copy, tfn, url_info = info, /no_download
                scurl, tfn, info = info
                if info.size gt 0 then asfflags[k] = 1b
            endfor
            if total(asfflags) ne 0 then begin
                ; read ae and dst.
                tmp = strjoin([yr,mo,dy],'-')+'/'+hr+':00'
                tet = stoepoch(tmp)
                load = 0
                if size(omni,/type) ne 8 then load = 1 else begin
                    if min(omni.epoch) gt tet then load = 1
                    if max(omni.epoch) lt tet then load = 1
                endelse
                tmp = strjoin([yr,mo,dy],'-')
                if load then omni = sread_omni(tmp)
                tmp = min(tet-omni.epoch, rec1, /absolute)
                rec2 = (rec1+59) < (n_elements(omni.epoch)-1)
                symh = min(omni.symh[rec1:rec2])
                ae = max(omni.ae[rec1:rec2])
                
                ; construct command.
                cmd = ''
                cmd+= yr+'_'+mo+dy+'_'+hr +spc
                cmd+= string(symh,format='(I4)')+spc
                cmd+= string(ae,format='(I4)')+spc
                cmd+= strjoin(sites[where(asfflags eq 1b)],',')
                printf, console, cmd
                openw, lun, ofn, /get_lun, /append
                printf, lun, cmd
                free_lun, lun
            endif
        endfor
    endfor
end