
pro cusp_check_files, locroot

    if strmid(locroot,strlen(locroot)-1) eq '/' then $      ; remove trailing '/'.
        locroot = strmid(locroot,0,strlen(locroot)-1)
    fns = file_search(locroot+'/*')
    for i = 0, n_elements(fns)-1 do begin
        tfn = fns[i]
        if ~file_test(tfn, /directory) then continue
        id = strmid(tfn, strlen(locroot)+1)
        tid = strmid(id,0,9)
        tfns = file_search(tfn+'/*')
        if tfns[0] eq '' then continue      ; no files.
        ; convert yyyy_mmdd_xxxx to yyyy_mmdd_hh_xxxx.
        for j = 0, n_elements(tfns)-1 do begin
            if file_test(tfns[j], /directory) then continue
            ttfn = file_basename(tfns[j])
            ttfns = strsplit(ttfn,'_',/extract)
            if strjoin(ttfns[0:2],'_') eq id then begin     ; already in yyyy_mmdd_hh.
                continue
            endif else if ttfns[2] eq 'polar' then begin    ; yyyy_mmdd_polar_xxxx.
                continue
            endif else if strjoin(ttfns[0:1],'_') eq tid then begin
                spawn, 'mv "'+tfn+'/'+ttfn+'" "'+tfn+'/'+strjoin([id,ttfns[2:*]],'_')+'"'
            endif
        endfor
        
        ; reread all files.
        tfns = file_search(tfn+'/*')
        if n_elements(tfns) ne 15 then begin
            print, id, n_elements(tfns)
            if n_elements(tfns) ne 3 then $
                cusp_read_plot_data, id, /save_data, /reload
            continue
        endif
        tfns = file_basename(tfns)
        tfns = tfns[sort(tfns)]
        fn0s = id+'_'+['efluxes.pdf','esa.pdf','field_and_poynt_freq_band_scaled.pdf',$
            'field_and_poynt_freq_band.pdf','field_and_poynt.pdf','footpoint.pdf','hydra.pdf',$
            'ke_special.pdf','ke_special.svg','mat_spec_fast.pdf','mat_spec_polar.pdf','overview.pdf']
        fn0s = [fn0s, $
            tid+'_'+['polar_b_despike.pdf','polar_e56dot0_dbcorr.pdf','polar_spc2fac.pdf']]
        fn0s = fn0s[sort(fn0s)]
        for j = 0, n_elements(tfns)-1 do begin
            if tfns[j] ne fn0s[j] then begin
                print, id
                print, tfns[j]
                print, fn0s[j]
                break
            endif
        endfor
    endfor
end
locroot = shomedir()+'/Google Drive/works/works/cusp/cusp list conjun'
;locroot = shomedir()+'/cusp/'
cusp_check_files, locroot
end