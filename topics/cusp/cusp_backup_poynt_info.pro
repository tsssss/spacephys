
pro cusp_backup_poynt_info

    type = 'south_imf'
    
    ifn = shomedir()+'/Google Drive/works/data/cusp/cusp_stats_'+type+'.tplot'
    tvar = 'cusp_stats_'+type
    
    if file_test(ifn) eq 0 then begin
        ids = cusp_id(type)
        fns = shomedir()+'/Google Drive/works/data/cusp/'+ids+'_all_data.tplot'
        nfn = n_elements(fns)
        tinfo = cusp_info_struct()
        infos = replicate(tinfo, nfn)
        ninfo = nfn
        
        for i = 0, nfn-1 do begin
            printf, -1, fns[i]
            tplot_restore, filename = fns[i]
            get_data, 'scidat', tmp, info
            infos[i] = info
        endfor
        
        store_data, tvar, ninfo, infos
        tplot_save, tvar, filename = ifn
    endif else begin
        tplot_restore, filename = ifn
        get_data, tvar, ninfo, infos
    endelse
    
    ; convert to from high latitude to low latitude, negative pperp means into cusp.
    po_dis = (infos.polar.cusp.entry.dis+infos.polar.cusp.exit.dis)*0.5
    po_pflux = infos.polar.sb.fh
    po_pperp = infos.polar.sv.fh
    flag = infos.polar.cusp.entry.ilat-infos.polar.cusp.exit.ilat
    idx = where(flag le 0, cnt)
    if cnt ne 0 then po_pperp[idx]*= -1
    fa_dis = (infos.fast.cusp.entry.dis+infos.fast.cusp.exit.dis)*0.5
    fa_pflux = infos.fast.sb.fh
    fa_pperp = infos.fast.sv.fh
    flag = infos.fast.cusp.entry.ilat-infos.fast.cusp.exit.ilat
    idx = where(flag le 0, cnt)
    if cnt ne 0 then po_pperp[idx]*= -1
    
    dis = [po_dis,fa_dis]
    pflux = [po_pflux,fa_pflux]
    pperp = [po_pperp,fa_pperp]
    pratio = pperp/pflux
    
    print, mean(pratio)
    
    stop
    

end