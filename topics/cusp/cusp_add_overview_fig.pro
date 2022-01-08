pro cusp_add_overview_fig, locroot

    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun.log'
    
    fns = file_search(locroot+'/*')
    
    for i = 0, n_elements(fns)-1 do begin
        tfn = fns[i]
        if ~file_test(tfn, /directory) then continue
        id = strmid(tfn, strlen(locroot)+1)
        ofn = tfn+'/'+strmid(id,0,9)+'_overview.pdf'
        if file_test(ofn) then file_delete, ofn
        ofn = tfn+'/'+id+'_overview.eps'
        if file_test(ofn) then file_delete, ofn
        ofn = tfn+'/'+id+'_overview.pdf'
        if file_test(ofn) then continue
        info = cusp_read_conjun_list(logfile, event = id)
        if size(info,/type) ne 8 then continue      ; no enough info.
        orootdir = locroot
        
        ; prepare.
        potr = info.polar.plot_time
        fatr = info.fast.plot_time
        potrcusp = info.polar.cusp_time
        fatrcusp = info.fast.cusp_time
        if n_elements(fatrcusp) eq 1 then continue  ; no fast cusp tr.
        if potr[1] lt potr[0] then potr[1]+= 86400d
        if fatr[1] lt fatr[0] then fatr[1]+= 86400d
        if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
        if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d
        
        tinfo = info.polar
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        if file_test(fn) eq 0 then message, 'file does not exist ...'
        sdt = ssdtread(fn)
        pre = 'po_'
        t1  = sdt.var.polarinvariantlatitude.depend_0
        dis  = sdt.var.polarspcraftdist.value
        store_data, pre+'dis', data = {x:t1, y:dis}, limits = {ytitle:'Dist (Re)'}
        sdt = 0d
        
        alltr = []
        
        dt = 600        ; 10 min.
        get_data, 'po_dis', t0, tmp
        if min(interpol(tmp, t0, potrcusp)) le 2.5 then dt = 120    ; 2 min.
        tmp = potrcusp
        ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
        ttr = ttr-(ttr mod dt) & ttr[1]+= dt
        alltr = [min([ttr,alltr]),max([ttr,alltr])]
        
        dt = 60        ; 1 min.
        tmp = fatrcusp
        ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
        ttr = ttr-(ttr mod dt) & ttr[1]+= dt
        alltr = [min([ttr,alltr]),max([ttr,alltr])]
        
        ofn = orootdir+'/'+id+'/'+id+'_overview.eps'
        sgpsopen, ofn, xsize = 6, ysize = 8, /inch
        plot_polar_fast_summary, stoepoch(alltr,'unix'), /no_delete
        ;    timebar, [potrcusp,fatrcusp], color = red, thick = 2
        sgpsclose, /pdf
    endfor
end

locroot = shomedir()+'/Google Drive/works/works/cusp/cusp list conjun'
cusp_add_overview_fig, locroot
end