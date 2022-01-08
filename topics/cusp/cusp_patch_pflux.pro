;+
; patch to reset poynting flux for given event id.
;-
pro cusp_patch_pflux, id, fast = fast, polar = polar

    if keyword_set(fast) then forpolar = 0 else forpolar = 1
    if keyword_set(polar) then forpolar = 1 else forpolar = 0

    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

    
; **** read basic info.
    ; log info, search on conjunction list, then on polar list.
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event = id)
    if size(loginfo,/type) ne 8 then begin
        logfile = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'
        loginfo = cusp_read_conjun_list(logfile, event = id, /nofast)
    endif

    if size(loginfo,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif

    print, 'processing event: '+id

    ; info, id, pre0, [de,db,pf]idx, plottr, cusptr.
    info = forpolar? loginfo.polar: loginfo.fast
    pre0 = forpolar? 'po_': 'fa_'

    deidx = 0   ; dEv.
    dbidx = 1   ; dBp.
    pfidx = 2   ; Sb.

    plottr = info.plot_time
    if plottr[1] lt plottr[0] then plottr[1]+= 86400d

    cusptr = info.cusp_time
    if cusptr[1] lt cusptr[0] then cusptr[1]+= 86400d

    ; filt0, filts, delims, facdelim, tscls.
    filt0 = info.filters    ; original filters.
    filts = filt0           ; filters, small > large.
    filts = filts[uniq(filts,sort(filts))]

    delims = info.noisedelim    ; below and above are noise.
    delims = minmax(delims)     ; sort.
    if delims[0] lt 0 then delims[0] = 0
    if delims[1] lt 0 then delims[1] = plottr[1]-plottr[0]

    facdelim = info.faclim      ; above is fac.
    if facdelim eq -1 then facdelim = max(filts)

    tmp = info.scaleinfo        ; preliminary, need datarate to modify.
    tscls = smkgmtrc(tmp[0],tmp[1],tmp[2],'n')

; **** load field data.
    ; [pre0]_[de,db]_fac, pre0_[ilat,mlt,dis], original field.
    ; [pre0]_tscl. scales in time in second.
    ; [pre0]_rscl. scales in # of records.

    if forpolar then begin
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        polar_sdt_prep_poynting_flux, fn, e56 = info.e56, $
            tr = plottr, eventid = id, /noplot
        vars = pre0+['b0_spc','spc2fac','db_spc','de_spc','mlat']
        store_data, vars, /delete
    endif else begin
        fn = rootdir+'/data/cusp/fa_sdt_fld_'+id+'_'+ $
            string(info.orbit,format='(I05)')+'.tplot'
        fast_sdt_prep_poynting_flux, fn, $  ; plot time add padding time.
            trplot = plottr+(plottr[1]-plottr[0])*[-1,1]
        vars = pre0+['pos','vel','alt','alt','b0_gei','b0_gei']
        store_data, vars, /delete
    endelse

    ; dr, rscls, tscls. data rate, scales in # of rec and time.
    get_data, pre0+'de_fac', t0
    dr = sdatarate(t0)      ; data rate.
    tmp = round(tscls/dr)
    rscls = double(tmp[uniq(tmp)])
    tscls = rscls*dr

; **** calc mat bands.
    ; mat filter.
    maxfilt = 1e31                      ; the last band is fac.
    matfilts = [0,filts,maxfilt]        ; filters for mat spectrogram.
    nmatband = n_elements(matfilts)-1
    matbandids = '_f'+string(indgen(nmatband),format='(I0)')

    dename = pre0+'de_fac'
    dbname = pre0+'db_fac'
    pfname = pre0+'pf_fac'
    stplot_calc_pflux_mat, dename, dbname, pfname, tscale = tscls, $
        filter = matfilts, ids = matbandids

    ; put all background into the fac band.
    vars = [dename,dbname]+'_mat'
    for i = 0, n_elements(vars)-1 do begin
        tmp = 0d
        for j = 0, nmatband-2 do begin
            get_data, vars[i]+matbandids[j], t0, dat
            tmp+= dat
        endfor
        get_data, vars[i], t0, dat
        store_data, vars[i]+matbandids[nmatband-1], t0, dat-tmp
    endfor

    remove_wake = 1
    ; remove background for the fac band.
    if keyword_set(remove_wake) then begin
        dename = pre0+'de_fac_mat'
        dbname = pre0+'db_fac_mat'
        vars = [dename,dbname]
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            idx = nmatband-1
            tvar = vars[i]+matbandids[idx]
            get_data, tvar, t0, dat
            bg = scalcbg(dat)
            store_data, tvar, t0, dat-bg
            tvar = vars[i]+'_nl'
            store_data, tvar, t0, bg, limits = {labels:['bg','','']}
        endfor
        nmatband+= 1
        matbandids = [matbandids,'_nl']
    endif

    ; update poynting flux.
    dename = pre0+'de_fac_mat'
    dbname = pre0+'db_fac_mat'
    pfname = pre0+'pf_fac_mat'
    vars = remove_wake? matbandids[nmatband-2:*]: matbandids[nmatband-1:*]
    for i = 0, n_elements(vars)-1 do $
        stplot_calc_pflux_mat, dename+vars[i], dbname+vars[i], pfname+vars[i]
    store_data, dename+'_*_comp'+['1','2','3'], /delete
    store_data, dbname+'_*_comp'+['1','2','3'], /delete
    ; total used field, exclude fac band and high freq noise.
    vars = pre0+['de','db','pf']+'_fac_mat'
    nvar = n_elements(vars)
    for i = 0, nvar-1 do begin
        tvar = tnames(vars[i]+'_f?')
        tvar = tvar[1:n_elements(tvar)-2]   ; exclude fac & high freq band.
        stplot_total, tvar, newname = vars[i]
    endfor


; **** to this point, we get the following tplot vars.
; [pre0]_[de,db,pf]_fac_mat_[f0,f1,f2,...,nl]. The components.
; [pre0]_[de,db,pf]_mat. The sum.
; [pre0]_[de,db]_fac.
; [pre0]_[b,ilat,mlt,dis].


    ; map and integrate.
    vars = pre0+'pf_fac_mat'+['',matbandids]
    get_data, vars[0], t0
    for j = 0, n_elements(vars)-1 do begin
        smap2iono, vars[j], pre0+'dis', b = pre0+'b', $
            newname = vars[j], coef = coef
        cusp_int_eflux, vars[j], pre0+'ilat', tr = cusptr
        tplot_rename, vars[j], vars[j]+'_map'
    endfor
    store_data, pre0+'map_coef', t0, coef
    

    ; summarize the results.
    ; [v,p,b] components are [0,1,2] in indices.
    ; the first band is high freq noise, the last band is background.
    ; the second last is fac.
    pfbands = dblarr(nmatband,3)        ; all the line-integrated S.
    pfstars = dblarr(nmatband)          ; coef for S*.
    ebratio = dblarr(nmatband,3)        ; e/b ratio.
    for i = 0, nmatband-1 do begin
        suf = matbandids[i]
        get_data, pre0+'pf_fac_mat'+suf+'_map', tmp, tmp, dat
        pfbands[i,*] = dat
        dename = pre0+'de_fac_mat'+suf
        dbname = pre0+'db_fac_mat'+suf

        ; calc e/b ratio and coef for S*.
        ebratio[i,*] = stplot_ebratio(dename, dbname, $
            deidx = deidx, dbidx = dbidx, trange = cusptr, method = 'minmax')
        pfstars[i] = stplot_pfstar(dename, dbname, trange = cusptr)
    endfor

    ; clean up. there are conflicts if do not delete the vars.
    tmp = ['ele_keflux','ion_keflux','pf_fac_mat_para']
    vars = ['ilat','mlt','dis','de_fac','db_fac','pf_fac_mat',tmp,tmp+'_map']
    scivars = ['po_'+vars,'fa_'+vars,'scidat']
    store_data, scivars, /delete
    
    ; apply change to event info. 
    datdir = rootdir+'/data/cusp'
    infofn = datdir+'/'+id+'_all_data.tplot'
    infovar = 'scidat'
    tplot_restore, filename = infofn
    
    get_data, infovar, scidatt0, scidat
    
    tmp = forpolar? 'POLAR': 'FAST'
    idx = where(tag_names(scidat) eq tmp)
    scidat.(idx).ebratio = ebratio
    scidat.(idx).pfstar = pfstars
    scidat.(idx).filters = matfilts
    scidat.(idx).nfilter = nmatband
    scidat.(idx).facdelim = facdelim
    scidat.(idx).noisedelim = delims

    scidat.(idx).sv.fs = pfbands[*,0]
    scidat.(idx).sp.fs = pfbands[*,1]
    scidat.(idx).sb.fs = pfbands[*,2]

    scidat.(idx).sv.fh = total(pfbands[1:nmatband-3,0])
    scidat.(idx).sp.fh = total(pfbands[1:nmatband-3,1])
    scidat.(idx).sb.fh = total(pfbands[1:nmatband-3,2])
    scidat.(idx).sv.fl = pfbands[nmatband-2,0]
    scidat.(idx).sp.fl = pfbands[nmatband-2,1] 
    scidat.(idx).sb.fl = pfbands[nmatband-2,2] 

    hem = scidat.(idx).hem
    if hem lt 0 then begin
        scidat.(idx).sb.fs*= -1
        scidat.(idx).sb.fh*= -1
        scidat.(idx).sb.fl*= -1
        scidat.(idx).sp.fs*= -1
        scidat.(idx).sp.fh*= -1
        scidat.(idx).sp.fl*= -1
    endif

    ; efluxes.
    scidat.(idx).pflux = scidat.(idx).sb.fh
    scidat.(idx).eflux = scidat.(idx).kei+scidat.(idx).kee+scidat.(idx).pflux
    scidat.ratio_pflux = scidat.polar.pflux/scidat.fast.pflux
    scidat.ratio_eflux = scidat.polar.eflux/scidat.fast.eflux

    store_data, infovar, scidatt0, scidat

    tmp = file_dirname(infofn)
    if ~file_test(tmp,/directory) then file_mkdir, tmp

    tmp = ['ele_keflux','ion_keflux','pf_fac_mat_para']
    tplot_save, scivars, filename = infofn
    store_data, scivars, /delete
    
end

ids = cusp_calc_id(cusp_id('south_imf'),cusp_id('polar_south_imf'))
;foreach id, ids do cusp_patch_pflux, id, /polar
ids = cusp_id('south_imf')
;foreach id, ids do cusp_patch_pflux, id, /fast

end
