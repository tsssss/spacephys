;+
; Use the way to estimate E/B ratio to get a measurement for E and B.
; Calc the ratio Rb = Bp/Bv and Re = Ev/Ep. In the idea case, they equal.
; The corrected poynting flux S_star = S*Rb/Re.
; We calculate the correction for each frequency bands.
;-

eventid = '1998_1001_03'
test = 1


    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    loginfo = cusp_read_conjun_list(logfile, event = eventid)
    if size(loginfo,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif
    
    
    ; directory to save output data and pdfs.
    figdir = rootdir+'/works/cusp/cusp list conjun'
    datdir = rootdir+'/data/cusp'

    if keyword_set(test) then begin
        figdir = shomedir()+'/cusp'
        datdir = shomedir()+'/cusp/data'
    endif
    
    if ~file_test(figdir,/directory) then file_mkdir, figdir
    

    ; reload and recalc data?
    if (tnames('*'))[0] eq '' then reload = 1   ; no data loaded.
    get_data, 'scidat', t0, info
    if size(dat,/type) eq 8 then if eventid ne info.id then reload = 1   ; new.




; **** basic info for both loading data and plotting.
    potr = loginfo.polar.plot_time
    fatr = loginfo.fast.plot_time
    potrcusp = loginfo.polar.cusp_time
    fatrcusp = loginfo.fast.cusp_time
    if potr[1] lt potr[0] then potr[1]+= 86400d
    if fatr[1] lt fatr[0] then fatr[1]+= 86400d
    if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
    if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d
    id = loginfo.id

    ; filters, min to max, uniq.
    ; Will take the minimum and maximum filters as noise delimiters,
    ; if delimiters are excplitly set, then use them.
    ; Will take the second large filter as delimiter for fac and wave,
    ; if fac delimiter is set, then use it.

    ; filt0, original filter: f1,f2,f3,f4.
    ; We also want to include 0 and infinity (maxf).
    ;
    ; delims, rule out noise: n1, n2.
    ;
    ; Say n1=f2, n2=f4, the wanted filters are f2, f3, f4, ie, f1',f2',f3'.
    ; The filters used to filter mat spectrogram are: 0,f1',f2',f3',maxf.
    ; There are 5 filters, 4 bands, among the bands, 3 are wanted.
    ; High freq noise band is 0-f1', low freq noise band is f3'-maxf.
    ;
    ; faclim, sets where fac bands are.
    ; Say faclim=f2', then f1'-faclim are wave bands, faclim-f3' are fac bands.

    maxfilt = 1e31

    ; polar.
    pofilt0 = loginfo.polar.filters     ; original filters.
    pofaclim = loginfo.polar.faclim     ; above is fac.
    podelims = loginfo.polar.noisedelim ; below and above are noise.
    if (where(podelims lt 0))[0] eq -1 then podelims = minmax(podelims) ; sort if no -1.

    ; sort from small to large.
    pofilts = pofilt0
    pofilts = pofilts[sort(pofilts)]
    pofilts = pofilts[uniq(pofilts)]
    if pofilts[0] ne 0 then pofilts = [0,pofilts]   ; auto fill 0.

    ; default setting: [0,max original filters].
    if podelims[0] lt 0 then podelims[0] = min(pofilts)
    if podelims[1] lt 0 then podelims[1] = maxfilt

    ; get the wanted filters, add noise delimiters, and fac limit.
    idx = where(pofilts ge podelims[0] and pofilts le podelims[1], ponfilt)
    pofilts = pofilts[idx]
    if min(pofilts) gt podelims[0] then pofilts = [podelims[0],pofilts]
    if max(pofilts) lt podelims[1] then pofilts = [pofilts,podelims[1]]
    
    ; use faclim to omit unecessary filters.
    if pofaclim lt 0 then pofaclim = pofilts[n_elements(pofilts)-2]
    if pofaclim gt podelims[1] then pohasfac = 0 else begin
        pohasfac = 1
        pofacidx = where(pofilts lt pofaclim)
        pofilts = [pofilts[pofacidx],pofaclim,max(pofilts)]
    endelse

    ; includes wave bands and fac (if pohasfac=1).
    ponfilt = n_elements(pofilts)

    ; the filters used in filtering in mat spectrogram.
    pomatfilts = pofilts
    if min(pomatfilts) gt 0 then begin
        pomatfilts = [0,pomatfilts]
        pohashighnoise = 1
    endif else pohashighnoise = 0
    if max(pomatfilts) lt maxfilt then begin
        pomatfilts = [pomatfilts,maxfilt]
        pohaslownoise = 1
    endif else pohaslownoise = 0
    ponmatband = n_elements(pomatfilts)-1
    pomatbandids = string(indgen(ponmatband),format='(I0)')

    ; the wanted bands are within the wanted filters.
    ponband = ponfilt-1     ; includes wave bands and fac, separate later.
    pobandids = 'b'+string(indgen(ponband),format='(I0)')    ; in var name.
    pobandlabels = strarr(ponband)
    for i = 0, ponband-1 do pobandlabels[i] = $
        sgnum2str(pofilts[i],msgn=3)+'-'+sgnum2str(pofilts[i+1],msgn=3)+'s'


    ; fast.
    fafilt0 = loginfo.fast.filters     ; original filters.
    fafaclim = loginfo.fast.faclim     ; above is fac.
    fadelims = loginfo.fast.noisedelim ; below and above are noise.
    if (where(fadelims lt 0))[0] eq -1 then fadelims = minmax(fadelims) ; sort if no -1.

    ; sort from small to large.
    fafilts = fafilt0
    fafilts = fafilts[sort(fafilts)]
    fafilts = fafilts[uniq(fafilts)]
    if fafilts[0] ne 0 then fafilts = [0,fafilts]   ; auto fill 0.

    ; default setting: [0,max original filters].
    if fadelims[0] lt 0 then fadelims[0] = min(fafilts)
    if fadelims[1] lt 0 then fadelims[1] = maxfilt

    ; get the wanted filters, add noise delimiters, and fac limit.
    idx = where(fafilts ge fadelims[0] and fafilts le fadelims[1], fanfilt)
    fafilts = fafilts[idx]
    if min(fafilts) gt fadelims[0] then fafilts = [fadelims[0],fafilts]
    if max(fafilts) lt fadelims[1] then fafilts = [fafilts,fadelims[1]]
    
    ; use faclim to omit unecessary filters.
    if fafaclim lt 0 then fafaclim = fafilts[n_elements(fafilts)-2]
    if fafaclim gt fadelims[1] then fahasfac = 0 else begin
        fahasfac = 1
        fafacidx = where(fafilts lt fafaclim)
        fafilts = [fafilts[fafacidx],fafaclim,max(fafilts)]
    endelse

    ; includes wave bands and fac (if fahasfac=1).
    fanfilt = n_elements(fafilts)

    ; the filters used in filtering in mat spectrogram.
    famatfilts = fafilts
    if min(famatfilts) gt 0 then begin
        famatfilts = [0,famatfilts]
        fahashighnoise = 1
    endif else fahashighnoise = 0
    if max(famatfilts) lt maxfilt then begin
        famatfilts = [famatfilts,maxfilt]
        fahaslownoise = 1
    endif else fahaslownoise = 0
    fanmatband = n_elements(famatfilts)-1
    famatbandids = string(indgen(fanmatband),format='(I0)')

    ; the wanted bands are within the wanted filters.
    fanband = fanfilt-1     ; includes wave bands and fac, separate later.
    fabandids = 'b'+string(indgen(fanband),format='(I0)')    ; in var name.
    fabandlabels = strarr(fanband)
    for i = 0, fanband-1 do fabandlabels[i] = $
        sgnum2str(fafilts[i],msgn=3)+'-'+sgnum2str(fafilts[i+1],msgn=3)+'s'



    ; the wanted component, 0-based counting.
    deidx = 0   ; dEv.
    dbidx = 1   ; dBp.
    pfidx = 2   ; Sb.
    dematsuf = '_comp'+string(deidx,format='(I0)')
    dbmatsuf = '_comp'+string(dbidx,format='(I0)')
    
    ndim = 3    ; 3-d fields.



; **** plot settings and generate plots to disk.
    posl = [0.1,0.10,0.40,0.90]
    posr = [0.6,0.10,0.90,0.90]
    faclabs = ['v','p','b']
    falabs = 'fa_'+['ilat','mlt','dis']
    polabs = 'po_'+['ilat','mlt','dis']
    rgb = [6,4,2] & red = 6
    charsz = 1
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'charsize', charsz*0.9
    tplot_options, 'xcharsize', charsz*0.7
    tplot_options, 'ycharsize', charsz*0.8
    tplot_options, 'zcharsize', charsz*0.6
    tplot_options, 'constant', 0
    time_stamp, /off


; **** calculate all data and save to tplot.
    ; poynting flux, KE flux, integrated fluxes.
    if keyword_set(reload) then begin
        store_data, '*', /delete

    ; **** data structure.
        info = cusp_info_struct()

        ; dst, ae, imf bz, by.
        conjtr = minmax([potrcusp,fatrcusp])
        dt = 0 & conjtr[0]-= dt & conjtr[1]+= dt    ; allow a padding time.
        tmp = sread_omni(conjtr)
        tmp2 = sfmepoch(tmp.epoch,'unix')
        idx = where(tmp2 ge conjtr[0] and tmp2 le conjtr[1])
        info.ae = max(tmp.ae[idx])
        info.dst = min(tmp.symh[idx])
        info.imfbz = minmax(tmp.bgse[*,2])
        info.imfby = minmax(tmp.bgse[*,1])

        info.id = id
        info.conjtr = conjtr

        info.polar.cusp.entry.ut = potrcusp[0]
        info.polar.cusp.exit.ut = potrcusp[1]
        info.polar.filters = pofilt0
        info.polar.nfilter = n_elements(pofilt0)
        info.polar.faclimit = pofaclim
        info.polar.noisedelim = podelims

        info.fast.cusp.entry.ut = fatrcusp[0]
        info.fast.cusp.exit.ut = fatrcusp[1]
        info.fast.filters = fafilt0
        info.fast.nfilter = n_elements(fafilt0)
        info.fast.faclimit = fafaclim
        info.fast.noisedelim = fadelims
        

    ; **** read polar data.
        poinfo = loginfo.polar

        ; kinetic energy flux.
        ; po_[ion,ele]_keflux.
        fn = rootdir+'/works/cusp/cusp list conjun/'+id+'/'+id+'_ke_special.svg'
        yrange = poinfo.keflux    ; [ele,ion].
        type = poinfo.ketype   ; etp,p2p,ptp.
        case type of
            'etp':ex = {epstopdf:1}
            'ptp':ex = {pstopdf:1}
            'p2p':ex = {ps2pdf:1}
        endcase
        polar_read_ke_flux, fn, potr, yrange[0]*[-1,1], yrange[1]*[-1,1], _extra = ex

        ; field data.
        ; po_[de,db]_fac, original field.
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        polar_sdt_prep_poynting_flux, fn, e56 = poinfo.e56, tr = potr, eventid = id, $
            cusptr = potrcusp, orootdir = figdir+'/'+id, titpre = 'Event ID: '+id+'    '
        vars = ['po_b','po_b0_spc','po_spc2fac','po_db_spc','po_de_spc','po_mlat']
        store_data, vars, /delete

        ; poynting flux.
        ; po_[de,db]_fac_mat, 3-D fields used in Poynting flux calc.
        ; po_[de,db]_fac_mat[1,2,3,...], 3-D fields in freq bands.
        dename = 'po_de_fac' & dbname = 'po_db_fac' & pfname = 'po_pf_fac'
        stplot_calc_pflux_mat, dename, dbname, pfname, tscale = potscl, $
            filter = pomatfilts, ids = pomatbandids, scaleinfo = poinfo.scaleinfo

        ; group poynting flux into bands.
        ; high freq noise band.
        vars = [dename,dbname,pfname]+'_mat'
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_nh' 
            if pohashighnoise then begin
                stplot_renew, vars[i]+pomatbandids[0], newname = tvar, /delete
            endif else store_data, tvar, potr, dblarr(2,3)
        endfor
        ; low freq noise band.
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_nl' 
            if pohaslownoise then begin
                stplot_renew, vars[i]+pomatbandids[ponmatband-1], newname = tvar, /delete
            endif else store_data, tvar, potr, dblarr(2,3)
        endfor
        ; wave bands.
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_'+pobandids
            for j = 0, ponband-1 do stplot_renew, $
                vars[i]+pomatbandids[pohashighnoise+j], newname = tvar[j], /delete
        endfor
        ; remove background for the lowest freq wave band.
        vars = [dename,dbname]+'_mat'
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_'+pobandids[ponband-1]
            get_data, tvar, t0, dat
            bg = scalcbg(dat)
            store_data, tvar, t0, dat-bg
            tvar = vars[i]+'_nl'
            if pohaslownoise then begin
                get_data, tvar, t0, dat
                store_data, tvar, t0, dat+bg
            endif else store_data, tvar, t0, bg
        endfor
        ; update poynting flux
        tmp = '_mat_'+pobandids[ponband-1]
        stplot_calc_pflux_mat, dename+tmp, dbname+tmp, pfname+tmp
        ; total used field, exclude low and high freq noise.
        vars = [dename,dbname,pfname]+'_mat'
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            tvar = vars[i]+['_'+pobandids]
            if pohasfac then tvar = vars[i]+['_'+pobandids[0:ponband-2]]    ; exclude fac band.
            stplot_total, tvar, newname = vars[i]
        endfor

        ; for all energy fluxes, map and integrate.
        ; flux -> flux_map.
        vars = ['po_'+['ele_keflux','ion_keflux','pf_fac_mat'+['','_'+pobandids]]]
        for j = 0, n_elements(vars)-1 do begin
            smap2iono, vars[j], 'po_dis', newname = vars[j]+'_map'
            cusp_int_eflux, vars[j]+'_map', 'po_ilat', tr = potrcusp
        endfor
        ; store parallel poynting flux.
        pfvars = ['po_pf_fac_mat'+['','_map']]
        for j = 0, n_elements(pfvars)-1 do begin
            get_data, pfvars[j], t0, tmp, tmp2
            idx = (n_elements(tmp2) eq 1)? 0: 2
            store_data, pfvars[j]+'_para', t0, tmp[*,2], tmp2[idx], $
                limits = {ytitle:'(mW/m!U2!N)',labels:'S!D||!N'}
        endfor
        stplot_renew, 'po_pf_fac_mat_map_para', $
            newname = 'po_pf_fac_mat_para_map', /delete


    ; **** read fast data.
        fainfo = loginfo.fast

        ; kinetic energy flux.
        ; fa_[ion,ele]_keflux.
        fn = rootdir+'/data/cusp/fa_sdt_esa_'+id+'_'+ $
            string(fainfo.orbit,format='(I05)')+'.tplot'
        tplot_restore, filename = fn
        stplot_renew, 'ele_eflux', newname = 'fa_ele_keflux', /delete
        stplot_renew, 'ion_eflux', newname = 'fa_ion_keflux', /delete
        vars = ['ion_*','ele_*']
        store_data, vars, /delete
        ; despike for ion ke. !!! need better despike algorithm.
        get_data, 'fa_ion_keflux', t0, dat
        idx = where(abs(dat) le 1, cnt)
        store_data, 'fa_ion_keflux', t0, interpol(dat[idx],t0[idx],t0)

        ; field data.
        ; fa_[de,db]_fac, original field.
        fn = rootdir+'/data/cusp/fa_sdt_fld_'+id+'_'+ $
            string(fainfo.orbit,format='(I05)')+'.tplot'
        if file_search(fn) eq '' then $
            fn = rootdir+'/data/cusp/fa_sdt_fld_'+id+'_'+$
            string(fainfo.orbit,format='(I05)')+'.tplot'
        fast_sdt_prep_poynting_flux, fn, $  ; plot time add padding time.
            trplot = fatr+(fatr[1]-fatr[0])*[-1,1]
        vars = ['fa_pos','fa_vel','fa_alt','alt','fa_b0_gei','fa_b0_gei','fa_b']
        store_data, vars, /delete

        ; poynting flux.
        ; fa_[de,db]_fac_mat, 3-D fields used in poynting flux calc.
        ; fa_[de,db]_fac_matf[1,2,3,...], 3-D fields in freq bands.
        dename = 'fa_de_fac' & dbname = 'fa_db_fac' & pfname = 'fa_pf_fac'
        stplot_calc_pflux_mat, dename, dbname, pfname, tscale = fatscl, $
            filter = famatfilts, ids = famatbandids, scaleinfo = fainfo.scaleinfo

        ; group poynting flux into bands.
        ; high freq noise band.
        vars = [dename,dbname,pfname]+'_mat'
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_nh' 
            if fahashighnoise then begin
                stplot_renew, vars[i]+famatbandids[0], newname = tvar, /delete
            endif else store_data, tvar, fatr, dblarr(2,3)
        endfor
        ; low freq noise band.
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_nl' 
            if fahaslownoise then begin
                stplot_renew, vars[i]+famatbandids[fanmatband-1], newname = tvar, /delete
            endif else store_data, tvar, fatr, dblarr(2,3)
        endfor
        ; wave bands.
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_'+fabandids
            for j = 0, fanband-1 do stplot_renew, $
                vars[i]+famatbandids[fahashighnoise+j], newname = tvar[j], /delete
        endfor
        ; remove background for the lowest freq wave band.
        vars = [dename,dbname]+'_mat'
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            tvar = vars[i]+'_'+fabandids[fanband-1]
            get_data, tvar, t0, dat
            bg = scalcbg(dat)
            store_data, tvar, t0, dat-bg
            tvar = vars[i]+'_nl'
            if fahaslownoise then begin
                get_data, tvar, t0, dat
                store_data, tvar, t0, dat+bg
            endif else store_data, tvar, t0, bg
        endfor
        ; update poynting flux
        tmp = '_mat_'+fabandids[fanband-1]
        stplot_calc_pflux_mat, dename+tmp, dbname+tmp, pfname+tmp
        ; total used field, exclude low and high freq noise.
        vars = [dename,dbname,pfname]+'_mat'
        nvar = n_elements(vars)
        for i = 0, nvar-1 do begin
            tvar = vars[i]+['_'+fabandids]
            if fahasfac then tvar = vars[i]+['_'+fabandids[0:fanband-2]]    ; exclude fac band.
            stplot_total, tvar, newname = vars[i]
        endfor

        ; for all energy fluxes, map and integrate.
        ; flux -> flux_map.
        vars = ['fa_'+['ele_keflux','ion_keflux','pf_fac_mat'+['','_'+fabandids]]]
        for j = 0, n_elements(vars)-1 do begin
            smap2iono, vars[j], 'fa_dis', newname = vars[j]+'_map'
            cusp_int_eflux, vars[j]+'_map', 'fa_ilat', tr = fatrcusp
        endfor
        ; store parallel poynting flux.
        pfvars = ['fa_pf_fac_mat'+['','_map']]
        for j = 0, n_elements(pfvars)-1 do begin
            get_data, pfvars[j], t0, tmp, tmp2
            idx = (n_elements(tmp2) eq 1)? 0: 2
            store_data, pfvars[j]+'_para', t0, tmp[*,2], tmp2[idx], $
                limits = {ytitle:'(mW/m!U2!N)',labels:'S!D||!N'}
        endfor
        stplot_renew, 'fa_pf_fac_mat_map_para', $
            newname = 'fa_pf_fac_mat_para_map', /delete


    ; **** fill info structure.
        ; polar.
        pre = 'po_'
        trcusp = potrcusp
        get_data, pre+'ilat', t0, dat
        info.polar.cusp.entry.ilat = interpol(dat,t0,trcusp[0])
        info.polar.cusp.exit.ilat = interpol(dat,t0,trcusp[1])
        get_data, pre+'mlt', t0, dat
        info.polar.cusp.entry.mlt = interpol(dat,t0,trcusp[0])
        info.polar.cusp.exit.mlt = interpol(dat,t0,trcusp[1])
        get_data, pre+'dis', t0, dat
        info.polar.cusp.entry.dis = interpol(dat,t0,trcusp[0])
        info.polar.cusp.exit.dis = interpol(dat,t0,trcusp[1])
    
        hem = info.polar.cusp.entry.ilat & hem = hem/abs(hem)
        info.polar.hem = hem
    
        get_data, pre+'ele_keflux_map', tmp, tmp2, dat
        info.polar.kee = dat
        get_data, pre+'ion_keflux_map', tmp, tmp2, dat
        info.polar.kei = dat
    
        pfbands = dblarr(ponband,3)
        for i = 0, ponband-1 do begin
            suf = '_'+pobandids[i]
            get_data, pre+'pf_fac_mat'+suf+'_map', tmp, tmp2, dat
            pfbands[i,*] = dat
            tmp = stplot_ebratio(pre+'de_fac_mat'+suf, pre+'db_fac_mat'+suf, $
                deidx = deidx, dbidx = dbidx, trange = potrcusp, method = 'minmax')
            info.polar.ebratio[i,*] = tmp[*]
            info.polar.pfstar[i] = stplot_pfstar(pre+'de_fac_mat'+suf, pre+'db_fac_mat'+suf, trange = potrcusp)
        endfor
        info.polar.sv.fs = pfbands[*,0]
        info.polar.sp.fs = pfbands[*,1]
        info.polar.sb.fs = pfbands[*,2]
        if pohasfac then begin
            info.polar.sv.fh = total(pfbands[0:ponband-2,0])
            info.polar.sp.fh = total(pfbands[0:ponband-2,1])
            info.polar.sb.fh = total(pfbands[0:ponband-2,2])
            info.polar.sv.fl = total(pfbands[ponband-1:*,0])
            info.polar.sp.fl = total(pfbands[ponband-1:*,1])
            info.polar.sb.fl = total(pfbands[ponband-1:*,2])
        endif else begin
            info.polar.sv.fh = total(pfbands[*,0])
            info.polar.sp.fh = total(pfbands[*,1])
            info.polar.sb.fh = total(pfbands[*,2])
            info.polar.sv.fl = 0d
            info.polar.sp.fl = 0d
            info.polar.sb.fl = 0d
        endelse
    
        ; ion ratio.
        get_data, pre+'ion_keflux_map', tmp, tmp2
        info.max_kei = abs(min(hem*tmp2[where(tmp ge potrcusp[0] and tmp le potrcusp[1])]))

        vars = pre+'ion_keflux_map_abs'
        store_data, vars, tmp, abs(tmp2)
        cusp_int_eflux, vars, pre+'ilat', tr = potrcusp
        get_data, vars, tmp, tmp2, dat
        info.ratio_ion = info.polar.kei/dat  ; -1 for all upward, 1 for all down.
    
        ; fast.
        pre = 'fa_'
        trcusp = fatrcusp
        get_data, pre+'ilat', t0, dat
        info.fast.cusp.entry.ilat = interpol(dat,t0,trcusp[0])
        info.fast.cusp.exit.ilat = interpol(dat,t0,trcusp[1])
        get_data, pre+'mlt', t0, dat
        info.fast.cusp.entry.mlt = interpol(dat,t0,trcusp[0])
        info.fast.cusp.exit.mlt = interpol(dat,t0,trcusp[1])
        get_data, pre+'dis', t0, dat
        info.fast.cusp.entry.dis = interpol(dat,t0,trcusp[0])
        info.fast.cusp.exit.dis = interpol(dat,t0,trcusp[1])
    
        hem = info.fast.cusp.entry.ilat & hem = hem/abs(hem)
        info.fast.hem = hem
    
        get_data, pre+'ele_keflux_map', tmp, tmp2, dat
        info.fast.kee = dat
        get_data, pre+'ion_keflux_map', tmp, tmp2, dat
        info.fast.kei = dat
    
        pfbands = dblarr(fanband,3)
        for i = 0, fanband-1 do begin
            suf = '_'+fabandids[i]
            get_data, pre+'pf_fac_mat'+suf+'_map', tmp, tmp2, dat
            pfbands[i,*] = dat
            tmp = stplot_ebratio(pre+'de_fac_mat'+suf, pre+'db_fac_mat'+suf, $
                deidx = deidx, dbidx = dbidx, trange = fatrcusp, method = 'minmax')
            info.fast.ebratio[i,*] = tmp[*]
            info.fast.pfstar[i] = stplot_pfstar(pre+'de_fac_mat'+suf, pre+'db_fac_mat'+suf, trange = fatrcusp)
        endfor

        info.fast.sv.fs = pfbands[*,0]
        info.fast.sp.fs = pfbands[*,1]
        info.fast.sb.fs = pfbands[*,2]
        if fahasfac then begin
            info.fast.sv.fh = total(pfbands[0:fanband-2,0])
            info.fast.sp.fh = total(pfbands[0:fanband-2,1])
            info.fast.sb.fh = total(pfbands[0:fanband-2,2])
            info.fast.sv.fl = total(pfbands[fanband-1:*,0])
            info.fast.sp.fl = total(pfbands[fanband-1:*,1])
            info.fast.sb.fl = total(pfbands[fanband-1:*,2])
        endif else begin
            info.fast.sv.fh = total(pfbands[*,0])
            info.fast.sp.fh = total(pfbands[*,1])
            info.fast.sb.fh = total(pfbands[*,2])
            info.fast.sv.fl = 0d
            info.fast.sp.fl = 0d
            info.fast.sb.fl = 0d
        endelse

        ; convert to up/down, + is down, - is up. It's a rotation around v-axis.
        if hem lt 0 then begin
            info.ratio_ion*= -1
    
            info.polar.kei*= -1
            info.polar.kee*= -1
            info.polar.sb.fs*= -1
            info.polar.sb.fh*= -1
            info.polar.sb.fl*= -1
            info.polar.sp.fs*= -1
            info.polar.sp.fh*= -1
            info.polar.sp.fl*= -1
    
            info.fast.kei*= -1
            info.fast.kee*= -1
            info.fast.sb.fs*= -1
            info.fast.sb.fh*= -1
            info.fast.sb.fl*= -1
            info.fast.sp.fs*= -1
            info.fast.sp.fh*= -1
            info.fast.sp.fl*= -1
        endif

        ; relative position.
        info.dt = (info.fast.cusp.entry.ut-info.polar.cusp.entry.ut)/3600d
        info.dmlt = (info.fast.cusp.entry.mlt+info.fast.cusp.exit.mlt)*0.5-$
            (info.polar.cusp.entry.mlt+info.polar.cusp.exit.mlt)*0.5
        info.dr = (info.polar.cusp.entry.dis+info.polar.cusp.exit.dis)*0.5-$
            (info.fast.cusp.entry.dis+info.fast.cusp.exit.dis)*0.5

        ; efluxes.
        info.polar.pflux = info.polar.sb.fh;+info.polar.sb.fl
        info.fast.pflux = info.fast.sb.fh;+info.fast.sb.fl
        info.polar.eflux = info.polar.kei+info.polar.kee+info.polar.pflux
        info.fast.eflux = info.fast.kei+info.fast.kee+info.fast.pflux
        info.ratio_pflux = info.polar.pflux/info.fast.pflux
        info.ratio_eflux = info.polar.eflux/info.fast.eflux
    
        ; Alfven speed.
        if info.polar.cusp.entry.ilat gt info.polar.cusp.exit.ilat then begin
            ilat = info.polar.cusp.exit.ilat
            dis = info.polar.cusp.exit.dis
        endif else begin
            ilat = info.polar.cusp.entry.ilat
            dis = info.polar.cusp.entry.dis
        endelse
        info.polar.va = smodelva(ilat, dis)
    
        if info.fast.cusp.entry.ilat gt info.fast.cusp.exit.ilat then begin
            ilat = info.fast.cusp.exit.ilat
            dis = info.fast.cusp.exit.dis
        endif else begin
            ilat = info.fast.cusp.entry.ilat
            dis = info.fast.cusp.entry.dis
        endelse
        info.fast.va = smodelva(ilat, dis)
    
        store_data, 'scidat', sfmdate(id,'%Y_%m%d_%H UT'), info

    endif


end