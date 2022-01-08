;+
; read info in the log file, use the info to calc Poynting flux, and to 
; plot other data. Save plots into pdf files, save data into tplot file.
;-
pro cusp_polar_save_plot_and_data, eventid, reload = reload, save_data = save_data, no_plot = no_plot, test = test
    
    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    logfile = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'
    loginfo = cusp_read_conjun_list(logfile, event = eventid, /nofast)
    poinfo = loginfo.polar
    if size(loginfo,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif
    
    
    ; directory to save output data and pdfs.
    figdir = rootdir+'/works/cusp/cusp list polar 2-4Re'
    datdir = rootdir+'/data/cusp'

    ; test version do not overwrite files in work dir.
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
; [de,db,pf]idx. the component to be plotted in MAT plots.
; potr. time range for plots.
; potrcusp. cusp entry and exit times.
; pofilt0. time filters from log file.
; pofilts. time filters, uniq and sorted from small to large.
; podelimes. delimiters for low and high freq noises.
; pofacdelim. delimiter for fac.

    ; the wanted component, 0-based counting.
    deidx = 0   ; dEv.
    dbidx = 1   ; dBp.
    pfidx = 2   ; Sb.

    potr = poinfo.plot_time
    potrcusp = poinfo.cusp_time
    if potr[1] lt potr[0] then potr[1]+= 86400d
    if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
    id = loginfo.id

    ; scales controls the MAT decomposition.
    ; filters are determined from the pattern in the MAT specotrogram.
    ; delimiters defines the scales between which poynting flux is counted.
    ; facdelimiter defines FAC scale, do not seperate freq bands for FAC.

    ; filters, min to max, uniq.
    ; if hashighfreqdelim = 1, then take 0-highfreqdelim to be noise.
    ; if haslowfreqdelim = 1, then take lowfreqdelim-maxscale to be background.
    ; if hasfacdelim = 1, then take facdelim-maxscale as fac signal.
    ; By default, take 0 as highfreqdelim and largest filter as lowfreqdelim
    ; and facdelim.

    ; filt0, original filter: f1,f2,f3,f4.
    ; matfilts: 0, f1, f2, f3, f4, maxfilt.
    ; there will be 5 bands.

    maxfilt = 1e31
    ; polar.
    ; filters, sort from small to large.
    pofilt0 = poinfo.filters     ; original filters.
    pofilts = pofilt0
    pofilts = pofilts[sort(pofilts)]
    pofilts = pofilts[uniq(pofilts)]

    ; delimiters.
    podelims = [0,poinfo.noisedelim] ; below and above are noise.
    tmp = where(podelims eq -1, cnt)
    if cnt eq 2 then podelims = minmax(podelims) ; sort if no -1.
    if podelims[0] lt 0 then podelims[0] = 0
    if podelims[1] lt 0 then podelims[1] = potr[1]-potr[0]

    ; fac limit.
    pofacdelim = poinfo.faclim     ; above is fac.
    if pofacdelim eq -1 then pofacdelim = max(pofilts)


    ; podelims. updated to match Polar and FAST phase space.
    autodelims = 0
    if keyword_set(autodelims) then begin
        podelims[0] = potr[1]-potr[0]
        podelims[1] = podelims[1]<pofacdelim
    endif else begin
        podelims = [0,max(pofilts)]
    endelse


; **** plot settings and generate plots to disk.
    posl = [0.1,0.10,0.40,0.90]
    posr = [0.6,0.10,0.90,0.90]
    faclabs = ['v','p','b']
    falabs = 'fa_'+['ilat','mlt','dis']
    polabs = 'po_'+['ilat','mlt','dis']
    rgb = [6,4,2] & red = 6
    ct = 43
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
        conjtr = minmax([potrcusp])
        dt = 0 & conjtr[0]-= dt & conjtr[1]+= dt    ; allow a padding time.
        tmp = sread_omni(conjtr)
        tmp2 = sfmepoch(tmp.epoch,'unix')
        idx = where(tmp2 ge conjtr[0] and tmp2 le conjtr[1])
        info.ae = max(tmp.ae[idx])
        info.dst = min(tmp.symh[idx])
        info.imfbz = minmax(tmp.bgse[*,2])
        info.imfby = minmax(tmp.bgse[*,1])

        info.id = id

        info.polar.cusp.entry.ut = potrcusp[0]
        info.polar.cusp.exit.ut = potrcusp[1]
        info.polar.filters = pofilt0
        info.polar.nfilter = n_elements(pofilt0)
        info.polar.facdelim = pofacdelim
        info.polar.noisedelim = podelims


    ; **** read kinetic energy flux.

        ; po_[ion,ele]_keflux.
        fn = rootdir+'/works/cusp/cusp list polar 2-4Re/'+id+'/'+id+'_ke_special.svg'
        yrange = poinfo.keflux    ; [ele,ion].
        type = strtrim(poinfo.ketype,2)   ; etp,p2p,ptp.
        case type of
            'etp':ex = {epstopdf:1}
            'ptp':ex = {pstopdf:1}
            'p2p':ex = {ps2pdf:1}
        endcase
        polar_read_ke_flux, fn, potr, yrange[0]*[-1,1], yrange[1]*[-1,1], _extra = ex


    ; **** read field data.
    ; [po,fa]_[de,db]_fac, original field.
    ; [po,fa]_tscl. scales in time in second.
    ; [po,fa]_rscl. scales in # of records.

        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        polar_sdt_prep_poynting_flux, fn, e56 = poinfo.e56, tr = potr, eventid = id, $
            cusptr = potrcusp, orootdir = figdir+'/'+id, titpre = 'Event ID: '+id+'    '
        vars = ['po_b','po_b0_spc','po_spc2fac','po_db_spc','po_de_spc','po_mlat']
        store_data, vars, /delete
        get_data, 'po_de_fac', t0
        podr = sdatarate(t0)
        tmp = poinfo.scaleinfo
        potscls = smkgmtrc(tmp[0],tmp[1],tmp[2],'n')
        tmp = round(potscls/podr)
        porscls = double(tmp[uniq(tmp)])
        potscls = porscls*podr

        
    ; **** calc Poynting flux.
    ; [po,fa]_[de,db,pf]_fac_mat. total dE/dB/S.
    ; [po,fa]_[de,db,pf]_fac_mat_f[0,1,2,...]. dE/dB/S in bands.
    ; [po,fa]matfilts. filters used in Poynting flux calc.
    ; [po,fa]nmatband. # of bands from filtering. # of matfilts-1.
    ; [po,fa]matbandids. The labels for mat bands. '_f[0,1,...]'.
    ; [po,fa]nband. # of wanted bands.
        pomatfilts = pofilts
        pomatfilts = [0,pomatfilts,maxfilt]
        ponmatband = n_elements(pomatfilts)-1
        pomatbandids = '_f'+string(indgen(ponmatband),format='(I0)')
        
        dename = 'po_de_fac' & dbname = 'po_db_fac' & pfname = 'po_pf_fac'
        stplot_calc_pflux_mat, dename, dbname, pfname, tscale = potscls, $
            filter = pomatfilts, ids = pomatbandids;, scaleinfo = poinfo.scaleinfo


        remove_wake = 1
        ; remove background for the lowest freq wave band.
        if keyword_set(remove_wake) then begin
            dename = 'po_de_fac_mat'
            dbname = 'po_db_fac_mat'
            pfname = 'po_pf_fac_mat'
            vars = [dename,dbname]
            nvar = n_elements(vars)
            for i = 0, nvar-1 do begin
                idx = ponmatband-1
                tvar = vars[i]+pomatbandids[idx]
                get_data, tvar, t0, dat
                bg = scalcbg(dat)
                store_data, tvar, t0, dat-bg
                tvar = vars[i]+'_nl'
                store_data, tvar, t0, bg, limits = {labels:'low freq bg'+strarr(3)}
            endfor
            ponmatband+= 1
            pomatbandids = [pomatbandids,'_nl']

            ; update poynting flux.
            vars = pomatbandids[ponmatband-2:*]
            for i = 0, n_elements(vars)-1 do $
                stplot_calc_pflux_mat, dename+vars[i], dbname+vars[i], pfname+vars[i]
            store_data, dename+'_f?_comp'+['1','2','3'], /delete
            store_data, dbname+'_f?_comp'+['1','2','3'], /delete
            ; total used field, exclude low and high freq noise.
            vars = [dename,dbname,pfname]+'_mat'
            nvar = n_elements(vars)
            for i = 0, nvar-1 do begin
                tvar = tnames(vars[i]+'_f?')
                stplot_total, tvar, newname = vars[i]
            endfor
        endif


        ponband = ponmatband
        pobandids = pomatbandids
        ponfilt = ponband
        pobandlabels = strarr(ponband)
        for i = 0, ponband-2 do pobandlabels[i] = $
            sgnum2str(sround(pomatfilts[i],error=1e-3))+'-'+sgnum2str(sround(pomatfilts[i+1],error=1e-3))+'s'
        pobandlabels[ponband-1] = 'low freq bg'
        

        ; for all energy fluxes, map and integrate.
        ; flux -> flux_map.
        vars = ['po_'+['ele_keflux','ion_keflux','pf_fac_mat'+['',pobandids]]]
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
            suf = pobandids[i]
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
        

        info.polar.sv.fh = total(pfbands[*,0])
        info.polar.sp.fh = total(pfbands[*,1])
        info.polar.sb.fh = total(pfbands[*,2])
        info.polar.sv.fl = 0d
        info.polar.sp.fl = 0d
        info.polar.sb.fl = 0d
    
        ; ion ratio.
        get_data, pre+'ion_keflux_map', tmp, tmp2
        info.max_kei = abs(min(hem*tmp2[where(tmp ge potrcusp[0] and tmp le potrcusp[1])]))

        vars = pre+'ion_keflux_map_abs'
        store_data, vars, tmp, abs(tmp2)
        cusp_int_eflux, vars, pre+'ilat', tr = potrcusp
        get_data, vars, tmp, tmp2, dat
        info.ratio_ion = info.polar.kei/dat  ; -1 for all upward, 1 for all down.
    

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
        endif

        ; efluxes.
        info.polar.pflux = info.polar.sb.fh;+info.polar.sb.fl
        info.polar.eflux = info.polar.kei+info.polar.kee+info.polar.pflux
    
        ; Alfven speed.
        if info.polar.cusp.entry.ilat gt info.polar.cusp.exit.ilat then begin
            ilat = info.polar.cusp.exit.ilat
            dis = info.polar.cusp.exit.dis
        endif else begin
            ilat = info.polar.cusp.entry.ilat
            dis = info.polar.cusp.entry.dis
        endelse
        info.polar.va = smodelva(ilat, dis)
    
        store_data, 'scidat', sfmdate(id,'%Y_%m%d_%H UT'), info

    endif


; **** save data to disk.
    ; save labels, de_fac, db_fac, pf_fac_mat,
    ; ele_keflux, ion_keflux, pf_fac_mat_para in situ and mapped.
    ; save event info to (1) include other info and (2) info for recovering
    ; quantities are not saved (e.g., mat spectrograms, fields in freq bands).
    if keyword_set(save_data) then begin
        tmp = ['ele_keflux','ion_keflux','pf_fac_mat_para']
        vars = ['ilat','mlt','dis','de_fac','db_fac','pf_fac_mat',tmp,tmp+'_map']
        vars = ['po_'+vars,'scidat']
        store_data, 'event_info', potr[0], loginfo
        vars = [vars,'event_info']
        ofn = datdir+'/'+id+'_all_data'
        tmp = file_dirname(ofn)
        if ~file_test(tmp,/directory) then file_mkdir, tmp
        tplot_save, vars, filename = ofn
    endif
    while !d.window ne -1 do wdelete, !d.window



    if keyword_set(no_plot) then return
    options, ['fa_','po_']+'dis', 'ytitle', 'R'
    

    ; **** plot 0: field preprocessing (several plots).
    ; done by polar_sdt_prep_poynting_flux.
    ofn = figdir+'/'+id+'/'+strmid(id,0,9)+'_polar_*.pdf'
    tmp = file_search(ofn)
    for i = 0, n_elements(tmp)-1 do if tmp[i] ne '' then file_delete, tmp[i]
    ct = 43
    red = 6

    ; **** plot 1: de, db, pf total.
    ; shows the final 3-D de, db, pf, each in common scale.
    ofn = figdir+'/'+id+'/'+id+'_field_and_poynt.pdf'
    sgopen, ofn
    device, decomposed = 0
    loadct2, ct
    erase
    ; polar.
    vars = 'po_'+['de_fac_mat','db_fac_mat','pf_fac_mat']
    options, vars[0], 'labels', 'dE'+faclabs
    options, vars[0], 'ytitle', '(mV/m)'
    options, vars[1], 'labels', 'dB'+faclabs
    options, vars[1], 'ytitle', '(nT)'
    options, vars[2], 'labels', 'S'+faclabs
    options, vars[2], 'ytitle', '(mW/m!U2!N)'
    for j = 0, n_elements(vars)-1 do begin
        stplot_split, vars[j]
        get_data, vars[j], t0, tmp
        idx = where(t0 ge potr[0] and t0 le potr[1])
        ylim, vars[j]+'_comp?', min(tmp[idx,*]), max(tmp[idx,*])
    endfor
    vars = tnames(vars+'_comp?')
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar)
    tplot, vars, var_label = polabs, trange = potr, position = pos, title = 'Polar', /noerase
    timebar, potrcusp, color = red, thick = 2
    store_data, vars, /delete

    xyouts, 0.5, 0.95, /normal, 'Event ID: '+id+'    Total 3-D field and Poynting flux', alignment = 0.5, charsize = 1.25
    store_data, vars, /delete
    sgclose
    

    ; **** plot 2: dE, dB, S in freq bands.
    ; decompose of 3-D fields.
    ofn = figdir+'/'+id+'/'+id+'_field_and_poynt_freq_band.pdf'
    sgopen, ofn, xsize = 9, ysize = 8, /inch
    device, decomposed = 0
    loadct2, ct
    ymax = 0.90
    ymin = 0.10
    nvar = ponfilt+1
    depos = sgcalcpos(nvar, position = [0.1,ymin,0.3,ymax])
    dbpos = sgcalcpos(nvar, position = [0.4,ymin,0.6,ymax])
    pfpos = sgcalcpos(nvar, position = [0.7,ymin,0.9,ymax])
    erase

    ; polar.
    ponewbandids = pobandids
    ponewbandlabels = pobandlabels

    podevars = 'po_de_fac_mat'+ponewbandids
    nvar = n_elements(podevars)
    options, podevars, 'ystyle', 1
    options, podevars, 'ytitle', 'dE'
    options, podevars, 'colors', [6,4,2]
    for i = 0, nvar-1 do begin
        get_data, podevars[i], t0, dat
        options, podevars[i], 'yrange', sg_autolim(dat, ntick = ntick)
        options, podevars[i], 'yticks', ntick
    endfor
    pos = depos
    tplot, podevars, var_label = polabs, trange = potr, position = pos, title = 'Polar dE', /noerase
    timebar, potrcusp, color = red, thick = 2

    podbvars = 'po_db_fac_mat'+ponewbandids
    nvar = n_elements(podbvars)
    options, podbvars, 'labels', faclabs
    options, podbvars, 'ystyle', 1
    options, podbvars, 'ytitle', 'dB'
    options, podbvars, 'colors', [6,4,2]
    for i = 0, nvar-1 do begin
        get_data, podbvars[i], t0, dat
        options, podbvars[i], 'yrange', sg_autolim(dat, ntick = ntick)
        options, podbvars[i], 'yticks', ntick
    endfor
    pos = dbpos
    tplot, podbvars, var_label = polabs, trange = potr, position = pos, title = 'Polar dB', /noerase
    timebar, potrcusp, color = red, thick = 2

    popfvars = 'po_pf_fac_mat'+ponewbandids
    nvar = n_elements(popfvars)
    tmp = info.polar.ebratio
    for i = 0, ponband-1 do begin
        labs = ['dE: '+sgnum2str(tmp[i,1],nsgn=2)+' mV/m', $
            'dB: '+sgnum2str(tmp[i,2],nsgn=2)+' nT', $
            'E/B: '+sgnum2str(tmp[i,0],nsgn=2)+' km/s']
        options, popfvars[i], 'labels', labs
    endfor
    options, popfvars, 'ystyle', 1
    options, popfvars, 'ytitle', 'S'
    options, popfvars, 'colors', [6,4,2]
    for i = 0, nvar-1 do begin
        get_data, popfvars[i], t0, dat
        options, popfvars[i], 'yrange', sg_autolim(dat, ntick = ntick)
        options, popfvars[i], 'yticks', ntick
    endfor
    pos = pfpos
    titl = 'Polar S, v!DA!N = '+sgnum2str(info.polar.va,nsgn=2)+' km/s'
    tplot, popfvars, var_label = polabs, trange = potr, position = pos, title = titl, /noerase
    timebar, potrcusp, color = red, thick = 2

    xyouts, 0.5, 0.95, /normal, 'Event ID: '+id+'    Freq bands of 3-D dE, dB, and Poynting flux', alignment = 0.5, charsize = 1.25
    sgclose
    

    ; **** plot 2': same quantities but in same scale.
    ofn = figdir+'/'+id+'/'+id+'_field_and_poynt_freq_band_scaled.pdf'
    sgopen, ofn, xsize = 9, ysize = 8, /inch
    device, decomposed = 0
    loadct2, ct
    erase
    
    stplot_minmax, podevars, /set, newname = podevars+'_scaled'
    stplot_minmax, podbvars, /set, newname = podbvars+'_scaled'
    stplot_minmax, popfvars, /set, newname = popfvars+'_scaled'
    
    nvar = n_elements(podevars)
    pos = depos
    tplot, podevars+'_scaled', var_label = polabs, trange = potr, position = pos, title = 'Polar dE', /noerase
    timebar, potrcusp, color = red, thick = 2
    pos = dbpos
    tplot, podbvars+'_scaled', var_label = polabs, trange = potr, position = pos, title = 'Polar dB', /noerase
    timebar, potrcusp, color = red, thick = 2
    pos = pfpos
    titl = 'Polar S, v!DA!N = '+sgnum2str(info.polar.va,nsgn=2)+' km/s'
    tplot, popfvars+'_scaled', var_label = polabs, trange = potr, position = pos, title = titl, /noerase
    timebar, potrcusp, color = red, thick = 2
    
    xyouts, 0.5, 0.95, /normal, 'Event ID: '+id+'    Scaled freq bands of 3-D dE, dB, and Poynting flux', alignment = 0.5, charsize = 1.25
    sgclose
    vars = [podevars,podbvars,popfvars]+'_scaled'
    store_data, vars, /delete


    ; **** plot 3: energy fluxes, in situ and mapped.
    ; final in situ and mapped KEi, KEe, Spara.
    ofn = figdir+'/'+id+'/'+id+'_efluxes.pdf'
    sgopen, ofn
    device, decomposed = 0
    loadct2, ct
    erase

    ; polar.
    vars = 'po_'+['ele_keflux','ion_keflux','pf_fac_mat_para']
    tmp = ['KEe','KEi','S!L||!N']+' map!C!C  INT = '
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i]+'_map', tmp1, tmp2, dat
        options, vars[i]+'_map', 'labels', tmp[i]+sgnum2str(dat,ndec=0)+' W/m'
    endfor
    vars = [vars,vars+'_map']
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar)
    options, vars, 'ytitle', '(mW/m!U2!N)'
    options, vars[0], 'labels', 'KEe'
    options, vars[1], 'labels', 'KEi'
    options, vars[2], 'labels', 'S!L||!N'
    tplot, vars, trange = potr, position = pos, /noerase, var_label = polabs, title = 'Polar'
    timebar, potrcusp, color = red, thick = 2

    xyouts, 0.5, 0.95, /normal, 'Event ID: '+id+'    Energy fluxes in situ and mapped', alignment = 0.5, charsize = 1.25
    sgclose
    

    ; **** plot 4: Polar de, db spectrograms, filtered field.

    dematsuf = '_comp'+string(deidx,format='(I0)')
    dbmatsuf = '_comp'+string(dbidx,format='(I0)')


    ; example of decomposition of de and db.
    posde = [0.1,0.10,0.4,0.90]
    posdb = [0.6,0.10,0.9,0.90]
    ofn = figdir+'/'+id+'/'+id+'_mat_spec_polar.pdf'

    sgopen, ofn, xsize = 9, ysize = 6, /inch
    device, decomposed = 0
    loadct2, ct
    erase

    ; dEv.
    dename = 'po_de_fac'

    ; mat spectrogram of the wanted component.
    stplot_index, dename, deidx, newname = dename+dematsuf
    stplot_mat, dename+dematsuf, scale = potscls, zrange = loginfo.polar.de

    vars = dename+'_mat'+ponewbandids
    nvar = n_elements(vars)
    for i = 0, nvar-1 do stplot_index, vars[i], deidx, newname = vars[i]+dematsuf

    vars = dename+'_mat'+ponewbandids
    vars = [dename+dematsuf+['','_mat'],vars+dematsuf]
    options, vars, 'ytitle', '(mV/m)'
    options, vars, 'yticks', 3

    tvar = dename+dematsuf+'_mat'
    options, tvar, 'labels', 'dEv!C  original'
    options, tvar, 'ytitle', 'Period (s)'
    options, tvar, 'ztitle', '(mV/m)'
    options, tvar, 'zticks', 2
    options, dename+dematsuf, 'labels', 'original!C  dEv'
    tvar = dename+'_mat'+ponewbandids+dematsuf
    for i = 0, n_elements(tvar)-1 do options, tvar[i], 'labels', ponewbandlabels[i]

    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posde)

    tplot, vars, var_label = polabs, position = pos, /noerase,$
        title = 'Polar dEv', trange = potr

    ; add filter lines to spec.
    plot, [0,1], poinfo.scaleinfo[0:1], /nodata, /noerase, /ylog, $
        position = pos[*,1], yrange = poinfo.scaleinfo[0:1], ystyle = 5, xstyle = 5
    for i = 0, n_elements(pofilts)-1 do $
        plots, [0,1], pofilts[i]*[1,1], color = 255
    timebar, potrcusp, color = red, thick = 2

    store_data, vars, /delete
    
    ; dBp
    dbname = 'po_db_fac'

    ; mat spectrogram of the wanted component.
    stplot_index, dbname, dbidx, newname = dbname+dbmatsuf
    stplot_mat, dbname+dbmatsuf, scale = potscls, zrange = loginfo.polar.db

    vars = dbname+'_mat'+ponewbandids
    nvar = n_elements(vars)
    for i = 0, nvar-1 do stplot_index, vars[i], dbidx, newname = vars[i]+dbmatsuf

    vars = dbname+'_mat'+ponewbandids
    vars = [dbname+dbmatsuf+['','_mat'],vars+dbmatsuf]
    options, vars, 'ytitle', '(nT)'
    options, vars, 'yticks', 3

    tvar = dbname+dbmatsuf+'_mat'
    options, tvar, 'labels', 'dBp!C  original'
    options, tvar, 'ytitle', 'Period (s)'
    options, tvar, 'ztitle', '(nT)'
    options, tvar, 'zticks', 2
    options, dbname+dbmatsuf, 'labels', 'original!C  dBp'
    tvar = dbname+'_mat'+ponewbandids+dbmatsuf
    for i = 0, n_elements(tvar)-1 do options, tvar[i], 'labels', ponewbandlabels[i]

    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posdb)

    tplot, vars, var_label = polabs, position = pos, /noerase,$
        title = 'Polar dBp', trange = potr

    ; add filter lines to spec.
    plot, [0,1], poinfo.scaleinfo[0:1], /nodata, /noerase, /ylog, $
        position = pos[*,1], yrange = poinfo.scaleinfo[0:1], ystyle = 5, xstyle = 5
    for i = 0, n_elements(pofilts)-1 do $
        plots, [0,1], pofilts[i]*[1,1], color = 255
    timebar, potrcusp, color = red, thick = 2

    store_data, vars, /delete

    xyouts, 0.5, 0.95, /normal, 'Event ID: '+id+'    Polar dEv, dBp decomposition', alignment = 0.5, charsize = 1.25
    sgclose
    

    ; **** plot 4': FAST de, db spectrograms, filtered field.

    ; **** plot 5, FAST ESA.

    ; **** plot 6, footpoint.
    alltr = []
    satnames = ['polar (*,b)']
    satcolors = [sgcolor('blue')]
    satsyms = [2]     ; [*,+].
    pos = [0.15,0.1,0.85,0.8]
    thick = 2
    ofn = figdir+'/'+id+'/'+id+'_footpoint.pdf'
    sgopen, ofn, xsize = 6, ysize = 3, /inch
    sgtruecolor
    ; determine hemisphere.
    get_data, 'po_ilat', t0, poilat
    if interpol(poilat, t0, potrcusp[0]) ge 0 then begin    ; north hem.
        sgset_map, xrange = [90,270], pos = pos, color = sgcolor('black'), $
            ytickv = [50,60,70,80], ytickpos = 225, yticknudge = [-1.4,-0.4], $
            xtickpos = 47
    endif else begin                                        ; south hem.
        sgset_map, xrange = [90,270], pos = pos, color = sgcolor('black'), $
            ytickv = -[50,60,70,80], ytickpos = 225, yticknudge = [-1.4,-0.4], $
            xtickpos = -47, yrange = [-50,-90]
    endelse
    ; use cusp time, expand 2.5 cusp time on both sides.
    ; polar ilat, mlt.
    j = 0
    dt = 600        ; 10 min.
    get_data, 'po_dis', t0, tmp
    if min(interpol(tmp, t0, potrcusp)) le 2.5 then dt = 120
    tmp = potrcusp
    ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
    ttr = ttr-(ttr mod dt) & ttr[1]+= dt
    alltr = [min([ttr,alltr]),max([ttr,alltr])]
    tuts = smkarthm(ttr[0], ttr[1], dt, 'dx')
    get_data, 'po_ilat', t0, tmp
    tilat = interpol(tmp, t0, tuts)
    get_data, 'po_mlt', t0, tmp
    tmlt = interpol(tmp, t0, tuts)*15
    plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j]
    plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j], $
        psym = satsyms[j], symsize = 0.4
    arrow, tmlt[-2], tilat[-2], tmlt[-1], tilat[-1], /data, $
        color = satcolors[j], /solid, thick = thick
    plots, interpol(tmlt, tuts, potrcusp), interpol(tilat, tuts, potrcusp), $
        color = satcolors[j], thick = thick*5
    get_data, 'po_dis', t0, tmp
    tmp = string(interpol(tmp, t0, potrcusp[0]), format='(F3.1)')
    xyouts, 0.6, 0.25, /normal, 'Polar: *, R = '+tmp+' Re', color = satcolors[j]

    xyouts, 0.15, 0.9, /normal, 'Event ID: '+id+'    Footprints of Polar', charsize = 1.25
    sgclose


    ; **** plot 7: omni, polar, fast overview.
    ofn = figdir+'/'+id+'/'+id+'_overview.pdf'
    sgopen, ofn, xsize = 6, ysize = 8, /inch
    device, decomposed = 0
    loadct2, ct
    titl = 'Event ID: '+id+'    Summary plot of OMNI, Polar, and FAST'
    plot_polar_fast_summary, stoepoch(alltr,'unix'), /no_delete, title = titl, /nofast
    sgclose
    
    
; **** output info to console.
;    cusp_gen_excel_form, id, /load, test = test
;    print, 'Polar ion ratio: ', info.ratio_ion
;    print, 'Poynting flux ratio: ', info.ratio_pflux
;    print, 'Energy flux ratio: ', info.ratio_eflux
;    
;    print, 'Polar integrated poynt in freq bands (v,p,b): '
;    print, dindgen(ponband)
;    print, info.polar.sv.fs[0:ponband-1]
;    print, info.polar.sp.fs[0:ponband-1]
;    print, info.polar.sb.fs[0:ponband-1]
;    print, 'FAST integrated poynt in freq bands (v,p,b): '
;    print, dindgen(fanband)
;    print, info.fast.sv.fs[0:fanband-1]
;    print, info.fast.sp.fs[0:fanband-1]
;    print, info.fast.sb.fs[0:fanband-1]
end

ids = cusp_id('polar_south_imf')
for i = 0, n_elements(ids)-1 do cusp_polar_save_plot_and_data, ids[i], /reload, /save_data, /no_plot
end