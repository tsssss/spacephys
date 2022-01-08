;+
; de0, db0 need to share uniform time.
; main inputs: tscinfo, tfts, dezr, dbzr.
;-

pro stplot_calc_poynting_flux_mat, de0, db0, pf0, trange = tr, $
    scale = tscs, scaleinfo = tscinfo, filter = tfts, $
    derange = dezr, dbrange = dbzr, deindex = deidx, dbindex = dbidx, $
    veclabs = veclabs, veccols = veccols, fninfo = fninfo
    
    if n_params() ne 3 then message, 'need 3 parameters: de0, db0, pf0 ...'
    
    ; **** preparation.
    pfname = pf0 & pfunit = '(mW/m!U2!N)'
    method = 'mat' & order = 2
    
    dename = tnames(de0) & deunit = '(mV/m)'
    dbname = tnames(db0) & dbunit = '(nT)'
    
    get_data, dbname, t0
    dr = sdatarate(t0) & dr1 = 1d/dr            ; data rate and frequency.
    nrec = n_elements(t0)                       ; number of record.
    
    if n_elements(veclabs) eq 0 then veclabs = ['v','p','b']    ; vector labels.
    if n_elements(veccols) eq 0 then veccols = [6,4,2]          ; vector colors.
    if n_elements(tr) eq 0 then tr = t0[[0,nrec-1]]             ; time range.
    tridx = where(t0 ge tr[0] and t0 le tr[1])  ; index for wanted time.
    
    if n_elements(ytitle) eq 0 then ytitle = pfunit
    if n_elements(labels) ne 3 then labels = $
        stagexist('labels',lim)? lim.labels: veclabs
    if n_elements(colors) ne 3 then colors = $
        stagexist('colors',lim)? lim.colors: [6,4,2]
    
    get_data, dename, t0, def0, limits = lim
    get_data, dbname, t0, dbf0, limits = lim
    deyr = minmax(def0[tridx,*])                ; dE yrange for wanted time.
    dbyr = minmax(dbf0[tridx,*])                ; dB yrange for wanted time.
        
    ; split components. 'de_fac' to 'de_facv', etc.
    stplot_split, dename, newname = dename+veclabs, labels = 'dE'+veclabs
    stplot_split, dbname, newname = dbname+veclabs, labels = 'dB'+veclabs
    
    ; find which component for mat example. choose the largest amplitude.
    if n_elements(deidx) eq 0 then begin
        vars = dename+veclabs
        for i = 0, n_elements(veclabs)-1 do begin
            get_data, vars[i], t0, tmp
            if abs(max(tmp[tridx])) eq abs(max(deyr)) then deidx = i
        endfor
    endif
    if n_elements(dbidx) eq 0 then begin
        vars = dbname+veclabs
        for i = 0, n_elements(veclabs)-1 do begin
            get_data, vars[i], t0, tmp
            if abs(max(tmp[tridx])) eq abs(max(dbyr)) then dbidx = i
        endfor
    endif
    deidxname = dename+veclabs[deidx]
    dbidxname = dbname+veclabs[dbidx]
    
    get_data, deidxname, t0, tmp, limits = lim
    store_data, deidxname+'_smooth', t0, $
        smooth(tmp, 30/dr, /edge_mirror, /nan), limits = lim
    
    ; flag for collecting time scale info, derange, dbrange, filters.
    nudge = 'n'
    
    ; check mat scale. want tscs, mintsc, maxtsc, nsc.
    if n_elements(tscs) ne 0 then begin             ; have time scale.
        mintsc = min(tscs, max = maxtsc)
        nsc = n_elements(tscs)
    endif else begin
        mintsc = 4*dr & maxtsc = 0.5*nrec*dr & nsc = 50
        case n_elements(tscinfo) of
            3: begin                                ; full time scale info.
                mintsc = tscinfo[0] & maxtsc = tscinfo[1]
                nsc = tscinfo[2] & end
            2: begin                                ; min,max scale in time.
                mintsc = tscinfo[0] & maxtsc = tscinfo[1] & end
            1: begin                                ; nscale.
                nsc = tscinfo[2] & end
            0: nudge = 'y'
        endcase
        tscs = smkgmtrc(mintsc,maxtsc,nsc,'n')
    endelse
    
    ; check de and db zrange.
    if n_elements(dezr) eq 0 then begin
        nudge = 'y'
        get_data, deidxname, t0, tmp
        dezr = minmax(tmp)/nsc/3
    endif
    if n_elements(dbzr) eq 0 then begin
        nudge = 'y'
        get_data, dbidxname, t0, tmp
        dbzr = minmax(tmp)/nsc
    endif
    dezr = abs(max(dezr))*[-1,1]
    dbzr = abs(max(dbzr))*[-1,1]
    
    ; do mat for indexed dE and dB component. 'de_facv' to 'de_facv_mat'.
    stplot_mat, deidxname, scale = tscs, zrange = dezr
    stplot_mat, dbidxname, scale = tscs, zrange = dbzr
    mintsc = min(tscs, max = maxtsc) & nsc = n_elements(tscs)
    
    ; **** plot settings.
    sgopen, 0, xsize = 7, ysize = 6, /inch
    sgindexcolor, ct = 43

    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.8
    time_stamp, /off
    
    posde = [0.1,0.15,0.45,0.95]
    posdb = [0.6,0.15,0.95,0.95]

    ; **** collect tscs, mintsc, maxtsc, nsc, dezr, dbzr, fts, nft, ftids.
    if nudge eq 'y' then begin
        erase
        vars = deidxname+['_smooth','_mat'] & nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posde)
        tplot, vars, position = pos, tr = tr, /noerase, title = ''
        vars = dbidxname+['','_mat'] & nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posdb)
        tplot, vars, position = pos, tr = tr, /noerase, title = ''
        while nudge eq 'y' do begin
            printf, -1, 'Current scale info: '
            printf, -1, [mintsc,maxtsc,nsc]
            tmp = ''
            read, tmp, prompt = 'Time scale info? (min, max, n): '
            if tmp ne '' then begin
                tscinfo = double(strsplit(tmp,', ',/extract))
                case n_elements(tscinfo) of
                    3: begin                            ; full time scale info.
                        mintsc = tscinfo[0] & maxtsc = tscinfo[1]
                        nsc = tscinfo[2] & end
                    2: begin                            ; min,max scale in time.
                        mintsc = tscinfo[0] & maxtsc = tscinfo[1] & end
                    1: if tscinfo[0] ne 0 then nsc = tscinfo[0]
                    else: ; do nothing.
                endcase
                tscs = smkgmtrc(mintsc,maxtsc,nsc,'n')
                stplot_mat, deidxname, scale = tscs, zrange = dezr
                stplot_mat, dbidxname, scale = tscs, zrange = dbzr
                mintsc = min(tscs, max = maxtsc) & nsc = n_elements(tscs)
            endif
            ; read zrange.
            printf, -1, 'Current dE zrange: '
            printf, -1, dezr
            printf, -1, 'Current dB zrange: '
            printf, -1, dbzr
            tmp = ''
            read, tmp, prompt = 'dE zrange? (min,max): '
            if tmp ne '' then begin
                dezr = double(strsplit(tmp,', ',/extract))
                if n_elements(dezr) eq 1 then dezr = max(dezr,/absolute)*[-1d,1]
            endif
            tmp = ''
            read, tmp, prompt = 'dB zrange? (min,max): '
            if tmp ne '' then begin
                dbzr = double(strsplit(tmp,', ',/extract))
                if n_elements(dbzr) eq 1 then dbzr = max(dbzr,/absolute)*[-1d,1]
            endif
            options, deidxname, 'zrange', dezr
            options, dbidxname, 'zrange', dbzr
            
            erase
            vars = deidxname+['_smooth','_mat'] & nvar = n_elements(vars)
            pos = sgcalcpos(nvar, position=posde)
            tplot, vars, position = pos, tr = tr, /noerase, title = ''
            vars = dbidxname+['','_mat'] & nvar = n_elements(vars)
            pos = sgcalcpos(nvar, position=posdb)
            tplot, vars, position = pos, tr = tr, /noerase, title = ''
            
            printf, -1, 'Current scale info: '
            printf, -1, [mintsc,maxtsc,nsc]
            printf, -1, 'Current dE zrange: '
            printf, -1, dezr
            printf, -1, 'Current dB zrange: '
            printf, -1, dbzr
            
            read, nudge, prompt = 'Tune up scale and zrange? (y/n): '
            if nudge eq '' then nudge = 'y'
        endwhile
    endif
    

    ; check filters. want tfts, nft, ftids.
    if n_elements(tfts) le 1 then begin
        nudge = 'y' & tfts = [mintsc,maxtsc]
    endif
    tfts = stofilter(tfts) & nft = n_elements(tfts)/2
    ftids = string(indgen(nft)+1,format='(I0)')
    
    ; do filter for indexed dE and dB component. 'de_facv' to 'de_facvf1'.
    for i = 0, nft-1 do begin
        stplot_filter, deidxname+'_mat', filter = tfts[i,*], ifilter = ftids[i]
        stplot_filter, dbidxname+'_mat', filter = tfts[i,*], ifilter = ftids[i]
    endfor
    
    if nudge eq 'y' then begin
        erase
        vars = deidxname+['_mat','_mat'+ftids] & nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posde)
        tplot, vars, position = pos, tr = tr, /noerase, title = ''
        vars = dbidxname+['_mat','_mat'+ftids] & nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posdb)
        tplot, vars, position = pos, tr = tr, /noerase, title = ''

        while nudge eq 'y' do begin
            ctime, tmp, tfts, /exact
            if n_elements(tfts) ge 2 then begin
                store_data, deidxname+'f'+ftids, /delete
                store_data, dbidxname+'f'+ftids, /delete
                tfts = stofilter(tfts) & nft = n_elements(tfts)/2
                ftids = string(indgen(nft)+1,format='(I0)')

                for i = 0, nft-1 do begin
                    stplot_filter, deidxname+'_mat', filter = tfts[i,*], ifilter = ftids[i]
                    stplot_filter, dbidxname+'_mat', filter = tfts[i,*], ifilter = ftids[i]
                endfor

                erase
                vars = deidxname+['_mat','_mat'+ftids] & nvar = n_elements(vars)
                pos = sgcalcpos(nvar, position=posde)
                tplot, vars, position = pos, tr = tr, /noerase, title = ''
                vars = dbidxname+['_mat','_mat'+ftids] & nvar = n_elements(vars)
                pos = sgcalcpos(nvar, position=posdb)
                tplot, vars, position = pos, tr = tr, /noerase, title = ''
            endif

            printf, -1, 'Current filter: '
            printf, -1, [tfts[*,0],tfts[nft*2-1]]
            read, nudge, prompt = 'Tune up filter? (y/n): '
            if nudge eq '' then nudge = 'y'
        endwhile
    endif
    
    ; **** should have tscs, mintsc, maxtsc, nsc, dezr, dbzr, tfts, nft, ftids.
    ; do mat for each component of dE and dB, eg, de_facv -> de_facv_mat.
    defnames = dename+veclabs
    dbfnames = dbname+veclabs
    ndim = n_elements(veclabs)
    fnames = [defnames, dbfnames]
    for i = 0, 2*ndim-1 do stplot_mat, fnames[i], scale = tscs
    sfts = strarr(nft)
    for i = 0, nft-1 do $
        sfts[i] = snum2str(tfts[i,0],/short)+'-'+snum2str(tfts[i,1],/short)+'s'
    
    ; do filter, eg, de_facv_mat -> de_facv_mat1, de_facv_mat2.
    for i = 0, ndim-1 do for j = 0, nft-1 do stplot_filter, $
        defnames[i]+'_mat', 'mat', filter = tfts[j,*], ifilter = ftids[j], $
        limits = {ytitle:deunit, labels:'dE'+veclabs[i]+'!C  '+sfts[j], colors:veccols[i]}
    for i = 0, ndim-1 do for j = 0, nft-1 do stplot_filter, $
        dbfnames[i]+'_mat', 'mat', filter = tfts[j,*], ifilter = ftids[j], $
        limits = {ytitle:dbunit, labels:'dB'+veclabs[i]+'!C  '+sfts[j], colors:veccols[i]}
    
    ; combine fields, eg, de_fac[v,p,b]_mat1 -> de_fac_mat1. order matters!!!
    for i = 0, nft-1 do begin
        tlabs = veclabs & tlabs[n_elements(tlabs)-1] += '!C  '+sfts[i]
        stplot_merge, defnames+'_mat'+ftids[i], newname = $
            dename+'_mat'+ftids[i], limits = {colors:veccols, labels:'dE'+tlabs, ytitle:deunit}
        stplot_merge, dbfnames+'_mat'+ftids[i], newname = $
            dbname+'_mat'+ftids[i], limits = {colors:veccols, labels:'dB'+tlabs, ytitle:dbunit}
    endfor
        
    ; calc poynting flux for each bands.
    for i = 0, nft-1 do begin
        get_data, dename+'_mat'+ftids[i], t0, de
        get_data, dbname+'_mat'+ftids[i], t0, db
        store_data, pfname+'_mat'+ftids[i], t0, spoynt(de,db), $
            limits = {ytitle:pfunit,colors:veccols, labels:'S'+veclabs+['','','!C  '+sfts[i]]}
    endfor
    
    ; sum up db, de, pf.
    tits = [deunit,dbunit,pfunit]
    labs = ['dE','dB','S']
    vars = [dename,dbname,pfname]
    for i = 0, n_elements(vars)-1 do begin
        stplot_total, vars[i]+'_mat'+ftids, newname = vars[i]+'_mat', $
            limits = {ytitle:tits[i], labels:labs[i]+veclabs, colors:veccols}
        get_data, vars[i], tmp, dat
        get_data, vars[i]+'_mat', tmp, dat1
        store_data, vars[i]+'_matd', tmp, dat-dat1, $
            limits = {ytitle:tits[i], labels:labs[i]+veclabs, colors:veccols}
    endfor
    store_data, pfname+'_matd', /delete
    
    ; **** make some plots.
    if ~n_elements(make_plot) then return
    tstr = time_string(t0[0],tformat='YYYY_MMDD')
    fninfo = n_elements(fninfo)? '_'+fninfo: ''
    fn = shomedir()+'/'+tstr+fninfo+'_mat_spectrogram.eps'
    sgopen, fn
    erase   ; mat and filters and filtered E/B fields.
    vars = pfname+'_mat_filter'
    store_data, vars, tr, [tfts[*,0],tfts[nft*2-1]]##[1,1], limits = {colors:intarr(nft+1)-1}
    store_data, deidxname+'_comb', data = [deidxname+'_mat',vars], limits = {yrange:[mintsc,maxtsc],zrange:dezr}
    store_data, dbidxname+'_comb', data = [dbidxname+'_mat',vars], limits = {yrange:[mintsc,maxtsc],zrange:dbzr}
    vars = deidxname+['_smooth','_comb','_mat'+ftids] & nvar = n_elements(vars)
    pos = transpose(sgcalcpos(nvar, position=posde))
    stplot_uniyr, deidxname+'_mat'+ftids, trange = tr
    tplot, vars, position = pos, tr = tr, /noerase, title = deidxname
    vars = dbidxname+['','_comb','_mat'+ftids] & nvar = n_elements(vars)
    pos = transpose(sgcalcpos(nvar, position=posdb))
    stplot_uniyr, dbidxname+'_mat'+ftids, trange = tr
    tplot, vars, position = pos, tr = tr, /noerase, title = dbidxname
    sgclose
;    sps2pdf, fn, /rm
    
    fn = shomedir()+'/'+tstr+fninfo+'_field_poynt_freq_band.eps'
    sgopen, fn, xsize = 6, ysize = 8, /inch
    vars = ['']
    for i = 0, nft-1 do vars = [vars,[dename,dbname,pfname]+'_mat'+ftids[i]]
    vars = vars[1:*]
    options, vars, 'yticks', 3
    tplot, vars, tr = tr, /noerase, title = tstr+' Field and Poynting flux in freq bands'
    sgclose
;    sps2pdf, fn, /rm
    
    fn = shomedir()+'/'+tstr+fninfo+'_field_poynt.eps'
    sgopen, fn, xsize = 6, ysize = 8, /inch
    vars = dename+veclabs+'_mat' & stplot_split, dename+'_mat', newname = vars, colors = veccols & stplot_uniyr, vars, trange = tr
    vars = dbname+veclabs+'_mat' & stplot_split, dbname+'_mat', newname = vars, colors = veccols & stplot_uniyr, vars, trange = tr
    vars = pfname+veclabs+'_mat' & stplot_split, pfname+'_mat', newname = vars, colors = veccols & stplot_uniyr, vars, trange = tr
    vars = [pfname+veclabs,dename+veclabs,dbname+veclabs]+'_mat'
    tplot, vars, tr = tr, /noerase, title = tstr+' Field and Poynting flux'
    sgclose
;    sps2pdf, fn, /rm
    
    while !d.window ne -1 do wdelete, !d.window
end

pofn = sdiskdir('Research')+'/data/cusp/po_sdt_fld_1998_1001_03.sdt'
potr = time_double(['1998-10-01/02:00','1998-10-01/04:30'])
;polar_sdt_prep_poynting_flux, pofn
stplot_calc_poynting_flux_mat, 'po_de_fac', 'po_db_fac', 'po_pf_fac', fninfo = 'polar', $
    filter = [2820d,860,236,52,13], scaleinfo = [6,3500,60], derange = 0.5, dbrange = 1.2
end