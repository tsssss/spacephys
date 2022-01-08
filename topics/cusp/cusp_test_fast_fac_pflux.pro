;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro cusp_test_fast_fac_pflux, type, ids = ids

    sgindexcolor, 43
    blue = 2
    red = 6
    
    get_data, 'fa_sb_fac_orig', tmp, fa_sb_fac_orig
    if n_elements(fa_sb_fac_orig) gt 1 then goto, do_plot
    
    ; assume type must be there if no id is set.
    if n_elements(ids) eq 0 then ids = cusp_id(type) else type = 'usrdef'

    nid = n_elements(ids)
    fa_sb_fac_orig = dblarr(nid)
    fa_sb_fac_nobg = dblarr(nid)
    
    maxfilt = 1e31

    ; read fac poynting flux.
    ; calc poynting flux using background remove algorithm.
    for k = 0, nid-1 do begin
        
        id = ids[k]
        
        ; read fac poynting flux.
        rootdir = shomedir()+'/Google Drive/works/data/cusp'
        fn = rootdir+'/'+id+'_all_data.tplot'
        tplot_restore, filename = fn
        get_data, 'scidat', tmp, info

        fa_sb_fac_orig[k] = info.fast.sb.fl


    ; **** calc fac poynting flux using background remove algorithm.

        rootdir = shomedir()+'/Google Drive/works'
        logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
        loginfo = cusp_read_conjun_list(logfile, event = id)
        fainfo = loginfo.fast

        fatr = loginfo.fast.plot_time
        fatrcusp = loginfo.fast.cusp_time
        if fatr[1] lt fatr[0] then fatr[1]+= 86400d
        if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d

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

        tvar = 'fa_db_fac'
        get_data, tvar, t0, dat
        bg = scalcbg(dat)
        store_data, tvar, t0, dat-bg
        tvar = 'fa_de_fac'
        get_data, tvar, t0, dat
        bg = scalcbg(dat)
        store_data, tvar, t0, dat-bg

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


        ; get the fac poynting flux, map and integrate.
        tvar = pfname+'_mat_'+fabandids[fanband-1]
        smap2iono, tvar, 'fa_dis', newname = tvar+'_map'
        cusp_int_eflux, tvar+'_map', 'fa_ilat', tr = fatrcusp

        get_data, tvar+'_map', t0, dat, val
        fa_sb_fac_nobg[k] = val[2]

        hem = info.fast.cusp.entry.ilat & hem = hem/abs(hem)
        fa_sb_fac_nobg[k]*= hem
        
        store_data, '*', /delete

    endfor

    store_data, 'fa_sb_fac_orig', 0, fa_sb_fac_orig
    store_data, 'fa_sb_fac_nobg', 0, fa_sb_fac_nobg

    do_plot:
    
    get_data, 'fa_sb_fac_orig', 0, fa_sb_fac_orig
    get_data, 'fa_sb_fac_nobg', 0, fa_sb_fac_nobg

    ofn = shomedir()+'/cusp/fig/'+'cusp_test_fig_fast_fac_pflux'+type+'.eps'
    sgpsopen, ofn, xsize = 4, ysize = 4, /inch
    sgtruecolor
    erase


    tpos = [0.15,0.15,0.85,0.85]
    xx = fa_sb_fac_orig
    yy = fa_sb_fac_nobg
    xr = [1e-1,1e5]
    yr = [1e-1,1e5]

    ex = {noerase:1, position:tpos, normal:1, xrange:xr, yrange:yr, xlog:1, ylog:1, ytitle:'FAST S!Drm bg', xtitle:'FAST S!Dorig'}
    plot, xx, yy, _extra = ex, xstyle = 1, ystyle = 1, nodata = 1, title = 'FAST FAC S, removed background vs original'
    plots, xr, yr, linestyle = 1
    tmp = findgen(11)*2*!dpi/10
    txs = cos(tmp)
    tys = sin(tmp)
    usersym, txs, tys, color = 0
    plots, xx, yy, psym = 8, symsize = 0.7
    
    sgpsclose, /pdf

end

cusp_test_fast_fac_pflux, 'default'
end
