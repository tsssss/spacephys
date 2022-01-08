;+
; Read and save third order data, including
;   filtered S, mapped and integrated KEi, KEe, S
;   finished scinfo
;-


pro cusp_save_3rd_order_data, eventid, test=test, $
    no_plot=noplot, save_data=save_data
    
    
    if n_elements(eventid) eq 0 then message, 'no event id ...'
    
;---constant.
    re = 6378d & re1 = 1d/re
    h0 = 100d   ; km.
    rad = !dpi/180d
    deg = 180d/!dpi
    
    
;---settings.
    ; rootdir to save the data and plots.
    rootdir = sdiskdir('GoogleDrive')+'/My Drive/works'
    if file_test(rootdir) eq 0 then rootdir = shomedir()

    ; dir to save plots and data.
    figdir = rootdir+'/works/cusp/cusp list conjun'
    datdir = rootdir+'/data/cusp'
    if ~file_test(figdir,/directory) then file_mkdir, figdir
    if ~file_test(datdir,/directory) then file_mkdir, datdir

    ; to prevent overwriting useful data without confirm.
    if keyword_set(test) then begin
        figdir = shomedir()+'/cusp'
        datdir = shomedir()+'/cusp/data'
    endif


    deidx = 0   ; v: north-south.
    dbidx = 1   ; p: east-west.
    pfidx = 2   ; b: parallel.

    sats = ['po','fa']
    
    ; read basic info for the conjunc event.

    
    
    posu = [0,0.4,1,1]
    posd = [0,0,1,0.4]
    
    labfac = ['n-s','e-w','b']
    
    ; plot.
    zchsz = 0.8d
    chsz0 = 1d
    !p.font = 1
    
    ; color.
    ct = 43
    red = 6
    blue = 2
    white = 255
    rgb = [6,4,2]


;---plot settings and generate plots to disk.
    posl = [0.1,0.10,0.40,0.90]
    posr = [0.6,0.10,0.90,0.90]
    faclabs = ['v','p','b']
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'zcharsize', zchsz
    time_stamp, /off
    
    

;---load 1st and 2nd order data.
    ifn1 = datdir+'/'+eventid+'_1st_order_data.tplot'
    if file_test(ifn1) eq 0 then begin
        cusp_save_1st_order_data, eventid, /save_data
    endif
    ifn2 = datdir+'/'+eventid+'_2nd_order_data.tplot'
    if file_test(ifn2) eq 0 then begin
        cusp_save_2nd_order_data, eventid, /save_data
    endif
    store_data, '*', /delete
    tplot_restore, filename=ifn1
    tplot_restore, filename=ifn2
    
    infofn = datdir+'/'+eventid+'_scinfo.tplot'
    tplot_restore, filename=infofn




;---calculate data and save to tplot.
; [po,fa]_[pf1,de1,db1]_fac. filtered field.
; [po,fa]_[pfb,,kei,kee]_map. mapped & filtered data.
; [po,fa]_[pfb_map_int]. This is the mapped and integrated parallel pflux as a function of periods.
; [po,fa]_ebratio. This is the E/B ratio as a function of periods.

    get_data, 'scinfo', tutr, scinfo
    foreach tsat, sats do begin
        pre0 = tsat+'_'
        tinfo = (tsat eq 'po')? scinfo.polar: scinfo.fast
        cusputr = tinfo.cusp_time

        ; filter.
        pmin = (tinfo.dr0)*4
        pmax = (tinfo.cusp_time[1]-tinfo.cusp_time[0])
        stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf1_fac', $
            scaleinfo=scaleinfo, filter=[pmin,pmax]


        ; map.
        vars = pre0+['ion_keflux','ele_keflux','pf1_fac']
        get_data, pre0+'map_coef', uts, mapc
        foreach tvar, vars do begin
            get_data, tvar, tuts, dat, limit=lims
            ndim = (tvar eq pre0+'pf1_fac')? 3: 1
            tmapc = interpol(mapc, uts, tuts)
            for i=0, ndim-1 do dat[*,i] = dat[*,i]*tmapc
            store_data, tvar+'_map', tuts, dat, limit=lims
        endforeach

        stplot_split, pre0+'pf1_fac_map', newnames=pre0+'pf'+faclabs+'_map'
        stplot_renew, pre0+'ion_keflux_map', newname=pre0+'kei_map'
        stplot_renew, pre0+'ele_keflux_map', newname=pre0+'kee_map'

        ; integrate.
        coef = (re+h0)*rad
        vars = pre0+['kei','kee','pf'+faclabs]+'_map'
        nvar = n_elements(vars)
        tval = dblarr(nvar)
        get_data, pre0+'ilat', uts, ilats
        for i=0, nvar-1 do begin
            tvar = vars[i]
            get_data, tvar, tuts, tdat
            idx = where(tuts ge cusputr[0] and tuts le cusputr[1], cnt)
            tuts = tuts[idx]
            tdat = tdat[idx]
            tilats = interpol(ilats, uts, tuts)
            v1 = 0.5*(tdat[1:cnt-1]+tdat[0:cnt-2])
            v2 = abs(tilats[1:cnt-1]-tilats[0:cnt-2])
            tval[i] = total(v1*v2*coef,/nan)
        endfor
        
        tinfo.int_kei = tval[0]
        tinfo.int_kee = tval[1]
        tinfo.int_pfb = tval[4]
        tinfo.int_pfs = tval[2:4]
        tinfo.filters = [pmin,pmax]
        
        ; ebratio and parallel pflux integrated over time.
        ; pre_pf1_fac_mor_spec_3 saves the pflux which directly sum to pf1_fac.
        get_data, pre0+'pf1_fac_mor_spec_3', uts, dat, ps
        ; map.
        get_data, pre0+'map_coef', tuts, mapc
        mapc = interpol(mapc, tuts, uts)
        np = n_elements(ps)
        for i=0, np-1 do dat[*,i] *= mapc
        ; integrate over time.
        coef = (re+h0)*rad
        get_data, pre0+'ilat', tuts, ilats
        ilats = interpol(ilats, tuts, uts)
        idx = where(uts ge cusputr[0] and uts le cusputr[1], cnt)
        dat = dat[idx,*]
        ilats = ilats[idx]
        v1 = 0.5*(dat[1:cnt-1,*]+dat[0:cnt-2,*])
        v2 = abs(ilats[1:cnt-1]-ilats[0:cnt-2])
        dat = fltarr(np)
        for i=0, np-1 do dat[i] = total(v1[*,i]*v2*coef,/nan)
        store_data, pre0+'pfb_spec', dat, ps
        
        ; ebratio.
        get_data, pre0+'de_mor_fft_info', 0, einfo
        get_data, pre0+'db_mor_fft_info', 0, binfo
        ebratio = sqrt(einfo.gws/binfo.gws)*1e3
        store_data, pre0+'ebratio_spec', ebratio, ps, limits={unit:'km/s'}       
        
        
    ;---ion and electron ratio.
        get_data, pre0+'ilat', uts, ilats
        vars = pre0+['kei','kee','pfb']+'_map'
        nvar = n_elements(vars)
        tval = dblarr(nvar)
        get_data, pre0+'ilat', uts, ilats
        for i=0, nvar-1 do begin
            tvar = vars[i]
            get_data, tvar, tuts, tdat
            idx = where(tuts ge cusputr[0] and tuts le cusputr[1], cnt)
            tuts = tuts[idx]
            tdat = abs(tdat[idx])
            tilats = interpol(ilats, uts, tuts)
            v1 = 0.5*(tdat[1:cnt-1]+tdat[0:cnt-2])
            v2 = abs(tilats[1:cnt-1]-tilats[0:cnt-2])
            tval[i] = total(v1*v2*coef,/nan)
        endfor
        tinfo.r_kei = tinfo.int_kei/tval[0]
        tinfo.r_kee = tinfo.int_kee/tval[1]
        tinfo.r_pfb = tinfo.int_pfb/tval[2]

        
    ;---save results.
        if tsat eq 'po' then begin
            scinfo.polar = tinfo
        endif else begin
            scinfo.fast = tinfo
        endelse
    endforeach
    
    
    
    
    store_data, 'scinfo', tutr, scinfo
    tplot_save, 'scinfo', filename = infofn

    
    
;---save data to disk.
    if keyword_set(save_data) then begin
        vars = [['pf1','de1','db1']+'_fac',['pfb','kei','kee']+'_map',['pfb','ebratio']+'_spec']
        vars = ['po_'+vars,'fa_'+vars]
        ofn = datdir+'/'+eventid+'_3rd_order_data.tplot'
        tplot_save, vars, filename = ofn
    endif
    
    
end

;id = '1998_1001_03'
;cusp_save_1st_order_data, id, /save_data
;cusp_save_2nd_order_data, id, /save_data
;cusp_save_3rd_order_data, id, /save_data
;end

ids = cusp_id('all')
ids = ['1998_0925_05']
foreach id, ids do begin
    print, '----processing '+id+' ...'
    cusp_save_1st_order_data, id, /save_data
    cusp_save_2nd_order_data, id, /save_data
    cusp_save_3rd_order_data, id, /save_data
endforeach

end
