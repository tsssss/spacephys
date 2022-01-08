;+
; Calc Poynting flux under temporal assumption.
;-

;---Constant.
    re = 6378d & re1 = 1d/re
    h0 = 100d   ; km.
    rad = !dpi/180d
    deg = 180d/!dpi


;---Settings.
    deidx = 0   ; v: north-south.
    dbidx = 1   ; p: east-west.
    pfidx = 2   ; b: parallel.
    faclabs = ['v','p','b']


;---Preparation.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then message, 'rootdir does not exist ...'

    datdir = rootdir+'/data/cusp'
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    datfn = datdir+'/cusp_int_pfb.tplot'
    
    ; if file exists, stop unless told to.
    if file_test(datfn) eq 1 then stop

    ids = cusp_id('south_imf')
    nid = n_elements(ids)


;---Calc Poynting flux for Polar.
    poinfos = []
    fainfos = []
    
    ; integrated pflux, fast, polar, and polar temporal.
    fpfb1s = []
    ppfb1s = []
    ppfb2s = []
    
    foreach tid, ids do begin
        ; read log file for each event.
        loginfo = cusp_read_conjun_list(logfile, event=tid)
        poinfo = loginfo.polar
        fainfo = loginfo.fast
        
        ; load all data.
        ifn1 = datdir+'/'+tid+'_1st_order_data.tplot'
        if file_test(ifn1) eq 0 then cusp_save_1st_order_data, tid, /save_data
        ifn2 = datdir+'/'+tid+'_2nd_order_data.tplot'
        if file_test(ifn2) eq 0 then cusp_save_2nd_order_data, tid, /save_data
        ifn3 = datdir+'/'+tid+'_3rd_order_data.tplot'
        if file_test(ifn3) eq 0 then cusp_save_3rd_order_data, tid, /save_data
        store_data, '*', /delete
        tplot_restore, filename=ifn1
        tplot_restore, filename=ifn2
        tplot_restore, filename=ifn3

        infofn = datdir+'/'+tid+'_scinfo.tplot'
        tplot_restore, filename=infofn
        
        
        get_data, 'scinfo', tutr, scinfo
        pre0 = 'po_'
        poinfo = scinfo.polar
        fainfo = scinfo.fast
        cusputr = poinfo.cusp_time

        ; filter. Change Polar max to Fast's.
        pmin = (poinfo.dr0)*4
        pmax = (fainfo.cusp_time[1]-fainfo.cusp_time[0])

        stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf1_fac', $
            scaleinfo=scaleinfo, filter=[pmin,pmax]

        ; map.
        tvar = pre0+'pf1_fac'
        get_data, pre0+'map_coef', uts, mapc
        get_data, tvar, tuts, dat, limit=lims
        ndim = 3
        tmapc = interpol(mapc, uts, tuts)
        for i=0, ndim-1 do dat[*,i] = dat[*,i]*tmapc
        store_data, tvar+'_map', tuts, dat, limit=lims
        stplot_split, tvar+'_map', newnames=pre0+'pf'+faclabs+'_map'

        ; integrate.
        coef = (re+h0)*rad
        vars = pre0+'pf'+faclabs+'_map'
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
        
        ppfb2s = [ppfb2s,tval[2]]
        ppfb1s = [ppfb1s,poinfo.int_pfb]
        fpfb1s = [fpfb1s,fainfo.int_pfb]
    endforeach

    int_pfb = {po:ppfb1s, fa:fpfb1s, po_temp:ppfb2s}
    tvar = 'int_pfb'
    store_data, tvar, 0, int_pfb
    tplot_save, tvar, filename=datfn
end
