;+
; Use Local Alfven speed as a threshold to calculate the Poynting flux power below and above the threshold.
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
    rootdir = sdiskdir('GoogleDrive')+'/My Drive/works'
    if file_test(rootdir) eq 0 then message, 'rootdir does not exist ...'

    datdir = rootdir+'/data/cusp'
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    datfn = datdir+'/cusp_int_pfb_ebratio.tplot'

    ids = cusp_id('south_imf')
    nid = n_elements(ids)
    
    int_kei = []
    hem = []


;---Calc Poynting flux for Polar.
    
    ; integrated pflux, for E/B ratio above and below Va.
    pres = ['po','fa']+'_'
    foreach pre0, pres do begin
        pfbls = []  ; E/B ratio < Va.
        pfbhs = []  ; E/B ratio > Va.
        foreach tid, ids do begin
            printf, -1, 'Processing '+tid+' ...'
            ; read log file for each event.
            loginfo = cusp_read_conjun_list(logfile, event=tid)
            
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
            
            
        ;---Get Alfven speed.
            get_data, 'scinfo', tutr, scinfo
            scinfo = (pre0 eq 'po_')? scinfo.polar: scinfo.fast
            cusputr = scinfo.cusp_time
            va = scinfo.va  ; Alfven speed for all H+.
            print, 'Va: ',va,' km/s'


        ;---Get E/B ratio as a function of periods, and the spectrum of parallel pflux.
            get_data, pre0+'ebratio_spec', ebratio, ps
            get_data, pre0+'pfb_spec', pfb_spec, ps

            filters = scinfo.filters
            idx = where(ps ge filters[0] and ps le filters[1])
            ps = ps[idx]
            pfb_spec = pfb_spec[idx]
            ebratio = ebratio[idx]

            print, 'E/B ratio: ',minmax(ebratio),' km/s'
            print, 'Period range: ',minmax(ps),' sec'

            
        ;---Divide total power into above and below Alfven speed.
            if pre0 eq 'po_' then begin
                int_kei = [int_kei, scinfo.int_kei]
                hem = [hem, scinfo.ilat0/abs(scinfo.ilat0)]
            endif
            int_pfb = scinfo.int_pfb
            
            idx = where(ebratio ge va, cnt)
            pva = (cnt eq 0)? ps[0]: ps[max(idx)] ; the smallest period below va.
            print, 'Period at Va: ',pva, ' sec'
            
            ; low frequency, or non-Alfvenic.
            idx = where(ps ge pva, cnt)
            tpfbl = (cnt eq 0)? 0: total(pfb_spec[idx])
            ; high frequency, or Alfvenic.
            idx = where(ps lt pva, cnt)
            tpfbh = (cnt eq 0)? 0: total(pfb_spec[idx])
            
            pfbls = [pfbls,tpfbl]
            pfbhs = [pfbhs,tpfbh]
            print, 'S total: ', int_pfb
            print, 'S Alfvenic + non-Alfvenic: ', tpfbl+tpfbh
            print, 'S non-Alfvenic and S Alfvenic): ', [tpfbl,tpfbh]
        endforeach

        tvar = pre0+'pfb_ebratio'
        if pre0 eq 'po_' then ppfb_ebratio = {non_alf:pfbls, alfven:pfbhs} $
        else fpfb_ebratio = {non_alf:pfbls, alfven:pfbhs}
    endforeach
    
    store_data, 'po_pfb_ebratio', 0, ppfb_ebratio
    store_data, 'fa_pfb_ebratio', 0, fpfb_ebratio
    vars = pres+'pfb_ebratio'
    tplot_save, vars, filename=datfn

end
