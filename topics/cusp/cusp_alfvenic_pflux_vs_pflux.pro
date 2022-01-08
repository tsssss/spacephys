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
            
            
            get_data, 'scinfo', tutr, scinfo
            scinfo = (pre0 eq 'po_')? scinfo.polar: scinfo.fast
            cusputr = scinfo.cusp_time
            va = scinfo.va
            va *= sqrt(2)  ; assume H+.
            print, 'Va: ',va,' km/s'


        ;---Calc E/B ratio as a function of periods.
            get_data, pre0+'de_mor_fft_info', 0, einfo
            get_data, pre0+'db_mor_fft_info', 0, binfo

            ; E in mV/m, B in nT, ebratio in km/s.
            ps = einfo.ps
            ebratio = sqrt(einfo.gws/binfo.gws)*1e3
            print, 'E/B ratio: ',minmax(ebratio),' km/s'
            print, 'Period range: ',minmax(ps),' sec'
            
        ;---Calculate the spectrogram of parallel pflux.
            stplot_index, pre0+'pf_fac', 2, newname=pre0+'pfb'
            get_data, pre0+'pfb', uts, pfb
            get_data, pre0+'map_coef', tuts, mapc
            mapc = interpol(mapc, tuts, uts)
            pfb = pfb*mapc
            store_data, pre0+'pfb', uts, pfb
            stplot_mor, pre0+'pfb', tscale = ps
            
            get_data, pre0+'pfb_mor', uts, dat, ps
            idx = where(uts ge cusputr[0] and uts le cusputr[1])
            intpflux = total(dat[idx,*], 1)
            
            
        ;---Divide total power into above and below Alfven speed.
            idx = where(ebratio ge va, cnt)
            pva = (cnt eq 0)? ps[0]: ps[max(idx)] ; the smallest period below va.
            print, 'Period at Va: ',pva, ' sec'
            
            idx = where(ps ge pva, cnt)
            tpfbl = (cnt eq 0)? 0: total(intpflux[idx])
            idx = where(ps lt pva, cnt)
            tpfbh = (cnt eq 0)? 0: total(intpflux[idx])
            
            pfbls = [pfbls,tpfbl]
            pfbhs = [pfbhs,tpfbh]
            print, 'S_FAC, S_Alfven: ', [tpfbl,tpfbh]
        endforeach

        tvar = pre0+'pfb_ebratio'
        if pre0 eq 'po_' then ppfb_ebratio = {lt_va:pfbls, gt_va:pfbhs} $
        else fpfb_ebratio = {lt_va:pfbls, gt_va:pfbhs}
    endforeach
    
    store_data, 'po_pfb_ebratio', 0, ppfb_ebratio
    store_data, 'fa_pfb_ebratio', 0, fpfb_ebratio
    vars = pres+'pfb_ebratio'
    tplot_save, vars, filename=datfn

end
