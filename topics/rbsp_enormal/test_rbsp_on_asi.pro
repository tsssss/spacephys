;+
; plot rbsp-a and -b footprint on asf. 
;-

    utr = time_double(['2012-11-14/02:00','2012-11-14/05:00'])

    ; load rbsp footprint.
    probes = ['a','b']
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        tvar = pre0+'vsc'
        dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'vsvy')
        
        tmp = sfmepoch(dat.epoch, 'unix')
        dr = sdatarate(tmp[sort(tmp)])
        uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
        dat = sinterpol(dat.vsvy, tmp, uts)
        
        cp0 = rbsp_efw_get_cal_params(utr[0])
        case strlowcase(tprobe) of
            'a': cp = cp0.a
            'b': cp = cp0.b
            else: dprint, 'Invalid probe name. Calibration aborted.'
        endcase
        gain = cp.adc_gain_vdc
        offset = cp.adc_offset_vdc
        for i = 0, 5 do dat[*,i] =  (dat[*,i]-offset[i])*gain[i]
        vsc = mean(dat[*,0:3], dimension = 2)
        store_data, tvar, uts, vsc, limits = {ytitle:'Vsc!C(V)', labels:'Vsc'}
        plot, vsc, /nodata & yr = !y.crange & dyr = abs(yr[1]-yr[0])
        if dyr le 0.1 then yr = mean(yr)+[-1,1]*dyr ; 0.2 V range.
        options, tvar, 'yrange', yr
        
        ; load position.
        dat = sread_rbsp_efw_l3(utr, probes = tprobe)
        tmp = sfmepoch(dat.epoch, 'unix')
        dr = sdatarate(tmp[sort(tmp)])
        uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
        tdat = sinterpol(dat.pos_gse, tmp, uts)/6378d
        tvar = pre0+'pos_gse'
        pre1 = 'RBSP-'+strupcase(tprobe)
        store_data, tvar+'x', uts, tdat[*,0], limits = {ytitle:pre1+' X GSE (Re)'}
        store_data, tvar+'y', uts, tdat[*,1], limits = {ytitle:pre1+' Y GSE (Re)'}
        store_data, tvar+'z', uts, tdat[*,2], limits = {ytitle:pre1+' Z GSE (Re)'}
    endforeach
    
    get_data, 'rbspa_vsc', uts, vsca
    get_data, 'rbspb_vsc', tmp, vscb & vscb = sinterpol(vscb, tmp, uts)
    
    store_data, 'rb_vsc', uts, [[vsca],[vscb]], limits = $
        {ytitle:'Vsc!C(V)', labels:['RBSP-A','RBSP-B'], colors:[6,2]}
    
    vars = 'rb_vsc'
    labs = ['rbspa_pos_gse'+['z','y','x'],'rbspb_pos_gse'+['z','y','x']]
    
    
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero', 1
    tplot_options, 'no_interp', 1
    
    ofn = shomedir()+'/ab_conj_vsc.pdf'
    sgopen, ofn, xsize = 10, ysize = 6, /inch
    device, decomposed = 0
    loadct2, 43
    
    tplot, vars, var_label = labs, trange = utr, position = [0.2,0.3,0.9,0.9]

    sgclose
stop

;    ; load asi.
;    minlat = 50
;    asipos = [0,0,1,1]
;    sites = ['kian','gako','fykn','whit','atha','tpas','gill','pina','gbay']
;    ; 'inuv','mgcr','pgeo','kapu','fsmi','rank','kuuj'
;
;    uts = smkarthm(utr[0],utr[1],3,'dx')
;    rootdir = shomedir()+'/psbl_de_asi'
;    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
;
;    nrec = n_elements(uts)
;    for i = 0, nrec-1 do begin
;        tut = uts[i]
;        ofn = rootdir+'/thg_asf_mosaic_'+time_string(tut,tformat='YYYY_MMDD_hhmm_ss')+'.png'
;        asimlt = sread_thg_mlt(uts[i], sites, minlat = minlat, type = 'ast', /nomoon)
;
;        timg = asimlt.mltimg
;        sgopen, ofn, xsize = 6, ysize = 6, /inch
;        sgtv, timg, position = asipos, /nan
;        sgset_map, pos = asipos, $
;            color = sgcolor('white'), charsize = 1, majorgrid = 1, minorgrid = 1
;        xyouts, asipos[0], asipos[1], color = blk, $
;            'ASI: '+time_string(tut,tformat='hh:mm:ss'), /normal
;        sgclose
;    endfor


end
