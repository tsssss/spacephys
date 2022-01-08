
utr = time_string(['2013-06-07/04:52','2013-06-07/05:02'])
utr = time_string(['2013-06-07/04:45','2013-06-07/05:15'])

probes = ['a','b']
tprobe = 'a'
pre0 = 'rbsp'+tprobe+'_'
nsec = 16
secs = smkarthm(0,360d/nsec,nsec,'x0')

;---load data.
    datfn = srootdir()+'/test_2013_0607_hope_flow_velocity.tplot'
    if file_test(datfn) eq 0 then load = 1 else load = 0
    if load then begin
        
        ;rbsp_load_spice_kernels, /all
        
        foreach tprobe, probes do begin
            pre0 = 'rbsp'+tprobe+'_'
            
            hopel2 = sread_rbsp_hope_l2(utr, probe=tprobe)
            store_data, pre0+'hopel2', utr, hopel2
            
            ets = hopel2.epoch
            uts = sfmepoch(ets, 'unix')
            nrec = n_elements(ets)
            uvws = dblarr(nrec,9)
            for i = 0, nrec-1 do begin
                tstr = time_string(uts[i], tformat = 'YYYY-MM-DDThh:mm:ss.ffffff')
                cspice_str2et, tstr, et
                cspice_pxform, 'RBSP'+strupcase(tprobe)+'_SCIENCE', 'GSM', et, pxform
                uvws[i,*] = pxform[*]   ; the colums are (u)(v)(w) in gsm.
            endfor
            
            tvar = pre0+'uvw_gsm'
            store_data, tvar, uts, uvws
        endforeach
        
        vars = ['hopel2','uvw_gsm']
        vars = ['rbspa_'+vars, 'rbspb_'+vars]
        tplot_save, vars, filename=datfn
    endif else begin
        store_data, '*', /delete
        tplot_restore, filename=datfn
    endelse


    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        get_data, pre0+'hopel2', utr, hopel2
        uts0 = sfmepoch(hopel2.epoch_ion,'unix')
        nrec = n_elements(uts0)
        dr0 = sdatarate(uts0)
        hflux0 = total(hopel2.fpdu[*,*,*,2],2,/nan)
        
        uts1 = sfmepoch(hopel2.epoch_ele,'unix')
        eflux0 = total(hopel2.fedu[*,*,*,2],2,/nan)
        
        
        tvar = pre0+'l2_p_flux'
        store_data, tvar, uts0, hflux0, secs, $
            limits={ylog:0, yrange:[0,360], ystyle:1, zlog:1, spec:1, no_interp:1}
            
        tvar = pre0+'l2_e_flux'
        store_data, tvar, uts1, eflux0, secs, $
            limits={ylog:0, yrange:[0,360], ystyle:1, zlog:1, spec:1, no_interp:1, zrange:[1e9,1e11]}
            
            
        get_data, pre0+'uvw_gsm', uts, uvws
        store_data, pre0+'u_gsm', uts, uvws[*,[0,3,6]], limits = {colors:[6,4,2], labels:['x','y','z']}
        store_data, pre0+'v_gsm', uts, uvws[*,[1,4,7]], limits = {colors:[6,4,2], labels:['x','y','z']}
        store_data, pre0+'w_gsm', uts, uvws[*,[2,5,8]], limits = {colors:[6,4,2], labels:['x','y','z']}
    endforeach
    

    ofn = 0
    ofn = shomedir()+'/hope_l2_px2_2013_0607.pdf'
    sgopen, ofn, xsize=8.5, ysize=11, /inch
    
    device, decomposed=0
    loadct2, 43
    
    vars = ['l2_e_flux','u_gsm']
    vars = ['rbspa_'+vars, 'rbspb_'+vars]
    tplot, vars, trange=utr

    sgclose


end
