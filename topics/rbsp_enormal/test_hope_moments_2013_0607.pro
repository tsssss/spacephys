;+
; test rbsp hope moments for 2013_0607 event.
;-



;---Settings.
    rootdir = '/Users/Sheng/Google Drive'
    utr0 = time_double(['2013-06-07/04:45','2013-06-07/05:15'])
    utr0 = time_double(['2013-06-07/04:52','2013-06-07/05:05'])
    utra = time_double(['2013-06-07/04:55','2013-06-07/05:00'])
    utrb = time_double(['2013-06-07/04:55','2013-06-07/05:00'])
    

    afn = rootdir+'/rbspa_hope_moments_2013_0607_v01.cdf'
    bfn = rootdir+'/rbspb_hope_moments_2013_0607_v01.cdf'
    
    ions = ['p','o','he']
    ionmass = [1d,16,4]
    allsp = ['e',ions]
    
    mapidx = 2  ; 't01'.
    mu0 = 4*!dpi*1e-7
    va0 = 22d

    probes = ['a','b']
    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        tplot_options, 'constant', 0
        tplot_options, 'labflag', -1
        tplot_options, 'ynozero', 1
        
        
        ;---Load hope moments, and B field.
        _2013_0607_load_data
        cdf2tplot, afn, prefix='rbspa_'
        cdf2tplot, bfn, prefix='rbspb_'
        
        
        ;---Test density.
        ; shows: all ions add up to electron density, o+ is 4 times of h+ sometimes.
        vars = pre0+['p','o','he']+'_density'
        stplot_merge, vars, newname=pre0+'i_density', limits={ytitle:'(cm!U-3!N)',colors:[6,4,2],labels:['H+','O+','He+']}
        
        get_data, pre0+ions[0]+'_density', uts
        nrec = n_elements(uts)
        mhdnden = dblarr(nrec)
        for i=0, n_elements(ions)-1 do begin
            get_data, pre0+ions[i]+'_density', uts, nden
            mhdnden += nden
        endfor
        store_data, pre0+'mhd_density', uts, mhdnden
        get_data, pre0+'e_density', tuts, eden
        store_data, pre0+'hope_density', uts, [[mhdnden],[interpol(eden,tuts,uts)]], $
            limits={ytitle:'(cm!U-3!N)', labels:['ion','ele'], colors:[6,2]}
            
        ;    vars = pre0+['e','mhd','i']+'_density'
        ;    tplot, vars, trange=utr0
        ;    stop
        
        ;---Test velocity.
        ; shows: oxygen velocity is almost the bulk velocity.
        vars = pre0+ions+'_vbulk'
        options, vars, 'colors', [6,4,2]
        options, vars, 'labels', ['x','y','z']
        
        get_data, pre0+ions[0]+'_vbulk', uts
        nrec = n_elements(uts)
        mhdmden = dblarr(nrec)
        mhdvbulk = dblarr(nrec,3)
        for i=0, n_elements(ions)-1 do begin
            get_data, pre0+ions[i]+'_density', uts, nden
            get_data, pre0+ions[i]+'_vbulk', uts, vbulk
            mhdmden += nden*ionmass[i]
            for j=0, 2 do mhdvbulk[*,j] += vbulk[*,j]*nden*ionmass[i]
        endfor
        for j=0, 2 do mhdvbulk[*,j] *= (1d/mhdmden)
        store_data, pre0+'mhd_vbulk', uts, mhdvbulk, limits={colors:[6,4,2],labels:['x','y','z']}
        store_data, pre0+'mhd_mden', uts, mhdmden
        
        vars = pre0+'mhd_vbulk'
        get_data, vars[0], uts
        get_data, pre0+'b0_gse', tuts, b0gse
        b0gse = sinterpol(b0gse, tuts, uts)
        b0gsm = sgse2gsm(b0gse, stoepoch(uts,'unix'))
        b0hat = sunitvec(b0gsm)
        
        
        
        ;---Test energy flux, enthalpy.
        ; h+ energy flux is usually much larger than o+.
        ; shows: enthalpy is the largest contribution to energy flux.
        vars = pre0+allsp+'_energy_flux'
        options, vars, 'colors', [6,4,2]
        options, vars, 'labels', ['x','y','z']
        
        vars = pre0+allsp+'_enthalpy'
        options, vars, 'colors', [6,4,2]
        options, vars, 'labels', ['x','y','z']
        
        
        get_data, pre0+ions[0]+'_energy_flux'
        nrec = n_elements(uts)
        ioneflux = dblarr(nrec,3)
        for i=0, n_elements(ions)-1 do begin
            get_data, pre0+ions[i]+'_energy_flux', uts, eflux
            ioneflux += eflux
        endfor
        store_data, pre0+'i_energy_flux', uts, ioneflux, limits={colors:[6,4,2],labels:['x','y','z']}
        
        
        get_data, pre0+ions[0]+'_enthalpy'
        nrec = n_elements(uts)
        ionenthalpy = dblarr(nrec,3)
        for i=0, n_elements(ions)-1 do begin
            get_data, pre0+ions[i]+'_enthalpy', uts, eflux
            ionenthalpy += eflux
        endfor
        store_data, pre0+'i_enthalpy', uts, ionenthalpy, limits={colors:[6,4,2],labels:['x','y','z']}
        
        
        
        ; rotate vectors to fac.
        vars = pre0+['i_enthalpy','i_energy_flux','mhd_vbulk']
        get_data, vars[0], uts
        get_data, pre0+'b0_gse', tuts, b0gse
        b0gse = sinterpol(b0gse, tuts, uts)
        b0gsm = sgse2gsm(b0gse, stoepoch(uts,'unix'))
        bhat = sunitvec(b0gsm)
        get_data, pre0+'pos_gse', tuts, rgse
        rgse = sinterpol(rgse, tuts, uts)
        rgsm = sgse2gsm(rgse, stoepoch(uts,'unix'))
        rhat = sunitvec(rgsm)
        phat = sunitvec(scross(rhat,bhat))
        vhat = scross(bhat,phat)
        
        foreach tvar, vars do begin
            get_data, tvar, uts, dat
            store_data, tvar+'_fac', uts, [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]], $
                limits={colors:[6,4,2], labels:['para','west','north']}
        endforeach
        
        vars = pre0+['e_enthalpy','e_energy_flux']
        get_data, vars[0], uts
        get_data, pre0+'b0_gse', tuts, b0gse
        b0gse = sinterpol(b0gse, tuts, uts)
        b0gsm = sgse2gsm(b0gse, stoepoch(uts,'unix'))
        bhat = sunitvec(b0gsm)
        get_data, pre0+'pos_gse', tuts, rgse
        rgse = sinterpol(rgse, tuts, uts)
        rgsm = sgse2gsm(rgse, stoepoch(uts,'unix'))
        rhat = sunitvec(rgsm)
        phat = sunitvec(scross(rhat,bhat))
        vhat = scross(bhat,phat)
        
        foreach tvar, vars do begin
            get_data, tvar, uts, dat
            store_data, tvar+'_fac', uts, [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]], $
                limits={colors:[6,4,2], labels:['para','west','north']}
        endforeach
        
        
        ; map.
        vars = pre0+['e_enthalpy','e_energy_flux','i_energy_flux','pf']+'_fac'
        foreach tvar, vars do begin
            get_data, tvar, uts, dat, limits=lim
            get_data, pre0+'map_coef', tuts, mcoef
            mcoef = sinterpol(mcoef[*,mapidx], tuts, uts)
            for i=0, 2 do dat[*,i] *= mcoef
            store_data, tvar+'_map', uts, dat, limits=lim
        endforeach
        
        
        
        ; |B|, V_alfven.
        get_data, pre0+'b_gse', uts, dat
        dat = snorm(dat) & dat = smooth(dat, 11d/sdatarate(uts), /edge_truncate)
        store_data, pre0+'bmag', uts, dat, limits={ytitle:'(nT)', labels:'|B|'}
        
        get_data, pre0+'mhd_mden', uts, mden
        get_data, pre0+'bmag', tuts, bmag
        bmag = interpol(bmag, tuts, uts)
        va = bmag/sqrt(mden)*va0
        store_data, pre0+'va', uts, va
        
        
        ; average westward velocity.
        get_data, pre0+'mhd_vbulk_fac', uts, dat
        if tprobe eq 'a' then tutr = utra else tutr = utrb
        idx = where(uts ge tutr[0] and uts le tutr[1])
        dat = dat[idx,1]
        print, 'mean, stddev Vbulk (km/s):', mean(dat), stddev(dat)
        ; RBSP-A: 72+/-50 km/s, or 22-122 km/s.
        ; RBSP-B: 118+/-77 km/s, or 41-195 km/s.
        
        ofn = shomedir()+'/fig_2013_0607_'+pre0+'moments.pdf'
        ofn = 0
        sgopen, ofn, xsize=8, ysize=10, /inch
        
        device, decomposed=0
        loadct2, 43
        
        xchsz = double(!d.x_ch_size)/!d.x_size
        ychsz = double(!d.y_ch_size)/!d.y_size
        
        
        tvar = pre0+['e','i']+'_energy_flux_fac_map'
        options, tvar, 'ytitle', '(mW/m!U2!N)'
        
        tvar = pre0+'pf_fac_map'
        options, tvar, 'labels', ['para','west','north']
        
        tvar = pre0+'mhd_vbulk_fac'
        options, tvar, 'ytitle', '(km/s)'
        
        
        vars = pre0+[['hope','i']+'_density',['e','i']+'_energy_flux_fac_map','pf_fac_map','bmag','mhd_vbulk_fac']
        nvar = n_elements(vars)
        
        
        figlab = ['a','b','c','d','e','f','g']+'. '+['HOPE density','Ion density',['Ele Eflux FAC','Ion Eflux FAC','S FAC']+'!C  @100km', '|B|', 'MHD Vbulk']
        
        poss = sgcalcpos(nvar, lmargin=20, rmargin=8)
        tplot, vars, trange=utr0, position=poss, title='RBSP-'+strupcase(tprobe)
        
        for i=0, nvar-1 do $
            xyouts, poss[0,i]-xchsz*15, poss[3,i]-ychsz*0.5, /normal, figlab[i]
        
        
        vars = pre0+[['hope','i']+'_density',['e','i']+'_energy_flux_fac_map','pf_fac_map','bmag','va','mhd_vbulk_fac']
        tplot, vars, trange=utr0
        stop
        sgclose
    endforeach    

end
