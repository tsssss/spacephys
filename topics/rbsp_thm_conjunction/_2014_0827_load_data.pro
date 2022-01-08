pro _2014_0827_load_data

;---Constant.
    deg = 180d/!dpi
    rad = !dpi/180d
    re = 6378d & re1 = 1d/re
    r0 = 100d/re+1
    secofday = 86400d


;---Results from looking data.
    ; when DF arrive the s/c.
    pre0s = ['thd','the','tha','rba']+'_'
    dfuts = time_double('2014-08-27/'+['09:28','09:35','09:42','09:48'])
    fmlts = [0.11d,1.13,1.90,4.64]
    fdmlts = [0.01d,0.09,0.12,0.20]
    fmlats = [63.4d,66.6,68.4,64.9]
    fdmlats = [2.3d,3.4,3.6,1.2]
    fvels = [0,2.2,1.6,6.8]
    
;    for i=0,2 do begin
;        dt = dfuts[i+1]-dfuts[i]
;        dmlt = fmlts[i+1]-fmlts[i]
;        vel = dmlt*15/dt*60
;        print, vel
;    endfor
;    stop



;---Settings
    utr0 = time_double(['2014-08-27/08:30','2014-08:27/10:30']) ; time to load data.
    utr1 = time_double(['2014-08-27/08:40','2014-08:27/10:10']) ; time for plot.
    utr1 = time_double(['2014-08-27/09:10','2014-08:27/10:10']) ; time for plot.
    utr2 = utr0+[-1,1]*3600d    ; a more extended time.
    
    
    reload = 0
    reload_pos = 0  ; in-situ positions.
    reload_map = 0  ; mapped positions.
    reload_esa = 0  ; particle data.
    reload_ebf = 0  ; E and B fields.
    reload_ewo = 0  ; ewogram for atha and pina at different mlat.
    
    
    rootdir = shomedir()
    datfn = rootdir+'/2014_0827_all_data.tplot'

    
    allvars = []
    
    
    rb_probes = ['a']
    th_probes = ['a','d','e']
    allpre0 = ['rb'+rb_probes,'th'+th_probes]+'_'
    
    
    ; setting for mapping.
    models = ['t89','t96','t01','t04s']
    nmodel = n_elements(models)
    dir = -1    ; map to northern hemisphere.

    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    
    tplot_options, 'labflag', -1
    tplot_options, 'constant', 0
    tplot_options, 'xticklen', -0.02
    tplot_options, 'yticklen', -0.01
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 6

    
    if file_test(datfn) eq 1 then begin
        tplot_restore, filename=datfn
        if ~keyword_set(reload) then return
    endif
    
    foreach pre0, pre0s, i do begin
        get_data, pre0+'fpt_mlon', uts, fmlts
        tmlts = sinterpol(fmlts, uts, dfuts[i])
        print, pre0
        print, minmax(tmlts)
        print, mean(tmlts), stddev(tmlts)
    endforeach
    


;---Load esa.
    load = 0
    tvar = ['density','vbulk','ele_enspec','ion_enspec']
    esavars = []
    foreach pre0, ['rb'+rb_probes,'th'+th_probes]+'_' do esavars = [esavars,pre0+tvar]
    allvars = [allvars,esavars]
    
    ion_en_yrng = [10,4e4]
    ele_en_yrng = [10,4e4]

    foreach tvar, esavars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_esa) then load = 1
    
    if load then begin
        probes = rb_probes
        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            hope = sread_rbsp_hope_moments(utr0, probe=tprobe)
            uts = hope.ut_ele
            store_data, pre0+'density', uts, hope.e_density, limits={ytitle:'(cm!U-3!N)', labels:'Density'}
            uts = hope.ut_ion
            nrec = n_elements(uts)
            vbulk = dblarr(nrec,3)
            mdens = dblarr(nrec)
            ions = ['p','o','he']
            mass = [1d,16,4]
            tnames = strlowcase(tag_names(hope))
            foreach tion, ions, i do begin
                idx = where(tnames eq tion+'_density')
                ndens = hope.(idx)
                mdens += mass[i]*ndens
                idx = where(tnames eq tion+'_vbulk')
                tvbulk = hope.(idx)
                for j=0,2 do vbulk[*,j] += mass[i]*ndens*tvbulk[*,j]
            endforeach
            for i=0,2 do vbulk[*,i] *= (1d/mdens)
            store_data, pre0+'vbulk', uts, vbulk, limits={ytitle:'(km/s)', colors:rgb, labels:'GSM V'+xyz}
            
            vars = ['PITCH_ANGLE',$
                'Epoch_Ele','HOPE_ENERGY_Ele',$
                'Epoch_Ion','HOPE_ENERGY_Ion',$
                'FEDU','FPDU','FODU']
            hopel3 = sread_rbsp_hope_l3(utr0, probe=tprobe, vars=vars)
            ens = hopel3.hope_energy_ion
            dat = total(hopel3.fpdu,/nan,3) ; in (1/cm^2-s-sr-keV).
            dat = dat*ens*1e-3   ; in (eV/cm^2-s-sr-eV).
            idx = where(finite(dat,/nan),cnt)
            if cnt ne 0 then dat[idx] = 0
            store_data, pre0+'ion_enspec', uts, dat, ens, limits=$
                {ytitle:'Energy (eV)', ztitle:'(eV/cm!U2!N-s-sr-eV)', spec:1, no_interp:1, $
                zlog:1, ylog:1, ystyle:1, yrange:ion_en_yrng}

            
            uts = sfmepoch(hopel3.epoch_ele,'unix')
            ens = hopel3.hope_energy_ele
            dat = total(hopel3.fedu,/nan,3)
            dat = dat*ens*1e-3   ; in (eV/cm^2-s-sr-eV).
            idx = where(finite(dat,/nan),cnt)
            if cnt ne 0 then dat[idx] = 0
            store_data, pre0+'ele_enspec', uts, dat, ens, limits=$
                {ytitle:'Energy (eV)', ztitle:'(eV/cm!U2!N-s-sr-eV)', spec:1, no_interp:1, $
                zlog:1, ylog:1, ystyle:1, yrange:ele_en_yrng}
        endforeach
        
        probes = th_probes
        foreach tprobe, probes do begin
            pre0 = 'th'+tprobe+'_'
            var1s = ['peir_'+['time','velocity_gsm','en_eflux','en_eflux_yaxis'],$
                'peer_'+['time','density','en_eflux','en_eflux_yaxis']]
            var0s = pre0+var1s
            esa = sread_thm_esa_l2(utr0, probe=tprobe, vars=var0s, newname=var1s)
            uts = esa.peir_time
            store_data, pre0+'vbulk', uts, esa.peir_velocity_gsm, limits={ytitle:'(km/s)',colors:rgb,labels:'GSM V'+xyz}
            store_data, pre0+'ion_enspec', uts, esa.peir_en_eflux, esa.peir_en_eflux_yaxis, $
                limits={ytitle:'Energy (eV)', ztitle:'(eV/cm!U2!N-s-sr-eV)', spec:1, no_interp:1, $
                zlog:1, ylog:1, ystyle:1, yrange:ion_en_yrng}
            uts = esa.peer_time
            store_data, pre0+'density', uts, esa.peer_density, limits={ytitle:'(cm!U-3!N)', labels:'Density'}
            store_data, pre0+'ele_enspec', uts, esa.peer_en_eflux, esa.peer_en_eflux_yaxis, $
                limits={ytitle:'Energy (eV)', ztitle:'(eV/cm!U2!N-s-sr-eV)', spec:1, no_interp:1, $
                zlog:1, ylog:1, ystyle:1, yrange:ele_en_yrng}
        endforeach
        
        options, 'rb'+rb_probes+'_ele_enspec', 'zrange', [1e5,1e9]
        options, 'th'+th_probes+'_ele_enspec', 'zrange', [1e4,1e8]
        
        options, 'rb'+rb_probes+'_ion_enspec', 'zrange', [1e5,1e8]
        options, 'th'+th_probes+'_ion_enspec', 'zrange', [1e4,1e7]

        tplot_save, allvars, filename=datfn
    endif
    
    
    

;---Load position.
    load = 0
    tvar = ['r_gsm','mlt','mlat','lshell']
    posvars = []
    foreach pre0, ['rb'+rb_probes,'th'+th_probes]+'_' do posvars = [posvars,pre0+tvar]
    allvars = [allvars,posvars]
    
    foreach tvar, posvars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_pos) then load = 1
    
    if load then begin
        probes = rb_probes
        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            pos = sread_rbsp_spice_product(utr0, probe=tprobe)
            uts = pos.ut_pos
            store_data, pre0+'r_gsm', uts, pos.pos_gsm, limits={ytitle:'(Re)',colors:rgb, labels:xyz}
            store_data, pre0+'mlt', uts, pos.mlt
            store_data, pre0+'mlat', uts, pos.mlat*deg
            store_data, pre0+'lshell', uts, pos.lshell
        endforeach

        probes = th_probes
        foreach tprobe, probes do begin
            pre0 = 'th'+tprobe+'_'
            var0s = ['Epoch','XYZ_GSM','SM_LCT_T','SM_LAT','L_VALUE']
            var1s = ['epoch','pos_gsm','mlt','mlat','lshell']
            pos = sread_thm_orbit(utr0, probe=tprobe, vars=var0s, newname=var1s)
            uts = sfmepoch(pos.epoch,'unix')
            idx = where(pos.mlt ge 12, cnt)
            if cnt ne 0 then pos.mlt[idx] -= 24
            idx = where(uts ge utr2[0] and uts le utr2[1])
            uts = uts[idx]
            store_data, pre0+'r_gsm', uts, pos.pos_gsm[idx,*], limits={ytitle:'(Re)',colors:rgb,labels:xyz}
            store_data, pre0+'mlt', uts, pos.mlt[idx]
            store_data, pre0+'mlat', uts, pos.mlat[idx]
            store_data, pre0+'lshell', uts, pos.lshell[idx]
        endforeach
        
        tplot_save, allvars, filename=datfn
    endif


;---Load mapped position.
    load = 0
    tvar = 'fpt_'+['mlt','mlat','mlon']
    mapvars = models+'_par'
    foreach pre0, ['rb'+rb_probes,'th'+th_probes]+'_' do mapvars = [mapvars,pre0+tvar]
    allvars = [allvars,mapvars]

    foreach tvar, mapvars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_map) then load = 1
    
    if load then begin
        ; load geopack par's. [models_par].
        foreach tmodel, models do sgeopack_par, utr0, tmodel
        
        foreach pre0, ['rb'+rb_probes,'th'+th_probes]+'_' do begin
            get_data, pre0+'r_gsm', tuts, posgsm
            tets = stoepoch(tuts, 'unix')
            tnrec = n_elements(tuts)
            hfmlats = dblarr(tnrec,nmodel)
            hfmlons = dblarr(tnrec,nmodel)
            hfmlts = dblarr(tnrec,nmodel)
            
            for j=0, nmodel-1 do begin
                tmodel = models[j]
                t89 = (tmodel eq 't89')? 1: 0
                t96 = (tmodel eq 't96')? 1: 0
                t01 = (tmodel eq 't01')? 1: 0
                t04s = (tmodel eq 't04s')? 1: 0
                storm = (tmodel eq 't04s')? 1: 0
                
                get_data, tmodel+'_par', tmp, pars
                pars = sinterpol(pars, tmp, tuts)
                tnrec = n_elements(tuts)
                
                for i=0, tnrec-1 do begin
                    tet = tets[i]
                    geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
                    geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
                    
                    xp = posgsm[i,0]
                    yp = posgsm[i,1]
                    zp = posgsm[i,2]
                    ; map to the ionosphere.
                    ; [xyz]f is the footpoint position.
                    par = reform(pars[i,*])
                    geopack_trace, xp,yp,zp, dir, par, xf,yf,zf, r0=r0, $
                        /refine, /ionosphere, $
                        t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm
                    ; convert footpoint from gsm to mag.
                    geopack_conv_coord, xf,yf,zf, /from_gsm, $
                        txf,tyf,tzf, /to_mag
                    hfmlats[i,j] = asin(tzf/r0)*deg
                    hfmlons[i,j] = atan(tyf,txf)*deg
                    hfmlts[i,j] = slon2lt(hfmlons[i,j], tet, /mag, /deg)/15
                endfor
            endfor
            store_data, pre0+'fpt_mlat', tuts, hfmlats
            store_data, pre0+'fpt_mlon', tuts, hfmlons
            store_data, pre0+'fpt_mlt', tuts, hfmlts
        endforeach
        
        tplot_save, allvars, filename=datfn
    endif


;---Load EWOgram.
    load = 0
    sites = ['atha','pina']
    dmlat = 1d  ; deg.
    dmlt = 0.05 ; hr.
    mlat0s = smkarthm(62,65,dmlat,'dx')
    ewovars = []
    foreach tsite,sites do ewovars = [ewovars,tsite+'_ewo_'+string(mlat0s,format='(I0)')]
    allvars = [allvars,ewovars]
    
    minelev = 5d
    height = 110d       ; km.
    hgtidx = where([90,110,150] eq height)

    foreach tvar, ewovars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_ewo) then load = 1

    if load then begin
        foreach tmlat, mlat0s do begin
            mlatrng = tmlat+[0,dmlat]
            tdatfn = rootdir+'/ewo_'+time_string(utr0[0],tformat='YYYY_MMDD')+'_'+sgnum2str(mlatrng[0])+'deg.tplot'
            if file_test(tdatfn) ne 0 then tplot_restore, filename=tdatfn else begin
                foreach tsite, sites do begin
                    asf = sread_thg_asi(utr0, tsite, type='asf')
                    asc = sread_thg_asc(0, tsite, type='asf', vars=['mlat','mlon','elev'])
                    
                    uts = asf.utsec
                    nrec = n_elements(uts)
                    mlats = reform(asc.(0).mlat[hgtidx,*,*])
                    mlons = reform(asc.(0).mlon[hgtidx,*,*])
                    elevs = asc.(0).elev
                    midns = -15.5179*(uts*(1d/secofday) mod 1)*24+72.5124   ; midnight mlon.
                    
                    ; converts corner position to center position.
                    mlats = (mlats[1:*,1:*]+mlats[1:*,0:255]+mlats[0:255,1:*]+mlats[0:255,0:255])*0.25
                    mlons = (mlons[1:*,1:*]+mlons[1:*,0:255]+mlons[0:255,1:*]+mlons[0:255,0:255])*0.25
                    
                    tmp = minmax(mlons) & tmp = [tmp-min(midns),tmp-max(midns)]/15
                    mltrng = minmax(tmp) & mltrng = mltrng-(mltrng mod dmlt)+[0,dmlt]
                    mlts = smkarthm(mltrng[0], mltrng[1], dmlt, 'dx')
                    nmlt = n_elements(mlts)
                    
                    ewos = dblarr(nrec,nmlt)
                    for i=0, nrec-1 do begin
                        img = double(reform(asf.img[i,*,*]))
                        ; remove edge.
                        edge = where(elevs lt minelev or ~finite(elevs))
                        img[edge] = mean(img[0:10,0:10], /nan)
                        ; crude corner average subtraction. use thg_asf_site_offset?
                        img = img - img[0] > 0
                        ; scale luminosity, adopted from thm_asi_merge_mosaic.
                        img *= 64d/(median(img) > 1)
                        tmlts = (mlons-midns[i])/15
                        ; map to ewogram.
                        for j=0, nmlt-1 do begin
                            idx = where(mlats ge mlatrng[0] and mlats le mlatrng[1] $
                                and tmlts ge mlts[j]-dmlt*0.5 and tmlts lt mlts[j]+dmlt*0.5, cnt, complement=idx2)
                            if cnt ne 0 then ewos[i,j] = mean(img[idx],/nan)
                        endfor
                    endfor
                    
                    tvar = tsite+'_ewo_'+string(tmlat,format='(I0)')
                    idx = where(ewos eq 0, cnt)
                    if cnt ne 0 then ewos[idx] = !values.d_nan
                    store_data, tvar, uts, ewos, mlts, limits={ytitle:'(hr)', spec:1, yrange:mltrng, no_interp:1}
                endforeach
                tplot_save, sites+'_ewo_'+string(tmlat,format='(I0)'), filename=tdatfn
            endelse
        endforeach

        tplot_save, allvars, filename=datfn
    endif
    

;---Load E/B field.
    load = 0
    tvar = ['e_gsm','b_gsm','b0_gsm','edot0_gsm']
    ebfvars = models+'_par'
    foreach pre0, ['rb'+rb_probes,'th'+th_probes]+'_' do ebfvars = [ebfvars,pre0+tvar]
    ebfvars = [ebfvars,'rb'+rb_probes+'_quvw2gsm']
    allvars = [allvars,ebfvars]
        
    foreach tvar, ebfvars do if tnames(tvar) eq '' then load = 1
    if keyword_set(reload_ebf) then load = 1

    bmin = 5d   ; nT.
    rmin = 0.25 ; percent.
    
    if load then begin
        probes = rb_probes
        spinrate = 12d  ; sec.
        dr0 = 1d/16
        blen0s = [100d,100,12]  ; m, twice of boom lengths.
        fshrts = [1d,1,1]       ; use 1 temperarily.
        blen1s = blen0s*fshrts
        
        ; uniform time for E/B.
        uts = smkarthm(utr0[0], utr0[1], dr0, 'dx')

        foreach tprobe, probes do begin
            pre0 = 'rb'+tprobe+'_'
            
            ; bfield.
            emfisis = sread_rbsp_emfisis_l3(utr0, type='hires', probe=tprobe, coord='gsm')
            tuts = sfmepoch(emfisis.epoch,'unix')
            bgsm = sinterpol(emfisis.mag,tuts,uts)
            store_data, pre0+'b_gsm', uts, bgsm, $
                limit = {ytitle:'(nT)', labels:'GSM B'+xyz, colors:rgb}
            b0gsm = bgsm
            for i=0,2 do b0gsm[*,i] = scalcbg(bgsm[*,i])
            store_data, pre0+'b0_gsm', uts, b0gsm, $
                limit = {ytitle:'(nT)', labels:'GSM B0'+xyz, colors:rgb}
            

            ; efield.
            spice = sread_rbsp_spice_product(utr1, probe=tprobe)
            tuts = spice.ut_cotran
            quvw2gsm = qslerp(spice.q_uvw2gsm, tuts, uts)
            store_data, pre0+'quvw2gsm', uts, quvw2gsm
            muvw2gsm = transpose(qtom(quvw2gsm))

            
            dat = sread_rbsp_efw_l2(utr1, probes=tprobe, type='vsvy')
            tuts = sfmepoch(dat.epoch, 'unix')
            nrec = n_elements(uts)
            vsvy = dat.vsvy
            vsvy = sinterpol(vsvy, tuts, uts)

            eu = (vsvy[*,0]-vsvy[*,1])/blen1s[0]*1e3   ; V -> V/m -> mV/m.
            ev = (vsvy[*,2]-vsvy[*,3])/blen1s[1]*1e3
            ew = dblarr(nrec)
            
            nspin = 1
            tnrec = nspin*spinrate/dr0
            eu = eu-smooth(eu, tnrec, /edge_truncate, /nan)
            ev = ev-smooth(ev, tnrec, /edge_truncate, /nan)
            euvw = [[eu],[ev],[ew]]

            ; rotate to gsm.
            ex = eu*muvw2gsm[0,0,*] + ev*muvw2gsm[1,0,*] + ew*muvw2gsm[2,0,*]
            ey = eu*muvw2gsm[0,1,*] + ev*muvw2gsm[1,1,*] + ew*muvw2gsm[2,1,*]
            ez = eu*muvw2gsm[0,2,*] + ev*muvw2gsm[1,2,*] + ew*muvw2gsm[2,2,*]
            
            tvar = pre0+'e_gsm'
            store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
                {ytitle:'(mV/m)', colors:rgb, labels:'GSM E'+xyz}


            ; calc e dot0.
            get_data, pre0+'b0_gsm', uts, b0gsm
            bx = b0gsm[*,0] & by = b0gsm[*,1] & bz = b0gsm[*,2]
            bu = bx*muvw2gsm[0,0,*] + by*muvw2gsm[0,1,*] + bz*muvw2gsm[0,2,*]
            bv = bx*muvw2gsm[1,0,*] + by*muvw2gsm[1,1,*] + bz*muvw2gsm[1,2,*]
            bw = bx*muvw2gsm[2,0,*] + by*muvw2gsm[2,1,*] + bz*muvw2gsm[2,2,*]
            
            buvw = [[bu],[bv],[bw]]
            bmag = snorm(buvw)
            euvw[*,2] = -(euvw[*,0]*buvw[*,0]+euvw[*,1]*buvw[*,1])/buvw[*,2]
            idx = where(bmag le bmin or abs(buvw[*,2])/bmag le rmin, cnt)
            if cnt ne 0 then euvw[idx,2] = !values.d_nan
            
            ex = euvw[*,0]*muvw2gsm[0,0,*] + euvw[*,1]*muvw2gsm[1,0,*] + euvw[*,2]*muvw2gsm[2,0,*]
            ey = euvw[*,0]*muvw2gsm[0,1,*] + euvw[*,1]*muvw2gsm[1,1,*] + euvw[*,2]*muvw2gsm[2,1,*]
            ez = euvw[*,0]*muvw2gsm[0,2,*] + euvw[*,1]*muvw2gsm[1,2,*] + euvw[*,2]*muvw2gsm[2,2,*]

            tvar = pre0+'edot0_gsm'
            store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
                {ytitle:'(mV/m)', colors:rgb, labels:'GSM Edot0'+xyz}
        endforeach
        
        probes = th_probes
        spinrate = 3d   ; sec.
        dr0 = 1d/16     ; what's themis data rate??
        
        foreach tprobe, probes do begin
            pre0 = 'th'+tprobe+'_'
            
            ; thx_e_dsl
            efsvar0s = pre0+['efs_dot0_dsl','efs_dot0_time']
            efsvar1s = ['efs_dsl','efs_utsec']
            efs = sread_thm_efi_l2(utr0, vars=efsvar0s, newname=efsvar1s, probes=tprobe)
            
            tvar = pre0+'e_dsl1'
            uts = efs.efs_utsec
            edsl = efs.efs_dsl
            edsl[*,2] = 0
            if tprobe eq 'd' then begin
                ; remove bad data for TH-D.
                tut = time_double('2014-08-27/09:25')
                idx = where(uts le tut, cnt)
                edsl[idx,*] = !values.d_nan
            endif
            ; remove offset.
            for i=0, 1 do edsl[*,i] = edsl[*,i]-scalcbg(edsl[*,i])
            ; despike.
            sdespike, edsl
            store_data, tvar, uts, edsl, limits={ytitle:'(mV/m)', colors:rgb, labels:'DSL E0'+xyz}
            emag = sqrt(edsl[*,0]^2+edsl[*,1]^2)
            store_data, pre0+'emag', uts, emag, limits={ytitle:'(mV/m)', labels:['|E|']}

            ; thx_b_dsl.
            fgmvar0s = pre0+['fgs_dsl','fgs_time']
            fgmvar1s = ['fgs_dsl','fgs_utsec']
            fgm = sread_thm_fgm_l2(utr0, vars=fgmvar0s, newname=fgmvar1s, probes=tprobe)
            tvar = pre0+'b_dsl'
            tuts = fgm.fgs_utsec
            bdsl = fgm.fgs_dsl
            ; remove bad data for TH-D.
            if tprobe eq 'd' then begin
                tut = time_double('2014-08-27/09:25')
                idx = where(uts le tut, cnt)
                if cnt ne 0 then bdsl[idx,*] = !values.d_nan
            endif
            store_data, tvar, tuts, bdsl, limits={ytitle:'(nT)', colors:rgb, labels:'DSL B'+xyz, ystyle:0}
            
            
            ; calc E dot0
            ; thx_[b_dsl1,b0,db_dsl,e_dsl].
            get_data, pre0+'b_dsl', tuts, bdsl
            get_data, pre0+'e_dsl1', uts, edsl
            bdsl = sinterpol(bdsl, tuts, uts)
            for i=0,2 do bdsl[*,i] = scalcbg(bdsl[*,i])
            bmag = snorm(bdsl)
            edsl[*,2] = -(edsl[*,0]*bdsl[*,0]+edsl[*,1]*bdsl[*,1])/bdsl[*,2]
            idx = where(bmag le bmin or abs(bdsl[*,2])/bmag le rmin, cnt)
            if cnt ne 0 then edsl[idx,2] = !values.d_nan

            store_data, pre0+'e_dsl', uts, edsl, limits={ytitle:'(mV/m)', colors:rgb, labels:'DSL E0'+xyz}
            store_data, pre0+'b0_dsl', uts, bdsl, limits={ytitle:'(nT)', colors:rgb, labels:'DSL B'+xyz}


            ; rotate to GSM.
            vars = pre0+['e_dsl','b_dsl']
            thm_load_state, /get_support, probe=tprobe, trange=utr0
            thm_cotrans, probe=tprobe, in_coord='dsl', out_coord='gsm', pre0+'e_dsl', pre0+'edot0_gsm'
            thm_cotrans, probe=tprobe, in_coord='dsl', out_coord='gsm', pre0+'e_dsl1', pre0+'e_gsm'
            thm_cotrans, probe=tprobe, in_coord='dsl', out_coord='gsm', pre0+'b_dsl', pre0+'b_gsm'
            thm_cotrans, probe=tprobe, in_coord='dsl', out_coord='gsm', pre0+'b0_dsl', pre0+'b0_gsm'
            store_data, pre0+'state*', /delete
            
            options, pre0+'e_gsm', 'colors', rgb
            options, pre0+'e_gsm', 'labels', 'GSM E'+xyz
            options, pre0+'edot0_gsm', 'colors', rgb
            options, pre0+'edot0_gsm', 'labels', 'GSM Edot0'+xyz
            options, pre0+'b_gsm', 'colors', rgb
            options, pre0+'b_gsm', 'labels', 'GSM B'+xyz
            options, pre0+'b0_gsm', 'colors', rgb
            options, pre0+'b0_gsm', 'labels', 'GSM B'+xyz
        endforeach
        
        tplot_save, allvars, filename=datfn
    endif


 
;    vars = [['thd','the','tha','rba']+'_b_gsm',[['thd','the','tha']+'_e_dsl','rba_e_mgse']]
;    vars = ['thd_'+['b_gsm','e_dsl'],'the_'+['b_gsm','e_dsl'],'tha_'+['b_gsm','e_dsl'],'rba_'+['b_gsm','e_mgse']]
;    vars = ['thd_'+['b_gsm','vbulk'],'the_'+['b_gsm','vbulk'],'tha_'+['b_gsm','vbulk'],'rba_'+['b_gsm','vbulk']]
;    options, 'rba_b_gsm', 'yrange', [0,120]
;    tplot, vars, trange=utr0
;    stop


    
end
