;+
; Load RBSP and ASI data.
; This version is different from _2013_0607_load_data and _2013_0607_load_data2
;
; I updated the way auroral images works, so it is now based on MLat/MLon,
; which shrinks the data size and runs faster. It is also much easier to calc
; keogram or ewogram.
;-
;

pro _2013_0607_load_event_info, id, sites=sites
    common constants_and_settings

    probes = ['a','b']
    pres = 'rb'+probes+'_'

    gamma = 5d/3    ; adiabatic constant.
    nboom = 6

    tinfo = {$
        n: 0d, $
        t: 0d, $
        mass: 0d, $
        cs: 0d, $       ; sound speed, km/s, sqrt(\gamma*kT/m).
        placeholder: 0b}

    sc_info = {$
        probe: '', $
        time_range: [0d,0], $       ; time range for averaging quantities.
        data_rate_burst: 0d, $      ; in sec, the time for burst data.
        fpt_mlon_range: [0d,0], $
        fpt_mlat_range: [0d,0], $
        e: tinfo, $
        h: tinfo, $
        o: tinfo, $
        he: tinfo, $
        o_ratio: 0d, $  ; n_o/(n_o+n_h)s.
        n_ratio: 0d, $  ; n_e/(n_o+n_h).
        mhd_rho: 0d, $  ; MHD mass density, in atomic mass*cc.
        mhd_p: 0d, $    ; MHD pressure, in eV*cc.
        r_gsm: [0d,0,0], $  ; averaged GSM position.
        b_gsm: [0d,0,0], $  ; averaged GSM B field.
        r_sm: [0d,0,0], $   ; converted from rgsm.
        bmag: 0d, $         ; |bgsm|.
        va: 0d, $       ; Alfven speed, km/s.
        cs: 0d, $       ; Sound speed, km/s.
        vexb_gsm: [0d,0,0], $   ; ExB velocity in GSM in km/s.
        vexb_fac: [0d,0,0], $   ; ExB velocity in FAC in km/s.
        time_compression: 0d, $ ; UT sec for the initial compression.
        placeholder: 0d}
    rbspa = sc_info
    rbspb = sc_info

    foreach probe, probes, ii do begin
        tinfo = (probe eq 'a')? rbspa: rbspb
        tinfo.probe = probe
        pre0 = pres[ii]
        time_range = (probe eq 'a')? time_range_rbspa: time_range_rbspb
        tinfo.time_range = time_range
        tinfo.data_rate_burst = (probe eq 'a')? 1d/4096: 1d/1024
        all_tags = strlowcase(tag_names(tinfo))

        ; The range of s/c footpoint predicted by models.
        vars = pre0+['fpt_mlon','fpt_mlat']
        tags = 'fpt_'+['mlon_range','mlat_range']
        foreach var, vars, jj do begin
            get_data, var, times, data
            index = where(times ge time_range[0] and times le time_range[1])
            tinfo.(where(all_tags eq tags[jj])) = minmax(data[index,*])
        endforeach

        ; The density and temperature meassured by HOPE for electrons.
        vars = pre0+['n','t']
        foreach var, vars, jj do begin
            get_data, var, times, data
            index = where(times ge time_range[0] and times le time_range[1])
            tinfo.e.(jj) = mean(data[index],/nan)
        endforeach

        ; The density and temperature meassured by HOPE for the ions.
        tags = ['h','o','he']
        pre1s = pre0+['p','o','he']+'_'
        foreach pre1, pre1s, jj do begin
            get_data, pre1+'density', times, data
            index = where(times ge time_range[0] and times le time_range[1])
            tinfo.(where(all_tags eq tags[jj])).n = mean(data[index],/nan)
            get_data, pre1+'t_avg', times, data
            index = where(times ge time_range[0] and times le time_range[1])
            tinfo.(where(all_tags eq tags[jj])).t = mean(data[index],/nan)
        endforeach
        tinfo.e.mass = 0.91e-30/1.67e-27
        tinfo.h.mass = 1d
        tinfo.o.mass = 16d
        tinfo.he.mass = 4d

        ; Assuming quasi-neutrality, neglecting He+.
        tinfo.n_ratio = tinfo.e.n/(tinfo.h.n+tinfo.o.n)   ; to scale ion density to electron density.
        ; The ratio of O+ in the number density, neglecting He+.
        tinfo.o_ratio = tinfo.o.n/(tinfo.h.n+tinfo.o.n)

        ; MHD rho = \sum_s (n_s * m_s)
        tinfo.mhd_rho = (tinfo.h.n*tinfo.h.mass+tinfo.o.n*tinfo.o.mass)*tinfo.n_ratio
        ; MHD P = \sum_s (n_s * T_s)
        tinfo.mhd_p = tinfo.e.n*tinfo.e.t + (tinfo.h.n*tinfo.h.t + tinfo.o.n*tinfo.o.t)*tinfo.n_ratio


        ; The average B GSM, R GSM, ExB velocity in GSM and FAC.
        tags = ['r_gsm','b_gsm','vexb_gsm','vexb_fac']
        vars = pre0+['r_gsm','b_gsm','vexb_gsm','vexb_fac']
        foreach var, vars, jj do begin
            get_data, var, times, data
            index = where(times ge time_range[0] and times le time_range[1], count)
            tinfo.(where(all_tags eq tags[jj])) = total(data[index,*],1)/count
        endforeach

        epoch = stoepoch(mean(time_range),'unix')
        tinfo.r_sm = sgsm2sm(tinfo.r_gsm, epoch)
        tinfo.bmag = snorm(tinfo.b_gsm)


        ; Alfven speed.
        va0 = 22.0d     ; km/s, for B in nT, n in cc, m in atomic mass.
        tinfo.va = va0*tinfo.bmag/sqrt(tinfo.mhd_rho)

        ; Sound speed.
        cs0 = 9.8d      ; km/s, for gamma=1, P in eV*cc and rho in atomic mass*cc, or T in eV and mass in atomic mass.
        tinfo.cs = cs0*sqrt(gamma*tinfo.mhd_p/tinfo.mhd_rho)
        tags = ['e','h','o','he']
        foreach tag, tags do begin
            index = where(all_tags eq tag)
            tinfo.(index).cs = cs0*sqrt(gamma*tinfo.(index).t/tinfo.(index).mass)
        endforeach


        ; Time of the initial compression.
        get_data, pre0+'bmag_smooth', times, bmag
        index = where(times ge time_range_compression[0] and times le time_range_compression[1])
        times = times[index]
        bmag = bmag[index]
        min_bmag = min(bmag, index)
        tinfo.time_compression = times[index]

        if probe eq 'a' then rbspa = tinfo else rbspb = tinfo
    endforeach


    event_info = {$
        id: id, $
        time_range_long: time_range_long, $
        time_range_long_pad: time_range_long_pad, $
        time_range_short: time_range_short, $
        time_range_compression: time_range_compression, $
        probes: probes, $
        pres: pres, $
        data_rate_pos: data_rate_pos, $             ; the data rate for pos, mapping vars.
        models: models, $                           ; a string array of the models.
        model_colors: model_colors, $               ; an array for true-colors.
        model_chosen: model_chosen, $               ; a string for the chosen model.
        dir: dir, $                                 ; a number of tracing direction, trace to the northern hemisphere.
        r0: r0, $                                   ; in Re, for the tracing height for 100 km altitude.
        data_rate_field: data_rate_field, $         ; the data rate for field, pflux vars.
        spin_rate: spin_rate, $                     ; the nominal spin rate in sec.
        boom_lengths: boom_lengths, $               ; the boom lengths in m.
        edot0_min_bmag: edot0_min_bmag, $           ; the min |B| for calc edot0.
        edot0_min_bratio: edot0_min_bratio, $       ; the min ratio of Bw/|B| for calc edot0.
        scale_info: scale_info, $                   ; the structure for calc pflux morlet wavelet.
        filter: filter, $                           ; [2] array, the min&max period in sec.
        use_edot0: use_edot0, $                     ; a boolean, 1: use edot0.
        hope_ion_min_energy: hope_ion_min_energy, $ ; the min energy in eV for the ion moments.
        hope_ele_min_energy: hope_ele_min_energy, $ ; the min energy in eV for the ele moments.
        ions: ions, $                               ; a string array of ion species.
        species: species, $                         ; a string array for the species that have moments available.
        data_rate_asf: data_rate_asf, $             ; the data rate for asf.
        time_range_asf: time_range_asf, $           ; the time range for calc asf.
        time_range_ewo: time_range_ewo, $           ; the time range for calc ewo.
        time_range_keo: time_range_keo, $           ; the time range for calc keo.
        rbspa: rbspa, $
        rbspb: rbspb, $
        sites: sites, $
        nboom: nboom, $
        boom_labels: string(indgen(nboom)+1,format='(I0)'), $
        boom_colors: indgen(nboom), $
        placeholder:1b}

    store_data, id+'_event_info', 0, event_info

end


pro _2013_0607_load_pos, time_range, probes=probes
    common constants_and_settings

    foreach probe, probes, ii do begin
        pre0 = pres[ii]
        spice = sread_rbsp_spice_product(time_range, probe=probe)

        uts = spice.ut_pos
        store_data, pre0+'r_gsm', uts, spice.pos_gsm, limits={ytitle:'(Re)',colors:rgb, labels:xyz}
        store_data, pre0+'mlt', uts, spice.mlt
        store_data, pre0+'mlat', uts, spice.mlat*deg
        store_data, pre0+'lshell', uts, spice.lshell
        store_data, pre0+'dis', uts, snorm(spice.pos_gsm)

        uts = spice.ut_cotran
        quvw2gsm = spice.q_uvw2gsm
        store_data, pre0+'q_uvw2gsm', uts, quvw2gsm

        ; Make sure that these variables are uniform in time.
        times = make_bins(time_range, data_rate_pos)
        vars = pre0+['r_gsm','mlt','mlat','lshell','dis']
        foreach var, vars do begin
            get_data, var, uts, dat
            store_data, var, times, sinterpol(dat,uts, times)
        endforeach
    endforeach
end


pro _2013_0607_load_mapping_info, time_range, probes=probes
    common constants_and_settings

    ; Prepare model parameters.
    foreach model, models do sgeopack_par, time_range, model

    foreach probe, probes, ii do begin
        pre0 = pres[ii]

        pos_var = pre0+'r_gsm'
        if tnames(pos_var) eq '' then _2013_0607_load_pos, time_range, probe=probe

        get_data, pos_var, times, rgsms
        epochs = stoepoch(times, 'unix')
        nrec = n_elements(times)
        nmodel = n_elements(models)
        hfmlats = fltarr(nrec,nmodel)       ; footpoint MLat, in deg.
        hfmlons = fltarr(nrec,nmodel)       ; footpoint MLon, in deg.
        hfmlts = fltarr(nrec,nmodel)        ; footpoint MLT, in hr.
        bmodgsms = fltarr(nrec,3)           ; in-situ B GSM.
        mapcoefs = fltarr(nrec,nmodel)      ; ratio of |B| between in-situ and footpoint.


        foreach model, models, jj do begin
            suf0 = '_'+model
            t89 = (model eq 't89')? 1: 0
            t96 = (model eq 't96')? 1: 0
            t01 = (model eq 't01')? 1: 0
            t04s = (model eq 't04s')? 1: 0
            storm = (model eq 't04s')? 1: 0

            get_data, model+'_par', tmp, pars
            pars = sinterpol(pars, tmp, times)

            for i=0, nrec-1 do begin
                epoch = epochs[i]
                geopack_epoch, epoch, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
                geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date

                ; In-situ position, in Re.
                xp = rgsms[i,0]
                yp = rgsms[i,1]
                zp = rgsms[i,2]

                ; Map to the ionosphere.
                ; [xyz]f is the footpoint position.
                par = reform(pars[i,*])
                geopack_trace, xp,yp,zp, dir, par, xf,yf,zf, r0=r0, /igrf, t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm
                ; Convert footpoint from gsm to mag.
                geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_mag
                hfmlats[i,jj] = asin(tzf/r0)*deg
                hfmlons[i,jj] = atan(tyf,txf)*deg
                hfmlts[i,jj] = mlon2mlt(hfmlons[i,jj], times[jj])

                ; Model B GSM, in-situ.
                ; Internal B field.
                geopack_igrf_gsm, xp,yp,zp, bxp,byp,bzp
                ; External B field.
                case model of
                    't89': geopack_t89, par, xp,yp,zp, tbx,tby,tbz
                    't96': geopack_t96, par, xp,yp,zp, tbx,tby,tbz
                    't01': geopack_t01, par, xp,yp,zp, tbx,tby,tbz
                    't01s': geopack_t01, par, xp,yp,zp, tbx,tby,tbz /storm
                    't04s': geopack_ts04, par, xp,yp,zp, tbx,tby,tbz
                    else: begin tbx = 0 & tby = 0 & tbz = 0 & end
                endcase
                bxp+= tbx & byp+= tby & bzp+= tbz
                bmodgsms[i,*] = [bxp,byp,bzp]

                ; Model B GSM, footpoint.
                geopack_igrf_gsm, xf,yf,zf, bxf,byf,bzf

                ; Mapping coefficent.
                mapcoefs[i,jj] = snorm([bxf,byf,bzf])/snorm([bxp,byp,bzp])
            endfor
            store_data, pre0+'bmod_gsm'+suf0, times, bmodgsms, limits=$
                {ytitle:'(nT)', labels:'Bmod GSM'+xyz, colors:rgb}
        endforeach


        ; Save to memory.
        store_data, pre0+'fpt_mlat', times, hfmlats, limits={ytitle:'(deg)', labels:models, colors:model_colors, ynozero:1}
        store_data, pre0+'fpt_mlon', times, hfmlons, limits={ytitle:'(deg)', labels:models, colors:model_colors, ynozero:1}
        store_data, pre0+'fpt_mlt', times, hfmlts, limits={ytitle:'(hr)', labels:models, colors:model_colors, ynozero:1}
        store_data, pre0+'cmap_model', times, mapcoefs, limits={ytitle:'#', labels:models, colors:model_colors, ynozero:1}

        ; Fix steps in MLat.
        tmp_vars = pre0+'fpt_mlat_'+models
        mlat_var = pre0+'fpt_mlat'
        get_data, mlat_var, limits=limits
        stplot_split, mlat_var, newnames=tmp_vars
        foreach tvar, tmp_vars do fix_model_mlat, tvar, to=tvar
        stplot_merge, tmp_vars, newname=mlat_var, /delete
        store_data, mlat_var, limits=limits
    endforeach
end


pro _2013_0607_load_field, time_range, probes=probes
    common constants_and_settings

    foreach probe, probes, ii do begin
        pre0 = pres[ii]

        ; Uniform time for E/B.
        times = make_bins(time_range, data_rate_field)


    ;---Load Q.
        q_var = pre0+'q_uvw2gsm'
        if tnames(q_var) eq '' then _2013_0607_load_pos, time_range, probes=probes
        get_data, q_var, uts, quvw2gsm
        quvw2gsm = qslerp(quvw2gsm, uts, times)
        muvw2gsm = transpose(qtom(quvw2gsm))

    ;---Load B.
        rbsp_read_bfield, time_range, probe=probe
        bvar = 'rbsp'+probe+'_b_gsm'
        get_data, bvar, tuts, bgsm
        bgsm = sinterpol(bgsm, tuts, times)
        store_data, pre0+'b_gsm', times, bgsm, limits={ytitle:'(nT)', labels:'GSM B'+xyz, colors:rgb}
        store_data, bvar, /delete

        ; Calc B0 and dB.
        b0gsm = bgsm
        ndim = (size(b0gsm,/dimensions))[1]
        for i=0,ndim-1 do b0gsm[*,i] = scalcbg(bgsm[*,i])
        store_data, pre0+'b0_gsm', times, b0gsm, limits={ytitle:'(nT)', labels:'GSM B0'+xyz, colors:rgb}
        store_data, pre0+'db_gsm', times, bgsm-b0gsm, limits={ytitle:'(nT)', labels:'GSM dB'+xyz, colors:rgb}

        ; Calc |B| and its smoothed version.
        bmag = snorm(bgsm)
        store_data, pre0+'bmag', times, bmag, limits = {ytitle:'(nT)', labels:'|B|'}
        ; calc |B| smoothed.
        order = 1
        width = round(spin_rate/data_rate_field)
        for i=0,order do bmag = smooth(bmag, width, /edge_truncate)
        store_data, pre0+'bmag_smooth', times, bmag, limits = {ytitle:'(nT)', labels:'|B| smth'}

        ; Modifies the mapping coefficient using the measured B field.
        tvar = pre0+'cmap_model'
        if tnames(tvar) eq '' then _2013_0607_load_mapping_info, time_range, probes=probes
        get_data, tvar, uts, cmaps, limits=limits
        get_data, pre0+'bmag_smooth', times, bmag
        bmag = interpol(bmag, times, uts)
        foreach model, models, jj do begin
            get_data, pre0+'bmod_gsm_'+model, uts, bmodgsm
            cmaps[*,jj] *= snorm(bmodgsm)/bmag
        endforeach
        store_data, pre0+'cmap', uts, cmaps, limits=limits


    ;---Load E.
        ; Load V[1-6].
        dat = sread_rbsp_efw_l2(time_range, probes=probe, type='vsvy')
        tuts = sfmepoch(dat.epoch, 'unix')
        vsvy = sinterpol(dat.vsvy, tuts, times)

        ; Calc Vsc.
        tvar = pre0+'vsc'
        vsc = mean(vsvy[*,0:1], dimension=2)
        store_data, tvar, times, vsc, limits={ytitle:'(V)', labels:'Vsc'}

        ; Calc E[uvw].
        tvar = pre0+'e_uvw'
        eu = (vsvy[*,0]-vsvy[*,1])/boom_lengths[0]*1e3   ; V -> V/m -> mV/m.
        ev = (vsvy[*,2]-vsvy[*,3])/boom_lengths[1]*1e3
        ew = dblarr(n_elements(eu))

        ; Remove dc-offset.
        nspin = 1
        width = nspin*spin_rate/data_rate_field
        eu = eu-smooth(eu, width, /edge_truncate, /nan)
        ev = ev-smooth(ev, width, /edge_truncate, /nan)
        store_data, tvar, times, [[eu],[ev],[ew]], limits={ytitle:'(mV/m)', labels:'E'+uvw, colors:rgb}

        ; Calc E GSM.
        get_data, pre0+'e_uvw', times, euvw
        euvw[*,2] = 0
        ex = euvw[*,0]*muvw2gsm[0,0,*] + euvw[*,1]*muvw2gsm[1,0,*] + euvw[*,2]*muvw2gsm[2,0,*]
        ey = euvw[*,0]*muvw2gsm[0,1,*] + euvw[*,1]*muvw2gsm[1,1,*] + euvw[*,2]*muvw2gsm[2,1,*]
        ez = euvw[*,0]*muvw2gsm[0,2,*] + euvw[*,1]*muvw2gsm[1,2,*] + euvw[*,2]*muvw2gsm[2,2,*]
        tvar = pre0+'e_gsm'
        store_data, tvar, times, [[ex],[ey],[ez]], limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM E'+xyz}

        ; Calc E dot0 GSM.
        get_data, pre0+'b0_gsm', times, b0gsm
        bx = b0gsm[*,0] & by = b0gsm[*,1] & bz = b0gsm[*,2]
        bu = bx*muvw2gsm[0,0,*] + by*muvw2gsm[0,1,*] + bz*muvw2gsm[0,2,*]
        bv = bx*muvw2gsm[1,0,*] + by*muvw2gsm[1,1,*] + bz*muvw2gsm[1,2,*]
        bw = bx*muvw2gsm[2,0,*] + by*muvw2gsm[2,1,*] + bz*muvw2gsm[2,2,*]

        buvw = [[bu],[bv],[bw]]
        bmag = snorm(buvw)
        euvw[*,2] = -(euvw[*,0]*buvw[*,0]+euvw[*,1]*buvw[*,1])/buvw[*,2]
        idx = where(bmag le edot0_min_bmag or abs(buvw[*,2])/bmag le edot0_min_bratio, cnt)
        if cnt ne 0 then euvw[idx,2] = !values.d_nan

        ex = euvw[*,0]*muvw2gsm[0,0,*] + euvw[*,1]*muvw2gsm[1,0,*] + euvw[*,2]*muvw2gsm[2,0,*]
        ey = euvw[*,0]*muvw2gsm[0,1,*] + euvw[*,1]*muvw2gsm[1,1,*] + euvw[*,2]*muvw2gsm[2,1,*]
        ez = euvw[*,0]*muvw2gsm[0,2,*] + euvw[*,1]*muvw2gsm[1,2,*] + euvw[*,2]*muvw2gsm[2,2,*]

        tvar = pre0+'edot0_gsm'
        store_data, tvar, times, [[ex],[ey],[ez]], limits={ytitle:'(mV/m)', colors:rgb, labels:'GSM Edot0'+xyz}
    endforeach

end


pro _2013_0607_load_pflux, time_range, probes=probes
    common constants_and_settings

    foreach probe, probes, ii do begin
        pre0 = pres[ii]

        ; Convert E/B field to fac.
        get_data, pre0+'b0_gsm', uts, b0gsm
        get_data, pre0+'r_gsm', tuts, rgsm
        rgsm = sinterpol(rgsm, tuts, uts)
        rhat = sunitvec(rgsm)
        bhat = sunitvec(b0gsm)
        phat = sunitvec(scross(rhat,bhat))
        vhat = scross(bhat,phat)

        tvar = use_edot0? pre0+'edot0_gsm': pre0+'e_gsm'
        get_data, tvar, uts, egsm
        defac = [[sdot(egsm,bhat)],[sdot(egsm,phat)],[sdot(egsm,vhat)]]
        store_data, pre0+'de_fac', uts, defac, limits={ytitle:'(mV/m)', labels:'FAC dE'+fac, colors:rgb}

        tvar = pre0+'db_gsm'
        get_data, tvar, uts, dbgsm
        dbfac = [[sdot(dbgsm,bhat)],[sdot(dbgsm,phat)],[sdot(dbgsm,vhat)]]
        store_data, pre0+'db_fac', uts, dbfac, limits={ytitle:'(nT)', labels:'FAC dB'+fac, colors:rgb}

        ; Calculate the instataneous pflux.
        stplot_calc_pflux_mor, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', scaleinfo=scale_info, filter=filter
        tvar = pre0+'pf_fac'
        options, tvar, 'colors', rgb
        options, tvar, 'labels', 'FAC S'+fac


        ; Calc ExB velocity.
        bmag = snorm(b0gsm)
        coef = 1d3/bmag^2   ; convert to km/s for E in mV/m and B in nT.
        ndim = (size(b0gsm,/dimensions))[1]

        tvar = pre0+'vexb_gsm'
        vexb = vec_cross(egsm, b0gsm)
        for jj=0, ndim-1 do vexb[*,jj] *= coef
        store_data, tvar, uts, vexb, limits={ytitle:'(km/s)', labels:'GSM ExB V'+xyz, colors:rgb}

        tvar = pre0+'vexb_fac'
        bmag = snorm(b0gsm)
        b0fac = b0gsm & b0fac[*,1:2] = 0 & b0fac[*,0] = bmag
        vexb = vec_cross(defac, b0fac)
        for jj=0, ndim-1 do vexb[*,jj] *= coef
        store_data, tvar, uts, vexb, limits={ytitle:'(km/s)', labels:'FAC ExB V'+fac, colors:rgb}
    endforeach

end


pro _2013_0607_load_hope, time_range, probes=probes
    common constants_and_settings

    foreach probe, probes, ii do begin
        pre0 =  pres[ii]

    ;---Load Mageis electron.
        mageis = sread_rbsp_mageis_l3(time_range, probes=probe)
        times = sfmepoch(mageis.epoch,'unix')
        fedu = total(mageis.fedu,3,/nan)
        store_data, pre0+'e_en_mageis', times, fedu, mageis.fedu_energy, limits = $
            {spec:1, no_interp:1, ylog:1, zlog:1, yrange:[4e1,4e3], ystyle:1, $
            ytitle:'MAGEIS!CElectron!CEenergy!C(keV)', $
            zticks:3, zrange:[1e-2,1e5], $
            ztitle:'(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'}

    ;---Load HOPE electron, proton, oxygen data.
        vars = ['PITCH_ANGLE',$
            'Epoch_Ele','HOPE_ENERGY_Ele',$
            'Epoch_Ion','HOPE_ENERGY_Ion',$
            'FEDU','FPDU','FODU']

        hopel3 = sread_rbsp_hope_l3(time_range, probes=probe, vars=vars)
        pas = hopel3.pitch_angle

        ; ion spec.
        uts = sfmepoch(hopel3.epoch_ion,'unix')
        nrec = n_elements(uts)
        ens = hopel3.hope_energy_ion

        dat = hopel3.fodu
        type = 'o'
        zr = [1e4,1e7]
        idx = where(finite(dat,/nan),cnt)
        if cnt ne 0 then dat[idx] = 0
        for i = 0, nrec-1 do begin
            idx = where(ens[i,*] le hope_ion_min_energy, cnt)
            if cnt ne 0 then dat[i,idx,*] = 0
        endfor

        paspec = total(dat, 2)
        idx = where(paspec eq 0)
        paspec[idx] = 7

        tvar = pre0+type+'_pa'
        store_data, tvar, uts, paspec, pas, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
            yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
            ytitle:'Pitch Angle (deg)'}

        enspec = total(dat, 3)
        idx = where(enspec eq 0)
        enspec[idx] = 7

        tvar = pre0+type+'_en'
        store_data, tvar, uts, enspec, ens, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr, $
            yrange:[3,3e4], ystyle:1, ylog:1, $
            ytitle:'Energy (eV)'}


        dat = hopel3.fpdu
        type = 'h'
        zr = [1e4,1e7]
        idx = where(finite(dat,/nan),cnt)
        if cnt ne 0 then dat[idx] = 0
        for i = 0, nrec-1 do begin
            idx = where(ens[i,*] le hope_ion_min_energy, cnt)
            if cnt ne 0 then dat[i,idx,*] = 0
        endfor

        paspec = total(dat, 2)
        idx = where(paspec eq 0)
        paspec[idx] = 7

        tvar = pre0+type+'_pa'
        store_data, tvar, uts, paspec, pas, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
            yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
            ytitle:'Pitch Angle (deg)'}

        enspec = total(dat, 3)
        idx = where(enspec eq 0)
        enspec[idx] = 7

        tvar = pre0+type+'_en'
        store_data, tvar, uts, enspec, ens, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr, $
            yrange:[30,3e4], ystyle:1, ylog:1, $
            ytitle:'Energy (eV)'}


        ; electron.
        uts = sfmepoch(hopel3.epoch_ele,'unix')
        nrec = n_elements(uts)
        ens = hopel3.hope_energy_ele

        dat = hopel3.fedu
        type = 'e'
        zr = [1e7,1e10]
        idx = where(finite(dat,/nan),cnt)
        if cnt ne 0 then dat[idx] = 0
        for i = 0, nrec-1 do begin
            idx = where(ens[i,*] le hope_ele_min_energy, cnt)
            if cnt ne 0 then dat[i,idx,*] = 0
        endfor

        paspec = total(dat, 2)
        idx = where(paspec eq 0)
        paspec[idx] = 7

        tvar = pre0+type+'_pa'
        store_data, tvar, uts, paspec, pas, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr*10, $
            yrange:[0,180], ystyle:1, ytickv:[0,90,180], yticks:2, yminor:9, $
            ytitle:'Pitch Angle (deg)'}

        enspec = total(dat, 3)
        idx = where(enspec eq 0)
        enspec[idx] = 7

        tvar = pre0+type+'_en'
        store_data, tvar, uts, enspec, ens, limits = $
            {spec:1, no_interp:1, zlog:1, zrange:zr, $
            yrange:[200,2e4], ylog:1, ystyle:1, $
            ytitle:'Energy (eV)'}

        vars = pre0+['h_en','h_pa','e_en','e_pa','o_en','o_pa']
        options, vars, 'spec', 1
        options, vars, 'no_interp', 1
        options, vars, 'zlog', 1
        options, vars, 'ztitle', '(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'

        vars = pre0+['h_pa','e_pa','o_pa']
        options, vars, 'yrange', [0,180]
        options, vars, 'ystyle', 1

        tvar = pre0+'e_en'
        options, tvar, 'ytitle', 'e- energy!C(eV)'
        tvar = pre0+'h_en'
        options, tvar, 'ytitle', 'H+ energy!C(eV)'

        tvar = pre0+'e_pa'
        options, tvar, 'ytitle', 'e- pitch!C(deg)'
        tvar = pre0+'h_pa'
        options, tvar, 'ytitle', 'H+ pitch!C(deg)'
        tvar = pre0+'o_pa'
        options, tvar, 'ytitle', 'O+ pitch!C(deg)'


    ;---Load HOPE even moments.
        hopemom = sread_rbsp_hope_l3(time_range, probes=probe, type='mom', $
            vars = ['Epoch_Ion','Dens_o_30','Dens_p_30'])
        if size(hopemom,/type) ne 8 then return
        uts = sfmepoch(hopemom.epoch_ion,'unix')
        store_data, pre0+'ni', uts, [[hopemom.dens_p_30],[hopemom.dens_o_30]], $
            limits = {ytitle:'n!Di!N!C(cm!U-3!N)', ylog:1, constant:1, labels:$
        ['n!DH!N>30eV','n!DO!N>30eV'], colors:[black,red]}


        hopemom = sread_rbsp_hope_l3(time_range, probes=probe, type='mom')
        if size(hopemom,/type) ne 8 then return
        uts = sfmepoch(hopemom.epoch_ele,'unix')
        store_data, pre0+'n', uts, hopemom.dens_e_200, $
            limits = {ytitle:'n!De!N!C(cm!U-3!N)', ylog:1, constant:1, labels:'n!De!N>200eV'}
        store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
            limits = {ytitle:'T!De!N!C(eV)', ylog:1, colors:[red,black], labels:['Tperp','Tpara'], constant:1000}
        tvar = pre0+'t'
        get_data, tvar, uts, dat
        idx = where(dat eq 1e20, cnt)
        if cnt ne 0 then dat[idx] = !values.d_nan
        store_data, tvar, uts, dat


    ;---Fit Vsc to density.
        get_data, pre0+'vsc', uts, vsc
        get_data, pre0+'n', tuts, hopen
        idx = where(finite(hopen))
        hopen = sinterpol(hopen[idx], tuts[idx], uts, /nan)

        ; hopen = a*exp(vsc*b)
        ; alog(hopen) = alog(a)+vsc*b
        tx0 = vsc
        ty0 = alog(hopen)

        a = regress(tx0, ty0, const = b)
        ty1 = a[0]*vsc+b[0]

        ybar = mean(ty0)
        sstot = total((ty0-ybar)^2)
        ssres = total((ty0-ty1)^2)

        r2 = 1-ssres/sstot
        fit_param = {a:a, b:b, r2:r2}

        efwn = exp(ty1)
        store_data, pre0+'n_efw', uts, efwn, fit_param, limits = $
            {ytitle:'Density!C(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
            labels:'n!De,EFW Vsc'}
        store_data, pre0+'n_combine', uts, [[efwn],[hopen]], limits = $
            {ytitle:'(cm!U-3!N)', ylog:1, yrange:[0.02,2], ystyle:1, $
            labels:['n!De,EFW Vsc!N', 'n!De,HOPE>'+string(hope_ele_min_energy,format='(I0)')+'eV'], colors:[red,black]}

    endforeach
end


pro _2013_0607_load_moments, time_range, probes=probes
    common constants_and_settings

    momvar0s = ['density','t_avg','vbulk','number_flux','energy_flux','enthalpy']
    units = ['(cm!U-3!N)','(eV)','(km/s)','(#/cm!U-2-s)','(mW/m!U2!N)','(mW/m!U2!N)']

    momvar1s = ['vbulk','number_flux','energy_flux','enthalpy']     ; vectors.

    foreach probe, probes, ii do begin
        pre0 =  pres[ii]

        ; Calculate the moments.
        ofn = shomedir()+'/'+pre0+'moments.cdf'
        rbsp_gen_hope_moments, probe=probe, time_range, ion_energy_range=[hope_ion_min_energy,1e4], filename=ofn
        moms = sread_rbsp_hope_moments(filename=ofn)
        file_delete, ofn

        ; Combine times of electrons and ions.
        uts = [moms.ut_ion,moms.ut_ele]
        uts = uts[sort(uts)]

        vnames = strlowcase(tag_names(moms))
        foreach tspecie, species do begin
            pre1 = tspecie+'_'
            utname = (tspecie eq 'e')? 'ut_ele': 'ut_ion'
            tuts = moms.(where(vnames eq utname))
            foreach tvar, momvar0s, i do begin
                dat = moms.(where(vnames eq pre1+tvar))
                dat = sinterpol(dat, tuts, uts)
                ttvar = pre0+pre1+tvar
                store_data, ttvar, uts, dat, limits={ytitle:units[i],labels:tvar}
            endforeach

            vecvars = pre0+pre1+momvar1s
            options, vecvars, 'colors', rgb

            ; Calculate FAC unit vectors.
            get_data, pre0+'b0_gsm', tuts, b0gsm
            b0gsm = sinterpol(b0gsm, tuts, uts)
            get_data, pre0+'r_gsm', tuts, rgsm
            rgsm = sinterpol(rgsm, tuts, uts)
            rhat = sunitvec(rgsm)
            bhat = sunitvec(b0gsm)
            phat = sunitvec(scross(rhat,bhat))
            vhat = scross(bhat,phat)

            foreach tvar, vecvars do begin
                get_data, tvar, uts, dat, limits=lim
                store_data, tvar+'_gsm', uts, dat, limits=lim
                options, tvar+'_gsm', 'labels', 'GSM '+xyz
                options, tvar+'_gsm', 'colors', rgb

                dat = [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]]
                store_data, tvar+'_fac', uts, dat, limits=lim
                options, tvar+'_fac', 'labels', 'FAC '+fac
                options, tvar+'_fac', 'colors', rgb
            endforeach
        endforeach
    endforeach
end


pro _2013_0607_load_mlonimg, time_range, sites=sites
    common constants_and_settings

    bg_count = 3000d
    scale_factor = 60d/(4000-3000)

    foreach site, sites do begin
        themis_read_mlonimg_per_site, time_range, site=site

        mlonimg_var = 'thg_'+site+'_asf_mlonimg'
        get_data, mlonimg_var, times, mlonimgs
;        ntime = n_elements(times)
;        min_counts = fltarr(ntime)
;        max_counts = fltarr(ntime)
;        avg_counts = fltarr(ntime)
;        med_counts = fltarr(ntime)
;        for ii=0, ntime-1 do begin
;            timg = reform(mlonimgs[ii,*,*])
;            index = where(timg ne 0, count)
;            if count eq 0 then continue
;            tdat = timg[index]
;            min_counts[ii] = min(tdat)
;            max_counts[ii] = max(tdat)
;            avg_counts[ii] = mean(tdat)
;            med_counts[ii] = median(tdat)
;        endfor
;        stop

        mlonimgs = (mlonimgs-bg_count)*scale_factor>0
        store_data, mlonimg_var, times, mlonimgs
    endforeach

end


pro _2013_0607_load_data3, filename=filename, $
    event_info=event_info, $    ; output, return the structure of event info.
    probes=probes, $            ; Default is ['a','b'].
    reload_pos=reload_pos, $
    reload_mapping_info=reload_mapping_info, $
    reload_field=reload_field, $
    reload_pflux=reload_pflux, $
    reload_hope=reload_hope, $
    reload_moments=reload_moments, $
    reload_mlonimg=reload_mlonimg, $
    reload_event_info=reload_event_info, $
    placeholder=placeholder


    common constants_and_settings, $
    ; constants.
        deg,rad,re,re1,secofday, $
        fac,xyz,uvw,rgb,red,black,white, $
    ; settings.
        time_range_long,time_range_long_pad,time_range_short, pres, $
        time_range_compression, $
        data_rate_pos,models,model_colors,model_chosen,dir,r0, $
        data_rate_field,spin_rate,boom_lengths,edot0_min_bmag,edot0_min_bratio, $
        scale_info,filter,use_edot0, $
        hope_ion_min_energy,hope_ele_min_energy, $
        ions, species, $
        data_rate_asf,time_range_asf,time_range_ewo,time_range_keo, $
        mlonimg_mlon_step, mlonimg_mlat_step, $
        time_range_rbspa, time_range_rbspb

;---Constants.
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    secofday = 86400d

    fac = ['b','w','n'] ; para, west, north.
    xyz = ['x','y','z']
    uvw = ['u','v','w']
    rgb = sgcolor(['red','green','blue'])
    black = sgcolor('black')
    red = sgcolor('red')
    white = sgcolor('white')


;---Plot settings.
    tplot_options, 'labflag', -1
    tplot_options, 'yticklen', -0.01
    tplot_options, 'xticklen', -0.02
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 3
    tplot_options, 'num_lab_min', 6
    tplot_options, 'zcharsize', 0.8


;---Settings.
    update_datfn = 0
    time_range_long = time_double(['2013-06-07/04:45','2013-06-07/05:15'])  ; long time range, used for overview plots.
    time_range_long_pad = time_range_long+[-1,1]*300                        ; 5 min pad time, used for loading data.
    time_range_short = time_double(['2013-06-07/04:52','2013-06-07/05:02'])  ; short time range, used for zoom-in plots.
    time_range_rbspa = time_double(['2013-06-07/04:54','2013-06-07/05:00'])  ; the time range for calc average values.
    time_range_rbspb = time_double(['2013-06-07/04:55','2013-06-07/05:01'])  ; the time range for calc average values.
    time_range_compression = time_double(['2013-06-07/04:54','2013-06-07/04:56'])   ; the time range for finding the foot of initial compression.


;---Check inputs.
    ; check datfn.
    load = 0
    if n_elements(filename) eq 0 then load = 1 else if file_test(filename) eq 0 then load = 1
    if load eq 1 then begin
        rootdir = googledir()+'/works/works/rbsp_psbl1'
        datdir = rootdir+'/data'
        if file_test(datdir,/directory) eq 0 then file_mkdir, datdir
        datfn = datdir+'/2013_0607_all_data_v3.tplot'
    endif else datfn = filename
    if file_test(datfn) eq 1 then tplot_restore, filename=datfn

    ; check probes.
    if n_elements(probes) eq 0 then probes = ['a','b']



;---Load position.
    pres = 'rb'+probes+'_'
    data_rate_pos = 60d             ; common data rate for pos, mapping.

    load = keyword_set(reload_pos)
    tvar = ['r_gsm','mlt','mlat','lshell','dis','q_uvw2gsm']
    pos_vars = []
    foreach pre0, pres do pos_vars = [pos_vars,pre0+tvar]
    foreach tvar, pos_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_pos, time_range_long_pad, probes=probes
        update_datfn = 1
    endif


;---Load mapping info.
    r0 = 100d/re+1
    models = ['t89','t96','t01','t04s']
    model_colors = [6,4,3,2]
    model_chosen = 't01'
    dir = -1    ; map to the northern hemisphere.

    load = keyword_set(reload_mapping_info)
    tvar = ['fpt_'+['mlt','mlat','mlon'],'bmod_gsm_'+models,'cmap_model']
    map_vars = models+'_par'
    foreach pre0, pres do map_vars = [map_vars,pre0+tvar]
    foreach tvar, map_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_mapping_info, time_range_long_pad, probes=probes
        update_datfn = 1
    endif


;---Load E/B field.
    data_rate_field = 1d/16         ; common data rate for fields, pflux.
    spin_rate = 12d                 ; sec.
    boom_lengths = [100d,100,12]    ; in m, twice of boom lengths.
    edot0_min_bmag = 5d ; nT.
    edot0_min_bratio = 0.25         ; percent.

    load = keyword_set(reload_field)
    tvar = ['b_gsm','b0_gsm','db_gsm','bmag','bmag_smooth','e_gsm','edot0_gsm','cmap','vsc']
    field_vars = []
    foreach pre0, pres do field_vars = [field_vars,pre0+tvar]
    foreach tvar, field_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_field, time_range_long_pad, probes=probes
        update_datfn = 1
    endif


;---Load pflux.
    scale_info = {s0:0.25d, s1:1200d, dj:1d/8, ns:0d}   ; for pflux calc.
    filter = [0.25,100]             ; in sec.
    use_edot0 = 1

    load = keyword_set(reload_pflux)
    tvar = [['de','db','pf']+'_fac','pf_fac_mor_spec_'+['1','2','3'],'vexb_'+['gsm','fac']]
    pflux_vars = []
    foreach pre0, pres do pflux_vars = [pflux_vars,pre0+tvar]
    foreach tvar, pflux_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_pflux, time_range_long, probes=probes
        update_datfn = 1
    endif


;---Load HOPE.
    hope_ion_min_energy = 30d       ; eV.
    hope_ele_min_energy = 200d      ; eV.

    load = keyword_set(reload_hope)
    tvar = ['e_en_mageis','e_en','e_pa','h_en','h_pa','o_en','o_pa', $
        'ni','n','t','n_efw','n_combine']
    hope_vars = []
    foreach pre0, pres do hope_vars = [hope_vars,pre0+tvar]
    foreach tvar, hope_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_hope, time_range_long_pad, probes=probes
        update_datfn = 1
    endif


;---Load moments.
    ions = ['p','o','he']
    species = ['e',ions]

    load = keyword_set(reload_moments)
    tvar = ['vbulk','number_flux','energy_flux','enthalpy']
    tvar = ['density','t_avg',tvar+'_gsm',tvar+'_fac']
    mom_vars = []
    foreach pre0, pres do foreach type, species do mom_vars = [mom_vars,pre0+type+'_'+tvar]
    foreach tvar, mom_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_moments, time_range_long_pad, probes=probes
        update_datfn = 1
    endif


;---Load auroral images.
    time_range_asf = time_double(['2013-06-07/04:45','2013-06-07/05:05'])  ; short time range.
    time_range_ewo = time_double(['2013-06-07/04:54','2013-06-07/04:57'])  ; ewo time range.
    time_range_keo = time_double(['2013-06-07/04:45','2013-06-07/05:05'])  ; keo time range.
    data_rate_asf = 3d              ; common data rate for asf.
    sites = 'pina'

    load = keyword_set(reload_mlonimg)
    mlonimg_vars = []
    tvar = ['mlonimg','new_image_size','new_mlon_bins','new_mlat_bins','new_elevs']
    foreach site, sites do mlonimg_vars = [mlonimg_vars,'thg_'+site+'_asf_'+tvar]
    foreach tvar, mlonimg_vars do if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_mlonimg, time_range_asf, sites=sites
        update_datfn = 1
    endif


;---Load event info.
    load = keyword_set(reload_event_info)
    id = time_string(time_range_long[0],tformat='YYYY_MMDD')
    tvar = id+'_event_info'
    if tnames(tvar) eq '' then load = 1
    if load then begin
        _2013_0607_load_event_info, id, sites=sites
        update_datfn = 1
    endif
    event_info_var = id+'_event_info'
    event_info = get_var_data(event_info_var)


;---Save changes.
    if update_datfn then begin
        vars_to_save = [pos_vars, map_vars, field_vars, pflux_vars, hope_vars, mom_vars, mlonimg_vars, event_info_var]
        lprmsg, 'Saving data to '+datfn+' ...'
        tplot_save, vars_to_save, filename=datfn
    endif

end

_2013_0607_load_data3, reload_event_info=1, reload_mapping_info=1
end
