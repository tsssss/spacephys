
pro rbsp_efw_plot_thm_ast, tr, model, scl, probe = probe

    ; prepare.
    if ~keyword_set(probe) then probe = ['a'] & nprob = n_elements(probe)
    if n_elements(model) eq 0 then model = 't89'
    re = 6374D & re1 = 1D/re
    
    ; initiate rbsp efw, emfisis, cdf_leap_second, istp.
    ; the trailing '/' is ANNOYING.
    rbsp_efw_init, local_data_dir = spreproot('rbsp')
    cdf_leap_second_init
    rbsp_emfisis_init, local_data_dir = spreproot('rbsp')
    istp_init
    timespan, tr[0], tr[1]-tr[0], /second
    rbsp_load_spice_kernels
    thm_init, local_data_dir = spreproot('themis')
    
    ; **** load rbsp position in gsm.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        rbsp_load_spice_state, probe = probe[i], $
            coord = 'gsm', /no_spice_load
        get_data, pre0+'state_pos_gsm', data = tmp
        tmp.y *= re1        ; from km to re.
        store_data, pre0+'pos_gsm', data = tmp
        ; delete vars.
        store_data, pre0+'state_*', /delete
    endfor
    
    ; **** load rbsp e and b field.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        ; ** load e in mgse.
        rbsp_load_efw_waveform_l2, probe = probe[i]
        get_data, pre0+'efw_e-spinfit-mgse_efield_spinfit_mgse', data = tmp
        tmp.y[*,0] = 0      ; throw e_spin.
        store_data, pre0+'e_mgse', data = tmp
        ; convert from mgse to gse.
        rbsp_mgse2gse, pre0+'e_mgse', newname = pre0+'e_gse', $
            probe = probe[i], /no_spice_load
        ; delete vars.
        store_data, pre0+'efw_*', /delete
        
        ; ** load b in gse.
        rbsp_load_emfisis, probe = probe[i], coord = 'gse'
        get_data, pre0+'emfisis_l3_4sec_gse_Mag', data = tmp
        store_data, pre0+'b_gse', data = tmp
        ; delete vars.
        store_data, pre0+'emfisis_*', /delete
    endfor
    
    ; **** get model b field.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        get_data, pre0+'b_gse', data = tmp
        bmod = tmp.y
        len = 1200D/sdatarate(tmp.x)    ; 20 min.
        for j = 0, 2 do bmod[*,j] = smooth(bmod[*,j], len)
        store_data, pre0+'bmod_gse', data = {x:tmp.x, y:bmod}
        store_data, pre0+'db_gse', data = {x:tmp.x, y:tmp.y-bmod}
    endfor
    
    ; *** decompose e field.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        get_data, pre0+'e_gse', data = de
        get_data, pre0+'bmod_gse', data = tmp
        bmod = sinterpol(tmp.y, tmp.x, de.x)
        bhat = sunitvec(bmod)
        
        p = atan(bhat[*,1],bhat[*,0])
        cosp = cos(p) & sint = bhat[*,2]
        cost = bhat[*,0]/cosp & sinp = bhat[*,1]/cost
        
        tmp = de.y
        ex = cost*cosp*tmp[*,0] + cost*sinp*tmp[*,1] - sint*tmp[*,2]
        ey = -sinp*tmp[*,0] + cosp*tmp[*,1]
        ez = sint*sinp*tmp[*,0] + sint*cosp*tmp[*,2] + cost*tmp[*,2]
        store_data, pre0+'de_fac', data = {x:de.x, y:[[ex],[ey],[ez]]}
    endfor
    
    ; **** get footprint.
    ; prepare mapping.
    dir = -1        ; always north.
    r0 = 1+110*re1  ; 110 km altitude.
    sgeopack_par, tr, model, /delete  ; get tplot var <model>_par.
    t89 = 0 & t96 = 0 & t01 = 0
    case model of
        't89': t89 = 1
        't96': t96 = 1
        't01': t01 = 1
    endcase
    ; map pos to fpt, convert to mag.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        get_data, pre0+'pos_gsm', data = tmp
        uts = tmp.x & ets = 1000*uts+62167219200000D
        pos0 = tmp.y & pos1 = pos0      ; in re.
        ; interpolate par.
        if model ne 't89' then begin
            get_data, model+'_par', data = tmp
            pars = sinterpol(tmp.y, tmp.x, uts)
        endif
        ; loop for each time.
        for j = 0, n_elements(uts)-1 do begin
            ; set geopack.
            geopack_epoch, ets[j], yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
            geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001D, /date
            ; pos in gsm, which is the mapping coord.
            x0 = pos0[j,0] & y0 = pos0[j,1] & z0 = pos0[j,2]
            if model ne 't89' then par = reform(pars[j,*]) else par = 2
            geopack_trace, x0, y0, z0, dir, par, xf, yf, zf, $
                epoch = ets[i], /refine, /ionosphere, $
                t89 = t89, t96 = t96, t01 = t01
            ; conver from gsm to geo.
            geopack_conv_coord, xf, yf, zf, /from_gsm, $
                x1, y1, z1, /to_geo
            pos1[j,*] = [x1,y1,z1]
        endfor
        rtod = 180D/!dpi
        mlat = asin(pos1[*,2]*(1/r0))*rtod
        mlon = atan(pos1[*,1],pos1[*,0])*rtod
        store_data, pre0+'fpt_geo', data = {x:uts, y:pos1}
        store_data, pre0+'fpt_lonlat', data = {x:uts, y:[[mlon],[mlat]]}
    endfor
    
    ; **** load thm_ast.
    device, decomposed = 0
    loadct, 1
    for i = 0, nprob-1 do begin
        get_data, pre0+'de_fac', data = tmp
        uts = tmp.x
        idx = where(uts ge tr[0] and uts le tr[1], nrec)
        uts = uts[idx]
        nrec = 1
        
        get_data, pre0+'fpt_lonlat', data = tmp
        for j = 0, nrec-1 do begin
            show_time = time_string(uts[j])
            thm_asi_create_mosaic, show_time, /thumb, /no_color, $
                track1 = tmp.y
        endfor
    endfor
        

end

tr = time_double(['2013-04-14/07:00','2013-04-14/10:00'])
scls = [1,10,100]
models = ['t01']
rbsp_efw_plot_thm_ast, tr, models, scls

end