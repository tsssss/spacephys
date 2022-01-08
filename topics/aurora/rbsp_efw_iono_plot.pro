
pro rbsp_efw_iono_plot, tr, model, scl, sites = sites, e_dot0 = e_dot0, $
    probe = probe, type = type

    ; prepare.
    if ~keyword_set(probe) then probe = ['a'] & nprob = n_elements(probe)
    if n_elements(model) eq 0 then model = 't89'
    re = 6374D & re1 = 1D/re
    if n_elements(type) eq 0 then type = 'ast'
    
    ; initiate rbsp efw, emfisis, cdf_leap_second, istp.
    ; the trailing '/' is ANNOYING.
    rbsp_efw_init
    cdf_leap_second_init
    rbsp_emfisis_init
    istp_init
    timespan, tr[0], tr[1]-tr[0], /second
    rbsp_load_spice_kernels
    thm_init
    
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
    
    ; *** decompose e field, ex is along b, ey is east-west, ez is up.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        get_data, pre0+'e_gse', data = de
        get_data, pre0+'db_gse', data = db
        get_data, pre0+'bmod_gse', data = tmp
        bmod = sinterpol(tmp.y, tmp.x, de.x)
        db = sinterpol(db.y, db.x, de.x)
        bhat = sunitvec(bmod)
        
        p = atan(bhat[*,1],bhat[*,0])
        cosp = cos(p) & sint = bhat[*,2]
        cost = bhat[*,0]/cosp & sinp = bhat[*,1]/cost
        
        tmp = de.y
        ex = cost*cosp*tmp[*,0] + cost*sinp*tmp[*,1] - sint*tmp[*,2]
        ey = -sinp*tmp[*,0] + cosp*tmp[*,1]
        ez = sint*sinp*tmp[*,0] + sint*cosp*tmp[*,2] + cost*tmp[*,2]
        store_data, pre0+'de_fac', data = {x:de.x, y:[[ex],[ey],[ez]]}
        ; de_dot0.
        if keyword_set(e_dot0) then begin
            de.y[*,0] = (de.y[*,1]*bmod[*,1]+de.y[*,2]*bmod[*,2])*$
                (-1/bmod[*,0])
            store_data, pre0+'de_gse_dot0', data = {x:de.x, y:de}
            tmp = de.y
            ex = cost*cosp*tmp[*,0] + cost*sinp*tmp[*,1] - sint*tmp[*,2]
            ey = -sinp*tmp[*,0] + cosp*tmp[*,1]
            ez = sint*sinp*tmp[*,0] + sint*cosp*tmp[*,2] + cost*tmp[*,2]
            store_data, pre0+'de_fac_dot0', data = {x:de.x, y:[[ex],[ey],[ez]]}
        endif
        tmp = db
        bx = cost*cosp*tmp[*,0] + cost*sinp*tmp[*,1] - sint*tmp[*,2]
        by = -sinp*tmp[*,0] + cosp*tmp[*,1]
        bz = sint*sinp*tmp[*,0] + sint*cosp*tmp[*,2] + cost*tmp[*,2]
        store_data, pre0+'db_fac', data = {x:de.x, y:[[bx],[by],[bz]]}
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
        uts = tmp.x & ets = 1000D*uts+62167219200000D
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
                epoch = ets[j], /refine, /ionosphere, $
                t89 = t89, t96 = t96, t01 = t01
            ; conver from gsm to mag.
            geopack_conv_coord, xf, yf, zf, /from_gsm, $
                x1, y1, z1, /to_mag
            pos1[j,*] = [x1,y1,z1]
        endfor
        rtod = 180D/!dpi
        mlat = asin(pos1[*,2]*(1/r0))*rtod
        mlon = atan(pos1[*,1],pos1[*,0])*rtod
        store_data, pre0+'fpt_mag', data = {x:uts, y:pos1}
        store_data, pre0+'fpt_lonlat', data = {x:uts, y:[[mlon],[mlat]]}
    endfor
    
    ; options.
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        pre1 = 'rbsp'+probe[i]+'!C'
        vars = pre0+['pos_gsm','e_mgse','e_gse','b_gse','bmod_gse','db_gse',$
            'de_fac','db_fac','fpt_mag','de_gse_dot0','de_fac_dot0']
        options, vars, 'colors', [6,4,2]
        options, vars, 'labels', ['x','y','z']
        
        options, pre0+'pos_gsm', 'ytitle', pre1+'pos gsm (Re)'
        options, pre0+'e_mgse', 'ytitle', pre1+'dE mgse (mV/m)'
        options, pre0+'e_gse', 'ytitle', pre1+'dE gse (mV/m)'
        options, pre0+'b_gse', 'ytitle', pre1+'B gse (nT)'
        options, pre0+'bmod_gse', 'ytitle', pre1+'B model gse (nT)'
        options, pre0+'db_gse', 'ytitle', pre1+'dB gse (nT)'
        options, pre0+'de_fac', 'ytitle', pre1+'dE fac (mV/m)'
        options, pre0+'db_fac', 'ytitle', pre1+'dB fac (nT)'
        options, pre0+'de_gse_dot0', 'ytitle', pre1+'dE gse dot0 (mV/m)'
        options, pre0+'de_fac_dot0', 'ytitle', pre1+'dE fac dot0 (mV/m)'
    endfor
    
           
    ; **** plot footprint and electric on aurora image.
    device, decomposed = 0
    loadct, 1
    picsz = 600
    minlat = 50
    scle = 0.05          ; 1mV/m ~ 10 deg.
    sclb = scle*0.25
    window, 0, xsize = picsz, ysize = picsz
    get_data, 'rbsp'+probe[0]+'_de_fac', data = tmp
    uts = tmp.x & idx = where(uts ge tr[0] and uts le tr[1], nrec)
    uts = uts[idx] & ets = stoepoch(uts, 'unix')
    
    ; trim footpoint, de and db data.
    mlon = dblarr(nrec,nprob) & mlat = dblarr(nrec,nprob)
    defac = dblarr(nrec,3,nprob) & dbfac = dblarr(nrec,3,nprob)
    for i = 0, nprob-1 do begin
        pre0 = 'rbsp'+probe[i]+'_'
        get_data, pre0+'fpt_lonlat', data = tmp
        tmp = sinterpol(tmp.y, tmp.x, uts)
        mlon[*,i] = tmp[*,0] & mlat[*,i] = tmp[*,1]
        ; de and db.
        get_data, pre0+'de_fac', data = tmp
        defac[*,*,i] = sinterpol(tmp.y, tmp.x, uts)
        get_data, pre0+'db_fac', data = tmp
        dbfac[*,*,i] = sinterpol(tmp.y, tmp.x, uts)
    endfor
    ; convert de_fac to position.
    ae = reform(atan(defac[*,2,*],-defac[*,1,*])) & sinae = sin(ae) & cosae = cos(ae)
    ab = reform(atan(dbfac[*,2,*],-dbfac[*,1,*])) & sinab = sin(ab) & cosab = cos(ab)
    de = reform(sqrt(defac[*,2,*]^2+defac[*,1,*]^2)*scle)      ; x deg/mV.
    db = reform(sqrt(dbfac[*,2,*]^2+dbfac[*,1,*]^2)*sclb)      ; x deg/mV.
    
    for j = 0, nrec-1 do begin
        ; aurora image.
        wicimg = read_thm_asi(ets[j], sites, type = type, $
            locroot = spreproot('themis'), midn = midn)
        tv, congrid(wicimg, picsz, picsz)
        map_set, name = 'AzimuthalEquidistant', 90, 0, 0, $
            /noborder, position = [0,0,1,1], /noerase, /isotropic, $
            limit = [minlat-0.01,-180,90,180], color = 255, $
            latdel = 10, londel = 45, glinestyle = 1
        xyouts, [0,88,180,272], minlat+[1,1,2,1], ['00','06','12','18'], $
            alignment = 0.5, /data, color = 255
        xyouts, 45+intarr(5), [5,6,7,8]*10, ['50','60','70','80'], $
            /data, color = 255
        xyouts, 10, 10, sfmepoch(ets[j]), /device, color = 255
        ; reset map coord to the real one.
        map_set, name = 'AzimuthalEquidistant', 90, 0, midn, $
            /noborder, position = [0,0,1,1], /noerase, /isotropic, $
            limit = [minlat-0.01,-180-midn,90,180-midn]
        
        ; plot footpoint, de, db.
        for i = 0, nprob-1 do begin
            device, decomposed = 0
            plots, mlon[*,i], mlat[*,i], color = 255, linestyle = 1
            xyouts, mlon[0,i], minlat+7, probe[i], color = 255, alignment = 0.5
            idx = j
            xy = map_proj_forward(mlon[idx,i], mlat[idx,i])
            t = atan(-xy[0,*],-xy[1,*]) & cost = cos(t) & sint = sin(t)
            x1 = de[idx,i]*(sinae[idx,i]*sint+cosae[idx,i]*cost)+xy[0,*]
            y1 = de[idx,i]*(sinae[idx,i]*cost-cosae[idx,i]*sint)+xy[1,*]
            x2 = db[idx,i]*(sinab[idx,i]*sint+cosab[idx,i]*cost)+xy[0,*]
            y2 = db[idx,i]*(sinab[idx,i]*cost-cosab[idx,i]*sint)+xy[1,*]
            xy = map_proj_inverse(x1, y1)
            x1 = rebin(xy[0,*], n_elements(xy)) & x1[0:*:2] = mlon[idx,i]
            y1 = rebin(xy[1,*], n_elements(xy)) & y1[0:*:2] = mlat[idx,i]
            xy = map_proj_inverse(x2, y2)
            x2 = rebin(xy[0,*], n_elements(xy)) & x2[0:*:2] = mlon[idx,i]
            y2 = rebin(xy[1,*], n_elements(xy)) & y2[0:*:2] = mlat[idx,i]
            device, decomposed = 1
            plots, x2, y2, color = '0000ff'x    ; db fac.
            plots, x1, y1, color = 'ffff00'x    ; de fac.
        endfor
        xyouts, 10, 25, 'de', color = 'ffff00'x, /device
        xyouts, 40, 25, 'db', color = '0000ff'x, /device
        rt = shomedir()+'/'+sfmepoch(ets[0], 'YYYY_MMDD')+'_'+type
        if file_test(rt) eq 0 then file_mkdir, rt
        fn = rt+'/thm_'+type+'_'+sfmepoch(ets[j], 'YYYY_MMDD_hhmmss')+'.png'
        write_png, fn, tvrd(0,0,picsz,picsz*0.5+1,/true)
        device, decomposed = 0
    endfor
    
    vars = pre0+['pos_gsm','bmod_gse','de_fac','de_fac_dot0','db_fac']
    tplot, vars, title = tl

    return
    
end

tr = time_double(['2013-04-14/06:50','2013-04-14/10:00'])
;tr = time_double(['2013-04-14/09:50','2013-04-14/10:00'])
scls = [1,10,100]
models = ['t89']
type = 'ast'
rbsp_efw_iono_plot, tr, models, scls, /e_dot0, probe = ['a','b'], type = type

; rbsp_efw_iono_plot, time_double(['2013-04-14/02:00','2013-04-14/11:00']), 't89', [1,10,100], /e_dot0, probe=['a','b'], type = 'asf'
; rbsp_efw_iono_plot, time_double(['2013-05-01/04:00','2013-05-01/10:00']), 't89', [1,10,100], /e_dot0, probe=['b'], type = 'asf'
; rbsp_efw_iono_plot, time_double(['2013-03-09/00:00','2013-03-09/12:00']), 't89', [1,10,100], /e_dot0, probe=['a','b'], type = 'asf'

vars = ['rbspa_'+['bmod_gse','de_fac','db_fac'],$
    'rbspb_'+['bmod_gse','de_fac','db_fac']]

;thm_init, local_data_dir = spreproot('themis')
;thm_asi_create_mosaic, '2013-04-14/07:00:00', /thumb, /no_color
;ttrace2iono, 'rbspa_pos_gsm', newname = 'rbspa_ifoot_geo', external_model = 't89', par = 2D, in_coord = 'gsm', out_coord = 'geo'
;get_data, 'rbspa_ifoot_geo', data = a
;lat2 = !radeg * atan(a.y[*,2],sqrt(a.y[*,0]^2+a.y[*,1]^2))
;lon2 = !radeg * atan(a.y[*,1],a.y[*,0])
;plots, lon2, lat2, color = 255
end