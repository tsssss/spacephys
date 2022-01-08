;+
; load rbsp pos, e, b, and calc de, db, bmod, fpoint.
; calc poynting flux, 1st time interactively, set tscinfo, tftrs, dezr, dbzr once decided.
; calc mosaic thg.
; save data to rbsp_efw_fld_yyyy_mmdd_hh.tplot.
;
; x:along b, y:east-west, z:perp in north sense.
;
; utr is time range in ut.
; probes, set 'a' and/or 'b'.
; sites, set the sites for thg, note do not provide exclude b/c encourage use only several sites.
; model is for footpoint mapping, 't89' by default.
; e_dot0, set to calc spin-axis e using e.b = 0.
; mosaic, set to calc mosaic for given sites, write to a cdf file in home directory.
; extra, transmit keywords for stplot_calc_poynting_flux_mat, including tscinfo, tftrs, dezr, dbzr.
;-
pro rbsp_thg_prep_data, utr, probes = probes, sites = sites, ofn = ofn, $
    model = model, e_dot0 = e_dot0, mosaic = mosaic, _extra = extra

    ; **** check inputs.
    if n_elements(utr) ne 2 then message, 'no time range ...'
    
    if ~keyword_set(probes) then probes = ['a','b']
    nprobe = n_elements(probes)

    nsite = n_elements(sites)
;    if nsite eq 0 then message, 'no thg site ...'
    
    if n_elements(model) eq 0 then model = 't89'
    
    ; **** settings.
    re = 6374d & re1 = 1d/re
    
    ; **** initiate rbsp efw, emfisis, cdf_leap_second, istp.
    rbsp_efw_init
    cdf_leap_second_init
    rbsp_emfisis_init
    istp_init
    timespan, utr[0], utr[1]-utr[0], /second
    rbsp_load_spice_kernels
    
    ; **** load rbsp position in gsm.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        rbsp_load_spice_state, probe = probes[i], $
            coord = 'gsm', /no_spice_load
        get_data, pre0+'state_pos_gsm', data = tmp
        tmp.y *= re1        ; from km to re.
        store_data, pre0+'pos_gsm', data = tmp
        ; delete vars.
        store_data, pre0+'state_*', /delete
    endfor
    
    ; **** load rbsp e and b field.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        ; ** load e in mgse.
        rbsp_load_efw_waveform_l2, probe = probes[i]
        get_data, pre0+'efw_e-spinfit-mgse_efield_spinfit_mgse', data = tmp
        tmp.y[*,0] = 0      ; throw e_spin.
        store_data, pre0+'e_mgse', data = tmp
        ; convert from mgse to gse.
        rbsp_mgse2gse, pre0+'e_mgse', newname = pre0+'e_gse', $
            probe = probes[i], /no_spice_load
        ; delete vars.
        store_data, pre0+'efw_*', /delete
        
        ; ** load b in gse.
        rbsp_load_emfisis, probe = probes[i], coord = 'gse'
        get_data, pre0+'emfisis_l3_4sec_gse_Mag', data = tmp
        store_data, pre0+'b_gse', data = tmp
        ; delete vars.
        store_data, pre0+'emfisis_*', /delete
    endfor
    
    ; **** get model b field.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        get_data, pre0+'b_gse', data = tmp
        bmod = tmp.y
        len = 1200D/sdatarate(tmp.x)    ; 20 min.
        for j = 0, 2 do bmod[*,j] = smooth(bmod[*,j], len)
        store_data, pre0+'bmod_gse', data = {x:tmp.x, y:bmod}
        store_data, pre0+'db_gse', data = {x:tmp.x, y:tmp.y-bmod}
    endfor
    
    ; *** decompose e field, ex is along b, ey is east-west, ez is up.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        get_data, pre0+'e_gse', data = de
        get_data, pre0+'db_gse', data = db
        get_data, pre0+'bmod_gse', data = tmp
        bmod = sinterpol(tmp.y, tmp.x, de.x)
        db = sinterpol(db.y, db.x, de.x)
        bhat = sunitvec(bmod)
        
        p = atan(bhat[*,1],bhat[*,0])
        cosp = cos(p) & sint = bhat[*,2]
        cost = bhat[*,0]/cosp & sinp = bhat[*,1]/cost
        
        vec = de.y
        x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
        y =       -sinp*vec[*,0] + cosp*vec[*,1]
        z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
        store_data, pre0+'de_fac', data = {x:de.x, y:[[x],[y],[z]]}
        ; de_dot0.
        if keyword_set(e_dot0) then begin
            de.y[*,0] = (de.y[*,1]*bmod[*,1]+de.y[*,2]*bmod[*,2])*$
                (-1/bmod[*,0])
            store_data, pre0+'de_gse_dot0', data = {x:de.x, y:de}
            vec = de.y
            x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
            y =       -sinp*vec[*,0] + cosp*vec[*,1]
            z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
            store_data, pre0+'de_fac_dot0', data = {x:de.x, y:[[x],[y],[z]]}
        endif
        vec = db
        x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
        y =       -sinp*vec[*,0] + cosp*vec[*,1]
        z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
        store_data, pre0+'db_fac', data = {x:de.x, y:[[x],[y],[z]]}
    endfor
    
    ; **** map footprint.
    ; prepare mapping.
    dir = -1        ; always north hem, b/c conjugate to thm_asi.
    r0 = 1+110*re1  ; 110 km altitude.
    sgeopack_par, utr, model, /delete  ; get tplot var <model>_par.
    t89 = 0 & t96 = 0 & t01 = 0
    case model of
        't89': t89 = 1
        't96': t96 = 1
        't01': t01 = 1
    endcase
    ; map pos to fpt, convert to mag.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
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
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        pre1 = 'rbsp'+probes[i]+'!C'
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
        
    ; **** calc poynting flux.
;    for i = 0, nprobe-1 do begin
;        pre0 = 'rbsp'+probes[i]+'_'
;        stplot_calc_poynting_flux_mat, pre0+'de_fac', pre0+'db_fac', pre0+'pf_fac', trange = utr, _extra = extra
;    endfor


    if nsite eq 0 then goto, save_exit

    ; **** calc count rate along s/c footpoint.
    etr = stoepoch(utr,'unix')
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        fpvars = pre0+'fpt_lonlat'
        del = 1
        thm_asi_countrate, etr, sites, fpvars, del = del, sc = 'rbsp'+probes[i]
    endfor
    
    ; **** generate thg mosaic.
    if keyword_set(mosaic) then fn = read_thg_asi_mosaic(utr, sites)
    
    save_exit:
    vars = []
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        tvar = ['pos_gsm','e_gse','b_gse','bmod_gse','de_fac','db_fac','fpt_lonlat*','pf_fac_mat']
        vars = tnames([vars, pre0+tvar])
    endfor
    
    if n_elements(ofn) eq 0 then $
        ofn = shomedir()+'/rbsp_efw_fld_'+time_string(utr[0],tformat='YYYY_MMDD_hh')
    
    tplot_save, vars, filename = ofn

end

;; 2013_0501.
utr = time_double(['2013-05-01/04:00','2013-05-01/10:00'])
probes = ['b']
sites = ['atha','tpas']
mosaic = 0
tscinfo = [10,5000,45]
tftrs = [21,65,262,1855]
deyr = [-1,1]
dbyr = [-1,1]

; 2015_0218.
utr = time_double(['2015-02-18/00:00','2015-02-18/02:00'])
probes = ['a']
sites = ['gbay']
mosaic = 0
tscinfo = [10,5000,45]
tftrs = [21,65,262,1855]
deyr = [-1,1]
dbyr = [-1,1]

;; 2013_0414.
;utr = time_double(['2013-04-14/04:00','2013-04-14/10:00'])
;probes = ['a','b']
;sites = ['gill','tpas']
;mosaic = 0
;tscinfo = [10,5000,45]
;deyr = [-1,1]*0.1
;dbyr = [-1,1]*1
;tftrs = [11,76,338,1616]

;; 2013_0501.
;tr = time_double(['2013-05-01/10:00','2013-05-01/12:00'])
;probes = ['b']
;tscinfo = [10,2000,30]
;deyr = [-1,1]*0.1
;dbyr = [-1,1]*1
;tftrs = [10,76.4,316.5,2000]

rbsp_thg_prep_data, utr, probes = probes, sites = sites, mosaic = mosaic, $
    scaleinfo = tscinfo, filter = tftrs, derange = deyr, dbrange = dbyr
end
