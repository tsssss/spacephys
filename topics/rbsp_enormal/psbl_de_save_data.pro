
pro psbl_de_save_data, id, trange = utr, rootdir = rootdir

    if n_elements(id) eq 0 then message, 'no id ...'
    if n_elements(probe) eq 0 then probe = strmid(id,strlen(id)-1)
    if n_elements(utr) eq 0 then $
        utr = time_double(strmid(id,0,14),tformat='YYYY_MMDD_hhmm')+[-2,8]*60
    if n_elements(rootdir) eq 0 then rootdir = shomedir()+'/psbl_de'
    
    pre0 = 'rbsp'+probe+'_'
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    minratio = 0.1  ; Bx/B ratio in MGSE>
    rgb = [6,4,2]
    dt0 = 86400d
    utr0 = utr-(utr mod dt0)+[0,dt0]
    timespan, utr0[0], dt0, /second
    
    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero', 1
    tplot_options, 'ymargin', [5,5]
    tplot_options, 'xmargin', [25,15]
    
    
; **** load E field. l3, mgse. And pos_gse, mlt/mlat/lshell.
; rbspx_de_[mgse,svy], rbspx_[mlt,lshell,ilat,mlat].
    
    efwl3 = sread_rbsp_efw_l3(utr0, probes = probe)
    if size(efwl3,/type) ne 8 then return
    uts = sfmepoch(efwl3.epoch,'unix',/epoch16)
    
    store_data, pre0+'de_mgse', uts, efwl3.efield_inertial_frame_mgse, $
        limits = {ytitle:'dE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
    store_data, pre0+'mlt', uts, (efwl3.mlt_lshell_mlat)[*,0], $
        limits = {ytitle:'MLT (hr)'}
    store_data, pre0+'lshell', uts, (efwl3.mlt_lshell_mlat)[*,1], $
        limits = {ytitle:'L (Re)'}
    store_data, pre0+'ilat', uts, acos(1/sqrt((efwl3.mlt_lshell_mlat)[*,1]))*deg, $
        limits = {ytitle:'ILat (deg)'}
    store_data, pre0+'mlat', uts, (efwl3.mlt_lshell_mlat)[*,2], $
        limits = {ytitle:'MLat (deg)'}
    store_data, pre0+'pos_gse', uts, efwl3.pos_gse, $
        limits = {colors:rgb, labels:['x','y','z']}
        
        
    vars = ['esvy']
    rbsp_load_efw_waveform, probe = probe, trange = utr0, datatype = vars, type = 'calibrated'
    tvar = pre0+'efw_esvy'
    get_data, tvar, tmp, dat
    dat[*,2] = 0
    store_data, pre0+'de_svy', tmp, dat, $
        limits = {ytitle:'dE survey!C(mV/m)', colors:rgb, labels:'UVW '+['x','y','z']}
    store_data, pre0+'efw_esvy*', /delete
    
    
; **** load B field. l3, gse, emfisis, 1sec.
; rbspx_b_gse, rbspx_b.
    
    emfisis = sread_rbsp_emfisis(utr0, probes = probe)
    if size(emfisis,/type) ne 8 then return
    uts = sfmepoch(emfisis.epoch,'unix',/tt2000)
    store_data, pre0+'b_gse', uts, emfisis.mag, $
        limits = {ytitle:'B!C(nT)',colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, pre0+'b', uts, sqrt(total(emfisis.mag^2,2)), $
        limits = {ytitle:'B magnitude!C(nT)'}
        
        
; **** map footprint.
; prepare mapping.
    model = 't89'
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
    get_data, pre0+'pos_gse', data = tmp
    uts = tmp.x & ets = 1000D*uts+62167219200000D
    pos0 = tmp.y*re1 & pos1 = pos0      ; in re.
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
        ; pos in gse, which is the mapping coord.
        x0 = pos0[j,0] & y0 = pos0[j,1] & z0 = pos0[j,2]
        ; convert from gse to gsm.
        geopack_conv_coord, x0, y0, z0, /from_gse, $
            x1, y1, z1, /to_gsm
        dir = (z1 gt 0)? -1: 1
        if model ne 't89' then par = reform(pars[j,*]) else par = 2
        geopack_trace, x1, y1, z1, dir, par, xf, yf, zf, $
            epoch = ets[j], /refine, /ionosphere, $
            t89 = t89, t96 = t96, t01 = t01
        ; convert from gse to mag.
        geopack_conv_coord, xf, yf, zf, /from_gsm, $
            x1, y1, z1, /to_mag
        pos1[j,*] = [x1,y1,z1]
    endfor
    mlat = asin(pos1[*,2]*(1/r0))*deg
    mlon = atan(pos1[*,1],pos1[*,0])*deg
    store_data, pre0+'fpt_mlon', uts, mlon, $
        limits = {ytitle:'MLon/fpt (deg)'}
    store_data, pre0+'fpt_mlat', uts, mlat, $
        limits = {ytitle:'MLat/fpt (deg)'}
    store_data, pre0+'fpt_mlt', uts, slon2lt(mlon, stoepoch(uts,'unix'), /mag, /deg)/15, $
        limits = {ytitle:'MLT/fpt (hr)'}
        
; **** calc b_mgse, de_dot0_[mgse,gse].
    get_data, pre0+'de_mgse', t0, de_mgse
    de_mgse[*,0] = 0
    store_data, pre0+'de_mgse', t0, de_mgse
    rbsp_mgse2gse, pre0+'de_mgse', newname = pre0+'de_gse'
    get_data, pre0+'de_gse', t0, de_gse
    get_data, pre0+'b_gse', tmp, db_gse
    db_gse = sinterpol(db_gse, tmp, t0)
    store_data, pre0+'b_gse', t0, db_gse
    rbsp_mgse2gse, pre0+'b_gse', newname = pre0+'b_mgse', /inverse
    get_data, pre0+'b_mgse', tmp, b_mgse
    de_mgse[*,0] = -(b_mgse[*,1]*de_mgse[*,1]+b_mgse[*,2]*de_mgse[*,2])/b_mgse[*,0]
    ratio = abs(b_mgse[*,0])/snorm(b_mgse)
    store_data, pre0+'bx_ratio', t0, ratio
    idx = where(ratio le minratio, cnt)
    if cnt ne 0 then de_mgse[idx,0] = !values.d_nan
    store_data, pre0+'de_dot0_mgse', t0, de_mgse
    rbsp_mgse2gse, pre0+'de_dot0_mgse', newname = pre0+'de_dot0_gse'
    
    tvar = pre0+'bx_ratio'
    options, tvar, 'ytitle', 'Bx/B MGSE'
    ylim, tvar, 0,1, 0
    options, tvar, 'constant', minratio
    
    tvar = pre0+'de_dot0_gse'
    options, tvar, 'colors', rgb
    options, tvar, 'labels', 'dot0 GSE '+['x','y','z']
    options, tvar, 'ytitle', 'dE dot0!C(mV/m)'
    tvar = pre0+'de_dot0_mgse'
    options, tvar, 'colors', rgb
    options, tvar, 'labels', 'dot0 MGSE '+['x','y','z']
    options, tvar, 'ytitle', 'dE dot0!C(mV/m)'
    tvar = pre0+'b_mgse'
    options, tvar, 'colors', rgb
    options, tvar, 'labels', 'MGSE '+['x','y','z']
    options, tvar, 'ytitle', 'B!C(nT)'
    tvar = pre0+'de_gse'
    options, tvar, 'colors', rgb
    options, tvar, 'labels', 'GSE '+['x','y','z']
    options, tvar, 'ytitle', 'dE!C(mV/m)'
    
    ; angle between e and b.
    tvar = pre0+'ebangle'
    store_data, tvar, t0, sang(de_gse, db_gse, /degree)
    ylim, tvar, 30,150, 0
    options, tvar, 'ytitle', 'E/B angle!C(deg)'
    options, tvar, 'labels', 'E GSE/B GSE'
    options, tvar, 'constant', 90
    
    
; **** calc de_fac, de_dot0_fac, db_fac.
; decompose background B and dB.
    get_data, pre0+'b_gse', t0, bgse
    b0gse = bgse
    for i = 0, 2 do b0gse[*,i] = scalcbg(bgse[*,i])
    store_data, pre0+'b0_gse', t0, b0gse
    store_data, pre0+'db_gse', t0, bgse-b0gse
    
    ; FAC dE. ex is along b, ey is east-west, ez is up.
    get_data, pre0+'de_gse', data = de
    get_data, pre0+'de_dot0_gse', data = dedot0
    get_data, pre0+'db_gse', data = db
    get_data, pre0+'b0_gse', data = tmp
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
    store_data, pre0+'de_fac', data = {x:de.x, y:[[x],[y],[z]]}, $
        limits = {colors:rgb, labels:'dE FAC '+['b','p','v'], ytitle:'dE!C(mV/m)'}
    
    vec = dedot0.y
    x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
    y =       -sinp*vec[*,0] + cosp*vec[*,1]
    z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
    store_data, pre0+'de_dot0_fac', data = {x:dedot0.x, y:[[x],[y],[z]]}, $
        limits = {colors:rgb, labels:'dE dot0 FAC '+['b','p','v'], ytitle:'dE!C(mV/m)'}

    vec = db
    x =  cost*(cosp*vec[*,0] + sinp*vec[*,1]) + sint*vec[*,2]
    y =       -sinp*vec[*,0] + cosp*vec[*,1]
    z = -sint*(cosp*vec[*,0] + sinp*vec[*,1]) + cost*vec[*,2]
    store_data, pre0+'db_fac', data = {x:dedot0.x, y:[[x],[y],[z]]}, $
        limits = {colors:rgb, labels:'dB FAC '+['b','p','v'], ytitle:'dB!C(nT)'}


; **** calc poynting flux.
    tfilters = [20,60,250]
    filterids = ['1','2']
    nfilter = 2
    coordlabs = ['b','p','v']
    rgb = [6,4,2]
    dt = 20*60
    sclinfo = [10,dt/2,40]
    
    denames = pre0+['de_fac','de_dot0_fac']
    dbnames = pre0+['db_fac','db_fac']
    pfnames = pre0+['pf_fac','pf_dot0_fac']
    suffixs = '_mat'+filterids
    
    ; calc poynting flux.
    for i = 0, n_elements(denames)-1 do begin
        dename = denames[i]
        dbname = dbnames[i]
        pfname = pfnames[i]
        
        stplot_calc_pflux_mat, dename, dbname, pfname, $
            filter = tfilters, scaleinfo = sclinfo
            
        vars = [dename,dbname,pfname]+'_mat'
        options, tvar, 'constant', 0
        
        ; map poynting flux.
        vars = pfname+[suffixs,'_mat']
        for j = 0, n_elements(vars)-1 do begin
            tvar = vars[j]
            get_data, pre0+'b0_gse', t0, dat
            store_data, pre0+'b0', t0, sqrt(total(dat^2,2))
            smap2iono, tvar, pre0+'fpt_mlat', b = pre0+'b0', $
                newname = tvar+'_map', coef = coef
            tvar = vars[j]+'_map'
            options, tvar, 'constant', 0
            options, tvar, 'colors', rgb
        endfor
        store_data, pre0+'_map_coef', t0, coef
    endfor
    
; **** HOPE moments.
    hopemom = sread_rbsp_hope_l3(utr0, probes = probe, type = 'mom')
    if size(hopemom,/type) ne 8 then return
    uts = sfmepoch(hopemom.epoch_ele,'unix')
    store_data, pre0+'n', uts, hopemom.dens_e_200, $
        limits = {ytitle:'N!Ie!N!C(cm!E-3!N)', ylog:1, constant:1}
    store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
        limits = {ytitle:'T!I!N!C(eV)', ylog:1, colors:[6,0], labels:['Tperp','Tpara'], constant:1000}
    tvar = pre0+'t'
    get_data, tvar, uts, dat
    idx = where(dat eq 1e20, cnt)
    if cnt ne 0 then dat[idx] = !values.d_nan
    store_data, tvar, uts, dat
    
    ; plasma beta. 2mu0*nkT/B^2
    tmp = 2*4*!dpi*1e-7*1e6*1.6e-19*1e18    ;
    get_data, pre0+'n', t0, dat & dat*= tmp
    get_data, pre0+'t', t0, tmp & dat*= tmp[*,0]
    get_data, pre0+'b', t1, tmp & tmp = interpol(tmp, t1, t0) & dat/= tmp^2
    store_data, pre0+'beta', t0, dat, $
        limits = {ytitle:'Beta'}
        
    
    
; **** omni Dst and AE.
    omni = sread_omni(utr0)
    uts = sfmepoch(omni.epoch,'unix')
    store_data, 'dst', uts, omni.symh, limits = {ytitle:'Dst (nT)'}
    store_data, 'ae', uts, omni.ae, limits = {ytitle:'AE (nT)'}
        
    
; save data to disk.
    vars = pre0+['de_mgse','de_gse','de_fac','de_svy', $
        'de_dot0_mgse','de_dot0_gse','de_dot0_fac', $
        'b_mgse','b_gse','db_gse','db_fac','b0_gse', $
        'bx_ratio','ebangle','ilat','mlt','lshell', $
        'fpt_mlat','fpt_mlon','fpt_mlt', $
        'de_fac_mat'+['','1','2'], $
        'de_dot0_fac_mat'+['','1','2'], $
        'db_fac_mat'+['','1','2'], $
        'pf_fac_mat'+['','1','2']+'_map', $
        'pf_dot0_fac_mat'+['','1','2']+'_map', 'map_coef']
    nvar = n_elements(vars)
    for i = 0, nvar-1 do begin
        get_data, vars[i], t0, dat
        idx = where(t0 ge utr[0] and t0 le utr[1])
        store_data, vars[i], t0[idx], dat[idx,*]
    endfor
    ofn = rootdir+'/psbl_de_'+id+'.tplot'
    vars = [vars,pre0+['n','t','beta'],'dst','ae']
    tplot_save, vars, filename = ofn
end

;ilogfn = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
;    'dipolarization/list_large_de_round3.log'
;nheader = 3
;headers = strarr(nheader)
;nline = file_lines(ilogfn)-nheader
;lines = strarr(nline)
;openr, lun, ilogfn, /get_lun
;readf, lun, headers
;readf, lun, lines
;free_lun, lun
;
;for i = 0, nline-1 do begin
;    tline = lines[i]
;    id = strmid(tline,0,16)
;    
;    t1 = time_double(strmid(id,0,9)+strmid(tline,20,5),tformat='YYYY_MMDDhh:mm')
;    t2 = time_double(strmid(id,0,9)+strmid(tline,29,5),tformat='YYYY_MMDDhh:mm')
;    if t2 lt t1 then t2+= 86400d
;
;    psbl_de_save_data, id, trange = [t1,t2]+[-2,2]*60
;endfor
;
;end
