;+
; generate survey plot for 30min. show E/B fields, H/O spectrograms, density, temperature.
;-
pro de_at_gradb_survey_plot3, utr, tprobe, id = id
    
    rootdir = shomedir()+'/de_at_gradb/30min_survey'
    pre0 = 'rbsp'+tprobe+'_'
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
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
    
    
; **** load hope data.
    vars = ['PITCH_ANGLE',$
            'Epoch_Ele','HOPE_ENERGY_Ele',$
            'Epoch_Ion','HOPE_ENERGY_Ion',$
            'FEDU','FPDU','FODU']
    hopel2 = sread_rbsp_hope_l3(utr0, probes = tprobe, vars = vars)
    if size(hopel2,/type) ne 8 then return
    uts = sfmepoch(hopel2.epoch_ion,'unix')
    store_data, pre0+'h_en', uts, total(hopel2.fpdu,3,/nan), hopel2.hope_energy_ion
    store_data, pre0+'h_pa', uts, total(hopel2.fpdu,2,/nan), hopel2.pitch_angle
    store_data, pre0+'o_en', uts, total(hopel2.fodu,3,/nan), hopel2.hope_energy_ion
    store_data, pre0+'o_pa', uts, total(hopel2.fodu,2,/nan), hopel2.pitch_angle
    uts = sfmepoch(hopel2.epoch_ele,'unix')
    store_data, pre0+'e_en', uts, total(hopel2.fedu,3,/nan), hopel2.hope_energy_ele
    store_data, pre0+'e_pa', uts, total(hopel2.fedu,2,/nan), hopel2.pitch_angle

    vars = pre0+['h_en','h_pa','o_en','o_pa','e_en','e_pa']
    options, vars, 'spec', 1
    options, vars, 'no_interp', 1
    options, vars, 'zlog', 1
    options, vars, 'ztitle', '(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'

    vars = pre0+['h_en','o_en','e_en']
    options, vars, 'ylog', 1
    options, vars, 'yrange', [10,4e4]
    
    vars = pre0+['e_en','e_pa']
    options, vars, 'zrange', [1e4,1e10]
    vars = pre0+['h_en','h_pa']
    options, vars, 'zrange', [1e4,1e7]
    vars = pre0+['o_en','o_pa']
    options, vars, 'zrange', [1e4,1e6]

    vars = pre0+['h_pa','o_pa','e_pa']
    options, vars, 'yrange', [0,180]
    
    tvar = pre0+'e_en'
    options, tvar, 'ytitle', 'e- energy!C(eV)'
    tvar = pre0+'h_en'
    options, tvar, 'ytitle', 'H+ energy!C(eV)'
    tvar = pre0+'o_en'
    options, tvar, 'ytitle', 'O+ energy!C(eV)'

    tvar = pre0+'e_pa'
    options, tvar, 'ytitle', 'e- pitch!C(deg)'
    tvar = pre0+'h_pa'
    options, tvar, 'ytitle', 'H+ pitch!C(deg)'
    tvar = pre0+'o_pa'
    options, tvar, 'ytitle', 'O+ pitch!C(deg)'
    
    hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom')
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



; **** load E field. l3, mgse. And pos_gse, mlt/mlat/lshell.
; rbspx_de_[mgse,svy], rbspx_[mlt,lshell,ilat,mlat].
    efwl3 = sread_rbsp_efw_l3(utr, probes = tprobe)
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
    rbsp_load_efw_waveform, probe = tprobe, trange = utr0, datatype = vars, type = 'calibrated'
    tvar = pre0+'efw_esvy'
    get_data, tvar, tmp, dat
    dat[*,2] = 0
    store_data, pre0+'de_svy', tmp, dat, $
        limits = {ytitle:'dE survey!C(mV/m)', colors:rgb, labels:'UVW '+['x','y','z']}
    store_data, pre0+'efw_esvy*', /delete

    
; **** load B field. l3, gse, emfisis, 1sec.
; rbspx_b_gse, rbspx_b.

    emfisis = sread_rbsp_emfisis(utr0, probes = tprobe)
    if size(emfisis,/type) ne 8 then return
    uts = sfmepoch(emfisis.epoch,'unix',/tt2000)
    store_data, pre0+'b_gse', uts, emfisis.mag, $
        limits = {ytitle:'B!C(nT)',colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, pre0+'b', uts, sqrt(total(emfisis.mag^2,2)), $
        limits = {ytitle:'B magnitude!C(nT)'}


; **** map footprint.
    ; prepare mapping.
    r0 = 1+110*re1  ; 110 km altitude.
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
    store_data, pre0+'fpt_mag', data = {x:uts, y:pos1}
    store_data, pre0+'fpt_mlat', data = {x:uts, y:mlat}, $
        limits = {ytitle:'MLat/fpt (deg)'}
    
    
    ; plasma beta. 2mu0*nkT/B^2
    tmp = 2*4*!dpi*1e-7*1e6*1.6e-19*1e18    ; 
    get_data, pre0+'n', t0, dat & dat*= tmp
    get_data, pre0+'t', t0, tmp & dat*= tmp[*,0]
    get_data, pre0+'b', t1, tmp & tmp = interpol(tmp, t1, t0) & dat/= tmp^2
    store_data, pre0+'beta', t0, dat, $
        limits = {ytitle:'Beta'}
        
    
    ; e_gse.
    get_data, pre0+'de_mgse', t0, dat
    dat[*,0] = 0
    store_data, pre0+'de_mgse', t0, dat
    rbsp_mgse2gse, pre0+'de_mgse', newname = pre0+'de_gse'
    get_data, pre0+'de_gse', t0, de_gse
    get_data, pre0+'b_gse', tmp, db_gse
    db_gse = sinterpol(db_gse, tmp, t0)
    
    tvar = pre0+'de_gse'
    options, tvar, 'colors', [6,4,2]
    options, tvar, 'labels', 'GSE '+['x','y','z']
    options, tvar, 'ytitle', 'dE GSE!C(mV/m)'
    
    ; angle between e and b.
    tvar = pre0+'ebangle'
    store_data, tvar, t0, sang(de_gse, db_gse, /degree)
    ylim, tvar, 70,110, 0
    options, tvar, 'ytitle', 'E/B angle!C(deg)'

    vars = pre0+['de_mgse','de_svy','de_gse','b_gse','b','e_en','h_en','o_en','n','t','beta','ebangle']
    labs = pre0+['mlt','lshell','ilat','fpt_mlat']
    titl = 'RBSP-'+strupcase(tprobe)+' survey plot, '+ $
        time_string(utr[0])+' to '+time_string(utr[1],tformat='hh:mm')
    
    
;    ofn = 0
;    sgopen, ofn, xsize = 8.5, ysize = 11, /inch
;    device, decomposed = 0
;    loadct2, 43
;    tplot, vars, var_label = labs, trange = utr, title = titl
;    sgclose
    
    
    
    ofn = rootdir+'/de_at_gradb_'+id+'_rbsp'+tprobe+'.pdf'
    sgopen, ofn, xsize = 8.5, ysize = 11, /inch
    device, decomposed = 0
    loadct2, 43

    titl = 'RBSP-'+strupcase(tprobe)+' survey plot, '+ $
        time_string(utr[0])+' to '+time_string(utr[1])

    tplot, vars, var_label = labs, trange = utr, title = titl
    sgclose
end