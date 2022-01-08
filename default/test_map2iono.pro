; 1, difference is small (0.98-1.00) b/w mappings use (1) distance and 
;   (2) ratio of B in situ and B 100 km (dipole model).
; 2, B in situ from goepack dipole model is close to B in situ measured
;   by Polar, but slightly smaller (~95%), b/c diamagnetic effect(?).
; 3, mapping uses distance should be fine, the error is 2% for 4Re.

pro test_map2iono

    re = 6378.137d  ; km.
    r0 = 100/re+1d  ; re.

    ifnfa = sdiskdir('Works')+'/works/polarcap/data/fa_sdt_fld_19980925_08278.tplot'
    trfa = time_double(['1998-09-25/04:27','1998-09-25/04:32'])
    trpo = time_double(['1998-09-25/05:00','1998-09-25/06:00'])
    
    ; load polar eflux and magnetic field.
    ; the model field returned is t96 from sdt.
    polar_sdt_prep_poynting_flux, $
        sdiskdir('Works')+'/works/polarcap/data/po_sdt_fld_1998092415.sdt'
    read_polar_ke, sdiskdir('Works')+'/agu/agu2013/data/polar_hyd_kei.dat', $
        'kei', time_double('1998-09-25/05:00')
    
    ; method 1, use distance.
    get_data, 'kei', uts, kei
    get_data, 'dis', tmp, dis
    dis = sinterpol(dis, tmp, uts)
    dat = kei*(dis*(1d/r0))^3
    store_data, 'kei_map1', uts, dat
    
    ; method 2, use extremely simple dipole model.
    get_data, 'bmod_spc', tmp, b1
    b1 = sinterpol(b1, tmp, uts)
    b1 = sqrt(total(b1^2, 2))
    get_data, 'mlat', tmp, mlat
    mlat = sinterpol(mlat, tmp, uts)
    coef = sqrt(3*sin(mlat*!dpi/180)^2+1d)
    b0 = 31025.2/(r0^3)*coef
    bdipole = 31025.2/(dis^3)*coef
    dat = kei*b0/b1
    store_data, 'kei_map2', uts, dat
    
    ; method 3, use geopack dipole model.
    ets = stoepoch(uts,'unix')
    tmp = sread_polar_pos(t = stoepoch(trpo,'unix'))
    pos = sgse2gsm(tmp.pos_gse, ets)*(1d/re)
    pos = sinterpol(pos, tmp.epoch, ets)
    b2 = pos
    for i = 0, n_elements(ets)-1 do begin
        cdf_epoch, ets[i], yr, mo, dy, hr, mi, sc, /breakdown_epoch
        geopack_recalc, yr, stodoy(yr,mo,dy), hr, mi, sc, tilt = tilt
        geopack_dip, pos[i,0], pos[i,1], pos[i,2], bx, by, bz, tilt = tilt
        b2[i,*] = [bx,by,bz]
    endfor
    b2 = sqrt(total(b2^2, 2))
    dat = kei*b0/b2
    store_data, 'kei_map3', uts, dat
    
    get_data, 'kei_map1', uts, map1
    get_data, 'kei_map2', uts, map2
    store_data, 'ratio2', uts, map2/map1
    store_data, 'bratio1', uts, bdipole/b1  ; simple dipole/t96.
    store_data, 'bratio2', uts, b2/b1       ; geopack dipole/t96.
    
    vars = ['kei_map1','kei_map2','kei_map3', 'ratio2','bratio1','bratio2']
    ylim, vars, -20, 15, 0
    ylim, ['ratio2','bratio1','bratio2'], 0.9, 1.1, 0
    tplot, vars, trange = trpo, var_label = ['ilat','mlt','dis']
   
end
