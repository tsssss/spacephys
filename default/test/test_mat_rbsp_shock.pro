
pro test_mat_rbsp_shock
    compile_opt idl2
    
    ; prepare.
    probes = ['a']
    nprobe = n_elements(probes)
    tr = time_double(['2013-10-08/20:00','2013-10-08/21:00'])
    timespan, tr, tr[1]-tr[0], /second
    rbsp_efw_init
    tplot_options, 'labflag', 1
;    rbsp_load_spice_kernels
    re = 6374D & re1 = 1D/re
    tpad = 1200d        ; 20 min.

    ; **** load rbsp position in gse.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        rbsp_load_spice_state, probe = probes[i], coord = 'gse', /no_spice_load
        get_data, pre0+'state_pos_gse', data = tmp
        tmp.y *= re1        ; from km to re.
        store_data, pre0+'pos_gse', data = tmp
        ; delete vars.
        store_data, pre0+'state_*', /delete
    endfor
    
    ; **** load rbsp e and b field.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        ; ** load e in mgse.
        rbsp_load_efw_waveform_l2, probe = probes[i]
        get_data, pre0+'efw_e-spinfit-mgse_e12_spinfit_mgse', data = tmp
        tmp.y[*,0] = 0      ; throw e_spin.
        store_data, pre0+'de_mgse', data = tmp
        ; convert from mgse to gse.
        rbsp_mgse2gse, pre0+'de_mgse', newname = pre0+'de_gse', $
            probe = probes[i], /no_spice_load
        ; delete vars.
        store_data, pre0+'efw_*', /delete
        
        ; ** load b in gse.
        rbsp_load_emfisis, probe = probes[i], coord = 'gse'
        get_data, pre0+'emfisis_l3_4sec_gse_Mag', data = tmp
        store_data, pre0+'b_gse', data = tmp
        ; delete vars.
        store_data, pre0+'emfisis_*', /delete
        dr = sdatarate(tmp.x)
    endfor
    
    ; **** get model b field.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        get_data, pre0+'b_gse', data = tmp
        bmod = tmp.y
        for j = 0, 2 do bmod[*,j] = smooth(bmod[*,j], tpad/dr)
        store_data, pre0+'bmod_gse', data = {x:tmp.x, y:bmod}
        store_data, pre0+'db_gse', data = {x:tmp.x, y:tmp.y-bmod}
    endfor
    
    ; **** uniform time.
    prea = 'rbspa_'
    timespan, tr[0], tr[1]-tr[0], /second
    vars = prea+['b_gse','bmod_gse','db_gse','de_gse']
    stplot_intpl, vars, t0, tinfo = [tr[0]-tpad,tr[1]+tpad,dr]
    vars = prea+['b_gse','bmod_gse','db_gse','de_gse','pos_gse']
    options, vars, 'colors', [6,4,2]
    options, vars, 'labels', ['x','y','z']
    vars = prea+['bmod_gse','db_gse','b_gse','de_gse','pos_gse']
    tplot, vars
end

stplot_calc_poynting_flux, 'rbspa_de_gse', 'rbspa_db_gse', 'rbspa_pf_gse', $
    scaleinfo = [5,500,20], filter = [10,40,200]
zlim, 'rbspa_db_gse_comp?_mat', -4, 4, 0
get_data, 'rbspa_pf_gse_mat', t0, pf
get_data, 'rbspa_bmod_gse', t0, bmod
pfpara = sdot(pf, sunitvec(bmod))
store_data, 'rbspa_pfpara_gse', t0, pfpara
tplot, 'rbspa_d?_gse_comp?_mat'
stop
tplot, ['rbspa_pf_gse_matf*','rbspa_pfpara_gse']
stop
tplot, ['rbspa_'+['de','db','pf']+'_gse_matf1',$
    'rbspa_'+['de','db','pf']+'_gse_matf2']
end