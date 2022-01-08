;+
; To check the moment calculated by Cristian Ferradas.
; Check UVW and GSM.
;-

;---Settings.
    time_range = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
    rgb = sgcolor(['red','green','blue'])
    abc = ['a','b','c']
    xyz = ['x','y','z']
    uvw = ['u','v','w']
    fac = ['north','west','para']
    ion_energy_range = [30.,60000]
    probe = 'a'
    test = 0
    tplot_options, 'labflag', -1


;---Load data.
    data_file = join_path([srootdir(),'RBSPa_20130607_mom_avg180_30eV_E_60000eV.tplot'])
    store_data, '*', /delete
    tplot_restore, filename=data_file
    vars = 'V_'+['xyz','uvw','gsm','fac']+'_ion'
    options, vars, 'colors', rgb
    options, 'V_xyz_ion', 'labels', 'ABC V'+abc
    options, 'V_gsm_ion', 'labels', 'GSM V'+xyz
    options, 'V_uvw_ion', 'labels', 'UVW V'+uvw
    options, 'V_fac_ion', 'labels', 'FAC V'+fac
    
    get_data, 'V_gsm_ion', times, vgsm
    store_data, 'V_mag_ion', times, snorm(vgsm), limits={ytitle:'(km/s)', labels:'|V| Cristian'}


    data_file = join_path([srootdir(),'rbspa_20130607_mom_sheng.tplot'])
    if file_test(data_file) eq 0 then begin
        cdf_file = join_path([srootdir(),'rbspa_20130607_mom_sheng.cdf'])
        if file_test(cdf_file) eq 0 then rbsp_gen_hope_moments, probe=probe, time_range, p_erng=ion_energy_range, ofn=cdf_file
        if keyword_set(test) then rbsp_gen_hope_moments, probe=probe, time_range, ion_energy_range=ion_energy_range, ofn=cdf_file

        moms = sread_rbsp_hope_moments(filename=cdf_file)
        times = moms.ut_ion
        vel = moms.o_vbulk
        store_data, 'rbspa_v_gsm_o', times, vel, limits={ytitle:'(km/s)', colors:rgb, labels:'GSM V'+xyz}
        store_data, 'rbspa_v_mag_o', times, snorm(vel), limits={ytitle:'(km/s)', labels:'|V| Sheng'}
        
        eflux = moms.o_energy_flux
        store_data, 'rbspa_eflux_gsm_o', times, eflux, limits={ytitle:'(mW/m!U2!N)', colors:rgb, labels:'GSM T'+xyz}
        
        eflux = moms.o_number_flux
        store_data, 'rbspa_nflux_gsm_o', times, eflux, limits={ytitle:'(#/cm!U2!N-s)', colors:rgb, labels:'GSM F'+xyz}

        stop
    endif

end
