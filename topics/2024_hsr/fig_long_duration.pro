function fig_long_duration, time_range, probes=probes


;---Load AE and Dst.
    omni_read_index, time_range
    dst_var = 'dst'
    ystep = 50
    yrange = minmax([0,minmax(get_var_data(dst_var))])
    yrange += [-1,1]*abs(total(yrange*[-1,1]))*0.05
    set_ytick, dst_var, ystep=ystep, yrange=yrange
    options, dst_var, 'constant', [-50,0]
    ae_var = 'ae'
    yrange = [0,max(get_var_data(ae_var))*1.05]
    ystep = 600
    set_ytick, ae_var, ystep=ystep, yrange=yrange

    
;---Load RBSP data.
    foreach probe, probes do begin
        mission_probe = 'rbsp'+probe
        prefix = 'rbsp'+probe+'_'
        r_gsm_var = rbsp_read_orbit(time_range, probe=probe)
        b_gsm_var = rbsp_read_bfield(time_range, probe=probe)
        external_models = ['t89','t96','t01','t04s']
        external_models = ['t89']
        nmodel = n_elements(external_models)
        bmod_vars = list()
        foreach external_model, external_models do begin
            bmod_vars.add, lets_read_geopack_bfield(orbit_var=r_gsm_var,internal_model='igrf',external_model=external_model)
        endforeach
        bmod_vars = bmod_vars.toarray()
        mlat_vars = lets_read_mlat_vars(orbit_var=r_gsm_var)
        
        
        bmod_var = bmod_vars[0]
        db_gsm_vars = lets_decompose_bfield(b_var=b_gsm_var, bmod_var=bmod_var)
        
        
    ;---B tilt.
        perigee_dis = 2.
        b_tilt_vars = list()
        foreach var, [b_gsm_var,bmod_vars] do begin
            var_out = streplace(var,'gsm','tilt')
            b_tilt_vars.add, lets_calc_vec_elev(var, coord='sm', var_info=var_out)
        endforeach
        b_tilt_vars = b_tilt_vars.toarray()
        
        b_tilt_var = b_tilt_vars[0]
        bmod_tilt_vars = b_tilt_vars[1:*]
        db_tilt_vars = bmod_tilt_vars
        foreach var_in, db_tilt_vars, vid do begin
            var_out = streplace(var_in, 'b_tilt','db_tilt')
            db_tilt_vars[vid] = lets_subtract_vars(b_tilt_var,var_in,save_to=var_out)
            yrange = [min(get_var_data(var_out))*1.05,20]
            set_ytick, var_out, ystep=25, yrange=yrange
            options, var_out, 'constant', 0
        endforeach
        
        labels = ['Obs',strupcase(external_models)]
        colors = sgcolor(['red','blue','green','orange','purple'])
        b_tilt_combo_var = prefix+'b_tilt_combo'
        b_tilt_combo_var = stplot_merge(b_tilt_vars, output=b_tilt_combo_var, labels=labels[0:nmodel], colors=colors[0:nmodel])
        set_ytick, b_tilt_combo_var, ystep=30, yrange=[0,90]
        options, b_tilt_combo_var, ytitle='(nT)'
        
        db_tilt_combo_var = prefix+'db_tilt_combo'
        db_tilt_combo_var = stplot_merge(db_tilt_vars, output=db_tilt_combo_var, labels=labels[1:nmodel], colors=colors[1:nmodel])
        set_ytick, db_tilt_combo_var, ystep=30, yrange=[0,90]
        options, db_tilt_combo_var, ytitle='(nT)'


    ;---B mag.
        b_mag_vars = list()
        foreach var, [b_gsm_var,bmod_vars] do begin
            var_out = streplace(var,'gsm','mag')
            b_mag_vars.add, lets_calc_vec_mag(var, var_info=var_out)
        endforeach
        b_mag_vars = b_mag_vars.toarray()
    
        b_mag_var = b_mag_vars[0]
        bmod_mag_vars = b_mag_vars[1:*]
        db_mag_vars = bmod_mag_vars
        foreach var_in, db_mag_vars, vid do begin
            var_out = streplace(var_in, 'b_mag','db_mag')
            db_mag_vars[vid] = lets_subtract_vars(b_mag_var,var_in,save_to=var_out)
        endforeach
    
        labels = ['Obs',strupcase(external_models)]
        colors = sgcolor(['red','blue','green','orange','purple'])
        ;b_mag_combo_var = prefix+'b_mag_combo'
        ;b_mag_combo_var = stplot_merge(b_mag_vars, output=b_mag_combo_var, labels=labels, colors=colors)
    
        db_mag_combo_var = prefix+'db_mag_combo'
        db_mag_combo_var = stplot_merge(db_mag_vars, output=db_mag_combo_var, labels=labels[1:nmodel], colors=colors[1:nmodel])

    
    ;---Thermal plasma
        foreach species, ['e','p','o'] do begin
            var = rbsp_read_en_spec(time_range, probe=probe, species=species)
        endforeach
        var = rbsp_read_en_spec_combo(time_range, probe=probe, species='o')
        

    endforeach
    
    labels = strupcase('rbsp-'+probes)
    colors = sgcolor(['red','blue'])
    
    vars = 'rbsp'+probes+'_db_tilt_igrf_t89'
    tilt_var = 'db_tilt_igrf_t89'
    tilt_var = stplot_merge(vars, output=tilt_var, colors=colors, labels=labels)
    yrange = [min(get_var_data(tilt_var))*1.05,10]
    set_ytick, tilt_var, ystep=20, yrange=yrange, yminor=4
    options, tilt_var, constant=0, ytitle='(nT)'
    
    
    vars = 'rbsp'+probes+'_dis'
    var = 'dis'
    var = stplot_merge(vars, output=var, colors=colors, labels=labels)
    yrange = [1,6.5]
    set_ytick, var, ystep=2, yrange=yrange, yminor=2
    options, var, ytitle='(Re)'
    
    stop
end


time_range = time_double(['2015-03-17','2015-03-18/12:00'])
time_range = time_double(['2013-06-06/00:00','2013-06-09'])
;time_range = time_double(['2013-04-30/24:00','2013-05-02/12:00'])
probes = ['a','b']
print, fig_long_duration(time_range, probes=probes)
end