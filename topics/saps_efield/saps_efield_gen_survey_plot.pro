function saps_efield_gen_survey_plot, input_time_range, probe=probe

    ; ion
    ; ele
    ; E field
    tr = time_double(input_time_range)
    prefix = 'rbsp'+probe+'_'

    ion_var = rbsp_read_en_spec(tr, probe=probe, species='p')
    ele_var = rbsp_read_en_spec(tr, probe=probe, species='e')
    o_var = rbsp_read_en_spec(tr, probe=probe, species='o')
    e_var = rbsp_read_efield(tr, probe=probe)
    b_var = rbsp_read_bfield(tr, probe=probe)
    edot0_mgse_var = lets_calc_edotb0(e_var=e_var, b_var=b_var, $
        var_info=prefix+'edot0_rbsp_mgse')
    mlt_var = rbsp_read_mlt(tr, probe=probe)
    r_var = rbsp_read_orbit(tr, probe=probe)
    options, r_var, 'mission_probe', 'rbsp'+probe
    
    q_gsm2fac_var = lets_define_fac(b_var=b_var, r_var=r_var)
    edot0_fac_var = prefix+'edot0_fac'
    foreach var, [edot0_mgse_var] do begin
        in_coord = strlowcase(get_setting(var,'coord'))
        out_var = streplace(var,in_coord,'fac')
        var_fac = lets_cotran([in_coord,'fac'], input=var, output=out_var, q_var=q_gsm2fac_var)
    endforeach
    
    options, e_var, 'yrange', [-1,1]*15
    stop


end


tr = ['2013-06-01','2013-06-02']
tr = ['2013-03-17','2013-03-18']
probes = ['a','b']
foreach probe, probes do begin
    print, saps_efield_gen_survey_plot(tr, probe=probe)
endforeach

end