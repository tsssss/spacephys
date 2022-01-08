;+
; Check ion_vel data and direction relative to B.
;-

    time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
    probes = ['a','d','e']
    foreach probe, probes do begin
        var = 'th'+probe+'_u_gsm'
        if ~check_if_update(var, time_range, dtime=600) then continue
        themis_read_ion_vel, time_range, probe=probe
    endforeach

    _2014_0828_10_load_data

    foreach probe, probes do begin
        prefix = 'th'+probe+'_'

        define_fac, prefix+'b_gsm', prefix+'r_gsm', time_var=prefix+'b_gsm'
        to_fac, prefix+'u_gsm', to=prefix+'u_fac'
        
        get_data, prefix+'u_gsm', times, u_gsm
        r_gsm = get_var_data(prefix+'r_gsm', at=times)
        ur = -vec_dot(u_gsm, sunitvec(r_gsm))
        store_data, prefix+'ur', times, ur
    endforeach
    
    vars = []
    foreach probe, probes do begin
        prefix = 'th'+probe+'_'
        vars = [vars, prefix+['b_gsm','ur']]
    endforeach
    tplot, vars

end
