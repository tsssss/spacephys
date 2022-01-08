;+
; Load preprocessed Poynting flux data.
; Have compared pflux vs pflux_dot0 for the first 3 years. They are comparable. Thus, will use pflux because it is available for all times.
;
; coord=. Default is gsm.
;-

pro rbsp_read_pflux, time, probe=probe, errmsg=errmsg, coord=coord

    errmsg = ''
    prefix = 'rbsp'+probe+'_'
    if n_elements(coord) eq 0 then coord = 'gsm'
    unit = 'mW/m!U2!N'

    ; read 'rbspx_pf_fac_norm'.
    pflux_grant_read_preprocessed_pflux, time, probe=probe
    pf_fac_norm_var = prefix+'pf_fac_norm'
    get_data, pf_fac_norm_var, common_times, pf_fac, limits=lim
    
    ; Convert earthward to parallel.
    ; because the latter follows the right-hand rule.
    rbsp_read_orbit, time, probe=probe
    r_gsm_var = prefix+'r_gsm'
    interp_time, r_gsm_var, common_times
    
    r_sm = cotran(get_var_data(r_gsm_var), common_times, 'gsm2sm')
    index = where(r_sm[*,2] lt 0, count)
    if count ne 0 then pf_fac[index,0] *= -1
    pf_fac_var = prefix+'pf_fac'
    lim.labels[0] = 'S!D||!N'
    store_data, pf_fac_var, common_times, pf_fac, limits=lim
    add_setting, pf_fac_var, {$
        coord: 'FAC'}
    
    
    ; Convert FAC to GSM.
    b0_gsm_var = prefix+'b0_gsm'
    pflux_grant_read_preprocessed_ebfield, time, probe=probe, id='b0'
    
    define_fac, b0_gsm_var, r_gsm_var
    pf_gsm_var = prefix+'pf_gsm_norm'
    q_var = prefix+'q_gsm2fac'
    from_fac, pf_fac_var, to=pf_gsm_var, q_var=q_var
    
end

time = time_double(['2014-07-31','2014-08-11'])
time = time_double(['2016-05-20','2016-06-01'])
time = time_double(['2018-02-20','2018-03-01'])
probe = 'b'
rbsp_read_pflux, time, probe=probe

prefix = 'rbsp'+probe+'_'
;to_fac, prefix+'pf_gsm_norm', to=prefix+'pf_fac_norm2', q_var=prefix+'q_gsm2fac'
end