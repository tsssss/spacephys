;+
; Load v_gsm and check the direction.
;-

pro check_themis_vel, time_range, probe=probe

    themis_read_ion_vel, time_range, probe=probe
    prefix = 'th'+probe+'_'
    v_gsm_var = prefix+'u_gsm'
    get_data, v_gsm_var, times, v_gsm
    v_sm = cotran(v_gsm, times, 'gsm2sm')
    

    themis_read_orbit, time_range, probe=probe
    r_gsm_var = prefix+'r_gsm'
    r_gsm = get_var_data(r_gsm_var, at=times)
    r_sm = cotran(r_gsm, times, 'gsm2sm')
    
    
    rhat = -r_sm & rhat[*,2] = 0
    rhat = sunitvec(rhat)
    zhat = r_sm & zhat[*,[0,1]] = 0 & zhat[*,2] = 1
    what = sunitvec(vec_cross(zhat, rhat))
    zhat = vec_cross(rhat, what)
    
    vr = vec_dot(v_sm, rhat)
    vz = vec_dot(v_sm, zhat)
    vw = vec_dot(v_sm, what)
    v_var = prefix+'v_fac'
    store_data, v_var, times, [[vr],[vw],[vz]]
    add_setting, v_var, /smart, dictionary($
        'display_type', 'vector', $
        'unit', 'km/s', $
        'short_name', 'V', $
        'coord', 'Cyn', $
        'coord_labels', ['r','west','z'] )
        
        
    themis_read_bfield, time_range, probe=probe
    b_gsm_var = prefix+'b_gsm'
    define_fac, b_gsm_var, r_gsm_var, time_var=v_gsm_var
    to_fac, v_gsm_var
    
end

time_range = time_double(['2014-08-28/09:50','2014-08-28/11:00'])
probes = ['a','d','e']


time_range = time_double(['2008-01-09/11:00','2008-01-09/12:10'])
probes = ['a','c','d','e']


time_range = time_double(['2017-03-28/02:30','2017-03-28/03:45'])
probes = ['a','d','e']


time_range = time_double(['2008-02-29/08:00','2008-02-29/09:00'])
probes = ['a','c','d','e']


time_range = time_double(['2008-03-20/11:30','2008-03-20/12:50'])
probes = ['c','b','d','e']


azim_dp_read_theta, time_range, probe='thd'


foreach probe, probes do begin
    azim_dp_read_theta, time_range, probe='th'+probe
    check_themis_vel, time_range, probe=probe
    stplot_split, 'th'+probe+'_v_fac', newnames='th'+probe+'_v'+['r','w','z']
    stplot_split, 'th'+probe+'_u_fac', newnames='th'+probe+'_u'+['r','w','z']
endforeach

foreach probe, probes do begin
    print, probe, mean(get_var_data('th'+probe+'_uw'),/nan)
endforeach

end