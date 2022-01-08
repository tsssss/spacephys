;+
; Remove background Efield, including E_coro and E_vxb.
;
; b_var=.
; e_var=.
; v_var=.
; r_var=.
; save_to=.
; probe=.
;-

pro rbsp_remove_efield_bg, e_var=e_var, b_var=b_var, $
    v_var=v_var, r_var=r_var, save_to=de_var, probe=probe

    if n_elements(de_var) eq 0 then message, 'No output e_var ...'
    prefix = 'rbsp'+probe+'_'
    get_data, e_var, times, e_vec, limits=limits
    coord = strlowcase(limits.coord)

    ecoro_var = prefix+'ecoro_'+coord
    calc_ecoro, b_var=b_var, r_var=r_var, save_to=ecoro_var, probe=probe
    evxb_var = prefix+'evxb_'+coord
    calc_evxb, b_var=b_var, v_var=v_var, save_to=evxb_var, probe=probe

    e_vec = get_var_data(e_var, times=times)
    ecoro_vec = get_var_data(ecoro_var, at=times)
    evxb_vec = get_var_data(evxb_var, at=times)
    de_vec = e_vec-ecoro_vec-evxb_vec

    store_data, de_var, times, de_vec, limits=limits

end