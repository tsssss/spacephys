;+
; Solve 3-sc timing.
;-

function azim_df_solve_triad_timing, tts, rrs

    retval = dictionary()
    ndim = 3

    tt = tts[1:2]-tts[0]
    index = where(tt eq 0, count)
    if count ne 0 then begin
        lprmsg, 'cannot solve triad, dt = 0, skip ...', log_file
        return, retval
    endif

    rr = transpose(rrs[1:2,0:1]-(rrs[0,0:1] ## [1,1]))
    index = where(snorm(transpose(rr)) eq 0, count)
    if count ne 0 then begin
        lprmsg, 'cannot solve triad, dr = 0, skip ...', log_file
        return, retval
    endif

    vv = la_linear_equation(rr,tt)
    vhat = sunitvec(vv)
    rr_normal = dblarr(ndim-1)
    for ii=0, ndim-2 do rr_normal[ii] = sdot(rr[*,ii],vhat)
    fit_result = linfit(tt, rr_normal)
    re = constant('re')
    vmag = fit_result[1]*re
    if ~finite(vmag) then begin
        lprmsg, 'Infinite vmag, skip ...', log_file
        return, retval
    endif

    km_per_s_2_deg_per_min = 1d/re*60*constant('deg')
    center_rr = total(rrs[*,0:1],1)/ndim
    omega_vector = vec_cross([vhat,0], sunitvec([center_rr,0]))*vmag/snorm(center_rr)*km_per_s_2_deg_per_min
    omega = omega_vector[2] ; negative is eastward.

    triad = dictionary($
        'vhat', vhat, $
        'vmag', vmag, $
        'omega', omega)
    return, triad

end
