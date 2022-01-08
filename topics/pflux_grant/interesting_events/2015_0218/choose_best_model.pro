;+
; Return the best model for the given B field and position.
;
; r_var=. var for r_gsm.
; b_var=. var for b_gsm.
;-

function choose_best_model, r_var=r_var, b_var=b_var

    errmsg = ''
    retval = !null
    if tnames(r_var) eq '' then begin
        errmsg = 'No input r_var ...'
        return, retval
    endif

    if tnames(b_var) eq '' then begin
        errmsg = 'No input b_var ...'
        return, retval
    endif

    models = ['t89','t96','t01','t04s']
    prefix = get_prefix(r_var)
    igrf = 1
    foreach model, models do begin
        read_geopack_bfield, r_var, model=model, prefix=prefix, igrf=igrf
    endforeach

    get_data, r_var, times
    b_gsm = get_var_data(b_var, at=times)
    nmodel = n_elements(models)
    diffs = fltarr(nmodel)
    foreach model, models, model_id do begin
        bmod_gsm = get_var_data(prefix+'bmod_gsm_'+model)
        db_gsm = b_gsm-bmod_gsm
        diffs[model_id] = stddev(snorm(db_gsm))
    endforeach

    min_stddev = min(diffs, min_index)
    return, models[min_index]

end


time_range = time_double(['2015-02-18/02:05','2015-02-18/02:15'])
probe = 'a'
rbsp_read_bfield, time_range, probe=probe
rbsp_read_orbit, time_range, probe=probe, coord='gsm'
prefix = 'rbsp'+probe+'_'
the_model = choose_best_model(r_var=prefix+'r_gsm', b_var=prefix+'b_gsm')
end
