;+
; Calculate the standard deviation between measured and model B at all s/c.
;-


time_range = time_double(['2014-08-28/10:00','2014-08-28/11:00'])
models = ['t89','t96','t01','t04s']
probes = ['rbspb','tha','thd','the','g13','g15']

_2014_0828_10_load_data

nprobe = n_elements(probes)
nmodel = n_elements(models)
b_diff = dblarr(nprobe,nmodel)
foreach probe, probes, probe_id do begin
    prefix = probe+'_'
    b_meas = get_var_data(prefix+'b_gsm', times=times)
    foreach model, models, model_id do begin
        b_mod = get_var_data(prefix+'bmod_gsm_'+model, at=times)
        db = b_mod-b_meas
        b_diff[probe_id,model_id] = total(snorm(db),/nan)
    endforeach
endforeach

model_b_diff = reform(total(b_diff,1))
tmp = min(model_b_diff, index)
best_model = models[index]
print, best_model
end
