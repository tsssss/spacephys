;+
; :Returns: double in [n,3], vector in TSCS.
; :Arguments:
;   vecs: in, required, double. Vector in instrument fixed frame.
;   times: in, optional, double. Dummy.
; :Keywords:
;   probe: in, required, string.
;   errmsg: out, optional, string.
;-

function ct_ts_aci2ts_tscs, vecs, times, probe=probe, errmsg=errmsg

    compile_opt idl2
    return, tracers_instr2ts_tscs(vecs, times, instr_str='aci', probe=probe, errmsg=errmsg)

end
