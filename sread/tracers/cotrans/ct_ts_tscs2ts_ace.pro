;+
; :Returns: double in [n,3], vector in ACE fixed frame.
; :Arguments:
;   vecs: in, required, double. Vector in TSCS.
;   times: in, optional, double. Dummy.
; :Keywords:
;   probe: in, required, string.
;   errmsg: out, optional, string.
;-

function ct_ts_tscs2ts_ace, vecs, times, probe=probe, errmsg=errmsg

    compile_opt idl2
    return, tracers_ts_tscs2instr(vecs, times, instr_str='ace', probe=probe, errmsg=errmsg)

end
