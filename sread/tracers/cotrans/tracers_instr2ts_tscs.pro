;+
; :Returns: double in [n,3], vector in TSCS.
; :Arguments:
;   vecs: in, required, double. Vector in instrument fixed frame.
;   times: in, optional, double. Dummy.
; :Keywords:
;   instr_str: in, required, string. Instrument name.
;   probe: in, required, string.
;   errmsg: out, optional, string.
;-

function tracers_instr2ts_tscs, vecs, times, instr_str=instr_str, probe=probe, errmsg=errmsg

    compile_opt idl2
    errmsg = ''
    retval = !null

    m_instr2tscs = tracers_get_m_instr2tscs(instr_str, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval

    vec_tscs = rotate_vector(vecs, m_instr2tscs)
    return, vec_tscs

end
