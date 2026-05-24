;+
; :Returns: double in [n,3], vector in instrument fixed frame.
; :Arguments:
;   vecs: in, required, double. Vector in TSCS.
;   times: in, optional, double. Dummy.
; :Keywords:
;   instr_str: in, required, string. Instrument name.
;   probe: in, required, string.
;   errmsg: out, optional, string.
;-

function tracers_ts_tscs2instr, vecs, times, instr_str=instr_str, probe=probe, errmsg=errmsg

    compile_opt idl2
    errmsg = ''
    retval = !null

    m_instr2tscs = tracers_get_m_instr2tscs(instr_str, probe=probe, errmsg=errmsg)
    if errmsg ne '' then return, retval

    m_tscs2instr = transpose(m_instr2tscs)
    vec_instr = rotate_vector(vecs, m_tscs2instr)
    return, vec_instr

end
